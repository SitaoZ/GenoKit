# -*- coding: utf-8 -*-
'''
filename: extract_hgvs.py
date: 2026/01/23 12:21
author: Sitao Zhu
mail: zhusitao1990@163.com
Download RefSeq records and visualize DNA-, RNA- and protein-level HGVS variants.
'''


import os
import re
import shutil
import logging
import subprocess
import matplotlib
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from dataclasses import dataclass
from urllib.parse import urlencode
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from typing import List, Optional, Sequence, Tuple, Union




matplotlib.use("Agg")
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42



_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_SVIEWER_URL = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi"
_PROXY_ENVIRONMENT_KEYS = (
    "http_proxy", "https_proxy", "all_proxy",
    "HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY",
)
_POSITION = r"(?:-\d+|\*\d+|\d+)"
_NUCLEOTIDES = r"ACGTUBDHKMNRSVWYacgtubdhkmnrsvwy"
_HGVS_RE = re.compile(
    r"^(?P<accession>[A-Za-z]{1,8}_[0-9]+(?:\.[0-9]+)?):"
    r"(?P<coordinate>[gcnrmop])\.(?P<description>.+)$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class HGVSVariant:
    original: str
    accession: str
    coordinate_type: str
    kind: str
    start_token: str
    end_token: str
    reference: str = ""
    alternate: str = ""
    predicted: bool = False
    repeat_count: int = 0
    repeat_units: Tuple[Tuple[str, int], ...] = ()


@dataclass(frozen=True)
class ProteinVariant:
    original: str
    accession: str
    coordinate_type: str
    kind: str
    position: int
    reference: str
    alternate: str
    predicted: bool = False
    end_position: int = 0
    end_reference: str = ""
    repeat_count: int = 0
    extension_length: Optional[int] = None


@dataclass(frozen=True)
class AlleleVariant:
    """Two or more HGVS changes described on one or two alleles."""
    original: str
    accession: str
    coordinate_type: str
    components: Tuple[Union[HGVSVariant, ProteinVariant], ...]
    groups: Tuple[Tuple[int, ...], ...]
    phase: str
    kind: str = "alleles"


@dataclass(frozen=True)
class LocatedVariant:
    variant: HGVSVariant
    start: int
    end: int
    reference: str
    alternate: str
    warning: str = ""


@dataclass(frozen=True)
class LocatedProteinVariant:
    variant: ProteinVariant
    position: int
    reference: str
    alternate: str
    warning: str = ""
    end_position: int = 0


@dataclass(frozen=True)
class LocatedAlleleVariant:
    variant: AlleleVariant
    components: Tuple[Union[LocatedVariant, LocatedProteinVariant], ...]


@dataclass(frozen=True)
class CDSConsequence:
    codon_number: int
    protein_label: str
    frameshift: bool
    reference_codon: Tuple[str, str, str]
    alternate_codon: Tuple[str, str, str]
    reference_amino_acid: Tuple[str, str, str]
    alternate_amino_acid: Tuple[str, str, str]


_AA_THREE = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "X": "Xaa",
}
_AA_ONE = {three: one for one, three in _AA_THREE.items()}
_AA_ONE.update({"Sec": "U", "Pyl": "O"})
_PROTEIN_SUBSTITUTION_RE = re.compile(
    r"^(?:(?P<ref_three>[A-Z][a-z]{2})|(?P<ref_one>[ARNDCQEGHILKMFPSTWYVXUO]))"
    r"(?P<position>[1-9][0-9]*)"
    r"(?:(?P<alt_three>[A-Z][a-z]{2})|(?P<alt_one>[ARNDCQEGHILKMFPSTWYVXUO*]))$"
)


Variant = Union[HGVSVariant, ProteinVariant, AlleleVariant]
Located = Union[LocatedVariant, LocatedProteinVariant]

_AA_TOKEN = r"(?:[A-Z][a-z]{2}|[ARNDCQEGHILKMFPSTWYVXUO*])"


def _amino_acid(token: str) -> str:
    """Convert one HGVS amino-acid token to its one-letter representation."""
    if len(token) == 1:
        return "*" if token == "*" else token
    value = _AA_ONE.get(token)
    if value is None:
        raise ValueError("Unknown HGVS amino-acid code: %s" % token)
    return value


def _amino_acid_sequence(text: str) -> str:
    """Convert a concatenated one- or three-letter HGVS peptide sequence."""
    if not text:
        return ""
    if re.fullmatch(r"[ARNDCQEGHILKMFPSTWYVXUO*]+", text):
        return text
    tokens = re.findall(r"[A-Z][a-z]{2}", text)
    if not tokens or "".join(tokens) != text:
        raise ValueError("Invalid HGVS amino-acid sequence: %s" % text)
    return "".join(_amino_acid(token) for token in tokens)


def _parse_repeat_description(
    description: str,
) -> Optional[Tuple[str, str, str, str, int, Tuple[Tuple[str, int], ...]]]:
    """Parse HGVS repeated-sequence forms with positions and optional units."""
    match = re.fullmatch(
        rf"(?P<start>{_POSITION})(?:_(?P<end>{_POSITION}))?"
        rf"(?P<body>(?:[{_NUCLEOTIDES}]+\[[1-9][0-9]*\])+|\[[1-9][0-9]*\])",
        description,
    )
    if not match:
        return None
    body = match.group("body")
    if body.startswith("["):
        units: Tuple[Tuple[str, int], ...] = ()
        count = int(body[1:-1])
    else:
        parsed_units = re.findall(
            rf"([{_NUCLEOTIDES}]+)\[([1-9][0-9]*)\]", body
        )
        if not parsed_units or "".join(
            "%s[%s]" % item for item in parsed_units
        ) != body:
            return None
        units = tuple((sequence, int(number)) for sequence, number in parsed_units)
        count = units[0][1] if len(units) == 1 else 0
    start = match.group("start")
    end = match.group("end") or start
    return "rpt", start, end, "", count, units


def _parse_description(description: str) -> Tuple[
    str, str, str, str, str, int, Tuple[Tuple[str, int], ...]
]:
    """Return kind, start, end, reference and alternate from an HGVS suffix."""
    repeated = _parse_repeat_description(description)
    if repeated:
        kind, start, end, reference, count, units = repeated
        return kind, start, end, reference, "", count, units
    patterns = (
        (
            "delins",
            rf"^(?P<start>{_POSITION})(?:_(?P<end>{_POSITION}))?"
            rf"(?:[{_NUCLEOTIDES}]+)?delins(?P<alt>[{_NUCLEOTIDES}]+)$",
        ),
        (
            "sub",
            rf"^(?P<start>{_POSITION})(?P<ref>[{_NUCLEOTIDES}]+)>(?P<alt>[{_NUCLEOTIDES}]+)$",
        ),
        (
            "ins",
            rf"^(?P<start>{_POSITION})_(?P<end>{_POSITION})ins(?P<alt>[{_NUCLEOTIDES}]+)$",
        ),
        (
            "del",
            rf"^(?P<start>{_POSITION})(?:_(?P<end>{_POSITION}))?del(?P<ref>[{_NUCLEOTIDES}]*)$",
        ),
        (
            "dup",
            rf"^(?P<start>{_POSITION})(?:_(?P<end>{_POSITION}))?dup(?P<ref>[{_NUCLEOTIDES}]*)$",
        ),
        (
            "inv",
            rf"^(?P<start>{_POSITION})_(?P<end>{_POSITION})inv$",
        ),
    )
    for kind, pattern in patterns:
        match = re.fullmatch(pattern, description)
        if match:
            values = match.groupdict(default="")
            start = values["start"]
            return (
                kind, start, values.get("end") or start,
                values.get("ref", ""), values.get("alt", ""), 0, (),
            )
    raise ValueError(
        "Unsupported HGVS description %r. Supported operations: substitution, "
        "del, ins, delins, dup, inv and repeated sequence." % description
    )


def _validate_nucleotide_style(
    coordinate_type: str,
    kind: str,
    reference: str,
    alternate: str,
) -> None:
    """Enforce HGVS DNA/RNA case and alphabet conventions."""
    provided = reference + alternate
    if coordinate_type == "r":
        if provided != provided.lower() or "t" in provided:
            raise ValueError(
                "HGVS RNA nucleotides must be lower-case and use u instead of t"
            )
    elif provided != provided.upper() or "U" in provided:
        raise ValueError(
            "HGVS DNA nucleotides must be upper-case and use T instead of U"
        )
    if kind == "sub" and (len(reference) != 1 or len(alternate) != 1):
        raise ValueError(
            "An HGVS substitution replaces exactly one nucleotide with one nucleotide; "
            "use delins for a multi-nucleotide replacement"
        )
    if kind in {"del", "dup"} and reference:
        raise ValueError(
            "Do not append the deleted/duplicated sequence in current HGVS syntax"
        )
    if kind == "delins" and len(reference) == 1 and len(alternate) == 1:
        raise ValueError("A one-nucleotide replacement must be described with >, not delins")


def _parse_protein_description(
    original: str,
    accession: str,
    description: str,
    predicted: bool,
) -> ProteinVariant:
    """Parse HGVS protein sub/del/ins/delins/dup/rpt/ext/fs descriptions."""
    substitution = _PROTEIN_SUBSTITUTION_RE.fullmatch(description)
    if substitution:
        reference = _amino_acid(
            substitution.group("ref_one") or substitution.group("ref_three")
        )
        alternate = _amino_acid(
            substitution.group("alt_one") or substitution.group("alt_three")
        )
        position = int(substitution.group("position"))
        if reference == "*":
            raise ValueError(
                "A protein substitution cannot use a termination codon as its reference"
            )
        if reference == alternate:
            raise ValueError(
                "Use p.%s%d= for an unchanged amino acid, not a substitution"
                % (_AA_THREE.get(reference, reference), position)
            )
        return ProteinVariant(
            original, accession, "p", "sub", position, reference, alternate,
            predicted=predicted,
        )

    frameshift = re.fullmatch(
        rf"(?P<ref>{_AA_TOKEN})(?P<position>[1-9][0-9]*)"
        rf"(?P<alt>{_AA_TOKEN})?fs(?:Ter(?P<distance>[1-9][0-9]*|\?))?",
        description,
    )
    if frameshift:
        distance = frameshift.group("distance")
        return ProteinVariant(
            original, accession, "p", "fs", int(frameshift.group("position")),
            _amino_acid(frameshift.group("ref")),
            _amino_acid(frameshift.group("alt")) if frameshift.group("alt") else "",
            predicted=predicted,
            extension_length=int(distance) if distance and distance.isdigit() else None,
        )

    extension = re.fullmatch(r"Met1ext-(?P<site>[1-9][0-9]*|\?)", description)
    if extension:
        site = extension.group("site")
        return ProteinVariant(
            original, accession, "p", "ext", 1, "M", "",
            predicted=predicted,
            extension_length=-int(site) if site.isdigit() else None,
        )
    extension = re.fullmatch(
        rf"Ter(?P<position>[1-9][0-9]*)(?P<alt>{_AA_TOKEN})extTer"
        rf"(?P<length>[1-9][0-9]*|\?)",
        description,
    )
    if extension:
        length = extension.group("length")
        return ProteinVariant(
            original, accession, "p", "ext", int(extension.group("position")),
            "*", _amino_acid(extension.group("alt")), predicted=predicted,
            extension_length=int(length) if length.isdigit() else None,
        )

    repeated = re.fullmatch(
        rf"(?P<ref1>{_AA_TOKEN})(?P<start>[1-9][0-9]*)"
        rf"(?:_(?P<ref2>{_AA_TOKEN})(?P<end>[1-9][0-9]*))?"
        rf"\[(?P<count>[1-9][0-9]*)\](?:;\[[1-9][0-9]*\])?",
        description,
    )
    if repeated:
        end = int(repeated.group("end") or repeated.group("start"))
        return ProteinVariant(
            original, accession, "p", "rpt", int(repeated.group("start")),
            _amino_acid(repeated.group("ref1")), "", predicted=predicted,
            end_position=end,
            end_reference=(
                _amino_acid(repeated.group("ref2"))
                if repeated.group("ref2") else ""
            ),
            repeat_count=int(repeated.group("count")),
        )

    edit = re.fullmatch(
        rf"(?P<ref1>{_AA_TOKEN})(?P<start>[1-9][0-9]*)"
        rf"(?:_(?P<ref2>{_AA_TOKEN})(?P<end>[1-9][0-9]*))?"
        rf"(?P<kind>delins|del|dup|ins)(?P<alt>(?:{_AA_TOKEN})*)",
        description,
    )
    if edit:
        kind = edit.group("kind")
        start = int(edit.group("start"))
        end = int(edit.group("end") or start)
        alternate_text = edit.group("alt")
        if kind in {"ins", "delins"} and not alternate_text:
            raise ValueError("Protein %s requires an inserted amino-acid sequence" % kind)
        if kind in {"del", "dup"} and alternate_text:
            raise ValueError("Protein %s must not append an amino-acid sequence" % kind)
        if kind == "ins" and (not edit.group("end") or end != start + 1):
            raise ValueError(
                "Protein insertion coordinates must identify two adjacent amino acids"
            )
        return ProteinVariant(
            original, accession, "p", kind, start,
            _amino_acid(edit.group("ref1")),
            _amino_acid_sequence(alternate_text), predicted=predicted,
            end_position=end,
            end_reference=(
                _amino_acid(edit.group("ref2")) if edit.group("ref2") else ""
            ),
        )

    raise ValueError(
        "Unsupported protein HGVS description %r. Supported operations: sub, del, "
        "ins, delins, dup, rpt, ext and fs." % description
    )


def _split_allele_edits(description: str) -> Optional[Tuple[str, List[List[str]]]]:
    """Return phase and edit groups for the three official allele syntaxes."""
    trans = re.fullmatch(r"\[(.+)\];\[(.+)\]", description)
    if trans:
        return "trans", [[trans.group(1)], [trans.group(2)]]
    cis = re.fullmatch(r"\[(.+)\]", description)
    if cis:
        content = cis.group(1)
        group_predicted = content.startswith("(") and content.endswith(")")
        if content.startswith("(") and content.endswith(")"):
            content = content[1:-1]
        edits = content.split(";")
        if len(edits) > 1:
            if group_predicted:
                edits = ["(%s)" % edit for edit in edits]
            return "cis", [edits]
    if "(;)" in description:
        edits = description.split("(;)")
        if len(edits) > 1:
            return "unknown", [[edit] for edit in edits]
    return None


def _parse_alleles(
    original: str, accession: str, coordinate_type: str, description: str
) -> Optional[AlleleVariant]:
    parsed = _split_allele_edits(description)
    if not parsed:
        return None
    phase, edit_groups = parsed
    components: List[Union[HGVSVariant, ProteinVariant]] = []
    groups: List[Tuple[int, ...]] = []
    for group in edit_groups:
        indices = []
        for edit in group:
            edit = edit.strip()
            component_predicted = edit.startswith("(") and edit.endswith(")")
            if component_predicted:
                edit = edit[1:-1]
            component_text = "%s:%s.%s" % (accession, coordinate_type, edit)
            component = (
                _parse_protein_description(
                    component_text, accession, edit, component_predicted
                )
                if coordinate_type == "p"
                else parse_hgvs(component_text)
            )
            if isinstance(component, AlleleVariant):
                raise ValueError("Nested allele descriptions are not supported")
            indices.append(len(components))
            components.append(component)
        groups.append(tuple(indices))
    return AlleleVariant(
        original, accession, coordinate_type, tuple(components), tuple(groups), phase
    )


def parse_hgvs(text: str) -> Variant:
    """Parse a DNA, RNA or protein HGVS expression without silently guessing."""
    text = text.strip()
    if any(character.isspace() for character in text):
        raise ValueError("Spaces are not permitted inside an HGVS description")
    match = _HGVS_RE.fullmatch(text)
    if not match:
        raise ValueError("Invalid HGVS expression: %s" % text)
    coordinate_type = match.group("coordinate").lower()
    if match.group("coordinate") != coordinate_type:
        raise ValueError("The HGVS coordinate prefix must be lower-case")
    if coordinate_type not in {"g", "c", "n", "m", "o", "r", "p"}:
        raise ValueError(
            "Coordinate type %r is not yet supported for structure plotting; "
            "use g./c./n./m./o. for DNA, r. for RNA or p. for protein." % coordinate_type
        )
    if "." not in match.group("accession"):
        raise ValueError("Use an accession.version reference identifier, for example NM_004006.3")
    description = match.group("description")
    allele = _parse_alleles(
        text, match.group("accession"), coordinate_type, description
    )
    if allele:
        return allele
    predicted = description.startswith("(") and description.endswith(")")
    if description.startswith("(") != description.endswith(")"):
        raise ValueError("Unbalanced parentheses in HGVS description")
    if predicted:
        description = description[1:-1]
    if coordinate_type == "p":
        return _parse_protein_description(
            text, match.group("accession"), description, predicted
        )
    kind, start, end, reference, alternate, repeat_count, repeat_units = (
        _parse_description(description)
    )
    _validate_nucleotide_style(coordinate_type, kind, reference, alternate)
    for repeat_sequence, _count in repeat_units:
        _validate_nucleotide_style(coordinate_type, "rpt", repeat_sequence, "")
    if "_" in description and start == end and kind in {"del", "dup", "inv"}:
        raise ValueError("An HGVS range must contain two different positions")
    return HGVSVariant(
        original=text,
        accession=match.group("accession"),
        coordinate_type=coordinate_type,
        kind=kind,
        start_token=start,
        end_token=end,
        reference=reference.upper().replace("U", "T"),
        alternate=alternate.upper().replace("U", "T"),
        predicted=predicted,
        repeat_count=repeat_count,
        repeat_units=tuple(
            (sequence.upper().replace("U", "T"), count)
            for sequence, count in repeat_units
        ),
    )


def _read_hgvs_inputs(args) -> List[str]:
    expressions = []
    if args.string:
        expressions.append(args.string.strip())
    if args.input:
        with open(args.input, encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # Accept a plain one-column file and the first field of CSV/TSV.
                candidate = re.split(r"[\t,]", line, maxsplit=1)[0].strip()
                if candidate.lower() in {"hgvs", "variant", "variation"}:
                    continue
                if ":" not in candidate:
                    raise ValueError("No HGVS expression on input line %d" % line_number)
                expressions.append(candidate)
    if not expressions:
        raise ValueError("Provide one HGVS expression with -s, or a file with -i")
    return expressions


def _request_ncbi_text(url: str, params: dict, use_environment: bool) -> str:
    """Request one NCBI record, optionally bypassing configured proxies."""
    retry = Retry(
        total=0,
        connect=0,
        read=0,
        status=0,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
    )
    session = requests.Session()
    session.trust_env = use_environment
    session.mount("https://", HTTPAdapter(max_retries=retry))
    try:
        response = session.get(
            url,
            params=params,
            headers={"User-Agent": "GenoKit HGVS viewer"},
            timeout=(8, 30),
        )
        response.raise_for_status()
        if not response.text.lstrip().startswith("LOCUS"):
            raise RuntimeError("NCBI response did not contain a flat-file record")
        return response.text
    finally:
        session.close()


def _curl_ncbi_text(url: str, params: dict) -> str:
    """Use the system TLS stack as a final proxy-free NCBI fallback."""
    executable = shutil.which("curl")
    if not executable:
        raise RuntimeError("curl is not installed")
    environment = os.environ.copy()
    for key in _PROXY_ENVIRONMENT_KEYS:
        environment.pop(key, None)
    command = [
        executable,
        "--location",
        "--fail",
        "--silent",
        "--show-error",
        "--noproxy", "*",
        "--connect-timeout", "8",
        "--max-time", "35",
        "%s?%s" % (url, urlencode(params)),
    ]
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=40,
        env=environment,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr.strip() or "curl exited with an error")
    if not completed.stdout.lstrip().startswith("LOCUS"):
        raise RuntimeError("NCBI response did not contain a flat-file record")
    return completed.stdout


def download_genbank(
    accession: str,
    cache_dir: str,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
    protein: bool = False,
) -> Path:
    """Download and cache a nucleotide GenBank or protein GenPept record."""
    cache_path = Path(cache_dir).expanduser() / (accession + (".gp" if protein else ".gb"))
    if cache_path.exists() and cache_path.stat().st_size:
        return cache_path
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    params = {
        "db": "protein" if protein else "nuccore",
        "id": accession,
        "rettype": "gp" if protein else "gbwithparts",
        "retmode": "text",
        "tool": "GenoKit",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    sviewer_params = {
        "db": "protein" if protein else "nuccore",
        "id": accession,
        "report": "genpept" if protein else "genbank",
        "retmode": "text",
    }
    attempts = (
        ("EFetch using configured network settings", _EFETCH_URL, params, True),
        ("EFetch without proxy settings", _EFETCH_URL, params, False),
        ("Sequence Viewer using configured network settings", _SVIEWER_URL, sviewer_params, True),
        ("Sequence Viewer without proxy settings", _SVIEWER_URL, sviewer_params, False),
    )
    errors = []
    record_text = None
    for label, url, request_params, use_environment in attempts:
        try:
            record_text = _request_ncbi_text(url, request_params, use_environment)
            break
        except (requests.RequestException, RuntimeError) as error:
            errors.append("%s: %s" % (label, error))
            logging.debug("%s failed for %s: %s", label, accession, error)
    if record_text is None:
        for label, url, request_params in (
            ("curl EFetch without proxy settings", _EFETCH_URL, params),
            ("curl Sequence Viewer without proxy settings", _SVIEWER_URL, sviewer_params),
        ):
            try:
                record_text = _curl_ncbi_text(url, request_params)
                break
            except (OSError, subprocess.SubprocessError, RuntimeError) as error:
                errors.append("%s: %s" % (label, error))
                logging.debug("%s failed for %s: %s", label, accession, error)
    if record_text is None:
        proxy_hint = ""
        if any(os.environ.get(key) for key in _PROXY_ENVIRONMENT_KEYS):
            proxy_hint = (
                " Proxy environment variables are configured; verify that the local proxy "
                "is running, or temporarily unset HTTP_PROXY/HTTPS_PROXY/http_proxy/https_proxy."
            )
        raise RuntimeError(
            "Unable to download %s from NCBI through EFetch, Sequence Viewer, direct "
            "connections, or curl.%s You can also provide the exact-version record with "
            "--genbank. Last error: %s"
            % (accession, proxy_hint, errors[-1] if errors else "unknown error")
        )
    temporary = cache_path.with_suffix(cache_path.suffix + ".tmp")
    temporary.write_text(record_text, encoding="utf-8")
    temporary.replace(cache_path)
    return cache_path


def load_genbank(accession: str, args, protein: bool = False) -> SeqRecord:
    """Read a user-supplied or cached NCBI GenBank/GenPept record."""
    if args.genbank:
        record_path = Path(args.genbank)
    else:
        record_path = download_genbank(
            accession,
            args.cache_dir,
            email=args.email or os.environ.get("NCBI_EMAIL"),
            api_key=args.api_key or os.environ.get("NCBI_API_KEY"),
            protein=protein,
        )
    with record_path.open(encoding="utf-8") as handle:
        record = SeqIO.read(handle, "genbank")
    requested = accession.upper()
    available = {record.id.upper(), record.name.upper()}
    if requested not in available:
        raise ValueError(
            "%s record %s does not match the exact requested accession.version %s"
            % ("GenPept" if protein else "GenBank", record.id, accession)
        )
    return record


def _cds_bounds(record: SeqRecord) -> Tuple[int, int]:
    cds_features = [feature for feature in record.features if feature.type == "CDS"]
    if not cds_features:
        raise ValueError("The GenBank record has no CDS feature required for c. coordinates")
    return min(int(feature.location.start) for feature in cds_features), max(
        int(feature.location.end) for feature in cds_features
    )


def _map_position(token: str, coordinate_type: str, record: SeqRecord) -> int:
    """Map one HGVS coordinate to a zero-based reference-sequence offset."""
    if coordinate_type in {"g", "n", "m", "o"}:
        if token.startswith(("-", "*")):
            raise ValueError("%s. coordinates must be positive reference positions" % coordinate_type)
        return int(token) - 1

    cds_features = [feature for feature in record.features if feature.type == "CDS"]
    if coordinate_type == "r" and not cds_features:
        # HGVS r. numbering follows n. numbering for a non-coding RNA.
        if token.startswith(("-", "*")):
            raise ValueError(
                "r. coordinates on a non-coding transcript must be positive positions"
            )
        return int(token) - 1

    # On a coding RNA, HGVS r. numbering follows the associated c. numbering:
    # r.1 is the first base of the translation-initiation codon, not the first
    # base in the transcript FASTA record.  The same mapping applies to c.
    cds_start, cds_end = _cds_bounds(record)
    if token.startswith("-"):
        return cds_start - int(token[1:])
    if token.startswith("*"):
        return cds_end + int(token[1:]) - 1
    return cds_start + int(token) - 1


def _display_nucleotides(sequence: str, coordinate_type: str) -> str:
    """Use HGVS RNA letters for r. output while keeping DNA internally."""
    if coordinate_type == "r":
        return sequence.replace("T", "U").lower()
    return sequence.upper()


def _inverse_position(position: int, coordinate_type: str, record: SeqRecord) -> str:
    """Convert a zero-based transcript base back to a simple HGVS coordinate."""
    has_cds = any(feature.type == "CDS" for feature in record.features)
    if coordinate_type in {"g", "n", "m", "o"} or (
        coordinate_type == "r" and not has_cds
    ):
        return str(position + 1)
    cds_start, cds_end = _cds_bounds(record)
    if position < cds_start:
        return "-%d" % (cds_start - position)
    if position >= cds_end:
        return "*%d" % (position - cds_end + 1)
    return str(position - cds_start + 1)


def _exon_boundaries(record: SeqRecord) -> set:
    exons = sorted(
        (int(part.start), int(part.end))
        for feature in record.features if feature.type == "exon"
        for part in feature.location.parts
    )
    return {end for _start, end in exons[:-1]}


def _three_prime_normalization(
    sequence: str,
    start: int,
    end: int,
    kind: str,
    alternate: str,
    record: SeqRecord,
) -> Tuple[int, int, str]:
    """Right-shift unambiguous simple variants without crossing exon junctions.

    Ambiguous reference symbols (for example the leading N gap in chromosome
    records) do not prove sequence identity and therefore must not be used to
    move an HGVS description to a different coordinate.
    """
    boundaries = _exon_boundaries(record)
    if kind in {"del", "dup"}:
        while (
            end < len(sequence)
            and end not in boundaries
            and sequence[start] in "ACGT"
            and sequence[end] in "ACGT"
            and sequence[start] == sequence[end]
        ):
            start += 1
            end += 1
    elif kind == "ins" and alternate:
        while (
            start < len(sequence)
            and start not in boundaries
            and sequence[start] in "ACGT"
            and alternate[0] in "ACGT"
            and sequence[start] == alternate[0]
        ):
            start += 1
            end = start
            alternate = alternate[1:] + alternate[0]
    return start, end, alternate


def _normalized_description(
    variant: HGVSVariant,
    start: int,
    end: int,
    alternate: str,
    record: SeqRecord,
) -> str:
    prefix = "%s:%s." % (variant.accession, variant.coordinate_type)
    if variant.kind == "ins":
        positions = "%s_%s" % (
            _inverse_position(start - 1, variant.coordinate_type, record),
            _inverse_position(start, variant.coordinate_type, record),
        )
        suffix = "ins%s" % _display_nucleotides(alternate, variant.coordinate_type)
    else:
        first = _inverse_position(start, variant.coordinate_type, record)
        last = _inverse_position(end - 1, variant.coordinate_type, record)
        positions = first if first == last else "%s_%s" % (first, last)
        suffix = variant.kind
    description = positions + suffix
    return prefix + ("(%s)" % description if variant.predicted else description)


def locate_variant(
    variant: HGVSVariant,
    record: SeqRecord,
    allow_reference_mismatch: bool = False,
    skip_three_prime_check: bool = False,
) -> LocatedVariant:
    """Map HGVS positions and validate the stated reference allele."""
    sequence = str(record.seq).upper().replace("U", "T")
    has_cds = any(feature.type == "CDS" for feature in record.features)
    if variant.coordinate_type == "n" and has_cds:
        raise ValueError(
            "n. is for a non-coding reference sequence; use c. or r. for this coding transcript"
        )
    start = _map_position(variant.start_token, variant.coordinate_type, record)
    if (
        variant.kind == "rpt"
        and variant.end_token == variant.start_token
        and variant.repeat_units
    ):
        # In forms such as g.123CAG[23], position 123 is the first base of
        # the repeat unit; its length determines the highlighted interval.
        end = start + len(variant.repeat_units[0][0])
    else:
        end = _map_position(variant.end_token, variant.coordinate_type, record) + 1
    if end <= start:
        raise ValueError("HGVS ranges must be listed from 5' to 3'")
    if variant.kind == "ins":
        if end - start != 2:
            raise ValueError("HGVS insertion coordinates must identify two adjacent bases")
        start = start + 1
        end = start
    if not 0 <= start <= end <= len(sequence):
        raise ValueError(
            "Variant coordinates map outside reference %s (length %d)" % (record.id, len(sequence))
        )

    observed = sequence[start:end]
    warning = ""
    if variant.reference and variant.kind not in {"dup"}:
        if observed != variant.reference:
            warning = "Reference mismatch: HGVS=%s, RefSeq=%s" % (
                _display_nucleotides(variant.reference, variant.coordinate_type),
                _display_nucleotides(observed, variant.coordinate_type),
            )
            if not allow_reference_mismatch:
                raise ValueError(
                    "%s: %s. If the mismatch is intentional and only needs "
                    "visualization, add --allow-reference-mismatch"
                    % (variant.original, warning)
                )
            logging.warning("%s: %s", variant.original, warning)

    if variant.kind == "inv" and end - start < 2:
        raise ValueError("An HGVS inversion must contain more than one nucleotide")

    if variant.kind == "rpt" and variant.repeat_units:
        stated_unit = "".join(sequence for sequence, _count in variant.repeat_units)
        if len(variant.repeat_units) == 1 and observed != stated_unit:
            warning = "Reference mismatch: HGVS repeat unit=%s, RefSeq=%s" % (
                _display_nucleotides(stated_unit, variant.coordinate_type),
                _display_nucleotides(observed, variant.coordinate_type),
            )
            if not allow_reference_mismatch:
                raise ValueError(
                    "%s: %s. If the mismatch is intentional and only needs "
                    "visualization, add --allow-reference-mismatch"
                    % (variant.original, warning)
                )
            logging.warning("%s: %s", variant.original, warning)

    normalized_start, normalized_end, normalized_alternate = _three_prime_normalization(
        sequence, start, end, variant.kind, variant.alternate, record
    )
    normalized_duplication_start = None
    if variant.kind == "ins" and normalized_alternate:
        candidate_start = normalized_start - len(normalized_alternate)
        if (
            candidate_start >= 0
            and sequence[candidate_start:normalized_start] == normalized_alternate
        ):
            normalized_duplication_start = candidate_start
    if not skip_three_prime_check and (
        normalized_start != start
        or normalized_end != end
        or normalized_alternate != variant.alternate
    ):
        if normalized_duplication_start is not None:
            first = _inverse_position(
                normalized_duplication_start, variant.coordinate_type, record
            )
            last = _inverse_position(
                normalized_start - 1, variant.coordinate_type, record
            )
            positions = first if first == last else "%s_%s" % (first, last)
            suggestion = "%s:%s.%sdup" % (
                variant.accession, variant.coordinate_type, positions
            )
            reason = "normalizes to a tandem copy and must be described as a duplication"
        else:
            suggestion = _normalized_description(
                variant, normalized_start, normalized_end,
                normalized_alternate, record
            )
            reason = "is not at the most 3' position permitted by HGVS"
        raise ValueError(
            "%s %s; use %s. "
            "To visualize the supplied, non-normalized description unchanged, "
            "add --skip-three-prime-check"
            % (variant.original, reason, suggestion)
        )

    if (
        not skip_three_prime_check
        and variant.kind == "ins"
        and variant.alternate
    ):
        upstream_start = start - len(variant.alternate)
        if upstream_start >= 0 and sequence[upstream_start:start] == variant.alternate:
            first = _inverse_position(
                upstream_start, variant.coordinate_type, record
            )
            last = _inverse_position(
                start - 1, variant.coordinate_type, record
            )
            positions = first if first == last else "%s_%s" % (first, last)
            raise ValueError(
                "This insertion is a tandem copy of the immediately 5' sequence and "
                "must be described as a duplication; use %s:%s.%sdup"
                % (variant.accession, variant.coordinate_type, positions)
            )

    if variant.kind == "sub":
        alternate = variant.alternate
    elif variant.kind == "del":
        alternate = ""
    elif variant.kind == "ins":
        alternate = variant.alternate
    elif variant.kind == "delins":
        alternate = variant.alternate
    elif variant.kind == "dup":
        alternate = observed + observed
    elif variant.kind == "inv":
        alternate = str(Seq(observed).reverse_complement())
    elif variant.kind == "rpt":
        if variant.repeat_units:
            alternate = "".join(
                sequence * count for sequence, count in variant.repeat_units
            )
        else:
            alternate = observed * variant.repeat_count
    else:
        raise ValueError("Unsupported variation type: %s" % variant.kind)
    return LocatedVariant(variant, start, end, observed, alternate, warning)


def locate_protein_variant(
    variant: ProteinVariant,
    record: SeqRecord,
    allow_reference_mismatch: bool = False,
) -> LocatedProteinVariant:
    """Map a protein HGVS edit and validate stated terminal amino acids."""
    sequence = str(record.seq).upper()
    position = variant.position - 1
    if variant.kind == "ext" and variant.reference == "*":
        valid_position = 0 <= position <= len(sequence)
    else:
        valid_position = 0 <= position < len(sequence)
    if not valid_position:
        raise ValueError(
            "Protein position %d is outside %s (length %d aa)"
            % (variant.position, record.id, len(sequence))
        )
    end_position = variant.end_position or variant.position
    maximum_end = (
        len(sequence) + 1
        if variant.kind == "ext" and variant.reference == "*"
        else len(sequence)
    )
    if end_position < variant.position or end_position > maximum_end:
        raise ValueError(
            "Protein range %d-%d is outside %s (length %d aa)"
            % (variant.position, end_position, record.id, len(sequence))
        )
    end_offset = end_position
    if variant.kind == "ins":
        observed = ""
        position = variant.position
        end_offset = position
    elif variant.kind == "ext" and variant.reference == "*":
        observed = "*"
        end_offset = position
    else:
        observed = sequence[position:end_offset]
    warning = ""
    if variant.kind == "ins":
        observed_first = sequence[variant.position - 1:variant.position]
        observed_last = sequence[end_position - 1:end_position]
    else:
        observed_first = observed[:1] or ""
        observed_last = observed[-1:] or ""
    if variant.reference and observed_first != variant.reference:
        warning = "Reference mismatch: HGVS=%s, RefSeq=%s" % (
            _AA_THREE.get(variant.reference, variant.reference),
            _AA_THREE.get(observed_first, observed_first or "empty"),
        )
        if not allow_reference_mismatch:
            raise ValueError("%s: %s" % (variant.original, warning))
        logging.warning("%s: %s", variant.original, warning)
    if variant.end_reference and observed_last != variant.end_reference:
        end_warning = "Range end mismatch: HGVS=%s, RefSeq=%s" % (
            _AA_THREE.get(variant.end_reference, variant.end_reference),
            _AA_THREE.get(observed_last, observed_last or "empty"),
        )
        if not allow_reference_mismatch:
            raise ValueError("%s: %s" % (variant.original, end_warning))
        warning = "; ".join(filter(None, (warning, end_warning)))
        logging.warning("%s: %s", variant.original, end_warning)

    if variant.kind == "sub":
        alternate = variant.alternate
    elif variant.kind == "del":
        alternate = ""
    elif variant.kind == "ins":
        alternate = variant.alternate
    elif variant.kind == "delins":
        alternate = variant.alternate
    elif variant.kind == "dup":
        alternate = observed + observed
    elif variant.kind == "rpt":
        alternate = observed * variant.repeat_count
    elif variant.kind in {"ext", "fs"}:
        # The changed peptide beyond the breakpoint cannot be reconstructed
        # from a protein reference alone; the context panel labels it symbolically.
        alternate = variant.alternate
    else:
        raise ValueError("Unsupported protein variation type: %s" % variant.kind)
    return LocatedProteinVariant(
        variant=variant,
        position=position,
        reference=observed,
        alternate=alternate,
        warning=warning,
        end_position=end_offset,
    )


def locate_allele_variant(
    variant: AlleleVariant,
    record: SeqRecord,
    allow_reference_mismatch: bool = False,
    skip_three_prime_check: bool = False,
) -> LocatedAlleleVariant:
    """Locate every component of an HGVS cis/trans/unknown-phase description."""
    located_components: List[Located] = []
    for component in variant.components:
        if isinstance(component, ProteinVariant):
            located_components.append(
                locate_protein_variant(
                    component, record,
                    allow_reference_mismatch=allow_reference_mismatch,
                )
            )
        else:
            located_components.append(
                locate_variant(
                    component, record,
                    allow_reference_mismatch=allow_reference_mismatch,
                    skip_three_prime_check=skip_three_prime_check,
                )
            )
    return LocatedAlleleVariant(variant, tuple(located_components))


def _feature_intervals(record: SeqRecord, feature_type: str) -> List[Tuple[int, int]]:
    intervals = []
    for feature in record.features:
        if feature.type == feature_type:
            intervals.extend((int(part.start), int(part.end)) for part in feature.location.parts)
    return sorted(set(intervals))


def _gene_label(record: SeqRecord) -> str:
    for feature in record.features:
        if feature.type == "gene" and feature.qualifiers.get("gene"):
            return feature.qualifiers["gene"][0]
    return record.annotations.get("source", "")


def _translation_table(record: SeqRecord) -> int:
    for feature in record.features:
        if feature.type == "CDS" and feature.qualifiers.get("transl_table"):
            return int(feature.qualifiers["transl_table"][0])
    return 1


def _translate(sequence: str, table: int) -> str:
    usable_length = len(sequence) - (len(sequence) % 3)
    if usable_length == 0:
        return ""
    return str(Seq(sequence[:usable_length]).translate(table=table))


def _codon_text(sequence: str) -> str:
    return " ".join(sequence[index:index + 3] for index in range(0, len(sequence), 3))


def _amino_text(sequence: str) -> str:
    # The compact one-letter context remains readable in a monospaced row;
    # the inferred consequence title uses unambiguous three-letter names.
    return sequence


def _surrounding_segments(sequence: str, focus_start: int, unit: int, flank_units: int = 2):
    context_start = max(0, focus_start - flank_units * unit)
    context_start -= context_start % unit
    context_end = min(len(sequence), focus_start + (flank_units + 1) * unit)
    left = sequence[context_start:focus_start]
    focus = sequence[focus_start:min(len(sequence), focus_start + unit)]
    right = sequence[min(len(sequence), focus_start + unit):context_end]
    return left, focus, right


def _with_separators(parts: Tuple[str, str, str], formatter):
    left, focus, right = (formatter(part) for part in parts)
    if left:
        left += " "
    if right:
        right = " " + right
    return left, focus or "∅", right


def cds_consequence(record: SeqRecord, located: LocatedVariant) -> Optional[CDSConsequence]:
    """Calculate local codon and amino-acid contexts for a CDS-contained variant."""
    if located.variant.coordinate_type not in {"c", "r"} or not any(
        feature.type == "CDS" for feature in record.features
    ):
        return None
    cds_start, cds_end = _cds_bounds(record)
    if located.start == located.end:
        inside_cds = cds_start <= located.start < cds_end
    else:
        inside_cds = cds_start <= located.start < located.end <= cds_end
    if not inside_cds:
        return None

    transcript = str(record.seq).upper().replace("U", "T")
    reference_cds = transcript[cds_start:cds_end]
    relative_start = located.start - cds_start
    relative_end = located.end - cds_start
    alternate_cds = (
        reference_cds[:relative_start] + located.alternate + reference_cds[relative_end:]
    )
    codon_start = (relative_start // 3) * 3
    codon_number = codon_start // 3 + 1
    frameshift = (len(located.alternate) - (located.end - located.start)) % 3 != 0

    reference_codon = _with_separators(
        _surrounding_segments(reference_cds, codon_start, 3), _codon_text
    )
    alternate_codon = _with_separators(
        _surrounding_segments(alternate_cds, codon_start, 3), _codon_text
    )

    table = _translation_table(record)
    reference_protein = _translate(reference_cds, table)
    alternate_protein = _translate(alternate_cds, table)
    amino_index = codon_number - 1
    reference_amino_acid = _with_separators(
        _surrounding_segments(reference_protein, amino_index, 1, flank_units=3), _amino_text
    )
    alternate_amino_acid = _with_separators(
        _surrounding_segments(alternate_protein, amino_index, 1, flank_units=3), _amino_text
    )

    reference_aa = reference_protein[amino_index] if amino_index < len(reference_protein) else "X"
    alternate_aa = alternate_protein[amino_index] if amino_index < len(alternate_protein) else "X"
    if frameshift:
        first_changed = amino_index
        comparison_length = min(len(reference_protein), len(alternate_protein))
        while (
            first_changed < comparison_length
            and reference_protein[first_changed] == alternate_protein[first_changed]
        ):
            first_changed += 1
        reference_aa = (
            reference_protein[first_changed] if first_changed < len(reference_protein) else "X"
        )
        alternate_aa = (
            alternate_protein[first_changed] if first_changed < len(alternate_protein) else "X"
        )
        protein_label = "p.(%s%d%sfs)" % (
            _AA_THREE.get(reference_aa, reference_aa),
            first_changed + 1,
            _AA_THREE.get(alternate_aa, alternate_aa),
        )
    elif located.variant.kind == "sub":
        if reference_aa == alternate_aa:
            protein_label = "p.(%s%d=)" % (
                _AA_THREE.get(reference_aa, reference_aa), codon_number
            )
        else:
            protein_label = "p.(%s%d%s)" % (
                _AA_THREE.get(reference_aa, reference_aa),
                codon_number,
                _AA_THREE.get(alternate_aa, alternate_aa),
            )
    else:
        # A formally normalized protein del/dup/ins/delins requires applying
        # protein-level comparison and its own 3' rule.  Keep the visualization
        # informative without emitting a misleading HGVS protein description.
        protein_label = "p.(?)"
    return CDSConsequence(
        codon_number=codon_number,
        protein_label=protein_label,
        frameshift=frameshift,
        reference_codon=reference_codon,
        alternate_codon=alternate_codon,
        reference_amino_acid=reference_amino_acid,
        alternate_amino_acid=alternate_amino_acid,
    )


def _draw_global_axis(ax, record: SeqRecord, located: LocatedVariant) -> None:
    length = len(record.seq)
    exons = _feature_intervals(record, "exon") or [(0, length)]
    cds = _feature_intervals(record, "CDS")
    transcript_level = located.variant.coordinate_type in {"c", "n", "r"}
    show_transcript_utrs = bool(cds) and transcript_level
    five_prime_utr_color = "#FEC84B"
    three_prime_utr_color = "#75E0A7"
    exon_color = "#EAECF0"
    cds_color = "#2E90FA"
    ax.hlines(0.5, 1, length, color="#667085", linewidth=2, zorder=1)
    for index, (start, end) in enumerate(exons):
        ax.add_patch(
            Rectangle(
                (start + 1, 0.33), end - start, 0.34,
                facecolor="#D0D5DD" if index % 2 else exon_color,
                edgecolor="#667085", linewidth=0.8, zorder=2,
            )
        )
    if show_transcript_utrs:
        cds_start = min(start for start, _end in cds)
        cds_end = max(end for _start, end in cds)
        for start, end in exons:
            five_prime_end = min(end, cds_start)
            if start < five_prime_end:
                ax.add_patch(
                    Rectangle(
                        (start + 1, 0.33), five_prime_end - start, 0.34,
                        facecolor=five_prime_utr_color, edgecolor="#B54708",
                        linewidth=0.8, zorder=3,
                    )
                )
            three_prime_start = max(start, cds_end)
            if three_prime_start < end:
                ax.add_patch(
                    Rectangle(
                        (three_prime_start + 1, 0.33), end - three_prime_start, 0.34,
                        facecolor=three_prime_utr_color, edgecolor="#067647",
                        linewidth=0.8, zorder=3,
                    )
                )
    for start, end in cds:
        ax.add_patch(
            Rectangle(
                (start + 1, 0.27), end - start, 0.46,
                facecolor=cds_color, edgecolor="#175CD3", linewidth=0.8, zorder=3,
            )
        )
    # Keep exon junctions visible even where the CDS overlay spans them.
    exon_boundaries = sorted({end for _start, end in exons[:-1]})
    for boundary in exon_boundaries:
        ax.vlines(boundary + 1, 0.27, 0.73, color="#475467", linewidth=1.0, zorder=4)
    marker = located.start + 0.5 if located.start == located.end else (located.start + located.end + 1) / 2
    ax.annotate(
        located.variant.original,
        xy=(marker, 0.76), xytext=(marker, 1.18),
        ha="center", va="bottom", fontsize=9, color="#D92D20",
        arrowprops={"arrowstyle": "-|>", "color": "#D92D20", "lw": 1.5},
    )
    ax.set_xlim(1, max(2, length, marker))
    ax.set_ylim(0, 1.55)
    ax.set_yticks([])
    axis_name = "Transcript" if transcript_level else "Reference sequence"
    if transcript_level:
        ax.annotate(
            "",
            xy=(length, 0.13), xytext=(1, 0.13),
            arrowprops={"arrowstyle": "-|>", "color": "#475467", "lw": 1.1},
        )
        ax.text(1, 0.08, "5′", ha="left", va="top", fontsize=9, color="#344054")
        ax.text(length, 0.08, "3′", ha="right", va="top", fontsize=9, color="#344054")
    ax.set_xlabel("%s coordinate (1-based)" % axis_name)
    ax.set_title(
        "Global %s structure: %s  %s"
        % (axis_name.lower(), record.id, _gene_label(record)),
        loc="left",
    )
    ax.spines[["left", "right", "top"]].set_visible(False)
    if show_transcript_utrs:
        legend_handles = [
            Patch(facecolor=five_prime_utr_color, edgecolor="#B54708", label="5′ UTR"),
            Patch(facecolor=cds_color, edgecolor="#175CD3", label="CDS"),
            Patch(facecolor=three_prime_utr_color, edgecolor="#067647", label="3′ UTR"),
            Patch(facecolor="#D92D20", label="Variant"),
        ]
    elif cds:
        legend_handles = [
            Patch(facecolor=exon_color, edgecolor="#667085", label="Exon"),
            Patch(facecolor=cds_color, edgecolor="#175CD3", label="CDS"),
            Patch(facecolor="#D92D20", label="Variant"),
        ]
    else:
        legend_handles = [
            Patch(facecolor=exon_color, edgecolor="#667085", label="Exon"),
            Patch(facecolor="#D92D20", label="Variant"),
        ]
    ax.legend(
        handles=legend_handles,
        loc="upper right", frameon=False, ncol=len(legend_handles),
    )


def _context_segments(record: SeqRecord, located: LocatedVariant, flank: int):
    sequence = str(record.seq).upper().replace("U", "T")
    left_bound = max(0, located.start - flank)
    right_bound = min(len(sequence), located.end + flank)
    coordinate_type = located.variant.coordinate_type
    left = _display_nucleotides(sequence[left_bound:located.start], coordinate_type)
    right = _display_nucleotides(sequence[located.end:right_bound], coordinate_type)
    reference = _compact_nucleotide_sequence(
        located.reference, coordinate_type
    )
    alternate = _compact_nucleotide_sequence(
        located.alternate, coordinate_type
    )
    return left_bound, right_bound, left, reference, alternate, right


def _compact_nucleotide_sequence(
    sequence: str, coordinate_type: str, maximum: int = 30
) -> str:
    """Keep long indel/repeat alleles legible in the local context panel."""
    displayed = _display_nucleotides(sequence, coordinate_type)
    if not displayed:
        return "∅"
    if len(displayed) <= maximum:
        return displayed
    return "%s...(%d nt)...%s" % (
        displayed[:8], len(displayed), displayed[-8:]
    )


def _draw_sequence_row(ax, y, label, left, allele, right, allele_color, fontsize=11):
    display = left + allele + right
    allele_start = len(left)
    allele_end = allele_start + len(allele)
    for index, base in enumerate(display):
        color = allele_color if allele_start <= index < allele_end else "#344054"
        weight = "bold" if allele_start <= index < allele_end else "normal"
        ax.text(index, y, base, ha="center", va="center", family="monospace", fontsize=fontsize,
                color=color, fontweight=weight)
    ax.text(-2.0, y, label, ha="right", va="center", fontsize=10, fontweight="bold")
    return len(display)


def _position_label(located: LocatedVariant) -> str:
    coordinate = "%s.%s" % (located.variant.coordinate_type, located.variant.start_token)
    if located.variant.end_token != located.variant.start_token:
        coordinate += "_%s" % located.variant.end_token
    if located.start == located.end:
        reference_position = "between reference positions %d and %d" % (
            located.start,
            located.start + 1,
        )
    elif located.end - located.start == 1:
        reference_position = "reference position %d" % (located.start + 1)
    else:
        reference_position = "reference positions %d-%d" % (located.start + 1, located.end)
    return "%s  |  %s" % (coordinate, reference_position)


def _draw_context_axis(
    ax,
    record: SeqRecord,
    located: LocatedVariant,
    flank: int,
    consequence: Optional[CDSConsequence],
) -> None:
    left_bound, right_bound, left, reference, alternate, right = _context_segments(record, located, flank)
    if consequence:
        reference_y, alternate_y, marker_y = 5.35, 4.55, 3.88
    else:
        reference_y, alternate_y, marker_y = 2.15, 1.35, 0.68
    reference_width = _draw_sequence_row(
        ax, reference_y, "Nucleotide Ref", left, reference, right, "#D92D20"
    )
    alternate_width = _draw_sequence_row(
        ax, alternate_y, "Nucleotide Alt", left, alternate, right, "#039855"
    )
    marker_x = len(left) + max(0, len(reference) - 1) / 2
    ax.annotate(
        _position_label(located),
        xy=(marker_x, reference_y - 0.14),
        xytext=(marker_x, marker_y),
        ha="center", va="top", fontsize=9, color="#B42318",
        arrowprops={"arrowstyle": "-|>", "color": "#D92D20", "lw": 1.1},
    )

    widths = [reference_width, alternate_width]
    if consequence:
        if located.variant.coordinate_type == "r":
            reference_codon = tuple(
                _display_nucleotides(part, "r") for part in consequence.reference_codon
            )
            alternate_codon = tuple(
                _display_nucleotides(part, "r") for part in consequence.alternate_codon
            )
        else:
            reference_codon = consequence.reference_codon
            alternate_codon = consequence.alternate_codon
        ax.text(
            0, 3.25,
            "CDS consequence: codon %d  |  %s%s"
            % (
                consequence.codon_number,
                consequence.protein_label,
                "  |  frameshift" if consequence.frameshift else "",
            ),
            ha="left", va="center", fontsize=10, fontweight="bold", color="#175CD3",
        )
        codon_reference_width = _draw_sequence_row(
            ax, 2.55, "Codon Ref", *reference_codon, "#D92D20", fontsize=10
        )
        codon_alternate_width = _draw_sequence_row(
            ax, 1.90, "Codon Alt", *alternate_codon, "#039855", fontsize=10
        )
        amino_reference_width = _draw_sequence_row(
            ax, 1.05, "Protein Ref", *consequence.reference_amino_acid, "#D92D20", fontsize=9
        )
        amino_alternate_width = _draw_sequence_row(
            ax, 0.38, "Protein Alt", *consequence.alternate_amino_acid, "#039855", fontsize=9
        )
        widths.extend(
            [codon_reference_width, codon_alternate_width, amino_reference_width, amino_alternate_width]
        )

    ax.set_xlim(-max(16, len("Nucleotide Reference") + 2), max(widths))
    ax.set_ylim((-0.35, 6.05) if consequence else (-0.35, 2.85))
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title(
        "Local sequence context: reference %d-%d | %s %s>%s"
        % (
            left_bound + 1,
            right_bound,
            located.variant.kind,
            reference,
            alternate,
        ),
        loc="left",
    )
    if located.warning:
        ax.text(0, -0.12, located.warning, color="#D92D20", fontsize=9, ha="left", va="top")
    for spine in ax.spines.values():
        spine.set_visible(False)


def _protein_feature_label(feature) -> str:
    for key in ("region_name", "product", "name", "function", "note"):
        values = feature.qualifiers.get(key)
        if values:
            return str(values[0]).split(";")[0]
    return feature.type


def _draw_protein_global_axis(
    ax,
    record: SeqRecord,
    located: LocatedProteinVariant,
) -> None:
    length = len(record.seq)
    ax.hlines(0.5, 1, max(2, length), color="#667085", linewidth=5, zorder=1)
    # GenPept ``FEATURES`` can contain many feature types (for example
    # ``Region``, ``Site``, ``Bond`` and ``mat_peptide``).  The global protein
    # overview deliberately visualizes only records whose feature key is
    # ``Region``; other feature types are annotations of a different semantic
    # level and previously made the overview unnecessarily dense.
    features = [
        feature for feature in record.features
        if feature.type.lower() == "region"
        and int(feature.location.end) > int(feature.location.start)
    ]
    region_color = "#84CAFF"
    labeled = 0
    for index, feature in enumerate(features):
        start = int(feature.location.start)
        end = int(feature.location.end)
        y = 0.31 + 0.12 * (index % 3)
        ax.add_patch(
            Rectangle(
                (start + 1, y), max(1, end - start), 0.26,
                facecolor=region_color, edgecolor="#475467",
                linewidth=0.7, alpha=0.9, zorder=2,
            )
        )
        if labeled < 8 and (end - start >= max(8, length * 0.025)):
            center = (start + end + 1) / 2
            # Distribute labels over three vertical levels and connect each
            # label to the upper edge of its Region block.  The connector is
            # particularly useful when several RefSeq regions overlap.
            label_y = 0.84 + 0.13 * (labeled % 3)
            ax.annotate(
                _protein_feature_label(feature)[:30],
                xy=(center, y + 0.26),
                xytext=(center, label_y),
                ha="center", va="bottom", fontsize=7, color="#344054",
                arrowprops={
                    "arrowstyle": "-",
                    "color": "#667085",
                    "linewidth": 0.7,
                    "shrinkA": 2,
                    "shrinkB": 0,
                },
                annotation_clip=True,
                zorder=4,
            )
            labeled += 1
    marker = (
        located.position + 1
        if located.end_position <= located.position + 1
        else (located.position + located.end_position + 1) / 2
    )
    ax.annotate(
        located.variant.original,
        xy=(marker, 0.75), xytext=(marker, 1.25),
        ha="center", va="bottom", fontsize=9, color="#D92D20",
        arrowprops={"arrowstyle": "-|>", "color": "#D92D20", "lw": 1.5},
    )
    ax.set_xlim(1, max(2, length, marker))
    ax.set_ylim(0, 1.62)
    ax.set_yticks([])
    ax.set_xlabel("Protein coordinate (amino acid, 1-based)")
    ax.set_title("Global protein structure: %s  %s" % (record.id, _gene_label(record)), loc="left")
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.legend(
        handles=[
            Patch(facecolor=region_color, edgecolor="#475467", label="Annotated region"),
            Patch(facecolor="#D92D20", label="Variant"),
        ],
        loc="upper right", frameon=False, ncol=2,
    )


def _draw_protein_context_axis(
    ax,
    record: SeqRecord,
    located: LocatedProteinVariant,
    flank: int,
) -> None:
    sequence = str(record.seq).upper()
    left_bound = max(0, located.position - flank)
    reference_end = max(located.position, located.end_position)
    right_bound = min(len(sequence), reference_end + flank)
    left = sequence[left_bound:located.position]
    right = sequence[reference_end:right_bound]
    reference = located.reference or "∅"
    alternate = located.alternate or "∅"
    alternate_right = right
    if located.variant.kind == "fs":
        alternate = (
            (_AA_THREE.get(located.variant.alternate, located.variant.alternate) or "?")
            + "fs"
            + (
                "Ter%d" % located.variant.extension_length
                if located.variant.extension_length is not None else ""
            )
        )
        alternate_right = ""
    elif located.variant.kind == "ext":
        if located.variant.reference == "M":
            alternate = "EXT(%s)+M" % (
                located.variant.extension_length
                if located.variant.extension_length is not None else "?"
            )
        else:
            alternate = "%sextTer%s" % (
                _AA_THREE.get(located.variant.alternate, located.variant.alternate),
                located.variant.extension_length
                if located.variant.extension_length is not None else "?",
            )
    reference_width = _draw_sequence_row(
        ax, 2.15, "Protein Ref", left, reference, right, "#D92D20"
    )
    alternate_width = _draw_sequence_row(
        ax, 1.35, "Protein Alt", left, alternate, alternate_right, "#039855"
    )
    marker_x = len(left)
    ax.annotate(
        "protein position %d  |  %s" % (
            located.variant.position, located.variant.original.split(":", 1)[1]
        ),
        xy=(marker_x, 2.0), xytext=(marker_x, 0.68),
        ha="center", va="top", fontsize=9, color="#B42318",
        arrowprops={"arrowstyle": "-|>", "color": "#D92D20", "lw": 1.1},
    )
    ax.set_xlim(-14, max(reference_width, alternate_width))
    ax.set_ylim(-0.35, 2.85)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title(
        "Local amino-acid context: residues %d-%d | %s | %s > %s"
        % (
            left_bound + 1,
            right_bound,
            located.variant.kind,
            located.reference or "empty",
            alternate,
        ),
        loc="left",
    )
    if located.warning:
        ax.text(0, -0.12, located.warning, color="#D92D20", fontsize=9, ha="left", va="top")
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_protein_hgvs(
    record: SeqRecord,
    located: LocatedProteinVariant,
    output: Path,
    context: int,
    dpi: int,
) -> Path:
    """Create a two-level protein structure and amino-acid context figure."""
    context_width = min(18, max(11, 0.22 * (context * 2 + 1)))
    figure, axes = plt.subplots(
        2, 1, figsize=(context_width, 7.6),
        gridspec_kw={"height_ratios": [1.25, 1], "hspace": 0.45},
    )
    _draw_protein_global_axis(axes[0], record, located)
    _draw_protein_context_axis(axes[1], record, located, context)
    figure.suptitle("GenoKit HGVS Protein Variant Viewer", fontsize=15, fontweight="bold", y=0.98)
    figure.text(
        0.01, 0.01,
        "Reference annotation: NCBI RefSeq GenPept record",
        fontsize=8, color="#667085",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(figure)
    return output


def plot_hgvs(record: SeqRecord, located: LocatedVariant, output: Path, context: int, dpi: int) -> Path:
    """Create one two-level HGVS structure/context figure."""
    consequence = cds_consequence(record, located)
    context_width = min(18, max(11, 0.22 * (context * 2 + max(1, len(located.alternate)))))
    figure, axes = plt.subplots(
        2, 1, figsize=(context_width, 10.2 if consequence else 7.6),
        gridspec_kw={
            "height_ratios": [1.15, 2.0] if consequence else [1.25, 1],
            "hspace": 0.45,
        },
    )
    _draw_global_axis(axes[0], record, located)
    _draw_context_axis(axes[1], record, located, context, consequence)
    figure.suptitle("GenoKit HGVS Variant Viewer", fontsize=15, fontweight="bold", y=0.98)
    figure.text(
        0.01, 0.01,
        "Reference annotation: NCBI RefSeq GenBank record",
        fontsize=8, color="#667085",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(figure)
    return output


def _add_global_component_markers(ax, components: Sequence[Located]) -> None:
    """Add markers after the first component to a shared allele overview."""
    for index, component in enumerate(components[1:], start=1):
        if isinstance(component, LocatedProteinVariant):
            marker = (
                component.position + 1
                if component.end_position <= component.position + 1
                else (component.position + component.end_position + 1) / 2
            )
        else:
            marker = (
                component.start + 0.5
                if component.start == component.end
                else (component.start + component.end + 1) / 2
            )
        ax.annotate(
            component.variant.original.split(":", 1)[-1],
            xy=(marker, 0.76), xytext=(marker, 1.18 + 0.18 * (index % 2)),
            ha="center", va="bottom", fontsize=8, color="#D92D20",
            arrowprops={"arrowstyle": "-|>", "color": "#D92D20", "lw": 1.2},
        )


def plot_allele_hgvs(
    record: SeqRecord,
    located: LocatedAlleleVariant,
    output: Path,
    context: int,
    dpi: int,
) -> Path:
    """Visualize all components of a cis, trans or unknown-phase HGVS allele."""
    components = located.components
    if not components:
        raise ValueError("An allele description must contain at least two changes")
    protein = isinstance(components[0], LocatedProteinVariant)
    figure_width = min(18, max(12, 0.22 * (context * 2 + 1)))
    figure, axes = plt.subplots(
        len(components) + 1, 1,
        figsize=(figure_width, 3.7 + 3.8 * len(components)),
        gridspec_kw={"hspace": 0.68},
    )
    if protein:
        _draw_protein_global_axis(axes[0], record, components[0])
    else:
        _draw_global_axis(axes[0], record, components[0])
    # The single-variant overview uses the full expression. In an allele
    # overview a compact component suffix avoids collisions between markers.
    for text_artist in axes[0].texts:
        if hasattr(text_artist, "xy"):
            label = text_artist.get_text()
            if ":" in label:
                text_artist.set_text(label.split(":", 1)[-1])
                text_artist.set_position((text_artist.xy[0], 1.18))
            break
    _add_global_component_markers(axes[0], components)
    for axis, component in zip(axes[1:], components):
        if isinstance(component, LocatedProteinVariant):
            _draw_protein_context_axis(axis, record, component, context)
        else:
            _draw_context_axis(
                axis, record, component, context,
                cds_consequence(record, component),
            )
    figure.suptitle(
        "GenoKit HGVS Allele Viewer | %s phase\n%s"
        % (located.variant.phase, located.variant.original),
        fontsize=14, fontweight="bold", y=0.975,
    )
    figure.text(
        0.01, 0.01,
        "Reference annotation: NCBI RefSeq %s record"
        % ("GenPept" if protein else "GenBank"),
        fontsize=8, color="#667085",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(figure)
    return output


def _output_paths(output: Optional[str], variants: Sequence[Variant]) -> List[Path]:
    if output and len(variants) == 1:
        path = Path(output)
        return [path if path.suffix.lower() in {".png", ".pdf", ".svg"} else path.with_suffix(".png")]
    if output:
        requested = Path(output)
        directory = requested if not requested.suffix else requested.parent
        prefix = "hgvs" if not requested.suffix else requested.stem
        suffix = requested.suffix.lower() if requested.suffix.lower() in {".png", ".pdf", ".svg"} else ".png"
    else:
        directory = Path.cwd()
        prefix = "hgvs"
        suffix = ".png"
    paths = []
    for index, variant in enumerate(variants, start=1):
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", variant.original)
        paths.append(directory / ("%s_%02d_%s%s" % (prefix, index, safe_name, suffix)))
    return paths


def get_hgvs(args):
    """CLI entry point for DNA, RNA and protein HGVS visualization."""

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S",
                        filename="HGVS.log"
                        )
    if args.context < 0:
        raise ValueError("--context must be zero or greater")
    if args.dpi < 50:
        raise ValueError("--dpi must be at least 50")
    variants = [parse_hgvs(item) for item in _read_hgvs_inputs(args)]
    if args.variation not in {None, "other"} and any(
        item.kind != args.variation for item in variants
    ):
        raise ValueError("--variation does not match the supplied HGVS expression")
    outputs = _output_paths(args.output, variants)
    record_cache = {}
    try:
        for variant, output in zip(variants, outputs):
            protein = variant.coordinate_type == "p"
            cache_key = (variant.accession, protein)
            if cache_key not in record_cache:
                record_cache[cache_key] = load_genbank(
                    variant.accession, args, protein=protein
                )
            record = record_cache[cache_key]
            if isinstance(variant, AlleleVariant):
                located_allele = locate_allele_variant(
                    variant,
                    record,
                    allow_reference_mismatch=args.allow_reference_mismatch,
                    skip_three_prime_check=args.skip_three_prime_check,
                )
                plot_allele_hgvs(
                    record, located_allele, output, args.context, args.dpi
                )
            elif protein:
                located_protein = locate_protein_variant(
                    variant,
                    record,
                    allow_reference_mismatch=args.allow_reference_mismatch,
                )
                plot_protein_hgvs(
                    record, located_protein, output, args.context, args.dpi
                )
            else:
                located = locate_variant(
                    variant,
                    record,
                    allow_reference_mismatch=args.allow_reference_mismatch,
                    skip_three_prime_check=args.skip_three_prime_check,
                )
                plot_hgvs(record, located, output, args.context, args.dpi)
            print("HGVS figure written: %s" % output)
    except (OSError, RuntimeError, ValueError) as error:
        raise SystemExit("GenoKit hgvs: error: %s" % error) from None
    return outputs


if __name__ == "__main__":
    # Parser-independent example for developing the HGVS grammar.
    example = parse_hgvs("NM_004006.3:r.123c>g")
    print(example)
