#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SQLite-backed annotation database with legacy-pickle compatibility.

The public reader intentionally behaves like GenoKit's historical nested
``dict`` database::

    db[gene_id]["gene"]
    db[gene_id][transcript_id]["mRNA"]
    db[gene_id][transcript_id]["exon"]
    t2g[transcript_id]

For SQLite databases only the requested gene model is materialized.  Whole-
genome ``items()``/``values()`` iteration uses one ordered SQL cursor and keeps
only one gene model alive at a time.  Existing pickle databases are returned
unchanged, so old files remain usable by all migrated command modules.
"""

from __future__ import annotations

import os
import re
import pickle
import sqlite3
from pathlib import Path
from collections import OrderedDict
from collections.abc import Mapping
from typing import Iterator, Optional, Tuple

from tqdm import tqdm

from GenoKit.utils.util import gene_symbol, open_auto, record, utr3_type, utr5_type


SQLITE_MAGIC = b"SQLite format 3\x00"
# Version 1 was also used by an earlier, structurally different SQLite layout.
# New databases use version 3; readers still inspect columns for compatibility
# with already-created version-1 files of either layout.
SCHEMA_VERSION = "3"
RECORD_COLUMNS = (
    "seqid", "source", "feature_type", "start", "end", "score",
    "strand", "phase", "attribute",
)
_GTF_GENE_ID = 'gene_id "'
_GTF_TRANSCRIPT_ID = 'transcript_id "'
_GFF_ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
_GFF_PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")


def is_sqlite_database(path: os.PathLike | str) -> bool:
    """Return ``True`` when *path* has the SQLite file signature."""
    try:
        with open(path, "rb") as handle:
            return handle.read(len(SQLITE_MAGIC)) == SQLITE_MAGIC
    except OSError:
        return False


def detect_annotation_style(
    path: os.PathLike | str,
    max_records: int = 200,
) -> Optional[str]:
    """Identify GTF or GFF3 from headers and ninth-column attributes.

    File suffixes are intentionally ignored because annotation files are often
    renamed or gzip-compressed.  ``None`` means that the inspected records did
    not contain enough format-specific evidence.
    """
    gtf_signals = 0
    gff_signals = 0
    inspected = 0
    with open_auto(path, "rb") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.lower().startswith(b"##gff-version"):
                return "gff"
            if stripped.startswith(b"#"):
                continue
            fields = raw_line.decode("utf-8").rstrip("\r\n").split("\t", 8)
            if len(fields) != 9:
                continue
            attributes = fields[8]
            if re.search(r'(?:^|;)\s*(?:gene_id|transcript_id)\s+"', attributes):
                gtf_signals += 1
            if re.search(r"(?:^|;)\s*(?:ID|Parent)=", attributes):
                gff_signals += 1
            inspected += 1
            if inspected >= max_records:
                break
    if gtf_signals and not gff_signals:
        return "gtf"
    if gff_signals and not gtf_signals:
        return "gff"
    if gtf_signals != gff_signals:
        return "gtf" if gtf_signals > gff_signals else "gff"
    return None


def _record_from_values(values) -> record:
    return record(
        values[0], values[1], values[2], int(values[3]), int(values[4]),
        values[5], values[6], values[7], values[8],
    )


def _record_values(item: record) -> tuple:
    return (
        item.chr, item.source, item.type, int(item.start), int(item.end),
        item.score, item.strand, item.phase, item.attribute,
    )


def _quoted_gtf_attribute(attributes: str, prefix: str) -> str:
    start = attributes.find(prefix)
    if start < 0:
        raise ValueError("Missing %s attribute" % prefix[:-2])
    start += len(prefix)
    end = attributes.find('"', start)
    if end < 0:
        raise ValueError("Unterminated %s attribute" % prefix[:-2])
    return attributes[start:end]


def _gff_attribute(attributes: str, pattern: re.Pattern, name: str) -> str:
    match = pattern.search(attributes)
    if not match:
        raise ValueError("Missing %s attribute in GFF record" % name)
    return match.group(1)


def _annotation_size(path: os.PathLike | str) -> Optional[int]:
    """Return a meaningful byte total for an uncompressed annotation.

    ``gzip.open`` reports decompressed bytes while ``getsize`` reports the
    compressed file size.  An indeterminate progress bar is therefore more
    accurate for ``.gz`` input than a bar that exceeds 100 percent.
    """
    if os.fspath(path).lower().endswith(".gz"):
        return None
    return os.path.getsize(path)


def _connect_for_build(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(str(path))
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA locking_mode=EXCLUSIVE")
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.execute("PRAGMA cache_size=-65536")
    connection.execute("PRAGMA page_size=16384")
    return connection


def _create_schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        ) WITHOUT ROWID;

        CREATE TABLE genes (
            gene_id TEXT PRIMARY KEY,
            symbol TEXT NOT NULL,
            ordinal INTEGER NOT NULL,
            seqid TEXT NOT NULL,
            source TEXT NOT NULL,
            feature_type TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            score TEXT NOT NULL,
            strand TEXT NOT NULL,
            phase TEXT NOT NULL,
            attribute TEXT NOT NULL,
            explicit_record INTEGER NOT NULL DEFAULT 0
        ) WITHOUT ROWID;

        CREATE TABLE transcripts (
            transcript_id TEXT PRIMARY KEY,
            gene_id TEXT NOT NULL,
            ordinal INTEGER NOT NULL,
            seqid TEXT NOT NULL,
            source TEXT NOT NULL,
            feature_type TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            score TEXT NOT NULL,
            strand TEXT NOT NULL,
            phase TEXT NOT NULL,
            attribute TEXT NOT NULL
        ) WITHOUT ROWID;

        CREATE TABLE features (
            feature_id INTEGER PRIMARY KEY,
            transcript_id TEXT NOT NULL,
            feature_key TEXT NOT NULL,
            ordinal INTEGER NOT NULL,
            seqid TEXT NOT NULL,
            source TEXT NOT NULL,
            feature_type TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            score TEXT NOT NULL,
            strand TEXT NOT NULL,
            phase TEXT NOT NULL,
            attribute TEXT NOT NULL
        );
        """
    )


_GENE_UPSERT = """
INSERT INTO genes(
    gene_id, symbol, ordinal, seqid, source, feature_type, start, end, score,
    strand, phase, attribute, explicit_record
) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT(gene_id) DO UPDATE SET
    symbol = CASE
        WHEN excluded.explicit_record = 1 THEN excluded.symbol
        WHEN genes.symbol = genes.gene_id AND excluded.symbol != excluded.gene_id
            THEN excluded.symbol
        ELSE genes.symbol END,
    ordinal = MIN(genes.ordinal, excluded.ordinal),
    seqid = CASE WHEN excluded.explicit_record = 1 THEN excluded.seqid ELSE genes.seqid END,
    source = CASE WHEN excluded.explicit_record = 1 THEN excluded.source ELSE genes.source END,
    feature_type = CASE WHEN excluded.explicit_record = 1 THEN excluded.feature_type ELSE genes.feature_type END,
    start = CASE
        WHEN excluded.explicit_record = 1 THEN excluded.start
        WHEN genes.explicit_record = 0 THEN MIN(genes.start, excluded.start)
        ELSE genes.start END,
    end = CASE
        WHEN excluded.explicit_record = 1 THEN excluded.end
        WHEN genes.explicit_record = 0 THEN MAX(genes.end, excluded.end)
        ELSE genes.end END,
    score = CASE WHEN excluded.explicit_record = 1 THEN excluded.score ELSE genes.score END,
    strand = CASE WHEN excluded.explicit_record = 1 THEN excluded.strand ELSE genes.strand END,
    phase = CASE WHEN excluded.explicit_record = 1 THEN excluded.phase ELSE genes.phase END,
    attribute = CASE WHEN excluded.explicit_record = 1 THEN excluded.attribute ELSE genes.attribute END,
    explicit_record = MAX(genes.explicit_record, excluded.explicit_record)
"""


_TRANSCRIPT_UPSERT = """
INSERT INTO transcripts(
    transcript_id, gene_id, ordinal, seqid, source, feature_type, start,
    end, score, strand, phase, attribute
) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT(transcript_id) DO UPDATE SET
    gene_id=excluded.gene_id,
    ordinal=MIN(transcripts.ordinal, excluded.ordinal),
    seqid=excluded.seqid,
    source=excluded.source,
    feature_type=excluded.feature_type,
    start=excluded.start,
    end=excluded.end,
    score=excluded.score,
    strand=excluded.strand,
    phase=excluded.phase,
    attribute=excluded.attribute
"""


_FEATURE_INSERT = """
INSERT INTO features(
    transcript_id, feature_key, ordinal, seqid, source, feature_type,
    start, end, score, strand, phase, attribute
) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


def _placeholder_gene_values(
    gene_id: str,
    ordinal: int,
    fields: list,
    symbol: Optional[str] = None,
) -> tuple:
    symbol = symbol or gene_id
    attribute = 'gene_id "%s";' % gene_id
    if symbol != gene_id:
        attribute += ' gene_name "%s";' % symbol
    return (
        gene_id, symbol, ordinal,
        fields[0], fields[1], "gene", int(fields[3]),
        int(fields[4]), "", fields[6], ".", attribute, 0,
    )


def _explicit_gene_values(gene_id: str, ordinal: int, fields: list) -> tuple:
    return (
        gene_id, gene_symbol(fields[8], gene_id), ordinal,
        fields[0], fields[1], fields[2], int(fields[3]),
        int(fields[4]), fields[5], fields[6], fields[7], fields[8], 1,
    )


def _transcript_values(transcript_id: str, gene_id: str, ordinal: int, fields: list) -> tuple:
    return (
        transcript_id, gene_id, ordinal, fields[0], fields[1], fields[2],
        int(fields[3]), int(fields[4]), fields[5], fields[6], fields[7], fields[8],
    )


def _feature_values(transcript_id: str, feature_key: str, ordinal: int, fields: list) -> tuple:
    return (
        transcript_id, feature_key, ordinal, fields[0], fields[1], fields[2],
        int(fields[3]), int(fields[4]), fields[5], fields[6], fields[7], fields[8],
    )


def _parse_gtf_to_sqlite(
    annotation_path: os.PathLike | str,
    connection: sqlite3.Connection,
) -> Tuple[int, int, int]:
    five_prime = utr5_type("gtf")
    three_prime = utr3_type("gtf")
    transcript_ids = set()
    feature_batch = []
    gene_count = transcript_count = feature_count = 0
    file_size = _annotation_size(annotation_path)
    with open_auto(annotation_path, "rb") as handle, tqdm(
        total=file_size, unit="B", unit_scale=True, mininterval=0.5,
        ncols=80, desc="Reading annotation", colour="GREEN",
    ) as progress:
        for ordinal, raw_line in enumerate(handle):
            progress.update(len(raw_line))
            if raw_line.startswith(b"#") or not raw_line.strip():
                continue
            fields = raw_line.decode("utf-8").rstrip("\r\n").split("\t", 8)
            if len(fields) != 9:
                raise ValueError("Invalid GTF line %d: expected 9 columns" % (ordinal + 1))
            feature_type = fields[2]
            if "gene" in feature_type:
                gene_id = _quoted_gtf_attribute(fields[8], _GTF_GENE_ID)
                connection.execute(_GENE_UPSERT, _explicit_gene_values(gene_id, ordinal, fields))
                gene_count += 1
                continue
            if feature_type == "miRNA":
                continue
            transcript_id = _quoted_gtf_attribute(fields[8], _GTF_TRANSCRIPT_ID)
            gene_id = _quoted_gtf_attribute(fields[8], _GTF_GENE_ID)
            if "RNA" in feature_type or "transcript" in feature_type:
                connection.execute(
                    _GENE_UPSERT,
                    _placeholder_gene_values(
                        gene_id, ordinal, fields,
                        gene_symbol(fields[8], gene_id),
                    ),
                )
                connection.execute(
                    _TRANSCRIPT_UPSERT,
                    _transcript_values(transcript_id, gene_id, ordinal, fields),
                )
                transcript_ids.add(transcript_id)
                transcript_count += 1
                continue
            if transcript_id not in transcript_ids:
                # Preserve the historical parser's behaviour: subfeatures that
                # precede their transcript record are ignored.
                continue
            feature_key = None
            if feature_type == "CDS":
                feature_key = "CDS"
            elif "exon" in feature_type:
                feature_key = "exon"
            elif feature_type == five_prime:
                feature_key = five_prime
            elif feature_type == three_prime:
                feature_key = three_prime
            elif feature_type in {"start_codon", "stop_codon"}:
                feature_key = feature_type
            if feature_key:
                feature_batch.append(
                    _feature_values(transcript_id, feature_key, ordinal, fields)
                )
                if len(feature_batch) >= 10000:
                    connection.executemany(_FEATURE_INSERT, feature_batch)
                    feature_batch.clear()
                feature_count += 1
        if feature_batch:
            connection.executemany(_FEATURE_INSERT, feature_batch)
    return gene_count, transcript_count, feature_count


def _parse_gff_to_sqlite(
    annotation_path: os.PathLike | str,
    connection: sqlite3.Connection,
) -> Tuple[int, int, int]:
    five_prime = utr5_type("gff")
    three_prime = utr3_type("gff")
    transcript_to_gene = {}
    feature_batch = []
    gene_count = transcript_count = feature_count = 0
    file_size = _annotation_size(annotation_path)
    with open_auto(annotation_path, "rb") as handle, tqdm(
        total=file_size, unit="B", unit_scale=True, mininterval=0.5,
        ncols=80, desc="Reading annotation", colour="GREEN",
    ) as progress:
        for ordinal, raw_line in enumerate(handle):
            progress.update(len(raw_line))
            if raw_line.startswith(b"#") or not raw_line.strip():
                continue
            fields = raw_line.decode("utf-8").rstrip("\r\n").split("\t", 8)
            if len(fields) != 9:
                raise ValueError("Invalid GFF line %d: expected 9 columns" % (ordinal + 1))
            if not fields[8].endswith(";"):
                fields[8] += ";"
            feature_type = fields[2]
            if "gene" in feature_type:
                gene_id = _gff_attribute(fields[8], _GFF_ID_RE, "ID")
                connection.execute(_GENE_UPSERT, _explicit_gene_values(gene_id, ordinal, fields))
                gene_count += 1
            if feature_type == "miRNA":
                continue
            if "RNA" in feature_type or "transcript" in feature_type:
                if "gene" in feature_type:
                    continue
                transcript_id = _gff_attribute(fields[8], _GFF_ID_RE, "ID")
                gene_id = _gff_attribute(fields[8], _GFF_PARENT_RE, "Parent").split(",", 1)[0]
                transcript_to_gene[transcript_id] = gene_id
                connection.execute(
                    _GENE_UPSERT,
                    _placeholder_gene_values(gene_id, ordinal, fields, gene_id),
                )
                connection.execute(
                    _TRANSCRIPT_UPSERT,
                    _transcript_values(transcript_id, gene_id, ordinal, fields),
                )
                transcript_count += 1
                continue
            feature_key = None
            if feature_type == "CDS":
                feature_key = "CDS"
            elif feature_type == "exon":
                feature_key = "exon"
            elif feature_type == five_prime:
                feature_key = five_prime
            elif feature_type == three_prime:
                feature_key = three_prime
            elif feature_type in {"start_codon", "stop_codon"}:
                feature_key = feature_type
            if not feature_key:
                continue
            parent_text = _gff_attribute(fields[8], _GFF_PARENT_RE, "Parent")
            for transcript_id in parent_text.split(","):
                if transcript_id not in transcript_to_gene:
                    continue
                feature_batch.append(
                    _feature_values(transcript_id, feature_key, ordinal, fields)
                )
                if len(feature_batch) >= 10000:
                    connection.executemany(_FEATURE_INSERT, feature_batch)
                    feature_batch.clear()
                feature_count += 1
        if feature_batch:
            connection.executemany(_FEATURE_INSERT, feature_batch)
    return gene_count, transcript_count, feature_count


def build_sqlite_database(
    annotation_path: os.PathLike | str,
    style: str,
    output_path: os.PathLike | str,
) -> Path:
    """Stream a GTF/GFF annotation into an indexed SQLite database."""
    style = style.lower()
    if style not in {"gtf", "gff"}:
        raise ValueError("Annotation style must be gtf or gff")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    if temporary_path.exists():
        temporary_path.unlink()
    connection = _connect_for_build(temporary_path)
    try:
        _create_schema(connection)
        with connection:
            if style == "gtf":
                _parse_gtf_to_sqlite(annotation_path, connection)
            else:
                _parse_gff_to_sqlite(annotation_path, connection)

        # Counts are derived from committed rows rather than parser event
        # counts because duplicate transcript/gene records are upserted.
        gene_count = connection.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
        transcript_count = connection.execute("SELECT COUNT(*) FROM transcripts").fetchone()[0]
        feature_count = connection.execute("SELECT COUNT(*) FROM features").fetchone()[0]
        with connection:
            connection.executemany(
                "INSERT INTO metadata(key, value) VALUES (?, ?)",
                [
                    ("schema_version", SCHEMA_VERSION),
                    ("style", style),
                    ("gene_count", str(gene_count)),
                    ("transcript_count", str(transcript_count)),
                    ("feature_count", str(feature_count)),
                ],
            )
            finalization_steps = (
                "CREATE UNIQUE INDEX genes_ordinal_idx ON genes(ordinal, gene_id)",
                "CREATE INDEX transcripts_gene_idx ON transcripts(gene_id, ordinal, transcript_id)",
                "CREATE UNIQUE INDEX transcripts_ordinal_idx ON transcripts(ordinal, transcript_id)",
                "CREATE INDEX features_transcript_idx ON features(transcript_id, feature_key, start, feature_id)",
                "CREATE INDEX features_key_idx ON features(feature_key, transcript_id)",
                "ANALYZE",
            )
            with tqdm(
                total=len(finalization_steps), unit="step", mininterval=0.2,
                ncols=80, desc="Finalizing database", colour="GREEN",
            ) as progress:
                for statement in finalization_steps:
                    connection.execute(statement)
                    progress.update(1)
        connection.execute("PRAGMA optimize")
    except Exception:
        connection.close()
        if temporary_path.exists():
            temporary_path.unlink()
        raise
    else:
        connection.close()
        os.replace(temporary_path, output_path)
    return output_path


class _SQLiteReader:
    """Process-local read-only SQLite connection owner."""

    def __init__(self, path: os.PathLike | str):
        self.path = os.fspath(Path(path).resolve())
        self._connection = None
        self._pid = None
        self._legacy_layout = None
        self._has_gene_symbol_column = None

    def _connect(self) -> sqlite3.Connection:
        pid = os.getpid()
        if self._connection is None or self._pid != pid:
            if self._connection is not None:
                self._connection.close()
            uri = Path(self.path).as_uri() + "?mode=ro"
            self._connection = sqlite3.connect(uri, uri=True)
            self._connection.execute("PRAGMA query_only=ON")
            self._connection.execute("PRAGMA temp_store=MEMORY")
            self._connection.execute("PRAGMA cache_size=-16384")
            # A bounded mapping keeps whole-genome scans fast without mapping
            # an entire multi-gigabyte human database into process RSS.
            self._connection.execute("PRAGMA mmap_size=67108864")
            self._pid = pid
        return self._connection

    def close(self) -> None:
        if self._connection is not None:
            self._connection.close()
            self._connection = None
            self._pid = None

    @property
    def legacy_layout(self) -> bool:
        """Detect the earlier SQLite layout by columns, not version text.

        Historical databases and the first current-layout release both used
        ``schema_version=1``, so the metadata value alone cannot distinguish
        them.  The old genes table has
        ``chr``/``type`` and no ``ordinal`` column.
        """
        if self._legacy_layout is None:
            columns = {
                row[1] for row in self._connect().execute("PRAGMA table_info(genes)")
            }
            self._legacy_layout = "ordinal" not in columns and "chr" in columns
        return self._legacy_layout

    @property
    def has_gene_symbol_column(self) -> bool:
        """Return whether this database stores symbols in ``genes.symbol``.

        Databases created before schema version 3 do not have this column.
        Their symbols are derived from the original ninth-column annotation,
        which keeps existing SQLite files usable without migration.
        """
        if self._has_gene_symbol_column is None:
            columns = {
                row[1] for row in self._connect().execute("PRAGMA table_info(genes)")
            }
            self._has_gene_symbol_column = "symbol" in columns
        return self._has_gene_symbol_column

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __getstate__(self):
        return {"path": self.path}

    def __setstate__(self, state):
        self.path = state["path"]
        self._connection = None
        self._pid = None
        self._legacy_layout = None
        self._has_gene_symbol_column = None


class SQLiteTranscriptToGene(_SQLiteReader, Mapping):
    """Read-only mapping compatible with the historical ``t2g`` dict."""

    def __len__(self) -> int:
        row = self._connect().execute(
            "SELECT value FROM metadata WHERE key='transcript_count'"
        ).fetchone()
        return int(row[0])

    def __iter__(self) -> Iterator[str]:
        order = "gene_id, transcript_id" if self.legacy_layout else "ordinal, transcript_id"
        cursor = self._connect().execute(
            f"SELECT transcript_id FROM transcripts ORDER BY {order}"
        )
        for row in cursor:
            yield row[0]

    def __getitem__(self, transcript_id: str) -> str:
        row = self._connect().execute(
            "SELECT gene_id FROM transcripts WHERE transcript_id=?",
            (transcript_id,),
        ).fetchone()
        if row is None:
            raise KeyError(transcript_id)
        return row[0]

    def __contains__(self, transcript_id: object) -> bool:
        if not isinstance(transcript_id, str):
            return False
        return self._connect().execute(
            "SELECT 1 FROM transcripts WHERE transcript_id=? LIMIT 1",
            (transcript_id,),
        ).fetchone() is not None


class SQLiteAnnotationDB(_SQLiteReader, Mapping):
    """Lazy, read-only mapping exposing GenoKit's legacy nested structure."""

    def __init__(self, path: os.PathLike | str, cache_genes: int = 8):
        super().__init__(path)
        self.cache_genes = max(0, int(cache_genes))
        self._gene_cache = OrderedDict()
        self._style = None

    @property
    def style(self) -> str:
        if self._style is None:
            row = self._connect().execute(
                "SELECT value FROM metadata WHERE key IN ('style', 'annotation_style') "
                "ORDER BY CASE key WHEN 'style' THEN 0 ELSE 1 END LIMIT 1"
            ).fetchone()
            if row is None:
                raise ValueError(f"Annotation style is missing from SQLite metadata: {self.path}")
            self._style = row[0]
        return self._style

    @property
    def transcript_count(self) -> int:
        row = self._connect().execute(
            "SELECT value FROM metadata WHERE key='transcript_count'"
        ).fetchone()
        return int(row[0])

    def __len__(self) -> int:
        row = self._connect().execute(
            "SELECT value FROM metadata WHERE key='gene_count'"
        ).fetchone()
        return int(row[0])

    def __iter__(self) -> Iterator[str]:
        order = "chr, start, end, gene_id" if self.legacy_layout else "ordinal, gene_id"
        cursor = self._connect().execute(
            f"SELECT gene_id FROM genes ORDER BY {order}"
        )
        for row in cursor:
            yield row[0]

    def __contains__(self, gene_id: object) -> bool:
        if not isinstance(gene_id, str):
            return False
        return self._connect().execute(
            "SELECT 1 FROM genes WHERE gene_id=? LIMIT 1", (gene_id,)
        ).fetchone() is not None

    @staticmethod
    def _empty_transcript(mrna: record, style: str) -> dict:
        return {
            "mRNA": mrna,
            "CDS": [],
            "exon": [],
            utr5_type(style): [],
            utr3_type(style): [],
            "start_codon": [],
            "stop_codon": [],
        }

    @staticmethod
    def _sort_transcript_features(transcript: dict) -> dict:
        """Restore the feature ordering used by the legacy pickle parser.

        GTF records for negative-strand transcripts are commonly written in
        transcript order (high coordinate to low coordinate).  Downstream
        GenoKit code expects every feature list to be in ascending genomic
        coordinate order and applies strand-aware numbering separately.
        """
        for feature_records in transcript.values():
            if isinstance(feature_records, list) and len(feature_records) > 1:
                feature_records.sort(key=lambda item: item.start)
        return transcript

    def _gene_record(self, row) -> record:
        return _record_from_values(row[2:11])

    def gene_symbol(self, gene_id: str) -> str:
        """Return a stored symbol or derive it for an older SQLite database."""
        if self.has_gene_symbol_column:
            row = self._connect().execute(
                "SELECT symbol, attribute FROM genes WHERE gene_id=?", (gene_id,),
            ).fetchone()
            if row is None:
                raise KeyError(gene_id)
            return row[0] or gene_symbol(row[1], gene_id)
        row = self._connect().execute(
            "SELECT attribute FROM genes WHERE gene_id=?", (gene_id,),
        ).fetchone()
        if row is None:
            raise KeyError(gene_id)
        symbol = gene_symbol(row[0], gene_id)
        if symbol == gene_id and self.style == "gtf":
            transcript_row = self._connect().execute(
                "SELECT attribute FROM transcripts WHERE gene_id=? "
                "ORDER BY transcript_id LIMIT 1",
                (gene_id,),
            ).fetchone()
            if transcript_row is not None:
                symbol = gene_symbol(transcript_row[0], gene_id)
        return symbol

    def _cache(self, gene_id: str, gene: dict) -> dict:
        if self.cache_genes:
            self._gene_cache[gene_id] = gene
            self._gene_cache.move_to_end(gene_id)
            while len(self._gene_cache) > self.cache_genes:
                self._gene_cache.popitem(last=False)
        return gene

    def __getitem__(self, gene_id: str) -> dict:
        cached = self._gene_cache.get(gene_id)
        if cached is not None:
            self._gene_cache.move_to_end(gene_id)
            return cached
        connection = self._connect()
        if self.legacy_layout:
            gene_sql = """
                SELECT gene_id, 0 AS ordinal, chr, source, type, start, end,
                       score, strand, phase, attribute, NULL AS symbol
                FROM genes WHERE gene_id=?
            """
        elif self.has_gene_symbol_column:
            gene_sql = """
                SELECT gene_id, ordinal, seqid, source, feature_type, start, end,
                       score, strand, phase, attribute, symbol
                FROM genes WHERE gene_id=?
            """
        else:
            gene_sql = """
                SELECT gene_id, ordinal, seqid, source, feature_type, start, end,
                       score, strand, phase, attribute, NULL AS symbol
                FROM genes WHERE gene_id=?
            """
        gene_row = connection.execute(gene_sql, (gene_id,)).fetchone()
        if gene_row is None:
            raise KeyError(gene_id)
        gene_record = self._gene_record(gene_row)
        gene = {
            "gene": gene_record,
            "symbol": gene_row[11] or gene_symbol(gene_record.attribute, gene_id),
        }
        style = self.style
        if self.legacy_layout:
            transcript_sql = """
                SELECT transcript_id, chr, source, type, start, end,
                       score, strand, phase, attribute
                FROM transcripts WHERE gene_id=?
                ORDER BY start, end, transcript_id
            """
        else:
            transcript_sql = """
                SELECT transcript_id, seqid, source, feature_type, start, end,
                       score, strand, phase, attribute
                FROM transcripts WHERE gene_id=?
                ORDER BY ordinal, transcript_id
            """
        transcript_rows = connection.execute(transcript_sql, (gene_id,)).fetchall()
        for transcript_row in transcript_rows:
            transcript_id = transcript_row[0]
            gene[transcript_id] = self._empty_transcript(
                _record_from_values(transcript_row[1:]), style,
            )
        if self.legacy_layout:
            feature_sql = """
                SELECT f.transcript_id, f.feature_type, f.chr, f.source,
                       f.type, f.start, f.end, f.score, f.strand,
                       f.phase, f.attribute
                FROM features AS f
                WHERE f.gene_id=?
                ORDER BY f.transcript_id, f.feature_type, f.start,
                         f.ordinal, f.feature_id
            """
        else:
            feature_sql = """
                SELECT f.transcript_id, f.feature_key, f.seqid, f.source,
                       f.feature_type, f.start, f.end, f.score, f.strand,
                       f.phase, f.attribute
                FROM features AS f
                JOIN transcripts AS t ON t.transcript_id=f.transcript_id
                WHERE t.gene_id=?
                ORDER BY t.ordinal, f.feature_key, f.start, f.feature_id
            """
        feature_rows = connection.execute(feature_sql, (gene_id,))
        for feature_row in feature_rows:
            transcript_id, feature_key = feature_row[:2]
            if transcript_id in gene:
                gene[transcript_id].setdefault(feature_key, []).append(
                    _record_from_values(feature_row[2:])
                )
        for transcript_id, transcript in gene.items():
            if transcript_id not in {"gene", "symbol"}:
                self._sort_transcript_features(transcript)
        return self._cache(gene_id, gene)

    def items(self):
        """Stream complete gene models using three merge-joined cursors.

        Returning one wide SQL join would repeat the full gene and transcript
        attributes for every exon/CDS/UTR row.  The independent ordered
        cursors below transfer each record exactly once and merge them while
        retaining only the current gene model in Python memory.
        """
        connection = self._connect()
        if self.legacy_layout:
            gene_sql = """
                SELECT gene_id, chr, source, type, start, end, score,
                       strand, phase, attribute, NULL AS symbol
                FROM genes ORDER BY gene_id
            """
            transcript_sql = """
                SELECT t.gene_id, t.transcript_id, t.chr, t.source,
                       t.type, t.start, t.end, t.score, t.strand,
                       t.phase, t.attribute
                FROM transcripts AS t
                ORDER BY t.gene_id, t.transcript_id
            """
            feature_sql = """
                SELECT f.gene_id, f.transcript_id, f.feature_type, f.chr,
                       f.source, f.type, f.start, f.end, f.score,
                       f.strand, f.phase, f.attribute
                FROM features AS f
                ORDER BY f.gene_id, f.transcript_id, f.feature_type,
                         f.start, f.ordinal, f.feature_id
            """
        elif self.has_gene_symbol_column:
            gene_sql = """
                SELECT gene_id, seqid, source, feature_type, start, end, score,
                       strand, phase, attribute, symbol
                FROM genes ORDER BY ordinal, gene_id
            """
            transcript_sql = """
                SELECT t.gene_id, t.transcript_id, t.seqid, t.source,
                       t.feature_type, t.start, t.end, t.score, t.strand,
                       t.phase, t.attribute
                FROM genes AS g
                JOIN transcripts AS t ON t.gene_id=g.gene_id
                ORDER BY g.ordinal, g.gene_id, t.ordinal, t.transcript_id
            """
            feature_sql = """
                SELECT t.gene_id, f.transcript_id, f.feature_key, f.seqid,
                       f.source, f.feature_type, f.start, f.end, f.score,
                       f.strand, f.phase, f.attribute
                FROM genes AS g
                JOIN transcripts AS t ON t.gene_id=g.gene_id
                JOIN features AS f ON f.transcript_id=t.transcript_id
                ORDER BY g.ordinal, g.gene_id, t.ordinal, t.transcript_id,
                         f.feature_key, f.start, f.feature_id
            """
        else:
            gene_sql = """
                SELECT gene_id, seqid, source, feature_type, start, end, score,
                       strand, phase, attribute, NULL AS symbol
                FROM genes ORDER BY ordinal, gene_id
            """
            transcript_sql = """
                SELECT t.gene_id, t.transcript_id, t.seqid, t.source,
                       t.feature_type, t.start, t.end, t.score, t.strand,
                       t.phase, t.attribute
                FROM genes AS g
                JOIN transcripts AS t ON t.gene_id=g.gene_id
                ORDER BY g.ordinal, g.gene_id, t.ordinal, t.transcript_id
            """
            feature_sql = """
                SELECT t.gene_id, f.transcript_id, f.feature_key, f.seqid,
                       f.source, f.feature_type, f.start, f.end, f.score,
                       f.strand, f.phase, f.attribute
                FROM genes AS g
                JOIN transcripts AS t ON t.gene_id=g.gene_id
                JOIN features AS f ON f.transcript_id=t.transcript_id
                ORDER BY g.ordinal, g.gene_id, t.ordinal, t.transcript_id,
                         f.feature_key, f.start, f.feature_id
            """
        gene_rows = connection.execute(gene_sql)
        transcript_rows = iter(connection.execute(transcript_sql))
        feature_rows = iter(connection.execute(feature_sql))
        transcript_row = next(transcript_rows, None)
        feature_row = next(feature_rows, None)
        style = self.style
        for gene_row in gene_rows:
            gene_id = gene_row[0]
            gene_record = _record_from_values(gene_row[1:10])
            gene = {
                "gene": gene_record,
                "symbol": gene_row[10] or gene_symbol(gene_record.attribute, gene_id),
            }
            while transcript_row is not None and transcript_row[0] == gene_id:
                transcript_id = transcript_row[1]
                transcript = self._empty_transcript(
                    _record_from_values(transcript_row[2:]), style,
                )
                while (
                    feature_row is not None
                    and feature_row[0] == gene_id
                    and feature_row[1] == transcript_id
                ):
                    transcript.setdefault(feature_row[2], []).append(
                        _record_from_values(feature_row[3:])
                    )
                    feature_row = next(feature_rows, None)
                gene[transcript_id] = self._sort_transcript_features(transcript)
                transcript_row = next(transcript_rows, None)
            yield gene_id, self._cache(gene_id, gene)

    def values(self):
        for _gene_id, gene in self.items():
            yield gene

    def iter_gene_records(self):
        if self.legacy_layout:
            sql = """
                SELECT gene_id, chr, source, type, start, end, score,
                       strand, phase, attribute
                FROM genes ORDER BY chr, start, end, gene_id
            """
        else:
            sql = """
                SELECT gene_id, seqid, source, feature_type, start, end, score,
                       strand, phase, attribute
                FROM genes ORDER BY ordinal, gene_id
            """
        cursor = self._connect().execute(sql)
        for row in cursor:
            yield row[0], _record_from_values(row[1:])

    def iter_transcript_models(
        self, feature_keys=None, include_feature_attributes=True,
    ):
        """Yield ``(gene_id, transcript_id, transcript_dict)`` lazily.

        ``feature_keys`` lets a command request only the annotation rows it
        actually consumes (for example, exon extraction does not read every
        CDS and UTR row).  All standard keys remain present in the returned
        transcript dictionary, preserving the legacy object contract.
        """
        connection = self._connect()
        if self.legacy_layout:
            transcript_sql = """
                SELECT gene_id, transcript_id, chr, source, type,
                       start, end, score, strand, phase, attribute
                FROM transcripts ORDER BY transcript_id
            """
        else:
            transcript_sql = """
                SELECT gene_id, transcript_id, seqid, source, feature_type,
                       start, end, score, strand, phase, attribute
                FROM transcripts ORDER BY ordinal, transcript_id
            """
        transcript_rows = connection.execute(transcript_sql)
        keys = None if feature_keys is None else tuple(dict.fromkeys(feature_keys))
        if keys == ():
            feature_rows = iter(())
        else:
            where = ""
            parameters = ()
            if keys is not None:
                placeholders = ",".join("?" for _key in keys)
                where = "WHERE f.feature_key IN (%s)" % placeholders
                parameters = keys
            attribute_column = "f.attribute" if include_feature_attributes else "''"
            if self.legacy_layout:
                where = where.replace("f.feature_key", "f.feature_type")
                feature_sql = """
                    SELECT f.transcript_id, f.feature_type, f.chr, f.source,
                           f.type, f.start, f.end, f.score, f.strand,
                           f.phase, %s
                    FROM features AS f
                    %s
                    ORDER BY f.transcript_id, f.feature_type, f.start,
                             f.ordinal, f.feature_id
                """ % (attribute_column, where)
            else:
                feature_sql = """
                    SELECT f.transcript_id, f.feature_key, f.seqid, f.source,
                           f.feature_type, f.start, f.end, f.score, f.strand,
                           f.phase, %s
                    FROM transcripts AS t
                    JOIN features AS f ON f.transcript_id=t.transcript_id
                    %s
                    ORDER BY t.ordinal, t.transcript_id,
                             f.feature_key, f.start, f.feature_id
                """ % (attribute_column, where)
            feature_rows = iter(connection.execute(feature_sql, parameters))
        feature_row = next(feature_rows, None)
        style = self.style
        for transcript_row in transcript_rows:
            gene_id, transcript_id = transcript_row[:2]
            transcript = self._empty_transcript(
                _record_from_values(transcript_row[2:]), style,
            )
            while feature_row is not None and feature_row[0] == transcript_id:
                transcript.setdefault(feature_row[1], []).append(
                    _record_from_values(feature_row[2:])
                )
                feature_row = next(feature_rows, None)
            yield (
                gene_id,
                transcript_id,
                self._sort_transcript_features(transcript),
            )

    def count_transcripts_with_feature(self, feature_key: str) -> int:
        """Return the number of transcripts containing *feature_key*."""
        column = "feature_type" if self.legacy_layout else "feature_key"
        return self._connect().execute(
            f"SELECT COUNT(DISTINCT transcript_id) FROM features WHERE {column}=?",
            (feature_key,),
        ).fetchone()[0]

    def statistics(self, style: Optional[str] = None) -> dict:
        connection = self._connect()
        utr5 = utr5_type(style or self.style)
        utr3 = utr3_type(style or self.style)
        feature_column = "feature_type" if self.legacy_layout else "feature_key"
        gene_count = len(self)
        transcript_count = self.transcript_count
        coding_transcripts = connection.execute(
            f"SELECT COUNT(DISTINCT transcript_id) FROM features WHERE {feature_column}='CDS'"
        ).fetchone()[0]
        coding_genes = connection.execute(
            """
            SELECT COUNT(DISTINCT t.gene_id)
            FROM transcripts AS t
            JOIN features AS f ON f.transcript_id=t.transcript_id
            WHERE f.%s='CDS'
            """ % feature_column
        ).fetchone()[0]
        utr5_count = connection.execute(
            f"SELECT COUNT(DISTINCT transcript_id) FROM features WHERE {feature_column}=?",
            (utr5,),
        ).fetchone()[0]
        utr3_count = connection.execute(
            f"SELECT COUNT(DISTINCT transcript_id) FROM features WHERE {feature_column}=?",
            (utr3,),
        ).fetchone()[0]
        multiple_exon_count = connection.execute(
            """
            SELECT COUNT(*) FROM (
                SELECT transcript_id FROM features
                WHERE %s='exon'
                GROUP BY transcript_id HAVING COUNT(*) >= 2
            )
            """ % feature_column
        ).fetchone()[0]
        return {
            "gene_nums": gene_count,
            "transcript_nums": transcript_count,
            "coding_gene_nums": coding_genes,
            "coding_transcript_nums": coding_transcripts,
            "utr5_transcript_nums": utr5_count,
            "utr3_transcript_nums": utr3_count,
            "multiple_exon_transcript_nums": multiple_exon_count,
        }

    def __getstate__(self):
        return {"path": self.path, "cache_genes": self.cache_genes}

    def __setstate__(self, state):
        self.path = state["path"]
        self.cache_genes = state.get("cache_genes", 8)
        self._connection = None
        self._pid = None
        self._legacy_layout = None
        self._has_gene_symbol_column = None
        self._gene_cache = OrderedDict()
        self._style = None


def load_annotation_database(path: os.PathLike | str):
    """Load either a lazy SQLite database or a legacy pickle database."""
    if is_sqlite_database(path):
        return SQLiteAnnotationDB(path), SQLiteTranscriptToGene(path)
    with open(path, "rb") as handle:
        return pickle.load(handle)


def iter_gene_records(database):
    """Yield ``(gene_id, gene_record)`` without loading transcript features."""
    if hasattr(database, "iter_gene_records"):
        yield from database.iter_gene_records()
    else:
        for gene_id, gene in database.items():
            yield gene_id, gene["gene"]


def iter_gene_records_with_symbols(database):
    """Yield ``(gene_id, symbol, gene_record)`` across both backends."""
    if hasattr(database, "iter_gene_records"):
        for gene_id, gene_record in database.iter_gene_records():
            yield gene_id, database.gene_symbol(gene_id), gene_record
    else:
        for gene_id, gene in database.items():
            yield (
                gene_id,
                gene.get("symbol") or gene_symbol(gene["gene"].attribute, gene_id),
                gene["gene"],
            )


def iter_transcript_models(
    database, feature_keys=None, include_feature_attributes=True,
):
    """Yield transcript models from SQLite or a legacy nested dictionary."""
    if hasattr(database, "iter_transcript_models"):
        yield from database.iter_transcript_models(
            feature_keys, include_feature_attributes,
        )
    else:
        for gene_id, gene in database.items():
            for transcript_id, transcript in gene.items():
                if transcript_id not in {"gene", "symbol"}:
                    yield gene_id, transcript_id, transcript


def get_gene_symbol(database, gene_id, gene=None) -> str:
    """Return a gene symbol across SQLite and legacy nested databases."""
    if hasattr(database, "gene_symbol"):
        return database.gene_symbol(gene_id)
    gene = database[gene_id] if gene is None else gene
    return gene.get("symbol") or gene_symbol(gene["gene"].attribute, gene_id)


def iter_transcript_models_with_symbols(
    database, feature_keys=None, include_feature_attributes=True,
):
    """Yield ``(gene_id, symbol, transcript_id, transcript)`` lazily."""
    if hasattr(database, "iter_transcript_models"):
        current_gene_id = None
        current_symbol = None
        for gene_id, transcript_id, transcript in database.iter_transcript_models(
            feature_keys, include_feature_attributes,
        ):
            if gene_id != current_gene_id:
                current_gene_id = gene_id
                current_symbol = database.gene_symbol(gene_id)
            yield gene_id, current_symbol, transcript_id, transcript
    else:
        for gene_id, gene in database.items():
            symbol = gene.get("symbol") or gene_symbol(
                gene["gene"].attribute, gene_id,
            )
            for transcript_id, transcript in gene.items():
                if transcript_id not in {"gene", "symbol"}:
                    yield gene_id, symbol, transcript_id, transcript


def transcript_count(database, t2g=None) -> int:
    """Return transcript count without walking an SQLite gene hierarchy."""
    if hasattr(database, "transcript_count"):
        return database.transcript_count
    if t2g is not None:
        return len(t2g)
    return sum(
        len(gene) - ("gene" in gene) - ("symbol" in gene)
        for gene in database.values()
    )
