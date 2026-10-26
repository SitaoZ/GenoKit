"""Unified k-mer databases used by GenoKit design commands.

The SQLite file is a lightweight catalogue.  Genome-scale counts stay in the
compact Jellyfish format and the catalogue records which file belongs to each
``(source_type, k)`` dataset.  Small Python-engine databases may be stored
directly in SQLite.  Keeping the large count payload in Jellyfish avoids
turning a several-GB binary hash into a much larger SQLite table while still
providing one stable database path to all design commands.
"""

from __future__ import annotations

import json
import os
import pickle
import shutil
import sqlite3
import subprocess
from datetime import datetime, timezone
from typing import Dict, Iterable, Mapping, Optional


SCHEMA_VERSION = 1
SOURCE_TYPES = {"genome", "transcriptome"}
_SQLITE_HEADER = b"SQLite format 3\x00"


def _normalise_sequence(sequence: str) -> str:
    return str(sequence).upper().replace("U", "T")


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _canonical_sequence(sequence: str) -> str:
    reverse = sequence.translate(_COMPLEMENT)[::-1]
    return min(sequence, reverse)


def _relative_payload_path(database_path: str, payload_path: str) -> str:
    database_dir = os.path.dirname(os.path.abspath(database_path))
    payload_path = os.path.abspath(payload_path)
    try:
        return os.path.relpath(payload_path, database_dir)
    except ValueError:
        return payload_path


def _resolve_payload_path(database_path: str, payload_path: str) -> str:
    if os.path.isabs(payload_path):
        return payload_path
    return os.path.abspath(os.path.join(os.path.dirname(database_path), payload_path))


def initialise_kmer_catalog(path: str) -> None:
    """Create or upgrade a GenoKit k-mer SQLite catalogue."""
    path = os.path.abspath(os.path.expanduser(path))
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with sqlite3.connect(path) as connection:
        connection.executescript(
            """
            PRAGMA journal_mode=WAL;
            CREATE TABLE IF NOT EXISTS metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            ) WITHOUT ROWID;
            CREATE TABLE IF NOT EXISTS datasets (
                id INTEGER PRIMARY KEY,
                source_type TEXT NOT NULL CHECK(source_type IN ('genome', 'transcriptome')),
                k INTEGER NOT NULL CHECK(k > 0),
                storage TEXT NOT NULL CHECK(storage IN ('jellyfish', 'sqlite')),
                payload_path TEXT,
                jellyfish_path TEXT,
                canonical INTEGER NOT NULL DEFAULT 1,
                source TEXT,
                unique_kmers INTEGER,
                created_at TEXT NOT NULL,
                UNIQUE(source_type, k)
            );
            CREATE TABLE IF NOT EXISTS counts (
                dataset_id INTEGER NOT NULL,
                kmer TEXT NOT NULL,
                count INTEGER NOT NULL,
                PRIMARY KEY(dataset_id, kmer),
                FOREIGN KEY(dataset_id) REFERENCES datasets(id) ON DELETE CASCADE
            ) WITHOUT ROWID;
            """
        )
        connection.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES('schema_version', ?)",
            (str(SCHEMA_VERSION),),
        )
        dataset_columns = {
            row[1] for row in connection.execute("PRAGMA table_info(datasets)")
        }
        if "jellyfish_path" not in dataset_columns:
            connection.execute("ALTER TABLE datasets ADD COLUMN jellyfish_path TEXT")


def register_jellyfish_dataset(
    catalog_path: str,
    payload_path: str,
    k: int,
    source_type: str,
    canonical: bool,
    source: Optional[str] = None,
    jellyfish_path: Optional[str] = None,
) -> None:
    """Register a Jellyfish payload, replacing the same source/k entry."""
    if source_type not in SOURCE_TYPES:
        raise ValueError(f"Unknown k-mer source type: {source_type!r}")
    initialise_kmer_catalog(catalog_path)
    stored_path = _relative_payload_path(catalog_path, payload_path)
    with sqlite3.connect(catalog_path) as connection:
        old = connection.execute(
            "SELECT id FROM datasets WHERE source_type=? AND k=?",
            (source_type, int(k)),
        ).fetchone()
        if old:
            connection.execute("DELETE FROM counts WHERE dataset_id=?", (old[0],))
        connection.execute(
            """
            INSERT INTO datasets(
                source_type, k, storage, payload_path, jellyfish_path, canonical, source,
                unique_kmers, created_at
            ) VALUES(?, ?, 'jellyfish', ?, ?, ?, ?, NULL, ?)
            ON CONFLICT(source_type, k) DO UPDATE SET
                storage=excluded.storage,
                payload_path=excluded.payload_path,
                jellyfish_path=excluded.jellyfish_path,
                canonical=excluded.canonical,
                source=excluded.source,
                unique_kmers=NULL,
                created_at=excluded.created_at
            """,
            (
                source_type,
                int(k),
                stored_path,
                os.path.abspath(os.path.expanduser(jellyfish_path)) if jellyfish_path else None,
                int(bool(canonical)),
                source,
                datetime.now(timezone.utc).isoformat(),
            ),
        )


def register_sqlite_counts(
    catalog_path: str,
    counts: Mapping[str, int],
    k: int,
    source_type: str,
    canonical: bool,
    source: Optional[str] = None,
) -> None:
    """Store a small count dictionary directly in the catalogue."""
    if source_type not in SOURCE_TYPES:
        raise ValueError(f"Unknown k-mer source type: {source_type!r}")
    initialise_kmer_catalog(catalog_path)
    with sqlite3.connect(catalog_path) as connection:
        connection.execute("PRAGMA foreign_keys=ON")
        old = connection.execute(
            "SELECT id FROM datasets WHERE source_type=? AND k=?",
            (source_type, int(k)),
        ).fetchone()
        if old:
            connection.execute("DELETE FROM counts WHERE dataset_id=?", (old[0],))
            connection.execute("DELETE FROM datasets WHERE id=?", (old[0],))
        cursor = connection.execute(
            """
            INSERT INTO datasets(
                source_type, k, storage, payload_path, canonical, source,
                unique_kmers, created_at
            ) VALUES(?, ?, 'sqlite', NULL, ?, ?, ?, ?)
            """,
            (
                source_type,
                int(k),
                int(bool(canonical)),
                source,
                len(counts),
                datetime.now(timezone.utc).isoformat(),
            ),
        )
        dataset_id = cursor.lastrowid
        rows = (
            (dataset_id, _normalise_sequence(kmer), int(count))
            for kmer, count in counts.items()
        )
        connection.executemany(
            "INSERT INTO counts(dataset_id, kmer, count) VALUES(?, ?, ?)", rows
        )


def is_sqlite_database(path: str) -> bool:
    try:
        with open(path, "rb") as handle:
            return handle.read(len(_SQLITE_HEADER)) == _SQLITE_HEADER
    except OSError:
        return False


def _resolve_executable(value: Optional[str], default: str) -> Optional[str]:
    if value:
        expanded = os.path.expanduser(value)
        if os.path.sep in expanded:
            return expanded if os.path.isfile(expanded) and os.access(expanded, os.X_OK) else None
        return shutil.which(expanded)
    return shutil.which(default)


class JellyfishKmerCounts:
    """Mapping-like adapter around a persistent ``jellyfish query`` process."""

    def __init__(self, path: str, k: Optional[int] = None, jellyfish: Optional[str] = None):
        self.path = os.path.abspath(os.path.expanduser(path))
        self.jellyfish = _resolve_executable(jellyfish, "jellyfish")
        if not self.jellyfish:
            raise ValueError(
                "Jellyfish is required to query this k-mer database; install it or "
                "provide --jellyfish PATH."
            )
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"K-mer payload does not exist: {self.path}")
        self.k = int(k) if k else self._read_k()
        self._process = None
        self._cache: Dict[str, int] = {}

    def _read_k(self) -> int:
        result = subprocess.run(
            [self.jellyfish, "info", "-j", self.path],
            text=True,
            capture_output=True,
            check=True,
        )
        try:
            header = json.loads(result.stdout)
        except json.JSONDecodeError as error:
            raise ValueError(f"Cannot read k from Jellyfish database {self.path}") from error
        for key in ("key_len", "mer_len", "length"):
            if key in header:
                value = int(header[key])
                return value // 2 if key == "key_len" else value
        command = str(header.get("cmdline", header.get("command", "")))
        fields = command.split()
        for flag in ("-m", "--mer-len"):
            if flag in fields and fields.index(flag) + 1 < len(fields):
                return int(fields[fields.index(flag) + 1])
        raise ValueError(f"Cannot determine k for Jellyfish database {self.path}; pass --kmer-size.")

    def _start(self) -> None:
        if self._process is not None and self._process.poll() is None:
            return
        self._process = subprocess.Popen(
            [self.jellyfish, "query", "-i", self.path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

    def get(self, sequence: str, default: int = 0) -> int:
        sequence = _normalise_sequence(sequence)
        if len(sequence) != self.k or set(sequence) - set("ACGT"):
            return default
        if sequence in self._cache:
            return self._cache[sequence]
        self._start()
        assert self._process is not None
        assert self._process.stdin is not None
        assert self._process.stdout is not None
        self._process.stdin.write(sequence + "\n")
        self._process.stdin.flush()
        response = self._process.stdout.readline().strip()
        if not response:
            message = ""
            if self._process.stderr is not None:
                message = self._process.stderr.read().strip()
            raise RuntimeError(f"Jellyfish query failed for {self.path}: {message}")
        count = int(response.split()[-1])
        if len(self._cache) >= 100_000:
            self._cache.clear()
        self._cache[sequence] = count
        return count

    def supports_length(self, length: int) -> bool:
        return int(length) == self.k

    def get_exact(self, sequence: str, default=None):
        sequence = _normalise_sequence(sequence)
        if not self.supports_length(len(sequence)):
            return default
        return self.get(sequence, 0)

    def get_many(self, sequences: Iterable[str]) -> Dict[str, int]:
        unique = list(dict.fromkeys(_normalise_sequence(seq) for seq in sequences))
        result = {
            sequence: self._cache.get(sequence, 0)
            for sequence in unique
            if sequence in self._cache
        }
        missing = [
            sequence for sequence in unique
            if sequence not in result and len(sequence) == self.k and not (set(sequence) - set("ACGT"))
        ]
        if missing:
            self._start()
            assert self._process is not None
            assert self._process.stdin is not None
            assert self._process.stdout is not None
            for start in range(0, len(missing), 1000):
                chunk = missing[start:start + 1000]
                self._process.stdin.write("\n".join(chunk) + "\n")
                self._process.stdin.flush()
                for sequence in chunk:
                    response = self._process.stdout.readline().strip()
                    if not response:
                        raise RuntimeError(f"Jellyfish batch query failed for {self.path}")
                    result[sequence] = int(response.split()[-1])
        for sequence in unique:
            result.setdefault(sequence, 0)
        if len(self._cache) + len(result) >= 100_000:
            self._cache.clear()
        self._cache.update(result)
        return result

    def __getitem__(self, sequence: str) -> int:
        return self.get(sequence, 0)

    def close(self) -> None:
        if self._process is None:
            return
        if self._process.stdin is not None:
            self._process.stdin.close()
        try:
            self._process.wait(timeout=2)
        except subprocess.TimeoutExpired:
            self._process.terminate()
        self._process = None

    def __bool__(self) -> bool:
        return True

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


class PickleKmerCounts:
    """Compatibility adapter for legacy GenoKit pickle dictionaries."""

    def __init__(self, path: str, k: Optional[int] = None):
        with open(path, "rb") as handle:
            counts = pickle.load(handle)
        if not isinstance(counts, dict):
            raise ValueError(f"Invalid legacy k-mer dictionary: {path}")
        self._counts = counts
        self.k = int(k) if k else (len(next(iter(counts))) if counts else 0)

    def get(self, sequence: str, default: int = 0) -> int:
        return int(self._counts.get(_normalise_sequence(sequence), default))

    def supports_length(self, length: int) -> bool:
        return int(length) == self.k

    def get_exact(self, sequence: str, default=None):
        sequence = _normalise_sequence(sequence)
        if not self.supports_length(len(sequence)):
            return default
        return self.get(sequence, 0)

    def get_many(self, sequences: Iterable[str]) -> Dict[str, int]:
        return {sequence: self.get(sequence, 0) for sequence in sequences}

    def __getitem__(self, sequence: str) -> int:
        return self.get(sequence, 0)

    def __bool__(self) -> bool:
        return bool(self._counts)


class KmerCatalog:
    """Read-only, source-aware view of one GenoKit k-mer catalogue."""

    def __init__(
        self,
        path: str,
        k: Optional[int] = None,
        expected_source: Optional[str] = None,
        jellyfish: Optional[str] = None,
    ):
        self.path = os.path.abspath(os.path.expanduser(path))
        if expected_source and expected_source not in SOURCE_TYPES:
            raise ValueError(f"Unknown expected k-mer source: {expected_source!r}")
        self.expected_source = expected_source
        self.jellyfish = jellyfish
        self._connection = sqlite3.connect(f"file:{self.path}?mode=ro", uri=True)
        self._connection.row_factory = sqlite3.Row
        self._adapters: Dict[int, object] = {}
        rows = self._connection.execute(
            "SELECT * FROM datasets ORDER BY source_type, k"
        ).fetchall()
        if not rows:
            raise ValueError(f"No k-mer datasets were found in {self.path}")
        if expected_source:
            rows = [row for row in rows if row["source_type"] == expected_source]
            if not rows:
                available = sorted(
                    {row[0] for row in self._connection.execute("SELECT DISTINCT source_type FROM datasets")}
                )
                raise ValueError(
                    f"This command requires a {expected_source} k-mer dataset, but "
                    f"{self.path} contains: {', '.join(available)}."
                )
        else:
            selected = [row for row in rows if k is None or int(row["k"]) == int(k)]
            sources_by_k: Dict[int, set[str]] = {}
            for row in selected:
                sources_by_k.setdefault(int(row["k"]), set()).add(row["source_type"])
            ambiguous = sorted(value for value, sources in sources_by_k.items() if len(sources) > 1)
            if ambiguous:
                values = ", ".join(str(value) for value in ambiguous)
                raise ValueError(
                    f"The k-mer source is ambiguous for k={values} in {self.path}. "
                    "Select genome or transcriptome explicitly (primer: --kmer-source)."
                )
        self._rows = {int(row["k"]): row for row in rows}
        if k is not None and int(k) not in self._rows:
            available = ", ".join(str(value) for value in sorted(self._rows))
            raise ValueError(
                f"k={k} is not available for {expected_source or 'the selected source'} "
                f"in {self.path}; available k values: {available}."
            )
        selected_rows = [self._rows[int(k)]] if k is not None else list(self._rows.values())
        if expected_source and any(not int(row["canonical"]) for row in selected_rows):
            raise ValueError(
                f"The {expected_source} k-mer dataset in {self.path} was built without "
                "reverse-complement canonical counting. Rebuild it with GenoKit kmer -r."
            )
        self.k = int(k) if k is not None else None

    @property
    def available_k(self):
        return tuple(sorted(self._rows))

    def supports_length(self, length: int) -> bool:
        return int(length) in self._rows

    def get_exact(self, sequence: str, default=None):
        """Query the complete sequence only; never substitute a shorter k."""
        sequence = _normalise_sequence(sequence)
        query_k = len(sequence)
        if not self.supports_length(query_k):
            return default
        return int(self._adapter(query_k).get(sequence, 0))

    def _adapter(self, k: int):
        if k in self._adapters:
            return self._adapters[k]
        row = self._rows[k]
        if row["storage"] == "jellyfish":
            payload = _resolve_payload_path(self.path, row["payload_path"])
            stored_executable = (
                row["jellyfish_path"] if "jellyfish_path" in row.keys() else None
            )
            executable = self.jellyfish or stored_executable
            adapter = JellyfishKmerCounts(payload, k=k, jellyfish=executable)
        else:
            adapter = _SQLiteDatasetCounts(
                self._connection, int(row["id"]), k, bool(row["canonical"])
            )
        self._adapters[k] = adapter
        return adapter

    def _query_k(self, sequence_length: int) -> Optional[int]:
        if sequence_length in self._rows:
            return sequence_length
        if self.k is not None and self.k <= sequence_length:
            return self.k
        choices = [value for value in self._rows if value <= sequence_length]
        return max(choices) if choices else None

    def get(self, sequence: str, default: int = 0) -> int:
        sequence = _normalise_sequence(sequence)
        query_k = self._query_k(len(sequence))
        if query_k is None:
            return default
        adapter = self._adapter(query_k)
        if len(sequence) == query_k:
            return int(adapter.get(sequence, default))
        return max(
            int(adapter.get(sequence[start : start + query_k], default))
            for start in range(len(sequence) - query_k + 1)
        )

    def get_many(self, sequences: Iterable[str]) -> Dict[str, int]:
        unique = list(dict.fromkeys(_normalise_sequence(seq) for seq in sequences))
        grouped: Dict[int, list[str]] = {}
        result: Dict[str, int] = {}
        for sequence in unique:
            query_k = self._query_k(len(sequence))
            if query_k is None:
                result[sequence] = 0
            elif len(sequence) == query_k:
                grouped.setdefault(query_k, []).append(sequence)
            else:
                result[sequence] = self.get(sequence, 0)
        for query_k, group in grouped.items():
            result.update(self._adapter(query_k).get_many(group))
        return result

    def __getitem__(self, sequence: str) -> int:
        return self.get(sequence, 0)

    def __bool__(self) -> bool:
        return True

    def close(self) -> None:
        for adapter in self._adapters.values():
            close = getattr(adapter, "close", None)
            if close:
                close()
        self._connection.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


class _SQLiteDatasetCounts:
    def __init__(
        self, connection: sqlite3.Connection, dataset_id: int, k: int, canonical: bool
    ):
        self.connection = connection
        self.dataset_id = dataset_id
        self.k = k
        self.canonical = canonical

    def _key(self, sequence: str) -> str:
        sequence = _normalise_sequence(sequence)
        return _canonical_sequence(sequence) if self.canonical else sequence

    def get(self, sequence: str, default: int = 0) -> int:
        sequence = self._key(sequence)
        row = self.connection.execute(
            "SELECT count FROM counts WHERE dataset_id=? AND kmer=?",
            (self.dataset_id, sequence),
        ).fetchone()
        return int(row[0]) if row else default

    def get_many(self, sequences: Iterable[str]) -> Dict[str, int]:
        requested = list(dict.fromkeys(_normalise_sequence(seq) for seq in sequences))
        key_to_requested = {}
        for sequence in requested:
            key_to_requested.setdefault(self._key(sequence), []).append(sequence)
        keys = list(key_to_requested)
        count_by_key = {key: 0 for key in keys}
        for start in range(0, len(keys), 900):
            chunk = keys[start : start + 900]
            placeholders = ",".join("?" for _ in chunk)
            rows = self.connection.execute(
                f"SELECT kmer, count FROM counts WHERE dataset_id=? AND kmer IN ({placeholders})",
                (self.dataset_id, *chunk),
            )
            count_by_key.update((row[0], int(row[1])) for row in rows)
        return {
            sequence: count_by_key[key]
            for key, sequences_for_key in key_to_requested.items()
            for sequence in sequences_for_key
        }


def open_kmer_database(
    path: Optional[str],
    k: Optional[int] = None,
    expected_source: Optional[str] = None,
    jellyfish: Optional[str] = None,
):
    """Open SQLite, Jellyfish, or legacy pickle counts through one API."""
    if not path:
        return {}
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.exists(path):
        raise FileNotFoundError(f"K-mer database does not exist: {path}")
    if is_sqlite_database(path):
        return KmerCatalog(path, k=k, expected_source=expected_source, jellyfish=jellyfish)

    # A raw Jellyfish file is accepted for backward compatibility, but it has
    # no source metadata; only a catalogue can enforce genome/transcriptome use.
    executable = _resolve_executable(jellyfish, "jellyfish")
    if executable:
        for info_args in (("info", "-j", path), ("info", path)):
            probe = subprocess.run(
                [executable, *info_args],
                text=True,
                capture_output=True,
            )
            if probe.returncode == 0:
                return JellyfishKmerCounts(path, k=k, jellyfish=executable)
    try:
        return PickleKmerCounts(path, k=k)
    except (pickle.UnpicklingError, EOFError, ValueError) as error:
        if not executable:
            raise ValueError(
                f"{path} is not SQLite or a legacy pickle database. If it is a "
                "Jellyfish file, install Jellyfish or provide --jellyfish PATH."
            ) from error
        raise ValueError(f"Unsupported k-mer database format: {path}") from error


def kmer_frequencies(database, sequences: Iterable[str]) -> Dict[str, int]:
    """Return frequencies in one batch when the backend supports it."""
    sequences = list(dict.fromkeys(_normalise_sequence(seq) for seq in sequences))
    if not database:
        return {sequence: 0 for sequence in sequences}
    get_many = getattr(database, "get_many", None)
    if get_many:
        return {sequence: int(count) for sequence, count in get_many(sequences).items()}
    return {sequence: int(database.get(sequence, 0)) for sequence in sequences}
