"""Optional k-mer and Bowtie2 specificity checks for primer candidates."""

import os
import shutil
import subprocess
import tempfile

from GenoKit.database.kmer_database import open_kmer_database


class PrimerSpecificityChecker:
    def __init__(self, kmer_database=None, jellyfish=None, bowtie2=None,
                 bowtie2_index=None, seed_size=12, max_offtargets=1,
                 kmer_size=None, expected_source=None):
        self.kmer_database = kmer_database
        self.jellyfish = self._executable(jellyfish, 'jellyfish')
        self.bowtie2 = self._executable(bowtie2, 'bowtie2')
        self.bowtie2_index = bowtie2_index
        self.seed_size = seed_size
        self.max_offtargets = max_offtargets
        self._counts = (
            open_kmer_database(
                kmer_database,
                k=kmer_size,
                expected_source=expected_source,
                jellyfish=jellyfish,
            )
            if isinstance(kmer_database, (str, os.PathLike))
            else (kmer_database or {})
        )

    @staticmethod
    def _executable(value, default):
        if value:
            return value if os.path.isfile(value) else shutil.which(value)
        return shutil.which(default)

    def count(self, sequence):
        sequence = sequence.upper()
        if not self._counts:
            return None
        get_exact = getattr(self._counts, 'get_exact', None)
        value = get_exact(sequence, None) if get_exact else self._counts.get(sequence, 0)
        return None if value is None else int(value)

    def _seed(self, primer):
        """Use the requested 3' seed or the nearest available catalogue k."""
        length = self.seed_size
        available = getattr(self._counts, 'available_k', ())
        if available and length not in available:
            eligible = [value for value in available if value <= len(primer)]
            if eligible:
                shorter = [value for value in eligible if value <= length]
                length = max(shorter) if shorter else min(eligible)
        return primer[-length:]

    def passes_kmer_filter(self, primer):
        if not self.kmer_database:
            return True
        if len(primer) < self.seed_size:
            return False
        full = self.count(primer)
        seed = self.count(self._seed(primer))
        return full == 1 and seed == 1

    def alignment_count(self, primer):
        if not self.bowtie2_index:
            return None
        if not self.bowtie2:
            raise ValueError('Bowtie2 index supplied but Bowtie2 is unavailable.')
        with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as query:
            query.write(f'>primer\n{primer}\n')
            query_path = query.name
        try:
            result = subprocess.run([self.bowtie2, '-x', self.bowtie2_index, '-f', '-U', query_path,
                                     '-a', '--end-to-end', '--very-sensitive', '-S', os.devnull],
                                    text=True, capture_output=True, check=True)
            return sum(1 for line in result.stderr.splitlines() if 'aligned exactly 1 time' in line)
        finally:
            os.unlink(query_path)

    def passes_primer(self, primer):
        if not self.passes_kmer_filter(primer):
            return False
        alignments = self.alignment_count(primer)
        return alignments is None or alignments <= self.max_offtargets
