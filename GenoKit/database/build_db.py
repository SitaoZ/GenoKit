#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: build_db.py
date: 2025/11/29 15:28
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import os
import sys
import gc
import time
import pickle
import logging
from tqdm import tqdm
from GenoKit.utils.util import gff_feature_dict
from GenoKit.utils.util import gtf_feature_dict
from GenoKit.database.annotation_db import build_sqlite_database, detect_annotation_style
from pyfaidx import Faidx, UnsupportedCompressionFormat

class _ProgressFileWriter:
    """Expose file.write while reporting the bytes emitted by pickle."""

    def __init__(self, file_handle, progress):
        self._file_handle = file_handle
        self._progress = progress
        self._pending_bytes = 0

    def write(self, data):
        written = self._file_handle.write(data)
        self._pending_bytes += written
        if self._pending_bytes >= 4 * 1024 * 1024:
            self.flush_progress()
        return written

    def flush_progress(self):
        if self._pending_bytes:
            self._progress.update(self._pending_bytes)
            self._pending_bytes = 0

    def __getattr__(self, name):
        return getattr(self._file_handle, name)


def _write_database(database_id, db, t2g):
    """Serialize the database while showing write activity and volume."""
    with open(database_id, 'wb', buffering=1024 * 1024) as file_handle, tqdm(
        unit='B',
        unit_scale=True,
        mininterval=0.5,
        ncols=80,
        desc='Writing database',
        colour='GREEN',
    ) as progress:
        writer = _ProgressFileWriter(file_handle, progress)
        try:
            pickle.dump((db, t2g), writer, protocol=pickle.HIGHEST_PROTOCOL)
        finally:
            writer.flush_progress()


def _ensure_fasta_index(genome_path):
    """检查或建立普通FASTA/BGZF FASTA所需的随机访问索引。

    普通FASTA只需要``.fai``；BGZF压缩的``.fa.gz``同时需要
    ``.fai``和``.gzi``。具体格式识别、索引缺失和索引过期检查交给
    pyfaidx处理，避免只检查``.fai``而漏掉BGZF的``.gzi``。
    """
    fasta_index = None
    try:
        fasta_index = Faidx(genome_path)
    except UnsupportedCompressionFormat as error:
        raise SystemExit(
            "GenoKit create: error: %s is compressed with ordinary gzip. "
            "GenoKit supports an uncompressed FASTA file or a BGZF-compressed "
            "FASTA (.fa.gz). Please decompress it and recompress it with bgzip "
            "before running GenoKit create again. pyfaidx: %s"
            % (genome_path, error)
        ) from None
    finally:
        if fasta_index is not None:
            fasta_index.close()



def create_db(args):
    '''
    This is a time- and memory-intense procedure,
     but it needs to be done only once for a given genome.
    parameters:
     args: arguments from argparse
     args.annotation gff/gtf file
     args.genome genome fasta file
     args.output_dir
     args.style gtf/gff
     args.prefix database prefix
    return:
     a dict object
    '''
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="Create.log"
    )
    logging.info("Parser database format")
    file_style = args.style.lower()
    if file_style not in {'gff', 'gtf'}:
        sys.exit('The file format should be gff or gtf.')
    detected_style = detect_annotation_style(args.annotation)
    if detected_style and detected_style != file_style:
        supplied_name = "GFF3" if file_style == "gff" else "GTF"
        detected_name = "GFF3" if detected_style == "gff" else "GTF"
        raise SystemExit(
            "GenoKit create: error: %s appears to be %s, but -s/--style %s "
            "tells GenoKit to parse it as %s. Use -s %s."
            % (
                args.annotation,
                detected_name,
                file_style,
                supplied_name,
                detected_style,
            )
        )

    database_format = getattr(args, 'database_format', 'sqlite').lower()
    os.makedirs(args.output_dir, exist_ok=True)
    if database_format == 'sqlite':
        # The SQLite writer parses and inserts records incrementally.  It never
        # constructs the genome-wide nested Python object, which is the main
        # source of peak memory in the historical pickle build.
        database_id = os.path.join(
            args.output_dir, args.prefix + '.' + file_style + '.sqlite'
        )
        build_started = time.perf_counter()
        build_sqlite_database(args.annotation, file_style, database_id)
        logging.info(
            'SQLite database built in %.1fs: %s',
            time.perf_counter() - build_started,
            database_id,
        )
    elif database_format == 'pickle':
        # Keep the old writer available for reproducibility and migration
        # checks.  Downstream commands auto-detect both formats.
        gc_was_enabled = gc.isenabled()
        if gc_was_enabled:
            gc.disable()
        try:
            parse_started = time.perf_counter()
            if file_style == 'gff':
                logging.info("Reading gff file")
                db, t2g = gff_feature_dict(args.annotation, file_style)
            else:
                logging.info("Reading gtf file")
                db, t2g = gtf_feature_dict(args.annotation, file_style)
            logging.info(
                'Annotation parsed in %.1fs.',
                time.perf_counter() - parse_started,
            )
            database_id = os.path.join(
                args.output_dir, args.prefix + '.' + file_style + '.pkl'
            )
            write_started = time.perf_counter()
            _write_database(database_id, db, t2g)
            logging.info(
                'Database written in %.1fs.',
                time.perf_counter() - write_started,
            )
        finally:
            if gc_was_enabled:
                gc.enable()
    else:
        raise ValueError('Database format must be sqlite or pickle.')
    print('Database written: %s' % database_id)
    # Faidx自动识别普通FASTA和BGZF FASTA。即使.fai已经存在也需要调用，
    # 因为BGZF还需要.gzi，并且现有索引可能已经落后于FASTA文件。
    if not getattr(args, 'skip_fasta_index', False):
        logging.info('Checking/building FASTA index...')
        _ensure_fasta_index(args.genome)
