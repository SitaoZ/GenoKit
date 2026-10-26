# -*- coding: utf-8 -*-
'''
filename: extract_view.py
date: 2025/11/28 14:46
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import sys
import pycircos
import numpy as np
from pathlib import Path
# from pycirclize import Circos
from GenoKit.database.database import genome_dict
from GenoKit.visualize.circos_genome import chrom_sector,cytoband
from GenoKit.visualize.circos_genome import scatter,line,bar,heatmap,link


np.random.seed(0)


_PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def _resolve_input_path(value):
    """Resolve a user path, falling back to GenoKit's packaged examples."""
    if not value:
        return None
    requested = Path(value).expanduser()
    if requested.is_file():
        return str(requested)
    if not requested.is_absolute() and requested.parts and requested.parts[0] == "example":
        packaged = _PACKAGE_ROOT.joinpath(*requested.parts)
        if packaged.is_file():
            return str(packaged)
    raise FileNotFoundError(
        "Circos input file not found: %s. Relative example/... paths are also "
        "searched under %s." % (value, _PACKAGE_ROOT / "example")
    )


def _track_type(path):
    """Distinguish a four-column line track from a six-column link track."""
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().strip()
    fields = header.split("\t") if "\t" in header and "," not in header else header.split(",")
    normalized = [field.strip().lower() for field in fields]
    if len(normalized) >= 6 or {"chr1", "start1", "end1", "chr2", "start2", "end2"}.issubset(normalized):
        return "link"
    if len(normalized) >= 4:
        return "line"
    raise ValueError(
        "Cannot detect the track type for %s: expected a four-column line "
        "track or a six-column link track" % path
    )


def _allocate_track_ranges(track_paths, outer_inner_radius=950, gap=10, width=80):
    """Place only the selected radial tracks in consecutive outer-to-inner slots."""
    ranges = {}
    upper = outer_inner_radius - gap
    for name in ("point", "line", "bar", "heatmap"):
        if not track_paths.get(name):
            continue
        lower = upper - width
        if lower <= 0:
            raise ValueError("Too many Circos tracks for the available radial space")
        ranges[name] = (lower, upper)
        upper = lower - gap
    return ranges


def chromosome_sort_key(chrom):
    """
    将染色体名称转换为排序键
    格式: chr[数字或X/Y]
    """
    # 去掉'chr'前缀
    if chrom.startswith('chr'):
        num = chrom[3:]
    else:
        num = chrom

    # 处理特殊情况
    if num == 'X':
        return (1, 23)  # X染色体在23号之后
    elif num == 'Y':
        return (1, 24)  # Y染色体在X之后
    elif num == 'M' or num == 'MT':
        return (1, 25)  # 线粒体染色体最后
    else:
        # 数字染色体按数值排序
        try:
            return (0, int(num))
        except ValueError:
            # 处理其他情况（如chr1_KI270706v1_random）
            return (2, chrom)  # 其他情况按字符串排序

# method 1 (not used)
def get_circos_pycirclize(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    #if args.style == 'GFF':
    #    db, t2g = gff_feature_dict(args.database, args.style)
    #else:
    #    db, t2g = gtf_feature_dict(args.database, args.style)
    # with open(args.database, 'rb') as f:
    #     db, t2g = cPickle.load(f)
    genome = genome_dict(args.genome)
    chroms = [ c for c in genome.keys() if c.startswith("chr") and ("_" not in c) and ("chrM" not in c) ]

    chroms_sort = sorted(chroms, key=chromosome_sort_key)
    print(chroms_sort)
    print(len(chroms))
    # Initialize Circos sectors
    sectors = {}
    for chrom in chroms_sort:
        sectors[chrom] = len(genome[chrom][1:20])
        # sectors = {"A": 10, "B": 15, "C": 12, "D": 20, "E": 15}

    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        # Plot sector name
        sector.text(f"{sector.name}", r=110, size=15)
        # Create x positions & random y values
        x = np.arange(sector.start, sector.end) + 0.5
        y = np.random.randint(0, 100, len(x))
        # Plot lines
        track1 = sector.add_track((80, 100), r_pad_ratio=0.1)
        track1.xticks_by_interval(interval=1)
        track1.axis()
        track1.line(x, y)
        # Plot points
        track2 = sector.add_track((55, 75), r_pad_ratio=0.1)
        track2.axis()
        track2.scatter(x, y)
        # Plot bars
        track3 = sector.add_track((30, 50), r_pad_ratio=0.1)
        track3.axis()
        track3.bar(x, y)

    # Plot links
    circos.link(("chr1", 0, 3), ("chr2", 15, 12))
    circos.link(("chr2", 0, 3), ("chr3", 7, 11), color="skyblue")
    circos.link(("chr3", 2, 5), ("chr5", 15, 12), color="chocolate", direction=1)
    circos.link(("chr4", 3, 5), ("chr4", 18, 15), color="lime", ec="black", lw=0.5, hatch="//", direction=2)
    circos.link(("chr5", 8, 10), ("chr5", 2, 8), color="violet", ec="red", lw=1.0, ls="dashed")

    circos.savefig(f"{args.output_dir}/example.{args.output_format}")


# PyCircos method used
def get_circos(args):
    '''
    PyCircos
    :param args:
    args: parse from argparse
    :return:
    '''

    Garc = pycircos.Garc
    Gcircle = pycircos.Gcircle
    chrom_path = _resolve_input_path(args.chrom)
    cytoband_path = _resolve_input_path(args.cytoband)
    point_path = _resolve_input_path(args.point)
    line_path = _resolve_input_path(args.line)
    bar_path = _resolve_input_path(args.bar)
    heatmap_path = _resolve_input_path(args.heatmap)
    link_path = _resolve_input_path(args.link)
    automatic_track_path = _resolve_input_path(args.track)
    if automatic_track_path:
        detected_type = _track_type(automatic_track_path)
        if detected_type == "line":
            if line_path:
                raise ValueError("Specify a line track only once with -l or --line")
            line_path = automatic_track_path
        else:
            if link_path:
                raise ValueError("Specify a link track only once with -k or --link")
            link_path = automatic_track_path

    track_ranges = _allocate_track_ranges({
        "point": point_path,
        "line": line_path,
        "bar": bar_path,
        "heatmap": heatmap_path,
    })
    link_radius = min(
        (radial_range[0] for radial_range in track_ranges.values()),
        default=950,
    ) - 10

    # Set chromosomes
    circle = Gcircle()
    if chrom_path:
        chrom_sector(circle, length_path=chrom_path)
    if cytoband_path:
        cytoband(circle, cytoband_path=cytoband_path)
    if point_path:
        scatter(circle, scatter_path=point_path, raxis_range=track_ranges["point"])
    if line_path:
        line(circle, line_path=line_path, raxis_range=track_ranges["line"])
    if bar_path:
        bar(circle, bar_path=bar_path, raxis_range=track_ranges["bar"])
    if heatmap_path:
        heatmap(circle, heatmap_path=heatmap_path, raxis_range=track_ranges["heatmap"])
    if link_path:
        link(circle, link_path=link_path, radius=link_radius)

    if args.print:
        circle.figure.savefig(sys.stdout.buffer, format=args.output_format)
        sys.stdout.buffer.flush()
    else:
        output_path = Path(args.output_file).expanduser()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        circle.figure.savefig(str(output_path))
