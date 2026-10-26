# -*- coding: utf-8 -*-
'''
filename: extract_view.py
date: 2025/11/28 14:46
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import os
import sys
import logging
from tqdm import tqdm
from Bio.Seq import Seq
from pathlib import Path
import matplotlib.pyplot as plt
from collections import defaultdict
from GenoKit.database.database import genome_dict
from GenoKit.database.annotation_db import load_annotation_database
from GenoKit.visualize.view_gene import GeneStructurePlotter
from GenoKit.extract.extract_uorf import GFF, uORF, get_mt_cds_seq
from GenoKit.extract.extract_dorf import dORF
from GenoKit.extract.utr_calculator import UTRCalculator


def parse_transcript(params):
    """
    parameters:
        params: argparse
    return:
        transcript_seq_list, a deque
    """
    gene_isoform = defaultdict(lambda: defaultdict(dict))
    g, ts, db, genome, style, output_format = params
    for t in ts:
        seq = ''
        mRNA = db[g][t]['mRNA']
        for e in db[g][t]['exon']:
            chrom, start, end, strand = e.chr, e.start, e.end, e.strand
            s = str(genome[chrom][start-1:end])
            seq += s
        seq = Seq(seq)
        if mRNA.strand == '-':
            seq = seq.reverse_complement()
        seq = str(seq)
        gene_isoform[g][t] = seq
    return gene_isoform
        

def feature_genome_location(db, g, t2g, genome, style, orf_length):
    '''
    the genome location of gene feature
    :param db:
    :param g:
    :param t2g:
    :param genome:
    :param style:
    :param orf_length:
    :return:
    '''
    ts = [i for i in db[g] if i not in {'gene', 'symbol'}]
    g = t2g[ts[0]]
    test_transcripts = []
    for t in ts:
        mRNA = db[g][t]['mRNA']
        # specify the transcript id; only one transcript
        mt, cds, exons_list, cds_list = get_mt_cds_seq(g, t, db, genome, style)
        if len(cds_list) == 0:
            logging.info(f"The transcript {t} of {g} gene does not have CDS region.")
        uORF_ = uORF(g, t, mRNA.chr, mRNA.strand, mt, cds)
        # figure the schametic
        if len(cds) % 3 != 0:
            logging.info(
                "Skip uORF visualization for %s: CDS length %d is not a multiple of 3.",
                t, len(cds)
            )
            uORF_dict = {}
        else:
            uORF_dict = uORF_.uorf_parse()
        # print('uORF_dict:', uORF_dict)
        # print("---------------------------exons_list", exons_list)
        # print("---------------------------uORF_dict:", uORF_dict)
        uorf_location_genome = []  # for schematic on genome # 3D list
        for key in sorted(uORF_dict.keys()):
            for it in uORF_dict[key]:
                # uORF length filter
                if len(it[9]) <= orf_length:
                    # default: 6
                    continue
                # for schematic on genome
                gff = GFF(exons_list, it, 'uORF')
                ex_locations = gff.uorf_exons_location()
                uorf_location_genome.append(ex_locations)
        # if t == 'ENST00000508070':
        #     # ENST00000508070是HMGCR 其中的一条没有编码CDS的转录本
        #     print(t, 'zhusitao', uorf_location_genome)
        dORF_ = dORF(t, orf_length, mRNA.chr, mRNA.strand, mt, cds)
        dORF_dict = dORF_.dorf_parse()
        dorf_location_genome = []  # for schematic on genome # 3D list
        for key in sorted(dORF_dict.keys()):
            for it in dORF_dict[key]:
                # dORF length filter；与 GenoKit dorf 的输出规则保持一致。
                if len(it[9]) < orf_length:
                    # default: 6
                    continue
                # for schematic on genome
                gff = GFF(exons_list, it, 'dORF')
                ex_locations = gff.uorf_exons_location()
                dorf_location_genome.append(ex_locations)

                # utr5 utr3 calculator
        # print(f"{t} exons_list:", exons_list)
        # print(f"{t} cds_list:", cds_list)
        # if t == 'ENST00000508070':
        #    print(t, 'fanfei', dorf_location_genome)
        calculator = UTRCalculator(exons_list, cds_list, mRNA.strand)
        utr_regions = calculator.get_utr_regions()
        # print(f"{t} utr:", utr_regions)
        # print(f"{t} uorf:", uorf_location_genome)
        # print(f"{t} dorf:", dorf_location_genome)
        # Dict type
        test_transcripts.append({'name': t,
                                 'strand': mRNA.strand,
                                 'exons': exons_list,
                                 'utr5': utr_regions['utr5'],
                                 'cds_exons': cds_list,
                                 'utr3': utr_regions['utr3'],
                                 'uorfs': uorf_location_genome,
                                 'dorfs': dorf_location_genome,
                                 })
        return test_transcripts

def get_view(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''

    output_dir = Path(args.output_dir).expanduser() if args.output_dir else None
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        log_path = output_dir / "View.log"
    else:
        log_path = Path("View.log")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        # filename="View.log"
        filename=str(log_path),
    )

    db, t2g = load_annotation_database(args.database)
    genome = genome_dict(args.genome)
    index = 0
    if not args.output_dir:
        print("No output directory specified. Please provide one and rerun")
        sys.exit(1)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if not args.gene:
        logging.info(f"Perform vision gene structure at the whole-genome scale")
        for g in tqdm(db, ncols=80, total=len(db),
                      desc = 'Vision processing:', colour="GREEN"):
            gene = db[g]['gene']
            test_transcripts = feature_genome_location(db, g, t2g, genome, args.style, args.length)
            plotter = GeneStructurePlotter(
                test_transcripts, g,
                orf_label_limit=getattr(args, 'orf_label_limit', 20)
            )
            fig = plotter.plot_transcripts()

            fig.savefig(f"{args.output_dir}/{g}.{args.output_format}")
            plt.close(fig)
            index += 1
            if index == 10:
                break
        logging.info(f"Finished vision gene structure at the whole-genome scale")
    else:
        logging.info(f"Perform vision gene structure at the gene scale")
        ts = [i for i in db[args.gene] if i not in {'gene', 'symbol'}]
        g = t2g[ts[0]]
        test_transcripts = []
        for t in ts:
            mRNA = db[g][t]['mRNA']
            # specify the transcript id; only one transcript
            params = (g, t, db, genome, args.style, args.length, args.output_format)
            mt, cds, exons_list, cds_list = get_mt_cds_seq(g, t, db, genome, args.style)
            # if t == 'ENST00000456328':
            #     print(t, mt, cds, exons_list, cds_list)
            if len(cds_list) == 0 :
                logging.info(f"The transcript {t} of {g} gene does not have CDS region.")
            uORF_ = uORF(g, t, mRNA.chr, mRNA.strand, mt, cds)
            # figure the schametic 
            if len(cds) % 3 != 0:
                logging.info(
                    "Skip uORF visualization for %s: CDS length %d is not a multiple of 3.",
                    t, len(cds)
                )
                uORF_dict = {}
            else:
                uORF_dict = uORF_.uorf_parse()
            # print('uORF_dict:', uORF_dict)
            # print("---------------------------exons_list", exons_list)
            # print("---------------------------uORF_dict:", uORF_dict)
            uorf_location_genome = [] # for schematic on genome # 3D list 
            for key in sorted(uORF_dict.keys()):
                for it in uORF_dict[key]:
                    # uORF length filter 
                    if len(it[9]) <= args.length:
                        # default: 6
                        continue
                    # for schematic on genome
                    gff = GFF(exons_list, it, 'uORF')
                    ex_locations = gff.uorf_exons_location()
                    uorf_location_genome.append(ex_locations)
            # if t == 'ENST00000508070':
            #     # ENST00000508070是HMGCR 其中的一条没有编码CDS的转录本
            #     print(t, 'zhusitao', uorf_location_genome)
            dORF_ = dORF(g, t, args.length, mRNA.chr, mRNA.strand, mt, cds)
            dORF_dict = dORF_.dorf_parse()
            dorf_location_genome = [] # for schematic on genome # 3D list 
            for key in sorted(dORF_dict.keys()):
                for it in dORF_dict[key]:
                    # dORF length filter；-l N 保留长度等于 N 的 dORF。
                    if len(it[9]) <= args.length:
                        # default: 6
                        continue
                    # for schematic on genome
                    gff = GFF(exons_list, it, 'dORF')
                    ex_locations = gff.uorf_exons_location()
                    dorf_location_genome.append(ex_locations) 
                
            # utr5 utr3 calculator
            # print(f"{t} exons_list:", exons_list)
            # print(f"{t} cds_list:", cds_list)
            #if t == 'ENST00000508070':
            #    print(t, 'fanfei', dorf_location_genome)
            calculator = UTRCalculator(exons_list, cds_list, mRNA.strand)
            utr_regions = calculator.get_utr_regions()
            # print(f"{t} utr:", utr_regions)
            # print(f"{t} uorf:", uorf_location_genome)
            # print(f"{t} dorf:", dorf_location_genome)
            # Dict type
            test_transcripts.append({'name': t, 
                                'strand': mRNA.strand, 
                                'exons': exons_list,
                                'utr5': utr_regions['utr5'],
                                'cds_exons': cds_list ,
                                'utr3': utr_regions['utr3'], 
                                'uorfs': uorf_location_genome,
                                'dorfs': dorf_location_genome,
                               })
        # plot gene structure
        # print("test_transcripts:", test_transcripts)
        plotter = GeneStructurePlotter(
            test_transcripts, g,
            orf_label_limit=getattr(args, 'orf_label_limit', 20)
        )
        fig = plotter.plot_transcripts() 
        fig.savefig(f"{args.output_dir}/{g}.{args.output_format}")
        plt.close(fig)
        logging.info(f"Finished vision gene structure at the gene scale")

if __name__ == "__main__":
    # 测试数据（包含6个转录本）
    test_transcripts = [
        {
            'name': 'Plus1',
            'strand': '+',
            'utr5': [(100, 300)],
            'cds_exons': [(320, 400), (420, 550)],
            'utr3': [(560, 700)],
            'uorfs': [[(120, 150)]],
            'dorfs': [[(600, 620)]]
        },
        {
            'name': 'Minus1',
            'strand': '-',
            'utr5': [(650, 800)],
            'cds_exons': [(570, 650), (420, 550)],
            'utr3': [(100, 300)],
            'uorfs': [[(750, 770)]],
            'dorfs': [[(150, 170)]]
        },
        {
            'name': 'Plus2',
            'strand': '+',
            'utr5': [(200, 350)],
            'cds_exons': [(370, 450), (470, 600)],
            'utr3': [(620, 750)],
            'uorfs': [[(220, 240)]],
            'dorfs': [[(700, 720)]]
        },
        {
            'name': 'Minus2',
            'strand': '-',
            'utr5': [(700, 850)],
            'cds_exons': [(620, 700), (470, 600)],
            'utr3': [(200, 350)],
            'uorfs': [[(800, 820)]],
            'dorfs': [[(250, 270)]]
        },
        {
            'name': 'Plus3',
            'strand': '+',
            'utr5': [(200, 290), (300, 450)],
            'cds_exons': [(470, 550), (570, 700)],
            'utr3': [(720, 850)],
            'uorfs': [[(260, 290), (320, 340)]],
            'dorfs': [[(800, 820)]]
        },
        {
            'name': 'Minus3',
            'strand': '-',
            'utr5': [(760, 900),(910,930),(950, 1000)],
            'cds_exons': [(670, 750), (570, 660)],
            'utr3': [(100, 200),(300, 450)],
            'uorfs': [[(850, 870),(880,890)], [(660, 888)]],
            'dorfs': [[(350, 370)]]
        }
    ]

    plotter = GeneStructurePlotter(test_transcripts)
    plotter.plot_transcripts()
