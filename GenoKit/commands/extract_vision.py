# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
import _pickle as cPickle
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from GenoKit.database.database import create_db, genome_dict
from GenoKit.utils.util import parse_output, gff_feature_dict, gtf_feature_dict
from GenoKit.commands.vision_gene import GeneStructurePlotter
from GenoKit.commands.extract_uorf import GFF, uORF, get_mt_cds_seq
from GenoKit.commands.extract_dorf import dORF
from GenoKit.commands.utr_calculator import UTRCalculator

_CSV_HEADER = ['GeneID','TranscriptID','siRNA','GC','Strand']

def parse_transcript(params):
    """
    parameters:
        params: argparse
    return:
        transcript_seq_list, a deque
    """
    gene_isoform = defaultdict(lambda: defaultdict(dict))
    g, ts, db, genome, style, output_format = params
    transcript_seq_list = deque()
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
        

def get_vision(args):
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
    with open(args.database, 'rb') as f:
        db, t2g = cPickle.load(f)
    genome = genome_dict(args.genome)
    sirna_seq_list = []
    index = 0
    if not args.gene:
        for g in tqdm(db, ncols=80, total=len(db), desc = 'Gene processing:'):
            gene = db[g]['gene']
            if isinstance(gene, dict):
                print(gene, g)
            chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
            # 
            ts = [i for i in db[g] if i != 'gene']
            params = (g, ts, db, genome, args.style, args.output_format)
            gene_isoform = parse_transcript(params)
            test_transcripts = gene_isoform[g]
            designer = SirnaDesigner(g, test_transcripts, length = 22)
            universal = designer.design(require_all=True)
            sirna_list = designer.get_top_designs(10)
            if args.output_format == 'fasta':
                # fasta (defalut)
                desc='strand:%s start:%d end:%d length=%d'%(gene.strand,
                                                              gene.start,
                                                              gene.end,
                                                              len(seq))
                geneRecord = SeqRecord(seq, id=g, description=desc)
                gene_seq_list.append(geneRecord)
            else:
                for seq,gc in sirna_list:
                    it = [g, "|".join(ts), seq, gc, gene.strand]
                    sirna_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
        parse_output(args, sirna_seq_list)
    else:
        ts =  [i for i in db[args.gene] if i != 'gene']
        g = t2g[ts[0]]
        test_transcripts = []
        for t in ts:
            mRNA = db[g][t]['mRNA']
            # specify the transcript id; only one transcript
            params = (g, t, db, genome, args.style, args.length, args.output_format)
            mt, cds, exons_list, cds_list = get_mt_cds_seq(g, t, db, genome, args.style)
            #if t == 'ENST00000508070':
            #    print(t, mt, cds, exons_list, cds_list)
            uORF_ = uORF(t, mRNA.chr, mRNA.strand, mt, cds)
            # figure the schametic 
            uORF_dict = uORF_.uorf_parse()
            # print("---------------------------exons_list", exons_list)
            # print("---------------------------uORF_dict:", uORF_dict)
            uorf_location_genome = [] # for schematic on genome # 3D list 
            for key in sorted(uORF_dict.keys()):
                for it in uORF_dict[key]:
                    # uORF length filter 
                    if len(it[8]) <= args.length:
                        # default: 6
                        continue
                    # for schematic on genome
                    gff = GFF(exons_list, it, 'uORF')
                    ex_locations = gff.uorf_exons_location()
                    uorf_location_genome.append(ex_locations)
            #if t == 'ENST00000508070':
            #    print(t, 'zhusitao', uorf_location_genome)
            dORF_ = dORF(t, args.length, mRNA.chr, mRNA.strand, mt, cds)
            dORF_dict = dORF_.dorf_parse()
            dorf_location_genome = [] # for schematic on genome # 3D list 
            for key in sorted(dORF_dict.keys()):
                for it in dORF_dict[key]:
                    # dORF length filter 
                    if len(it[8]) <= args.length:
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
        plotter = GeneStructurePlotter(test_transcripts)
        fig = plotter.plot_transcripts() 
        fig.savefig(args.output)

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

