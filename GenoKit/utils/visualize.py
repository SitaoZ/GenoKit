# -*- coding: utf-8 -*-
'''
filename: extract_uorf.py
date: 2025/12/12 15:18
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.ticker import FuncFormatter, MaxNLocator

mpl.rcParams['axes.formatter.useoffset'] = False

class visual_transcript(object):
    '''
    Visualize the uORF/dORF in transcript location (relative position) 
    Parameters:
        output: output file path, str
        transcript_id: str
        transcript_location: a 1D list postion of transcript
        cds_location: a 1D list CDS postion of transcript
        xorf_location: a 2D list position of uORFs/dORFs
        feature_type: str dORF or uORF
    Return:
     - a figure object
    '''
    def __init__(self, output, transcript_id, transcript_location,
                 cds_location, xorf_location, feature_type):
        # init
        self.output = output                  # output prefix of figure
        self.t_id = transcript_id             
        self.transcript = transcript_location # 1D list
        self.CDS = cds_location               # 1D list
        self.xORFs = xorf_location            # 2D list dORF[start, end]
        self.feature_type = feature_type      # uorf types: type 1 2 3
        # secondary 
        self.xt = self.transcript[0] # real position 
        self.xc = self.CDS[0]
        # figure global variable
        self.xticks = transcript_location
        self.yt = 10
        self.yc = 10
        self.height = 1
        self.scale = 0.5 
        self.tip = 1
        
    def draw(self):
        if len(self.transcript) != 2 or self.transcript[1] <= self.transcript[0]:
            raise ValueError('transcript_location must contain an increasing start/end pair')

        # Transcript、5'UTR、CDS和3'UTR属于同一条转录本结构轨道；
        # uORF/dORF才使用下面的独立轨道。
        track_count = 1 + len(self.xORFs)
        figure_height = max(4.2, 1.8 + track_count * 0.62)
        fig, ax = plt.subplots(figsize=(12, figure_height), dpi=150)

        transcript_y = track_count - 1
        xorf_y_values = [track_count - 2 - index
                         for index in range(len(self.xORFs))]
        block_height = 0.26

        # 先绘制完整的紫色转录本，再在同一轨道叠加橙色CDS。
        # 因此CDS两侧未被覆盖的紫色区段分别表示5'UTR和3'UTR。
        ax.add_patch(Rectangle(
            (self.transcript[0], transcript_y - block_height / 2),
            self.transcript[1] - self.transcript[0], block_height,
            facecolor='#DA70D6', edgecolor='#A64CAE', linewidth=0.9,
            zorder=2
        ))
        ax.add_patch(Rectangle(
            (self.CDS[0], transcript_y - block_height / 2),
            self.CDS[1] - self.CDS[0], block_height,
            facecolor='#FF8C00', edgecolor='#D87300', linewidth=0.9,
            zorder=3
        ))

        for xorf, y_position in zip(self.xORFs, xorf_y_values):
            start, end = xorf
            ax.add_patch(Rectangle(
                (start, y_position - block_height / 2),
                end - start, block_height,
                facecolor='#43CD80', edgecolor='#2D9E5F', linewidth=0.8,
                zorder=2
            ))

        # 5'UTR和3'UTR标签放在转录本轨道附近，不与坐标刻度混排。
        if self.CDS[0] > self.transcript[0]:
            ax.text((self.transcript[0] + self.CDS[0]) / 2,
                    transcript_y + 0.28, "5'UTR",
                    ha='center', va='bottom', fontsize=8, color='#555555')
        if self.CDS[1] < self.transcript[1]:
            ax.text((self.CDS[1] + self.transcript[1]) / 2,
                    transcript_y + 0.28, "3'UTR",
                    ha='center', va='bottom', fontsize=8, color='#555555')
        ax.text((self.CDS[0] + self.CDS[1]) / 2,
                transcript_y + 0.28, 'CDS',
                ha='center', va='bottom', fontsize=8, color='#555555')

        y_positions = [transcript_y] + xorf_y_values
        y_labels = [
            f'Transcript\n{self.t_id}'
        ] + [f'{self.feature_type} {index + 1}'
             for index in range(len(self.xORFs))]
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=9)
        ax.tick_params(axis='y', length=0, pad=10)

        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_major_locator(MaxNLocator(nbins=7, integer=True, min_n_ticks=4))
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda value, position: f'{int(round(value)):,}')
        )
        ax.set_xlabel('Transcript position (nt)', fontsize=10, labelpad=12)
        ax.tick_params(axis='x', labelsize=9, pad=6, length=5,
                       color='#707070', labelcolor='#404040')

        coordinate_span = self.transcript[1] - self.transcript[0]
        padding = max(coordinate_span * 0.025, 5)
        ax.set_xlim(self.transcript[0] - padding,
                    self.transcript[1] + padding)
        ax.set_ylim(-0.7, track_count - 0.15)

        for spine_name in ('left', 'right', 'bottom'):
            ax.spines[spine_name].set_visible(False)
        ax.spines['top'].set_color('#A0A0A0')
        ax.grid(axis='x', color='#E6E6E6', linewidth=0.7, linestyle='--')
        ax.set_axisbelow(True)

        fig.tight_layout(pad=1.4)
        output_path = f'{self.output}.pdf'
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        return fig



class visual_transcript_genome(object):
    '''
    Visualize a transcript on genome, (absolute position)
    Parameters:
        strand: str
        output: output file path, str
        transcript_id: transcript id, str
        exon_location: exon location blocks in genome, 2D list
        cds_location: CDS location blocks in genome, 2D list
        xorf_location: uORF/dORF location blocks in genome, 2D list
        feature_type: dORF or uORF
    Return:
        a figure object
    '''
    def __init__(self, strand, output, transcript_id,
                 exon_location, cds_location,
                 xorf_location, feature_type):
        # init
        self.strand = strand               # strand in genome
        self.output = output               # output prefix of a figure
        self.transcript_id = transcript_id
        self.exon_location = exon_location # 2D list 
        # print('exon_location=',exon_location)
        self.cds_location = cds_location   # 2D list 
        # print('cds_location=',cds_location)
        self.xorf_location = xorf_location # 3D list  
        # print("xorf_location=", xorf_location)
        self.feature_type = feature_type   # uorf types: type 1 2 3
        # figure global variable
        self.height = 1
        self.yt = 10                       # relative position
        self.yc = 10
        self.scale = 0.8
        self.tip = 1
        self.xticks = exon_location
        # intron from exon location 
        self.transcript_intron_location = self.transcript_intron_locate()
    
    def transcript_intron_locate(self):
        '''
        intron of transcript in genome coordinate, same to genome location 
        -1 
        '''
        intron_list = []
        intron_start = 0
        intron_end = 0
        # no intron 
        if len(self.exon_location) == 1:
            return None
        # intron exist 
        else:
            for i,exon in enumerate(self.exon_location):
                i = i + 1
                if i == 1 :
                    # first exon 
                    intron_start = exon[1]-1 # genome location should add 1 exon[1]+1; visualization ignore this
                    intron_list.append(intron_start)
                    continue 
                elif i == len(self.exon_location):
                    # last exon 
                    intron_end = exon[0]-1
                    intron_list.append(intron_end)
                else:
                    # inter exon
                    intron_list.append(exon[0]-1)
                    intron_list.append(exon[1]-1)
            array = np.array(intron_list)
            intron_location_2d = array.reshape(-1, 2)
            #if self.strand == '+':
            #    genome_intron_location_2d = intron_location_2d
            #else:
            #    genome_intron_location_2d = intron_location_2d[::-1]
            return intron_location_2d

    def cds_intron_location(self):
        '''
        intron of cds in genome coordinate, same to genome location 
        '''
        intron_list = []
        intron_start = 0
        intron_end = 0
        # no intron 
        if len(self.cds_location) == 1:
            return None
        # intron exist 
        else:
            for i,exon in enumerate(self.cds_location):
                i = i + 1
                if i == 1 :
                    # first exon 
                    intron_start = exon[1]-1 # genome location should add 1 exon[1]+1; visualization ignore this
                    intron_list.append(intron_start)
                    continue
                elif i == len(self.cds_location):
                    # last exon 
                    intron_end = exon[0]-1
                    intron_list.append(intron_end)
                else:
                    # inter exon
                    intron_list.append(exon[0]-1)
                    intron_list.append(exon[1]-1)
            array = np.array(intron_list)
            intron_location_2d = array.reshape(-1, 2)
            #if self.strand == '+':
            #    genome_intron_location_2d = intron_location_2d
            #else:
            #    genome_intron_location_2d = intron_location_2d[::-1]
            return intron_location_2d 
    
    def xorf_intron_location(self):
        '''
         uorf/dorf intron location  
        '''
        intron_list_all = []
        intron_start = 0
        intron_end = 0
        # 2D list 
        for xorf in self.xorf_location:
            # wotyout intron 
            # print('xorf:', xorf)
            if len(xorf) == 1:
                intron_list_all.append([])
            else:
                # with intron 
                if self.strand == '-':
                    xorf = xorf[::-1] # reverse list 
                intron_list = []
                for i,exon in enumerate(xorf):
                    i = i + 1
                    if i == 1 :
                        # first exon 
                        intron_start = exon[1]-1 # genome location should add 1 exon[1]+1; visualization ignore this
                        intron_list.append(intron_start)
                        continue
                    elif i == len(xorf):
                        # last exon 
                        intron_end = exon[0]-1
                        intron_list.append(intron_end)
                    else:
                        # inter exon
                        intron_list.append(exon[0]-1)
                        intron_list.append(exon[1]-1)
                array = np.array(intron_list)
                intron_location_2d = array.reshape(-1, 2)
                intron_list_all.append(intron_location_2d)    
                #if self.strand == '+':
                #    genome_intron_location_2d = intron_location_2d
                #else:
                #    genome_intron_location_2d = intron_location_2d[::-1]
        return intron_list_all
        
    def get_arrow_ab(self, exon):
        """get arrow posA and posB"""
        if self.strand == '+':
            pos_a = exon[0]
            pos_b = exon[1]
        else:
            pos_a = exon[1]
            pos_b = exon[0]
        return pos_a, pos_b
         
    def draw(self):
        if not self.exon_location:
            raise ValueError('exon_location is empty; cannot draw transcript structure')

        # Transcript外显子、UTR和CDS绘制在同一轨道；uORF/dORF使用独立轨道。
        # 画布高度随uORF/dORF数量增长，宽度固定为12英寸。
        track_count = 1 + len(self.xorf_location)
        figure_height = max(4.2, 1.8 + track_count * 0.62)
        fig, ax = plt.subplots(figsize=(12, figure_height), dpi=150)

        transcript_y = track_count - 1
        xorf_y_values = [track_count - 2 - index
                         for index in range(len(self.xorf_location))]
        block_height = 0.26

        def draw_blocks(blocks, y_position, color, edge_color):
            """绘制uORF/dORF区块及跨外显子的内含子连接线。"""
            sorted_blocks = sorted(blocks, key=lambda block: (block[0], block[1]))
            for left, right in zip(sorted_blocks, sorted_blocks[1:]):
                ax.plot([left[1], right[0]], [y_position, y_position],
                        color='#808080', linewidth=0.9, zorder=1)
            for block in sorted_blocks:
                ax.add_patch(Rectangle(
                    (block[0], y_position - block_height / 2),
                    block[1] - block[0], block_height,
                    facecolor=color, edgecolor=edge_color,
                    linewidth=0.8, zorder=2
                ))

        # 先画紫色外显子和灰色内含子连接线，再将橙色CDS区块叠加到
        # 同一条转录本轨道。未被CDS覆盖的紫色部分即为5'UTR和3'UTR。
        sorted_exons = sorted(
            self.exon_location, key=lambda block: (block[0], block[1])
        )
        for left, right in zip(sorted_exons, sorted_exons[1:]):
            ax.plot([left[1], right[0]], [transcript_y, transcript_y],
                    color='#808080', linewidth=1.0, zorder=1)
        for start, end in sorted_exons:
            ax.add_patch(Rectangle(
                (start, transcript_y - block_height / 2),
                end - start, block_height,
                facecolor='#DA70D6', edgecolor='#A64CAE',
                linewidth=0.9, zorder=2
            ))
        for start, end in sorted(self.cds_location,
                                 key=lambda block: (block[0], block[1])):
            ax.add_patch(Rectangle(
                (start, transcript_y - block_height / 2),
                end - start, block_height,
                facecolor='#FF8C00', edgecolor='#D87300',
                linewidth=0.9, zorder=3
            ))

        for blocks, y_position in zip(self.xorf_location, xorf_y_values):
            draw_blocks(blocks, y_position, '#43CD80', '#2D9E5F')

        genome_start = min(block[0] for block in self.exon_location)
        genome_end = max(block[1] for block in self.exon_location)
        if self.cds_location:
            cds_start = min(block[0] for block in self.cds_location)
            cds_end = max(block[1] for block in self.cds_location)
            if self.strand == '+':
                utr5_range = (genome_start, cds_start)
                utr3_range = (cds_end, genome_end)
            else:
                utr5_range = (cds_end, genome_end)
                utr3_range = (genome_start, cds_start)

            for label, (start, end) in (
                    ("5'UTR", utr5_range), ("3'UTR", utr3_range)):
                if end > start:
                    ax.text((start + end) / 2, transcript_y + 0.28, label,
                            ha='center', va='bottom', fontsize=8,
                            color='#555555')
            ax.text((cds_start + cds_end) / 2, transcript_y + 0.28, 'CDS',
                    ha='center', va='bottom', fontsize=8,
                    color='#555555')

        # 左侧轨道名称替代原来挤在每个区块末端的小字体标签。
        y_positions = [transcript_y] + xorf_y_values
        y_labels = [
            f'Transcript\n{self.transcript_id} ({self.strand})'
        ] + [f'{self.feature_type} {index + 1}'
             for index in range(len(self.xorf_location))]
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=9)
        ax.tick_params(axis='y', length=0, pad=10)

        # 使用真正的顶部坐标轴。MaxNLocator限制主刻度数量，千位分隔符和
        # 30度旋转保证八到十位的基因组坐标仍可区分。
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True, min_n_ticks=4))
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda value, position: f'{int(round(value)):,}')
        )
        ax.set_xlabel('Genomic position (bp)', fontsize=10, labelpad=12)
        ax.tick_params(axis='x', labelsize=9, pad=6, length=5,
                       color='#707070', labelcolor='#404040')
        for tick_label in ax.get_xticklabels():
            tick_label.set_rotation(30)
            tick_label.set_horizontalalignment('left')

        coordinate_span = max(genome_end - genome_start, 1)
        padding = max(coordinate_span * 0.025, 20)
        ax.set_xlim(genome_start - padding, genome_end + padding)
        ax.set_ylim(-0.7, track_count - 0.15)

        for spine_name in ('left', 'right', 'bottom'):
            ax.spines[spine_name].set_visible(False)
        ax.spines['top'].set_color('#A0A0A0')
        ax.grid(axis='x', color='#E6E6E6', linewidth=0.7, linestyle='--')
        ax.set_axisbelow(True)

        fig.tight_layout(pad=1.4)
        output_path = f'{self.output}.pdf'
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        return fig



if __name__ == "__main__":
    # a = visual_transcript('out', 'aaa',[1,1000],\
    #                        [100,800],\
    #                        [[30,90],[60,440]], 'type1')
    # a.draw()

    b = visual_transcript_genome('+', 'out', 'aaa',[[1,100],[200,300],[340,400],[440,550]], \
                                 [[50,100],[200,300],[340,400],[440,500]], \
                                 [[[20,100],[200,240]],[[20,100],[200,300],[340,360]],[[50,100],[200,300]]],'type')
    b.draw()
