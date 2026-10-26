# -*- coding: utf-8 -*-

'''
filename: design_sirna.py
date: 2025/12/8 10:22
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

from typing import List
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

class GeneStructurePlotter:
    def __init__(self, transcripts: List, gene: str, orf_label_limit: int = 20):
        self.transcripts = transcripts
        self.gene = gene
        self.orf_label_limit = max(0, orf_label_limit)
        self.min_pos = None
        self.max_pos = None

    def _calculate_global_range(self):
        """计算基因组坐标范围"""
        self.min_pos = None
        self.max_pos = None
        all_pos = []
        for transcript in self.transcripts:
            features = []
            for region in ['exons']:
                if transcript.get(region):
                    for exon in transcript[region]:
                        features.extend(exon)
            for region in ['utr5', 'utr3']:
                if transcript.get(region):
                    for exon in transcript[region]:
                        features.extend(exon)
            for exon in transcript['cds_exons']:
                features.extend(exon)
            for orf_type in ['uorfs', 'dorfs']:
                for orf in transcript[orf_type]:
                    for exon in orf:
                        features.extend(exon)
            
            if features:
                current_min = min(features)
                current_max = max(features)
                self.min_pos = min(self.min_pos, current_min) if self.min_pos else current_min
                self.max_pos = max(self.max_pos, current_max) if self.max_pos else current_max

    def _draw_directional_block(self, ax, start, end, y_center, height, color, strand):
        """绘制带方向箭头的区块"""
        block_width = abs(end - start)
        if block_width == 0:
            return

        # x 轴为基因组坐标，不能用 y 轴高度直接作为箭头宽度。
        # 箭头最多占区块的 35%，短 ORF 也不会出现负宽度矩形。
        genomic_span = max((self.max_pos or end) - (self.min_pos or start), 1)
        arrow_size = min(block_width * 0.35,
                         max(block_width * 0.15, genomic_span * 0.0015))
        
        # 绘制主体矩形
        if strand == '+':
            rect_start = start
            arrow_start = end - arrow_size
        else:
            rect_start = start + arrow_size
            arrow_start = start
        
        rect = Rectangle((rect_start, y_center - height/2),
                        block_width - arrow_size,
                        height,
                        facecolor=color['fill'],
                        edgecolor=color['edge'],
                        lw=0.8)
        ax.add_patch(rect)
        
        # 绘制方向箭头
        if strand == '+':
            arrow_points = [
                (end - arrow_size, y_center - height/2),
                (end, y_center),
                (end - arrow_size, y_center + height/2)
            ]
        else:
            arrow_points = [
                (start + arrow_size, y_center - height/2),
                (start, y_center),
                (start + arrow_size, y_center + height/2)
            ]
        
        arrow = Polygon(arrow_points,
                       closed=True,
                       facecolor=color['fill'],
                       edgecolor=color['edge'])
        ax.add_patch(arrow)

    def _adjust_coordinates(self, coords, strand):
        """根据链方向调整坐标顺序"""
        # print("coords, strand:", self.transcripts, coords, strand)
        return sorted(coords, reverse=(strand == '-'))

    @staticmethod
    def _assign_feature_lanes(features):
        """按基因组区间分配不重叠轨道，返回绘图任务和所需轨道数。"""
        indexed_features = []
        for feature_index, exons in enumerate(features):
            if not exons:
                continue
            start = min(block[0] for block in exons)
            end = max(block[1] for block in exons)
            indexed_features.append((start, end, feature_index, exons))

        indexed_features.sort(key=lambda item: (item[0], item[1], item[2]))
        lane_ends = []
        assignments = []
        for start, end, feature_index, exons in indexed_features:
            lane_index = None
            for candidate, lane_end in enumerate(lane_ends):
                if start > lane_end:
                    lane_index = candidate
                    lane_ends[candidate] = end
                    break
            if lane_index is None:
                lane_index = len(lane_ends)
                lane_ends.append(end)
            assignments.append((feature_index, exons, lane_index))
        return assignments, len(lane_ends)

    def _draw_multi_exon_orf(self, ax, orf_exons, y_level, color, strand,
                             orf_index, show_label=True):
        """绘制多外显子ORF"""
        sorted_exons = self._adjust_coordinates(orf_exons, strand)
        # orf_height = 0.25
        orf_height = 0.15
        
        # 绘制所有外显子
        for idx, (start, end) in enumerate(sorted_exons):
            self._draw_directional_block(ax, start, end, y_level, orf_height, color, strand)
            
            # 绘制连接线
            if idx < len(sorted_exons) - 1:
                if strand == "+":
                    next_start = sorted_exons[idx+1][0]
                    ax.hlines(y_level, end, next_start, 
                         colors=color['edge'], 
                         linestyles=':', 
                         linewidth=0.8)
                else:
                    next_start = sorted_exons[idx+1][1]
                    ax.hlines(y_level, next_start, start,
                         colors=color['edge'], 
                         linestyles=':', 
                         linewidth=0.8)
        
        # 添加ORF编号
        if sorted_exons and show_label:
            first_start, first_end = sorted_exons[0]
            ax.text((first_start + first_end)/2, y_level, 
                   str(orf_index+1), 
                   ha='center', va='center',
                   fontsize=7, color=color['edge'], weight='bold')

    def _draw_stranded_features(self, ax, features, base_y, region_height, orf_type, strand, vertical_direction):
        """绘制带方向的多外显子ORF"""
        color = self.color_palette['uorf'] if orf_type == 'uorf' else self.color_palette['dorf']
        assignments, lane_count = self._assign_feature_lanes(features)
        if lane_count == 0:
            return

        layer_spacing = region_height / lane_count
        show_labels = len(features) <= self.orf_label_limit
        for orf_idx, orf_exons, lane in assignments:
            y_level = base_y + vertical_direction * (lane + 0.5) * layer_spacing
            self._draw_multi_exon_orf(
                ax, orf_exons, y_level, color, strand, orf_idx,
                show_label=show_labels
            )

    def _draw_intron(self, ax, start, end, y, style, strand):
        """绘制内含子连接线"""
        # sorted_coords = sorted([start, end], reverse=(strand == '-'))
        sorted_coords = sorted([start, end])
        if sorted_coords[0] < sorted_coords[1]:
            ax.hlines(y, sorted_coords[0], sorted_coords[1], 
                     colors=style['color'],
                     linestyles=style['ls'],
                     linewidth=style['lw'])

    def _create_legend(self, ax):
        """创建外置图例"""
        legend_elements = [
            Rectangle((0,0),1,1, 
                     facecolor=self.color_palette['cds']['fill'], 
                     edgecolor=self.color_palette['cds']['edge'],
                     label='CDS'),
            Rectangle((0,0),1,1,
                     facecolor=self.color_palette['utr']['fill'],
                     edgecolor=self.color_palette['utr']['edge'],
                     label='UTR'),
            Rectangle((0,0),1,1,
                     facecolor=self.color_palette['noncoding']['fill'],
                     edgecolor=self.color_palette['noncoding']['edge'],
                     label='Non-coding Exon'),  # 新增图例项 Non-coding
            Rectangle((0,0),1,1,
                     facecolor=self.color_palette['uorf']['fill'],
                     edgecolor=self.color_palette['uorf']['edge'],
                     label='uORF'),
            Rectangle((0,0),1,1,
                     facecolor=self.color_palette['dorf']['fill'],
                     edgecolor=self.color_palette['dorf']['edge'],
                     label='dORF'),
            plt.Line2D([0], [0], 
                      color=self.color_palette['intron']['color'],
                      linestyle=self.color_palette['intron']['ls'],
                      linewidth=self.color_palette['intron']['lw'],
                      label='Intron'),
            plt.Line2D([0], [0], 
                      color=self.color_palette['uorf']['edge'],
                      linestyle=':',
                      linewidth=0.8,
                      label='Connection')
        ]
        # 将图例放在右侧外部
        leg = ax.legend(handles=legend_elements, 
                      loc='upper left',
                      bbox_to_anchor=(1.02, 1),
                      frameon=True,
                      framealpha=0.95,
                      edgecolor='#FFFFFF',
                      fontsize=9,
                      title='Legend:',
                      title_fontsize=10,
                      alignment = 'left')
        return leg

    def plot_transcripts(self):
        """主绘图方法"""
        if not self.transcripts:
            raise ValueError("没有可用于绘图的转录本。")

        fig_width = 14
        lane_height = 0.26
        transcript_layouts = []
        for transcript in self.transcripts:
            _, uorf_lanes = self._assign_feature_lanes(transcript.get('uorfs', []))
            _, dorf_lanes = self._assign_feature_lanes(transcript.get('dorfs', []))
            upper_region = max(0.26, uorf_lanes * lane_height)
            lower_region = max(0.26, dorf_lanes * lane_height)
            transcript_layouts.append({
                'upper_region': upper_region,
                'lower_region': lower_region,
                'height': 1.0 + upper_region + lower_region
            })

        fig_height = max(4.0, 1.8 + sum(item['height'] for item in transcript_layouts))

        # fig_height = 2 + 1.8 * len(self.transcripts) + extra_height  # 动态计算高度
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=120)
        self._calculate_global_range()
        
        # 专业配色方案
        self.color_palette = {
            'cds': {'fill': '#4C72B0', 'edge': '#2F4F7F'},
            'utr': {'fill': '#E0E0E0', 'edge': '#A0A0A0'},
            'uorf': {'fill': '#59A14F', 'edge': '#36752F'},
            'dorf': {'fill': '#E15759', 'edge': '#B83A3C'},
            'intron': {'color': '#606060', 'ls': (0, (5, 2)), 'lw': 1.2},
            'noncoding': {'fill': '#9467BD', 'edge': '#6A3D9A'}  # 新增非编码转录本颜色
        }
        
        # 计算每个转录本的y轴位置
        y_positions = []
        current_y = 0
        for transcript_layout in transcript_layouts:
            y_positions.append(current_y + transcript_layout['lower_region'] + 0.5)
            current_y += transcript_layout['height']

        for y_idx, (transcript,y_center) in enumerate(zip(self.transcripts, y_positions)):
            strand = transcript.get('strand', '+')
            transcript_layout = transcript_layouts[y_idx]
            has_cds = 'cds_exons' in transcript and transcript['cds_exons']
            layout = {
                'cds_height': 0.9,
                'utr_height': 0.6,
                'uorf_region': transcript_layout['upper_region'],
                'dorf_region': transcript_layout['lower_region'],
                'noncoding_height': 0.4  # 非编码外显子高度
            }
            
            # 处理5'UTR
            if transcript.get('utr5'):
                utr5 = self._adjust_coordinates([transcript['utr5']], strand)[0]
                for start, end in utr5:
                    self._draw_directional_block(ax, start, end, y_center, layout['utr_height'], 
                                            self.color_palette['utr'], strand)
                # 处理uORF
                self._draw_stranded_features(ax, transcript['uorfs'], 
                                            y_center + layout['utr_height']/2,
                                            layout['uorf_region'], 'uorf', strand, 1)

            # 处理CDS
            if has_cds:
                # 编码转录本：绘制CDS
                sorted_exons = self._adjust_coordinates(transcript['cds_exons'], strand)
                for start, end in sorted_exons:
                    self._draw_directional_block(ax, start, end, y_center, 
                                                layout['cds_height'],
                                                self.color_palette['cds'], strand)
            elif 'exons' in transcript:
                # 非编码转录本：绘制所有外显子
                sorted_exons = self._adjust_coordinates(transcript['exons'], strand)
                for start, end in sorted_exons:
                    self._draw_directional_block(ax, start, end, y_center, 
                                                 layout['noncoding_height'],
                                                 self.color_palette['noncoding'], strand)
                
            # 处理3'UTR
            if transcript.get('utr3'):
                utr3 = self._adjust_coordinates([transcript['utr3']], strand)[0]
                for start, end in utr3:
                    self._draw_directional_block(ax, start, end, y_center, layout['utr_height'],
                                            self.color_palette['utr'], strand)
                # 处理dORF
                self._draw_stranded_features(ax, transcript['dorfs'],
                                            y_center - layout['utr_height']/2,
                                            layout['dorf_region'], 'dorf', strand, -1)
            # 绘制内含子 utr5 + cds + utr3
            
            # intron_region = utr5 + sorted_exons + utr3
            # intron_region_sort = self._adjust_coordinates(intron_region, strand)
            # if strand == '-':
            #    # 负链也正向连线
            #    intron_region_sort = intron_region_sort[::-1]
            # for i in range(len(intron_region_sort)-1):
            #     self._draw_intron(ax, intron_region_sort[i][1], intron_region_sort[i+1][0],
            #                      y_center, self.color_palette['intron'], strand)
            # 绘制内含子
            intron_region = []
            if transcript.get('utr5'):
                intron_region.extend(utr5)
            if has_cds:
                intron_region.extend(sorted_exons)
            elif 'exons' in transcript:
                intron_region.extend(sorted_exons)
            if transcript.get('utr3'):
                intron_region.extend(utr3)
                
            if intron_region:  # 只有在外显子/UTR存在时才绘制内含子
                intron_region_sort = self._adjust_coordinates(intron_region, strand)
                if strand == '-':
                    intron_region_sort = intron_region_sort[::-1]
                for i in range(len(intron_region_sort)-1):
                    self._draw_intron(ax, intron_region_sort[i][1], intron_region_sort[i+1][0],
                                    y_center, self.color_palette['intron'], strand)

        # 坐标轴设置
        coordinate_padding = max((self.max_pos - self.min_pos) * 0.01, 50)
        ax.set_xlim(self.min_pos - coordinate_padding, self.max_pos + coordinate_padding)
        ax.set_ylim(-0.15, current_y + 0.15)
        ax.set_yticks(y_positions)
        ax.set_yticklabels([
            (f"{t.get('name', 'Transcript')} ({t.get('strand', '+')})"
             f"  uORF:{len(t.get('uorfs', []))}  dORF:{len(t.get('dorfs', []))}")
            for t in self.transcripts],
                          fontsize=9)  # 减小字体
        ax.set_xlabel('Genomic Position (bp)', fontsize=10)
        
        # 添加图例并调整布局
        self._create_legend(ax)
        plt.subplots_adjust(right=0.8)  # 为图例留出空间
        ax.set_title(f'Gene Structure of {self.gene} with Adaptive Layout\n',
                     fontsize=12, pad=15)
        fig.tight_layout()
        # plt.show()
        # fig.savefig("foo.pdf")
        return fig

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
        },
        {
            'name': 'Minus4',
            'strand': '-',
            'utr5': [],
            'cds_exons': [],
            'exons': [(670, 750), (570, 660)],
            'utr3': [],
            'uorfs': [],
            'dorfs': []
        },
    ]
    
    plotter = GeneStructurePlotter(test_transcripts)
    plotter.plot_transcripts()
