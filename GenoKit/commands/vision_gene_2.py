import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

class GeneStructurePlotter:
    def __init__(self, transcripts):
        self.transcripts = transcripts
        self.min_pos = None
        self.max_pos = None

    def _calculate_global_range(self):
        """计算基因组坐标范围"""
        all_pos = []
        for transcript in self.transcripts:
            features = []
            for region in ['utr5', 'utr3']:
                if transcript.get(region):
                    features.extend(transcript[region])
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
        arrow_size = height * 0.8
        block_width = abs(end - start)
        
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
        return sorted(coords, reverse=(strand == '-'))

    def _draw_multi_exon_orf(self, ax, orf_exons, y_level, color, strand, orf_index):
        """绘制多外显子ORF"""
        sorted_exons = self._adjust_coordinates(orf_exons, strand)
        orf_height = 0.25
        
        # 绘制所有外显子
        for idx, (start, end) in enumerate(sorted_exons):
            self._draw_directional_block(ax, start, end, y_level, orf_height, color, strand)
            
            # 绘制连接线
            if idx < len(sorted_exons) - 1:
                next_start = sorted_exons[idx+1][0]
                ax.hlines(y_level, end, next_start, 
                         colors=color['edge'], 
                         linestyles=':', 
                         linewidth=0.8)
        
        # 添加ORF编号
        if sorted_exons:
            first_start, first_end = sorted_exons[0]
            ax.text((first_start + first_end)/2, y_level, 
                   str(orf_index+1), 
                   ha='center', va='center',
                   fontsize=7, color=color['edge'], weight='bold')

    def _draw_stranded_features(self, ax, features, base_y, region_height, orf_type, strand, vertical_direction):
        """绘制带方向的多外显子ORF"""
        color = self.color_palette['uorf'] if orf_type == 'uorf' else self.color_palette['dorf']
        max_layers = 3
        layer_spacing = region_height / max_layers
        
        for orf_idx, orf_exons in enumerate(features):
            layer = orf_idx % max_layers
            y_level = base_y + vertical_direction * (layer + 0.5) * layer_spacing
            self._draw_multi_exon_orf(ax, orf_exons, y_level, color, strand, orf_idx)

    def _draw_intron(self, ax, start, end, y, style, strand):
        """绘制内含子连接线"""
        sorted_coords = sorted([start, end], reverse=(strand == '-'))
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
                      label='ORF Connection')
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
                      title_fontsize=10)
        return leg

    def plot_transcripts(self):
        """主绘图方法"""
        plt.style.use('seaborn-whitegrid')
        fig_width = 14
        fig_height = 2 + 1.8 * len(self.transcripts)  # 动态计算高度
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=120)
        self._calculate_global_range()
        
        # 专业配色方案
        self.color_palette = {
            'cds': {'fill': '#4C72B0', 'edge': '#2F4F7F'},
            'utr': {'fill': '#E0E0E0', 'edge': '#A0A0A0'},
            'uorf': {'fill': '#59A14F', 'edge': '#36752F'},
            'dorf': {'fill': '#E15759', 'edge': '#B83A3C'},
            'intron': {'color': '#606060', 'ls': (0, (5, 2)), 'lw': 1.2}
        }
        
        # 尺寸参数
        layout = {
            'cds_height': 0.9,
            'utr_height': 0.6,
            'track_spacing': 1.8,  # 减小轨道间距
            'orf_region': 1.0
        }

        for y_idx, transcript in enumerate(self.transcripts):
            strand = transcript.get('strand', '+')
            y_center = y_idx * layout['track_spacing']
            
            # 处理5'UTR
            if transcript.get('utr5'):
                utr5 = self._adjust_coordinates([transcript['utr5']], strand)[0]
                self._draw_directional_block(ax, *utr5, y_center, layout['utr_height'], 
                                            self.color_palette['utr'], strand)
                self._draw_stranded_features(ax, transcript['uorfs'], 
                                            y_center + layout['utr_height']/2,
                                            layout['orf_region'], 'uorf', strand, 1)

            # 处理CDS
            if transcript['cds_exons']:
                sorted_exons = self._adjust_coordinates(transcript['cds_exons'], strand)
                for start, end in sorted_exons:
                    self._draw_directional_block(ax, start, end, y_center, 
                                                layout['cds_height'],
                                                self.color_palette['cds'], strand)
                
                # 绘制内含子
                for i in range(len(sorted_exons)-1):
                    self._draw_intron(ax, sorted_exons[i][1], sorted_exons[i+1][0],
                                     y_center, self.color_palette['intron'], strand)

            # 处理3'UTR
            if transcript.get('utr3'):
                utr3 = self._adjust_coordinates([transcript['utr3']], strand)[0]
                self._draw_directional_block(ax, *utr3, y_center, layout['utr_height'],
                                            self.color_palette['utr'], strand)
                self._draw_stranded_features(ax, transcript['dorfs'],
                                            y_center - layout['utr_height']/2,
                                            layout['orf_region'], 'dorf', strand, -1)

        # 坐标轴设置
        ax.set_xlim(self.min_pos - 50, self.max_pos + 50)
        ax.set_ylim(-layout['track_spacing']/2, 
                   len(self.transcripts)*layout['track_spacing'])
        ax.set_yticks([i*layout['track_spacing'] for i in range(len(self.transcripts))])
        ax.set_yticklabels([f"{t.get('name', 'Transcript')} {i+1} ({t.get('strand', '+')})" 
                          for i, t in enumerate(self.transcripts)],
                          fontsize=9)  # 减小字体
        ax.set_xlabel('Genomic Position (bp)', fontsize=10)
        
        # 添加图例并调整布局
        self._create_legend(ax)
        plt.subplots_adjust(right=0.8)  # 为图例留出空间
        plt.title('Directional Gene Structure with Adaptive Layout\n', 
                fontsize=12, pad=15)
        plt.tight_layout()
        plt.show()

# 测试数据（包含6个转录本）
test_transcripts = [
    {
        'name': 'Plus1',
        'strand': '+',
        'utr5': (100, 300),
        'cds_exons': [(320, 400), (420, 550)],
        'utr3': (560, 700),
        'uorfs': [[(120, 150)]],
        'dorfs': [[(600, 620)]]
    },
    {
        'name': 'Minus1',
        'strand': '-',
        'utr5': (650, 800),
        'cds_exons': [(570, 650), (420, 550)],
        'utr3': (100, 300),
        'uorfs': [[(750, 770)]],
        'dorfs': [[(150, 170)]]
    },
    {
        'name': 'Plus2',
        'strand': '+',
        'utr5': (200, 350),
        'cds_exons': [(370, 450), (470, 600)],
        'utr3': (620, 750),
        'uorfs': [[(220, 240)]],
        'dorfs': [[(700, 720)]]
    },
    {
        'name': 'Minus2',
        'strand': '-',
        'utr5': (700, 850),
        'cds_exons': [(620, 700), (470, 600)],
        'utr3': (200, 350),
        'uorfs': [[(800, 820)]],
        'dorfs': [[(250, 270)]]
    },
    {
        'name': 'Plus3',
        'strand': '+',
        'utr5': (300, 450),
        'cds_exons': [(470, 550), (570, 700)],
        'utr3': (720, 850),
        'uorfs': [[(320, 340)]],
        'dorfs': [[(800, 820)]]
    },
    {
        'name': 'Minus3',
        'strand': '-',
        'utr5': (750, 900),
        'cds_exons': [(670, 750), (570, 700)],
        'utr3': (300, 450),
        'uorfs': [[(850, 870),(880,890)]],
        'dorfs': [[(350, 370)]]
    },
]

plotter = GeneStructurePlotter(test_transcripts)
plotter.plot_transcripts()

