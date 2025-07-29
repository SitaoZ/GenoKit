import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

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
            # 收集所有特征位置
            if transcript['utr5']:
                features.extend(transcript['utr5'])
            if transcript['utr3']:
                features.extend(transcript['utr3'])
            for exon in transcript['cds_exons']:
                features.extend(exon)
            for uorf in transcript['uorfs']:
                features.extend(uorf)
            for dorf in transcript['dorfs']:
                features.extend(dorf)
            
            if features:
                current_min = min(features)
                current_max = max(features)
                self.min_pos = min(self.min_pos, current_min) if self.min_pos else current_min
                self.max_pos = max(self.max_pos, current_max) if self.max_pos else current_max

    def _draw_stranded_features(self, ax, features, y_base, region_height, orf_type):
        """绘制带方向的ORF特征"""
        max_layers = 3  # 最大分层数
        layer_height = region_height * 0.8 / max_layers
        color = self.color_palette['uorf'] if orf_type == 'uorf' else self.color_palette['dorf']

        for idx, (start, end) in enumerate(features):
            # 计算分层位置
            layer = idx % max_layers
            vertical_offset = y_base + (layer_height * layer * (-1 if orf_type == 'dorf' else 1))
            
            # 设置透明度渐变
            alpha = 0.7 + 0.3 * (1 - layer/max_layers)
            
            # 绘制ORF框
            ax.add_patch(Rectangle(
                (start, vertical_offset),
                end - start,
                layer_height*0.8,
                facecolor=color['fill'],
                edgecolor=color['edge'],
                alpha=alpha,
                lw=0.6
            ))
            
            # 添加编号标注（长度>15bp才显示）
            if (end - start) > 15:
                ax.text((start+end)/2, vertical_offset + layer_height*0.4,
                       str(idx+1), 
                       ha='center', 
                       va='center',
                       fontsize=7,
                       color=color['edge'],
                       weight='bold')

    def _draw_intron(self, ax, start, end, y, style):
        """绘制内含子连接线"""
        if start < end:
            ax.hlines(y, start, end, 
                     colors=style['color'],
                     linestyles=style['ls'],
                     linewidth=style['lw'])

    def _create_legend(self, ax):
        """创建专业图例"""
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
                      label='Intron')
        ]
        ax.legend(handles=legend_elements, 
                loc='upper right',
                frameon=True,
                framealpha=0.95,
                edgecolor='#FFFFFF',
                fontsize=9,
                title='Legend:',
                title_fontsize=10)

    def plot_transcripts(self):
        """主绘图方法"""
        plt.style.use('seaborn-whitegrid')  # 使用专业样式
        fig, ax = plt.subplots(figsize=(12, 6), dpi=120)
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
            'cds_height': 0.8,
            'utr_height': 0.5,
            'track_spacing': 1.8,
            'orf_region': 0.6
        }

        for y_idx, transcript in enumerate(self.transcripts):
            y_center = y_idx * layout['track_spacing']
            
            # ===== 绘制5'UTR区域 =====
            if transcript['utr5']:
                utr5_start, utr5_end = transcript['utr5']
                # UTR主体
                ax.add_patch(Rectangle(
                    (utr5_start, y_center - layout['utr_height']/2),
                    utr5_end - utr5_start,
                    layout['utr_height'],
                    facecolor=self.color_palette['utr']['fill'],
                    edgecolor=self.color_palette['utr']['edge'],
                    lw=0.8
                ))
                # 绘制uORFs（上方分层）
                self._draw_stranded_features(
                    ax=ax,
                    features=transcript['uorfs'],
                    y_base=y_center + layout['utr_height']/2,
                    region_height=layout['orf_region'],
                    orf_type='uorf'
                )
                
            # ===== 绘制CDS区域 =====
            if transcript['cds_exons']:
                # 连接UTR到CDS
                if transcript.get('utr5'):
                    self._draw_intron(
                        ax=ax,
                        start=transcript['utr5'][1],
                        end=transcript['cds_exons'][0][0],
                        y=y_center,
                        style=self.color_palette['intron']
                    )
                
                # 绘制CDS外显子
                for i, (start, end) in enumerate(transcript['cds_exons']):
                    ax.add_patch(Rectangle(
                        (start, y_center - layout['cds_height']/2),
                        end - start,
                        layout['cds_height'],
                        facecolor=self.color_palette['cds']['fill'],
                        edgecolor=self.color_palette['cds']['edge'],
                        lw=1.2
                    ))
                    
                    # 绘制内含子
                    if i < len(transcript['cds_exons']) - 1:
                        self._draw_intron(
                            ax=ax,
                            start=end,
                            end=transcript['cds_exons'][i+1][0],
                            y=y_center,
                            style=self.color_palette['intron']
                        )
                
                # 连接CDS到3'UTR
                if transcript.get('utr3'):
                    self._draw_intron(
                        ax=ax,
                        start=transcript['cds_exons'][-1][1],
                        end=transcript['utr3'][0],
                        y=y_center,
                        style=self.color_palette['intron']
                    )
            
            # ===== 绘制3'UTR区域 =====
            if transcript['utr3']:
                utr3_start, utr3_end = transcript['utr3']
                ax.add_patch(Rectangle(
                    (utr3_start, y_center - layout['utr_height']/2),
                    utr3_end - utr3_start,
                    layout['utr_height'],
                    facecolor=self.color_palette['utr']['fill'],
                    edgecolor=self.color_palette['utr']['edge'],
                    lw=0.8
                ))
                # 绘制dORFs（下方分层）
                self._draw_stranded_features(
                    ax=ax,
                    features=transcript['dorfs'],
                    y_base=y_center - layout['utr_height']/2,
                    region_height=layout['orf_region'],
                    orf_type='dorf'
                )

        # ===== 坐标轴设置 =====
        ax.set_xlim(self.min_pos - 50, self.max_pos + 50)
        ax.set_ylim(-layout['track_spacing']/2, 
                   len(self.transcripts)*layout['track_spacing'])
        ax.set_yticks([i*layout['track_spacing'] for i in range(len(self.transcripts))])
        ax.set_yticklabels([f'Isoform {i+1}' for i in range(len(self.transcripts))],
                          fontsize=10)
        ax.set_xlabel('Genomic Position (bp)', fontsize=11, labelpad=10)
        ax.tick_params(axis='x', which='major', labelsize=9)
        
        # ===== 添加图例和标题 =====
        self._create_legend(ax)
        plt.title('Gene Structure with Alternative ORFs\n', 
                fontsize=13, pad=20, weight='semibold')
        
        plt.tight_layout()
        plt.show()

# 测试数据（包含密集ORF）
test_transcripts = [
    {
        'utr5': (100, 300),
        'cds_exons': [(320, 400), (420, 550), (570, 650)],
        'utr3': (670, 800),
        'uorfs': [
            (110, 130), (125, 150), 
            (140, 180), (170, 200),
            (210, 240), (230, 260)
        ],
        'dorfs': [
            (680, 700), (690, 720),
            (710, 740), (730, 760),
            (750, 780), (770, 790)
        ]
    },
    {
        'utr5': (150, 280),
        'cds_exons': [(300, 450), (480, 600)],
        'utr3': (620, 750),
        'uorfs': [
            (160, 190), (180, 210),
            (200, 230), (220, 250)
        ],
        'dorfs': [
            (630, 660), (650, 680),
            (670, 700), (690, 720)
        ]
    }
]

# 执行绘图
plotter = GeneStructurePlotter(test_transcripts)
plotter.plot_transcripts()
