#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: extract_ppi.py
date: 2025/12/10 09:53
author: Sitao Zhu
mail: zhusitao1990@163.com
'''


import logging
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Any, Optional



_PACKAGE_ROOT = Path(__file__).resolve().parents[1]
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def _resolve_input_path(value):
    """Resolve a user path, falling back to GenoKit's packaged examples."""
    requested = Path(value).expanduser()
    if requested.is_file():
        return requested
    if not requested.is_absolute() and requested.parts and requested.parts[0] == "example":
        packaged = _PACKAGE_ROOT.joinpath(*requested.parts)
        if packaged.is_file():
            return packaged
    raise FileNotFoundError(
        "PPI input file not found: %s. Relative example/... paths are also "
        "searched under %s." % (value, _PACKAGE_ROOT / "example")
    )


class PPINetworkAnalyzer:
    """
    蛋白质-蛋白质相互作用网络分析类

    功能：
    1. 构建和操作PPI网络
    2. 计算网络拓扑特征
    3. 可视化网络
    4. 分析和导出网络数据
    """

    def __init__(self, name: str = "PPI_Network"):
        """
        初始化PPI网络分析器

        参数：
        name: 网络名称
        """
        self.name = name
        self.network = nx.Graph(name=name)
        self.proteins = []
        self.centrality_measures = {}
        self.last_layout = None

    def _invalidate_analysis(self) -> None:
        """网络结构变化后清除派生结果，避免返回过期的中心性数值。"""
        self.centrality_measures.clear()

    def _sync_proteins(self) -> None:
        """让兼容旧代码的 proteins 列表始终与 NetworkX 节点保持一致。"""
        self.proteins = list(self.network.nodes())

    def add_proteins(self, proteins: List[str]) -> None:
        """
        添加蛋白质节点到网络

        参数：
        proteins: 蛋白质名称列表
        """
        proteins = list(proteins)
        if not proteins:
            return
        self.network.add_nodes_from(proteins)
        self._sync_proteins()
        self._invalidate_analysis()

    def add_interactions(self, interactions: List[Tuple[str, str]],
                         weights: Optional[List[float]] = None) -> None:
        """
        添加相互作用到网络

        参数：
        interactions: 相互作用的元组列表 [(protein1, protein2), ...]
        weights: 可选的权重列表
        """
        interactions = list(interactions)
        if not interactions:
            return

        if weights is None:
            self.network.add_edges_from(interactions)
        else:
            weights = list(weights)
            # 兼容原有调用方式：无权网络会传入空列表。
            if not weights:
                self.network.add_edges_from(interactions)
            elif len(weights) != len(interactions):
                raise ValueError(
                    "weights 的数量必须与 interactions 相同: "
                    f"{len(weights)} != {len(interactions)}"
                )
            else:
                weighted_edges = []
                for (protein1, protein2), weight in zip(interactions, weights):
                    try:
                        numeric_weight = float(weight)
                    except (TypeError, ValueError) as exc:
                        raise ValueError(f"无效的互作权重: {weight!r}") from exc
                    weighted_edges.append((protein1, protein2, numeric_weight))
                self.network.add_weighted_edges_from(weighted_edges)

        self._sync_proteins()
        self._invalidate_analysis()

    def remove_protein(self, protein: str) -> None:
        """
        从网络中移除蛋白质节点

        参数：
        protein: 要移除的蛋白质名称
        """
        if protein in self.network:
            self.network.remove_node(protein)
            self._sync_proteins()
            self._invalidate_analysis()

    def remove_interaction(self, protein1: str, protein2: str) -> None:
        """
        从网络中移除相互作用

        参数：
        protein1, protein2: 相互作用的两个蛋白质
        """
        if self.network.has_edge(protein1, protein2):
            self.network.remove_edge(protein1, protein2)
            self._invalidate_analysis()

    def get_network_statistics(self) -> Dict[str, Any]:
        """
        获取网络基本统计信息

        返回：
        包含网络统计信息的字典
        """
        num_nodes = self.network.number_of_nodes()
        is_connected = nx.is_connected(self.network) if num_nodes else False
        stats = {
            'name': self.name,
            'num_nodes': num_nodes,
            'num_edges': self.network.number_of_edges(),
            'density': nx.density(self.network),
            'is_connected': is_connected,
            'num_connected_components': nx.number_connected_components(self.network) if num_nodes else 0,
            'average_clustering': nx.average_clustering(self.network) if num_nodes else 0.0,
            'average_shortest_path_length': None,
            'diameter': None
        }

        # 只有在连通图中才计算这些指标
        if is_connected:
            stats['average_shortest_path_length'] = nx.average_shortest_path_length(self.network)
            stats['diameter'] = nx.diameter(self.network)

        return stats

    def calculate_centrality_measures(self, include: List[str] = None) -> Dict[str, Dict[str, float]]:
        """
        计算各种中心性指标

        参数：
        include: 要计算的中心性指标列表，可选 ['degree', 'betweenness', 'closeness', 'eigenvector']

        返回：
        包含各种中心性指标的字典
        """
        if include is None:
            include = ['degree', 'betweenness', 'closeness', 'eigenvector']

        include = list(dict.fromkeys(include))
        supported = {'degree', 'betweenness', 'closeness', 'eigenvector', 'pagerank'}
        unsupported = [metric for metric in include if metric not in supported]
        if unsupported:
            raise ValueError(f"不支持的中心性指标: {', '.join(unsupported)}")

        # 每种指标独立缓存。plot_centrality_comparison 再次调用时不会重复
        # 执行代价较高的 betweenness/closeness 计算。
        missing = [metric for metric in include if metric not in self.centrality_measures]
        if not self.network:
            for metric in missing:
                self.centrality_measures[metric] = {}
        else:
            if 'degree' in missing:
                self.centrality_measures['degree'] = nx.degree_centrality(self.network)
            if 'betweenness' in missing:
                self.centrality_measures['betweenness'] = nx.betweenness_centrality(self.network)
            if 'closeness' in missing:
                self.centrality_measures['closeness'] = nx.closeness_centrality(self.network)
            if 'eigenvector' in missing:
                try:
                    values = nx.eigenvector_centrality(self.network, max_iter=1000)
                except nx.PowerIterationFailedConvergence:
                    logging.warning("特征向量中心性未收敛，结果使用 0 填充。")
                    values = {node: 0.0 for node in self.network}
                self.centrality_measures['eigenvector'] = values
            if 'pagerank' in missing:
                self.centrality_measures['pagerank'] = nx.pagerank(self.network)

        return {metric: self.centrality_measures[metric] for metric in include}

    def get_top_hubs(self, n: int = 10, centrality_type: str = 'degree') -> List[Tuple[str, float]]:
        """
        获取最重要的中心节点

        参数：
        n: 返回的节点数量
        centrality_type: 使用的中心性类型

        返回：
        排序后的(节点, 中心性值)列表
        """
        if centrality_type not in self.centrality_measures:
            self.calculate_centrality_measures([centrality_type])

        centrality = self.centrality_measures[centrality_type]
        if n < 0:
            raise ValueError("n 必须大于或等于 0")
        return sorted(centrality.items(), key=lambda item: (-item[1], str(item[0])))[:n]

    def detect_communities(self, method: str = 'louvain') -> Dict[int, List[str]]:
        """
        检测网络中的社区结构

        参数：
        method: 社区检测方法，支持 'louvain', 'girvan_newman', 'label_propagation'

        返回：
        社区ID到节点列表的映射字典
        """
        if not self.network:
            return {}

        if method == 'louvain':
            louvain = getattr(nx.algorithms.community, 'louvain_communities', None)
            if louvain is not None:
                groups = louvain(self.network, weight='weight', seed=42)
            else:
                # 兼容较旧的 NetworkX，无需强制安装 python-louvain。
                logging.warning("当前 NetworkX 不支持 Louvain，改用 greedy modularity。")
                groups = nx.algorithms.community.greedy_modularity_communities(
                    self.network, weight='weight'
                )
        elif method == 'girvan_newman':
            if self.network.number_of_edges() == 0:
                groups = [{node} for node in self.network]
            else:
                groups = next(nx.algorithms.community.girvan_newman(self.network))
        elif method == 'label_propagation':
            groups = nx.algorithms.community.asyn_lpa_communities(
                self.network, weight='weight', seed=42
            )
        else:
            raise ValueError(f"不支持的社区检测方法: {method}")

        normalized = [sorted(group, key=str) for group in groups]
        normalized.sort(key=lambda group: (-len(group), str(group[0]) if group else ''))
        return {index: group for index, group in enumerate(normalized)}

    def get_degree_node_sizes(self,
                              base_size: float = 1000,
                              min_size: Optional[float] = None,
                              max_size: Optional[float] = None) -> Dict[Any, float]:
        """根据节点 degree 线性计算绘图大小。

        ``base_size`` 是所有节点 degree 相同时使用的大小，也作为自动缩放
        范围的基准。默认情况下，最低和最高 degree 分别映射到
        ``0.5 * base_size`` 和 ``2.0 * base_size``。

        返回值以节点为键，可确保节点大小与 NetworkX 的绘图顺序正确对应。
        """
        if base_size <= 0:
            raise ValueError("base_size 必须大于 0")

        min_size = base_size * 0.5 if min_size is None else min_size
        max_size = base_size * 2.0 if max_size is None else max_size
        if min_size <= 0 or max_size <= 0:
            raise ValueError("min_size 和 max_size 必须大于 0")
        if min_size > max_size:
            raise ValueError("min_size 不能大于 max_size")

        degrees = dict(self.network.degree())
        if not degrees:
            return {}

        min_degree = min(degrees.values())
        max_degree = max(degrees.values())
        if min_degree == max_degree:
            return {node: float(base_size) for node in degrees}

        degree_span = max_degree - min_degree
        size_span = max_size - min_size
        return {
            node: float(min_size + (degree - min_degree) / degree_span * size_span)
            for node, degree in degrees.items()
        }

    def _tree_layout(self) -> Dict[Any, np.ndarray]:
        """生成类似 igraph ``layout_as_tree`` 的分层布局。

        每个连通分量选择 degree 最大的节点作为根；不同分量沿 x 轴分开，
        因此森林或包含孤立节点的网络也能正常绘制。
        """
        positions = {}
        x_offset = 0.0
        components = sorted(
            nx.connected_components(self.network),
            key=lambda nodes: (-len(nodes), min(map(str, nodes)))
        )
        for component in components:
            subgraph = self.network.subgraph(component)
            root = min(
                component,
                key=lambda node: (-subgraph.degree(node), str(node))
            )
            distances = nx.single_source_shortest_path_length(subgraph, root)
            layers = {}
            for node, depth in distances.items():
                layers.setdefault(depth, []).append(node)

            component_width = max(len(nodes) for nodes in layers.values())
            for depth, nodes in sorted(layers.items()):
                nodes = sorted(nodes, key=str)
                start_x = x_offset + (component_width - len(nodes)) / 2.0
                for index, node in enumerate(nodes):
                    positions[node] = np.array([start_x + index, -float(depth)])
            x_offset += component_width + 1.0

        return nx.rescale_layout_dict(positions) if positions else {}

    def _grid_layout(self) -> Dict[Any, np.ndarray]:
        """将节点按接近正方形的规则网格排列。"""
        nodes = sorted(self.network.nodes(), key=str)
        if not nodes:
            return {}
        columns = int(np.ceil(np.sqrt(len(nodes))))
        positions = {
            node: np.array([float(index % columns), -float(index // columns)])
            for index, node in enumerate(nodes)
        }
        return nx.rescale_layout_dict(positions)

    def _shell_layout(self) -> Dict[Any, np.ndarray]:
        """按 degree 从高到低将节点放入同心层。"""
        degree_groups = {}
        for node, degree in self.network.degree():
            degree_groups.setdefault(degree, []).append(node)
        shells = [
            sorted(degree_groups[degree], key=str)
            for degree in sorted(degree_groups, reverse=True)
        ]
        return nx.shell_layout(self.network, nlist=shells) if shells else {}

    def get_layout(self, layout: str = 'auto', seed: int = 42) -> Dict[Any, np.ndarray]:
        """计算网络节点坐标，布局名称兼容常见 igraph 表达。

        ``auto`` 对森林采用 tree，小型一般网络采用 Kamada-Kawai，
        中型网络采用 Fruchterman-Reingold，大型网络采用 ForceAtlas2。
        所有含随机过程的布局均使用固定 seed，保证重复运行结果一致。
        """
        aliases = {
            'spring': 'fr',
            'fruchterman_reingold': 'fr',
            'kamada_kawai': 'kk',
            'circular': 'circle',
            'as_star': 'star',
            'reingold_tilford': 'tree',
            'drl': 'forceatlas2'
        }
        layout = aliases.get(layout.lower(), layout.lower())
        supported = {
            'auto', 'fr', 'kk', 'forceatlas2', 'circle', 'star', 'tree',
            'grid', 'shell', 'spectral', 'spiral', 'random'
        }
        if layout not in supported:
            raise ValueError(
                f"不支持的网络布局: {layout}；可用布局: {', '.join(sorted(supported))}"
            )
        if not self.network:
            self.last_layout = layout
            return {}

        if layout == 'auto':
            node_count = self.network.number_of_nodes()
            if nx.is_forest(self.network):
                layout = 'tree'
            elif node_count <= 40:
                layout = 'kk'
            elif node_count <= 500:
                layout = 'fr'
            else:
                layout = 'forceatlas2'

        if layout == 'fr':
            positions = nx.spring_layout(
                self.network, k=None, iterations=200, weight='weight', seed=seed
            )
        elif layout == 'kk':
            # PPI 权重通常代表置信度而不是几何距离，因此 KK 布局不把权重
            # 当作最短路径距离使用。
            positions = nx.kamada_kawai_layout(self.network, weight=None)
        elif layout == 'forceatlas2':
            forceatlas2 = getattr(nx, 'forceatlas2_layout', None)
            if forceatlas2 is None:
                logging.warning("当前 NetworkX 不支持 ForceAtlas2，改用 FR 布局。")
                layout = 'fr'
                positions = nx.spring_layout(
                    self.network, k=None, iterations=200, weight='weight', seed=seed
                )
            else:
                positions = forceatlas2(
                    self.network, max_iter=200, scaling_ratio=2.0,
                    gravity=1.0, weight='weight', seed=seed
                )
        elif layout == 'circle':
            positions = nx.circular_layout(self.network)
        elif layout == 'star':
            hub = min(
                self.network.nodes(),
                key=lambda node: (-self.network.degree(node), str(node))
            )
            outer_nodes = sorted(
                (node for node in self.network if node != hub), key=str
            )
            positions = nx.shell_layout(self.network, nlist=[[hub], outer_nodes])
        elif layout == 'tree':
            positions = self._tree_layout()
        elif layout == 'grid':
            positions = self._grid_layout()
        elif layout == 'shell':
            positions = self._shell_layout()
        elif layout == 'spectral':
            positions = nx.spectral_layout(self.network, weight='weight')
        elif layout == 'spiral':
            positions = nx.spiral_layout(self.network, equidistant=True)
        else:
            positions = nx.random_layout(self.network, seed=seed)

        self.last_layout = layout
        logging.info("PPI network layout: %s", layout)
        return positions

    def plot_network(self,
                     layout: str = 'auto',
                     node_size: int = 1000,
                     node_color: str = 'lightblue',
                     with_labels: bool = True,
                     figsize: Tuple[int, int] = (10, 8),
                     title: Optional[str] = None,
                     save_path: Optional[str] = None,
                     scale_node_by_degree: bool = True,
                     min_node_size: Optional[float] = None,
                     max_node_size: Optional[float] = None,
                     show_degree_legend: bool = True) -> None:
        """
        可视化网络

        参数：
        layout: 布局算法，支持 auto/fr/kk/forceatlas2/circle/star/tree/
                grid/shell/spectral/spiral/random，并兼容原有布局名称
        node_size: 节点基准大小
        node_color: 节点颜色
        scale_node_by_degree: 是否根据 degree 动态缩放节点
        min_node_size: 最低 degree 对应的节点大小，默认 node_size * 0.5
        max_node_size: 最高 degree 对应的节点大小，默认 node_size * 2.0
        show_degree_legend: 是否显示节点大小与 degree 的对应图例
        with_labels: 是否显示标签
        figsize: 图形大小
        title: 图形标题
        save_path: 保存路径
        """
        pos = self.get_layout(layout)

        fig, ax = plt.subplots(figsize=figsize)

        if scale_node_by_degree:
            degree_sizes = self.get_degree_node_sizes(
                base_size=node_size,
                min_size=min_node_size,
                max_size=max_node_size
            )
            plotted_node_sizes = [degree_sizes[node] for node in self.network.nodes()]
        else:
            degree_sizes = {node: float(node_size) for node in self.network.nodes()}
            plotted_node_sizes = node_size

        # 绘制网络
        nx.draw_networkx_nodes(self.network, pos,
                               node_size=plotted_node_sizes,
                               node_color=node_color,
                               alpha=0.8,
                               ax=ax)

        nx.draw_networkx_edges(self.network, pos,
                               width=1.0,
                               alpha=0.5,
                               edge_color='gray',
                               ax=ax)

        if with_labels:
            nx.draw_networkx_labels(self.network, pos,
                                    font_size=10,
                                    font_weight='bold',
                                    ax=ax)

        if scale_node_by_degree and show_degree_legend and self.network:
            unique_degrees = sorted(set(dict(self.network.degree()).values()))
            if len(unique_degrees) > 3:
                unique_degrees = [
                    unique_degrees[0],
                    unique_degrees[len(unique_degrees) // 2],
                    unique_degrees[-1]
                ]
            degree_to_size = {
                degree: degree_sizes[next(
                    node for node, node_degree in self.network.degree()
                    if node_degree == degree
                )]
                for degree in unique_degrees
            }
            # NetworkX/Matplotlib 的 node_size 是面积（points^2）。图例使用
            # 较小的直径示意，既保留相对大小，又不会遮挡实际网络节点。
            legend_handles = [
                plt.Line2D(
                    [], [], linestyle='', marker='o',
                    markersize=max(5.0, np.sqrt(size) * 0.35),
                    markerfacecolor=node_color, markeredgecolor=node_color,
                    alpha=0.8, label=f'Degree {degree}'
                )
                for degree, size in degree_to_size.items()
            ]
            ax.legend(handles=legend_handles, title='Node size',
                      loc='upper left', bbox_to_anchor=(1.01, 1.0),
                      borderaxespad=0.0, frameon=True, framealpha=0.9)

        if title is None:
            title = f"{self.name} - PPI Network"

        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.margins(0.15)
        ax.axis('off')
        if scale_node_by_degree and show_degree_legend and self.network:
            fig.tight_layout(rect=(0, 0, 0.84, 1))
        else:
            fig.tight_layout()

        if save_path:
            save_path = Path(save_path).expanduser()
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(str(save_path), dpi=300, bbox_inches='tight')
            logging.info(f"图形已保存至: {save_path}")
        else:
            plt.show()
        plt.close(fig)

    def plot_centrality_comparison(self,
                                   centrality_types: List[str] = None,
                                   top_n: int = 15,
                                   figsize: Tuple[int, int] = (12, 6),
                                   save_path: Optional[str] = None) -> None:
        """
        绘制不同中心性指标的对比图

        参数：
        centrality_types: 要对比的中心性类型
        top_n: 显示前N个节点
        figsize: 图形大小
        """
        if centrality_types is None:
            centrality_types = ['degree', 'betweenness', 'closeness']

        if top_n < 1:
            raise ValueError("top_n 必须大于 0")
        self.calculate_centrality_measures(centrality_types)

        # 获取所有中心性数据
        centrality_data = {}
        for c_type in centrality_types:
            centrality = self.centrality_measures[c_type]
            sorted_items = sorted(
                centrality.items(), key=lambda item: (-item[1], str(item[0]))
            )[:top_n]
            nodes = [item[0] for item in sorted_items]
            values = [item[1] for item in sorted_items]
            centrality_data[c_type] = {'nodes': nodes, 'values': values}

        # 绘制子图
        fig, axes = plt.subplots(1, len(centrality_types), figsize=figsize)

        if len(centrality_types) == 1:
            axes = [axes]

        for c_type, ax in zip(centrality_types, axes):
            data = centrality_data[c_type]
            nodes = data['nodes']
            values = data['values']

            colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(nodes)))
            ax.barh(nodes, values, color=colors)
            ax.set_xlabel(f'{c_type.capitalize()} Centrality', fontweight='bold')
            ax.set_title(f'Top {top_n} Nodes by {c_type.capitalize()} Centrality')
            ax.invert_yaxis()

        fig.suptitle('Centrality Measures Comparison', fontsize=14, fontweight='bold')
        fig.tight_layout()
        if save_path:
            save_path = Path(save_path).expanduser()
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(str(save_path), dpi=300, bbox_inches='tight')
            logging.info(f"图形已保存至: {save_path}")
        else:
            plt.show()
        plt.close(fig)

    def export_to_dataframe(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        将网络数据导出为DataFrame

        返回：
        包含网络信息的DataFrame
        """
        required_metrics = ['degree', 'betweenness', 'closeness']
        self.calculate_centrality_measures(required_metrics)

        # 节点数据
        nodes_data = []
        for node in sorted(self.network.nodes(), key=str):
            degree = self.network.degree(node)
            nodes_data.append({
                'protein': node,
                'degree': degree,
                'degree_centrality': self.centrality_measures.get('degree', {}).get(node, 0),
                'betweenness_centrality': self.centrality_measures.get('betweenness', {}).get(node, 0),
                'closeness_centrality': self.centrality_measures.get('closeness', {}).get(node, 0)
            })

        # 边数据
        edges_data = []
        edges = sorted(
            self.network.edges(data=True),
            key=lambda edge: (str(edge[0]), str(edge[1]))
        )
        for edge in edges:
            edges_data.append({
                'protein1': edge[0],
                'protein2': edge[1],
                'weight': edge[2].get('weight', 1.0)
            })

        nodes_df = pd.DataFrame(nodes_data)
        edges_df = pd.DataFrame(edges_data)

        return nodes_df, edges_df

    def save_network(self, filepath: str, format: str = 'gexf') -> None:
        """
        保存网络到文件

        参数：
        filepath: 文件路径
        format: 文件格式 ('gexf', 'graphml', 'gml')
        """
        writers = {
            'gexf': nx.write_gexf,
            'graphml': nx.write_graphml,
            'gml': nx.write_gml
        }
        format = format.lower()
        if format not in writers:
            raise ValueError(f"不支持的文件格式: {format}")

        filepath = Path(filepath).expanduser()
        filepath.parent.mkdir(parents=True, exist_ok=True)
        writers[format](self.network, str(filepath))

        # print(f"网络已保存至: {filepath}")
        logging.info(f"网络已保存至: {filepath}")

    def load_network(self, filepath: str, format: str = 'gexf') -> None:
        """
        从文件加载网络

        参数：
        filepath: 文件路径
        format: 文件格式 ('gexf', 'graphml', 'gml')
        """
        readers = {
            'gexf': nx.read_gexf,
            'graphml': nx.read_graphml,
            'gml': nx.read_gml
        }
        format = format.lower()
        if format not in readers:
            raise ValueError(f"不支持的文件格式: {format}")

        filepath = Path(filepath).expanduser()
        self.network = readers[format](str(filepath))
        self.name = self.network.graph.get('name', self.name)
        self._sync_proteins()
        self._invalidate_analysis()
        logging.info(f"网络已从 {filepath} 加载，包含 {len(self.proteins)} 个节点")


def get_ppi(args):
    '''

    :param args:
    :return:
    '''
    input_path = _resolve_input_path(args.input)
    output_dir = Path(args.output_dir).expanduser() if args.output_dir else None
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        log_path = output_dir / "PPI.log"
    else:
        log_path = Path("PPI.log")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename=str(log_path)
    )

    # 互作文件分隔符
    SEP = {'tab': '\t', 'comma': ','}
    separator = SEP[args.separator]
    logging.info(f"Perform PPI interaction network")
    # 创建PPI网络分析器实例
    ppi_analyzer = PPINetworkAnalyzer(name="Interaction_Network")

    # 添加蛋白质节点
    proteins = set()
    interactions = []
    weights = []
    logging.info(f"Processing PPI interaction file")
    with input_path.open(encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            array = line.split(separator)
            if len(array) != 2:
                print("PPI file contains missing values.")
            proteinA, proteinB = array[0], array[1]
            proteins.add(proteinA)
            proteins.add(proteinB)
            interactions.append((proteinA, proteinB))
            if len(array) == 3:
                weight = array[2]
                weights.append(weight)
    proteins = list(proteins)
    ppi_analyzer.add_proteins(proteins)

    # 添加相互作用（可带权重）
    ppi_analyzer.add_interactions(interactions, weights)

    # 获取网络统计信息
    stats = ppi_analyzer.get_network_statistics()
    logging.info("网络统计信息:")

    for key, value in stats.items():
        # print(f"  {key}: {value}")
        logging.info(f"  {key}: {value}")

    # 计算中心性指标
    centrality = ppi_analyzer.calculate_centrality_measures()

    # print("\n中心性指标:")
    logging.info("\n中心性指标:")
    for protein in proteins[:5]:  # 只显示前5个蛋白质
        # print(f"{protein}:")
        logging.info(f"{protein}")
        for c_type, values in centrality.items():
            # print(f"  {c_type}: {values.get(protein, 0):.3f}")
            logging.info(f"  {c_type}: {values.get(protein, 0):.3f}")

    # 获取最重要的中心节点
    top_hubs = ppi_analyzer.get_top_hubs(n=5, centrality_type='degree')
    # print(f"\nTop 5 Hub Proteins (by degree centrality):")
    logging.info(f"\nTop 5 Hub Proteins (by degree centrality):")
    for protein, centrality_value in top_hubs:
        # print(f"  {protein}: {centrality_value:.3f}")
        logging.info(f"  {protein}: {centrality_value:.3f}")

    # 检测社区
    # pip install python-louvain 0.16
    communities = ppi_analyzer.detect_communities(method='louvain')
    # print(f"\n检测到 {len(communities)} 个社区:")
    logging.info(f"\n检测到 {len(communities)} 个社区:")
    for comm_id, nodes in communities.items():
        # print(f"  社区 {comm_id}: {', '.join(nodes)}")
        logging.info(f"  社区 {comm_id}: {', '.join(nodes)}")

    # 可视化网络
    if output_dir:
        network_pdf_path = str(output_dir / "ppi_network.pdf")
    else:
        network_pdf_path = None
    ppi_analyzer.plot_network(
        layout=getattr(args, 'layout', 'auto'),
        node_size=1500,
        node_color='skyblue',
        figsize=(10, 8),
        title="Protein-Protein Interaction Network",
        save_path=network_pdf_path # 可以设置为 'ppi_network.png' 来保存
    )

    # 对比不同中心性指标
    if output_dir:
        centrality_pdf_path = str(output_dir / "ppi_centrality.pdf")
    else:
        centrality_pdf_path = None
    ppi_analyzer.plot_centrality_comparison(
        centrality_types=['degree', 'betweenness', 'closeness'],
        top_n=8,
        figsize=(15, 5),
        save_path=centrality_pdf_path
    )
    # 导出数据
    if output_dir:
        centrality_path = str(output_dir / "ppi_centrality.csv")
    else:
        centrality_path = None
    nodes_df, edges_df = ppi_analyzer.export_to_dataframe()
    nodes_df.to_csv(centrality_path, sep='\t', index=False)
    # print(f"\n节点数据 (前5行):")
    logging.info(f"\n节点数据 (前5行):")
    # print(nodes_df.head())
    logging.info(nodes_df.head())
    # print(f"\n边数据 (前5行):")
    logging.info(f"\n边数据 (前5行):")
    # print(edges_df.head())
    logging.info(edges_df.head())

    # 保存网络
    if output_dir:
        network_gexf_path = str(output_dir / "ppi_network.gexf")
    else:
        network_gexf_path = None
    ppi_analyzer.save_network(network_gexf_path, format='gexf')



if __name__ == "__main__":
    a = 0
    get_ppi(a)
