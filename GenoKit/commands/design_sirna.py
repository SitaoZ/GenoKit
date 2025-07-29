import re
from typing import Dict, List, Tuple
from collections import defaultdict

class SirnaDesigner:
    """
    siRNA设计类，支持多转录本处理
    
    参数:
        gene_name (str): 基因名称
        transcripts (Dict[str, str]): 转录本字典 {transcript_id: sequence}
        min_gc (float): 最小GC含量 (默认0.3)
        max_gc (float): 最大GC含量 (默认0.6)
        length (int): siRNA长度 (默认21)
    """
    def __init__(self, 
                 gene_name: str, 
                 transcripts: Dict[str, str],
                 min_gc: float = 0.3,
                 max_gc: float = 0.6,
                 length: int = 21):
        self.gene_name = gene_name
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.length = length
        self.candidates = defaultdict(set)
        self.sirnas = []
        self.transcripts = self._preprocess_sequences(transcripts)
    def _preprocess_sequences(self, transcripts: Dict[str, str]) -> Dict[str, str]:
        """清洗序列并验证有效性"""
        valid_seqs = {}
        for tid, seq in transcripts.items():
            clean_seq = re.sub(r'[^ATCGatcg]', '', seq).upper()
            try:
                if len(clean_seq) < self.length:
                    raise ValueError(f"转录本 {tid} 长度不足{self.length}nt")
            except Exception as e:
                print(f"\033[93mWARNING:\033[0m {e}, 该转录本将不会设计siRNA")
                continue
            valid_seqs[tid] = clean_seq
        return valid_seqs

    def generate_candidates(self) -> None:
        """生成所有可能的候选siRNA"""
        for tid, seq in self.transcripts.items():
            for i in range(len(seq) - self.length + 1):
                candidate = seq[i:i+self.length]
                self.candidates[candidate].add(tid)

    def filter_sirnas(self, 
                     require_all: bool = True, 
                     avoid_repeats: int = 3) -> None:
        """
        过滤候选siRNA
        
        参数:
            require_all (bool): 是否必须靶向所有转录本
            avoid_repeats (int): 禁止连续重复的碱基数
        """
        target_transcripts = set(self.transcripts.keys())
        
        filtered = []
        for seq, tids in self.candidates.items():
            # 检查靶向要求
            if require_all and (tids != target_transcripts):
                continue
                
            # GC含量检查
            gc = (seq.count('G') + seq.count('C')) / len(seq)
            if not (self.min_gc <= gc <= self.max_gc):
                continue
                
            # 排除连续重复
            if re.search(r'(.)\1{%d}' % (avoid_repeats-1), seq):
                continue
                
            # 优先选择AA/T开头
            priority = 1 if seq.startswith(('AA', 'T')) else 0
            
            filtered.append({
                'sequence': seq,
                'gc_content': round(gc, 2),
                'targets': len(tids),
                'start_pos': self._find_positions(seq),
                'priority': priority
            })

        # 按优先级和GC含量排序
        self.sirnas = sorted(filtered, 
                            key=lambda x: (-x['priority'], abs(x['gc_content']-0.45)))

    def _find_positions(self, seq: str) -> Dict[str, List[int]]:
        """找出序列在所有转录本中的出现位置"""
        positions = {}
        for tid, tseq in self.transcripts.items():
            positions[tid] = [m.start() for m in re.finditer(seq, tseq)]
        return positions

    def design(self, require_all: bool = True) -> List[Dict]:
        """完整设计流程"""
        self.generate_candidates()
        self.filter_sirnas(require_all)
        return self.sirnas

    def get_top_designs(self, n: int = 5) -> List[Tuple[str, float]]:
        """获取前N个最佳设计"""
        return [(s['sequence'], s['gc_content']) for s in self.sirnas[:n]]


if __name__ == "__main__":
    # 使用示例
    test_transcripts = {
        "NM_001": "ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCGCGATCGATCGATCGATCGA",
        "NM_002": "ATCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
    }

    designer = SirnaDesigner("TP53", test_transcripts, length = 22)
    designer = SirnaDesigner("TP53", test_transcripts, length = 22)
    
    # 设计靶向所有转录本的siRNA
    universal = designer.design(require_all=True)
    print(f"找到 {len(universal)} 条通用siRNA:")
    print(designer.get_top_designs(3))
    
    # 设计任意靶向siRNA
    designer.filter_sirnas(require_all=False)
    print(f"\n找到 {len(designer.sirnas)} 条任意靶向siRNA:")
    print(designer.get_top_designs(3))

