class UTRCalculator:
    def __init__(self, exons_list, cds_list, strand):
        """
        初始化UTR计算器
        
        参数:
        exons_list: 外显子位置列表，格式为[(start1, end1), (start2, end2), ...]
        cds_list: CDS位置列表，格式同上
        strand: str +/-
        """
        self.exons = self._merge_intervals(sorted(exons_list))
        self.cds = self._merge_intervals(sorted(cds_list))
        self.strand = strand
        
    def _merge_intervals(self, intervals):
        """合并重叠或相邻的区间"""
        if not intervals:
            return []
            
        merged = [intervals[0]]
        for current in intervals[1:]:
            last = merged[-1]
            if current[0] <= last[1] + 1:  # 重叠或相邻
                merged[-1] = (last[0], max(last[1], current[1]))
            else:
                merged.append(current)
        return merged
    
    def calculate_utr5(self):
        """计算5'UTR区域"""
        if not self.cds:
            # return self.exons
            return []
            
        if self.strand == '+':
            # 正链: 5'UTR在第一个CDS之前
            first_cds_start = self.cds[0][0]
            utr5 = []
            for exon in self.exons:
                if exon[1] < first_cds_start:
                    utr5.append(exon)
                elif exon[0] < first_cds_start:
                    utr5.append((exon[0], first_cds_start - 1))
                    break
                else:
                    break
            return utr5
        else:
            # 负链: 5'UTR在最后一个CDS之后
            last_cds_end = self.cds[-1][1]
            utr5 = []
            for exon in reversed(self.exons):
                if exon[0] > last_cds_end:
                    utr5.insert(0, exon)
                elif exon[1] > last_cds_end:
                    utr5.insert(0, (last_cds_end + 1, exon[1]))
                    break
                else:
                    break
            return utr5
    
    def calculate_utr3(self):
        """计算3'UTR区域"""
        if not self.cds:
            return []
            
        if self.strand == '+':
            # 正链: 3'UTR在最后一个CDS之后
            last_cds_end = self.cds[-1][1]
            utr3 = []
            for exon in reversed(self.exons):
                if exon[0] > last_cds_end:
                    utr3.insert(0, exon)
                elif exon[1] > last_cds_end:
                    utr3.insert(0, (last_cds_end + 1, exon[1]))
                    break
                else:
                    break
            return utr3
        else:
            # 负链: 3'UTR在第一个CDS之前
            first_cds_start = self.cds[0][0]
            utr3 = []
            for exon in self.exons:
                if exon[1] < first_cds_start:
                    utr3.append(exon)
                elif exon[0] < first_cds_start:
                    utr3.append((exon[0], first_cds_start - 1))
                    break
                else:
                    break
            return utr3
    
    def get_utr_regions(self):
        """获取所有UTR区域"""
        return {
            'utr5': self.calculate_utr5(),
            'utr3': self.calculate_utr3(),
            'strand': self.strand
        }


# 使用示例
if __name__ == "__main__":
    exons = [(88898, 89081), (89173, 89263), (89405, 89745)]
    cds = [(88977, 89081), (89173, 89263), (89405, 89529)]
    
    calculator = UTRCalculator(exons, cds, "+")
    utr_regions = calculator.get_utr_regions()
    
    print(f"链方向: {utr_regions['strand']}")
    print(f"5'UTR: {utr_regions['utr5']}")
    print(f"3'UTR: {utr_regions['utr3']}")
