import sys


class Peak(object):
    def __init__(self, chrom, pos0, pos1, score, strand) -> None:
        self.chrom = chrom
        self.pos0 = int(pos0)
        self.pos1 = int(pos1)
        self.score = float(score)
        self.strand = strand


class DivergentPeakPair(object):
    def __init__(self, line) -> None:
        (
            m_chr, m_pos0, m_pos1, m_name, m_score, m_strand, 
            p_chr, p_pos0, p_pos1, p_name, p_score, p_strand, 
        ) = line.strip().split('\t')[:12]
        self.m_peak = Peak(m_chr, m_pos0, m_pos1, m_score, m_strand)
        self.p_peak = Peak(p_chr, p_pos0, p_pos1, p_score, p_strand)


class PairSet:
    def __init__(self, pair: DivergentPeakPair) -> None:
        self.chrom = pair.m_peak.chrom
        self.range_m = [pair.m_peak.pos0, pair.m_peak.pos1]
        self.range_p = [pair.p_peak.pos0, pair.p_peak.pos1]
        self.m_peaks = [pair.m_peak]
        self.p_peaks = [pair.p_peak]

    def merge_if_intersect(self, pair: DivergentPeakPair):
        if pair.m_peak.chrom != self.chrom or pair.m_peak.pos0 >= self.range_p[1]:
            return False
        self.m_peaks.append(pair.m_peak)
        self.p_peaks.append(pair.p_peak)
        
        if pair.m_peak.pos0 < self.range_m[0]:
            self.range_m[0] = pair.m_peak.pos0
        if pair.m_peak.pos1 > self.range_m[1]:
            self.range_m[1] = pair.m_peak.pos1
        if pair.p_peak.pos0 < self.range_p[0]:
            self.range_p[0] = pair.p_peak.pos0
        if pair.p_peak.pos1 > self.range_p[1]:
            self.range_p[1] = pair.p_peak.pos1
        return True

    def export(self):
        self.m_peaks.sort(key=lambda x: x.pos1)
        self.p_peaks.sort(key=lambda x: x.pos0, reverse=True)
        while self.range_p[0] < self.range_m[1]:
            if self.m_peaks[-1].score > self.p_peaks[-1].score:
                if len(self.p_peaks) < 2 and len(self.m_peaks) > 1:
                    self.m_peaks.pop()
                elif len(self.p_peaks) > 1:
                    self.p_peaks.pop()
            else:
                if len(self.m_peaks) < 2 and len(self.p_peaks) > 1:
                    self.p_peaks.pop()
                elif len(self.m_peaks) > 1:
                    self.m_peaks.pop()
            self.range_p[0] = self.p_peaks[-1].pos0
            self.range_m[1] = self.m_peaks[-1].pos1
        
        highest_score = [
            max([x.score for x in self.m_peaks]),
            max([x.score for x in self.p_peaks])
        ]
        center = int(0.5 * (self.range_m[1] - 1 + self.range_p[0]))
        return "\t".join([
            self.chrom, 
            str(self.range_m[0]), 
            str(self.range_p[1]), 
            f"{self.chrom}:{self.range_m[0]}-{self.range_p[1]}", 
            str(int(min(highest_score))), 
            ".", 
            str(center), 
            str(center + 1), 
            "0,0,0", 
            "2",
            f"{self.range_m[1] - self.range_m[0]},{self.range_p[1] - self.range_p[0]}",
            f"0,{self.range_p[0] - self.range_m[0]}"
        ])

current_pair_set = None

for line in sys.stdin:
    pair = DivergentPeakPair(line)

    if not current_pair_set:
        current_pair_set = PairSet(pair)
        continue
    result = current_pair_set.merge_if_intersect(pair)

    if not result:
        print(current_pair_set.export())
        current_pair_set = PairSet(pair)

if current_pair_set:
    print(current_pair_set.export())    
