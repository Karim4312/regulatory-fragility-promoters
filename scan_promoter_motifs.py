import re
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd


MOTIF_PATTERNS = {
    'CCAAT': re.compile(r'CCAAT'),
    'TATA': re.compile(r'TATA[AT]A'),
    'GC_box': re.compile(r'GGGCGG|CCGCCC')
}


def find_motifs(sequence, patterns):
    results = {}
    for motif_name, pattern in patterns.items():
        results[motif_name] = [m.start() for m in pattern.finditer(sequence)]
    return results


def collapse_positions(positions, window=5):
    if not positions:
        return []
    clusters = [[positions[0]]]
    for pos in positions[1:]:
        if pos - clusters[-1][-1] <= window:
            clusters[-1].append(pos)
        else:
            clusters.append([pos])
    return [min(cluster) for cluster in clusters]


def scan_promoter_motifs(fasta_path, output_tsv):
    results = []

    for record in SeqIO.parse(fasta_path, 'fasta'):
        seq = str(record.seq)
        pid = record.id

        forward = find_motifs(seq, MOTIF_PATTERNS)
        reverse_seq = str(Seq(seq).reverse_complement())
        reverse = find_motifs(reverse_seq, MOTIF_PATTERNS)

        combined = {}
        for name, pattern in MOTIF_PATTERNS.items():
            fwd = forward[name]
            rev = [len(seq) - p - len(pattern.pattern) for p in reverse[name]]
            combined[name] = collapse_positions(sorted(fwd + rev))

        row = {'promoter_id': pid}
        for name, pos in combined.items():
            row[f'{name}_count'] = len(pos)
            row[f'{name}_positions'] = ';'.join(map(str, pos)) if pos else ''

        results.append(row)

    pd.DataFrame(results).to_csv(output_tsv, sep='\t', index=False)


if __name__ == "__main__":
    scan_promoter_motifs(
        fasta_path="promoters.fa",
        output_tsv="promoter_motifs.tsv"
    )
