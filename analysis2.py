"""
pcap_alignment_analysis_v2.py

Script to analyse two PCAP files (attacker and target captures) using multiple sequence alignment and distance algorithms.
Algorithms implemented:
  - Dynamic Time Warping (DTW)
  - Euclidean matching (direct Euclidean distance)
  - Needleman-Wunsch global alignment
  - Smith-Waterman local alignment

Only TCP payload lengths are analysed, and all alignments operate on normalized values in [0,1].

Usage:
  python3 pcap_alignment_analysis_v2.py --attacker <attacker.pcap> --target <target.pcap>
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from dtaidistance import dtw
from dtaidistance import dtw_visualisation as dtwvis
from dtaidistance import alignment
from matplotlib.patches import ConnectionPatch
from scapy.all import rdpcap
from scapy.layers.inet import TCP


def extract_sequence(pcap_file):
    """
    Extract sequence of TCP payload lengths for all non-SSH TCP traffic.
    """
    packets = rdpcap(pcap_file)
    lengths = []
    for pkt in packets:
        if TCP in pkt:
            if pkt[TCP].dport == 22 or pkt[TCP].sport == 22:
                continue
            data = bytes(pkt[TCP].payload)
            lengths.append(len(data) if data else 0)
    return lengths


def build_alignment_path(aligned1, aligned2):
    """
    Given two gap-padded aligned sequences of equal length,
    return a list of (i,j) pairs mapping positions in each original sequence.
    Gaps ('-') map to None.
    """
    path = []
    i = j = 0
    for a, b in zip(aligned1, aligned2):
        ai = i if a != '-' else None
        bj = j if b != '-' else None
        path.append((ai, bj))
        if a != '-':
            i += 1
        if b != '-':
            j += 1
    return path


def compute_dtw(seq1, seq2):
    """
    Compute DTW distance and best warping path on normalized sequences.
    Returns (distance, path).
    """
    distance = dtw.distance(seq1, seq2)
    _, paths = dtw.warping_paths(seq1, seq2, use_c=False)
    best_path = dtw.best_path(paths)
    return distance, best_path


def needleman_global(seq1, seq2, gap=0):
    """
    Global alignment (Needleman-Wunsch). Returns (score, aligned1, aligned2, path).
    Uses dtaidistance.aligners for scoring, then builds a uniform path.
    """
    score, scores, pointers = alignment.needleman_wunsch(seq1, seq2)
    _, aligned1, aligned2 = alignment.best_alignment(pointers, seq1, seq2, gap=gap)
    path = build_alignment_path(aligned1, aligned2)
    return score, aligned1, aligned2, path


def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Local alignment (Smith-Waterman). Returns (score, aligned1, aligned2, path).
    Operates on normalized input sequences.
    """
    n, m = len(seq1), len(seq2)
    H = np.zeros((n+1, m+1), dtype=int)
    ptr = np.zeros((n+1, m+1), dtype=object)

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n+1):
        for j in range(1, m+1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            diag = H[i-1, j-1] + s
            up   = H[i-1, j]   + gap
            left = H[i, j-1]   + gap
            H[i, j] = max(0, diag, up, left)
            if H[i, j] == diag:
                ptr[i, j] = 'd'
            elif H[i, j] == up:
                ptr[i, j] = 'u'
            elif H[i, j] == left:
                ptr[i, j] = 'l'
            if H[i, j] > max_score:
                max_score = H[i, j]
                max_pos = (i, j)

    # Traceback local region
    aligned1, aligned2 = [], []
    i, j = max_pos
    while i > 0 and j > 0 and H[i, j] > 0:
        move = ptr[i, j]
        if move == 'd':
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i, j = i-1, j-1
        elif move == 'u':
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        elif move == 'l':
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
        else:
            break

    i, j = max_pos
    path = []
    while i>0 and j>0 and H[i,j] > 0:
        path.append((i-1, j-1))
        score_diag   = H[i-1,j-1]
        score_up     = H[i-1,j]
        score_left   = H[i,j-1]
        s = match if seq1[i-1] == seq2[j-1] else mismatch

        if H[i,j] == score_diag + s:
            i, j = i-1, j-1
        elif H[i,j] == score_up + gap:
            i -= 1
        else:
            j -= 1

    aligned1.reverse()
    aligned2.reverse()
    path.reverse()
    # path = build_alignment_path(aligned1, aligned2)
    return max_score, aligned1, aligned2, path


def plot_alignment(seq1, seq2, path=None,
                   label1='Seq1', label2='Seq2',
                   ylabel='Value', max_lines=200, filename='alignment.png'):
    """
    Plot two gap-padded sequences and overlay alignment lines.
    Gaps ('-') are converted to np.nan.
    """
    plot1 = [x if x != '-' else np.nan for x in seq1]
    plot2 = [x if x != '-' else np.nan for x in seq2]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=False)
    ax1.plot(plot1, linewidth=2)
    ax1.set_title(f'{label1} (aligned)')
    ax1.set_ylabel(ylabel)
    ax1.grid(True, axis='x', linestyle='--', alpha=0.5)

    ax2.plot(plot2, linewidth=2)
    ax2.set_title(f'{label2} (aligned)')
    ax2.set_ylabel(ylabel)
    ax2.set_xlabel('Alignment index')
    ax2.grid(True, axis='x', linestyle='--', alpha=0.5)

    if path is not None:
        for idx, (i, j) in enumerate(path):
            if idx >= max_lines:
                break
            if i is None or j is None:
                continue
            # clamp indices
            if i < 0 or i >= len(plot1) or j < 0 or j >= len(plot2):
                continue
            con = ConnectionPatch(
                xyA=(i, plot1[i]), coordsA=ax1.transData,
                xyB=(j, plot2[j]), coordsB=ax2.transData,
                color='orange', alpha=0.5
            )
            fig.add_artist(con)

    plt.tight_layout()
    plt.savefig(f'figures/{filename}', dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='PCAP alignment analysis')
    parser.add_argument('--attacker', required=True)
    parser.add_argument('--target', required=True)
    args = parser.parse_args()

    seq1 = extract_sequence(args.attacker)
    seq2 = extract_sequence(args.target)
    print(f"Packets: attacker={len(seq1)}, target={len(seq2)}")

    # Normalize both sequences by overall max
    global_max = max(max(seq1) if seq1 else 1, max(seq2) if seq2 else 1)
    norm1 = [x/ global_max for x in seq1]
    norm2 = [x/ global_max for x in seq2]

    # NW global
    nw_score, nw1, nw2, nw_path = needleman_global(seq1, seq2)
    print(f"NW best score      = {nw_score:.3f}")
    plot_alignment(nw1, nw2, nw_path,
                   label1='Attacker (NW)', label2='Target (NW)',
                   ylabel='Norm payload', filename='nw_alignment.png')

    # SW local
    sw_score, sw1, sw2, sw_path = smith_waterman(seq1, seq2)
    print(f"SW best score      = {sw_score:.3f}")
    plot_alignment(sw1, sw2, sw_path,
                   label1='Attacker (SW)', label2='Target (SW)',
                   ylabel='Norm payload', filename='sw_alignment.png')

    # DTW
    dtw_dist, dtw_path = compute_dtw(seq1, seq2)
    print(f"DTW distance       = {dtw_dist:.3f}")

    dtwvis.plot_warping(norm1, norm2, dtw_path, filename="figures/dtw_alignment.png")

if __name__ == '__main__':
    main()
