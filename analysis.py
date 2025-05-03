"""
pcap_alignment_analysis.py

Script to analyse two PCAP files (attacker and target captures) using multiple sequence alignment and distance algorithms.
Algorithms implemented:
  - Dynamic Time Warping (DTW)
  - Euclidean matching (direct Euclidean distance)
  - Needleman-Wunsch global alignment
  - Smith-Waterman local alignment

Only TCP payload lengths are analysed across the full trace.

Outputs:
  - Alignment / distance scores for each algorithm
  - Two-panel alignment plot (attacker vs target) with matching lines

Dependencies:
  pip install scapy numpy matplotlib

Usage:
  python3 pcap_alignment_analysis.py --attacker attacker.pcap --target target.pcap
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from dtaidistance import dtw
from dtaidistance import dtw_visualisation as dtwvis
from matplotlib.patches import ConnectionPatch
from scapy.all import rdpcap
from scapy.layers.inet import TCP, IP


def extract_sequence(pcap_file):
    """
    Extract sequence of TCP payload lengths and inter-arrival times for all TCP traffic.
    """
    packets = rdpcap(pcap_file)
    lengths = []
    times = []
    for pkt in packets:
        if TCP in pkt:
            if pkt[TCP].dport == 22 or pkt[TCP].sport == 22:
                # Ignore SSH traffic
                continue
            data = bytes(pkt[TCP].payload)
            lengths.append(len(data) if data else 0)
    return lengths

def euclidean_distance(seq1, seq2):
    n = min(len(seq1), len(seq2))
    diff = np.array(seq1[:n]) - np.array(seq2[:n])
    return np.sqrt(np.sum(diff ** 2))

def compute_dtw(seq1, seq2):
    """
    Compute DTW distance and obtain the warping paths.
    """
    distance = dtw.distance(seq1, seq2)
    dtw_matrix, paths = dtw.warping_paths(seq1, seq2, use_c=False)
    dtwvis.plot_warpingpaths(seq1, seq2, paths, filename="figures/dtw_warping_paths.png")
    best_path = dtw.best_path(paths)
    return distance, best_path

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    n, m = len(seq1), len(seq2)
    score = np.zeros((n+1, m+1), dtype=int)
    pointer = np.empty((n+1, m+1), dtype=object)

    # initialize first column and first row
    for i in range(1, n+1):
        score[i, 0] = score[i-1, 0] + gap
        pointer[i, 0] = 'u'  # up
    for j in range(1, m+1):
        score[0, j] = score[0, j-1] + gap
        pointer[0, j] = 'l'  # left

    pointer[0, 0] = None

    # fill in the score matrix and pointer matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            diag = score[i-1, j-1] + s
            up   = score[i-1, j] + gap
            left = score[i, j-1] + gap
            max_score = max(diag, up, left)
            score[i, j] = max_score
            if max_score == diag:
                pointer[i, j] = 'd'
            elif max_score == up:
                pointer[i, j] = 'u'
            else:
                pointer[i, j] = 'l'

    # traceback to build the alignment path
    i, j = n, m
    path = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and pointer[i, j] == 'd':
            path.append((i-1, j-1))
            i, j = i-1, j-1
        elif i > 0 and pointer[i, j] == 'u':
            path.append((i-1, None))
            i -= 1
        elif j > 0 and pointer[i, j] == 'l':
            path.append((None, j-1))
            j -= 1
    path.reverse()

    return score[n, m], path

def smith_waterman(seq1, seq2, match=5, mismatch=-10, gap=-1):
    n, m = len(seq1), len(seq2)
    H = np.zeros((n+1, m+1), dtype=int)
    max_score, max_pos = 0, (0,0)
    for i in range(1, n+1):
        for j in range(1, m+1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            H[i,j] = max(
                0,
                H[i-1,j-1] + s,
                H[i-1,j]   + gap,
                H[i,j-1]   + gap
            )
            if H[i,j] > max_score:
                max_score, max_pos = H[i,j], (i,j)

    # Backtrack from max_pos
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

    path.reverse()
    i_start, i_end = path[0][0], path[-1][0]
    j_start, j_end = path[0][1], path[-1][1]
    return max_score, path, (i_start, i_end), (j_start, j_end)


def plot_alignment(seq1, seq2, path, label1='Seq1', label2='Seq2', ylabel='Value', max_lines=200, filename='alignment_plot.png'):

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=False)
    ax1.plot(seq1, linewidth=2)
    ax1.set_ylabel(ylabel)
    ax1.set_title(f'{label1} (packet index)')
    ax1.grid(True, axis='x', linestyle='--', alpha=0.5)

    ax2.plot(seq2, linewidth=2)
    ax2.set_ylabel(ylabel)
    ax2.set_xlabel('Series Index')
    ax2.set_title(f'{label2} (packet index)')
    ax2.grid(True, axis='x', linestyle='--', alpha=0.5)
    

    total = len(path)
    idxs = np.linspace(0, total-1, min(total, max_lines), dtype=int)
    for k in idxs:
        i, j = path[k]
        if i is None or j is None:
            continue
        con = ConnectionPatch(
            xyA=(i, seq1[i]), coordsA=ax1.transData,
            xyB=(j, seq2[j]), coordsB=ax2.transData,
            color='orange', alpha=0.6
        )
        fig.add_artist(con)
    plt.tight_layout()
    plt.savefig(f'figures/{filename}', dpi=300)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description='PCAP alignment analysis')
    parser.add_argument('--attacker', required=True)
    parser.add_argument('--target', required=True)
    args = parser.parse_args()

    seq1 = extract_sequence(args.attacker)
    print(seq1)
    seq2 = extract_sequence(args.target)
    print(seq2)
    print(f"Packets: attacker={len(seq1)}, target={len(seq2)}")
    
    max_seq1 = max(seq1)
    max_seq2 = max(seq2)
    norm_seq1 = np.array(seq1) / max_seq1
    norm_seq2 = np.array(seq2) / max_seq2
    
    # --- PAYLOAD LENGTHS ---
    nc, nc_path = needleman_wunsch(norm_seq1, norm_seq2)
    score = smith_waterman(norm_seq1, norm_seq2)
    distance, path = compute_dtw(norm_seq1, norm_seq2)
    
    print(f"NW best score: {nc:.2f}")
    print(f"SW best score: {score[0]:.2f}")
    print(f"DTW distance: {distance:.2f}")
    print(f"DTW alignment path length: {len(path)} pairs aligned")
    plot_alignment(seq1, seq2, nc_path, label1='Attacker Payloads (NW)', label2='Target Payloads (NW)', ylabel='Payload Length', filename="nw_alignment_plot.png")
    plot_alignment(seq1, seq2, path, label1='Attacker Payloads (DTW)', label2='Target Payloads (DTW)', ylabel='Payload Length', filename="dtw_alignment_plot.png")
    dtwvis.plot_warping(seq1, seq2, path, filename="figures/dtw_alignment_plot(dtai).png")


    print('-----SUBSEQUENCE-----')
    print('Following subsequences were extracted using Smith-Waterman that have the best match between them')
    # Extract subsequence using Smith-Waterman
    score, sw_path, (i0, i1), (j0, j1) = smith_waterman(seq1, seq2)
    print(f"SW best score: {score}")  
    print(f"Seq1 match indices: {i0}-{i1}")  
    print(f"Seq2 match indices: {j0}-{j1}")
    # Extract subsequence using the indices from Smith-Waterman
    # Note: i0, i1, j0, j1 are inclusive  
    sw_seq1 = seq1[i0:i1+1]
    sw_seq2 = seq2[j0:j1+1]
    max_sw1 = max(sw_seq1)
    max_sw2 = max(sw_seq2)
    norm_sw_seq1 = np.array(sw_seq1) / max_sw1
    norm_sw_seq2 = np.array(sw_seq2) / max_sw2

    dtw_matrix_sub, dtw_path_sub = compute_dtw(norm_sw_seq1, norm_sw_seq2)
    nw_sub, nw_sub_path = needleman_wunsch(norm_sw_seq1, norm_sw_seq2)

    print(f"NW best score (sub): {nw_sub:.2f}")
    print(f"DTW distance (sub): {dtw_matrix_sub:.2f}")
    print(f"DTW alignment path length: {len(dtw_path_sub)} pairs aligned")
    plot_alignment(sw_seq1, sw_seq2, nw_sub_path, label1='Attacker Subseq (NW)', label2='Target Subseq (NW)', ylabel='Payload Length', filename="nw_subseq_alignment_plot.png")
    plot_alignment(sw_seq1, sw_seq2, dtw_path_sub, label1='Attacker Subseq (DTW)', label2='Target Subseq (DTW)', ylabel='Payload Length', filename="dtw_subseq_alignment_plot.png")
    dtwvis.plot_warping(sw_seq1, sw_seq2, dtw_path_sub, filename="figures/dtw_subseq_alignment_plot(dtai).png")

if __name__ == '__main__':
    main()
