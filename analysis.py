import argparse
import numpy as np
import matplotlib.pyplot as plt
from dtaidistance import dtw
from dtaidistance import dtw_visualisation as dtwvis
from dtaidistance import alignment
from matplotlib.patches import ConnectionPatch
from scapy.all import rdpcap
from scapy.layers.inet import TCP, IP


def extract_sequence(pcap_file):
    """
    Extract sequence of TCP payload lengths and inter-arrival times for all TCP traffic.
    """
    packets = rdpcap(pcap_file)
    lengths = []
    for pkt in packets:
        if TCP in pkt:
            if pkt[TCP].dport == 22 or pkt[TCP].sport == 22:
                # Ignore SSH traffic
                continue
            data = bytes(pkt[TCP].payload)
            lengths.append(len(data) if data else 0)
    return lengths

def compute_dtw(seq1, seq2):
    """
    Compute DTW distance and obtain the warping paths.
    """
    distance = dtw.distance(seq1, seq2)
    _, paths = dtw.warping_paths(seq1, seq2, use_c=False)
    best_path = dtw.best_path(paths)
    return distance, best_path

def needleman_wunsch(seq1, seq2):
    values, scores, paths = alignment.needleman_wunsch(seq1, seq2)
    algn, s1a, s2a, = alignment.best_alignment(paths, seq1, seq2, gap=0)
    return values, s1a, s2a, algn


def smith_waterman(attacker_sequence, victim_sequence):
    from enum import IntEnum
    class Trace(IntEnum):
        STOP = 0
        LEFT = 1 
        UP = 2
        DIAGONAL = 3

    class Score(IntEnum):
        MATCH = 1
        MISMATCH = -1
        GAP = -1


    row, col = len(attacker_sequence) + 1, len(victim_sequence) + 1
    matrix = np.zeros(shape=(row, col), dtype=int)
    tracing_matrix = np.zeros(shape=(row, col), dtype=int)

    max_score = -1
    max_index = (-1, -1)

    # Fill the matrix
    for i in range(1, row):
        for j in range(1, col):
            val1 = attacker_sequence[i - 1]
            val2 = victim_sequence[j - 1]

            match_value = Score.MATCH if val1 == val2 else Score.MISMATCH
            diagonal = matrix[i - 1, j - 1] + match_value
            vertical = matrix[i - 1, j] + Score.GAP
            horizontal = matrix[i, j - 1] + Score.GAP

            matrix[i, j] = max(0, diagonal, vertical, horizontal)

            if matrix[i, j] == 0:
                tracing_matrix[i, j] = Trace.STOP
            elif matrix[i, j] == diagonal:
                tracing_matrix[i, j] = Trace.DIAGONAL
            elif matrix[i, j] == vertical:
                tracing_matrix[i, j] = Trace.UP
            else:
                tracing_matrix[i, j] = Trace.LEFT

            if matrix[i, j] > max_score:
                max_score = matrix[i, j]
                max_index = (i, j)

    # Traceback
    align1, align2 = [], []
    i, j = max_index
    while tracing_matrix[i, j] != Trace.STOP:
        if tracing_matrix[i, j] == Trace.DIAGONAL:
            align1.append(attacker_sequence[i - 1])
            align2.append(victim_sequence[j - 1])
            i -= 1
            j -= 1
        elif tracing_matrix[i, j] == Trace.UP:
            align1.append(attacker_sequence[i - 1])
            align2.append(0)
            i -= 1
        else:
            align1.append(0)
            align2.append(victim_sequence[j - 1])
            j -= 1

    return max_score, align1[::-1], align2[::-1]


def plot_alignment(seq1, seq2, label1='Seq1', label2='Seq2', ylabel='Value', max_lines=200, filename='alignment_plot.png'):

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
    
    plt.tight_layout()
    plt.savefig(f'figures/{filename}', dpi=300)

def plot_dtw(seq1, seq2, path=None, label1='Seq1', label2='Seq2', filename='dtw_original.png'):
    fig1, axs1 = plt.subplots(2, 1, figsize=(12, 8))
    
    axs1[0].set_title(label1)
    axs1[1].set_title(label2)
    dtwvis.plot_warping(seq1, seq2, path, fig=fig1, axs=axs1, filename=f"figures/{filename}")
    print(f"Saved visualization to figures/{filename}")

def main():
    parser = argparse.ArgumentParser(description='PCAP alignment analysis')
    parser.add_argument('--attacker', required=True)
    parser.add_argument('--target', required=True)
    args = parser.parse_args()

    seq1 = extract_sequence(args.attacker)
    seq2 = extract_sequence(args.target)
    print(f"Packets: attacker={len(seq1)}, target={len(seq2)}")
    
    # --- PAYLOAD LENGTHS ---
    nw, nw_seq1, nw_seq2, nw_path = needleman_wunsch(seq1, seq2)

    distance, path = compute_dtw(seq1, seq2)
    
    print(f"NW best score: {nw:.2f}")
    print(f"DTW distance: {distance:.2f}")
    print(f"DTW alignment path length: {len(path)} pairs aligned")
    plot_alignment(nw_seq1, nw_seq2, label1='Attacker Payloads (NW)', label2='Target Payloads (NW)', ylabel='Payload Length', filename="nw_alignment_plot_streched.png")
    plot_dtw(seq1, seq2, nw_path, label1='Attacker Payloads (NW)', label2='Target Payloads (NW)', filename="nw_alignment_plot.png")
    plot_dtw(seq1, seq2, path, label1='Attacker Payloads (DTW)', label2='Target Payloads (DTW)', filename="dtw_original.png")


    print('-----SUBSEQUENCE-----')
    print('Following subsequences were extracted using Smith-Waterman that have the best match between them')
    
    # Extract subsequence using Smith-Waterman
    score, sw_seq1, sw_seq2 = smith_waterman(seq1, seq2)
    print(f"SW best score: {score}")
    print(f"Subseq lengths: attacker={len(sw_seq1)}, target={len(sw_seq2)}")    

    sub_distance, dtw_path_sub = compute_dtw(sw_seq1, sw_seq2)
    nw_score, _, _, nw_path = needleman_wunsch(sw_seq1, sw_seq2)

    print(f"NW best score (sub): {nw_score:.2f}")
    print(f"DTW distance (sub): {sub_distance:.2f}")
    print(f"DTW alignment path length: {len(dtw_path_sub)} pairs aligned")
    plot_alignment(sw_seq1, sw_seq2, label1='Attacker Subseq (SW)', label2='Target Subseq (SW)', ylabel='Payload Length', filename="sw_subseq_alignment_plot.png")
    plot_dtw(sw_seq1, sw_seq2, nw_path, label1='Attacker Subseq (SW)', label2='Target Subseq (SW)', filename="sw_subseq_alignment_plot.png")
    plot_dtw(sw_seq1, sw_seq2, dtw_path_sub, label1='Attacker Subseq (DTW)', label2='Target Subseq (DTW)', filename="sw_subseq_dtw_plot.png")
if __name__ == '__main__':
    main()
