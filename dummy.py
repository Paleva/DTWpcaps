def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
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
    print(f"max_pos: {max_pos}")
    print(H)    
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
    print(f"path: {path}")
    i_start, i_end = path[0][0], path[-1][0]
    j_start, j_end = path[0][1], path[-1][1]
    print(f"i_start: {i_start}, i_end: {i_end}")
    print(f"j_start: {j_start}, j_end: {j_end}")
    return max_score, path, (i_start, i_end), (j_start, j_end)

import numpy as np
from dtw import *


random_nums1 = np.random.randint(0, 10, 10)
random_nums2 = np.random.randint(0, 10, 10)

alignemnt = dtw(random_nums2, random_nums1, keep_internals=True)
alignemnt.plot(type="twoway")

import matplotlib.pyplot as plt
plt.show()