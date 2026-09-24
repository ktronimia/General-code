#!/usr/bin/env python3
# heatmap_from_mafft.py

from Bio import AlignIO
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def percent_identity(seq1, seq2, ignore_gaps=True):
    """Compute percent identity between two aligned sequences."""
    matches, length = 0, 0
    for a, b in zip(str(seq1), str(seq2)):
        if ignore_gaps and (a == "-" or b == "-"):
            continue
        length += 1
        if a == b:
            matches += 1
    return 100.0 * matches / length if length > 0 else 0

def compute_identity_matrix(alignment, ignore_gaps=True):
    """Compute full pairwise identity matrix from alignment."""
    n = len(alignment)
    names = [rec.id for rec in alignment]
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                M[i, j] = 100.0
            elif i < j:
                val = percent_identity(alignment[i].seq, alignment[j].seq, ignore_gaps)
                M[i, j] = M[j, i] = round(val, 2)
    return pd.DataFrame(M, index=names, columns=names)

def plot_heatmap(df, out_png="identity_heatmap.png", title="Percent Identity Matrix"):
    plt.figure(figsize=(20, 10))
    ax = sns.heatmap(
        df, 
        annot=True, fmt=".2f",
        cmap="coolwarm",
        cbar_kws={"label": "Percent Identity"},
        annot_kws={"size": 12}  # 🔹 font for numbers
    )

    # Titles and labels
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=12, rotation=45, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=12, rotation=0)

    plt.tight_layout()
    plt.savefig(out_png, dpi=1000)
    plt.show()
    print(f"✅ Heatmap saved as {out_png}")

if __name__ == "__main__":
    # Change "aligned.fasta" to your MAFFT alignment file name
    alignment_file = "aligned.fasta"
    aln = AlignIO.read(alignment_file, "fasta")  # or "clustal" if .aln
    df = compute_identity_matrix(aln, ignore_gaps=True)
    df.to_csv("identity_matrix.csv")
    print("✅ Identity matrix saved as identity_matrix.csv")
    plot_heatmap(df, out_png="identity_heatmap.png")

