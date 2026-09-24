#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import AlignIO
import math


# -------------------------------
# Percent identity
# -------------------------------
def percent_identity(seq1, seq2):
    matches, valid = 0, 0

    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            valid += 1
            if a == b:
                matches += 1

    return (matches / valid) * 100 if valid > 0 else 0


# -------------------------------
# p-distance
# -------------------------------
def p_distance(seq1, seq2):
    mismatches, valid = 0, 0

    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            valid += 1
            if a != b:
                mismatches += 1

    return (mismatches / valid) * 100 if valid > 0 else 0


# -------------------------------
# Kimura 2-parameter distance
# -------------------------------
def kimura_2p(seq1, seq2):
    transitions = 0
    transversions = 0
    valid = 0

    transitions_pairs = [('A','G'), ('G','A'), ('C','T'), ('T','C')]

    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            valid += 1
            if a != b:
                if (a, b) in transitions_pairs:
                    transitions += 1
                else:
                    transversions += 1

    if valid == 0:
        return 0

    P = transitions / valid
    Q = transversions / valid

    try:
        distance = -0.5 * np.log(1 - 2*P - Q) - 0.25 * np.log(1 - 2*Q)
        return distance * 100
    except:
        return np.nan


# -------------------------------
# Compute matrices
# -------------------------------
def compute_matrices(seqs, model):
    n = len(seqs)

    identity = np.zeros((n, n))
    divergence = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            pid = percent_identity(seqs[i], seqs[j])
            identity[i, j] = pid

            if model == "p":
                divergence[i, j] = p_distance(seqs[i], seqs[j])
            elif model == "k2p":
                divergence[i, j] = kimura_2p(seqs[i], seqs[j])
            else:
                divergence[i, j] = 100 - pid

    return identity, divergence


# -------------------------------
# Plot heatmap with triangles
# -------------------------------
def plot_heatmap(identity, divergence, names, args):

    n = len(names)

    # Masks
    mask_upper = np.tril(np.ones_like(identity, dtype=bool))
    mask_lower = np.triu(np.ones_like(divergence, dtype=bool))

    # Dynamic font size
    if n <= 10:
        annot_size = 10
    elif n <= 20:
        annot_size = 8
    else:
        annot_size = 6

    fig, ax = plt.subplots(figsize=(10, 8))

    # Create separate axes for colorbars
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)

    cax1 = divider.append_axes("right", size="3%", pad=0.1)   # Identity bar
    cax2 = divider.append_axes("right", size="3%", pad=0.6)   # Divergence bar

    # -------------------------------
    # Upper triangle → Identity (RED)
    # -------------------------------
    hm1 = sns.heatmap(
        identity,
        mask=mask_upper,
        cmap=sns.color_palette("Reds", as_cmap=True),
        vmin=np.nanmin(identity),
        vmax=100,
        annot=True,
        fmt=".1f",
        annot_kws={"size": annot_size},
        linewidths=0.5,
        linecolor='black',
        cbar=True,
        cbar_ax=cax1,
        ax=ax
    )

    hm1.collections[0].colorbar.set_label("Identity (%)")

    # -------------------------------
    # Lower triangle → Divergence (BLUE)
    # -------------------------------
    hm2 = sns.heatmap(
        divergence,
        mask=mask_lower,
        cmap=sns.color_palette("Blues", as_cmap=True),
        vmin=0,
        vmax=np.nanmax(divergence),
        annot=True,
        fmt=".1f",
        annot_kws={"size": annot_size},
        linewidths=0.5,
        linecolor='black',
        cbar=True,
        cbar_ax=cax2,
        ax=ax
    )

    hm2.collections[0].colorbar.set_label("k2p Distance (%)")

    # Axis labels
    ax.set_xticks(np.arange(n)+0.5)
    ax.set_yticks(np.arange(n)+0.5)
    ax.set_xticklabels(names, rotation=45, ha='right')
    ax.set_yticklabels(names, rotation=0)

    # Diagonal stars
    for i in range(n):
        ax.text(i+0.5, i+0.5, "*",
                ha='center', va='center',
                fontsize=annot_size+2, color='black')

    ax.set_title("Percent Identity (Red, Upper) & Distance (Blue, Lower)", fontsize=14)

    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi, format=args.format)
    plt.close()



# -------------------------------
# Main
# -------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Advanced heatmap: identity (upper) + divergence (lower)"
    )

    parser.add_argument("-i", "--input", required=True, help="Aligned FASTA")
    parser.add_argument("-o", "--output", default="heatmap.png", help="Output image")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--format", default="png", choices=["png", "pdf", "svg"])
    parser.add_argument("--model", default="simple",
                        choices=["simple", "p", "k2p"],
                        help="Divergence model")

    args = parser.parse_args()

    alignment = AlignIO.read(args.input, "fasta")

    names = [rec.id for rec in alignment]
    seqs = [str(rec.seq) for rec in alignment]

    identity, divergence = compute_matrices(seqs, args.model)

    # Save combined matrix
    combined = np.full_like(identity, np.nan)

    for i in range(len(seqs)):
        for j in range(len(seqs)):
            if i < j:
                combined[i, j] = identity[i, j]
            elif i > j:
                combined[i, j] = divergence[i, j]

    df = pd.DataFrame(combined, index=names, columns=names)
    df.to_csv("matrix_combined.csv")

    plot_heatmap(identity, divergence, names, args)

    print("✅ Done!")
    print(f"📊 Saved: {args.output}")
    print("📄 Saved: matrix_combined.csv")


if __name__ == "__main__":
    main()
