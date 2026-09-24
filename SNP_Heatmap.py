#!/usr/bin/env python3

# ============================================================
# INSTALL DEPENDENCIES
# ============================================================
# conda install -c conda-forge biopython pandas numpy matplotlib seaborn


import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import AlignIO
from matplotlib.colors import LinearSegmentedColormap


# ============================================================
# USER SETTINGS
# ============================================================

ALIGNMENT_FILE = "alignment.fasta"

OUTPUT_MATRIX_CSV = "SNP_similarity_matrix.csv"
OUTPUT_LONG_CSV = "SNP_pairwise_similarity.csv"
OUTPUT_SNP_POSITIONS_CSV = "SNP_positions.csv"

OUTPUT_HEATMAP_PNG = "SNP_similarity_heatmap.png"
OUTPUT_HEATMAP_PDF = "SNP_similarity_heatmap.pdf"

FIGURE_DPI = 600

LABEL_FONT_SIZE = 7

# ------------------------------------------------------------
# FIRST SEQUENCE = REFERENCE STRAIN
# ------------------------------------------------------------

REFERENCE_INDEX = 0

# ------------------------------------------------------------
# VALID NUCLEOTIDES
# ------------------------------------------------------------

VALID_BASES = set("ACGT")


# ============================================================
# READ ALIGNMENT
# ============================================================

print("=" * 70)
print("Reading alignment...")
print("=" * 70)

try:

    alignment = AlignIO.read(
        ALIGNMENT_FILE,
        "fasta"
    )

except Exception as e:

    print(
        f"ERROR: Cannot read alignment file: {e}"
    )

    sys.exit(1)


n_sequences = len(
    alignment
)

alignment_length = (
    alignment.get_alignment_length()
)

print(
    f"Number of sequences: "
    f"{n_sequences}"
)

print(
    f"Alignment length: "
    f"{alignment_length:,} bp"
)


# ============================================================
# EXTRACT SEQUENCE IDS AND SEQUENCES
# ============================================================

sequence_names = [

    record.id

    for record in alignment

]

sequences = [

    str(record.seq).upper()

    for record in alignment

]

seq_array = np.array([

    list(seq)

    for seq in sequences

])


reference_name = sequence_names[
    REFERENCE_INDEX
]

print(
    f"Reference / first strain: "
    f"{reference_name}"
)


# ============================================================
# IDENTIFY SNP POSITIONS
# ============================================================

print(
    "\nIdentifying SNP positions..."
)

snp_positions = []

snp_information = []


for pos in range(
    alignment_length
):

    column = seq_array[
        :,
        pos
    ]

    # --------------------------------------------------------
    # RETAIN ONLY STANDARD NUCLEOTIDES
    # --------------------------------------------------------

    valid_column = [

        base

        for base in column

        if base in VALID_BASES

    ]

    # --------------------------------------------------------
    # SNP = AT LEAST TWO DIFFERENT NUCLEOTIDES
    # --------------------------------------------------------

    unique_bases = sorted(
        set(valid_column)
    )

    if len(
        unique_bases
    ) >= 2:

        snp_positions.append(
            pos
        )

        base_counts = {

            base: int(
                np.sum(
                    column == base
                )
            )

            for base in [
                "A",
                "C",
                "G",
                "T"
            ]

        }

        snp_information.append({

            "Alignment_Position":
                pos + 1,

            "Reference_Position":
                pos + 1,

            "Bases_Observed":
                ",".join(
                    unique_bases
                ),

            "A_Count":
                base_counts["A"],

            "C_Count":
                base_counts["C"],

            "G_Count":
                base_counts["G"],

            "T_Count":
                base_counts["T"]

        })


n_snps = len(
    snp_positions
)

print(
    f"Total SNP sites identified: "
    f"{n_snps:,}"
)


# ============================================================
# SAVE SNP POSITION INFORMATION
# ============================================================

snp_info_df = pd.DataFrame(
    snp_information
)


if not snp_info_df.empty:

    snp_info_df.to_csv(

        OUTPUT_SNP_POSITIONS_CSV,

        index=False

    )

else:

    pd.DataFrame(

        columns=[

            "Alignment_Position",
            "Reference_Position",
            "Bases_Observed",
            "A_Count",
            "C_Count",
            "G_Count",
            "T_Count"

        ]

    ).to_csv(

        OUTPUT_SNP_POSITIONS_CSV,

        index=False

    )


print(
    f"SNP position file saved: "
    f"{OUTPUT_SNP_POSITIONS_CSV}"
)


# ============================================================
# CHECK SNPs
# ============================================================

if n_snps == 0:

    print(
        "\nNo SNP sites were detected."
    )

    print(
        "All sequences are identical at valid nucleotide positions."
    )

    sys.exit(0)


# ============================================================
# EXTRACT SNP-ONLY ALIGNMENT
# ============================================================

print(
    "\nExtracting SNP-site alignment..."
)

snp_array = seq_array[
    :,
    snp_positions
]


print(

    f"SNP alignment dimensions: "

    f"{snp_array.shape[0]} sequences × "

    f"{snp_array.shape[1]} SNP sites"

)


# ============================================================
# CALCULATE PAIRWISE SNP SIMILARITY
# ============================================================

print(
    "\nCalculating pairwise SNP similarity..."
)


similarity_matrix = np.zeros(

    (
        n_sequences,
        n_sequences
    ),

    dtype=float

)


pairwise_results = []


for i in range(
    n_sequences
):

    seq1 = snp_array[i]


    for j in range(
        i,
        n_sequences
    ):

        seq2 = snp_array[j]


        # ----------------------------------------------------
        # COMPARE ONLY A/C/G/T POSITIONS
        # ----------------------------------------------------

        valid_mask = np.array([

            (a in VALID_BASES)

            and

            (b in VALID_BASES)

            for a, b

            in zip(
                seq1,
                seq2
            )

        ])


        comparable_positions = int(

            np.sum(
                valid_mask
            )

        )


        if comparable_positions == 0:

            similarity = np.nan

            identical_sites = 0

            different_sites = 0


        else:

            seq1_valid = seq1[
                valid_mask
            ]

            seq2_valid = seq2[
                valid_mask
            ]


            identical_sites = int(

                np.sum(

                    seq1_valid
                    ==
                    seq2_valid

                )

            )


            different_sites = (

                comparable_positions

                -

                identical_sites

            )


            similarity = (

                identical_sites

                /

                comparable_positions

            ) * 100


            # ------------------------------------------------
            # PREVENT NON-IDENTICAL PAIRS FROM
            # DISPLAYING AS 100.00%
            # ------------------------------------------------

            if (

                different_sites > 0

                and

                similarity >= 99.995

            ):

                similarity = 99.99


        similarity_matrix[
            i,
            j
        ] = similarity


        similarity_matrix[
            j,
            i
        ] = similarity


        # ----------------------------------------------------
        # LONG-FORMAT RESULTS
        # ----------------------------------------------------

        if i != j:

            pairwise_results.append({

                "Sequence_1":
                    sequence_names[i],

                "Sequence_2":
                    sequence_names[j],

                "Total_SNP_Sites":
                    n_snps,

                "Comparable_SNP_Sites":
                    comparable_positions,

                "Identical_SNP_Sites":
                    identical_sites,

                "Different_SNP_Sites":
                    different_sites,

                "SNP_Similarity_Percent":
                    similarity

            })


# ============================================================
# CREATE SIMILARITY DATAFRAME
# ============================================================

similarity_df = pd.DataFrame(

    similarity_matrix,

    index=sequence_names,

    columns=sequence_names

)


similarity_df.index.name = "Sequence"


# ============================================================
# SAVE SUPPLEMENTARY MATRIX CSV
# ============================================================

similarity_df.to_csv(

    OUTPUT_MATRIX_CSV,

    float_format="%.4f"

)


print(

    f"SNP similarity matrix saved: "

    f"{OUTPUT_MATRIX_CSV}"

)


# ============================================================
# SAVE LONG-FORMAT CSV
# ============================================================

pairwise_df = pd.DataFrame(
    pairwise_results
)


pairwise_df.to_csv(

    OUTPUT_LONG_CSV,

    index=False,

    float_format="%.4f"

)


print(

    f"Pairwise SNP similarity file saved: "

    f"{OUTPUT_LONG_CSV}"

)


# ============================================================
# SUMMARY STATISTICS
# ============================================================

upper_triangle_values = similarity_matrix[

    np.triu_indices(

        n_sequences,

        k=1

    )

]


upper_triangle_values = (

    upper_triangle_values[

        ~np.isnan(

            upper_triangle_values

        )

    ]

)


if len(
    upper_triangle_values
) > 0:

    print(
        "\n"
        + "=" * 70
    )

    print(
        "SNP SIMILARITY SUMMARY"
    )

    print(
        "=" * 70
    )

    print(

        f"Minimum similarity: "

        f"{np.min(upper_triangle_values):.4f}%"

    )

    print(

        f"Maximum similarity: "

        f"{np.max(upper_triangle_values):.4f}%"

    )

    print(

        f"Mean similarity:    "

        f"{np.mean(upper_triangle_values):.4f}%"

    )

    print(

        f"Median similarity:  "

        f"{np.median(upper_triangle_values):.4f}%"

    )

    print(
        "=" * 70
    )


# ============================================================
# DARK GREEN → GREEN COLOR MAP
#
# LOW SIMILARITY
#       ↓
# Very pale green
# Pale green
# Light green
# Medium green
# Green
# Dark green
# Very dark green
#       ↑
# HIGH SIMILARITY
# ============================================================

green_colors = [

    "#F0FFF0",

    "#D9F7D6",

    "#A8E6A3",

    "#66C266",

    "#2E9E2E",

    "#147A14",

    "#004D00"

]


custom_green_cmap = (

    LinearSegmentedColormap.from_list(

        "SNP_Similarity_Green",

        green_colors,

        N=256

    )

)


# ============================================================
# DETERMINE HEATMAP COLOR RANGE
# ============================================================

min_similarity = np.nanmin(

    upper_triangle_values

)


max_similarity = np.nanmax(

    similarity_matrix

)


if (

    max_similarity
    -
    min_similarity

    <

    0.01

):

    vmin = max(

        0,

        min_similarity - 0.01

    )

else:

    vmin = min_similarity


vmax = 100.0


print(
    "\nGenerating SNP similarity heatmap..."
)


print(

    f"Color scale range: "

    f"{vmin:.4f}% "

    f"to "

    f"{vmax:.2f}%"

)


# ============================================================
# MASK REFERENCE ROW AND COLUMN
#
# FIRST STRAIN REMAINS COMPLETELY BLANK
# ============================================================

heatmap_mask = np.zeros_like(

    similarity_df,

    dtype=bool

)


# Reference row
heatmap_mask[
    REFERENCE_INDEX,
    :
] = True


# Reference column
heatmap_mask[
    :,
    REFERENCE_INDEX
] = True


# Missing values
heatmap_mask |= np.isnan(

    similarity_matrix

)


# ============================================================
# FIGURE SIZE
# ============================================================

figure_size = max(

    10,

    n_sequences * 0.22

)


# ============================================================
# CREATE FIGURE
# ============================================================

fig, ax = plt.subplots(

    figsize=(

        figure_size,

        figure_size

    )

)


# ============================================================
# HEATMAP
#
# IMPORTANT:
# linewidths = 0
# linecolor = none
# rasterized = True
#
# These prevent white separation lines.
# ============================================================

sns.heatmap(

    similarity_df,

    mask=heatmap_mask,

    cmap=custom_green_cmap,

    vmin=vmin,

    vmax=vmax,

    square=True,

    linewidths=0,

    linecolor="none",

    rasterized=True,

    xticklabels=True,

    yticklabels=True,

    cbar=True,

    cbar_kws={

        "label":
            "SNP similarity (%)",

        "shrink":
            1.0,

        "aspect":
            30,

        "pad":
            0.02

    },

    ax=ax

)


# ============================================================
# FORCE ALL GRIDLINES OFF
# ============================================================

ax.grid(

    False

)

ax.xaxis.grid(

    False,

    which="both"

)

ax.yaxis.grid(

    False,

    which="both"

)


# ============================================================
# REMOVE AXES SPINES
# ============================================================

for spine in ax.spines.values():

    spine.set_visible(False)


# ============================================================
# REMOVE TICK MARKS
# ============================================================

ax.tick_params(

    axis="both",

    which="both",

    length=0

)


# ============================================================
# AXIS LABELS
# ============================================================

plt.xticks(

    rotation=90,

    fontsize=LABEL_FONT_SIZE

)

plt.yticks(

    rotation=0,

    fontsize=LABEL_FONT_SIZE

)


plt.xlabel(

    "Genome sequences",

    fontsize=10

)


plt.ylabel(

    "Genome sequences",

    fontsize=10

)


# ============================================================
# TITLE
# ============================================================

plt.title(

    f"Pairwise SNP Similarity Heatmap\n"

    f"Based on {n_snps:,} variable nucleotide positions",

    fontsize=12,

    pad=15

)


# ============================================================
# COLORBAR
# ============================================================

cbar = ax.collections[0].colorbar


# Remove colorbar border
cbar.outline.set_visible(False)


# Remove colorbar tick marks
cbar.ax.tick_params(

    length=0,

    labelsize=8

)


cbar.set_label(

    "SNP similarity (%)",

    fontsize=9

)


# ============================================================
# FORCE COLORBAR TO MATCH ONLY HEATMAP HEIGHT
# ============================================================

# Draw the figure once so Matplotlib has finalized
# the axes positions.

fig.canvas.draw()


# Get heatmap position
heatmap_position = ax.get_position()


# Get colorbar position
cbar_position = cbar.ax.get_position()


# ------------------------------------------------------------
# IMPORTANT:
#
# Colorbar bottom = heatmap bottom
# Colorbar top    = heatmap top
#
# Therefore the colorbar height is exactly the
# heatmap matrix height.
# ------------------------------------------------------------

cbar.ax.set_position([

    cbar_position.x0,

    heatmap_position.y0,

    cbar_position.width,

    heatmap_position.height

])


# ============================================================
# HEATMAP-ONLY BORDER
#
# The border surrounds ONLY the colored heatmap matrix.
# It does NOT surround:
# - title
# - axis labels
# - sequence names
# - colorbar
# ============================================================

heatmap_border = plt.Rectangle(

    (0, 0),

    n_sequences,

    n_sequences,

    fill=False,

    edgecolor="black",

    linewidth=0.8,

    zorder=20,

    clip_on=False

)


ax.add_patch(

    heatmap_border

)


# ============================================================
# DO NOT USE tight_layout() HERE
#
# tight_layout() can change the heatmap/colorbar
# positions after we matched their heights.
# ============================================================

# Instead, use controlled subplot spacing.

fig.subplots_adjust(

    left=0.20,

    right=0.88,

    bottom=0.20,

    top=0.90

)


# ============================================================
# RE-MATCH COLORBAR AFTER FINAL FIGURE POSITIONING
# ============================================================

fig.canvas.draw()


heatmap_position = ax.get_position()

cbar_position = cbar.ax.get_position()


cbar.ax.set_position([

    cbar_position.x0,

    heatmap_position.y0,

    cbar_position.width,

    heatmap_position.height

])


# ============================================================
# SAVE HIGH-RESOLUTION PNG
# ============================================================

plt.savefig(

    OUTPUT_HEATMAP_PNG,

    dpi=FIGURE_DPI,

    bbox_inches="tight",

    pad_inches=0.05

)


# ============================================================
# SAVE VECTOR PDF
# ============================================================

plt.savefig(

    OUTPUT_HEATMAP_PDF,

    bbox_inches="tight",

    pad_inches=0.05

)


plt.close()


# ============================================================
# FINAL MESSAGE
# ============================================================

print(
    "\n"
    + "=" * 70
)

print(
    "ANALYSIS COMPLETED SUCCESSFULLY"
)

print(
    "=" * 70
)

print(

    f"Alignment file:                 "
    f"{ALIGNMENT_FILE}"

)

print(

    f"Number of sequences:            "
    f"{n_sequences}"

)

print(

    f"Alignment length:               "
    f"{alignment_length:,} bp"

)

print(

    f"Reference / first strain:       "
    f"{reference_name}"

)

print(

    f"Total SNP sites:                "
    f"{n_snps:,}"

)

print(
    "\nOutput files:"
)

print(

    f"1. {OUTPUT_MATRIX_CSV}"

)

print(

    f"2. {OUTPUT_LONG_CSV}"

)

print(

    f"3. {OUTPUT_SNP_POSITIONS_CSV}"

)

print(

    f"4. {OUTPUT_HEATMAP_PNG}"

)

print(

    f"5. {OUTPUT_HEATMAP_PDF}"

)

print(
    "=" * 70
)
