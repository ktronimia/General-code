#!/usr/bin/env python3

"""
================================================================
LSDV PAIRWISE SEQUENCE IDENTITY HEATMAP + CSV
WITH TOP PHYLOGENETIC TREE
================================================================

INPUT
-----
1. Multiple sequence alignment
2. Newick phylogenetic tree

OUTPUT
------
1. LSDV_pairwise_identity_heatmap.png
2. LSDV_pairwise_identity_heatmap.pdf
3. LSDV_pairwise_identity_matrix.csv
4. LSDV_pairwise_identity_long_format.csv

COLOR SCALE
-----------
Lowest identity  -> Slightly sky blue
                   Light blue
                   Blue (~98%)
                   Deep blue
99.99% identity   -> Dark blue
100.00% identity  -> Black

100.00% is reserved ONLY for truly identical sequences.
================================================================
"""

# ================================================================
# USER CONFIGURATION
# ================================================================

alignment_path = "alignment.fasta"
tree_path = "phylogenetic_tree.nwk"

output_prefix = "LSDV_pairwise_identity"

dpi = 600

font_size = 12
leaf_font_size = 6.5

ignore_double_gaps = True

# Annotate cells only if number of genomes <= this value
annotate_max_sequences = 40


# ================================================================
# IMPORTS
# ================================================================

import os
import gc

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from matplotlib import gridspec

from matplotlib.colors import (
    LinearSegmentedColormap,
    Normalize
)

from Bio import AlignIO
from Bio import Phylo


# ================================================================
# MATPLOTLIB SETTINGS
# ================================================================

plt.rcParams.update({

    "font.size": font_size,

    "axes.titlesize": font_size,

    "axes.labelsize": font_size,

    "xtick.labelsize": leaf_font_size,

    "ytick.labelsize": leaf_font_size,

    "legend.fontsize": font_size,

    "pdf.fonttype": 42,

    "ps.fonttype": 42

})


# ================================================================
# READ ALIGNMENT
# ================================================================

def read_alignment(file_path):

    formats_to_try = [

        "fasta",
        "clustal",
        "phylip-relaxed",
        "phylip",
        "nexus"

    ]

    for fmt in formats_to_try:

        try:

            alignment = AlignIO.read(
                file_path,
                fmt
            )

            if len(alignment) < 2:

                raise ValueError(
                    "At least two sequences are required."
                )

            print(
                f"Alignment successfully read using: {fmt}"
            )

            print(
                f"Number of sequences: {len(alignment)}"
            )

            print(
                f"Alignment length: "
                f"{alignment.get_alignment_length()} bp"
            )

            return alignment

        except Exception:

            continue

    raise ValueError(
        "\nCould not read alignment file:\n"
        f"{file_path}"
    )


# ================================================================
# CALCULATE PAIRWISE IDENTITY
# ================================================================

def calculate_pairwise_identity(
    alignment,
    ignore_double_gaps=True
):

    """
    Calculate pairwise sequence identity.

    100.00%:
        Only completely identical sequences.

    Non-identical sequences:
        Maximum = 99.99%
    """

    n_sequences = len(alignment)

    print(
        "\nConverting alignment to NumPy arrays..."
    )

    sequences = [

        np.frombuffer(
            str(record.seq).upper().encode("ascii"),
            dtype=np.uint8
        )

        for record in alignment

    ]

    matrix = np.zeros(
        (n_sequences, n_sequences),
        dtype=np.float64
    )

    print(
        "\nCalculating pairwise sequence identity..."
    )

    for i in range(n_sequences):

        print(
            f"Processing sequence "
            f"{i + 1}/{n_sequences}",
            end="\r",
            flush=True
        )

        seq1 = sequences[i]

        # --------------------------------------------------------
        # Self identity
        # --------------------------------------------------------

        matrix[i, i] = 100.00

        for j in range(i + 1, n_sequences):

            seq2 = sequences[j]

            # ----------------------------------------------------
            # Ignore double gaps
            # ----------------------------------------------------

            if ignore_double_gaps:

                valid_mask = ~(
                    (seq1 == ord("-")) &
                    (seq2 == ord("-"))
                )

            else:

                valid_mask = np.ones(
                    len(seq1),
                    dtype=bool
                )

            valid_positions = np.count_nonzero(
                valid_mask
            )

            # ----------------------------------------------------
            # No valid positions
            # ----------------------------------------------------

            if valid_positions == 0:

                identity_percent = 0.0

            else:

                matches = np.count_nonzero(

                    (seq1 == seq2) &
                    valid_mask

                )

                identity_fraction = (
                    matches /
                    valid_positions
                )

                identity_percent = (
                    identity_fraction * 100.0
                )

                # ------------------------------------------------
                # Exact identity
                # ------------------------------------------------

                if identity_fraction == 1.0:

                    identity_percent = 100.00

                else:

                    # ------------------------------------------------
                    # Non-identical sequences cannot become 100%
                    # ------------------------------------------------

                    identity_percent = min(
                        identity_percent,
                        99.99
                    )

                    identity_percent = round(
                        identity_percent,
                        2
                    )

                    if identity_percent >= 100.00:

                        identity_percent = 99.99

            matrix[i, j] = identity_percent

            matrix[j, i] = identity_percent

    print(
        "\nPairwise identity calculation completed."
    )

    return matrix


# ================================================================
# READ PHYLOGENETIC TREE
# ================================================================

def read_and_prepare_tree(
    tree_file,
    sequence_ids
):

    print(
        "\nReading phylogenetic tree..."
    )

    tree = Phylo.read(
        tree_file,
        "newick"
    )

    alignment_ids = set(sequence_ids)

    terminals_to_remove = [

        terminal

        for terminal in tree.get_terminals()

        if terminal.name not in alignment_ids

    ]

    for terminal in terminals_to_remove:

        tree.prune(terminal)

    tree_ids = set(

        terminal.name

        for terminal in tree.get_terminals()

    )

    missing_in_tree = (
        alignment_ids - tree_ids
    )

    if missing_in_tree:

        print("\nWARNING:")

        print(
            f"{len(missing_in_tree)} alignment sequence(s) "
            "were not found in the phylogenetic tree."
        )

    return tree


# ================================================================
# GET TREE ORDER
# ================================================================

def get_tree_order(
    tree,
    sequence_ids
):

    alignment_ids = set(sequence_ids)

    ordered_ids = [

        terminal.name

        for terminal in tree.get_terminals()

        if terminal.name in alignment_ids

    ]

    return ordered_ids


# ================================================================
# REORDER MATRIX
# ================================================================

def reorder_matrix(
    matrix,
    sequence_ids,
    ordered_ids
):

    index_map = {

        sequence_id: index

        for index, sequence_id
        in enumerate(sequence_ids)

    }

    indices = [

        index_map[sequence_id]

        for sequence_id in ordered_ids

    ]

    return matrix[
        np.ix_(
            indices,
            indices
        )
    ]


# ================================================================
# CREATE BLUE -> BLACK COLORMAP
# ================================================================

def create_blue_black_colormap():

    """
    Color scheme:

    Lowest identity
        -> Slightly sky blue

    Intermediate
        -> Light blue

    ~98%
        -> Blue

    Higher identity
        -> Deep blue

    99.99%
        -> Dark blue

    100%
        -> BLACK
    """

    colors = [

        "#BFEFFF",   # sky blue
        "#87CEEB",   # light sky blue
        "#5BC0EB",   # light blue
        "#2196F3",   # blue
        "#1976D2",   # medium blue
        "#1565C0",   # deep blue
        "#0D47A1",   # dark blue
        "#061A40",   # very dark blue
        "#000000"    # black = 100%
    ]

    positions = [

        0.00,
        0.10,
        0.20,
        0.50,
        0.65,
        0.78,
        0.90,
        0.985,
        1.00

    ]

    return LinearSegmentedColormap.from_list(

        "LSDV_Blue_Black",

        list(
            zip(
                positions,
                colors
            )
        ),

        N=512

    )


# ================================================================
# CUSTOM NORMALIZATION
# ================================================================

class LSDVIdentityNormalize(Normalize):

    """
    Nonlinear normalization designed to clearly separate:

    ~98%
    99%
    99.5%
    99.9%
    99.99%
    100%

    99.99% remains dark blue.

    100% becomes black.
    """

    def __init__(
        self,
        vmin=90.0,
        vmax=100.0,
        clip=False
    ):

        super().__init__(
            vmin=vmin,
            vmax=vmax,
            clip=clip
        )

    def __call__(
        self,
        value,
        clip=None
    ):

        value = np.asarray(
            value,
            dtype=float
        )

        result = np.zeros_like(
            value,
            dtype=float
        )

        # --------------------------------------------------------
        # Lowest -> 97%
        # Sky blue region
        # --------------------------------------------------------

        mask1 = value < 97.0

        if self.vmin < 97.0:

            result[mask1] = (

                (
                    value[mask1] -
                    self.vmin
                )
                /
                (
                    97.0 -
                    self.vmin
                )

            ) * 0.20

        else:

            result[mask1] = 0.0

        # --------------------------------------------------------
        # 97 -> 98%
        # Sky blue -> blue
        # --------------------------------------------------------

        mask2 = (

            (value >= 97.0) &
            (value < 98.0)

        )

        result[mask2] = (

            0.20 +

            (
                (
                    value[mask2] -
                    97.0
                ) /
                1.0
            ) * 0.35

        )

        # --------------------------------------------------------
        # 98 -> 99%
        # Blue -> deep blue
        # --------------------------------------------------------

        mask3 = (

            (value >= 98.0) &
            (value < 99.0)

        )

        result[mask3] = (

            0.55 +

            (
                (
                    value[mask3] -
                    98.0
                ) /
                1.0
            ) * 0.20

        )

        # --------------------------------------------------------
        # 99 -> 99.5%
        # Deep blue
        # --------------------------------------------------------

        mask4 = (

            (value >= 99.0) &
            (value < 99.5)

        )

        result[mask4] = (

            0.75 +

            (
                (
                    value[mask4] -
                    99.0
                ) /
                0.5
            ) * 0.10

        )

        # --------------------------------------------------------
        # 99.5 -> 99.9%
        # Dark blue
        # --------------------------------------------------------

        mask5 = (

            (value >= 99.5) &
            (value < 99.9)

        )

        result[mask5] = (

            0.85 +

            (
                (
                    value[mask5] -
                    99.5
                ) /
                0.4
            ) * 0.10

        )

        # --------------------------------------------------------
        # 99.9 -> 99.99%
        # Very dark blue
        # --------------------------------------------------------

        mask6 = (

            (value >= 99.9) &
            (value < 99.99)

        )

        result[mask6] = (

            0.95 +

            (
                (
                    value[mask6] -
                    99.9
                ) /
                0.09
            ) * 0.035

        )

        # --------------------------------------------------------
        # 99.99 -> 100%
        #
        # Dark blue -> BLACK
        #
        # This creates a visible distinction between
        # 99.99% and 100%.
        # --------------------------------------------------------

        mask7 = (

            (value >= 99.99) &
            (value <= 100.0)

        )

        result[mask7] = (

            0.985 +

            (
                (
                    value[mask7] -
                    99.99
                ) /
                0.01
            ) * 0.015

        )

        return np.ma.masked_array(

            np.clip(
                result,
                0.0,
                1.0
            )

        )


# ================================================================
# PREPARE TREE COORDINATES
# ================================================================

def prepare_tree_coordinates(tree):

    depths = tree.depths()

    terminals = tree.get_terminals()

    max_depth = max(

        depths.get(
            terminal,
            0.0
        )

        for terminal in terminals

    )

    if max_depth == 0.0:

        depths = tree.depths(
            unit_branch_lengths=True
        )

        max_depth = max(

            depths.get(
                terminal,
                0.0
            )

            for terminal in terminals

        )

    terminal_positions = {

        terminal: index

        for index, terminal
        in enumerate(terminals)

    }

    positions = {}

    def assign_position(clade):

        if clade.is_terminal():

            position = (
                terminal_positions[clade]
            )

        else:

            child_positions = [

                assign_position(child)

                for child in clade.clades

            ]

            position = np.mean(
                child_positions
            )

        positions[clade] = position

        return position

    assign_position(
        tree.root
    )

    return (
        depths,
        positions,
        max_depth
    )


# ================================================================
# DRAW TOP TREE
# ================================================================

def draw_top_tree(
    tree,
    ax,
    depths,
    positions,
    max_depth,
    n
):

    def draw_clade(clade):

        y_parent = depths.get(
            clade,
            0.0
        )

        children = clade.clades

        if children:

            child_x = [

                positions[child]

                for child in children

            ]

            ax.plot(

                [
                    min(child_x),
                    max(child_x)
                ],

                [
                    y_parent,
                    y_parent
                ],

                color="black",

                linewidth=0.6

            )

            for child in children:

                y_child = depths.get(
                    child,
                    y_parent
                )

                x_child = positions[child]

                ax.plot(

                    [
                        x_child,
                        x_child
                    ],

                    [
                        y_parent,
                        y_child
                    ],

                    color="black",

                    linewidth=0.6

                )

                draw_clade(
                    child
                )

    draw_clade(
        tree.root
    )

    ax.set_xlim(
        -0.5,
        n - 0.5
    )

    ax.set_ylim(
        max_depth * 1.08,
        0
    )

    ax.set_xticks([])
    ax.set_yticks([])

    ax.axis("off")


# ================================================================
# EXPORT PAIRWISE IDENTITY CSV
# ================================================================

def export_identity_csv(
    matrix,
    ordered_ids,
    output_prefix
):

    print(
        "\nExporting pairwise identity CSV..."
    )

    # ============================================================
    # MATRIX FORMAT
    # ============================================================

    matrix_df = pd.DataFrame(

        matrix,

        index=ordered_ids,

        columns=ordered_ids

    )

    matrix_df.index.name = (
        "Genome"
    )

    matrix_csv = (
        f"{output_prefix}_matrix.csv"
    )

    matrix_df.to_csv(
        matrix_csv,
        float_format="%.2f"
    )

    print(
        "Matrix CSV saved:"
    )

    print(
        os.path.abspath(
            matrix_csv
        )
    )

    # ============================================================
    # LONG FORMAT
    # ============================================================

    print(
        "\nCreating long-format supplementary table..."
    )

    rows = []

    n = len(ordered_ids)

    for i in range(n):

        for j in range(i + 1, n):

            rows.append({

                "Genome_1":
                    ordered_ids[i],

                "Genome_2":
                    ordered_ids[j],

                "Pairwise_identity_percent":
                    round(
                        float(
                            matrix[i, j]
                        ),
                        2
                    )

            })

    long_df = pd.DataFrame(
        rows
    )

    long_csv = (
        f"{output_prefix}_long_format.csv"
    )

    long_df.to_csv(
        long_csv,
        index=False,
        float_format="%.2f"
    )

    print(
        "Long-format CSV saved:"
    )

    print(
        os.path.abspath(
            long_csv
        )
    )

    return matrix_df, long_df


# ================================================================
# CREATE FIGURE
# ================================================================

def create_figure(
    matrix,
    ordered_ids,
    tree,
    output_prefix,
    dpi=600
):

    n = len(
        ordered_ids
    )

    print(
        f"\nPreparing figure for {n} genomes..."
    )

    # ============================================================
    # COLOR MAP
    # ============================================================

    cmap = (
        create_blue_black_colormap()
    )

    # ============================================================
    # COLOR RANGE
    # ============================================================

    matrix_min = float(
        np.min(matrix)
    )

    vmin = min(
        matrix_min,
        97.0
    )

    norm = (
        LSDVIdentityNormalize(
            vmin=vmin,
            vmax=100.0
        )
    )

    # ============================================================
    # FIGURE
    # ============================================================

    fig = plt.figure(

        figsize=(
            18,
            20
        ),

        facecolor="white"

    )

    gs = gridspec.GridSpec(

        2,
        2,

        figure=fig,

        width_ratios=[
            1.0,
            0.04
        ],

        height_ratios=[
            0.22,
            1.0
        ],

        wspace=0.05,

        hspace=0.03

    )

    ax_top = fig.add_subplot(
        gs[0, 0]
    )

    ax_heatmap = fig.add_subplot(
        gs[1, 0]
    )

    ax_colorbar = fig.add_subplot(
        gs[1, 1]
    )

    # ============================================================
    # TREE
    # ============================================================

    print(
        "Preparing phylogenetic tree..."
    )

    depths, positions, max_depth = (
        prepare_tree_coordinates(
            tree
        )
    )

    draw_top_tree(

        tree,

        ax_top,

        depths,

        positions,

        max_depth,

        n

    )

    # ============================================================
    # HEATMAP
    # ============================================================

    print(
        "Drawing heatmap..."
    )

    image = ax_heatmap.imshow(

        matrix,

        cmap=cmap,

        norm=norm,

        interpolation="nearest",

        aspect="equal",

        rasterized=True

    )

    ax_heatmap.set_xlim(
        -0.5,
        n - 0.5
    )

    ax_heatmap.set_ylim(
        n - 0.5,
        -0.5
    )

    ax_heatmap.set_aspect(
        "equal",
        adjustable="box"
    )

    # ============================================================
    # LABELS
    # ============================================================

    ax_heatmap.set_xticks(
        np.arange(n)
    )

    ax_heatmap.set_yticks(
        np.arange(n)
    )

    ax_heatmap.set_xticklabels(

        ordered_ids,

        rotation=90,

        ha="center",

        va="top",

        fontsize=leaf_font_size

    )

    ax_heatmap.set_yticklabels(

        ordered_ids,

        fontsize=leaf_font_size

    )

    ax_heatmap.tick_params(

        axis="both",

        which="major",

        length=0,

        pad=3

    )

    # ============================================================
    # AXIS LABELS
    # ============================================================

    ax_heatmap.set_xlabel(

        "LSDV complete genomes",

        fontsize=font_size,

        labelpad=15

    )

    ax_heatmap.set_ylabel(

        "LSDV complete genomes",

        fontsize=font_size,

        labelpad=15

    )

    # ============================================================
    # BORDER
    # ============================================================

    for spine in (
        ax_heatmap.spines.values()
    ):

        spine.set_visible(True)

        spine.set_linewidth(
            0.8
        )

    # ============================================================
    # ANNOTATIONS
    # ============================================================

    if n <= annotate_max_sequences:

        print(
            "Adding identity annotations..."
        )

        for i in range(n):

            for j in range(n):

                value = (
                    matrix[i, j]
                )

                if value == 100.00:

                    label = "100"

                else:

                    label = (
                        f"{value:.2f}"
                    )

                ax_heatmap.text(

                    j,

                    i,

                    label,

                    ha="center",

                    va="center",

                    fontsize=6

                )

    # ============================================================
    # COLORBAR
    # ============================================================

    cbar = fig.colorbar(

        image,

        cax=ax_colorbar

    )

    cbar.set_label(

        "Pairwise sequence identity (%)",

        fontsize=font_size,

        labelpad=10

    )

    cbar.ax.tick_params(

        labelsize=font_size

    )

    # ------------------------------------------------------------
    # IMPORTANT TICKS
    #
    # 99.99 and 100 are deliberately included
    # ------------------------------------------------------------

    colorbar_ticks = [

        97.00,

        97.50,

        98.00,

        98.50,

        99.00,

        99.50,

        99.90,

        99.95,

        99.99,

        100.00

    ]

    if vmin < 97.0:

        colorbar_ticks.insert(

            0,

            round(
                vmin,
                2
            )

        )

    colorbar_ticks = sorted(
        set(
            tick

            for tick in colorbar_ticks

            if vmin <= tick <= 100.0
        )
    )

    cbar.set_ticks(
        colorbar_ticks
    )

    cbar.set_ticklabels(

        [

            f"{tick:.2f}"

            for tick
            in colorbar_ticks

        ]

    )

    # ============================================================
    # TITLE
    # ============================================================

    fig.suptitle(

        "Pairwise Genome Sequence Identity and "
        "Phylogenetic Relationships of LSDV Complete Genomes",

        fontsize=font_size,

        fontweight="bold",

        y=0.985

    )

    # ============================================================
    # SAVE PNG
    # ============================================================

    png_file = (

        f"{output_prefix}.png"

    )

    print(
        f"\nSaving PNG at {dpi} DPI..."
    )

    fig.savefig(

        png_file,

        dpi=dpi,

        facecolor="white",

        bbox_inches="tight"

    )

    print(
        "PNG saved:"
    )

    print(
        os.path.abspath(
            png_file
        )
    )

    # ============================================================
    # SAVE PDF
    # ============================================================

    pdf_file = (

        f"{output_prefix}.pdf"

    )

    print(
        "\nSaving vector PDF..."
    )

    fig.savefig(

        pdf_file,

        facecolor="white",

        bbox_inches="tight"

    )

    print(
        "PDF saved:"
    )

    print(
        os.path.abspath(
            pdf_file
        )
    )

    plt.close(
        fig
    )

    gc.collect()


# ================================================================
# MAIN
# ================================================================

def main():

    print(
        "=" * 70
    )

    print(
        "LSDV PAIRWISE SEQUENCE IDENTITY ANALYSIS"
    )

    print(
        "=" * 70
    )

    # ============================================================
    # CHECK INPUTS
    # ============================================================

    if not os.path.isfile(
        alignment_path
    ):

        raise FileNotFoundError(

            "\nAlignment file not found:\n"
            f"{alignment_path}"

        )

    if not os.path.isfile(
        tree_path
    ):

        raise FileNotFoundError(

            "\nPhylogenetic tree not found:\n"
            f"{tree_path}"

        )

    # ============================================================
    # READ ALIGNMENT
    # ============================================================

    print(
        "\nSTEP 1: Reading alignment..."
    )

    alignment = read_alignment(
        alignment_path
    )

    sequence_ids = [

        record.id

        for record in alignment

    ]

    # ============================================================
    # DUPLICATE CHECK
    # ============================================================

    if len(sequence_ids) != len(
        set(sequence_ids)
    ):

        raise ValueError(

            "\nDuplicate sequence IDs detected.\n"
            "Every sequence must have a unique ID."

        )

    # ============================================================
    # CALCULATE IDENTITY
    # ============================================================

    print(
        "\nSTEP 2: Calculating pairwise identity..."
    )

    identity_matrix = (
        calculate_pairwise_identity(

            alignment,

            ignore_double_gaps=
                ignore_double_gaps

        )
    )

    del alignment

    gc.collect()

    # ============================================================
    # READ TREE
    # ============================================================

    print(
        "\nSTEP 3: Reading phylogenetic tree..."
    )

    tree = (
        read_and_prepare_tree(

            tree_path,

            sequence_ids

        )
    )

    # ============================================================
    # TREE ORDER
    # ============================================================

    print(
        "\nSTEP 4: Determining phylogenetic order..."
    )

    ordered_ids = (
        get_tree_order(

            tree,

            sequence_ids

        )
    )

    if len(ordered_ids) < 2:

        raise ValueError(

            "\nFewer than two sequences are shared "
            "between alignment and tree."

        )

    print(

        "Number of genomes included: "
        f"{len(ordered_ids)}"

    )

    # ============================================================
    # REORDER MATRIX
    # ============================================================

    print(
        "\nSTEP 5: Reordering identity matrix..."
    )

    ordered_matrix = (
        reorder_matrix(

            identity_matrix,

            sequence_ids,

            ordered_ids

        )
    )

    del identity_matrix

    gc.collect()

    # ============================================================
    # EXPORT CSV
    # ============================================================

    print(
        "\nSTEP 6: Exporting pairwise identity tables..."
    )

    export_identity_csv(

        ordered_matrix,

        ordered_ids,

        output_prefix

    )

    # ============================================================
    # CREATE FIGURE
    # ============================================================

    print(
        "\nSTEP 7: Creating publication-ready heatmap..."
    )

    create_figure(

        matrix=ordered_matrix,

        ordered_ids=ordered_ids,

        tree=tree,

        output_prefix=output_prefix,

        dpi=dpi

    )

    # ============================================================
    # FINISHED
    # ============================================================

    print(
        "\n" + "=" * 70
    )

    print(
        "ANALYSIS COMPLETED SUCCESSFULLY"
    )

    print(
        "=" * 70
    )

    print(
        "\nGenerated files:"
    )

    print(
        f"  {output_prefix}.png"
    )

    print(
        f"  {output_prefix}.pdf"
    )

    print(
        f"  {output_prefix}_matrix.csv"
    )

    print(
        f"  {output_prefix}_long_format.csv"
    )


# ================================================================
# RUN
# ================================================================

if __name__ == "__main__":

    main()
