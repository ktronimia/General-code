#!/usr/bin/env python3

"""
===============================================================
LSDV PAIRWISE SEQUENCE IDENTITY HEATMAP
WITH TOP PHYLOGENETIC TREE
===============================================================

Input:
    1. Multiple sequence alignment
    2. Newick phylogenetic tree

Output:
    1. PNG - 600 DPI
    2. PDF - vector format

COLOR SCHEME:
    Lowest identity  -> Slightly sky blue
    Intermediate     -> Light blue
    ~98.00%          -> Blue
    99.00%           -> Deep blue
    99.50%           -> Dark blue
    99.90%           -> Very dark blue
    99.99%           -> Extremely dark blue
    100.00%          -> Darkest blue

IMPORTANT:
    99.99% and 100.00% are intentionally separated
    to make highly similar genomes visually distinguishable.

Features:
    - Pairwise genome identity calculated directly from alignment
    - Exact identical sequences = 100.00%
    - Non-identical sequences capped at 99.99%
    - Top phylogenetic tree
    - No left-side phylogenetic tree
    - Square heatmap
    - Leaf labels = 6.5 pt
    - Publication-ready PNG and PDF
===============================================================
"""


# =============================================================
# USER CONFIGURATION
# =============================================================

alignment_path = "alignment.fasta"
tree_path = "phylogenetic_tree.nwk"
output_path = "LSDV_pairwise_identity_heatmap"

dpi = 600
font_size = 12
leaf_font_size = 6.5

# Ignore positions where BOTH sequences contain gaps
ignore_double_gaps = True

# Annotate heatmap cells only when number of genomes <= this
annotate_max_sequences = 40


# =============================================================
# IMPORTS
# =============================================================

import os
import gc

import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap, Normalize

from Bio import AlignIO
from Bio import Phylo


# =============================================================
# GLOBAL MATPLOTLIB SETTINGS
# =============================================================

plt.rcParams.update({

    "font.size": font_size,

    "axes.titlesize": font_size,

    "axes.labelsize": font_size,

    "xtick.labelsize": leaf_font_size,

    "ytick.labelsize": leaf_font_size,

    "legend.fontsize": font_size,

    # Better font compatibility in PDF
    "pdf.fonttype": 42,
    "ps.fonttype": 42,

})


# =============================================================
# READ ALIGNMENT
# =============================================================

def read_alignment(file_path):

    """
    Read a multiple sequence alignment.

    Supported formats:
        FASTA
        CLUSTAL
        PHYLIP
        NEXUS
    """

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
                f"Alignment successfully read using format: {fmt}"
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


# =============================================================
# PAIRWISE SEQUENCE IDENTITY
# =============================================================

def calculate_pairwise_identity(
    alignment,
    ignore_double_gaps=True
):

    """
    Calculate pairwise sequence identity.

    Exact identical sequences:
        100.00%

    Non-identical sequences:
        Maximum displayed value = 99.99%
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
        dtype=np.float32
    )

    print(
        "\nCalculating pairwise sequence identity..."
    )

    for i in range(n_sequences):

        print(
            f"Processing sequence {i + 1}/{n_sequences}",
            end="\r",
            flush=True
        )

        seq1 = sequences[i]

        # Exact self identity
        matrix[i, i] = 100.00

        for j in range(i + 1, n_sequences):

            seq2 = sequences[j]

            # -------------------------------------------------
            # Ignore positions where BOTH sequences are gaps
            # -------------------------------------------------

            if ignore_double_gaps:

                valid_mask = ~(
                    (seq1 == ord("-"))
                    &
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

            # -------------------------------------------------
            # No valid positions
            # -------------------------------------------------

            if valid_positions == 0:

                identity_percent = 0.0

            else:

                matches = np.count_nonzero(
                    (seq1 == seq2)
                    &
                    valid_mask
                )

                identity_fraction = (
                    matches /
                    valid_positions
                )

                # -------------------------------------------------
                # Exact identical sequences
                # -------------------------------------------------

                if identity_fraction == 1.0:

                    identity_percent = 100.00

                else:

                    identity_percent = (
                        identity_fraction *
                        100.0
                    )

                    # -------------------------------------------------
                    # Non-identical sequences cannot be 100%
                    # -------------------------------------------------

                    identity_percent = min(
                        identity_percent,
                        99.99
                    )

                    identity_percent = round(
                        identity_percent,
                        2
                    )

                    # Final protection
                    if identity_percent >= 100.00:

                        identity_percent = 99.99

            matrix[i, j] = identity_percent

            matrix[j, i] = identity_percent

    print(
        "\nPairwise identity calculation completed."
    )

    return matrix


# =============================================================
# READ PHYLOGENETIC TREE
# =============================================================

def read_and_prepare_tree(
    tree_file,
    sequence_ids
):

    """
    Read Newick tree and retain only leaves
    present in the alignment.
    """

    print(
        "\nReading phylogenetic tree..."
    )

    tree = Phylo.read(
        tree_file,
        "newick"
    )

    alignment_ids = set(
        sequence_ids
    )

    terminals_to_remove = [

        terminal

        for terminal in tree.get_terminals()

        if terminal.name not in alignment_ids

    ]

    for terminal in terminals_to_remove:

        tree.prune(
            terminal
        )

    tree_ids = set(

        terminal.name

        for terminal in tree.get_terminals()

    )

    missing_in_tree = (
        alignment_ids -
        tree_ids
    )

    if missing_in_tree:

        print(
            "\nWARNING:"
        )

        print(
            f"{len(missing_in_tree)} alignment sequence(s) "
            "were not found in the phylogenetic tree."
        )

    return tree


# =============================================================
# GET TREE ORDER
# =============================================================

def get_tree_order(
    tree,
    sequence_ids
):

    """
    Get sequence order according to
    phylogenetic tree topology.
    """

    alignment_ids = set(
        sequence_ids
    )

    ordered_ids = [

        terminal.name

        for terminal in tree.get_terminals()

        if terminal.name in alignment_ids

    ]

    return ordered_ids


# =============================================================
# REORDER MATRIX
# =============================================================

def reorder_matrix(
    matrix,
    sequence_ids,
    ordered_ids
):

    """
    Reorder identity matrix according to
    phylogenetic tree leaf order.
    """

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


# =============================================================
# CUSTOM SKY BLUE -> BLUE -> DARK BLUE COLORMAP
# =============================================================

def create_blue_colormap():

    """
    Identity color scheme:

    Lowest identity
        -> Slightly sky blue

    Intermediate identity
        -> Light blue

    ~98%
        -> Blue

    99.00%
        -> Deep blue

    99.50%
        -> Dark blue

    99.90%
        -> Very dark blue

    99.99%
        -> Extremely dark blue

    EXACT 100%
        -> Darkest blue
    """

    colors = [

        # -----------------------------------------------------
        # Lowest identity
        # -----------------------------------------------------

        "#BFEFFF",

        # -----------------------------------------------------
        # Sky blue
        # -----------------------------------------------------

        "#87CEEB",

        "#5BC0EB",

        # -----------------------------------------------------
        # Light blue
        # -----------------------------------------------------

        "#3FA9F5",

        "#2196F3",

        # -----------------------------------------------------
        # ~98% identity
        # -----------------------------------------------------

        "#1976D2",

        # -----------------------------------------------------
        # Higher identities
        # -----------------------------------------------------

        "#1565C0",

        "#0D47A1",

        "#083B8A",

        # -----------------------------------------------------
        # 99.90%
        # -----------------------------------------------------

        "#062F73",

        # -----------------------------------------------------
        # 99.99%
        # -----------------------------------------------------

        "#041F4A",

        # -----------------------------------------------------
        # EXACT 100%
        # -----------------------------------------------------

        "#000000"

    ]

    positions = [

        0.00,

        0.10,

        0.20,

        0.35,

        0.50,

        # ~98%
        0.65,

        # Higher identities
        0.73,

        0.80,

        0.86,

        # 99.90%
        0.91,

        # 99.99%
        0.96,

        # 100%
        1.00

    ]

    cmap = LinearSegmentedColormap.from_list(

        "LSDV_SkyBlue_Blue_DarkBlue",

        list(
            zip(
                positions,
                colors
            )
        ),

        N=512

    )

    return cmap


# =============================================================
# CUSTOM NONLINEAR COLOR NORMALIZATION
# =============================================================

class LSDVIdentityNormalize(Normalize):

    """
    Nonlinear identity normalization.

    Color relationship:

        Lowest identity
            -> Sky blue

        97%
            -> Light/medium blue

        98%
            -> Blue

        99%
            -> Deep blue

        99.50%
            -> Dark blue

        99.90%
            -> Very dark blue

        99.99%
            -> Extremely dark blue

        100.00%
            -> Darkest blue

    The upper identity range is deliberately expanded
    so that 99.99% and 100.00% remain distinguishable.
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

        # =====================================================
        # BELOW 97%
        # SKY BLUE
        # =====================================================

        mask1 = value < 97.00

        if self.vmin < 97.00:

            result[mask1] = (

                (
                    value[mask1] -
                    self.vmin
                )

                /

                (
                    97.00 -
                    self.vmin
                )

            ) * 0.20

        else:

            result[mask1] = 0.0


        # =====================================================
        # 97.00–98.00%
        # SKY BLUE -> BLUE
        # =====================================================

        mask2 = (

            (value >= 97.00)

            &

            (value < 98.00)

        )

        result[mask2] = (

            0.20

            +

            (

                (
                    value[mask2] -
                    97.00
                )

                /

                1.00

            ) * 0.40

        )


        # =====================================================
        # 98.00–99.00%
        # BLUE -> DEEP BLUE
        # =====================================================

        mask3 = (

            (value >= 98.00)

            &

            (value < 99.00)

        )

        result[mask3] = (

            0.60

            +

            (

                (
                    value[mask3] -
                    98.00
                )

                /

                1.00

            ) * 0.15

        )


        # =====================================================
        # 99.00–99.50%
        # DEEP BLUE
        # =====================================================

        mask4 = (

            (value >= 99.00)

            &

            (value < 99.50)

        )

        result[mask4] = (

            0.75

            +

            (

                (
                    value[mask4] -
                    99.00
                )

                /

                0.50

            ) * 0.10

        )


        # =====================================================
        # 99.50–99.90%
        # DARK BLUE
        # =====================================================

        mask5 = (

            (value >= 99.50)

            &

            (value < 99.90)

        )

        result[mask5] = (

            0.85

            +

            (

                (
                    value[mask5] -
                    99.50
                )

                /

                0.40

            ) * 0.06

        )


        # =====================================================
        # 99.90–99.99%
        # VERY DARK BLUE
        # =====================================================

        mask6 = (

            (value >= 99.90)

            &

            (value < 99.99)

        )

        result[mask6] = (

            0.91

            +

            (

                (
                    value[mask6] -
                    99.90
                )

                /

                0.09

            ) * 0.05

        )


        # =====================================================
        # 99.99%
        # EXTREMELY DARK BLUE
        # =====================================================

        mask7 = (

            (value >= 99.99)

            &

            (value < 100.00)

        )

        result[mask7] = 0.96


        # =====================================================
        # EXACT 100.00%
        # DARKEST BLUE
        # =====================================================

        mask8 = (
            value >= 100.00
        )

        result[mask8] = 1.00


        return np.ma.masked_array(

            np.clip(
                result,
                0.0,
                1.0
            )

        )


# =============================================================
# PREPARE TREE COORDINATES
# =============================================================

def prepare_tree_coordinates(tree):

    """
    Calculate coordinates for drawing the
    top phylogenetic tree.
    """

    depths = tree.depths()

    terminals = tree.get_terminals()

    max_depth = max(

        depths.get(
            terminal,
            0.0
        )

        for terminal in terminals

    )

    # ---------------------------------------------------------
    # Handle trees without branch lengths
    # ---------------------------------------------------------

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


# =============================================================
# DRAW TOP PHYLOGENETIC TREE
# =============================================================

def draw_top_tree(
    tree,
    ax,
    depths,
    positions,
    max_depth,
    n
):

    """
    Draw top phylogenetic tree.

    Terminal positions align exactly with
    heatmap columns.
    """

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

            # -------------------------------------------------
            # Horizontal internal connector
            # -------------------------------------------------

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

            # -------------------------------------------------
            # Draw branches
            # -------------------------------------------------

            for child in children:

                y_child = depths.get(
                    child,
                    y_parent
                )

                x_child = positions[
                    child
                ]

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

    # ---------------------------------------------------------
    # Match heatmap width exactly
    # ---------------------------------------------------------

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

    ax.axis(
        "off"
    )


# =============================================================
# CREATE FIGURE
# =============================================================

def create_figure(
    matrix,
    ordered_ids,
    tree,
    output_prefix,
    dpi=600
):

    """
    Create final publication-ready figure.
    """

    n = len(
        ordered_ids
    )

    print(
        f"\nPreparing figure for {n} genomes..."
    )


    # =========================================================
    # COLOR MAP
    # =========================================================

    cmap = create_blue_colormap()


    # =========================================================
    # COLOR SCALE RANGE
    # =========================================================

    matrix_min = float(
        np.min(matrix)
    )

    # Include informative 97–100% range
    vmin = min(
        matrix_min,
        97.0
    )

    norm = LSDVIdentityNormalize(
        vmin=vmin,
        vmax=100.0
    )


    # =========================================================
    # FIGURE
    # =========================================================

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


    # =========================================================
    # PREPARE TREE
    # =========================================================

    print(
        "Preparing phylogenetic tree..."
    )

    depths, positions, max_depth = (

        prepare_tree_coordinates(
            tree
        )

    )


    # =========================================================
    # DRAW TOP TREE
    # =========================================================

    print(
        "Drawing top phylogenetic tree..."
    )

    draw_top_tree(

        tree,

        ax_top,

        depths,

        positions,

        max_depth,

        n

    )


    # =========================================================
    # DRAW HEATMAP
    # =========================================================

    print(
        "Drawing square heatmap..."
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


    # =========================================================
    # TICKS AND LABELS
    # =========================================================

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


    # =========================================================
    # AXIS LABELS
    # =========================================================

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


    # =========================================================
    # HEATMAP BORDER
    # =========================================================

    for spine in ax_heatmap.spines.values():

        spine.set_visible(True)

        spine.set_linewidth(
            0.8
        )


    # =========================================================
    # OPTIONAL CELL ANNOTATIONS
    # =========================================================

    if n <= annotate_max_sequences:

        print(
            "Adding identity annotations..."
        )

        for i in range(n):

            for j in range(n):

                value = matrix[
                    i,
                    j
                ]

                if value == 100.00:

                    label = "100"

                else:

                    label = (
                        f"{value:.2f}"
                    )

                # -------------------------------------------------
                # Automatically choose text color for readability
                # -------------------------------------------------

                if value >= 99.90:

                    text_color = "white"

                else:

                    text_color = "black"

                ax_heatmap.text(

                    j,

                    i,

                    label,

                    ha="center",

                    va="center",

                    fontsize=6,

                    color=text_color

                )

    else:

        print(

            f"Identity annotations disabled "
            f"for {n} genomes."

        )


    # =========================================================
    # COLORBAR
    # =========================================================

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


    # ---------------------------------------------------------
    # Colorbar ticks
    # ---------------------------------------------------------

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


    if vmin < 97.00:

        colorbar_ticks.insert(

            0,

            round(
                vmin,
                2
            )

        )


    colorbar_ticks = sorted(
        set(

            [

                tick

                for tick in colorbar_ticks

                if vmin <= tick <= 100.00

            ]

        )
    )


    cbar.set_ticks(
        colorbar_ticks
    )


    cbar.set_ticklabels(

        [

            f"{tick:.2f}"

            for tick in colorbar_ticks

        ]

    )


    # =========================================================
    # TITLE
    # =========================================================

    fig.suptitle(

        "Pairwise Genome Sequence Identity and "
        "Phylogenetic Relationships of LSDV Complete Genomes",

        fontsize=font_size,

        fontweight="bold",

        y=0.985

    )


    # =========================================================
    # SAVE PNG
    # =========================================================

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

        "PNG saved:\n"

        f"{os.path.abspath(png_file)}"

    )


    # =========================================================
    # SAVE VECTOR PDF
    # =========================================================

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

        "PDF saved:\n"

        f"{os.path.abspath(pdf_file)}"

    )


    # =========================================================
    # CLEAN MEMORY
    # =========================================================

    plt.close(
        fig
    )

    gc.collect()


# =============================================================
# MAIN
# =============================================================

def main():

    print(
        "=" * 70
    )

    print(
        "LSDV PAIRWISE IDENTITY HEATMAP"
    )

    print(
        "WITH TOP PHYLOGENETIC TREE"
    )

    print(
        "=" * 70
    )


    # =========================================================
    # CHECK INPUT FILES
    # =========================================================

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

            "\nPhylogenetic tree file not found:\n"

            f"{tree_path}"

        )


    # =========================================================
    # STEP 1: READ ALIGNMENT
    # =========================================================

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


    # =========================================================
    # CHECK DUPLICATE IDS
    # =========================================================

    if len(sequence_ids) != len(
        set(sequence_ids)
    ):

        raise ValueError(

            "\nDuplicate sequence IDs detected."

            "\nEvery sequence must have a unique ID."

        )


    # =========================================================
    # STEP 2: CALCULATE PAIRWISE IDENTITY
    # =========================================================

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


    # Release alignment memory
    del alignment

    gc.collect()


    # =========================================================
    # STEP 3: READ TREE
    # =========================================================

    print(
        "\nSTEP 3: Reading phylogenetic tree..."
    )


    tree = read_and_prepare_tree(

        tree_path,

        sequence_ids

    )


    # =========================================================
    # STEP 4: GET PHYLOGENETIC ORDER
    # =========================================================

    print(
        "\nSTEP 4: Determining phylogenetic order..."
    )


    ordered_ids = get_tree_order(

        tree,

        sequence_ids

    )


    if len(ordered_ids) < 2:

        raise ValueError(

            "\nFewer than two sequences are shared "

            "between the alignment and tree."

        )


    print(

        "Number of sequences included in final matrix: "

        f"{len(ordered_ids)}"

    )


    # =========================================================
    # STEP 5: REORDER MATRIX
    # =========================================================

    print(
        "\nSTEP 5: Reordering identity matrix..."
    )


    ordered_matrix = reorder_matrix(

        identity_matrix,

        sequence_ids,

        ordered_ids

    )


    del identity_matrix

    gc.collect()


    # =========================================================
    # STEP 6: CREATE FIGURE
    # =========================================================

    print(

        "\nSTEP 6: Generating publication-ready heatmap..."

    )


    create_figure(

        matrix=ordered_matrix,

        ordered_ids=ordered_ids,

        tree=tree,

        output_prefix=output_path,

        dpi=dpi

    )


    # =========================================================
    # FINISHED
    # =========================================================

    print(
        "\n" + "=" * 70
    )

    print(
        "COMPLETED SUCCESSFULLY"
    )

    print(
        "=" * 70
    )

    print(
        "\nOutput files:"
    )

    print(
        f"  {output_path}.png"
    )

    print(
        f"  {output_path}.pdf"
    )


# =============================================================
# EXECUTE
# =============================================================

if __name__ == "__main__":

    main()
