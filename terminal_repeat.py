#!/usr/bin/env python3
"""
Code below 
python3 terminal_repeat.py sequence.fasta

auto_itr.py
Automatically determine the longest inverted terminal repeat (ITR)
from FASTA sequences by scanning the entire sequence.
"""

import argparse


# ----------------------- Utility functions -----------------------

def read_fasta(path: str):
    """Read FASTA file and return dict {header: sequence}."""
    seqs = {}
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq_lines).upper()
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            seqs[header] = "".join(seq_lines).upper()
    return seqs


def revcomp(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def hamming(a: str, b: str) -> int:
    """Count mismatches between two same-length strings."""
    return sum(x != y for x, y in zip(a, b))


# ----------------------- ITR detection -----------------------

def find_longest_itr_auto(seq: str, mismatches: int = 0):
    """
    Automatically find the longest inverted terminal repeat (ITR)
    by scanning all possible repeat lengths up to half the sequence.
    Returns (length, repeat_seq, mismatches) or None.
    """
    n = len(seq)
    best = None
    # Scan from largest possible repeat down to 1
    for L in range(n // 2, 0, -1):
        left = seq[:L]
        right = revcomp(seq[-L:])
        mis = hamming(left, right)
        if mis <= mismatches:
            best = (L, left, mis)
            break  # stop at first (longest) valid repeat
    return best


# ----------------------- Main -----------------------

def main():
    parser = argparse.ArgumentParser(description="Automatically find the longest inverted terminal repeat (ITR).")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("--mismatches", type=int, default=0, help="Allowed mismatches (default: 0)")
    parser.add_argument("--out", default="auto_itr_results.txt", help="Output results file")
    args = parser.parse_args()

    seqs = read_fasta(args.fasta)
    results = []
    header_line = "Sequence_ID\tLength\tMismatches\tRepeat_Sequence\n"
    results.append(header_line)

    for name, seq in seqs.items():
        print("=" * 80)
        print(f"> {name}")
        print(f"Sequence length: {len(seq)} bp")

        best = find_longest_itr_auto(seq, mismatches=args.mismatches)
        if best:
            L, rep, mis = best
            print(f"\nLongest ITR found:")
            print(f"  Length: {L}")
            print(f"  Mismatches: {mis}")
            print(f"  Repeat: {rep}")
            results.append(f"{name}\t{L}\t{mis}\t{rep}\n")
        else:
            print("\nNo inverted terminal repeat found.")
            results.append(f"{name}\t0\tNA\tNo_repeat_found\n")

    with open(args.out, "w") as f:
        f.writelines(results)

    print("\nResults saved to:", args.out)


if __name__ == "__main__":
    main()
