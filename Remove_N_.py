#!/usr/bin/env python3
"""
Write Code WSL
python3 Remove_N_.py

Remove N from the consensus sequences
Combine reference + consensus, align with MAFFT, replace N in consensus using reference.
Works in WSL and with Ubuntu apt MAFFT.
For mafft install
conda install -c bioconda mafft

"""

import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import sys

# ------------------------
# User inputs
# ------------------------
REFERENCE = "reference.fasta"
CONSENSUS = "consensus.fasta"
COMBINED = "combined.fasta"
ALIGNED = "aligned.fasta"
OUTPUT = "consensus_N_replaced.fasta"

# ------------------------
# Check files
# ------------------------
for f in (REFERENCE, CONSENSUS):
    if not Path(f).is_file():
        sys.exit(f"ERROR: Missing file → {f}")

# ------------------------
# Step 1: Combine FASTA
# ------------------------
print("Combining FASTA files...")
with open(COMBINED, "w") as out:
    for seq in SeqIO.parse(REFERENCE, "fasta"):
        SeqIO.write(seq, out, "fasta")
    for seq in SeqIO.parse(CONSENSUS, "fasta"):
        SeqIO.write(seq, out, "fasta")

# ------------------------
# Step 2: Run MAFFT
# ------------------------
print("Running MAFFT alignment...")

try:
    # Simple, WSL-safe MAFFT call
    subprocess.run(
        ["mafft", COMBINED],
        stdout=open(ALIGNED, "w"),
        stderr=sys.stderr,
        check=True
    )
except subprocess.CalledProcessError:
    sys.exit(
        "ERROR: MAFFT failed. Make sure MAFFT is installed and in PATH. "
        "Alternatively, run alignment manually:\n"
        "  mafft combined.fasta > aligned.fasta"
    )

# ------------------------
# Step 3: Replace N
# ------------------------
records = list(SeqIO.parse(ALIGNED, "fasta"))
if len(records) != 2:
    sys.exit("ERROR: Alignment must contain exactly 2 sequences (reference + consensus)")

ref, cons = records
ref_seq = str(ref.seq).upper()
cons_seq = str(cons.seq).upper()

corrected = []
replaced_count = 0
for r, c in zip(ref_seq, cons_seq):
    if c == "N" and r in "ATGC":
        corrected.append(r)
        replaced_count += 1
    else:
        corrected.append(c)

# Remove alignment gaps
final_seq = "".join(corrected).replace("-", "")
cons.seq = Seq(final_seq)

# ------------------------
# Step 4: Write final corrected consensus
# ------------------------
SeqIO.write(cons, OUTPUT, "fasta")

# ------------------------
# Step 5: Cleanup temp file
# ------------------------
Path(COMBINED).unlink(missing_ok=True)

print("Done ✔")
print(f"N replaced: {replaced_count}")
print(f"Output file: {OUTPUT}")
