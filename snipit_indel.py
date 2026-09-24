#!/usr/bin/env python3

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Input/output files
input_file = "alignment_326_152358.fasta"
output_file = "alignment_SNPmasked.fasta"

# Read alignment
alignment = AlignIO.read(input_file, "fasta")
alignment_length = alignment.get_alignment_length()

# Find SNP positions
snp_positions = []

for i in range(alignment_length):
    column = alignment[:, i]

    # Ignore gaps and Ns when determining SNPs
    bases = {b.upper() for b in column if b.upper() not in {"N", "-"}}

    if len(bases) > 1:
        snp_positions.append(i)

print(f"Total SNP positions found: {len(snp_positions)}")

# Mask SNPs without changing gaps or existing Ns
masked_records = []

for record in alignment:
    seq = list(str(record.seq))

    for pos in snp_positions:
        # Replace only A/C/G/T with N
        if seq[pos].upper() in {"A", "C", "G", "T"}:
            seq[pos] = "N"
        # Leave '-' and existing 'N' unchanged

    masked_records.append(
        SeqRecord(
            Seq("".join(seq)),
            id=record.id,
            description=record.description
        )
    )

# Save alignment
masked_alignment = MultipleSeqAlignment(masked_records)
AlignIO.write(masked_alignment, output_file, "fasta")

print(f"Masked alignment written to: {output_file}")
