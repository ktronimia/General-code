from Bio import SeqIO

# Input and output files
input_fasta = "input_file.fasta"
output_fasta = "hbv_C_genes.fasta"

# Coordinates of S gene (1-based positions as in GenBank)
start1 = 2845
end1   = None  # means until end of genome
start2 = 1
end2   = 835

records_out = []

for record in SeqIO.parse(input_fasta, "fasta"):
    seq = record.seq
    length = len(seq)

    # First fragment: from 2848 to end of genome
    frag1 = seq[start1-1:length]  # python is 0-based

    # Second fragment: from 1 to 835
    frag2 = seq[start2-1:end2]

    # Concatenate fragments
    sgene_seq = frag1 + frag2

    # Create new record
    new_record = record[:0]  # copy metadata
    new_record.id = record.id
    new_record.seq = sgene_seq

    records_out.append(new_record)

# Write output fasta
SeqIO.write(records_out, output_fasta, "fasta")

print(f"Extracted S genes written to {output_fasta}")

