from Bio import SeqIO

input_fasta = "input.fasta"
output_fastq = "output.fastq"

with open(output_fastq, "w") as out_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Create fake quality scores (Phred score = 40, represented by 'I')
        record.letter_annotations["phred_quality"] = [40] * len(record.seq)
        SeqIO.write(record, out_handle, "fastq")
