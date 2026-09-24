from Bio import SeqIO
import pandas as pd

# Load aligned FASTA file
fasta_file = "aligned.fasta"
records = list(SeqIO.parse(fasta_file, "fasta"))

# Reference sequence (first sequence)
reference = str(records[0].seq)
reference_id = records[0].id

mutation_results = []

# Compare each sample against reference
for record in records[1:]:
    
    # Reference-based position counter
    ref_position = 0
    sample_id = record.id
    sample_seq = str(record.seq)

    mutations = []

    for ref_nt, sample_nt in zip(reference, sample_seq):

        # Skip gaps in reference sequence
        if ref_nt == "-":
            continue

        # Count only reference positions
        ref_position += 1

        # Ignore gaps in query sequence if needed
        if sample_nt == "-":
            continue

        # Identify mutation using reference-based numbering
        if ref_nt != sample_nt:
            mutation = f"{ref_nt}{ref_position}{sample_nt}"
            mutations.append(mutation)

    mutation_results.append({
        "Sample": sample_id,
        "Total_Mutations": len(mutations),
        "Mutations": ", ".join(mutations)
    })

# Convert to DataFrame
mutation_df = pd.DataFrame(mutation_results)

# Save results
mutation_df.to_csv("mutation_results.tsv", sep="	", index=False)

# Print results
print(mutation_df)
