from Bio import AlignIO
import matplotlib.pyplot as plt

# Load the alignment file (FASTA format)
input_file = "aligned.fasta"
alignment = AlignIO.read(input_file, "fasta")

query = alignment[0]  # Use the first sequence as reference/query
seq_count = len(alignment)
seq_length = alignment.get_alignment_length()

# Store alignment ranges for plotting
plot_data = []

for record in alignment[1:]:  # skip the query itself
    start = None
    for i, (q_base, s_base) in enumerate(zip(query.seq, record.seq)):
        if s_base != "-" and start is None:
            start = i
        elif s_base == "-" and start is not None:
            plot_data.append((record.id, start, i))
            start = None
    if start is not None:
        plot_data.append((record.id, start, seq_length))

# Plot
fig, ax = plt.subplots(figsize=(12, (seq_count+1)*0.5))  # +1 for query
y_positions = {record.id: i+1 for i, record in enumerate(alignment[1:], 1)}

# Add the query sequence as a full-length bar
ax.broken_barh([(0, seq_length)], (0.6, 0.8), facecolors='lightgreen')
ax.text(-5, 1, query.id, ha='right', va='center', fontsize=8)  # use real ID

# Add alignments
for seq_id, start, end in plot_data:
    ax.broken_barh([(start, end-start)], (y_positions[seq_id]-0.4, 0.8),
                   facecolors='skyblue')

# Add labels for other sequences
for seq_id, y in y_positions.items():
    ax.text(-5, y, seq_id, ha='right', va='center', fontsize=8)

ax.set_xlim(0, seq_length)
ax.set_xlabel("Query Sequence Position")
ax.set_yticks([])
ax.set_title("Alignment Summary (BLAST-style)")
plt.tight_layout()

# ---------- Auto-save with high resolution ----------
output_file = input_file.replace(".fasta", "_summary.png")
plt.savefig(output_file, dpi=1000, bbox_inches="tight")  # high-res save
print(f"Figure saved as {output_file}")
# ----------------------------------------------------

plt.show()
