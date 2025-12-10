import pandas as pd
import argparse

# -----------------------------
# Parse arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Filter high-confidence lncRNAs")
parser.add_argument("--in", required=True, help="Input merged lncRNA CSV file")
parser.add_argument("--out", required=True, help="Output filtered lncRNA CSV")
args = parser.parse_args()

input_file = getattr(args, "in")
output_file = args.out

# -----------------------------
# Load data
# -----------------------------
print("[INFO] Reading input file:", input_file)
df = pd.read_csv(input_file)

print("[INFO] Initial transcripts:", df.shape[0])

# -----------------------------
# Apply Filtering Rules
# -----------------------------
filtered = df[
    (df['label'].str.lower() == 'noncoding') &
    (df['PLEK_label'].str.lower().str.replace("-", "") == 'noncoding') &
    (df['coding_probability'] <= 0.1) &
    (df['PLEK_score'] <= -1) &
    (df['peptide_length'] < 50) &
    (df['Fickett_score'] < 0.5)
]

# Optional strict rule
filtered = filtered[filtered['ORF_integrity'] == -1]

# Remove duplicate IDs
filtered = filtered.drop_duplicates(subset=['ID'])

print("[INFO] After filtering:", filtered.shape[0])

# -----------------------------
# Save Output
# -----------------------------
filtered.to_csv(output_file, index=False)
print("[INFO] Saved:", output_file)

