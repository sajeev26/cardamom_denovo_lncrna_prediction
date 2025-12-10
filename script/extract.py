import argparse
from Bio import SeqIO
import pandas as pd

# ------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Extract lncRNA FASTA sequences using ID list and reference FASTA"
)

parser.add_argument("--lnc", required=True, help="High confidence lncRNA CSV file")
parser.add_argument("--fasta", required=True, help="Reference transcripts FASTA file")
parser.add_argument("--out", required=True, help="Output FASTA file for lncRNAs")

args = parser.parse_args()

# ------------------------------------------------------------
# Step 1: Load lncRNA ID list
# ------------------------------------------------------------
print("[INFO] Loading lncRNA CSV:", args.lnc)
df = pd.read_csv(args.lnc)

if "ID" not in df.columns:
    raise ValueError("ERROR: CSV must contain a column named 'ID'")

lnc_ids = set(df["ID"].astype(str).str.strip())

print(f"[INFO] Total lncRNA IDs loaded: {len(lnc_ids)}")

# ------------------------------------------------------------
# Step 2: Parse reference FASTA and extract matching sequences
# ------------------------------------------------------------
print("[INFO] Reading reference FASTA:", args.fasta)

records = SeqIO.parse(args.fasta, "fasta")
selected_records = []

for rec in records:
    if rec.id in lnc_ids:
        selected_records.append(rec)

print(f"[INFO] lncRNA sequences found: {len(selected_records)}")

# ------------------------------------------------------------
# Step 3: Write output FASTA
# ------------------------------------------------------------
SeqIO.write(selected_records, args.out, "fasta")
print("[INFO] lncRNA FASTA saved as:", args.out)

