import pandas as pd
import argparse

# -----------------------------
# Parse command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Filter and deduplicate non-coding RNAs")
parser.add_argument("--cpc2", required=True, help="Input CPC2 result file")
parser.add_argument("--plek", required=True, help="Input PLEK result file")
parser.add_argument("--out", required=True, help="Output file for final lncRNA list")

args = parser.parse_args()

# -----------------------------
# 1. Load CPC2 file
# -----------------------------
print("[INFO] Reading CPC2 file:", args.cpc2)

cpc2_df = pd.read_csv(args.cpc2, sep="\t")

# Detect ID column
id_col = "ID"
if "ID" not in cpc2_df.columns:
    if "#ID" in cpc2_df.columns:
        id_col = "#ID"
    else:
        raise ValueError("ERROR: No 'ID' or '#ID' column found in CPC2 file.")

# Normalize column name
cpc2_df[id_col] = cpc2_df[id_col].astype(str).str.strip()
cpc2_df.rename(columns={id_col: "ID"}, inplace=True)

# -----------------------------
# 2. Load PLEK file
# -----------------------------
print("[INFO] Reading PLEK file:", args.plek)

plek_data = []
with open(args.plek) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 3:
            label = parts[0].lower()
            score = float(parts[1])
            tid = parts[2].replace(">", "")
            plek_data.append([tid, label, score])

plek_df = pd.DataFrame(plek_data, columns=["ID", "PLEK_label", "PLEK_score"])

# -----------------------------
# 3. Merge CPC2 + PLEK
# -----------------------------
merged = cpc2_df.merge(plek_df, on="ID", how="outer")
print("[INFO] Total transcripts after merging:", merged.shape[0])

# -----------------------------
# 4. Identify ALL noncoding transcripts
#
# Keep transcripts where:
# - CPC2 says noncoding OR
# - PLEK says noncoding
# -----------------------------
noncoding = merged[
    (merged["label"].str.lower() == "noncoding") |
    (merged["PLEK_label"] == "noncoding")
]

print("[INFO] Number of noncoding transcripts before deduplication:", noncoding.shape[0])

# -----------------------------
# 5. Remove duplicate entries (keep ID once)
# -----------------------------
noncoding_unique = noncoding.drop_duplicates(subset=["ID"])

print("[INFO] Number of unique noncoding transcripts:", noncoding_unique.shape[0])

# -----------------------------
# 6. Save final lncRNA list
# -----------------------------
noncoding_unique.to_csv(args.out, index=False)
print("[INFO] Final noncoding (unique) transcripts saved to:", args.out)

