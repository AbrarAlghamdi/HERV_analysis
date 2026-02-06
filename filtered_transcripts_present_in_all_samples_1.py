#!/usr/bin/env python3
import os
import glob
import pandas as pd

# ------------ SETTINGS ------------
INPUT_DIR = "."           # directory with Telescope TSV files
OUTPUT_DIR = "./finalcount_common"
MIN_COUNT = 10
# ----------------------------------

os.makedirs(OUTPUT_DIR, exist_ok=True)

def is_header_line(line: str) -> bool:
    """Detect Telescope table header (quoted or unquoted)."""
    s = line.strip().replace('"', '').replace("'", "")
    return s.lower().startswith("transcript\t") and "final_count" in s.lower()

def read_telescope_tsv(path: str) -> pd.DataFrame:
    """Read Telescope TSV with RunInfo lines + quoted headers."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    header_idx = None
    for i, line in enumerate(lines):
        if is_header_line(line):
            header_idx = i
            break

    if header_idx is None:
        raise RuntimeError(f"No Telescope header found in {path}")

    df = pd.read_csv(
        path,
        sep="\t",
        header=header_idx,
        engine="python",
        quotechar='"'
    )

    df.columns = df.columns.astype(str).str.strip().str.replace('"', '')
    df["transcript"] = df["transcript"].astype(str).str.strip().str.replace('"', '')
    df["final_count"] = pd.to_numeric(df["final_count"], errors="coerce").fillna(0)

    return df[["transcript", "final_count"]]

# ------------ MAIN ------------
tsv_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.tsv")))
if len(tsv_files) < 2:
    raise RuntimeError("Need at least two TSV files.")

print(f"Found {len(tsv_files)} TSV files\n")

dfs = {}

for f in tsv_files:
    print(f"Reading: {f}")
    df = read_telescope_tsv(f)

    # Filter final_count ≥ MIN_COUNT
    df = df[df["final_count"] >= MIN_COUNT].copy()
    print(f"  → {df.shape[0]} transcripts after final_count ≥ {MIN_COUNT}")

    dfs[f] = df.set_index("transcript")

# ---- Transcripts present in ALL files ----
common_transcripts = set.intersection(*[set(d.index) for d in dfs.values()])
common_transcripts = sorted(common_transcripts)

print(f"\nTranscripts present in ALL files: {len(common_transcripts)}\n")

# ---- Build merged final_count matrix ----
merged = pd.DataFrame(index=common_transcripts)

for f, df in dfs.items():
    sample = os.path.basename(f).replace(".tsv", "")
    merged[sample] = df.loc[common_transcripts, "final_count"].astype(int)

out_matrix = os.path.join(
    OUTPUT_DIR,
    f"merged_final_count_matrix_common_finalcount{MIN_COUNT}.tsv"
)
merged.to_csv(out_matrix, sep="\t")

print(f"Saved merged final_count matrix:\n{out_matrix}")
print("\nDone ✅")
