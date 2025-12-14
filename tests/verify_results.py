import pandas as pd
import sys
import os

# Paths (passed as arguments)
if len(sys.argv) < 3:
    print("Usage: python verify_results.py <new_result> <expected_result>")
    sys.exit(1)

new_file = sys.argv[1]
expected_file = sys.argv[2]

if not os.path.exists(new_file):
    print(f"FAILED: Output file {new_file} was not created.")
    sys.exit(1)

# Load data
print(f"Loading files:\n New: {new_file}\n Expected: {expected_file}")
df_new = pd.read_csv(new_file)
df_expected = pd.read_csv(expected_file)

# --- PRE-PROCESSING FOR ROBUST COMPARISON ---

# 1. Sort Columns (ensure column order doesn't matter)
df_new = df_new.reindex(sorted(df_new.columns), axis=1)
df_expected = df_expected.reindex(sorted(df_expected.columns), axis=1)

# 2. Sort Rows (ensure row order doesn't matter)
# We assume the first column is the 'key' (e.g., Cell Barcode)
key_column = df_new.columns[0]
print(f"Sorting data by column: '{key_column}' to align rows...")

df_new = df_new.sort_values(by=key_column).reset_index(drop=True)
df_expected = df_expected.sort_values(by=key_column).reset_index(drop=True)

# --------------------------------------------

# Compare
try:
    # check_exact=False allows for tiny floating point differences
    # check_like=True is redundant if we sorted, but good safety
    pd.testing.assert_frame_equal(df_new, df_expected, check_exact=False, rtol=1e-5)
    print("SUCCESS: Results match expected data.")
except AssertionError as e:
    print("FAILED: Results differ from expected data!")
    print(e)
    # Optional: Print the first few differences to help debug
    print("\n--- Top of New Data ---")
    print(df_new.head())
    print("\n--- Top of Expected Data ---")
    print(df_expected.head())
    sys.exit(1)