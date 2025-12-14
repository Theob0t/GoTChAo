# tests/verify_results.py
import pandas as pd
import sys
import os

# Paths (passed as arguments)
new_file = sys.argv[1]
expected_file = sys.argv[2]

if not os.path.exists(new_file):
    print(f"FAILED: Output file {new_file} was not created.")
    sys.exit(1)

# Load data
df_new = pd.read_csv(new_file)
df_expected = pd.read_csv(expected_file)

# Compare
try:
    # check_exact=False allows for tiny floating point differences
    pd.testing.assert_frame_equal(df_new, df_expected, check_exact=False, rtol=1e-5)
    print("SUCCESS: Results match expected data.")
except AssertionError as e:
    print("FAILED: Results differ from expected data!")
    print(e)
    sys.exit(1)