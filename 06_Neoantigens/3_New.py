import pandas as pd
from mhcflurry import Class1PresentationPredictor

# 1. CONFIGURATION
# Ensure this matches your patient's HLA (from OptiType)
hla_allele = "HLA-A*02:06" 
input_file = 'batch_peptides.csv'
output_file = 'high_affinity_candidates.csv'

print("Loading MHCflurry model (this takes a few seconds)...")
predictor = Class1PresentationPredictor.load()

# 2. LOAD PEPTIDES
try:
    df = pd.read_csv(input_file)
    print(f"Loaded {len(df)} rows.")
except FileNotFoundError:
    print(f"Error: {input_file} not found. Run Step 1 first.")
    exit()

# === [NEW] 2.5 SAFETY FILTER ===
# Remove peptides with non-standard amino acids (like 'U', 'X', '*', etc.)
valid_aas = set("ACDEFGHIKLMNPQRSTVWY")

def is_valid(seq):
    # Returns True only if ALL letters in seq are standard
    return set(seq).issubset(valid_aas)

# Apply filter
initial_count = len(df)
df = df[df['Peptide'].apply(is_valid)].copy()
dropped_count = initial_count - len(df)

if dropped_count > 0:
    print(f"WARNING: Removed {dropped_count} peptides containing non-standard amino acids (e.g., U, X, *).")

unique_peptides = df['Peptide'].unique()
print(f"Testing {len(unique_peptides)} unique valid peptides.")

if len(unique_peptides) == 0:
    print("Error: No valid peptides left to test.")
    exit()
# ===============================

# 3. PREDICT BINDING
print(f"Predicting binding affinities against {hla_allele}...")
results = predictor.predict(
    peptides=unique_peptides,
    alleles=[hla_allele],
    verbose=1
)

# 4. MERGE & RANK
preds_df = pd.DataFrame(results)
final_df = df.merge(preds_df, left_on='Peptide', right_on='peptide')

# Filter for Strong Binders (<500 nM)
strong_binders = final_df[final_df['affinity'] < 500].copy()
strong_binders.sort_values(by='affinity', inplace=True)

# 5. SAVE
strong_binders.to_csv(output_file, index=False)

print("\n" + "="*40)
print(f"ANALYSIS COMPLETE")
print(f"Total Peptides Tested: {len(df)}")
print(f"High Binding Candidates Found (<500nM): {len(strong_binders)}")
print("="*40)

if not strong_binders.empty:
    print("\nTOP 5 NEOANTIGENS:")
    print(strong_binders[['Gene_ID', 'Mutation', 'peptide', 'affinity']].head(5))
    print(f"\nFull list saved to: {output_file}")
else:
    print("No strong binders found. Try a different HLA allele.")
