import pandas as pd
from mhcflurry import Class1PresentationPredictor

# 1. CONFIGURATION
# Ensure this matches your patient's HLA (from OptiType)
# If unknown, use 'HLA-A*02:01' as a standard placeholder
hla_allele = "HLA-A*02:06" 
input_file = 'batch_peptides.csv'
output_file = 'high_affinity_candidates.csv'

print("Loading MHCflurry model (this takes a few seconds)...")
predictor = Class1PresentationPredictor.load()

# 2. LOAD PEPTIDES
try:
    df = pd.read_csv(input_file)
    # Deduplicate to speed up prediction
    unique_peptides = df['Peptide'].unique()
    print(f"Loaded {len(df)} rows. Found {len(unique_peptides)} unique peptides to test.")
except FileNotFoundError:
    print(f"Error: {input_file} not found. Run Step 1 first.")
    exit()

if len(unique_peptides) == 0:
    print("Error: No peptides found in file.")
    exit()

# 3. PREDICT BINDING
print(f"Predicting binding affinities against {hla_allele}...")
results = predictor.predict(
    peptides=unique_peptides,
    alleles=[hla_allele],
    verbose=1
)

# 4. MERGE & RANK
# Convert results to DataFrame
preds_df = pd.DataFrame(results)
# Merge back with original data to keep Gene IDs
final_df = df.merge(preds_df, left_on='Peptide', right_on='peptide')

# Filter for Strong Binders (<500 nM)
# <50 nM = High Affinity (Strongest candidates)
# <500 nM = Binder (Good candidates)
strong_binders = final_df[final_df['affinity'] < 500].copy()

# Sort: Lowest affinity value is best
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
    # Display key columns
    print(strong_binders[['Gene_ID', 'Mutation', 'peptide', 'affinity']].head(5))
    print(f"\nFull list saved to: {output_file}")
else:
    print("No strong binders found. Try a different HLA allele or check peptide generation.")
