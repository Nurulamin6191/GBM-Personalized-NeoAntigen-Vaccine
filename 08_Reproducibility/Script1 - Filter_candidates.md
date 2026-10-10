**Script1-To Identify the candidates**
        nano 1_Filter_candidates.py
______________________________________________________________________________________________________________
    import pandas as pd
    import sys

# 1. SETUP
# Ensure this matches your filename exactly (Case Sensitive!)
input_file = 'SRR_annotated.txt'

print(f"Scanning {input_file}...")

# 2. COUNT METADATA LINES
# We need to calculate how many lines start with "##" so we can skip them.
skip_count = 0
try:
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                skip_count += 1
            else:
                break # We found the header line (starts with single #)
except FileNotFoundError:
    print(f"ERROR: Could not find file '{input_file}'.")
    print("Please check if the file is in this folder and the name matches exactly.")
    sys.exit()

print(f"Detected {skip_count} metadata lines. Reading data...")

# 3. LOAD DATA
# We skip the "##" lines. We do NOT use the comment parameter.
df = pd.read_csv(input_file, sep='\t', skiprows=skip_count)

# 4. CLEANUP HEADERS
# The first column usually has a '#' (e.g., #Uploaded_variation). We remove it.
df.rename(columns={col: col.replace('#', '') for col in df.columns}, inplace=True)

# 5. FILTER FOR MISSENSE (Neoantigen Candidates)
if 'Consequence' not in df.columns:
    print("ERROR: Column 'Consequence' not found. Check your input file format.")
    sys.exit()

# Filter for rows containing 'missense_variant'
missense_df = df[df['Consequence'].str.contains('missense_variant', na=False)]

# Filter for rows that have Amino Acid changes (remove empty ones)
if 'Amino_acids' in df.columns:
    # Use copy() to avoid SettingWithCopyWarning
    filtered_df = missense_df.dropna(subset=['Amino_acids']).copy()
else:
    filtered_df = missense_df.copy()

# 6. SAVE & REPORT
output_file = 'candidates.csv'
filtered_df.to_csv(output_file, index=False)

print(f"------------------------------------------------")
print(f"Total Variants:    {len(df)}")
print(f"Candidates Found:  {len(filtered_df)}")
print(f"Saved to:          {output_file}")
print(f"------------------------------------------------")

if not filtered_df.empty:
    print("TOP CANDIDATE INFO (Use for Step 2):")
    # Display the first row
    first_row = filtered_df.iloc[0]
    print(f"Gene: {first_row.get('Extra', 'Unknown')}") 
    print(f"Mutation: {first_row.get('Amino_acids', 'Unknown')}")
    print(f"Location: {first_row.get('Location', 'Unknown')}")
else:
    print("No missense variants found.")

#####------------------------------------------------------------------------------------------#####


