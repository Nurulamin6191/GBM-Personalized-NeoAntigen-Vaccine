import pandas as pd
import sys
import re

# ================= CONFIGURATION =================
input_file = 'SRR_annotated.vcf'   # Your VCF file
output_file = 'candidates.csv'     # Output for Step 2
# =================================================

print(f"Scanning {input_file}...")

# 1. READ VCF AND SKIP METADATA
# We need to find where the header line (#CHROM) starts
skip_count = 0
try:
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                skip_count += 1
            else:
                break
except FileNotFoundError:
    print(f"ERROR: Could not find '{input_file}'.")
    sys.exit()

print(f"Reading data (Skipping {skip_count} metadata lines)...")

# Load VCF into DataFrame
df = pd.read_csv(input_file, sep='\t', skiprows=skip_count)
# Clean up header (#CHROM -> CHROM)
df.rename(columns={col: col.replace('#', '') for col in df.columns}, inplace=True)

# 2. DEFINE PARSING LOGIC
# We need to extract: Gene ID (ENSG), Protein Position (123), and Amino Acid Change (L/P)
# This function handles both VEP (CSQ) and SnpEff (ANN) formats.

# 3-letter to 1-letter amino acid map (for SnpEff conversion)
aa_map = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E',
    'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
    'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', '*': '*'
}

def parse_info(info_str):
    gene_id = None
    pos = None
    aa_change = None
    
    # --- STRATEGY A: VEP (CSQ) ---
    if 'CSQ=' in info_str:
        try:
            # Extract the CSQ block
            csq = info_str.split('CSQ=')[1].split(';')[0]
            # Take first transcript
            first_trans = csq.split(',')[0]
            fields = first_trans.split('|')
            
            # VEP Standard Indices: 
            # 4=Gene(ENSG), 14=Protein_pos, 15=Amino_acids(L/P)
            if len(fields) >= 16:
                if 'missense_variant' in fields[1]: # Index 1 is Consequence
                    gene_id = fields[4]
                    pos = fields[14]
                    aa_change = fields[15] # Already in L/P format
        except Exception:
            pass

    # --- STRATEGY B: SnpEff (ANN) ---
    elif 'ANN=' in info_str:
        try:
            # Extract ANN block
            ann = info_str.split('ANN=')[1].split(';')[0]
            first_trans = ann.split(',')[0]
            fields = first_trans.split('|')
            
            # SnpEff Indices: 
            # 1=Annotation, 4=Gene_ID(ENSG), 10=HGVS.p(p.Ala12Val), 13=AA.pos
            if len(fields) >= 14:
                if 'missense_variant' in fields[1]:
                    gene_id = fields[4]
                    
                    # Clean Position (120/400 -> 120)
                    raw_pos = fields[13]
                    pos = raw_pos.split('/')[0] if raw_pos else None
                    
                    # Convert p.Ala12Val -> A/V
                    hgvsp = fields[10]
                    if hgvsp.startswith('p.'):
                        # Regex to capture p.Ala12Val -> Ala, Val
                        match = re.match(r'p\.([A-Z][a-z]{2})\d+([A-Z][a-z]{2})', hgvsp)
                        if match:
                            ref = aa_map.get(match.group(1), '?')
                            alt = aa_map.get(match.group(2), '?')
                            aa_change = f"{ref}/{alt}"
        except Exception:
            pass
            
    return pd.Series([gene_id, pos, aa_change])

# 3. APPLY PARSING
if 'INFO' not in df.columns:
    print("CRITICAL ERROR: No 'INFO' column found in VCF.")
    sys.exit()

print("Parsing VCF annotations (Gene, Position, Amino Acids)...")
df[['Gene', 'Protein_position', 'Amino_acids']] = df['INFO'].apply(parse_info)

# 4. FILTERING
# Remove rows where we couldn't find a Gene or Amino Acid change (i.e., not missense)
filtered_df = df.dropna(subset=['Gene', 'Amino_acids']).copy()

# Remove rows where Amino_acids format is wrong (must contain '/')
filtered_df = filtered_df[filtered_df['Amino_acids'].str.contains('/')]

# 5. SAVE
# Keep only necessary columns for the next script
final_cols = ['Gene', 'Protein_position', 'Amino_acids']
filtered_df[final_cols].to_csv(output_file, index=False)

print("-" * 40)
print(f"Total Variants in VCF: {len(df)}")
print(f"Missense Candidates:   {len(filtered_df)}")
print(f"Saved to:              {output_file}")
print("-" * 40)

if not filtered_df.empty:
    print("Top Candidate Preview:")
    print(filtered_df.iloc[0][final_cols])
else:
    print("No missense candidates found! Check your VCF annotation.")
