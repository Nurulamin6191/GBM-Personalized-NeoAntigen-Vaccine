#Script2 - For Peptide sequences
nano  2_peptide_seq.py
#####-------------------------------------------------------------------------------------------#####
import pandas as pd
import requests
import time
import sys

# ================= SETTINGS =================
input_file = 'candidates.csv'
output_file = 'batch_peptides.csv'
# Set to None to process ALL 2874 variants (Takes ~30-60 mins)
# Set to 50 to test first
LIMIT = None 
# ============================================

def get_sequence_by_id(gene_id):
    """Fetches sequence using Ensembl Stable ID (ENSG...)"""
    server = "https://rest.ensembl.org"
    # usage of the /sequence/id endpoint is more robust than symbol
    ext = f"/sequence/id/{gene_id}?type=protein;multiple_sequences=1;content-type=text/x-fasta"
    
    headers = { "Content-Type" : "application/json"}
    r = requests.get(server+ext, headers=headers)
    
    if r.status_code == 200:
        # The API might return multiple isoforms. We take the longest one.
        data = r.json()
        if len(data) > 0:
            # Sort by length, return the longest sequence
            data.sort(key=lambda x: len(x['seq']), reverse=True)
            return data[0]['seq']
    return None

def generate_mutant_peptides(seq, pos, new_aa):
    """Creates 9-mers around the mutation"""
    # Python index is pos-1
    idx = pos - 1
    
    if idx >= len(seq): return [] # Safety check
    
    # Create Mutated Full Sequence
    mutated_seq = seq[:idx] + new_aa + seq[idx+1:]
    
    # Extract 9-mers
    start = max(0, idx - 8)
    end = min(len(mutated_seq), idx + 9)
    region = mutated_seq[start:end]
    
    peptides = []
    for i in range(len(region) - 8):
        pep = region[i : i+9]
        if new_aa in pep: # Only keep peptides containing the mutation
            peptides.append(pep)
    return peptides

# 1. LOAD DATA
print(f"Reading {input_file}...")
df = pd.read_csv(input_file)

# Check for 'Gene' column (Ensembl ID)
if 'Gene' not in df.columns:
    print("CRITICAL ERROR: Column 'Gene' not found.")
    print("Columns available:", df.columns.tolist())
    sys.exit()

results = []
processed_count = 0

print(f"Starting batch processing (Limit: {LIMIT})...")

for index, row in df.iterrows():
    if LIMIT and processed_count >= LIMIT: break
    
    # A. Parse Mutation Info
    try:
        # Position
        if pd.isna(row['Protein_position']) or '-' in str(row['Protein_position']): continue
        pos = int(str(row['Protein_position']).split('/')[0])
        
        # Amino Acids (e.g. L/F)
        aa_change = row['Amino_acids']
        if '/' not in str(aa_change): continue
        old_aa, new_aa = aa_change.split('/')
        
        # Gene ID (ENSG...)
        gene_id = row['Gene']
        
    except (ValueError, IndexError):
        continue # Skip complex/malformed rows

    # B. Fetch & Generate
    print(f"[{processed_count+1}] Processing {gene_id} (Pos {pos})...")
    
    seq = get_sequence_by_id(gene_id)
    if seq:
        peps = generate_mutant_peptides(seq, pos, new_aa)
        for p in peps:
            results.append({
                'Gene_ID': gene_id,
                'Mutation': f"{old_aa}{pos}{new_aa}",
                'Peptide': p
            })
    else:
        print(f"   -> Failed to fetch sequence for {gene_id}")
        
    processed_count += 1
    # Be polite to the server
    time.sleep(0.1)

# Save
out_df = pd.DataFrame(results)
out_df.to_csv(output_file, index=False)
print(f"\nCOMPLETED. Generated {len(out_df)} peptides.")
print(f"Saved to {output_file}")
#####-----------------------------------------------------------------------------------------------------------#####


