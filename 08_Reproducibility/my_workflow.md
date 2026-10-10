conda create -n Gbm

##PHASE1-ACQUIRE DATA
#Install SRA toolkit(Downloads SRR files-under conda)
conda install -c bioconda sra-tools
#Setup SRA toolkit (Location of File Download)
vdb-config --interactive
#Download Raw Data
prefetch -p SRR30677899
#Uncompress (dump=unzip, gzip=compressed output, split=Forward&Reverse Reads)
fastq-dump --gzip --split-3 SRR30677899.sra

##PHASE2-QUALITY CHECK & TRIMMING
#Install fastqc & trim-galore
conda install -c bioconda fastqc trim-galore
#Quality check (Imp factor(Per Base Sequence Quality: Good = 28+, Bad <20.)
fastqc SRR30677899_1.fastq.gz SRR30677899_2.fastq.gz
#Trimming (trim=cut/remove, quality 20=quality less than 20)
trim_galore --quality 20 --fastqc SRR30677899_1.fastq.gz SRR30677899_2.fastq.gz

##PHASE3-ALIGNMENT WITH STAR ALIGNER & REFERENCE GENOME + GTF (Gene Map) DOWNLOAD 
#Install Star Aligner 
conda install -c bioconda star
#Download reference Genome from Ensembl (most Chromosome 19 used) As I know mutation is in Chromosome1
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
#unzip Chromosome & GTP (gunzip=only keep unzipped, *.gz=unzip all in .gz format)
gunzip *.gz
#Indexing(Creates star_index dir for further analysis)
STAR --runMode genomeGenerate \
--genomeDir ./star_index \
--genomeFastaFiles Homo_sapiens.GRCh38.dna.chromosome.1.fa \
--sjdbGTFfile Homo_sapiens.GRCh38.115.gtf \
--sjdbOverhang 100
#Alignment script(Creates .bam=gold, log.final.out=details, .out.tab=bases alignment, log.progress.out)
STAR --runMode alignReads \
--genomeDir ./star_index \
--readFilesIn SRR30677899_1_trimmed.fq.gz SRR30677899_2_trimmed.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix SRR_star \
--outSAMtype BAM SortedByCoordinate
#Output(log.final.out, uniquely mapped reads >80%=Excellent, 50%-80%=Okay, <50%=Poor
------------------------------------------------------------------------------------------------------------
#on error on Alignmnet (REPAIR USING BBMAP)
#Check error using 
rm SRR30677899_1_trimmed.fq.gz SRR30677899_2_trimmed.fq.gz
#Install Bioconda
conda install -c bioconda bbmap
#Repair command
repair.sh in1=SRR30677899_1_trimmed.fq.gz \
          in2=SRR30677899_2_trimmed.fq.gz \
          out1=SRR30677899_1_fixed.fq.gz \
          out2=SRR30677899_2_fixed.fq.gz \
          outs=singletons.fq.gz \
          repair=t
-------------------------------------------------------------------------------------------------------------

##PHASE4-VARIANT CALLING USING GATK (Gatk requires Picard(add reads group-Stamp)& Samtools for Indexing)
#Install Picard,Samtools and GATK
conda install -c bioconda picard gatk4 samtools
#Picard Add groups
picard AddOrReplaceReadGroups \
I=SRR_starAligned.sortedByCoord.out.bam \
O=SRR_picard_tagged.bam \
RGLB=illumina \
RGPU=unit1 \
RGSM=SRR30677899 \
RGPL=ILLUMINA
#Index using Samtools (Output=SRR_picard_tagged.bam.bai)####
samtools index SRR_picard_tagged.bam 
#Marks duplicate using Picard (OP=SRR_marked_dups.bam)
picard MarkDuplicates \
I=SRR_picard_tagged.bam \
O=SRR_marked_dups.bam \
M=marked_dup_metrics.txt
#Index using Samtools (Output=SRR_marked_dups.bam.bai)####
samtools index SRR_marked_dups.bam 
#Create reference Genome into Dictionary using Picard
picard CreateSequenceDictionary \
R=Homo_sapiens.GRCh38.dna.chromosome.1.fa \
O=Homo_sapiens.GRCh38.dna.chromosome.1.dict 
#Index reference Genome using samtools
samtools faidx Homo_sapiens.GRCh38.dna.chromosome.1.fa 
#Spltting Introns using gatk
gatk SplitNCigarReads \
-R Homo_sapiens.GRCh38.dna.chromosome.1.fa \
-I SRR_marked_dups.bam \
-O SRR_split.bam
#Variant calling (SNPs, Indels) using gatk
gatk HaplotypeCaller \
-R Homo_sapiens.GRCh38.dna.chromosome.1.fa \
-I SRR_split.bam \
-O SRR_variants.vcf
#Variant Filtering (FS=Fisher Strand Bias, Phred Score >30, QD (Quality by Depth)
gatk VariantFiltration \
-R Homo_sapiens.GRCh38.dna.chromosome.1.fa \
-V SRR_variants.vcf \
-O SRR_filtered.vcf \
--filter-name "FS_filter" --filter-expression "FS > 30.0" \
--filter-name "QD_filter" --filter-expression "QD < 2.0"
------------------------------------------------------------------------------------------------------------------------------------------------------------
##PHASE5-HLA Typing
#Install Optitype (HLA Typing)
conda install -c bioconda optitype
#HLA Typing (Folder hla_result)
OptiTypePipeline.py --input SRR30677899_1_fixed.fq.gz SRR30677899_2_fixed.fq.gz --rna --verbose --outdir hla_result

#Use the Web Interface Go to the Ensembl VEP Website. "https://www.ensembl.org/Multi/Tools/VEP" Upload your SRR_filtered.vcf file. Run it and download the result as a VCF/TXT file.
#Install Python dependencies
pip install pandas mhcflurry openpyxl requests
mhcflurry-downloads fetch

#Script1-To Identify the candidates
nano 1_Filter_candidates.py
#####---------------------------------------------------------------------------#####
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

#Script3- High Affinity Candidates 
nano 3_rank_candidates.py
#####-----------------------------------------------------------------------------------------------------------#####
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

#####-----------------------------------------------------------------------------------------------------------#####

#High Affinity (Findings in CMD)
head -n 11 high_affinity_candidates.csv | column -t -s,



---------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------
https://www.ensembl.org/Homo_sapiens/Info/Index  (Gene ID ===> Gene Name)
https://www.cbioportal.org/ (Gene contribution towards a cancer + Somatic Mutations)
https://portal.gdc.cancer.gov/ (Gene ===> Frequent Somatic Mutation list)
https://www.iedb.org/ (Immuno response prediction 
===> T-cell Epitope prediction ===> MHC I ====> Add Peptide sequence + Prediction Model) 
===> Immunogenic Score (>0 = good, 0.1 to 0.5 considered good)

----------------------------------------------------------------------------------------------------------------------------
##To check gene expression level
#Install subread
conda install -c bioconda subread

#Run the Counting Command (-T 4=4 threads for faster result, -p: Check for paired-end fragments)
featureCounts -T 4 -p -a Homo_sapiens.GRCh38.115.gtf -o gene_counts.txt SRR_starAligned.sortedByCoord.out.bam

#count the gene (End = Gene length, Read count) (Count < 5: Not Expressed, Count 5 - 50: Low Expression, Count > 50: High Expression. (Pass! This is a valid target.)
grep "ENSG00000XXXX" gene_counts.txt

