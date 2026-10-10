#Create Conda Environment
        conda create -n Gbm
________________________________________________________________________________________________________________________________
**PHASE1-ACQUIRE DATA**
#Install SRA toolkit(Downloads SRR files-under conda)
        conda install -c bioconda sra-tools
#Setup SRA toolkit (Location of File Download)
        vdb-config --interactive
#Download Raw Data
        prefetch -p SRR30677899
#Uncompress (dump=unzip, gzip=compressed output, split=Forward&Reverse Reads)
        fastq-dump --gzip --split-3 SRR30677899.sra
_________________________________________________________________________________________________________________________________
**PHASE2-QUALITY CHECK & TRIMMING**
#Install fastqc & trim-galore
        conda install -c bioconda fastqc trim-galore
#Quality check (Imp factor(Per Base Sequence Quality: Good = 28+, Bad <20.)
        fastqc SRR30677899_1.fastq.gz SRR30677899_2.fastq.gz
#Trimming (trim=cut/remove, quality 20=quality less than 20)
        trim_galore --quality 20 --fastqc SRR30677899_1.fastq.gz SRR30677899_2.fastq.gz
_________________________________________________________________________________________________________________________________
**PHASE3-ALIGNMENT WITH STAR ALIGNER & REFERENCE GENOME + GTF (Gene Map) DOWNLOAD** 
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
- Output(log.final.out, uniquely mapped reads >80%=Excellent, 50%-80%=Okay, <50%=Poor
__________________________________________________________________________________________________________________________________
**on error on Alignmnet (REPAIR USING BBMAP**)
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
___________________________________________________________________________________________________________________________________
**PHASE4-VARIANT CALLING USING GATK (Gatk requires Picard(add reads group-Stamp)& Samtools for Indexing)**
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
____________________________________________________________________________________________________________________________________
**PHASE5-HLA Typing**
#Install Optitype (HLA Typing)
        conda install -c bioconda optitype
#HLA Typing (Folder hla_result)
        OptiTypePipeline.py --input SRR30677899_1_fixed.fq.gz SRR30677899_2_fixed.fq.gz --rna --verbose --outdir hla_result

#Use the Web Interface Go to the Ensembl VEP Website.
- "https://www.ensembl.org/Multi/Tools/VEP" Upload your SRR_filtered.vcf file. Run it and download the result as a VCF/TXT file.
#Install Python dependencies
        pip install pandas mhcflurry openpyxl requests
        mhcflurry-downloads fetch
        
---> **Customized Python Script 1, 2, 3**
_____________________________________________________________________________________________________________________________________

#High Affinity (Findings in CMD)
        head -n 11 high_affinity_candidates.csv | column -t -s,
_____________________________________________________________________________________________________________________________________
- https://www.ensembl.org/Homo_sapiens/Info/Index  (Gene ID ===> Gene Name)
- https://www.cbioportal.org/ (Gene contribution towards a cancer + Somatic Mutations)
- https://portal.gdc.cancer.gov/ (Gene ===> Frequent Somatic Mutation list)
- https://www.iedb.org/ (Immuno response prediction 
        T-cell Epitope prediction ===> MHC I ====> Add Peptide sequence + Prediction Model) 
        Immunogenic Score (>0 = good, 0.1 to 0.5 considered good)

