# BMI 503 - Genomic Annotations and Motif Finding
## From GENCODE to p53 Binding Sites

**Date:** December 2, 2024  
**Instructor:** Pratik Dutta  
**Duration:** 90 minutes

---

## 📚 Table of Contents

1. [Introduction to GENCODE](#section1)
2. [Understanding GTF and GFF Formats](#section2)
3. [Downloading and Parsing GENCODE Data](#section3)
4. [Creating a DataFrame from GTF](#section4)
5. [Extracting TSS (Transcription Start Sites)](#section5)
6. [Finding Splice Sites](#section6)
7. [Extracting Promoter Regions (-1kb to TSS)](#section7)
8. [Getting FASTA Sequences](#section8)
9. [Finding p53 Binding Sites with Regex](#section9)
10. [Real-World Example: Complete Pipeline](#section10)

---

<a id="section1"></a>
## 1️⃣ Introduction to GENCODE

### What is GENCODE?

**GENCODE** = **GEN**ome en**CODE** - A comprehensive annotation of the human genome

```
🧬 GENCODE provides:
├── Gene locations (where genes are on chromosomes)
├── Transcript structures (exons, introns)
├── Protein-coding genes
├── Non-coding RNAs (lncRNA, miRNA)
└── Regulatory elements
```

### Why is GENCODE important?

- **Gold standard** annotation for human and mouse genomes
- Used by major consortia (ENCODE, GTEx, TCGA)
- Updated regularly with new discoveries
- Essential for RNA-seq, variant interpretation, regulatory analysis

### Key Concepts

```python
# Let's understand the hierarchy

Genome
  └── Chromosome (e.g., chr17)
        └── Gene (e.g., TP53)
              └── Transcript (e.g., TP53-201)
                    └── Exon 1, Exon 2, ... Exon 11
                          └── CDS (Coding Sequence)
```

**Example:**
```
Gene: TP53 (tumor protein p53)
Location: Chromosome 17
Transcripts: 12 different isoforms
Function: Tumor suppressor - "Guardian of the Genome"
```

---

<a id="section2"></a>
## 2️⃣ Understanding GTF and GFF Formats

### Setup

```python
# Install required packages
!pip install pandas biopython requests -q

import pandas as pd
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

print("✅ Setup complete!")
```

### What are GTF and GFF?

Both are **tab-delimited text files** that describe genomic features.

#### GTF (Gene Transfer Format)
- Standard format for GENCODE
- 9 columns
- Human/mouse focused

#### GFF (General Feature Format)
- More general, used for all organisms
- Version 3 (GFF3) is most common
- More flexible attribute syntax

### GTF File Structure

```
Format: 9 tab-separated columns

1. seqname     - Chromosome name
2. source      - Who created this annotation
3. feature     - Type (gene, transcript, exon, CDS)
4. start       - Start position (1-based, inclusive)
5. end         - End position (1-based, inclusive)
6. score       - Quality score (. if not applicable)
7. strand      - + (forward) or - (reverse)
8. frame       - Reading frame (0, 1, 2, or .)
9. attribute   - Semicolon-separated key-value pairs
```

### Let's Look at Real GTF Data

```python
# Create a sample GTF to understand the format
sample_gtf = """chr17	HAVANA	gene	7661779	7687550	.	-	.	gene_id "ENSG00000141510"; gene_type "protein_coding"; gene_name "TP53";
chr17	HAVANA	transcript	7661779	7687550	.	-	.	gene_id "ENSG00000141510"; transcript_id "ENST00000269305"; gene_name "TP53"; transcript_name "TP53-201";
chr17	HAVANA	exon	7687490	7687550	.	-	.	gene_id "ENSG00000141510"; transcript_id "ENST00000269305"; exon_number "1";
chr17	HAVANA	CDS	7687490	7687538	.	-	0	gene_id "ENSG00000141510"; transcript_id "ENST00000269305";
chr17	HAVANA	exon	7676592	7676736	.	-	.	gene_id "ENSG00000141510"; transcript_id "ENST00000269305"; exon_number "2";
chr17	HAVANA	CDS	7676594	7676736	.	-	1	gene_id "ENSG00000141510"; transcript_id "ENST00000269305";
"""

# Save it
with open('sample_tp53.gtf', 'w') as f:
    f.write(sample_gtf)

print("✅ Sample GTF created: sample_tp53.gtf")
print("\nLet's view it:")
print(sample_gtf)
```

### Understanding Each Line

```python
# Let's parse and explain one line
example_line = "chr17	HAVANA	gene	7661779	7687550	.	-	.	gene_id \"ENSG00000141510\"; gene_name \"TP53\";"

fields = example_line.split('\t')

print("📖 Breaking down the GTF line:\n")
print(f"1. Chromosome:  {fields[0]}")
print(f"2. Source:      {fields[1]}")
print(f"3. Feature:     {fields[2]}")
print(f"4. Start:       {fields[3]}")
print(f"5. End:         {fields[4]}")
print(f"6. Score:       {fields[5]}")
print(f"7. Strand:      {fields[6]}")
print(f"8. Frame:       {fields[7]}")
print(f"9. Attributes:  {fields[8]}")

print("\n🔍 This means:")
print("   TP53 gene is on chromosome 17")
print("   Located from position 7,661,779 to 7,687,550")
print("   On the MINUS strand (reads backward)")
print("   Total length: ~26 kb")
```

### Key Differences: GTF vs GFF3

```python
comparison = pd.DataFrame({
    'Feature': ['File Extension', 'Attribute Format', 'Hierarchy', 'Use Case'],
    'GTF': ['.gtf', 'key "value";', 'Implicit', 'Human/Mouse (GENCODE/ENSEMBL)'],
    'GFF3': ['.gff or .gff3', 'key=value;', 'Explicit with ID/Parent', 'All organisms']
})

print(comparison.to_string(index=False))
```

---

<a id="section3"></a>
## 3️⃣ Downloading and Parsing GENCODE Data

### Where to Get GENCODE Data?

**GENCODE Website:** https://www.gencodegenes.org/

```
Available files:
├── Comprehensive gene annotation (GTF) - ALL genes
├── Basic gene annotation (GTF) - Filtered subset
├── Transcript sequences (FASTA)
├── Protein sequences (FASTA)
└── Genome sequence (FASTA)
```

### Downloading GTF File

```python
import urllib.request
import gzip

# For this class, we'll download chromosome 22 (smaller, faster)
# Full genome GTF is ~1.5 GB!

print("📥 Downloading GENCODE GTF for chromosome 22...")
print("   (This may take 1-2 minutes)")

# URL for human GRCh38 release 45, chromosome 22
url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz"

# For class, let's use a smaller file
# In practice, download the full file

# Let's create a realistic sample instead
print("\n✅ Using sample data for demonstration")
print("   (In real analysis, download full GTF)")
```

### Creating a Sample Dataset

```python
# Let's create a more comprehensive sample GTF
# Based on real GENCODE data for chromosome 22

sample_gencode_gtf = """chr22	HAVANA	gene	10736171	10755301	.	-	.	gene_id "ENSG00000100320"; gene_type "protein_coding"; gene_name "RBFOX2";
chr22	HAVANA	transcript	10736171	10755301	.	-	.	gene_id "ENSG00000100320"; transcript_id "ENST00000216146"; gene_name "RBFOX2"; transcript_name "RBFOX2-201"; transcript_type "protein_coding";
chr22	HAVANA	exon	10754774	10755301	.	-	.	gene_id "ENSG00000100320"; transcript_id "ENST00000216146"; exon_number "1"; exon_id "ENSE00003654571";
chr22	HAVANA	CDS	10754774	10755220	.	-	0	gene_id "ENSG00000100320"; transcript_id "ENST00000216146";
chr22	HAVANA	start_codon	10755218	10755220	.	-	0	gene_id "ENSG00000100320"; transcript_id "ENST00000216146";
chr22	HAVANA	exon	10749600	10749854	.	-	.	gene_id "ENSG00000100320"; transcript_id "ENST00000216146"; exon_number "2";
chr22	HAVANA	CDS	10749600	10749854	.	-	0	gene_id "ENSG00000100320"; transcript_id "ENST00000216146";
chr22	HAVANA	exon	10748153	10748267	.	-	.	gene_id "ENSG00000100320"; transcript_id "ENST00000216146"; exon_number "3";
chr22	HAVANA	stop_codon	10748153	10748155	.	-	0	gene_id "ENSG00000100320"; transcript_id "ENST00000216146";
chr22	HAVANA	gene	17662484	17690409	.	+	.	gene_id "ENSG00000073756"; gene_type "protein_coding"; gene_name "PTPN1";
chr22	HAVANA	transcript	17662484	17690409	.	+	.	gene_id "ENSG00000073756"; transcript_id "ENST00000378103"; gene_name "PTPN1"; transcript_name "PTPN1-201"; transcript_type "protein_coding";
chr22	HAVANA	exon	17662484	17662621	.	+	.	gene_id "ENSG00000073756"; transcript_id "ENST00000378103"; exon_number "1";
chr22	HAVANA	start_codon	17662596	17662598	.	+	0	gene_id "ENSG00000073756"; transcript_id "ENST00000378103";
chr22	HAVANA	CDS	17662596	17662621	.	+	0	gene_id "ENSG00000073756"; transcript_id "ENST00000378103";
chr22	HAVANA	exon	17670685	17670815	.	+	.	gene_id "ENSG00000073756"; transcript_id "ENST00000378103"; exon_number "2";
chr22	HAVANA	CDS	17670685	17670815	.	+	2	gene_id "ENSG00000073756"; transcript_id "ENST00000378103";
chr22	HAVANA	gene	23523149	23703817	.	+	.	gene_id "ENSG00000100288"; gene_type "protein_coding"; gene_name "CPT1B";
chr22	HAVANA	transcript	23523149	23703817	.	+	.	gene_id "ENSG00000100288"; transcript_id "ENST00000359453"; gene_name "CPT1B"; transcript_type "protein_coding";
chr22	HAVANA	exon	23523149	23523277	.	+	.	gene_id "ENSG00000100288"; transcript_id "ENST00000359453"; exon_number "1";
"""

with open('gencode_sample.gtf', 'w') as f:
    f.write(sample_gencode_gtf)

print("✅ Created gencode_sample.gtf")
print(f"   Contains: {len(sample_gencode_gtf.split(chr(10)))} lines")
print(f"   Genes: RBFOX2, PTPN1, CPT1B")
```

---

<a id="section4"></a>
## 4️⃣ Creating a DataFrame from GTF

### Basic GTF Parser

```python
def parse_gtf_to_dataframe(gtf_file):
    """
    Parse GTF file into a pandas DataFrame
    
    Returns DataFrame with columns:
    - seqname, source, feature, start, end, score, strand, frame
    - Plus all attributes as separate columns
    """
    
    data = []
    
    with open(gtf_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Split by tab
            fields = line.strip().split('\t')
            
            if len(fields) < 9:
                continue
            
            # Basic fields
            record = {
                'seqname': fields[0],
                'source': fields[1],
                'feature': fields[2],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'score': fields[5],
                'strand': fields[6],
                'frame': fields[7],
                'attribute': fields[8]
            }
            
            # Parse attributes (key "value"; format)
            attributes = fields[8]
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr:
                    # Split on first space
                    parts = attr.split(' ', 1)
                    if len(parts) == 2:
                        key = parts[0]
                        value = parts[1].strip('"')
                        record[key] = value
            
            data.append(record)
    
    return pd.DataFrame(data)

# Parse our sample GTF
df = parse_gtf_to_dataframe('gencode_sample.gtf')

print("✅ GTF parsed into DataFrame!")
print(f"\nShape: {df.shape[0]} rows × {df.shape[1]} columns")
print(f"\nColumns: {list(df.columns)}")
print(f"\nFirst few rows:")
print(df.head())
```

### Exploring the DataFrame

```python
# Let's explore our parsed data

print("="*70)
print("📊 EXPLORING THE GENCODE DATA")
print("="*70)

print(f"\n1. Feature types:")
print(df['feature'].value_counts())

print(f"\n2. Genes in our sample:")
genes = df[df['feature'] == 'gene']
print(genes[['seqname', 'gene_name', 'start', 'end', 'strand']])

print(f"\n3. Transcripts per gene:")
transcripts = df[df['feature'] == 'transcript']
print(transcripts.groupby('gene_name').size())

print(f"\n4. Exons per transcript:")
exons = df[df['feature'] == 'exon']
if 'transcript_id' in exons.columns:
    print(exons.groupby('transcript_id').size())
```

### Filter for Specific Features

```python
# Create separate DataFrames for different features

genes_df = df[df['feature'] == 'gene'].copy()
transcripts_df = df[df['feature'] == 'transcript'].copy()
exons_df = df[df['feature'] == 'exon'].copy()
cds_df = df[df['feature'] == 'CDS'].copy()

print("✅ Created feature-specific DataFrames:")
print(f"   Genes:       {len(genes_df)}")
print(f"   Transcripts: {len(transcripts_df)}")
print(f"   Exons:       {len(exons_df)}")
print(f"   CDS:         {len(cds_df)}")
```

---

<a id="section5"></a>
## 5️⃣ Extracting TSS (Transcription Start Sites)

### What is TSS?

```
TSS = Transcription Start Site
     = Where RNA polymerase begins transcribing

For + strand genes: TSS = start position
For - strand genes: TSS = end position

Example:

+ strand:  ======gene======>
           ^
           TSS (at start)

- strand:  <======gene======
                           ^
                           TSS (at end)
```

### Extract TSS from Transcripts

```python
def extract_tss(transcripts_df):
    """
    Extract TSS from transcript annotations
    
    TSS depends on strand:
    - Plus strand (+): TSS is at start position
    - Minus strand (-): TSS is at end position
    """
    
    tss_data = []
    
    for idx, row in transcripts_df.iterrows():
        # Calculate TSS based on strand
        if row['strand'] == '+':
            tss = row['start']
        else:  # minus strand
            tss = row['end']
        
        tss_data.append({
            'gene_id': row.get('gene_id', ''),
            'gene_name': row.get('gene_name', ''),
            'transcript_id': row.get('transcript_id', ''),
            'chromosome': row['seqname'],
            'tss': tss,
            'strand': row['strand'],
            'gene_start': row['start'],
            'gene_end': row['end']
        })
    
    return pd.DataFrame(tss_data)

# Extract TSS
tss_df = extract_tss(transcripts_df)

print("✅ Extracted TSS for all transcripts!")
print(f"\nTotal TSS: {len(tss_df)}")
print(f"\nTSS data:")
print(tss_df.to_string(index=False))
```

### Visualize TSS

```python
# Let's visualize where TSS is for each gene

fig, ax = plt.subplots(figsize=(12, 4))

for idx, row in tss_df.iterrows():
    gene = row['gene_name']
    start = row['gene_start']
    end = row['gene_end']
    tss = row['tss']
    strand = row['strand']
    
    # Draw gene
    y = idx
    ax.plot([start, end], [y, y], 'b-', linewidth=5, alpha=0.5)
    
    # Mark TSS
    ax.plot(tss, y, 'ro', markersize=10, label='TSS' if idx == 0 else '')
    
    # Add gene name
    ax.text(start - 50000, y, gene, ha='right', va='center')
    
    # Add strand arrow
    arrow_pos = end if strand == '+' else start
    ax.annotate('', xy=(arrow_pos, y), 
                xytext=(arrow_pos - 20000 if strand == '+' else arrow_pos + 20000, y),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))

ax.set_xlabel('Chromosome Position (bp)')
ax.set_ylabel('Gene')
ax.set_title('TSS Locations for Sample Genes')
ax.set_yticks(range(len(tss_df)))
ax.set_yticklabels([])
ax.legend()
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig('tss_visualization.png', dpi=300, bbox_inches='tight')
print("\n✅ Visualization saved: tss_visualization.png")
plt.show()
```

---

<a id="section6"></a>
## 6️⃣ Finding Splice Sites

### What are Splice Sites?

```
Splice sites = Junctions between exons and introns

Pre-mRNA:  exon1 - intron - exon2 - intron - exon3
                  ^        ^        ^        ^
                  splice sites

Donor site (5'): End of exon (usually GT)
Acceptor site (3'): Start of exon (usually AG)
```

### Extract Splice Sites

```python
def extract_splice_sites(exons_df):
    """
    Extract splice sites from exon annotations
    
    For each transcript:
    - Find junctions between consecutive exons
    - Donor = end of upstream exon
    - Acceptor = start of downstream exon
    """
    
    splice_sites = []
    
    # Group by transcript
    for transcript_id, group in exons_df.groupby('transcript_id'):
        # Sort exons by position
        group = group.sort_values('start')
        
        # Get strand
        strand = group.iloc[0]['strand']
        chrom = group.iloc[0]['seqname']
        gene_name = group.iloc[0].get('gene_name', '')
        
        # Iterate through consecutive exons
        exon_list = group.sort_values('exon_number' if 'exon_number' in group.columns else 'start')
        
        for i in range(len(exon_list) - 1):
            exon1 = exon_list.iloc[i]
            exon2 = exon_list.iloc[i + 1]
            
            if strand == '+':
                donor = exon1['end']
                acceptor = exon2['start']
            else:
                donor = exon2['start']
                acceptor = exon1['end']
            
            splice_sites.append({
                'transcript_id': transcript_id,
                'gene_name': gene_name,
                'chromosome': chrom,
                'strand': strand,
                'donor_site': donor,
                'acceptor_site': acceptor,
                'intron_start': min(donor, acceptor),
                'intron_end': max(donor, acceptor),
                'intron_length': abs(acceptor - donor)
            })
    
    return pd.DataFrame(splice_sites)

# Extract splice sites
if len(exons_df) > 0 and 'transcript_id' in exons_df.columns:
    splice_sites_df = extract_splice_sites(exons_df)
    
    print("✅ Extracted splice sites!")
    print(f"\nTotal splice junctions: {len(splice_sites_df)}")
    print(f"\nSplice sites data:")
    print(splice_sites_df.to_string(index=False))
    
    print(f"\n📊 Intron length statistics:")
    print(f"   Mean: {splice_sites_df['intron_length'].mean():.0f} bp")
    print(f"   Min:  {splice_sites_df['intron_length'].min()} bp")
    print(f"   Max:  {splice_sites_df['intron_length'].max()} bp")
else:
    print("⚠️ Not enough exon data for splice site extraction")
```

---

<a id="section7"></a>
## 7️⃣ Extracting Promoter Regions (-1kb to TSS)

### What is a Promoter Region?

```
Promoter = Regulatory region upstream of gene
           Contains transcription factor binding sites

Common definition: -1000 bp to TSS

    -1000bp        TSS      gene
       |------------|=========>
       promoter
```

### Extract Promoter Coordinates

```python
def extract_promoter_regions(tss_df, upstream=1000, downstream=0):
    """
    Extract promoter regions around TSS
    
    Parameters:
    - upstream: bp upstream of TSS (default 1000)
    - downstream: bp downstream of TSS (default 0)
    
    Returns DataFrame with promoter coordinates
    """
    
    promoter_data = []
    
    for idx, row in tss_df.iterrows():
        tss = row['tss']
        strand = row['strand']
        
        if strand == '+':
            # For + strand: promoter is upstream (smaller coordinates)
            prom_start = tss - upstream
            prom_end = tss + downstream
        else:
            # For - strand: promoter is upstream (larger coordinates)
            prom_start = tss - downstream
            prom_end = tss + upstream
        
        # Ensure start < end
        prom_start, prom_end = min(prom_start, prom_end), max(prom_start, prom_end)
        
        # Don't allow negative coordinates
        prom_start = max(1, prom_start)
        
        promoter_data.append({
            'gene_name': row['gene_name'],
            'transcript_id': row['transcript_id'],
            'chromosome': row['chromosome'],
            'strand': row['strand'],
            'tss': tss,
            'promoter_start': prom_start,
            'promoter_end': prom_end,
            'promoter_length': prom_end - prom_start
        })
    
    return pd.DataFrame(promoter_data)

# Extract promoter regions (-1kb to TSS)
promoter_df = extract_promoter_regions(tss_df, upstream=1000, downstream=0)

print("✅ Extracted promoter regions!")
print(f"\nPromoter definition: -1000 bp to TSS")
print(f"\nPromoter data:")
print(promoter_df.to_string(index=False))

# Save promoter coordinates as BED file
bed_file = 'promoter_regions.bed'
with open(bed_file, 'w') as f:
    for idx, row in promoter_df.iterrows():
        f.write(f"{row['chromosome']}\t{row['promoter_start']}\t{row['promoter_end']}\t")
        f.write(f"{row['gene_name']}\t.\t{row['strand']}\n")

print(f"\n✅ Saved promoter coordinates: {bed_file}")
print("   (BED format - can be used with bedtools or genome browsers)")
```

---

<a id="section8"></a>
## 8️⃣ Getting FASTA Sequences

### Understanding FASTA Format

```
FASTA format:
>header (starts with >)
SEQUENCE
SEQUENCE
...

Example:
>chr22:10735171-10736171 RBFOX2 promoter
ATGCGATCGATCGATCG...
```

### Option 1: Using pyfaidx (with real genome)

```python
# This is how you'd extract real sequences from genome FASTA

"""
# First, download genome FASTA (one-time, ~1GB per chromosome)
# Example for chromosome 22:
!wget http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
!gunzip Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

# Then use pyfaidx:
from pyfaidx import Fasta

genome = Fasta('Homo_sapiens.GRCh38.dna.chromosome.22.fa')

def extract_sequence(chrom, start, end, strand):
    # Extract sequence
    seq = genome[chrom][start:end].seq
    
    # Reverse complement if minus strand
    if strand == '-':
        seq = reverse_complement(seq)
    
    return seq

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))
"""

print("📝 For real analysis:")
print("   1. Download genome FASTA from Ensembl")
print("   2. Install pyfaidx: pip install pyfaidx")
print("   3. Use code above to extract sequences")
```

### Option 2: Generate Sample Sequences (for demonstration)

```python
import random

def generate_sample_sequence(length):
    """Generate random DNA sequence for demonstration"""
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choices(bases, k=length))

def create_sample_fasta(promoter_df, output_file='promoter_sequences.fasta'):
    """
    Create FASTA file with sample sequences for promoters
    In real analysis, use actual genome sequences
    """
    
    with open(output_file, 'w') as f:
        for idx, row in promoter_df.iterrows():
            # Create header
            header = f">{row['chromosome']}:{row['promoter_start']}-{row['promoter_end']}"
            header += f" {row['gene_name']} {row['strand']}"
            
            # Generate sample sequence
            seq_length = row['promoter_end'] - row['promoter_start']
            sequence = generate_sample_sequence(seq_length)
            
            # Write to file
            f.write(header + '\n')
            # Write sequence in 60 bp chunks (standard FASTA format)
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')
    
    print(f"✅ Created FASTA file: {output_file}")
    print(f"   Contains {len(promoter_df)} sequences")
    return output_file

# Create sample FASTA
fasta_file = create_sample_fasta(promoter_df)

# Read and display first sequence
from Bio import SeqIO

print("\n📖 First sequence in FASTA file:")
for record in SeqIO.parse(fasta_file, 'fasta'):
    print(f"\nHeader: {record.id} {record.description}")
    print(f"Length: {len(record.seq)} bp")
    print(f"Sequence preview: {str(record.seq)[:80]}...")
    break  # Just show first one
```

---

<a id="section9"></a>
## 9️⃣ Finding p53 Binding Sites with Regular Expression

### Understanding p53 Motif

```
p53 binds to DNA with this pattern:
RRRCWWGYYY {spacer 0-10 bp} RRRCWWGYYY

Where:
R = Purine (A or G)
W = Weak bond (A or T)
Y = Pyrimidine (C or T)

Example matches:
GGACATGCCCGGGCATGTCC        ← no spacer
AGACATGCCCATGGGGCATGTCT     ← 2 bp spacer
GGGCATGTCCATAGGGACTTGCCT    ← 5 bp spacer
```

### Build the Regular Expression

```python
# Step-by-step regex construction

print("🔨 Building p53 motif regular expression:\n")

# Define IUPAC ambiguity codes
R = '[AG]'      # Purine
W = '[AT]'      # Weak bond
Y = '[CT]'      # Pyrimidine

print(f"R (purine) = {R}")
print(f"W (weak)   = {W}")
print(f"Y (pyrimidine) = {Y}")

# Build decamer pattern: RRRCWWGYYY
decamer = f"{R}{R}{R}C{W}{W}G{Y}{Y}{Y}"
print(f"\nDecamer pattern: {decamer}")
print("   Matches: AAACAAGCCC, GGACTTGCTT, etc.")

# Add spacer (0 to 10 bp)
spacer = ".{0,10}"
print(f"\nSpacer: {spacer}")
print("   Allows 0 to 10 nucleotides between decamers")

# Complete p53 motif
p53_pattern = f"{decamer}{spacer}{decamer}"

print(f"\n✅ Complete p53 regex pattern:")
print(f"   {p53_pattern}")
```

### Test the Pattern

```python
# Let's test with known p53 binding sites

test_sequences = [
    ("GGACATGCCCGGGCATGTCC", "Classic p53 site (no spacer)"),
    ("AGACATGCCCATGGGGCATGTCT", "With 2 bp spacer"),
    ("GGGCATGTCCATAGGGACTTGCCT", "With 5 bp spacer"),
    ("ATCGATCGATCG", "Not a p53 site"),
]

print("\n🧪 Testing p53 pattern on known sequences:\n")

for seq, description in test_sequences:
    match = re.search(p53_pattern, seq)
    if match:
        print(f"✓ MATCH: {description}")
        print(f"  Sequence: {seq}")
        print(f"  Found at position {match.start()}: {match.group()}")
        spacer_len = len(match.group()) - 20
        print(f"  Spacer length: {spacer_len} bp\n")
    else:
        print(f"✗ NO MATCH: {description}")
        print(f"  Sequence: {seq}\n")
```

### Search Promoter Sequences for p53 Sites

```python
def find_p53_sites_in_sequences(fasta_file, pattern):
    """
    Search all sequences in FASTA for p53 binding sites
    """
    
    results = []
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequence = str(record.seq).upper()
        
        # Find all matches
        for match in re.finditer(pattern, sequence):
            results.append({
                'sequence_id': record.id,
                'description': record.description,
                'position': match.start(),
                'p53_site': match.group(),
                'site_length': len(match.group()),
                'spacer_length': len(match.group()) - 20
            })
    
    return pd.DataFrame(results)

# Search for p53 sites
print("🔍 Searching for p53 binding sites in promoter sequences...\n")

# First, let's add some known p53 sites to our sequences
# (Since we generated random sequences, they won't have real sites)

# Let's create a new FASTA with embedded p53 sites
def create_fasta_with_p53_sites(promoter_df, output_file='promoters_with_p53.fasta'):
    """Create FASTA with embedded known p53 sites"""
    
    known_p53_sites = [
        "GGACATGCCCGGGCATGTCC",
        "AGACATGCCCATGGGGCATGTCT",
        "GGGCATGTCCATAGGGACTTGCCT",
        "AGGCAAGTTCAGGACAAGCCT",
    ]
    
    with open(output_file, 'w') as f:
        for idx, row in promoter_df.iterrows():
            header = f">{row['chromosome']}:{row['promoter_start']}-{row['promoter_end']}"
            header += f" {row['gene_name']} {row['strand']}"
            
            # Generate random sequence
            seq_length = row['promoter_end'] - row['promoter_start']
            sequence = generate_sample_sequence(seq_length)
            
            # Insert p53 site at random position
            if random.random() > 0.3:  # 70% chance of having a site
                site = random.choice(known_p53_sites)
                insert_pos = random.randint(100, len(sequence) - len(site) - 100)
                sequence = sequence[:insert_pos] + site + sequence[insert_pos+len(site):]
            
            f.write(header + '\n')
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')
    
    return output_file

# Create FASTA with p53 sites
fasta_with_p53 = create_fasta_with_p53_sites(promoter_df)
print(f"✅ Created: {fasta_with_p53}")

# Now search
p53_results = find_p53_sites_in_sequences(fasta_with_p53, p53_pattern)

print(f"\n✅ Search complete!")
print(f"   Found {len(p53_results)} p53 binding sites")
print(f"   In {p53_results['sequence_id'].nunique()} promoter sequences")

if len(p53_results) > 0:
    print(f"\n📊 Results:")
    print(p53_results.to_string(index=False))
    
    # Spacer distribution
    print(f"\n📏 Spacer length distribution:")
    print(p53_results['spacer_length'].value_counts().sort_index())
```

### Visualize Results

```python
if len(p53_results) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Number of sites per gene
    sites_per_gene = p53_results.groupby('sequence_id').size()
    axes[0].bar(range(len(sites_per_gene)), sites_per_gene.values, color='steelblue', alpha=0.7)
    axes[0].set_xlabel('Gene (index)', fontsize=12)
    axes[0].set_ylabel('Number of p53 Sites', fontsize=12)
    axes[0].set_title('p53 Sites per Promoter', fontsize=14, fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3)
    
    # Plot 2: Spacer length distribution
    spacer_counts = p53_results['spacer_length'].value_counts().sort_index()
    axes[1].bar(spacer_counts.index, spacer_counts.values, color='coral', alpha=0.7)
    axes[1].set_xlabel('Spacer Length (bp)', fontsize=12)
    axes[1].set_ylabel('Number of Sites', fontsize=12)
    axes[1].set_title('Spacer Length Distribution', fontsize=14, fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('p53_results.png', dpi=300, bbox_inches='tight')
    print("\n✅ Visualization saved: p53_results.png")
    plt.show()
else:
    print("⚠️ No p53 sites found (try running again - random generation)")
```

---

<a id="section10"></a>
## 🎯 Real-World Example: Complete Pipeline

### Putting It All Together

```python
print("="*70)
print("COMPLETE PIPELINE: GENCODE → p53 BINDING SITES")
print("="*70)

# Step 1: Parse GTF
print("\n[1/7] Parsing GENCODE GTF...")
df = parse_gtf_to_dataframe('gencode_sample.gtf')
print(f"      ✓ Parsed {len(df)} features")

# Step 2: Extract transcripts
print("\n[2/7] Extracting transcripts...")
transcripts_df = df[df['feature'] == 'transcript'].copy()
print(f"      ✓ Found {len(transcripts_df)} transcripts")

# Step 3: Extract TSS
print("\n[3/7] Calculating TSS positions...")
tss_df = extract_tss(transcripts_df)
print(f"      ✓ Extracted {len(tss_df)} TSS")

# Step 4: Define promoter regions
print("\n[4/7] Defining promoter regions (-1kb)...")
promoter_df = extract_promoter_regions(tss_df, upstream=1000)
print(f"      ✓ Defined {len(promoter_df)} promoter regions")

# Step 5: Get sequences
print("\n[5/7] Extracting promoter sequences...")
fasta_file = create_fasta_with_p53_sites(promoter_df, 'final_promoters.fasta')
print(f"      ✓ Created FASTA: {fasta_file}")

# Step 6: Build p53 pattern
print("\n[6/7] Building p53 motif pattern...")
R = '[AG]'
W = '[AT]'
Y = '[CT]'
p53_pattern = f"{R}{R}{R}C{W}{W}G{Y}{Y}{Y}.{{0,10}}{R}{R}{R}C{W}{W}G{Y}{Y}{Y}"
print(f"      ✓ Pattern: {p53_pattern}")

# Step 7: Search for p53 sites
print("\n[7/7] Searching for p53 binding sites...")
results = find_p53_sites_in_sequences(fasta_file, p53_pattern)
print(f"      ✓ Found {len(results)} p53 sites")

# Summary
print("\n" + "="*70)
print("📊 FINAL SUMMARY")
print("="*70)
print(f"Input:    GENCODE GTF annotation")
print(f"Genes:    {len(df[df['feature'] == 'gene'])}")
print(f"Promoters analyzed: {len(promoter_df)}")
print(f"p53 sites found:    {len(results)}")
print(f"Genes with p53:     {results['sequence_id'].nunique() if len(results) > 0 else 0}")

if len(results) > 0:
    print(f"\nTop genes with p53 sites:")
    gene_counts = results.groupby('description').size().sort_values(ascending=False).head()
    for gene, count in gene_counts.items():
        gene_name = gene.split()[0] if ' ' in gene else gene
        print(f"  • {gene_name}: {count} site(s)")

# Save final results
results.to_csv('p53_binding_sites_final.csv', index=False)
print(f"\n✅ Results saved: p53_binding_sites_final.csv")

print("\n" + "="*70)
print("🎉 PIPELINE COMPLETE!")
print("="*70)
```

---

## 📚 Summary and Key Takeaways

### What We Learned Today

```python
summary = {
    '1. GENCODE': 'Gold standard genome annotation',
    '2. GTF Format': '9 columns: chr, source, feature, start, end, score, strand, frame, attributes',
    '3. Parsing': 'Convert GTF to DataFrame for analysis',
    '4. TSS': 'Transcription start site - depends on strand direction',
    '5. Splice Sites': 'Junctions between exons (donor/acceptor)',
    '6. Promoters': 'Regulatory regions upstream of TSS (-1kb)',
    '7. FASTA': 'Sequence format for extracting DNA',
    '8. p53 Motif': 'RRRCWWGYYY {0-10bp} RRRCWWGYYY',
    '9. Regex': 'Pattern matching for motif finding'
}

print("=" * 70)
print("📖 KEY CONCEPTS")
print("=" * 70)
for topic, description in summary.items():
    print(f"{topic:20s} → {description}")
```

### Files Generated

```python
import os

print("\n" + "=" * 70)
print("📂 FILES CREATED")
print("=" * 70)

files = [
    'gencode_sample.gtf',
    'promoter_regions.bed',
    'promoter_sequences.fasta',
    'promoters_with_p53.fasta',
    'final_promoters.fasta',
    'p53_binding_sites_final.csv',
    'tss_visualization.png',
    'p53_results.png'
]

for f in files:
    if os.path.exists(f):
        size = os.path.getsize(f)
        print(f"✓ {f:40s} ({size:,} bytes)")
```

### Next Steps / Homework

```python
homework = """
🏠 OPTIONAL HOMEWORK:

1. Download real GENCODE data:
   - Go to https://www.gencodegenes.org/
   - Download full GTF for GRCh38
   - Analyze chromosome 22

2. Try different motifs:
   - NF-κB: GGGRNWYYCC
   - STAT3: TTCN2-4GAA
   - SP1: GGGGCGGGG

3. Compare + and - strands:
   - Search both DNA strands
   - Use reverse complement

4. Validate findings:
   - Check against ChIP-seq data
   - Use JASPAR database
   - Compare with known p53 targets

5. Advanced:
   - Calculate motif enrichment
   - Compare to random sequences
   - Build position weight matrix
"""

print(homework)
```

---

## 🎓 Class Exercise

**Try this now!**

```python
# EXERCISE: Find NF-κB binding sites
# NF-κB motif: GGGRNWYYCC
# R=[AG], W=[AT], Y=[CT], N=[ATGC]

print("="*70)
print("👨‍💻 CLASS EXERCISE: Find NF-κB Binding Sites")
print("="*70)

# TODO: Build the NF-κB regex pattern
# Hint: Similar to p53, but simpler (single motif, no spacer)

nfkb_pattern = ""  # YOUR CODE HERE

# Test it
test_seq = "GGGAGATTCC"  # Should match!

if nfkb_pattern:
    match = re.search(nfkb_pattern, test_seq)
    if match:
        print(f"✅ Correct! Found: {match.group()}")
    else:
        print("❌ Pattern doesn't match test sequence. Try again!")
else:
    print("💡 Hint: Start with [GGG]...")
```

---

## 🎉 Congratulations!

You've learned how to:
- ✅ Understand GENCODE annotations
- ✅ Parse GTF files
- ✅ Extract genomic coordinates (TSS, splice sites, promoters)
- ✅ Work with FASTA sequences
- ✅ Find DNA motifs using regular expressions
- ✅ Build a complete bioinformatics pipeline

**This is real bioinformatics work!** 🧬🔬

---

## 📧 Questions?

Post on Ed Discussion or come to office hours!

**Next class:** [Next topic]

---

**END OF NOTEBOOK**
