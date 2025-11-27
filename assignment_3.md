# Assignment 3: Bioinformatics File Formats and Analysis
**BMI 503 – Introduction to Programming for Biomedical Informatics**

**Instructors:** Prof. Ramana Davuluri, Prof. Fusheng Wang, and Rekha Sathian  
**Institution:** Department of Biomedical Informatics, Stony Brook University  
**Total Points:** 15 points

---

## 📋 Overview
This assignment introduces you to common bioinformatics file formats (FASTA, FASTQ, VCF) and teaches you how to parse and analyze genomic data using Python libraries including Biopython, PyVCF, and pandas.

---

## 📁 Datasets

All datasets are available in the course GitHub repository: `https://github.com/duttaprat/BMI_503/tree/main/BMI503-Data`

**Dataset 1: Gene Sequences (FASTA)**
- File: `gene_sequences.fasta`
- Contains: TP53 and BRCA1 gene sequences
- URL: `https://raw.githubusercontent.com/duttaprat/BMI_503/main/BMI503-Data/gene_sequences.fasta`

**Dataset 2: Sequencing Reads (FASTQ)**
- File: `sample_reads.fastq`
- Contains: Illumina sequencing reads with quality scores
- URL: `https://raw.githubusercontent.com/duttaprat/BMI_503/main/BMI503-Data/sample_reads.fastq`

**Dataset 3: Genetic Variants (VCF)**
- File: `variants.vcf`
- Contains: Somatic mutations from cancer samples
- URL: `https://raw.githubusercontent.com/duttaprat/BMI_503/main/BMI503-Data/variants.vcf`

---

## 📝 Assignment Tasks

### **Part I: Working with FASTA Files (5 points)**

FASTA format stores biological sequences with a header line (starting with '>') followed by sequence data.

#### Task 1.1: Parse FASTA with Biopython (2 points)

Write code to parse `gene_sequences.fasta` and extract basic information:

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# TODO: Parse the FASTA file
for record in SeqIO.parse("gene_sequences.fasta", "fasta"):
    # Extract and print:
    # - Sequence ID
    # - Sequence length
    # - GC content (%)
    pass
```

**Expected Output for each sequence:**
- Sequence ID (gene name)
- Total length in base pairs
- GC content as percentage

#### Task 1.2: Sequence Analysis (3 points)

For each gene sequence, perform the following analyses:

1. **Find the start codon (ATG)** - Report the position of the first ATG
2. **Count total codons** - How many complete 3-base codons are present?
3. **Translate to protein** - Use Biopython's translate function to convert DNA to amino acid sequence

**Create a visualization:**
- Bar chart comparing the lengths of TP53 and BRCA1
- Save as `sequence_comparison.png`

```python
import matplotlib.pyplot as plt

# Hint: Store lengths in a dictionary or list
genes = ['TP53', 'BRCA1']
lengths = [len1, len2]  # Replace with actual lengths

plt.bar(genes, lengths)
plt.xlabel('Gene')
plt.ylabel('Sequence Length (bp)')
plt.title('Gene Sequence Length Comparison')
plt.savefig('sequence_comparison.png', dpi=300, bbox_inches='tight')
```

---

### **Part II: Working with FASTQ Files (4 points)**

FASTQ format stores sequences with quality scores. Each read consists of 4 lines:
1. Header (starts with @)
2. Sequence
3. Separator (+)
4. Quality scores

#### Task 2.1: Parse FASTQ and Calculate Statistics (2 points)

```python
from Bio import SeqIO

# Parse FASTQ file
reads = list(SeqIO.parse("sample_reads.fastq", "fastq"))

# Calculate:
# 1. Total number of reads
# 2. Average read length
# 3. Average quality score per read
```

**Quality Score Conversion:**
Phred quality score Q = ASCII value - 33

```python
def calculate_average_quality(quality_string):
    """Convert ASCII quality scores to Phred scores and calculate average"""
    phred_scores = [ord(q) - 33 for q in quality_string]
    return sum(phred_scores) / len(phred_scores)
```

#### Task 2.2: Quality Control Analysis (2 points)

Analyze the quality of sequencing reads:

1. Count reads with average quality < Q20 (low quality)
2. Count reads with average quality ≥ Q30 (high quality)
3. Calculate the percentage of high-quality reads

**Create a visualization:**
- Histogram showing distribution of average quality scores across all reads
- Save as `quality_distribution.png`

```python
import matplotlib.pyplot as plt

# Collect average quality scores
avg_qualities = [calculate_average_quality(str(read.letter_annotations["phred_quality"])) 
                 for read in reads]

plt.hist(avg_qualities, bins=20, edgecolor='black')
plt.xlabel('Average Quality Score (Phred)')
plt.ylabel('Number of Reads')
plt.title('Distribution of Read Quality Scores')
plt.axvline(x=20, color='orange', linestyle='--', label='Q20 threshold')
plt.axvline(x=30, color='green', linestyle='--', label='Q30 threshold')
plt.legend()
plt.savefig('quality_distribution.png', dpi=300, bbox_inches='tight')
```

---

### **Part III: Working with VCF Files (6 points)**

VCF (Variant Call Format) stores information about genetic variants. It contains header lines (starting with ##) followed by a table with columns including CHROM, POS, REF, ALT, QUAL, FILTER, and INFO.

#### Task 3.1: Parse VCF with pandas (2 points)

```python
import pandas as pd

# Load VCF file (skip header lines starting with ##)
vcf_data = pd.read_csv("variants.vcf", sep="\t", comment="##", 
                       header=0)  # First non-## line is the header

# Calculate basic statistics:
# 1. Total number of variants
# 2. Number of variants per chromosome
# 3. Which chromosome has the most variants?
```

**Report:**
- Total variants in the file
- Variant count for each chromosome
- Chromosome with highest variant density

#### Task 3.2: Quality Filtering (2 points)

Apply standard quality control filters:

```python
# Filter 1: Quality score ≥ 30
high_qual_variants = vcf_data[vcf_data['QUAL'] >= 30]

# Filter 2: FILTER column = "PASS"
passed_variants = high_qual_variants[high_qual_variants['FILTER'] == 'PASS']

# Calculate percentage of variants passing both filters
pass_rate = (len(passed_variants) / len(vcf_data)) * 100
```

**Report:**
- Number of variants passing quality threshold
- Number of variants with FILTER="PASS"
- Overall pass rate as percentage

#### Task 3.3: Cancer Gene Analysis (2 points)

Focus on clinically significant cancer genes: **TP53, BRCA1, BRCA2, EGFR, KRAS**

The INFO column contains gene names and impact predictions. Parse this information:

```python
def extract_gene_from_info(info_string):
    """Extract gene name from INFO field"""
    # INFO format example: "GENE=TP53;IMPACT=HIGH;..."
    for item in info_string.split(';'):
        if item.startswith('GENE='):
            return item.split('=')[1]
    return None

def extract_impact_from_info(info_string):
    """Extract impact level from INFO field"""
    for item in info_string.split(';'):
        if item.startswith('IMPACT='):
            return item.split('=')[1]
    return None

# Add parsed columns
vcf_data['GENE'] = vcf_data['INFO'].apply(extract_gene_from_info)
vcf_data['IMPACT'] = vcf_data['INFO'].apply(extract_impact_from_info)

# Filter for cancer genes
cancer_genes = ['TP53', 'BRCA1', 'BRCA2', 'EGFR', 'KRAS']
cancer_variants = vcf_data[vcf_data['GENE'].isin(cancer_genes)]
```

**Analysis requirements:**
1. Count variants per gene
2. Count variants by IMPACT level (HIGH, MODERATE, LOW, MODIFIER)
3. Which gene has the most HIGH impact variants?

**Create visualizations:**

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Plot 1: Variants per gene
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

gene_counts = cancer_variants['GENE'].value_counts()
axes[0].bar(gene_counts.index, gene_counts.values)
axes[0].set_xlabel('Gene')
axes[0].set_ylabel('Number of Variants')
axes[0].set_title('Variants per Cancer Gene')
axes[0].tick_params(axis='x', rotation=45)

# Plot 2: Variants by impact
impact_counts = cancer_variants['IMPACT'].value_counts()
axes[1].bar(impact_counts.index, impact_counts.values, color=['red', 'orange', 'yellow', 'lightblue'])
axes[1].set_xlabel('Impact Level')
axes[1].set_ylabel('Number of Variants')
axes[1].set_title('Variants by Predicted Impact')

plt.tight_layout()
plt.savefig('cancer_variant_analysis.png', dpi=300, bbox_inches='tight')
```

---

## 💻 Setup and Data Download

**Step 1: Install required packages**
```python
# Run in Jupyter notebook or terminal
!pip install biopython pandas matplotlib seaborn
```

**Step 2: Download datasets**
```python
import urllib.request

base_url = "https://raw.githubusercontent.com/duttaprat/BMI_503/main/BMI503-Data/"

files = [
    "gene_sequences.fasta",
    "sample_reads.fastq",
    "variants.vcf"
]

for file in files:
    url = base_url + file
    urllib.request.urlretrieve(url, file)
    print(f"Downloaded: {file}")
```

**Alternative: Using wget (Linux/Mac)**
```bash
wget https://raw.githubusercontent.com/duttaprat/BMI_503/main/BMI503-Data/gene_sequences.fasta
wget https://raw.githubusercontent.com/duttaprat/BMI_503/main/BMI503-Data/sample_reads.fastq
wget https://raw.githubusercontent.com/duttaprat/BMI_503/main/BMI503-Data/variants.vcf
```

---

## 📤 Submission Requirements

**Due Date:** [To be announced on Brightspace]

**Submit a ZIP file named: `Assignment3_FirstnameLastname.zip`**

Your submission should contain:

```
Assignment3_FirstnameLastname/
├── assignment3.ipynb          # Jupyter notebook with all code
├── sequence_comparison.png    # Part I visualization
├── quality_distribution.png   # Part II visualization
├── cancer_variant_analysis.png # Part III visualization
└── RESULTS.txt                # Summary of your findings
```

**RESULTS.txt format:**
```
Student Name: [Your Full Name]
Student ID: [Your SBU ID]
Date: [Submission Date]

==========================================
PART I: FASTA ANALYSIS
==========================================
TP53:
  - Length: [X] bp
  - GC Content: [X.XX]%
  - First ATG position: [X]
  - Number of codons: [X]
  - Protein length: [X] amino acids

BRCA1:
  - Length: [X] bp
  - GC Content: [X.XX]%
  - First ATG position: [X]
  - Number of codons: [X]
  - Protein length: [X] amino acids

==========================================
PART II: FASTQ QUALITY CONTROL
==========================================
Total reads: [X]
Average read length: [X.X] bp
Overall average quality: Q[XX.X]

Quality Distribution:
  - Low quality reads (Q < 20): [X] ([XX.X]%)
  - High quality reads (Q ≥ 30): [X] ([XX.X]%)

==========================================
PART III: VCF VARIANT ANALYSIS
==========================================
Total variants: [X]
High-quality variants (QUAL ≥ 30, FILTER=PASS): [X] ([XX.X]%)

Chromosome with most variants: chr[X] ([X] variants)

Cancer Gene Analysis:
  - TP53: [X] variants
  - BRCA1: [X] variants
  - BRCA2: [X] variants
  - EGFR: [X] variants
  - KRAS: [X] variants

Impact Distribution:
  - HIGH: [X] variants
  - MODERATE: [X] variants
  - LOW: [X] variants
  - MODIFIER: [X] variants

Gene with most HIGH impact variants: [GENE] ([X] variants)

==========================================
TIME SPENT: [X] hours
==========================================
```

---

## 🎯 Grading Rubric (15 points total)

| Component | Points | Criteria |
|-----------|--------|----------|
| **Part I: FASTA** | | |
| Task 1.1: Parsing & GC content | 2 | Correct parsing, accurate GC calculation |
| Task 1.2: Sequence analysis | 3 | Start codon, codons, translation, visualization |
| **Part II: FASTQ** | | |
| Task 2.1: Basic statistics | 2 | Correct read count, length, quality calculation |
| Task 2.2: Quality control | 2 | Accurate filtering, proper visualization |
| **Part III: VCF** | | |
| Task 3.1: VCF parsing | 2 | Correct loading, chromosome statistics |
| Task 3.2: Quality filtering | 2 | Proper filter implementation |
| Task 3.3: Gene analysis | 2 | INFO parsing, gene/impact analysis, visualizations |
| **Code Quality** | | |
| Code organization | 1 | Clean, commented, follows best practices |
| Documentation | 1 | Complete RESULTS.txt with accurate values |
| **Total** | **15** | |

**Deductions:**
- Missing or incorrect visualization: -0.5 points each
- Code does not run: -2 points
- Missing RESULTS.txt: -1 point
- Late submission: -1 point per day

---

## 💡 Helpful Resources

**Python Libraries:**
- **Biopython documentation:** https://biopython.org/wiki/Documentation
- **Biopython Tutorial (PDF):** http://biopython.org/DIST/docs/tutorial/Tutorial.pdf
- **Pandas guide:** https://pandas.pydata.org/docs/user_guide/index.html

**File Format Specifications:**
- **FASTA format:** https://zhanggroup.org/FASTA/
- **FASTQ format:** https://en.wikipedia.org/wiki/FASTQ_format
- **VCF format:** https://samtools.github.io/hts-specs/VCFv4.2.pd