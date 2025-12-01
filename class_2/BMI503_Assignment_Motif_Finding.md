# BMI 503 - Assignment: Motif Finding in Human Genome
## Finding Regulatory Elements Using Regular Expressions

**Due Dates:**
- Assignment 3 (p53 Sites)



---

## 📚 Learning Objectives

By completing these assignments, you will:
- Apply regular expressions to real genomic sequences
- Understand splice site consensus sequences
- Discover transcription factor binding sites
- Connect computational patterns to biological function
- Practice the complete bioinformatics workflow: download → parse → extract → search → analyze

---

## 🔬 Background

In class, we learned how to:
1. Parse GTF files from GENCODE
2. Extract genomic coordinates (TSS, splice sites, promoters)
3. Retrieve DNA sequences using pyfaidx
4. Build regular expressions for DNA motifs
5. Search sequences for patterns

Now you'll apply this workflow to find **real regulatory elements** in the human genome!

---

---

## 🧬 ASSIGNMENT 3: p53 Binding Sites 

### Background

**What is p53?**
- Tumor suppressor protein (TP53 gene)
- Called the "**Guardian of the Genome**"
- Activates genes for DNA repair, cell cycle arrest, apoptosis
- **Mutated in ~50% of all human cancers!**

**p53 Binding Site Pattern:**
```
RRRCWWGYYY {spacer 0-10 bp} RRRCWWGYYY

Two "half-sites" (decamers):
- Each half-site: RRRCWWGYYY
- Separated by 0-10 bp spacer (variable length)
- p53 binds as a tetramer (4 proteins)

Where:
R = Purine (A or G)
W = Weak bond (A or T)
Y = Pyrimidine (C or T)

Example matches:
GGACATGCCCGGGCATGTCC          ← 0 bp spacer
AGACATGCCCATGGGGCATGTCT       ← 2 bp spacer
GGGCATGTCCATAGGGACTTGCCT      ← 5 bp spacer
```

**Why this matters:** Understanding where p53 binds helps us:
- Identify p53 target genes
- Understand cancer biology
- Design therapies

---

### Your Tasks

#### Task 3.1: Build the Complex Regular Expression 

This is the **most challenging** pattern!

**Step-by-step approach:**

1. Build one half-site first:
   - R (purine) = `[AG]`
   - W (weak) = `[AT]`
   - Y (pyrimidine) = `[CT]`
   - Pattern: `RRRCWWGYYY` = `[AG][AG][AG]C[AT][AT]G[CT][CT][CT]`

2. Add the spacer:
   - 0-10 bp spacer = `.{0,10}`
   - `.` means any character
   - `{0,10}` means 0 to 10 times

3. Add the second half-site (same as first)

**What to submit:**
```python
# Build the pattern:
half_site = r"[__][__][__]C[__][__]G[__][__][__]"
spacer = r".{_,__}"
p53_pattern = half_site + spacer + half_site

# Test your pattern on these known p53 sites:
test_sequences = [
    "GGACATGCCCGGGCATGTCC",        # Should match
    "AGACATGCCCATGGGGCATGTCT",     # Should match
    "GGGCATGTCCATAGGGACTTGCCT",    # Should match
    "ATCGATCGATCGATCG",             # Should NOT match
]

# Verify each one matches or doesn't match as expected
```

---

#### Task 3.2: Extract Promoter Sequences 

You already have `promoter_df` from class with promoter coordinates (-1kb to TSS).

**Step-by-step:**
1. For each promoter, extract the 1000 bp sequence
2. Handle strand correctly (reverse complement for minus strand)
3. Store: gene_name, sequence, chromosome, start, end, strand

**What to submit:**
- Code outline for extraction
- First 3 promoter sequences (just first 100 bp of each)

---

#### Task 3.3: Search for p53 Binding Sites 

Search each promoter sequence for p53 motifs.

**Hints:**
```python
import re

results = []
for each promoter:
    matches = re.finditer(p53_pattern, sequence)
    for match in matches:
        # Record: gene_name, position, matched_sequence, spacer_length
        # Spacer length = total_match_length - 20 (since each half-site is 10 bp)
```

**What to submit:**
- Code showing search logic
- Example of how you calculate spacer length
- List of first 10 matches found

---

#### Task 3.4: Analyze Spacer Length Distribution 

**Calculate:**
- How many matches have 0 bp spacer?
- How many have 1 bp? 2 bp? ... 10 bp?
- Which spacer length is most common?

**Create a histogram:**
```
Spacer Length Distribution
==========================
0 bp:  ████████ (23 sites)
1 bp:  ██████ (18 sites)
2 bp:  ████████████ (34 sites)
3 bp:  ██████████ (28 sites)
4 bp:  ████████████████ (42 sites)
5 bp:  ██████████████ (38 sites)
6 bp:  ████████ (22 sites)
7 bp:  ██████ (16 sites)
8 bp:  ████ (11 sites)
9 bp:  ██ (7 sites)
10 bp: ██ (5 sites)
```

You can create this using matplotlib or just ASCII art like above.

**What to submit:**
- Spacer distribution (table or plot)
- Which spacer length is most common?

---

#### Task 3.5: Comprehensive Report 

**Create a final report:**

```
p53 BINDING SITE ANALYSIS
=========================

Dataset: Chromosome 22, GENCODE v49
Region: Promoters (-1000 bp to TSS)

SUMMARY STATISTICS:
- Total promoters analyzed: _____
- Promoters with p53 sites: _____
- Percentage with sites: ____%
- Total p53 sites found: _____
- Average sites per promoter: ____

TOP 10 GENES WITH MOST p53 SITES:
1. ________ : __ sites
2. ________ : __ sites
3. ________ : __ sites
... (continue to 10)

SPACER LENGTH ANALYSIS:
- Most common spacer: __ bp (__ sites)
- Least common spacer: __ bp (__ sites)
- Average spacer length: __ bp

EXAMPLE p53 SITES:
Gene: ________ | Position: chr22:_______ | Strand: _ | Spacer: __ bp
Sequence: _________________________________

Gene: ________ | Position: chr22:_______ | Strand: _ | Spacer: __ bp
Sequence: _________________________________

(Show 3-5 examples)
```

---


---

## 📚 Resources

**In-Class Notebook:**
- GENCODE GTF parsing
- TSS and splice site extraction
- Promoter region definition
- pyfaidx sequence extraction
- Regular expression examples

**Python Documentation:**
- `re` module: https://docs.python.org/3/library/re.html
- `pandas` DataFrames: https://pandas.pydata.org/docs/

**Biology References:**
- GENCODE: https://www.gencodegenes.org/
- p53 Database: https://p53.fr/
- GeneCards: https://www.genecards.org/

---



**Good luck! 🧬🔬**




