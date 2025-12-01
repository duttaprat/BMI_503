# BMI 503 - Assignment: Motif Finding in Human Genome
## Finding Regulatory Elements Using Regular Expressions

**Due Dates:**
- Assignment 1 (Acceptor Sites): Wednesday, December 4, 2024
- Assignment 2 (Donor Sites): Friday, December 6, 2024  
- Assignment 3 (p53 Sites): Monday, December 9, 2024

**Total Points:** 100 (Acceptor: 25pts, Donor: 30pts, p53: 45pts)

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

## 🧬 ASSIGNMENT 1: Splice Acceptor Sites (25 points)

### Background

**What are splice acceptor sites?**
- Mark where **introns end** and **exons begin**
- Recognized by the spliceosome during pre-mRNA splicing
- Critical for proper gene expression

**Consensus Pattern: Y₁₀NCAG**
```
Y₁₀ = 10 pyrimidines (C or T)
N   = Any nucleotide (A, T, G, or C)
CAG = Nearly invariant (the splice point is before C)

Example sequences that match:
CCTTCTCTCTNCAG
TTTTTTTTTTACAG
CTCTCTCTCTGCAG
```

**Why this matters:** Mutations in splice sites cause ~15% of genetic diseases!

---

### Your Tasks

#### Task 1.1: Build the Regular Expression (5 pts)

Convert the consensus pattern `Y₁₀NCAG` into a Python regular expression.

**Hints:**
- Y (pyrimidine) = `[CT]`
- N (any base) = `[ATGC]`
- Use `{n}` for exact repetitions

**What to submit:**
```python
# Your pattern (fill in the blanks):
acceptor_pattern = r"[CT]{__}[____]CAG"
```

---

#### Task 1.2: Extract Sequences Around Acceptor Sites (8 pts)

From your class notebook, you already have `splice_sites_df` with donor and acceptor coordinates.

**Step-by-step approach:**
1. Filter for acceptor sites only
2. For each acceptor site, extract sequence from `(acceptor_position - 10)` to `(acceptor_position + 10)`
3. Handle strand correctly (use reverse complement for minus strand)
4. Store sequences in a list or DataFrame

**Hints:**
```python
# You already have the genome loaded:
genome = Fasta('Homo_sapiens.GRCh38.dna.chromosome.22.fa')

# For each acceptor site:
# - Get coordinates
# - Extract 21 bp sequence (10 bp upstream, acceptor site, 10 bp downstream)
# - Apply reverse complement if on minus strand
```

**What to submit:**
- Code showing extraction (can be pseudocode or outline)
- First 5 sequences as examples

---

#### Task 1.3: Search for Pattern Matches (7 pts)

Search each extracted sequence for your regex pattern.

**Hints:**
```python
import re

# For each sequence:
#   if re.search(pattern, sequence):
#       # Count as match
```

**What to submit:**
- Code showing how you search
- Count of total matches

---

#### Task 1.4: Calculate Statistics & Report (5 pts)

**Calculate:**
- Total acceptor sites analyzed
- Number matching consensus
- Percentage match rate
- Show 3 examples of matches and 2 examples of non-matches

**What to submit:**

Create a report with:

```
SPLICE ACCEPTOR SITE ANALYSIS
=============================

Dataset: Chromosome 22, GENCODE v49

Results:
- Total acceptor sites: ___
- Sites matching Y₁₀NCAG: ___
- Match rate: ___%

Example Matches:
1. Position chr22:_______ | Sequence: ____________ | Strand: _
2. Position chr22:_______ | Sequence: ____________ | Strand: _
3. Position chr22:_______ | Sequence: ____________ | Strand: _

Example Non-Matches:
1. Position chr22:_______ | Sequence: ____________ | Strand: _
2. Position chr22:_______ | Sequence: ____________ | Strand: _

Discussion (3-4 sentences):
Why do you think some acceptor sites don't match the consensus?
What percentage did you expect? Is your result reasonable?
```

---

## 🧬 ASSIGNMENT 2: Splice Donor Sites (30 points)

### Background

**What are splice donor sites?**
- Mark where **exons end** and **introns begin**
- The first step in splicing: spliceosome recognizes donor site first
- More conserved than acceptor sites

**Consensus Pattern: MAG/GTRAGT**
```
M = A or C (aMino)
A = Adenine
G = Guanine
/ = conceptual splice point (not in sequence)
GT = Nearly invariant (GT-AG rule)
R = Purine (A or G)
AGT = Variable but preferred

The actual sequence is: MAGGTRAGT (remove the /)

Example sequences that match:
AAGGTAGGT
CAGGTAGAT
AAGGTAAGT
```

**Biology:** The GT-AG rule states that most introns start with GT (donor) and end with AG (acceptor).

---

### Your Tasks

#### Task 2.1: Build the Regular Expression (5 pts)

Convert `MAGGTRAGT` into a Python regular expression.

**Hints:**
- M (amino) = `[AC]`
- R (purine) = `[AG]`

**What to submit:**
```python
donor_pattern = r"[__]AGG[__][__]AGT"
```

---

#### Task 2.2: Extract Sequences Around Donor Sites (8 pts)

Similar to Assignment 1, but for donor sites.

**Step-by-step:**
1. Filter `splice_sites_df` for donor sites
2. Extract sequence from `(donor_position - 10)` to `(donor_position + 10)`
3. Handle strand correctly
4. Store sequences

**What to submit:**
- Code outline showing extraction
- First 5 sequences as examples

---

#### Task 2.3: Search for Pattern Matches (7 pts)

Search each donor sequence for your regex pattern.

**What to submit:**
- Code showing search logic
- Count of matches

---

#### Task 2.4: Compare with Acceptor Sites (10 pts)

**Calculate:**
- Total donor sites analyzed
- Number matching consensus
- Percentage match rate
- **COMPARISON**: Create a table comparing donor vs acceptor match rates

**What to submit:**

```
SPLICE DONOR SITE ANALYSIS
==========================

Dataset: Chromosome 22, GENCODE v49

Results:
- Total donor sites: ___
- Sites matching MAGGTRAGT: ___
- Match rate: ___%

COMPARISON TABLE:
┌─────────────────┬───────────┬─────────┬────────────┐
│ Site Type       │ Total     │ Matches │ Match Rate │
├─────────────────┼───────────┼─────────┼────────────┤
│ Acceptor (Y₁₀N) │ _____     │ _____   │ ____%      │
│ Donor (MAGGTR)  │ _____     │ _____   │ ____%      │
└─────────────────┴───────────┴─────────┴────────────┘

Discussion (1 paragraph, 5-6 sentences):
1. Which site (donor or acceptor) has a higher match rate?
2. Why might one be more conserved than the other?
3. What does this tell us about splicing biology?
4. Are your results consistent with the GT-AG rule?
```

---

## 🧬 ASSIGNMENT 3: p53 Binding Sites (45 points)

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

#### Task 3.1: Build the Complex Regular Expression (10 pts)

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

#### Task 3.2: Extract Promoter Sequences (8 pts)

You already have `promoter_df` from class with promoter coordinates (-1kb to TSS).

**Step-by-step:**
1. For each promoter, extract the 1000 bp sequence
2. Handle strand correctly (reverse complement for minus strand)
3. Store: gene_name, sequence, chromosome, start, end, strand

**What to submit:**
- Code outline for extraction
- First 3 promoter sequences (just first 100 bp of each)

---

#### Task 3.3: Search for p53 Binding Sites (12 pts)

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

#### Task 3.4: Analyze Spacer Length Distribution (8 pts)

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

#### Task 3.5: Comprehensive Report (7 pts)

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

#### Task 3.6: Biological Interpretation (BONUS +10 pts)

**Research Challenge:**

Pick **5 genes** from your top results and look them up on [GeneCards](https://www.genecards.org).

**Answer these questions:**
1. What is each gene's function?
2. Is it related to cancer, DNA repair, cell cycle, or apoptosis?
3. Does it make biological sense that p53 would regulate this gene?

**What to submit:**

```
BIOLOGICAL INTERPRETATION
========================

Gene 1: ________
- Function: ____________
- Cancer-related? Yes/No - Explain
- Expected p53 target? Yes/No - Why?

Gene 2: ________
... (continue for all 5)

CONCLUSION (1 paragraph):
Overall, do your computational results make biological sense?
Are most of your top hits actually involved in processes p53 regulates?
```

---

## 📤 Submission Instructions

### For Each Assignment, Submit:

1. **Jupyter Notebook (.ipynb)** or **Python Script (.py)**
   - Well-commented code
   - Clear section headers
   - All outputs visible (if notebook)

2. **Report (PDF or Markdown)**
   - All required statistics
   - Examples
   - Discussion/interpretation

3. **Data Files:**
   - Assignment 1: `acceptor_sequences.csv`
   - Assignment 2: `donor_sequences.csv` + comparison table
   - Assignment 3: `p53_sites.csv` (all sites found)

### File Naming Convention:
```
LastName_FirstName_Assignment1.ipynb
LastName_FirstName_Assignment1_Report.pdf
LastName_FirstName_acceptor_sequences.csv
```

### Submit via:
[Specify: Blackboard / Email / GitHub / etc.]

---

## 🎯 Grading Rubric

### Assignment 1: Acceptor Sites (25 pts)
- Regular expression pattern (5 pts)
- Sequence extraction code (8 pts)
- Pattern search implementation (7 pts)
- Report with statistics and examples (5 pts)

### Assignment 2: Donor Sites (30 pts)
- Regular expression pattern (5 pts)
- Sequence extraction code (8 pts)
- Pattern search implementation (7 pts)
- Comparison table and discussion (10 pts)

### Assignment 3: p53 Sites (45 pts)
- Complex regex pattern with spacer (10 pts)
- Promoter sequence extraction (8 pts)
- Search implementation (12 pts)
- Spacer analysis and visualization (8 pts)
- Comprehensive report (7 pts)
- **BONUS:** Biological interpretation (+10 pts)

---

## 💡 Tips for Success

1. **Start early** - Test your regex on small examples first
2. **Use the class notebook** - Reference your in-class work
3. **Test incrementally** - Don't write all code at once
4. **Ask questions** - Office hours, Piazza, email
5. **Verify your patterns** - Test on known sequences before full analysis
6. **Comment your code** - Explain your logic
7. **Check strand orientation** - Most common source of errors!

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

## ❓ FAQ

**Q: Can I work in groups?**
A: You can discuss approaches, but submit your own code and report.

**Q: My match rates seem low. Is that OK?**
A: Yes! Real biology is messy. Discuss why in your report.

**Q: Can I use a different chromosome?**
A: Stick with chr22 for consistency, unless you have a good reason.

**Q: What if I can't extract some sequences (negative coordinates)?**
A: Skip those. Document how many you skipped.

**Q: Can I use libraries other than pandas and re?**
A: Yes, but justify why. Keep it simple when possible.

---

## 🎓 Learning Goals Achieved

By completing these assignments, you will have:
- ✅ Applied computational pattern matching to real genomic data
- ✅ Connected DNA sequences to biological function
- ✅ Practiced the complete bioinformatics workflow
- ✅ Understood splice sites and transcription factor binding
- ✅ Analyzed real regulatory elements in the human genome
- ✅ Developed skills for future genomics research

**Good luck! 🧬🔬**

---

## 📧 Contact

**Instructor:** Pratik Dutta  
**Course:** BMI 503 - Introduction to Programming for Biomedical Informatics  
**Institution:** Stony Brook University  
**Office Hours:** [TBD]  
**Email:** [TBD]

---

*Last Updated: December 1, 2024*
