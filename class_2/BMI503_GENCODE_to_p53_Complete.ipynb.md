# BMI 503 - Genomic Annotations and Motif Finding
## From GENCODE to p53 Binding Sites

**Date:** December 2, 2024  
**Instructor:** Pratik Dutta  
**Duration:** 90 minutes

---

## 📚 Table of Contents


 [Real-World Example: Complete Pipeline](#section10)


---

## Assignment 

```
Find p53 binding sites in human promoters using:
1. GTF parsing (get TSS locations)
2. Sequence extraction (promoter regions)
3. Regular expressions (find motifs)

p53 motif: RRRCWWGYYY {0-10bp} RRRCWWGYYY
```



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

```

### Search Promoter Sequences for p53 Sites







**END OF NOTEBOOK**
