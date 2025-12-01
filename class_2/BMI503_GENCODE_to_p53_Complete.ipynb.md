# BMI 503 - Genomic Annotations and Motif Finding
## From GENCODE to p53 Binding Sites

**Date:** December 2, 2024  
**Instructor:** Pratik Dutta  
**Duration:** 90 minutes

---

## 📚 Table of Contents


10. [Real-World Example: Complete Pipeline](#section10)


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
