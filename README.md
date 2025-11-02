# Python for Omics Data Analysis in Bioinformatics

A comprehensive two-class course teaching students how to access, analyze, and visualize omics data using Python. This repository contains all materials including presentations, Jupyter notebooks, and datasets for hands-on learning.

## 🎯 Course Overview

This course introduces students to Python-based tools and workflows for analyzing three major types of omics data:
- **Genomics**: DNA sequences, variants, and structural genomics
- **Transcriptomics**: Gene expression and RNA-seq analysis
- **Imagomics**: Medical imaging and spatial biology


**Duration**: Two 3-hour classes  
**Prerequisites**: Basic Python programming, familiarity with biological concepts

## 📚 Course Structure

### Class 1: Python for Omics Data Access and Basic Analysis
**Learning Objectives:**
- Understand different types of omics data and their applications
- Learn to access public omics databases using Python
- Master basic data manipulation and analysis for each omics type
- Visualize omics data effectively

**Topics Covered:**
1. **Introduction to Omics** (30 min)
   - Overview of major omics types (genomics, transcriptomics, proteomics, metabolomics, epigenomics, metagenomics, imagomics)
   - Focus on genomics, transcriptomics, and imagomics
   - Data formats and biological applications

2. **Genomics with Python** (70 min)
   - Data types: FASTA, FASTQ, VCF, BED, BAM
   - Data sources: NCBI, Ensembl, UCSC Genome Browser, 1000 Genomes
   - Python libraries: Biopython, pysam, PyVCF, cyvcf2
   - Hands-on: Sequence manipulation, variant analysis

3. **Transcriptomics with Python** (70 min)
   - Data types: Count matrices, TPM, FPKM, h5ad (AnnData format)
   - Data sources: GEO, SRA, GTEx, TCGA
   - Python libraries: GEOparse, scanpy, pandas, PyDESeq2, anndata
   - Hands-on: Download RNA-seq data, quality control, differential expression

4. **Imagomics with Python** (30 min)
   - Data types: Whole slide images, H&E staining, multiplex imaging
   - Data sources: TCGA imaging, Human Protein Atlas, IDC
   - Python libraries: OpenSlide, scikit-image, Pillow, napari
   - Demo: Loading and visualizing pathology images

### Class 2: Advanced Python Workflows and Multi-Omics Integration
**Learning Objectives:**
- Perform advanced analysis on genomics and transcriptomics data
- Understand spatial transcriptomics and its applications
- Integrate multiple omics data types
- Build reproducible analysis pipelines

**Topics Covered:**
1. **Advanced Genomics Analysis** (40 min)
   - Variant annotation and functional prediction
   - Population genetics analysis
   - Structural variant detection
   - Hands-on: Annotate variants, calculate allele frequencies

2. **Advanced Transcriptomics Analysis** (40 min)
   - Single-cell RNA-seq analysis
   - Pathway and gene set enrichment analysis
   - Trajectory analysis and pseudotime
   - Hands-on: scRNA-seq clustering, marker gene identification

3. **Spatial Transcriptomics** (60 min)
   - Integration of imagomics and transcriptomics
   - Spatial gene expression analysis
   - Tools: scanpy, squidpy, spatialdata
   - Hands-on: Analyze 10x Visium spatial transcriptomics data

4. **Multi-Omics Integration & Project Ideas** (20 min)
   - Strategies for integrating multiple omics layers
   - Real-world applications and case studies
   - Final project suggestions

## 🗂️ Repository Structure

```
.
├── README.md
├── requirements.txt
├── environment.yml
│
├── Class1_Omics_Foundations/
│   ├── presentations/
│   │   ├── 01_Introduction_to_Omics.pdf
│   │   ├── 02_Genomics_Overview.pdf
│   │   ├── 03_Transcriptomics_Overview.pdf
│   │   └── 04_Imagomics_Overview.pdf
│   │
│   ├── notebooks/
│   │   ├── 01_genomics_sequence_analysis.ipynb
│   │   ├── 02_genomics_variant_analysis.ipynb
│   │   ├── 03_transcriptomics_data_access.ipynb
│   │   ├── 04_transcriptomics_differential_expression.ipynb
│   │   └── 05_imagomics_visualization_demo.ipynb
│   │
│   └── data/
│       └── README.md │
├── Class2_Advanced_Analysis/
│   ├── presentations/
│   │   ├── 01_Advanced_Genomics.pdf
│   │   ├── 02_Advanced_Transcriptomics.pdf
│   │   ├── 03_Spatial_Transcriptomics.pdf
│   │   └── 04_Multi_Omics_Integration.pdf
│   │
│   ├── notebooks/
│   │   ├── 01_variant_annotation.ipynb
│   │   ├── 02_population_genetics.ipynb
│   │   ├── 03_single_cell_rnaseq.ipynb
│   │   ├── 04_pathway_enrichment.ipynb
│   │   ├── 05_spatial_transcriptomics_analysis.ipynb
│   │   └── 06_multi_omics_integration.ipynb
│   │
│   └── data/
│       └── README.md (links to public datasets)

```

## 🔧 Installation & Setup

### Option 1: Using Conda (Recommended)
```bash
# Clone the repository
git clone https://github.com/yourusername/python-omics-analysis.git
cd python-omics-analysis

# Create conda environment
conda env create -f environment.yml
conda activate omics-analysis

# Launch Jupyter
jupyter lab
```

### Option 2: Using pip
```bash
# Clone the repository
git clone https://github.com/yourusername/python-omics-analysis.git
cd python-omics-analysis

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Launch Jupyter
jupyter lab
```

### Option 3: Google Colab (No Installation Required)
All notebooks are designed to run on Google Colab. Simply click the "Open in Colab" badge at the top of each notebook.

## 📦 Required Python Packages

### Core Libraries
- `numpy >= 1.24.0`
- `pandas >= 2.0.0`
- `matplotlib >= 3.7.0`
- `seaborn >= 0.12.0`
- `jupyter >= 1.0.0`

### Genomics
- `biopython >= 1.81`
- `pysam >= 0.21.0`
- `cyvcf2 >= 0.30.0`
- `pybedtools >= 0.9.0`

### Transcriptomics
- `scanpy >= 1.9.0`
- `anndata >= 0.9.0`
- `GEOparse >= 2.0.3`
- `PyDESeq2 >= 0.4.0`
- `gseapy >= 1.0.0`

### Imagomics
- `scikit-image >= 0.21.0`
- `Pillow >= 10.0.0`
- `openslide-python >= 1.3.0`
- `napari >= 0.4.18`

### Spatial Transcriptomics
- `squidpy >= 1.3.0`
- `spatialdata >= 0.0.14`

### Utilities
- `requests >= 2.31.0`
- `tqdm >= 4.65.0`


## 📊 Datasets Used

All datasets are publicly available. Links and download instructions are provided in the `data/README.md` files within each class folder.

### Genomics Datasets
- Human reference genome (GRCh38)
- 1000 Genomes Project variants
- Example VCF files from ClinVar

### Transcriptomics Datasets
- GEO: GSE164073 (bulk RNA-seq example)
- 10x Genomics PBMC dataset (scRNA-seq)
- GTEx tissue expression data

### Imagomics Datasets
- TCGA whole slide images
- Human Protein Atlas images
- Example H&E stained sections

### Spatial Transcriptomics Datasets
- 10x Visium spatial transcriptomics data
- Human DLPFC spatial dataset

## 🤝 Contributing

Contributions are welcome! If you find bugs, have suggestions, or want to add new tutorials:

1. Fork the repository
2. Create a new branch (`git checkout -b feature/improvement`)
3. Make your changes
4. Commit your changes (`git commit -am 'Add new feature'`)
5. Push to the branch (`git push origin feature/improvement`)
6. Create a Pull Request

## 📖 Additional Resources

### Online Tutorials
- [Biopython Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
- [Scanpy Tutorials](https://scanpy-tutorials.readthedocs.io/)
- [scikit-image Documentation](https://scikit-image.org/docs/stable/)

### Books
- *Python for Biologists* by Martin Jones
- *Bioinformatics Data Skills* by Vince Buffalo
- *Python for Data Analysis* by Wes McKinney

### Databases & Repositories
- [NCBI](https://www.ncbi.nlm.nih.gov/)
- [Ensembl](https://www.ensembl.org/)
- [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/)
- [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/tcga)
- [Human Protein Atlas](https://www.proteinatlas.org/)

### Community
- [Biostars Forum](https://www.biostars.org/)
- [SEQanswers Forum](http://seqanswers.com/)
- [r/bioinformatics Reddit](https://www.reddit.com/r/bioinformatics/)



## 📧 Contact

For questions about the course materials:
- Create an issue in this repository
- Email: [pratik.dutta@stonybrook.edu]

## 📝 License

This work is licensed under a Creative Commons Attribution 4.0 International License (CC BY 4.0).

You are free to:
- Share — copy and redistribute the material
- Adapt — remix, transform, and build upon the material



## 🙏 Acknowledgments

- Stony Brook University Department of Biomedical Informatics
- BMI 503 students for their feedback and contributions
- Open-source bioinformatics community
- Public data providers (NCBI, TCGA, GEO, etc.)

---

**Course**: BMI 503 - Introduction to Computer Science for Biomedical Informatics  
**Semester**: Fall 2025  
**Last Updated**: November 2025
