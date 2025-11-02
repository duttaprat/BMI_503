{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Using `pybedtools` in Google Colab to Extract DNA Sequences\n",
    "\n",
    "This notebook demonstrates the complete, correct workflow for installing `pybedtools` and its dependencies in Google Colab and using it to extract specific sequences from a FASTA file based on genomic coordinates."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Step 1: Install `pybedtools` and `bedtools` using Conda\n",
    "\n",
    "The most reliable way to install bioinformatics tools in Colab is with `conda`, as it handles complex dependencies correctly. The default tools installed with `!apt-get` are often too old and cause errors.\n",
    "\n",
    "First, we install `condacolab`, which sets up a conda environment for us."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "!pip install -q condacolab\n",
    "import condacolab\n",
    "condacolab.install()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "After running the cell above, the Colab runtime will automatically restart. **This is normal and necessary.**\n",
    "\n",
    "Now we can use `conda` to install the correct, compatible versions of our tools from the trusted `bioconda` channel."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "!conda install -c bioconda bedtools pybedtools"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Step 2: Import Libraries and Download Data\n",
    "\n",
    "Now that our tools are installed, we can import `pybedtools`.\n",
    "\n",
    "We also need a reference genome to work with. We'll download the same sample FASTA file we used before, which contains sequences for `chr1` and `chr2`."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "import pybedtools\n",
    "import os\n",
    "\n",
    "# Download the sample FASTA file\n",
    "!wget -q -O sample_genome.fa https://raw.githubusercontent.com/duttaprat/BMI_Bootcamp_2023/master/Data_set/Genome_dataset/test_FASTA.fa\n",
    "\n",
    "# Define the path to our downloaded file\n",
    "fasta_path = \"sample_genome.fa\"\n",
    "\n",
    "print(f\"FASTA file downloaded to: {fasta_path}\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Step 3: Define Regions and Create a BedTool Object\n",
    "\n",
    "Next, we need to specify which genomic regions we want to extract. We can define these regions directly in a string using the standard BED format (`chrom  start  end`). We can then load this string directly into a `pybedtools.BedTool` object.\n",
    "\n",
    "Let's define two regions from `chr1` and one from `chr2`."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Define the coordinates for the regions we want to extract\n",
    "# Format: tab-separated string of \"chromosome start end name\"\n",
    "bed_definitions = \"\"\"\n",
    "chr1\\t10\\t40\\tRegionA_GeneX\n",
    "chr1\\t100\\t150\\tRegionA_PromoterY\n",
    "chr2\\t50\\t90\\tRegionB_GeneZ\n",
    "\"\"\"\n",
    "\n",
    "# Create a BedTool object from our string definition\n",
    "regions_to_extract = pybedtools.BedTool(bed_definitions, from_string=True)\n",
    "\n",
    "print(\"Created a BedTool object with the following regions:\")\n",
    "print(regions_to_extract)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Step 4: Extract the Sequences\n",
    "\n",
    "This is the key step. We use the `.sequence()` method on our `BedTool` object, telling it which FASTA file to use as a reference.\n",
    "\n",
    "The `name=True` argument is very useful: it ensures that the output FASTA headers match the names we provided in our BED definitions (e.g., `RegionA_GeneX`)."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Use the .sequence() method to extract the DNA for our defined regions\n",
    "extracted_sequences = regions_to_extract.sequence(\n",
    "    fi=fasta_path,      # Path to the reference FASTA file\n",
    "    name=True           # Use the name from the BED file as the FASTA header\n",
    ")\n",
    "\n",
    "# The output is a new BedTool object pointing to a temporary file with the results\n",
    "print(f\"Sequences extracted to temporary file: {extracted_sequences.fn}\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Step 5: View the Results\n",
    "\n",
    "Finally, let's print the contents of the output file to see our extracted DNA sequences. The output is in FASTA format, with the correct names and sequences for each region we requested."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# The BedTool object can be treated like a file to read its contents\n",
    "print(\"--- Extracted Sequences (FASTA format) ---\")\n",
    "print(open(extracted_sequences.fn).read())"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.10.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
