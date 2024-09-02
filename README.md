# Repository for the manuscript "A deep mutational scanning platform to characterize the fitness landscape of anti-CRISPR proteins"

This is the companion repository to the manuscript "A deep mutational scanning platform to characterize the fitness landscape of anti-CRISPR proteins".
It contains a Jupyter Notebook to reproduce the analyses and figures in the manuscript.

### Installation

*Note: If you are working on Windows, please use WSL (Windows Subsystem for Linux) for the following steps.*

If not already installed, install `conda` (e.g. using [miniconda](https://docs.conda.io/en/latest/miniconda.html)).

Use `git` to clone this repository:

```
git clone https://github.com/mjendrusch/acr-dms.git
```

Navigate into the newly created directory:

```bash
cd acr-dms
```

Source the install-script `setup.sh`:

```bash
source setup.sh
```

This will create and activate a new conda-environment `acrdms` containing all dependencies needed to reproduce our analyses.


Finally, download `acrdms_data_v2.tar.gz` from [our dataset deposited at Zenodo](https://zenodo.org/records/13374667) and unpack it in this directory.

```bash
tar -xzf acrdms_data_v2.tar.gz
```

This will create the data directory containing all raw data needed for analysis.

### Reproduce

To reproduce the anaĺyses in our manuscript, follow the installation instructions and start jupyter:

```
jupyter notebook
```

Open the notebook `AcrDMS.ipynb` and run the cells you want to evaluate from top to bottom, or click on "Cell" in the toolbar and select "Run All".
This will run all analyses and reproduce all plots used in our manuscript.

### Run InDel count
`ngs_indelcount.py` processes a set of multiplexed amplicon-seq results to extract the fraction of reads containing indels.
It demultiplexes the data using provided barcodes in `barcodes.csv` and compares the length of reads to the ground-truth sequence
in `sequence.fa` to detect reads containing insertions or deletions.
All `.csv` files for `ngs_indelcount.py` need to be ","-separated.

Steps:
1. Set up a directory containing all paired-end read fastq files for all replicates. In that directory:
2. Create a file named `readfiles.csv` in this directory, with one row per replicate and columns: `replicate name,relative path to read file 1,relative path to read file 2`
3. Create a file named `barcodes.csv` with format `condition name,forward barcode,reverse barcode`, with one header row and one row per barcoded condition.
   The barcodes need to be provided in the orientation they would appear in for a single merged paired-end read.
4. Create a file named `sequence.fa` with the expected sequence of the region being sequenced
5. Run `python ngs_indelcount.py path/to/directory/`.
This will create a file named `outputs.csv` in that directory, which counts the fraction of paired-end reads with indels for each replicate with format:
`condition name,combined barcode,fraction replicate 1,fraction replicate 2,...`
