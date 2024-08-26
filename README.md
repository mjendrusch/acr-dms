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
tar -xzf acr-dms-data.tar.gz
```

This will create the data directory containing all raw data needed for analysis.

### Reproduce

To reproduce the anaĺyses in our manuscript, follow the installation instructions and start jupyter:

```
jupyter notebook
```

Open the notebook `AcrDMS.ipynb` and run the cells you want to evaluate from top to bottom, or click on "Cell" in the toolbar and select "Run All".
This will run all analyses and reproduce all plots used in our manuscript.

