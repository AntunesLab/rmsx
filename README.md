# RMSX Trajectory Analysis Tool Tutorial

Welcome to the tutorial for setting up and running the **RMSX Trajectory Analysis Tool**. This guide will walk you through the installation process, setup, and usage of the tool for analyzing molecular dynamics (MD) trajectories.

---

## Table of Contents
- [Quick Start with Google Colab](#quick-start)
- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
    - [Software Requirements](#software-requirements)
    - [Python Dependencies](#python-dependencies)
    - [R Dependencies](#r-dependencies)
- [Installation](#installation)
    - [1. Clone the Repository](#1-clone-the-repository)
    - [2. Install the Python Package](#2-install-the-python-package)
    - [3. Install R and Required Packages](#3-install-r-and-required-packages)
- [Usage](#usage)
    - [Command-Line Interface](#command-line-interface)
    - [Python Module](#python-module)
- [Examples](#examples)
    - [Example 1: Basic Usage](#example-1-basic-usage)
    - [Example 2: Advanced Options](#example-2-advanced-options)
    - [Example 3: Selecting a Specific Chain](#example-3-selecting-a-specific-chain)
- [Troubleshooting](#troubleshooting)
- [Additional Resources](#additional-resources)
- [Contact Information](#contact-information)
- [Contributing](#contributing)
- [License](#license)

---
## Quick Start with Collab

```python
!pip install rpy2
%load_ext rpy2.ipython
```

```r
%%R
install.packages(c('ggplot2', 'gridExtra', 'reshape2', 'scales'), repos='http://cran.r-project.org/')

```


```python
# Clone the RMSX repository (replace with the actual repository URL)
!git clone https://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@github.com/AntunesLab/rmsx.git
!pip install -e ./rmsx/.
```

    

```python
import os, sys

# Add the parent directory of 'rmsx' to sys.path
sys.path.append('/content/rmsx/rmsx')
print(sys.path)
```

```python
from rmsx import run_rmsx
```

```python
psf_file = "/content/rmsx/test_files/1UBQ.psf"
dcd_file = "/content/rmsx/test_files/mon_sys.dcd"
pdb_file = "/content/rmsx/test_files/1UBQ.pdb"
output_dir = "/content/rmsx/test_files/colab_demo_data2"


run_rmsx(psf_file, dcd_file, pdb_file, output_dir, 15, 'Rscript', verbose=True, interpolate=False, triple=False)
```

## Introduction

The **RMSX Trajectory Analysis Tool** is a Python package designed to facilitate the analysis of molecular dynamics simulations. It automates the processing of trajectory data, calculation of RMSD (Root Mean Square Deviation) and RMSF (Root Mean Square Fluctuation), and generates informative plots using R scripts. The tool supports various file formats compatible with MDAnalysis and provides both a command-line interface (CLI) and a Python API for flexible use.

---

## Prerequisites

### Software Requirements

- **Operating System**: Linux, macOS, or Windows
- **Python**: Version 3.6 or higher
- **R**: Version 3.6 or higher

### Python Dependencies

The following Python packages are required:

- **MDAnalysis**
- **numpy**
- **pandas**
- **matplotlib** (optional, for plotting within Python)
- **pkg_resources**
- **argparse**

These dependencies will be installed automatically when you install the package using `pip`.

### R Dependencies

The R script used for plotting requires the following R packages:

- **ggplot2**
- **gridExtra**
- **reshape2**
- **scales**

---

## Installation

### 1. Clone the Repository

First, clone the repository from GitHub to your local machine:
```
git clone https://github.com/finn2400/rmsx.git cd rmsx
```
### 2. Install the Python Package

Install the package using `pip`. It's recommended to use a virtual environment.

#### Option A: Using a Virtual Environment

Create and activate a virtual environment:
```
# For virtualenv
python -m venv venv source venv/bin/activate
# On Windows: venv\Scripts\activate
# For conda conda create -n rmsx_env python=3.8 conda activate rmsx_env
```
Install the package in editable mode (useful for development):

```
pip install -e .
```
#### Option B: Install Globally

```
pip install .
```


### 3. Install R and Required Packages

#### Install R

If you don't have R installed, download and install it from CRAN.

#### Install Required R Packages

Open R or RStudio and install the required packages:

```
install.packages(c("ggplot2", "gridExtra", "reshape2", "scales"))
```
Alternatively, you can install the packages from the command line:

```
Rscript -e "install.packages(c('ggplot2', 'gridExtra', 'reshape2', 'scales'), repos='http://cran.rstudio.com/')"
```
---

## Usage

### Command-Line Interface

The package provides a command-line tool `rmsx_cli.py` for easy usage.

#### Syntax
```
rmsx_cli.py [options] psf_file dcd_file pdb_file
```
#### Positional Arguments

- `psf_file`: Path to the topology file (e.g., `.psf`, `.gro`)
- `dcd_file`: Path to the trajectory file (e.g., `.dcd`, `.trr`)
- `pdb_file`: Path to the PDB file containing the structure

#### Optional Arguments

- `--output_dir OUTPUT_DIR`: Output directory to save results (default: current directory)
- `--slice_size SLICE_SIZE`: Number of frames per slice for trajectory processing (default: 5)
- `--verbose`: Enable verbose output
- `--interpolate`: Enable interpolation in plots
- `--triple`: Generate triple plots
- `--help`: Show help message and exit

#### Accessing Help

To view the help message and see all options:

```
rmsx_cli.py --help
```
### Python Module

You can also use the tool directly within your Python scripts or Jupyter notebooks.

#### Importing the Module
```
from rmsx.analysis import run_rmsx
```
#### Function Syntax

```
run_rmsx(
    psf_file,
    dcd_file,
    pdb_file,
    r_script_name='triple_plot_rmsx.R',
    output_dir=None,
    slice_size=5,
    rscript_executable='Rscript',
    verbose=True,
    interpolate=False,
    triple=False)
```
#### Parameters

- `psf_file`: Path to the topology file
- `dcd_file`: Path to the trajectory file
- `pdb_file`: Path to the PDB file
- `r_script_name`: Name of the R script for plotting (default: `'triple_plot_rmsx.R'`)
- `output_dir`: Directory to save outputs (default: current directory)
- `slice_size`: Number of frames per slice (default: 5)
- `rscript_executable`: Command to execute R scripts (default: `'Rscript'`)
- `verbose`: Enable verbose output (default: `True`)
- `interpolate`: Enable interpolation in plots (default: `False`)
- `triple`: Generate triple plots (default: `False`)

---

## Examples

### Example 1: Basic Usage

Process a trajectory with default settings and generate plots.

#### Command Line

```
rmsx_cli.py structure.psf trajectory.dcd structure.pdb
```
#### Python Script

```
from rmsx.analysis import run_rmsx  run_rmsx(psf_file='structure.psf', dcd_file='trajectory.dcd', pdb_file='structure.pdb' )
```
### Example 2: Advanced Options

Process a trajectory, specify output directory and slice size, enable verbose output and triple plotting.

#### Command Line
```
rmsx_cli.py structure.psf trajectory.dcd structure.pdb \ --output_dir analysis_results \ --slice_size 20 \ --verbose \ --triple 
```
#### Python Script

```
from rmsx.analysis import run_rmsx
run_rmsx(psf_file='structure.psf', dcd_file='trajectory.dcd', pdb_file='structure.pdb', output_dir='analysis_results', slice_size=20, verbose=True,     triple=True)
```

### Example 3: Selecting a Specific Chain

When you run the tool, it will prompt you to select a chain from your structure for analysis.


`Available chains and their lengths (in residues): Chain A: 153 residues Chain B: 160 residues Please enter the chain ID you would like to analyze from the following options: A (153 residues), B (160 residues) Chain ID:`

Enter the chain ID (e.g., `A`) to proceed with the analysis on that specific chain.

---

## Troubleshooting

### Common Issues and Solutions

#### 1. R Script Not Found

**Error Message:**


`Fatal error: cannot open file '.../rmsx/r_scripts/triple_plot_rmsx.R': No such file or directory`

**Solution:**

- Ensure that the R scripts are included in your installation.
- Reinstall the package using `pip install -e .` in the package directory.
- Verify that the `r_scripts` directory exists within the installed package.

#### 2. MDAnalysis Errors Loading Files

**Error Message:**

```
ValueError: Universe.load_new(): File ... cannot be read
```
**Solution:**

- Verify that your topology and trajectory files are in compatible formats supported by MDAnalysis.
- Check that the files are not corrupted and correspond to the same simulation.

#### 3. R Packages Not Installed

**Error Message:**

```
Error in library(ggplot2) : there is no package called ‘ggplot2’
```
**Solution:**

- Install the required R packages using the instructions in the [R Dependencies](#r-dependencies) section.

#### 4. Chain Not Found

**Error Message:**

```
Chain 'C' is not available in the topology file.
```
**Solution:**

- Ensure that you enter a valid chain ID as displayed in the prompt.
- Chain IDs are case-sensitive; make sure to match the exact ID.

### Getting Help

If you encounter issues not covered here:

- Check the documentation and examples.
- Open an issue on the GitHub repository.
- Contact the maintainer via email.

---

## Additional Resources

- **MDAnalysis Documentation**: MDAnalysis User Guide
- **R Documentation**: [R Project](https://www.r-project.org/)

---

## Contact Information

If you have any questions, suggestions, or need assistance, please feel free to contact the maintainer:

- **Email**: fpberuld@cougarnet.uh.edu
- **GitHub Issues**: [GitHub Repository Issues](https://github.com/finn2400/rmsx/issues)

---

## License

This project is licensed under the MIT License. 
---

Thank you for using the RMSX Trajectory Analysis Tool! We hope this tutorial has helped you set up and run the tool effectively. Happy analyzing!
