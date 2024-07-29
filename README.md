# Longitudinal Principal Manifold Estimation (LPME)

This repository includes the documents needed to reproduce the simulations and applied analysis demonstrating potential uses of the Longitudinal Principal Manifold Estimation (LPME) algorithm developed in the following paper by Robert Zielinski, Kun Meng, and Ani Eloyan.

[Longitudinal Principal Manifold Estimation](https://arxiv.org/pdf/2407.17450)

This project is intended to create a longitudinal extension of the Principal Manifold Estimation (PME) algorithm introduced in [Meng and Eloyan, 2021](https://pubmed.ncbi.nlm.nih.gov/35813449/). The proposed LPME algorithm allows for the smooth estimation of manifolds over time. The files included here contain notebooks allowing the complete replication of the comparisons between LPME and naive applications of the PME and principal curve algorithms on simulated data. Notebooks also describe the analysis of hippocampus and thalamus data from structural MRI data from the Alzheimer's Disease Neuroimaging Initiative (ADNI) dataset. The relevant files needed to reproduce results can be found in the `code/` directory at the following locations:

- `code/01_lpme_simulations.qmd`
- `code/02_lpme_simulation_analysis.qmd`
- `code/03_adni_mri_preprocess.qmd`
- `code/04_lpme_adni_application.qmd`
- `code/05_lpme_adni_results.qmd`
- `code/06_generate_lpme_figures.qmd`
- `code/07_generate_misc_figures.qmd`

The code in these files often depend on helper functions, which can be found in the `code/functions` directory.

## Dependency installation:

We used the [renv package](https://rstudio.github.io/renv/) to manage this project's dependency environment. To install most of the necessary dependencies, use the following steps:

1. Clone this repository using `git clone [LINK]`.
2. Install the `renv` package using `install.packages()`.
3. Call `renv::restore()` to install dependencies.

Please see the [renv documentation](https://rstudio.github.io/renv/) for more information about this process.

The [fslr](https://github.com/muschellij2/fslr) package, used to preprocess MRI data, relies on a separate installation of [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/). Instructions for installation of the FSL software can be found [here](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallationhttps://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation). For more information about configuring your FSL installation to work with the `fslr` package, see the "fslr setup" section in [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4911193/).

Also, note that for this analysis, the PME and LPME algorithms were implemented in the [pme package](https://github.com/rjzielinski/pme). This package is still in the early stages of development, and may undergo substantial syntax changes in the future. Efforts will be made to maintain compatibility between the code in this repository and the updated versions of the `pme` package, but if you do experience errors related to this package, please create an issue in this repository so the problem can be addressed.

## ADNI Data

To demonstrate the use of the LPME algorithm for a real dataset, we considered the surfaces of the thalamuses and hippocampuses of participants of the ADNI study, which offers access to high quality longitudinal imaging data for individuals with Alzheimer's disease, as well as a number of healthy control participants. Use of the ADNI data is governed by a Data Use Agreement, so we cannot make this data publicly available. However, those interested can easily apply for access to the ADNI data as described [here](https://adni.loni.usc.edu/data-samples/access-data/). Upon receiving access to the data, the process used to download the data for this analysis is described below.

- Log in to the Data Archive, accessed [here](https://ida.loni.usc.edu/login.jsp?project=ADNI).
- Under the "Download" tab at the top of the page, select "Image Collections".
- Select "Advanced Search", and check the boxes for the ADNI 1 project phase, MRI modality, and T1 weighting. Then search for the relevant images.
- Select all of the images returned in the search results and add to a collection.
- Navigate to your collection under the "Data Collections" tab and select all the images. Then hit the "1-Click Download" button, and click the "Zip File 1" link on the subsequent window to initiate the download process. Note that this will be a large download and will likely take a long time to complete.
- Move your downloaded dataset to the `data/` directory and unzip the file to `data/adni`. Now you are ready to proceed to the image preprocessing steps in `code/03_adni_mri_preprocess.qmd`.

