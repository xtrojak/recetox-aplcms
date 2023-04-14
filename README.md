# recetox-aplcms
![GitHub R package version](https://img.shields.io/github/r-package/v/RECETOX/recetox-aplcms)
[![Conda](https://img.shields.io/conda/v/bioconda/r-recetox-aplcms)](https://anaconda.org/bioconda/r-recetox-aplcms)
[![codecov](https://codecov.io/gh/RECETOX/recetox-aplcms/branch/master/graph/badge.svg?token=1B664Z8DHT)](https://codecov.io/gh/RECETOX/recetox-aplcms)

apLCMS is a software which generates a feature table from a batch of LC/MS spectra. The m/z and retention time tolerance levels are estimated from the data. A run-filter is used to detect peaks and remove noise. Non-parametric statistical methods are used to find-tune peak selection and grouping. After retention time correction, a feature table is generated by aligning peaks across spectra. For further information on apLCMS please refer to https://mypage.cuhk.edu.cn/academics/yutianwei/apLCMS/.

This is a fork of the official [aplcms repo](https://github.com/tianwei-yu/apLCMS) that takes the project towards large-scale MS analyses.

## Installation

The newest version of the package can be installed through conda from [bioconda](https://anaconda.org/bioconda/r-recetox-aplcms) channel:

```
conda install -c bioconda r-recetox-aplcms
```

Alternatively, a series of Galaxy tools is available at https://github.com/RECETOX/galaxytools/tree/master/tools/recetox_aplcms.

## Usage

The tool can generate a feature table from a batch of LC/MS spectra given as a collection of mzML files. mzML is a well tested open-source format for mass spectrometer output files that can be readily utilized by the community and easily adapted for incremental advances in mass spectrometry technology.

In contrast to well-known XCMS tool, apLCMS can process profile mode data and fits a bi-Gaussian peak shape model to the data, resulting in better peak detection than XCMS. Drifts in retention time are also corrected in the tool, outputting an aligned feature table.

It operates in two modes - `unsupervised` and `hybrid`. `Unsupervised` mode of apLCMS is not relying on any existing knowledge about metabolites or any historically detected features. On the other hand, `Hybrid` version of apLCMS is incorporating the knowledge of known metabolites and historically detected features on the same machinery to help detect and quantify lower-intensity peaks. To use such knowledge, especially historical data, you must keep using the same chromatography system (otherwise the retention time will not match), and the same type of samples with similar extraction technique, such as human serum. For both modes, an equally-named function is exposed, parametrised with multiple arguments.

## Testing
Before being able to run the tests, it is necessary to fetch the required data using the following commands:

```
wget -P tests/testdata/adjusted -i tests/remote-files/adjusted.txt
wget -P tests/testdata/aligned -i tests/remote-files/aligned.txt
wget -P tests/testdata/extracted -i tests/remote-files/extracted.txt
wget -P tests/testdata/input -i tests/remote-files/input.txt
wget -P tests/testdata/recovered -i tests/remote-files/recovered.txt
wget -P tests/testdata/recovered/recovered-extracted -i tests/remote-files/recovered-extracted.txt
wget -P tests/testdata/recovered/recovered-corrected -i tests/remote-files/recovered-corrected.txt
wget -P tests/testdata/filtered -i tests/remote-files/filtered.txt
wget -P tests/testdata/filtered/run_filter -i tests/remote-files/run_filter.txt
wget -P tests/testdata/features -i tests/remote-files/features.txt
wget -P tests/testdata/clusters -i tests/remote-files/clusters.txt
wget -P tests/testdata/hybrid -i tests/remote-files/hybrid.txt
wget -P tests/testdata/template -i tests/remote-files/template.txt
wget -P tests/testdata/unsupervised -i tests/remote-files/unsupervised.txt
```

The `hybrid` and `unsupervised` tests are [reported](https://github.com/RECETOX/recetox-aplcms/issues/24) to be OS specific and may fail depending on the platform they are run on. To ensure reproducibility during development process you can run the tests in a designated Docker container as follows:
```
# from the repository root run
$ docker build -t recetox-aplcms .
```
After `docker-build` has built the image run:
```
$ docker run --rm -t -v $(pwd):/usr/src/recetox-aplcms recetox-aplcms
```
This will create a container and automatically run all the tests from the **tests** folder.

# Documentation for developers

## Setting up a development environment
The development environment can be set up in two ways, either via **VSCode's devcontainer** extension or a **docker container**.

### Devcontainer
To use a devcontainer you need VSCode with [Remote - Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension and docker installed on your machine:
- Clone your fork of the repository and open the folder in VSCode;
- From VSCode's command palette run `Remote-Containers: Open Folder in Container`. VSCode may take a few minutes building a container;
- After container is ready, open a **new** terminal and type `conda activate recetox-aplcms-dev` to activate Conda environment;
- Run `R` or `radian` to enter R terminal (we recommend `radian` due to its ease of use);
- A good starting point would be fetching the test data as described above, running `devtools::test()` and waiting until all tests pass to ensure the environment is set correctly.

### Docker container
To use a docker development environment you need **Docker** installed on your machine. If you don't have **Docker** you can follow installation instructions on Docker's [web](https://docs.docker.com/engine/install/).
- Clone your fork of the repository;
- From the package root folder run `docker build -t recetox-aplcms .` to build an image. This may take a few minutes.
- After the image is build start the container:
    ```bash
    $ docker run -it \
    -v $(pwd):/usr/src/recetox-aplcms \
    --entrypoint '/bin/bash' \
    recetox-aplcms
    ```
- Once in container, finish setting up the environment by running:
    ```bash
    $ apt update && apt upgrade
    ```
    ```shell
    $ apt install git && git config --global --add safe.directory /usr/src/recetox-aplcms
    ```
- Enter a Conda environment by running `conda activate recetox-aplcms-dev`
- Run `R` or `radian` to enter R terminal (we recommend `radian` due to its ease of use);
- A good starting point would be fetching the test data as described above, running `devtools::test()` and waiting until all tests pass to ensure the environment is set correctly.


## References
Yu, T., Park, Y., Johnson, J. M. & Jones, D. P. apLCMS—adaptive processing of high-resolution LC/MS data. Bioinformatics 25, 1930–1936 (2009). DOI: [10.1093/bioinformatics/btp291](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btp291).

Yu, T., Park, Y., Li, S. & Jones, D. P. Hybrid Feature Detection and Information Accumulation Using High-Resolution LC–MS Metabolomics Data. J. Proteome Res. 12, 1419–1427 (2013). DOI: [10.1021/pr301053d](https://pubs.acs.org/doi/10.1021/pr301053d).

Yu, T. & Jones, D. P. Improving peak detection in high-resolution LC/MS metabolomics data using preexisting knowledge and machine learning approach. Bioinformatics 30, 2941–2948 (2014). DOI: [10.1093/bioinformatics/btu430](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu430).

Yu, T. & Peng, H. Quantification and deconvolution of asymmetric LC-MS peaks using the bi-Gaussian mixture model and statistical model selection. BMC Bioinformatics 11, 559 (2010). DOI: [10.1186/1471-2105-11-559](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-559).

Liu, Q. et al. Addressing the batch effect issue for LC/MS metabolomics data in data preprocessing. Sci. Rep. 10, 13856 (2020). DOI: [10.1038/s41598-020-70850-0](https://doi.org/10.1038/s41598-020-70850-0).
