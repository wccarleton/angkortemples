# Project
## Overview
This repo contains the data and code used for the study presented in the following paper:

[*Bayesian regression versus machine learning for rapid dating of archaeological features identified in lidar scans using Angkor, Cambodia, as a case study*]()

### Note
The paper describing the anlayses conducted with the code in this repo has been submitted for review and this page will be upadated with more details once the paper has been accepted for publication. This repo may, therefore, be changed in response to reviewer comments. This version of the repo has, nevertheless, been archived with Zenodo.

## Abstract
Lidar (light-detection and ranging) has revolutionized archaeology. We are now able to produce high-resolution maps of archaeological surface features over vast areas, allowing us to see ancient land-use and anthropogenic landscape modification at previously un-imagined scales. In the tropics, this has enabled documentation of previously archaeologically unrecorded cities in various tropical regions, igniting scientific and popular interest ancient tropical urbanism. An emerging challenge, however, is to add temporal depth to this torrent of new spatial data because traditional archaeological investigations are time consuming and inherently destructive. So far, we are aware of only one attempt to apply statistics and machine learning to lidar data in order to add time-depth to newly acquired spatial data. Using temples at the well-known massive urban complex of Angkor in Cambodia as a case study, a predictive model was developed combining standard regression with novel machine learning methods to estimate temple foundation dates for undated Angkorian temples identified by lidar scans. The model's predictions were used to produce an historical population curve for Angkor and study urban expansion at this important ancient tropical urban centre. The approach, however, has certain limitations. Importantly, its handling of uncertainties leaves room for improvement, and like many machine learning approaches it is opaque regarding which predictor variables are most relevant. Here we describe a new study in which we investigated an alternative Bayesian regression approach applied to the same case study. We compare the two models in terms of their inner workings, results, and interpretive utility. We also use an updated database of Angkorian temples as the training dataset, allowing us to produce the most current estimate for temple foundations and historic spatiotemporal urban growth patterns at Angkor. Our results demonstrate that, in principle, predictive statistical and machine learning methods can be used to rapidly add chronological information to large lidar datasets and a Bayesian paradigm makes it possible to incorporate important uncertainties---especially chronological---into modelled temporal estimates.

## Software
The R scripts contained in this repository are intended for replication efforts and to improve the transparency of research. They are, of course, provided without warranty or technical support. That said, questions about the code can be directed to me, Chris Carleton, at ccarleton@protonmail.com.

### R
This analysis described in the associated manuscript was performed in R. Thus, you may need to download the latest version of [R](https://www.r-project.org/) in order to make use of the scripts described below.

### Nimble
This project made use of a Bayesian Analysis package called [Nimble](https://r-nimble.org/). See the Nimble website for documentation and a tutorial. Then, refer to the R scripts in this repo.

## Contact

[ORCID](https://orcid.org/0000-0001-7463-8638) |
[Google Scholar](https://scholar.google.com/citations?hl=en&user=0ZG-6CsAAAAJ) |
[Website](https://wccarleton.me)

## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
