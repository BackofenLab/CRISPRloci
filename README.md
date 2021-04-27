# CRISPRloci
####CRISPRloci provides an automated and comprehensive in silico characterization of CRISPR-Cas system on bacterial and archaeal genomes. It is a full suite for CRISPR locus characteriztion that includes CRISPR array orientation, detection of conserved leaders, cas gene annotation and subtype classification.

The web server interface of CRISPRloci is freely available at: rna.informatik.uni-freiburg.de/trunk/CRISPRloci
 
 
![GitHub Logo](/images/Webserver_flowchart.png)

### Installation and requirements

CRISPRloci_standalone.py has been tested with Python 3.7 To run it, we recommend installing the same library versions we used. Since we exported our classifiers following the [model persistence guideline from scikit-learn](https://scikit-learn.org/stable/modules/model_persistence.html), it is not guaranteed that they will work properly if loaded using other Python and/or library versions. For such, we recommend the use of our docker image or a conda virtual environment. They make it easy to install the correct Python and library dependencies without affecting the whole operating system (see below).

### First step: download the last version of the tool and extract it


```
wget https://github.com/BackofenLab/CRISPRloci/archive/v1.0.0.tar.gz
tar -xzf v1.0.0.tar.gz
```

### Second step: download the Hidden Markov (HMM) and Machine Learning (ML) models

Due to GitHub's file size constraints, we made our HMM and ML models available in Google Drive. You can download them [here](https://drive.google.com/file/d/1YbTxkn9KuJP2D7U1-6kL1Yimu_4RqSl1/view?usp=sharing) and [here](https://drive.google.com/file/d/1Nc5o6QVB6QxMxpQjmLQcbwQwkRLk-thM/view?usp=sharing). 
Save both tar.gz files inside CRISPRcasIdentifier's directory. It is not necessary to extract them, since the tool will do that the first time it is run.


### Third step: download the Hidden Markov (HMM) and Machine Learning (ML) models

We made our HMM and ML models available in Google Drive. You can download them from the following links:

* [Machine Learning Models](https://drive.google.com/file/d/1SkyS03hQG0P7bO7KvaQXgQ_J2Vfng85K/view?usp=sharing)
* [General HMM Models](https://drive.google.com/file/d/1yZ3rl0LPIk-LDLKl2KHb3ra5O2h9_jpA/view?usp=sharing)
* [Signature HMM Models](https://drive.google.com/file/d/1itiqV8djmrfwgwsByZLGe9olv7Uri5hV/view?usp=sharing)
* [Cas HMM Models](https://drive.google.com/file/d/1zRJNlgqAC6A8BEDXPNKmtN4AfnVhjJIR/view?usp=sharing)

Save all tar.gz files inside Casboundary's folder. It is not necessary to extract them, since the tool will do that the first time it is run.


### Fourth step (conda)

First we install Miniconda for python 3.
Miniconda can be downloaded from here: [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Install Miniconda.

``
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
``

Create and activate environment for CRISPRloci.

```
conda env create -f CRISPRloci-env.yml -n CRISPRloci-env
conda activate CRISPRloci-env
```

After using CRISPRloci_standalone.py you can deactivate the environment.

```
conda deactivate
```

### Quick run with the default parameters

#### DNA mode
In order to test the dna mode please execute the following the command 

```
python3.7 CRISPRloci_standalone.py -f Example/NC_005230.fasta -st dna
```

#### Protein mode
In order to test the protein mode tun the following command:
```
python3.7 CRISPRloci_standalone.py -f Example/NC_005230_proteins.fasta -st  protein
```


#### Viral mode
In order to test the viral more execute the following command
```
python3.7 CRISPRloci_standalone.py -f Example/Input3.fa -st repeat
```


## DNA mode

#### In order to pick this mode the flag `-st` should be set to `dna` 

### General flags

* `-f` : input DNA fasta file path. 

* `-output` folder where results will be stored

* `-cpu` number of CPUs to use

### Cas Protein search related flags

* `-r` : list of regressors to use. Available options: CART, ERT or SVM (default: ERT).

* `-c` : list of classifiers to use. Available options: CART, ERT or SVM (default: ERT).

* `-s` : list of HMM models to use, available options: HMM1 to HMM5 and HMM2019 (default: HMM2019). The models HMM1 to HMM5 are the ones that were originally used in our paper. HMM2019 consists on the HMM models that were obtained from the most recent dataset by [Makarova (2019)](https://www.nature.com/articles/s41579-019-0299-x). Setting this parameter is enough for the tool to know which ML models should be used.

* `-sc`: sequence completeness (used only when `-st` is set to `dna`). Available options: `complete` or `partial` (default: `complete`).

* `-m`: run mode. Available options: `classification`, `regression` or `combined` (default: `combined`).

* `-cg` : maximum number of contiguous gaps allowed in a cassette (default: 1)

* `-cm` :  which ML models to use. Available options: `ERT` or `DNN` (default: `ERT`).

### CRISPR array search related flag

* `--model` Model for the CRISPR array classification. Takes values: 8, 9, 10, ALL and specifies the classification model. The default value is `ALL`

* `--strand` Specifies if the array orientation should be predicted. Available options `True/False`. The default value is `True`

* `--is_element` Specifies if IS-Elements should be predicted. Available options `True/False`. The default value is `False`

* `--fast_run` option to skip the candidate enhancement. Available options `True/False`. The default value is `False`

* `--degenerated` allows search for degenerated repeat candidates on both ends of the CRISPR array candidate. Available options `True/False`. The default value: `True`

* `--min_len_rep` specifies the minimum length of repeats in a CRISPR array. The default value: 21

* `--max_len_rep` specifies the maximum length of repeats in a CRISPR array. The default value: 55

* `--min_len_spacer` specifies the minimum average length of spacers in a CRISPR array. The default value: 18

* `--max_len_spacer` specifies the maximum average length of spacers in a CRISPR array. The default value: 78

* `--min_repeats` specifies the minimum number of repeats in a CRISPR array. The default value: 3

* `--enhancement_max_min` specifies if the filter approximation based on the max. and min. elements should be built
The default value is True 

* `--enhancement_start_end` specifies if the filter approximation based on the max. and min. elements should be built
The default value is True 

* `--max_identical_spacers` specifies the number of maximum identical spacers in a CRISPR array. The default value: 4

* `--max_identical_cluster_spacers` specifies the number of maximum identical consequent spacers in a CRISPR array. The default value: 3

* `--margin_degenerated` specifies the maximum length difference between a new spacer sequence (obtained with the search of degenerated repeats) and the average value of spacer length in the array. The default value: 30

* `--max_edit_distance_enhanced` specifies the number of editing operations for candidate enhancement. The default value: 6

## Protein mode
#### In order to pick this mode the flag `-st` should be set to `protein` 


### General flags

* `-f` : input proteins fasta file path. 

* `-output` folder where results will be stored

* `-cpu` number of CPUs to use

### Cas Protein search related flags

* `-r` : list of regressors to use. Available options: CART, ERT or SVM (default: ERT).

* `-c` : list of classifiers to use. Available options: CART, ERT or SVM (default: ERT).

* `-s` : list of HMM models to use, available options: HMM1 to HMM5 and HMM2019 (default: HMM2019). The models HMM1 to HMM5 are the ones that were originally used in our paper. HMM2019 consists on the HMM models that were obtained from the most recent dataset by [Makarova (2019)](https://www.nature.com/articles/s41579-019-0299-x). Setting this parameter is enough for the tool to know which ML models should be used.

* `-sc`: sequence completeness (used only when `-st` is set to `dna`). Available options: `complete` or `partial` (default: `complete`).

* `-m`: run mode. Available options: `classification`, `regression` or `combined` (default: `combined`).

* `-cg` : maximum number of contiguous gaps allowed in a cassette (default: 1)

* `-cm` :  which ML models to use. Available options: `ERT` or `DNN` (default: `ERT`).

## Viral mode
#### In order to pick this mode the flag `-st` should be set to `virus` 


### General flags

* `-f` : input proteins fasta file path. 

* `-output` folder where results will be stored

* `-cpu` number of CPUs to use

### Search parameter

* `evalue_s` the number of expected hits with spacer database. The default value: 1e-7
