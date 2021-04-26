# CRISPRloci
####CRISPRloci provides an automated and comprehensive in silico characterization of CRISPR-Cas system on bacterial and archaeal genomes. It is a full suite for CRISPR locus characteriztion that includes CRISPR array orientation, detection of conserved leaders, cas gene annotation and subtype classification.

The web server interface of CRISPRloci is freely available at: rna.informatik.uni-freiburg.de/trunk/CRISPRloci
 
 
![GitHub Logo](/images/Webserver_flowchart.png)


### Quick run with the default parameters

#### DNA mode
In order to test the dna mode please execute the following the command 

```
python3.7 CRISPRloci_standalone.py -f Example/NC_005230.fasta -st dna
```

#### Protein mode
In order to test the protein mode tun the following command:
```
python3.7 master_script.py -f Example/NC_005230_proteins.fasta -st  protein
```


#### Viral mode
In order to test the viral more execute the following command
```
python3.7 master_script.py -f Example/Input3.fa -st repeat
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


#All flags (tmp remove)


* `-f` : input DNA fasta file path. 

* `-r` : list of regressors to use. Available options: CART, ERT or SVM (default: ERT).

* `-c` : list of classifiers to use. Available options: CART, ERT or SVM (default: ERT).

* `-s` : list of HMM models to use, available options: HMM1 to HMM5 and HMM2019 (default: HMM2019). The models HMM1 to HMM5 are the ones that were originally used in our paper. HMM2019 consists on the HMM models that were obtained from the most recent dataset by [Makarova (2019)](https://www.nature.com/articles/s41579-019-0299-x). Setting this parameter is enough for the tool to know which ML models should be used.

* `-sc`: sequence completeness (used only when `-st` is set to `dna`). Available options: `complete` or `partial` (default: `complete`).

* `-m`: run mode. Available options: `classification`, `regression` or `combined` (default: `combined`).

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

* `-cg` : maximum number of contiguous gaps allowed in a cassette (default: 1)

* `-cm` :  which ML models to use. Available options: `ERT` or `DNN` (default: `ERT`).

* `-output` folder where results will be stored

* `-cpu` number of CPUs to use

* `evalue_r` the number of expected hits with repeat database

* `evalue_s` the number of expected hits with spacer database