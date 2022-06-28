# RRMScorer documentation
RRMScorer allows the user to easily predict how likely a single RRM is to bind ssRNA using a carefully generated alignment for the RRM structures in complex with RNA, from which we analyzed the interaction patterns and derived the scores (Please address to the publication for more details on the method REF)

---

## Installation

#### Clone this repository to your working environment:
```console
$ git clone git@bitbucket.org:bio2byte/rrmscorer.git && cd rrmscorer
```

#### The following packages are required:

```python
python==3.10.4
numpy==1.21.5
pandas==1.4.2
biopython==1.79
matplotlib==3.5.2
scikit-learn==1.1.1
hmmer==3.3.2
```
#### Via [Conda](https://docs.conda.io/en/latest/):

```console
$ conda create --yes --name rrmscorer python==3.10.4
$ conda activate rrmscorer
$ conda install --yes --file requirements.txt
```

#### Via [Virtual Environment](https://docs.python.org/3/tutorial/venv.html):

```console
$ python3 -m venv rrmscorer-venv
$ source ./rrmscorer-venv/bin/activate
$ python -m pip install numpy==1.21.5 pandas==1.4.2 biopython==1.79 matplotlib==3.5.2 scikit-learn==1.1.1 hmmer
```

## How to run it:
Either you are using Conda or Virtual Environments for your installation, before executing this software features, you need to setup the Python environment.
Using Conda:

```console
$ conda activate rrmscorer
```
Using Virtual Environment:

```console
$ source ./rrmscorer-venv/bin/activate
```

Continue reading the next section to find further details about the available features.
In case you need to deactivate this Python environment:

Using Conda:

```console
$ conda deactivate
```

Using Virtual Environment:

```console
$ deactivate
```

## Features
RRMScorer has several features to either calculate the binding score for a specific RRM and RNA sequences, for a set of RRM sequences in a fasta file, or to explore which are the best RNA binders according to our scoring method.

### i) Single RRM vs RNA
To use this feature the user needs to input:

1. `-RRM` The identifier of an RRM already included in the original alignment following the format UniProtId_RRMnumber (e.g., P19339_RRM1)
2. `-RNA` The RNA sequence to test
3. `-ws` The window size to test (**Only 3 and 5 nucleotide windows are accepted**)
4. `-plot` [Optional] To generate score plots for all the RNA possible windows


```console
$ python rrm_rna_wrapper.py -RRM P19339_RRM1 -RNA UAUAUUAGUAGUA -ws 5 [-plot]
```

Example output:
```console
UAUAU -1.08
AUAUU -0.99
UAUUA -1.33
AUUAG -0.90
UUAGU -1.07
```

### ii) Fasta file with RRM sequences vs RNA
To use this feature the user needs to input:

1. `-fasta` Fasta file with 1 or more single RRM sequences. The sequences are aligned to the master alignment HMM
1. `-RNA` The RNA sequence to test
1. `-ws` The window size to test (**Only 3 and 5 nucleotide windows are accepted**)
1. `-plot` [Optional] To generate score plots for all the RNA possible windows

```console
$ python rrm_rna_wrapper.py -fasta input_files/rrm_seq.fasta -RNA UAUAUUAGUAGUA -ws 5 [-plot]
```


### iii) Fasta file / RRM id to find top-scoring RNAs
To use this feature the user needs to input:

1. `-fasta` Fasta file or RRM is as described in the previous cases.
1. `-ws` The window size to test (**Only 3 and 5 nucleotide windows are accepted**)
1. `-top` To find the top-scoring RNA for the specified RRM/s

```console
$ python rrm_rna_wrapper.py -fasta input_files/rrm_seq.fasta -ws 5 -top
```



