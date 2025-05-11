# WSL Setup Guide

## Install Windows Subsystem for Linux

1. Open PowerShell as Administrator (right-click in search -> open as administrator)
2. Run the following command:
   ```powershell
   dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
   ```
3. Run the following command:
   ```powershell
   dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
   ```
4. Restart your PC

## Install Ubuntu

1. Open Microsoft Store and find Ubuntu
2. Install -> Launch

## How to Use Ubuntu (WSL)

1. Open PowerShell 
2. Run: `wsl`
3. Your command line should start looking like this: 
   ```
   derfe@DESKTOP-D4CC88M:/mnt/c/Users/derfe$
   ```

## Install Conda

Run the following commands:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b -p $HOME/miniconda
eval "$($HOME/miniconda/bin/conda shell.bash hook)"
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc
echo 'eval "$(conda shell.bash hook)"' >> ~/.bashrc
source ~/.bashrc
```

## Setting Up RrmScorer Tool

1. Navigate to working folder:
   ```bash
   cd ~
   ```

2. Get RrmScorer:
   ```bash
   git clone https://bitbucket.org/bio2byte/rrmscorer.git
   cd rrmscorer
   ```

3. Install environment and all libraries:
   ```bash
   conda create --yes --name rrmscorer python=3.10.4
   conda activate rrmscorer
   conda config --add channels conda-forge
   conda install --yes numpy=1.21.5 pandas=1.4.2 biopython=1.79 matplotlib=3.5.2 scikit-learn=1.1.1
   conda install --yes -c bioconda hmmer=3.3.2
   ```

4. Test that everything works:
   ```bash
   python rrm_rna_wrapper.py -RRM P19339_RRM1 -RNA UAUAUUAGUAGUA -ws 5
   ```

5. Open in explorer `\\wsl.localhost\Ubuntu\home` and find rrmScorer folder. Add analyzer and data files:
   - `rrm_analyzer.py`
   - `extract_data.py`
   - `journal.pcbi.1010859.s010.txt`
   - `journal.pcbi.1010859.s011.txt`
   - `proten_cds.fasta`

## How to Use Analyzer

1. Generate extracted_data.txt:
   ```bash
   python extract_data.py journal.pcbi.1010859.s011.txt proten_cds.fasta
   ```

2. Run the analyzer with different options:

   ```bash
   # Basic analysis
   python rrm_analyzer.py extracted_data.txt output
   
   # For generating images connections vs random
   python rrm_analyzer.py extracted_data.txt output_for_plots --batch-plot-distributions --iterations 5000
   
   # For generating images connections vs random INCLUDING random RRMs and HOMO SAPIENS ORGANISMS
   python random_distribution.py extracted_data.txt
   
   # Images for scores + full - dataset mean
   python rrm_analyzer.py extracted_data.txt output_dir --batch-binding-comparison --window-size 5
   ```

   Where:
   - `extracted_data.txt` - generated file
   - `output` - name of output folder for score files

## Results

1. plot_output - All plots with random distributions for current RRM inside the range + all p-values plot
2. random_distribution_plots - All plots with random distrtbutions for random RRMs + Homo sapiens filter + all p-vales plots
3. 1024_random_distribution - All plots with distributions for random windows from 1024 possible + all p-values plot 
