# Transmission_DataAnalysis
Project to analyse data of n_TOF natCu Transmission campaign.


## File structure
This project is organized into the following directories and key files:

- `config/` — Files to enable reading data from input files:
  - `ConfigReader.h`: Class to access the input parameters in the different types
  - `configreader_bindings.cpp`: C++ Python binding to allow Python codes to access class ConfigReader
  
- `setup.py`: Builder of the bindings

- `input_files/` — Input files:
  - `Efficiency_plot_runlists.cmnd`: Input parameters for code Efficiency_plot_runlists.py
  - `Histograms_AllRuns.cmnd`: Input parameters for code Histograms_AllRuns.py
  - `Transmission_ratio_final.cmnd`: Input parameters for code Transmission_final.C

- `OUTPUT/` — Folder to store the outputs

- `Histograms_AllRuns.py` : Builds amplitude histogram for each run and allows the plot of a few if needed for inspection
- `STABILITY_allrunsMaxBin.py` : Imports the module that creates the amplitude histograms from Histograms_AllRuns.py and performs the detector stability analysis
- `EFFICIENCY_AllRuns.py` : Imports the module that creates the amplitude histograms from Histograms_AllRuns.py and performs an efficiency study
- `Transmission_final.C` : Creates ToF histograms and Transmission ratio
- `tof_to_E.C` : Converts Transmission graph from time to energy domain and computes cross section<br>
*Note:* in this code you have to change the path to your TTOFSORT<br>
*Note:* this code works only if you have first created a transmission histogram (by running the code Transmission_final.C) with the same number of bins

- `README.md` : Documentation and usage instructions (this file)
 



## Usage

### Importing the repository
To get this repository and enter in it:

```bash
git clone https://github.com/AliBerardi/Transmission_DataAnalysis.git
cd Transmission_DataAnalysis
``` 

### Bash

Setup bindings with:

```bash
python setup.py build_ext --inplace
```
Or:

```bash
pip install .
```

### Examples of running the codes

Run python codes with:
```bash
python3 EFFICIENCY_AllRuns.py
```

Run C++ macros with:
```bash
root -l -b -q 'Transmission_final.C(200)'
```
```bash
root -l -b -q 'tof_to_E(182.1, 200)'
```

