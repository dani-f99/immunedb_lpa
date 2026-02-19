--------------------------------------------------------------------------------
                              SYSTEM IMMUNOLOGY LAB
                                Haifa University
                                 Daniel Fridman 
                                      2026
--------------------------------------------------------------------------------


# Identification of Genetic Motifs in SARS-CoV-2 Specific B-Cell Receptor Repertoires Repertoiry

## 1. OVERVIEW
LPA-based detection of enriched and depleted BCR sequence properties. This tool utilizes the LPA (latent 
personal analysis) algorithm to analyze substitution mutation survival fractions and trimer (3-mer) usage 
patterns within the heavy chain variable (V-H) gene of B-cell receptors.


## 2. PREREQUISITES
Please ensure the following python modules are installed:
- `Pandas`
- `NumPy` 
- `SciPy`
- `tqdm`
- `sqlalchemy`


## 3. USAGE GUIDE
1. Congifgure the `congif.json` file (see section 4 - config) by:
   - rename `congif.json.example` to `congif.json`.
   - Fill the sql connection and database information as illusrated.
   - Once the custom python modules will be loaded the script will initiate
   the required folders and import the `congif.json` information into `config`
   variable.
2. Follow the steps in `notebook_tutorial.ipynb` file.

`Note: Run report will be saved in the reports folder`


## 4. CONFIG.JSON CONFIGURATION
The `config.json` is in a json format and it's purpose is to configure the ImmuneDB MySQL connection
and database on which the process will be performed.
before program usage:
- `sql`: Configure the sql connection information.
  - `adress` - ip adress of the MySQL server.
  - `port` - port of the MySQL server.
  - `username` - username credentials used to connect to the MySQL server.
  - `password` - password credentials used to connect to the MySQL server.
  - `database`- database name in the MySQL server. 
  - `metadata_columns` - Required metadata columns.


## 5. Pipeline Steps & Components
1. MySQL tables import.
2. Creation of required dataframes.
3. filtring of out unwatned data (non-functional clones, clones with 1 sequence etc).
4. Creating LPA appropriate input dataset.
5. Performing LPA analysis.
6. Performing PCA and on the LPA domain.
7. Visualization of the results (heatmap and scatterplot).


## 5. DIRECTORY STRUCTURE
 * data_raw -> Folder which contains the raw tables downloaded from the sql server.
 * data_processed -> Processed tables used for the analysis.
   * data_processed/raw_processed -> Modified raw tables, to be turned into the LPA input.
   * data_processed/lpa_input -> LPA input file.
 * reports -> Run reports, consult incase of error.
 * results_figures -> Figure output folder.
 * results_tables -> Tables results folder.
 * source -> Source code.


## 6. RESOUCES
* ImmuneDB GitHub - https://github.com/arosenfeld/immunedb
* LPA GitHub - https://github.com/ScanLab-ossi/LPA


