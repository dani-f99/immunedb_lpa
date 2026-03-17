#########
# Imports
from sqlalchemy import create_engine
import pandas as pd
import numpy as np
import json
import os


####################################################################
# DNA codon table - used for the translation from nt to aa sequences
protein = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
            "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
            "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
            "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
            "TCT" : "s", "CCT" : "P", "ACT" : "T", "GCT" : "A",
            "TCC" : "s", "CCC" : "P", "ACC" : "T", "GCC" : "A",
            "TCA" : "s", "CCA" : "P", "ACA" : "T", "GCA" : "A",
            "TCG" : "s", "CCG" : "P", "ACG" : "T", "GCG" : "A",
            "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
            "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
            "TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
            "TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
            "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
            "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
            "TGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
            "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G", 
            "---" : "-"
            }

########################################################
# Codons dictionary used for the nt sequence translation
codon_dic_updated = {
                'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                'TCT': "S'", 'TCC': "S'", 'TCA': "S'", 'TCG': "S'",
                'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',  # * for STOP
                'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',

                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGT': 'S"', 'AGC': 'S"', 'AGA': 'R', 'AGG': 'R',

                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
                }


#########################################
# custom function to round numbers upward
def round_up(number):
    num_dec = number
    num_round = round(number)
    
    if num_round < num_dec:
        value = num_round + 1
    else:
        value = num_round
    return value


##################################################################
# Custom function that deletes selected files from specific folder
def pruge_file(dir_path : str,
               files_2keep : list):
    """
    dir_path : str -> absulote path of the dir to be pruged
    files_2keep : list -> list of files to keep in the dir_path folder 
    """
    dir_list = os.listdir(dir_path)
    files_2delete = [os.path.join(dir_path, i) for i in dir_list if i not in files_2keep]

    for i in files_2delete:
        os.remove(i)

    print(f"> folder {dir_path} purged from the following files: {files_2delete}")


################################################################################$$##########
# Reading information from json file. Used to extract the parameters from the `config.json`.
def read_json(path:str = "config.json") -> dict:
    """
    path : str -> path of the json file
    """

    with open(path) as config:
        config_f = json.load(config)

    return config_f


###################################################################################################
# Custom function that cheecks if dit exists and if not create it (require absulote path as input).
def create_folders(dir_list : list):
    """dir_list : list -> list of dir to be created, if dirs exists than pass"""

    if isinstance(dir_list, (list, np.ndarray)) is False:
        raise Exception("> 'dir_list' argument isnt array-like please inpur list argument.")

    dirs_exists = [os.path.exists(i) for i in dir_list]

    for i in range(len(dirs_exists)):
        dir = dir_list[i]

        if dirs_exists[i] is False:
            os.mkdir(dir_list[i])
            print(f"> Dir `{dir}` was created.")
        else:
            print(f"> Dir `{dir}` already exists.")


######################################################################################################
# A custom function that connects to MySQL server, execute query and returns the results as dataframe.
class mysql_qry():
    def __init__(self,
                username:str = read_json()["sql"]["username"],
                password:str= read_json()["sql"]["password"],
                adress:str= read_json()["sql"]["adress"],
                port:str= read_json()["sql"]["port"],
                db_name:str = read_json()["sql"]["database"]):

        """
        username:str -> Username credentials for the MySQL server.
        password:str -> Password credentials for the MySQL server.
        adress:str -> IP adress of the MySQL server.
        port:str -> Port of the MySQL server.
        """

        # Setting up MySQL connenction
        connection_mysql = f"mysql+pymysql://{username}:{password}@{adress}:{port}/{db_name}"
        self.engine = create_engine(connection_mysql)
        print(f"> Established connecntion to the {db_name} database.")

    def run_qry(self,
                qry:str) -> pd.DataFrame:
        """qry : str -> SQL query to be executed.)"""

        # Executing the query
        qry_df = pd.read_sql(qry, self.engine)
        print(f"> Query executed successfully: \n{qry}")

        # Returing the table as pd.dataframe
        return qry_df

    def close_conn(self):
        # Closing the connenction
        self.engine.dispose()
        print("> MySQL connenction terminated.")


########################################################
# Translating NT sequense and creating trimers dataframe
def nt_transalte_104(cdr_seq:str, 
                     aa_start:int = 1, 
                     aa_end:int = None):
                     #returns:str = "trimers",
                     #output_type:str = "pandas"):
    
    """
    cdr_seq: str -> cdr_seq.dna_seq format.
    aa_start: int -> first amino acid position.
    aa_end: int -> last aa_position.
    #returns: str -> non-functional argument.
    #output_type: str -> non functional argument.
    """
    
    # heavy chain varaible regions, see IMGT documentation for more information.
    regions_dict = {"fw1": np.arange(1,27),
                    "cdr1": np.arange(27,39),
                    "fw2": np.arange(39,56),
                    "cdr2": np.arange(56,66),
                    "fw3": np.arange(66,105)}
    
    # the first column should be the cdr3_aa and the second germline\seq
    cdr_seq = cdr_seq.split(".") #sequence of the cdr, changes between sequences.
    cdr3_len = len(cdr_seq[0]) #cdr3 length
    nt_seq = cdr_seq[1] #1-104 sequence

    translated = []
    t_length = int((len(nt_seq)-len(nt_seq)%3)/3)

    #if the number of NT spacers in the NT sequence isn't equal to 3 there is a change in 
    #Reading frame and need to cheek the sequence 
    n_spaces = nt_seq.count("-")
    if  n_spaces % 3 != 0:
        raise Exception("Number of NT spacers dosent divide by 3, cheek sequence") 

    # Translating the NT sequence
    for i in range(1,t_length):
        codon = nt_seq[i*3-3:i*3]
      
        if codon in list(codon_dic_updated.keys()):
            aa = codon_dic_updated[codon]  
        else:
            aa = "-"
         
        translated.append(aa)

    # If no end was decided (aa_end argument) the code will itirate over all the dna sequence
    if aa_end is None:
        aa_end = t_length

    results_aa = pd.Series(data=translated[aa_start-1:aa_end], index=range(aa_start,aa_end))
    results_aa_cleaned = results_aa[results_aa != "-"]
    
    # output dataframe creation
    r_dict = {"first_aa":[], "trimer":[], "ncount":1}
    j_count = 1
    j_positions = []
    for j in range(1, len(results_aa_cleaned)-1):
        r_dict["trimer"].append("".join(results_aa_cleaned.values[j-1:j+2]))

        # giving different numring scheme to the J region beacuse the CDR3 is highly variable in it's length
        first_aa = results_aa_cleaned.index[j-1]
        seq_up2_j = 104 + cdr3_len

        if first_aa > seq_up2_j:
            r_dict["first_aa"].append(f"{j_count}j")
            j_positions.append(f"{j_count}j")
            j_count += 1

        else:
            r_dict["first_aa"].append(first_aa)

    trimer_result = pd.DataFrame(r_dict)

    # determining the region of cdr3 and j based on the cdr3 length
    regions_dict["cdr3"] =  np.arange(105,105+cdr3_len)
    regions_dict["j"] = j_positions

    # assigning each trimer it's location
    reg_list = []
    for value in trimer_result.first_aa.values:
        for key in regions_dict.keys():
            if value in regions_dict[key]:
                reg_list.append(key)
    
    trimer_result.insert(loc=1, 
                         column="region", 
                         value= reg_list)

    return trimer_result #, reg_list # results_aa_cleaned