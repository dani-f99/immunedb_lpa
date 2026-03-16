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

        