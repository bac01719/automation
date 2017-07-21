# bola coker
# 9th June, 2017
#
# utilities.py
# includes general functions
import urllib.request
import library.definitions as definitions
from Bio import Entrez
import os
import re

""" function to return a string given a dictionary """
def return_string(dictionary_list):
    dictionary_string=""
    for key_value in dictionary_list:
        if dictionary_list[key_value]!="":
            dictionary_string+="{} ".format(key_value)+dictionary_list[key_value]
        elif dictionary_list[key_value]=="":
            dictionary_string+=key_value
        dictionary_string+=" "
    return dictionary_string

""" function to check opal job is still running """
def opal_job_running(opal_status_response):
    if opal_status_response["code"] == definitions.OPAL_SUCCESS or \
                    opal_status_response["code"] == definitions.OPAL_FAILURE:
        return True
    else:
        return False

""" get opal server returned error message """
def get_error(opal_response):
    try:
        # get opal error link and read contents
        opal_error=urllib.request.urlopen(opal_response["stdErr"])
        opal_error_contents=opal_error.read()
        opal_error.close()
        return opal_error_contents.decode("utf-8")
    except Exception as Argument:
        print("An error has occurred \n%s" % Argument)
        raise

""" store stdOut and stdErr files """
def save_std_files(opal_response,parent_folder):
    try:
        # check parent folder exist. if not create
        if not os.path.exists(parent_folder):
            os.mkdir(parent_folder)
        # read stdout
        opal_stdout = urllib.request.urlopen(opal_response['stdOut'])
        opal_stdout_contents = opal_stdout.read()
        opal_stdout.close()
        # save stdout results to file
        stdout_file_path = parent_folder + definitions.FILE_SEPARATOR + "stdOut.txt"
        stdout_file = open(stdout_file_path, "w")
        stdout_file.write(opal_stdout_contents.decode("utf-8"))
        stdout_file.close()
        # read stderr
        opal_stderr = urllib.request.urlopen(opal_response['stdErr'])
        opal_stderr_contents = opal_stderr.read()
        opal_stderr.close()
        # save stderr results to file
        stderr_file_path = parent_folder + definitions.FILE_SEPARATOR + "stdErr.txt"
        stderr_file = open(stderr_file_path, "w")
        stderr_file.write(opal_stderr_contents.decode("utf-8"))
        stderr_file.close()
    except Exception as Argument:
        print("An error has occurred \n%s" % Argument)
        raise

""" function to download protein given pdb number and folder to download into """
def download_protein(protein,protein_folder_path):
    try:
        print(protein)
        protein_url = "https://files.rcsb.org/view/" + protein + ".pdb"
        print(protein_url)
        protein_res = urllib.request.urlopen(protein_url)
        protein_pdb = protein_res.read().decode("utf-8")
        protein_file_path = protein_folder_path + definitions.FILE_SEPARATOR + protein + ".pdb"
        protein_file = open(protein_file_path, "w")
        protein_file.write(protein_pdb)
        protein_file.close()
        return protein_file_path
    except Exception as Argument:
        print("An error has occurred \n%s" % Argument)
        raise

""" function to check organism exist """
def verify_entrez_option(database,option):
    verified=False
    try:
        Entrez.email = definitions.EMAIL_ADDRESS
        handle = Entrez.esearch(db=database, term=option)
        record = Entrez.read(handle)
        for id_number in record["IdList"]:
            # get full record detail in xml format for each id
            id_handle = Entrez.efetch(db=database, id=id_number, retmode="xml")
            if id_handle.read()!="":
                verified=True
    except Exception as Argument:
        print("An error has occurred \n%s" % Argument)
        raise
    finally:
        return verified

""" function to check SMILES string exist """
def verify_zinc_smiles(smiles):
    verified=False
    try:
        smiles_URL=definitions.ZINC_URL+\
                   "?structure.smiles="+re.escape(smiles)+"&structure.similarity=1.0"
        smiles_response = urllib.request.urlopen(smiles_URL)
        #print(smiles_response.read().decode("utf-8"))
        if smiles_response.read().decode("utf-8")!="":
            verified=True
    except Exception as Argument:
        return verified
        #print("An error has occurred \n%s" % Argument)
        #raise
    finally:
        return verified

""" function to parse SNP """
def verify_snp(snp):
    verified=False
    try:
        parse=re.compile(definitions.SNP_REG_EXP)
        if parse.match(snp)==None:
           verified=True
    except Exception as Argument:
        print("An error has occurred \n%s" % Argument)
        raise
    finally:
        return verified

""" function to verify mutation file """
def verify_mutation_file(mutation_file_path):
    parse=re.compile("^.\,.\d+.$")
    verified=False
    try:
        with open(mutation_file_path,"r") as mutations:
            for line in mutations:
                if parse.match(line)==None:
                    verified=False
                    break
                else:
                    verified=True
    except Exception as Argument:
        print("An error has occurred \n%s" % Argument)
        raise
    finally:
        return verified




