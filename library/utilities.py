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
import xml.dom.minidom
from pathlib import Path

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

""" functiom to call WhatIF web service and store file """
""" modified from WhatIF http://swift.cmbi.ru.nl/whatif/HTML/what-if-web-service-sample-scripts.tgz """
def save_whatif_ws(whatif_service,complex_file_path):
    print("\nWhatIf service: %s" % whatif_service)
    try:
        # read complex
        complex_file=open(complex_file_path,"r")
        complex_file_content=complex_file.read().encode("utf-8")
        complex_file.close()
        # upload PDB file to What IF webservice
        whatif_pdb = urllib.request.urlopen(definitions.WHATIF_UPLOADPDB, complex_file_content)
        whatif_pdb_xml = xml.dom.minidom.parse(whatif_pdb)
        whatif_pdb_id = whatif_pdb_xml.getElementsByTagName("response")[0].childNodes[0].data
        print("WhatIf protein id: %s" % whatif_pdb_id)
        # call web service
        whatif_service_url=definitions.WHATIF_REST+whatif_service+"/id/"+whatif_pdb_id
        print("WhatIf service url: %s\n" % whatif_service_url)
        whatif_response = urllib.request.urlopen(whatif_service_url)
        # store response in file
        whatif_folder=os.path.dirname(complex_file_path)+definitions.FILE_SEPARATOR+"whatif"
        if not os.path.exists(whatif_folder):
            os.mkdir(whatif_folder)
        whatif_file_name=os.path.basename(complex_file_path).replace(".pdb","_"+whatif_service+".xml")
        whatif_file_path=whatif_folder+definitions.FILE_SEPARATOR+whatif_file_name
        whatif_file=open(whatif_file_path,"w")
        whatif_file.write(whatif_response.read().decode("utf-8"))
        whatif_file.close()
        return whatif_file_path
    except Exception as Argument:
        print("An error has occurred \n%s" % Argument)
        raise

""" function to create pdb  file from vina's pdbqt protein and files """
def create_complex_pdb(poses_file_path,receptor_pdbqt_file_path):
    try:
        # obtain docking folder
        docking_folder=os.path.dirname(poses_file_path)
        # open poses file
        poses_file = open(poses_file_path, 'r')
        # write complex pdbqt
        # by copying receptor pdbqt and first pose to a new ligand_receptor pdbqt
        receptor_ligand_pdbqt_file_name = Path(receptor_pdbqt_file_path).name. \
            replace(".pdbqt", "_ligand.pdbqt")
        receptor_ligand_pdbqt_file_path = docking_folder + definitions.FILE_SEPARATOR + \
                                           receptor_ligand_pdbqt_file_name
        receptor_ligand_pdbqt = open(receptor_ligand_pdbqt_file_path, "w")
        with open(receptor_pdbqt_file_path, "r") as receptor_pdbqt:
            for line in receptor_pdbqt:
                receptor_ligand_pdbqt.write(line)
        poses_file.seek(0)
        # change ATOM in ligand.pdbqt file to HETATM
        # assign ligand to chain A
        while True:
            line = poses_file.readline()
            line = line.replace("LIG    1", "LIG A  1")
            "LIG    1"
            if line != "MODEL 2\n":
                if line.count("ATOM      ") == 1:
                    receptor_ligand_pdbqt.write(line.replace("ATOM      ", "HETATM    "))
                elif line.count("ATOM     ") == 1:
                    receptor_ligand_pdbqt.write(line.replace("ATOM     ", "HETATM   "))
                else:
                    receptor_ligand_pdbqt.write(line)
            else:
                break
        # close files
        poses_file.close()
        receptor_ligand_pdbqt.close()
        return receptor_ligand_pdbqt_file_path
    except Exception as Argument:
        print("An error has occurred \n%s" % Argument)
        raise




