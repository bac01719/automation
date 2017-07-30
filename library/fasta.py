# bola coker
# 21st May, 2017
#
# class to get fasta given NCBI number
from Bio import Entrez
import library.definitions as definitions
import os
import xml.etree.ElementTree
import re
import time
import urllib.request
from pathlib import Path
from Bio.PDB.PDBParser import PDBParser
import Bio.PDB.Selection as Selection
import library.utilities as utilities

class Fasta:
    """ class to obtain FASTA from protien NCBI number"""

    """ method to tranlate 3 letter aa to 1 letter """
    @staticmethod
    def __get_aa_code(amino_acid):
        aa_code = ""
        if str(amino_acid).upper() == "GLY":
            aa_code = "G"
        elif str(amino_acid).upper() == "PRO":
            aa_code = "P"
        elif str(amino_acid).upper() == "ALA":
            aa_code = "A"
        elif str(amino_acid).upper() == "VAL":
            aa_code = "V"
        elif str(amino_acid).upper() == "LEU":
            aa_code = "L"
        elif str(amino_acid).upper() == "ILE":
            aa_code = "I"
        elif str(amino_acid).upper() == "MET":
            aa_code = "M"
        elif str(amino_acid).upper() == "CYS":
            aa_code = "C"
        elif str(amino_acid).upper() == "PHE":
            aa_code = "F"
        elif str(amino_acid).upper() == "TYR":
            aa_code = "Y"
        elif str(amino_acid).upper() == "TRP":
            aa_code = "W"
        elif str(amino_acid).upper() == "HIS":
            aa_code = "H"
        elif str(amino_acid).upper() == "LYS":
            aa_code = "K"
        elif str(amino_acid).upper() == "ARG":
            aa_code = "R"
        elif str(amino_acid).upper() == "GLN":
            aa_code = "Q"
        elif str(amino_acid).upper() == "ASN":
            aa_code = "N"
        elif str(amino_acid).upper() == "GLU":
            aa_code = "E"
        elif str(amino_acid).upper() == "ASP":
            aa_code = "D"
        elif str(amino_acid).upper() == "SER":
            aa_code = "S"
        elif str(amino_acid).upper() == "THR":
            aa_code = "T"
        else:
            aa_code = "-"
        return aa_code

    """ method to return amino acid sequence with coords in pdf file """
    @staticmethod
    def __get_pdbseq(protein, protein_file_path):
        pdb_object = PDBParser(PERMISSIVE=True, QUIET=True)
        pdb_structure = pdb_object.get_structure(protein, protein_file_path)
        pdb_fasta = ""
        for residue in Selection.unfold_entities(pdb_structure, "R"):
            # print(residue.get_full_id())
            if residue.get_full_id()[3][0] == " " and residue.get_full_id()[2] == "A":
                pdb_fasta += Fasta.__get_aa_code(residue.get_resname())
                # print("-- %s , %s \n" % (residue.get_resname(),get_aa_code(residue.get_resname())))
        return pdb_fasta

    """method to obtain fasta given list of proteins and working directory"""
    @staticmethod
    def get_fasta(protein_list,home_folder,snp=""):
        try:
            # check if home folder exist if not give error
            os.chdir(home_folder)
            local_time=time.asctime( time.localtime(time.time()))
            local_time=re.sub(r'\W+', '_', local_time)
            snp_extension=""
            if snp!="":
                snp_extension=snp+"_"
            fasta_folder=home_folder.rstrip(definitions.FILE_SEPARATOR)+\
                         definitions.FILE_SEPARATOR+snp_extension+'analysis_'+local_time
            fasta_sequence=''
            if not os.path.exists(fasta_folder):
                os.mkdir(fasta_folder)
            fasta_list=[]
            # get fasta for each protein
            for protein in protein_list:
                Entrez.email = definitions.EMAIL_ADDRESS
                handle = Entrez.esearch(db=definitions.PROTEIN_DATABASE,
                                        term=protein[definitions.DICT_PROTEIN_NCBI_ID])
                record = Entrez.read(handle)
                for id_number in record["IdList"]:
                    # get full record detail in xml format for each id
                    id_handle = Entrez.efetch(db="protein", id=id_number, rettype="fasta", retmode="xml")
                    fasta_XML = id_handle.read()
                    # parse xml produced
                    root = xml.etree.ElementTree.fromstring(fasta_XML)
                    for tseq in root.findall('TSeq'):
                        fasta_sequence=tseq.find('TSeq_sequence').text
                    # add single point mutation if given
                    if snp!="":
                        parse=re.compile("(.)(\d+)(.)")
                        snp_items=parse.match(snp)
                        fasta_sequence=fasta_sequence.replace(snp_items.groups()[0],\
                                                              snp_items.groups()[2],\
                                                              int(snp_items.groups()[1]))
                    # create fasta_protein subfolder if not exist
                    protein_subfolder =re.sub(r'\W+', '_', protein[definitions.DICT_PROTEIN_NAME])
                    fasta_protein_folder=fasta_folder+definitions.FILE_SEPARATOR+protein_subfolder
                    if not os.path.exists(fasta_protein_folder):
                        os.mkdir(fasta_protein_folder)
                    # create fasta_protein_sequence subfolder if not exist
                    sequence_subfolder=protein[definitions.DICT_PROTEIN_NCBI_ID]+'_'+str(id_number)
                    fasta_protein_sequence_folder=fasta_protein_folder+definitions.FILE_SEPARATOR\
                                                  +sequence_subfolder
                    if not os.path.exists(fasta_protein_sequence_folder):
                        os.mkdir(fasta_protein_sequence_folder)
                    # store protein sequence in file
                    protein_sequence_file_path=fasta_protein_sequence_folder+definitions.FILE_SEPARATOR+ \
                                               protein[definitions.DICT_PROTEIN_NCBI_ID]+'_'+str(id_number)\
                                               +definitions.FASTA_FILE_EXTENSION
                    protein_sequence_file=open(protein_sequence_file_path,'w')
                    protein_sequence_file_header='>'+protein[definitions.DICT_PROTEIN_NCBI_ID]+'_'+str(id_number)
                    protein_sequence_file.write(protein_sequence_file_header+"\n")
                    protein_sequence_file.write(fasta_sequence+"\n")
                    protein_sequence_file.close()
                    handle.close()
                    fasta_list.append({definitions.DICT_FASTA_SEQUENCE:fasta_sequence,\
                                       definitions.DICT_FASTA_FILE_PATH:protein_sequence_file_path})
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise
        else:
            return {definitions.DICT_FASTA_LIST:fasta_list,\
                    definitions.DICT_FASTA_PROTEIN_FOLDER:fasta_protein_sequence_folder}

    """method to obtain template for homology modelling given a set of PDB names
       in comma separated string"""
    @staticmethod
    def get_templates(template_protein_names,fasta_item):
        try:
            # check if fast folder exist otherwise give error
            path_object=Path(fasta_item[definitions.DICT_FASTA_FILE_PATH])
            fasta_folder=str(path_object.parent)
            os.chdir(fasta_folder)
            # make homology folder
            homology_folder=fasta_folder+definitions.FILE_SEPARATOR+definitions.HOMOLOGY_FOLDER
            os.mkdir(homology_folder)
            # open template file
            template_file_name = path_object.name. \
                replace(definitions.FASTA_FILE_EXTENSION, definitions.HOMOLOGY_TEMPLATE_FILE_EXTENSION)
            template_file_path = homology_folder + definitions.FILE_SEPARATOR + template_file_name
            template_file = open(template_file_path, "w")
            # download proteins in template_protein_names list and update template file
            print("\nDownloading PDBs for template\n")
            print(template_protein_names)
            for protein in template_protein_names:
                #print(protein)
                #protein_url="https://files.rcsb.org/view/"+protein+".pdb"
                #print(protein_url)
                #protein_res=urllib.request.urlopen(protein_url)
                #protein_pdb=protein_res.read().decode("utf-8")
                #protein_file_path=homology_folder+definitions.FILE_SEPARATOR+protein+".pdb"
                #protein_file=open(protein_file_path,"w")
                #protein_file.write(protein_pdb)
                #protein_file.close()
                protein_file_path=utilities.download_protein(protein,homology_folder)
                template_sequence_title = ">" + protein + "|A"
                template_sequence_body = Fasta.__get_pdbseq(protein,protein_file_path)
                template_file.write(template_sequence_title + "\n")
                template_file.write(template_sequence_body + "\n")
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise
        else:
            template_file.close()
            return template_file_path






