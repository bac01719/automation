# bola coker
# 21st May, 2017
#
# class to get fasta given NCBI number
from Bio import Entrez
import library.definitions
import os
import xml.etree.ElementTree
import re

class Fasta:
    """ class to obtain FAST from protien NCBI number"""

    """method to obtain fasta given protein dictionary and working directory"""
    @staticmethod
    def get_fasta(protein_list,home_folder):
        try:
            # check if home folder exist if not give error
            os.chdir(home_folder)
            fasta_folder=home_folder.rstrip(library.definitions.FILE_SEPARATOR)+\
                         library.definitions.FILE_SEPARATOR+'fasta'
            fasta_sequence=''
            if not os.path.exists(fasta_folder):
                os.mkdir(fasta_folder)
            for protein in protein_list:
                Entrez.email = library.definitions.EMAIL_ADDRESS
                handle = Entrez.esearch(db=library.definitions.PROTEIN_DATABASE,
                                        term=protein[library.definitions.DICT_PROTEIN_NCBI_ID])
                record = Entrez.read(handle)
                for id_number in record["IdList"]:
                    # get full record detail in xml format for each id
                    id_handle = Entrez.efetch(db="protein", id=id_number, rettype="fasta", retmode="xml")
                    fasta_XML = id_handle.read()
                    # parse xml produced
                    root = xml.etree.ElementTree.fromstring(fasta_XML)
                    for tseq in root.findall('TSeq'):
                        fasta_sequence=tseq.find('TSeq_sequence').text
                    # create fasta_protein subfolder if not exist
                    protein_subfolder =re.sub(r'\W+', '_', protein[library.definitions.DICT_PROTEIN_NAME])
                    fasta_protein_folder=fasta_folder+library.definitions.FILE_SEPARATOR+protein_subfolder
                    if not os.path.exists(fasta_protein_folder):
                        os.mkdir(fasta_protein_folder)
                    # create fasta_protein_sequence subfolder if not exist
                    sequence_subfolder=protein[library.definitions.DICT_PROTEIN_NCBI_ID]+'_'+str(id_number)
                    fasta_protein_sequence_folder=fasta_protein_folder+library.definitions.FILE_SEPARATOR\
                                                  +sequence_subfolder
                    if not os.path.exists(fasta_protein_sequence_folder):
                        os.mkdir(fasta_protein_sequence_folder)
                    # store protein sequence in file
                    protein_sequence_file_path=fasta_protein_sequence_folder+library.definitions.FILE_SEPARATOR+ \
                                               protein[library.definitions.DICT_PROTEIN_NCBI_ID]+'_'+str(id_number)\
                                               +'.fasta'
                    protein_sequence_file=open(protein_sequence_file_path,'w')
                    protein_sequence_file.write(fasta_sequence)
                    protein_sequence_file.close()
                    handle.close()
                    return protein_sequence_file_path
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise