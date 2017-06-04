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

class Fasta:
    """ class to obtain FAST from protien NCBI number"""

    """method to obtain fasta given list of proteins and working directory"""
    @staticmethod
    def get_fasta(protein_list,home_folder):
        try:
            # check if home folder exist if not give error
            os.chdir(home_folder)
            local_time=time.asctime( time.localtime(time.time()))
            local_time=re.sub(r'\W+', '_', local_time)
            fasta_folder=home_folder.rstrip(definitions.FILE_SEPARATOR)+\
                         definitions.FILE_SEPARATOR+'analysis_'+local_time
            fasta_sequence=''
            if not os.path.exists(fasta_folder):
                os.mkdir(fasta_folder)
            fasta_file_path_list=[]
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
                    protein_sequence_file.write(fasta_sequence)
                    protein_sequence_file.close()
                    handle.close()
                    fasta_file_path_list.append(protein_sequence_file_path)
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise
        else:
            return fasta_file_path_list

    """method to obtain template for homology modelling given a set of PDB names
       in comma separated string"""
    @staticmethod
    def get_templates(template_protein_names,fasta_file_path):
        try:
            # check if fast folder exist otherwise give error
            path_object=Path(fasta_file_path)
            fasta_folder=str(path_object.parent)
            #fasta_folder=definitions.FILE_SEPARATOR.join(fasta_file_path.split(definitions.FILE_SEPARATOR)[0:10])
            os.chdir(fasta_folder)
            # make homology folder
            homology_folder=fasta_folder+definitions.FILE_SEPARATOR+definitions.HOMOLOGY_FOLDER
            os.mkdir(homology_folder)
            # populate temp
            # late file with proteins
            #template_file_name=fasta_file_path.split(definitions.FILE_SEPARATOR)[11].\
            #    replace(definitions.FASTA_FILE_EXTENSION,definitions.HOMOLOGY_TEMPLATE_FILE_EXTENSION)
            template_file_name=path_object.name.\
                replace(definitions.FASTA_FILE_EXTENSION,definitions.HOMOLOGY_TEMPLATE_FILE_EXTENSION)
            template_file_path=homology_folder+definitions.FILE_SEPARATOR+template_file_name
            template_file=open(template_file_path,"w")
            # call RCSB custom report web service
            pdb_url_sequence_string = "http://www.rcsb.org/pdb/rest/customReport.xml?pdbids=" + \
                                    template_protein_names + "&customReportColumns=sequence"
            pdb_response = urllib.request.urlopen(pdb_url_sequence_string)
            sequence_XML = pdb_response.read().decode("utf-8")
            #print(sequence_XML)
            root = xml.etree.ElementTree.fromstring(sequence_XML)
            templates_string=""
            for record in root.findall("record"):
                #print(record.find("dimEntity.chainId").text)
                if record.find("dimEntity.chainId").text==definitions.HOMOLOGY_MONOMER_CHAIN:
                    template_sequence_title=">" + record.find("dimEntity.structureId").text \
                                                + "|" + record.find("dimEntity.chainId").text
                    template_sequence_body=record.find("dimEntity.sequence").text
                    template_file.write(template_sequence_title+"\n")
                    template_file.write(template_sequence_body+"\n")
                    templates_string=templates_string+template_sequence_title+"\n"+template_sequence_body+"\n"
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise
        else:
            template_file.close()
            return template_file_path
