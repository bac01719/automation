# bola coker
# class to handle blast requests
# date: 12st May, 2017
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
import urllib.request
import urllib.parse
import xml.etree.ElementTree
import library.definitions as definitions
import os

class Blast:
    """ class to handle blast function"""

    """
    method to return dictionary of blast results as well as save to disk.
    input is fasta file path, percentage identity cutoff for proteins to be docked
    or complexes found that can go straight to data analysis, InChIKey of ligand    
    """
    @staticmethod
    def do_blast(protein_sequence_file_path,percentage_identity_cutoff,inchi_key):
        try:
            # open fasta file and read it
            fasta_file=open(protein_sequence_file_path,'r')
            fasta_sequence_header=fasta_file.readline()
            fasta_sequence=fasta_file.readline()
            fasta_file.close()
            # create file to include blast results
            blast_file_path=protein_sequence_file_path.replace(definitions.FASTA_FILE_EXTENSION,\
                                                               definitions.BLAST_FILE_EXTENSION)
            blast_file=open(blast_file_path,'w')
            # do blast
            result_handle = NCBIWWW.qblast(definitions.BLAST_PROGRAM, 'pdb', fasta_sequence,
                                           perc_ident=percentage_identity_cutoff)
            blast_records = NCBIXML.parse(result_handle)
            # initialise blast list to return results
            blast_list = []
            # write out blast records to disk and populate blast_list
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        percentage_identity=hsp.identities*100/len(fasta_sequence)
                        if percentage_identity > percentage_identity_cutoff:
                            # write blast into file
                            blast_file.write('****Blast results****\n\n')
                            blast_file.write('sequence:%s\n' % alignment.title)
                            blast_file.write('number of identities:%s\n' % hsp.identities)
                            blast_file.write('percentage identity:%s\n' % str(percentage_identity))
                            blast_file.write('percentage identity cutoff:%s\n' % str(percentage_identity_cutoff))
                            blast_file.write('number of positive:%s\n' % hsp.positives)
                            blast_file.write('number of gaps:%s\n' % hsp.gaps)
                            blast_protein=alignment.title.split("|")[3]
                            blast_file.write('protein:%s\n' % blast_protein)
                            #pdbl = PDBList()
                            blast_file.write('length:%s\n' % alignment.length)
                            blast_file.write('e value:%s\n' % hsp.expect)
                            blast_file.write(hsp.query[0:75] + '...\n')
                            blast_file.write(hsp.match[0:75] + '...\n')
                            blast_file.write(hsp.sbjct[0:75] + '...\n\n')
                            # write ligand data of pdb into file
                            blast_file.write("->Ligand information\n")
                            pdb_url_ligand_string = "http://www.rcsb.org//pdb/rest/ligandInfo?structureId=" + \
                                                 blast_protein
                            pdb_response = urllib.request.urlopen(pdb_url_ligand_string)
                            # print("%s\n" % pdbResponse.read().decode("utf-8"))
                            # print out ligand information
                            ligand_XML = pdb_response.read().decode("utf-8")
                            root = xml.etree.ElementTree.fromstring(ligand_XML)
                            blast_file.write("Protein Id:%s\n" % root.get("id"))
                            for ligand_info in root.findall('ligandInfo'):
                                for ligand in ligand_info.findall('ligand'):
                                    blast_file.write("Chemical id:%s\n" % ligand.get("chemicalID"))
                                    blast_file.write("type:%s\n" % ligand.get("type"))
                                    blast_file.write("molecular weight:%s\n" % ligand.get("molecularWeight"))
                                    blast_file.write("Chemical name:%s\n" % ligand.find("chemicalName").text)
                                    blast_file.write("Formula:%s\n" % ligand.find("formula").text)
                                    ligand_inchikey=ligand.find("InChIKey").text
                                    blast_file.write("InChIKey:%s\n" % ligand_inchikey)
                                    blast_file.write("InChI:%s\n" % ligand.find("InChI").text)
                                    blast_file.write("smiles:%s\n\n" % ligand.find("smiles").text)
                            # write pdb experimental details to file
                            blast_file.write("->Protein information\n")
                            pdb_url_description_string = "http://www.rcsb.org/pdb/rest/describePDB?structureId=" + \
                                                    blast_protein
                            pdb_response = urllib.request.urlopen(pdb_url_description_string)
                            description_XML= pdb_response.read().decode("utf-8")
                            root = xml.etree.ElementTree.fromstring(description_XML)
                            for PDB_info in root.findall("PDB"):
                                blast_file.write("Protein Id:%s\n" % PDB_info.get("structureId"))
                                blast_file.write("Experimental method:%s\n" % PDB_info.get("expMethod"))
                                pdb_resolution=PDB_info.get("resolution")
                                blast_file.write("Resolution:%s\n" % pdb_resolution)
                                blast_file.write("Number of entities:%s\n" % PDB_info.get("nr_entities"))
                                blast_file.write("Number of residues:%s\n" % PDB_info.get("nr_residues"))
                                blast_file.write("Number of atoms:%s\n\n" % PDB_info.get("nr_atoms"))
                            # populate blast_list
                            blast_list.append({definitions.DICT_BLAST_PROTEIN:blast_protein,
                            definitions.DICT_BLAST_PERCENTAGE_IDENTITY:round(percentage_identity,2),
                            definitions.DICT_BLAST_PROTEIN_RESOLUTION:float(pdb_resolution),
                            definitions.DICT_BLAST_INCHIKEY_FOUND:(ligand_inchikey==inchi_key)})
            blast_file.close()
            result_handle.close()
            return blast_list
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise



