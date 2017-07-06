# bola coker
# 21st May 2017
# main class to automate analysis of mutations in proteins for a given organism
#
import library.definitions as definitions
from library.gene import Gene
from library.protein import Protein
from library.fasta import Fasta
from library.blast import Blast
from library.clustalomega import ClustalOmega
from library.homology import Homology
from library.openbabel import OpenBabel
from library.pdb2pqr import PDB2PQR
from library.prepare_receptor import PrepareReceptor
from library.vina import Vina
import time

# define global paths
data_analysis_complex_file_path=""
docking_receptor_file_path=""
# get input - gene, organism, SMILES, analysis_folder
gene_list=Gene.get_genes("pncA","mycobacterium tuberculosis")
ligand_smiles="NC(=O)c1cnccn1"
home_folder="/home/bola/Documents/Private/BBK/project/other"
# print disclaimer
print("Pipeline to aid in research of effects of mutation on protein-ligand binding affinity")
print("Developed by Bola Coker")
print("References:")
print("RCSB website")
print("OPAL website - http://webservices.rbvi.ucsf.edu/opal2/dashboard?command=serviceList")
print("OPAL website - http://nbcr-222.ucsd.edu/opal2/dashboard?command=serviceList\n")
# print start time
print("Start : %s " % time.strftime("%Y-%m-%d %H:%M:%S"))
#gene_list=Gene.get_genes("gid","mycobacterium tuberculosis")
# get protein of gene
protein_list=Protein.get_proteins(gene_list)
# Get fasta of protein.
fasta_list=Fasta.get_fasta(protein_list,home_folder)
for fasta_item in fasta_list[definitions.DICT_FASTA_LIST]:
    print("processing %s\n" % fasta_item)
    # blast fast to get templates for homologue processing also check if pdb found includes SMILES string of ligand
    blast_list=Blast.do_blast(fasta_item[definitions.DICT_FASTA_FILE_PATH],\
                              definitions.BLAST_PERCENTAGE_IDENTITY_CUTOFF,ligand_smiles)
    print("blast results")
    print("lower percentage identity for homology %s" % str(definitions.HOMOLOGY_PERCENT_IDENTITY_LOWER_RANGE))
    print("upper percentage identity for homology %s" % str(definitions.HOMOLOGY_PERCENT_IDENTITY_UPPER_RANGE))
    print("minimum resoultion for homology %s" % str(definitions.HOMOLOGY_MINIMUM_RESOLUTION))
    print("minimum percentage identity to bypass homology %s" % str(definitions.MINIMUM_PERCENT_IDENTITY_FOR_OTHER))
    print("minimum resolution to bypass homology %s\n" % str(definitions.OTHER_MINIMUM_RESOLUTION))
    if len(blast_list)>0:
        template_protein_names=[]
        for blast in blast_list:
            if (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]>=definitions.HOMOLOGY_PERCENT_IDENTITY_LOWER_RANGE \
                    and blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]<=\
                        definitions.HOMOLOGY_PERCENT_IDENTITY_UPPER_RANGE) \
                    and blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]<=definitions.HOMOLOGY_MINIMUM_RESOLUTION:
                template_protein_names.append(blast[definitions.DICT_BLAST_PROTEIN])
                print("%s :" % blast[definitions.DICT_BLAST_PROTEIN], end="")
                print("********Homology modelling next (%s%%,%s,%s) **************" % \
                      (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],
                       blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION],
                       blast[definitions.DICT_BLAST_SMILES_FOUND]))
            elif blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]>=definitions.MINIMUM_PERCENT_IDENTITY_FOR_OTHER and \
                    blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]<=definitions.OTHER_MINIMUM_RESOLUTION:
                print("%s :" % blast[definitions.DICT_BLAST_PROTEIN],end="")
                if blast[definitions.DICT_BLAST_SMILES_FOUND]==True:
                   # TODO download blast protein and set data_analysis_complex_file_path
                   print("*********** data analysis next (%s%%,%s,%s) *************" % \
                         (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],
                          blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION],
                          blast[definitions.DICT_BLAST_SMILES_FOUND]))
                elif blast[definitions.DICT_BLAST_SMILES_FOUND]==False:
                   # TODO download blast protein and set docking_receptor_file_path
                   print("*********** docking next (%s%%,%s,%s) ****************" % \
                         (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],
                          blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION],
                          blast[definitions.DICT_BLAST_SMILES_FOUND]))
            else:
                print("%s :" % blast[definitions.DICT_BLAST_PROTEIN], end="")
                print("*********** low resolution (%s%%,%s,%s) ****************" % \
                      (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],
                       blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION],
                       blast[definitions.DICT_BLAST_SMILES_FOUND]))
    else:
        print("*********** blast did not find sequences with percentage identities above %s%% ****************" % \
              definitions.BLAST_PERCENTAGE_IDENTITY_CUTOFF)
    # homology processing
    if len(template_protein_names)>0:
        # save templates into fasta file path
        template_file_path=Fasta.get_templates(template_protein_names,fasta_item)
        # do multiple sequence alignment and save file
        clustal_omega_object=ClustalOmega(fasta_item[definitions.DICT_FASTA_FILE_PATH],template_file_path)
        # if no error store msa file
        if clustal_omega_object.get_status()==definitions.OPAL_SUCCESS:
            # save msa file
            msa_file_path=clustal_omega_object.save_msa()
            # do homology
            homology_object=Homology(msa_file_path)
            if homology_object.get_status()==definitions.OPAL_SUCCESS:
                homology_save_output=homology_object.save_output()
                print(homology_save_output)
                docking_receptor_file_path=homology_save_output[definitions.DICT_HOMOLOGY_BEST_PATH]
            else:
                # an error occurred
                print("Error with Modeller Opal web service")
                print(homology_object.get_error())
        else:
            # an error occurred
            print("Error with Clustal Omega Opal web service")
            print(clustal_omega_object.get_error())
    # docking processing
    if docking_receptor_file_path!="":
        # create pdbqt file ligand with given SMILEs string
        docking_folder=fasta_list[definitions.DICT_FASTA_PROTEIN_FOLDER]+definitions.FILE_SEPARATOR+"docking"
        open_babel_object=OpenBabel(ligand_smiles,docking_folder)
        if open_babel_object.get_status()==definitions.OPAL_SUCCESS:
            ligand_pdbqt_file_path=open_babel_object.save_ligand()
            # prepare protein
            # first protonate protein with pqr2pdb server
            pdb2pqr_object=PDB2PQR(docking_receptor_file_path,docking_folder)
            if pdb2pqr_object.get_status()==definitions.OPAL_SUCCESS:
                receptor_pqr_file_path=pdb2pqr_object.save_output()
                # second prepare receptor and convert pqr file to pdbqt autodock file
                prepare_receptor_object=PrepareReceptor(receptor_pqr_file_path,docking_folder)
                if prepare_receptor_object.get_status()==definitions.OPAL_SUCCESS:
                    receptor_pdbqt_file_path=prepare_receptor_object.save_output()
                    # third do docking with autodock vina opal web service
                    vina_object=Vina(receptor_pdbqt_file_path,ligand_pdbqt_file_path,docking_folder)
                    if vina_object.get_status()==definitions.OPAL_SUCCESS:
                        data_analysis_complex_file_path=vina_object.save_output()
                    else:
                        # an error occurred
                        print("Error with vina Opal web service")
                        print(vina_object.get_error())
                else:
                    # an error occurred
                    print("Error with prepare receptor Opal web service")
                    print(prepare_receptor_object.get_error())
            else:
                # an error occurred
                print("Error with PDB2PQR Opal web service")
                print(pdb2pqr_object.get_error())
        else:
            # an error occurred
            print("Error with Open babel SMILES converison Opal web service")
            print(open_babel_object.get_error())

# print end
print("End  : %s " % time.strftime("%Y-%m-%d %H:%M:%S"))



