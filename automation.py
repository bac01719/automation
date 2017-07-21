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
from library.duet import Duet
from library.mcsm_lig import MCSMLig
from library.replace_bfactor import ReplaceBFactor
from library.plip import PLIP
import library.utilities as utilities
import os
import time
from optparse import OptionParser
import sys
import re

__version__="automation 0.1 (c) 2017"
# set up command line parser
parser=OptionParser(version=__version__)
parser.add_option("-o","--organism",type="string",dest="organism",\
                  help="organism as defined in Entrex ((http://www.ncbi.nlm.nih.gov/Entrez) taxonomy database")
parser.add_option("-g","--gene",type="string",dest="gene",help="gene of organism as defined in Entrez gene database")
parser.add_option("-s","--SMILES",type="string",dest="smiles",\
                  help="SMILES string for ligand as defined in ZINC (http://zinc.docking.org) database")
parser.add_option("-f","--folder",type="string",dest="home_folder",\
                  help="full path to folder to store analysis results")
parser.add_option("-m","--mutation_file_path",type="string",dest="mutation_file_path",\
                  help="""Full path to SNP file required for DUET and mCSM-Lig analysis. 
                  CSV file with no header row in format of [chain,mutation] 
                  e.g, [A,V183L] or [K,Y99*]. Required for data analysis")""", default="")
parser.add_option("-v","--snp",type="string",dest="snp",\
                  help="single point mutation for structural analysis in form of [wild type,position,snp],e.g., V183L.",\
                  default="")
(options,args)=parser.parse_args()

# check parameters
if not options.organism:
    print("Organism required. For help use --help option")
    sys.exit()
elif not utilities.verify_entrez_option(definitions.ORGANISM_DATABASE,options.organism):
    print("Organism not found in Entrez taxonomy database. For help use --help option")
    sys.exit()
elif not options.gene:
    print("Gene required. For help use --help option")
    sys.exit()
elif not utilities.verify_entrez_option(definitions.GENE_DATABASE,options.gene):
    print("Gene not found i.mutation_filen Entrez gene database. For help use --help option")
    sys.exit()
elif not options.smiles:
    print("SMILES string for ligand required. For help use --help option")
    sys.exit()
elif not utilities.verify_zinc_smiles(options.smiles):
    print("SMILES string not found in ZINC database. For help use --help option")
    sys.exit()
elif not options.home_folder:
    print("full path to folder to store analysis results required. For help use --help option")
    sys.exit()
elif not os.path.exists(options.home_folder):
    print("full path to folder to store analysis results not found. For help use --help option")
    sys.exit()
elif options.mutation_file_path:
    if not os.path.exists(options.mutation_file_path):
        print("full path to folder to store analysis results not found. For help use --help option")
        sys.exit()
    elif not utilities.verify_mutation_file(options.mutation_file_path):
        print("error in mutation file. For help use --help option")
        sys.exit()
elif options.snp:
    parse=re.compile(definitions.SNP_REG_EXP)
    if parse.match(options.snp)==None:
        print("error in snp format. For help use --help option")
        sys.exit()

# capture arguments and set program global variables
gene_list=Gene.get_genes(options.gene,options.organism)
ligand_smiles=options.smiles
home_folder=options.home_folder
mutation_file_path=options.mutation_file_path
snp=options.snp

# define global paths
template_protein_names=[]
data_analysis_complex_names=[]
docking_protein_names=[]
# print disclaimer
print(__version__)
print("Pipeline to aid in research of effects of mutation on protein-ligand binding affinity")
print("Developed by Bola Coker\n\n")
print("References:")
print("Zinc website - %s" % definitions.ZINC_URL)
print("RCSB website - https://www.rcsb.org")
print("OPAL website - http://webservices.rbvi.ucsf.edu/opal2/dashboard?command=serviceList")
print("PLIP website - https://projects.biotec.tu-dresden.de/plip-web/plip/")
print("OPAL website - http://nbcr-222.ucsd.edu/opal2/dashboard?command=serviceList\n")

# print start time
print("Start : %s " % time.strftime("%Y-%m-%d %H:%M:%S"))
# print inputs
print("Command line arguments")
print("Gene %s : " % options.gene)
print("Organism %s : " % options.organism)
print("Ligand SMILES %s : " % ligand_smiles)
print("Home folder %s : " % home_folder)
print("mutation_file_path %s : " % mutation_file_path)
print("snp %s : \n" % snp)
# get protein of gene
protein_list=Protein.get_proteins(gene_list)
# Get fasta of protein.
fasta_list=Fasta.get_fasta(protein_list,home_folder,snp)
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
        for blast in blast_list:
            if (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]>=definitions.HOMOLOGY_PERCENT_IDENTITY_LOWER_RANGE \
                    and blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]<=\
                        definitions.HOMOLOGY_PERCENT_IDENTITY_UPPER_RANGE) \
                    and blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]<=definitions.HOMOLOGY_MINIMUM_RESOLUTION:
                template_protein_names.append({ definitions.DICT_BLAST_PROTEIN: blast[definitions.DICT_BLAST_PROTEIN],\
                                               definitions.DICT_BLAST_PERCENTAGE_IDENTITY :\
                                                   blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]})
                print("%s :" % blast[definitions.DICT_BLAST_PROTEIN], end="")
                print("********Homology modelling next (%s%%,%s,%s) **************" % \
                      (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],
                       blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION],
                       blast[definitions.DICT_BLAST_SMILES_FOUND]))
            elif blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]>=definitions.MINIMUM_PERCENT_IDENTITY_FOR_OTHER and \
                    blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]<=definitions.OTHER_MINIMUM_RESOLUTION:
                print("%s :" % blast[definitions.DICT_BLAST_PROTEIN],end="")
                if blast[definitions.DICT_BLAST_SMILES_FOUND]==True:
                   protein_file_path = utilities.download_protein(blast[definitions.DICT_BLAST_PROTEIN], \
                                                                   fasta_list[definitions.DICT_FASTA_PROTEIN_FOLDER])
                   data_analysis_complex_names.append(protein_file_path)
                   print("*********** data analysis next (%s%%,%s,%s) *************" % \
                         (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],
                          blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION],
                          blast[definitions.DICT_BLAST_SMILES_FOUND]))
                elif blast[definitions.DICT_BLAST_SMILES_FOUND]==False:
                   protein_file_path=utilities.download_protein(blast[definitions.DICT_BLAST_PROTEIN],\
                                                                fasta_list[definitions.DICT_FASTA_PROTEIN_FOLDER])
                   docking_protein_names.append(protein_file_path)
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
        # single template homology required
        # select protein with max alignment
        max_identity=0
        template_protein_names_2=template_protein_names
        if definitions.SINGLE_HOMOLOGY_TEMPLATE:    # otherwise use multiple segments for homology
            for template in template_protein_names:
                if template[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]>max_identity:
                    max_identity=template[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]
                    template_protein_names_2=[template[definitions.DICT_BLAST_PROTEIN]]
        # save templates into fasta file path
        template_file_path=Fasta.get_templates(template_protein_names_2,fasta_item)
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
                docking_protein_names.append(homology_save_output[definitions.DICT_HOMOLOGY_BEST_PATH])
            else:
                # an error occurred
                print("Error with Modeller Opal web service")
                print(homology_object.get_error())
        else:
            # an error occurred
            print("Error with Clustal Omega Opal web service")
            print(clustal_omega_object.get_error())
    # docking processing
    if len(docking_protein_names)>0:
        # create pdbqt file ligand with given SMILEs string
        docking_folder=fasta_list[definitions.DICT_FASTA_PROTEIN_FOLDER]+definitions.FILE_SEPARATOR+"docking"
        open_babel_object=OpenBabel(ligand_smiles,docking_folder)
        if open_babel_object.get_status()==definitions.OPAL_SUCCESS:
            ligand_pdbqt_file_path=open_babel_object.save_ligand()
            # prepare protein
            for docking_receptor_file_path in docking_protein_names:
                # first protonate protein with pqr2pdb server
                protein_pdb_name=os.path.basename(docking_receptor_file_path).replace(".pdb","")
                protein_docking_folder=docking_folder+definitions.FILE_SEPARATOR+protein_pdb_name
                pdb2pqr_object=PDB2PQR(docking_receptor_file_path,protein_docking_folder)
                if pdb2pqr_object.get_status()==definitions.OPAL_SUCCESS:
                    receptor_pqr_file_path=pdb2pqr_object.save_output()
                    # second prepare receptor and convert pqr file to pdbqt autodock file
                    prepare_receptor_object=PrepareReceptor(receptor_pqr_file_path,protein_docking_folder)
                    if prepare_receptor_object.get_status()==definitions.OPAL_SUCCESS:
                        receptor_pdbqt_file_path=prepare_receptor_object.save_output()
                        # third do docking with autodock vina opal web service
                        vina_object=Vina(receptor_pdbqt_file_path,ligand_pdbqt_file_path,protein_docking_folder)
                        if vina_object.get_status()==definitions.OPAL_SUCCESS:
                            data_analysis_complex_file_path=vina_object.save_output()
                            data_analysis_complex_names.append(data_analysis_complex_file_path)
                            # fourth produce PLIP file
                            plip_object=PLIP(data_analysis_complex_file_path)
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
    # data analysis processing
    if len(data_analysis_complex_names)>0 and mutation_file_path!="":
        data_analysis_folder = fasta_list[definitions.DICT_FASTA_PROTEIN_FOLDER] + definitions.FILE_SEPARATOR + \
                               "data_analysis"
        for data_analysis_complex_file_path in data_analysis_complex_names:
            complex_pdb_name=os.path.basename(data_analysis_complex_file_path).replace(".pdb","")
            complex_data_analysis_folder=data_analysis_folder+definitions.FILE_SEPARATOR+complex_pdb_name
            # perform duet analysis
            duet_analysis_object=Duet(complex_data_analysis_folder,data_analysis_complex_file_path)
            duet_csv_file_path=duet_analysis_object.get_file_path()
            # perform mcsm_lig analysis
            mcsm_lig_analysis_object=MCSMLig(complex_data_analysis_folder,data_analysis_complex_file_path,\
                                              mutation_file_path,"LIG",ligand_smiles)
            mcsm_lig_csv_file_path=mcsm_lig_analysis_object.get_file_path()
            # replace bfactor pdb field with values from mCSM,SDM,DUET and mCSM_Lig
            replace_bfactor_analysis_object=ReplaceBFactor(data_analysis_complex_file_path,duet_csv_file_path,\
                                                           mcsm_lig_csv_file_path)

# print end
print("End  : %s " % time.strftime("%Y-%m-%d %H:%M:%S"))



