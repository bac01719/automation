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

# define global paths
data_analysis_complex_file_path=""
docking_receptor_file_path=""
# get input - gene and organism
gene_list=Gene.get_genes("pncA","mycobacterium tuberculosis")
#gene_list=Gene.get_genes("gid","mycobacterium tuberculosis")
protein_list=Protein.get_proteins(gene_list)
fasta_list=Fasta.get_fasta(protein_list,'/home/bola/Documents/Private/BBK/project/other')
for fasta_item in fasta_list:
    print("processing %s\n\n" % fasta_item)
    blast_list=Blast.do_blast(fasta_item[definitions.DICT_FASTA_FILE_PATH],\
                              definitions.BLAST_PERCENTAGE_IDENTITY_CUTOFF,\
                              "IPEHBUMCGVEMRF-UHFFFAOYSA-N")
    if len(blast_list)>0:
        template_protein_names=[]
        for blast in blast_list:
            if (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]>=definitions.HOMOLOGY_PERCENT_IDENTITY_LOWER_RANGE \
                    and blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]<=\
                        definitions.HOMOLOGY_PERCENT_IDENTITY_UPPER_RANGE) \
                    and blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]<=definitions.HOMOLOGY_MINIMUM_RESOLUTION:
                template_protein_names.append(blast[definitions.DICT_BLAST_PROTEIN])
                print("%s :" % blast[definitions.DICT_BLAST_PROTEIN], end="")
                print("********Homology model next (%s%%,%s) **************" % \
                      (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],\
                       blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]))
            elif blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]>=definitions.MINIMUM_PERCENT_IDENTITY_FOR_OTHER and \
                    blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]<=definitions.OTHER_MINIMUM_RESOLUTION:
                print("%s :" % blast[definitions.DICT_BLAST_PROTEIN],end="")
                if blast[definitions.DICT_BLAST_INCHIKEY_FOUND]==True:
                   # TODO download blast protein and set data_analysis_complex_file_path
                   print("*********** data analysis next (%s%%,%s) *************" % \
                         (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],
                          blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]))
                elif blast[definitions.DICT_BLAST_INCHIKEY_FOUND]==False:
                   # TODO download blast protein and set docking_receptor_file_path
                   print("*********** docking next (%s%%,%s) ****************" % \
                         (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],
                          blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]))
            else:
                print("%s :" % blast[definitions.DICT_BLAST_PROTEIN], end="")
                print("*********** low resolution (%s%%,%s) ****************" % \
                      (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]))
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


