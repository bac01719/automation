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


gene_list=Gene.get_genes("pncA","mycobacterium tuberculosis")
#gene_list=Gene.get_genes("gid","mycobacterium tuberculosis")
protein_list=Protein.get_proteins(gene_list)
fasta_file_path_list=Fasta.get_fasta(protein_list,'/home/bola/Documents/Private/BBK/project/other')
for fasta_file_path in fasta_file_path_list:
    print("processing %s\n\n" % fasta_file_path)
    blast_list=Blast.do_blast(fasta_file_path,definitions.BLAST_PERCENTAGE_IDENTITY_CUTOFF,
                          "IPEHBUMCGVEMRF-UHFFFAOYSA-N")
    #print(blast_list)
    if len(blast_list)>0:
        template_protein_names=""
        for blast in blast_list:
            if (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]>=definitions.HOMOLOGY_PERCENT_IDENTITY_LOWER_RANGE and \
                    blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]<=definitions.HOMOLOGY_PERCENT_IDENTITY_UPPER_RANGE) \
                    and blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]<=definitions.HOMOLOGY_MINIMUM_RESOLUTION:
                template_protein_names=template_protein_names+blast[definitions.DICT_BLAST_PROTEIN]+","
                print("%s :" % blast[definitions.DICT_BLAST_PROTEIN], end="")
                print("********Homology model next (%s%%,%s) **************" % \
                      (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]))
            elif blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY]>=definitions.MINIMUM_PERCENT_IDENTITY_FOR_OTHER and \
                    blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]<=definitions.OTHER_MINIMUM_RESOLUTION:
                print("%s :" % blast[definitions.DICT_BLAST_PROTEIN],end="")
                if blast[definitions.DICT_BLAST_INCHIKEY_FOUND]==True:
                   print("*********** direct analysis next (%s%%,%s) *************" % \
                         (blast[definitions.DICT_BLAST_PERCENTAGE_IDENTITY],
                          blast[definitions.DICT_BLAST_PROTEIN_RESOLUTION]))
                elif blast[definitions.DICT_BLAST_INCHIKEY_FOUND]==False:
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
    # save templates into fasta file
    template_file_path=Fasta.get_templates(template_protein_names.strip(","),fasta_file_path)
    # do multiple sequence alignment and save file
    clustal_omega_object=ClustalOmega(template_file_path)
    # if no error store msa file
    if clustal_omega_object.get_status()==definitions.OPAL_SUCCESS:
        # save msa file
        clustal_omega_object.save_msa()
    else:
        # an error occurred
        print(clustal_omega_object.get_error())


