# bola coker
# 21st May 2017
# main class to automate analysis of mutations in proteins for a given organism
#
import library.definitions

from library.gene import Gene
gene_list=Gene.get_genes("pncA","mycobacterium tuberculosis")
#gene_list=Gene.get_genes("gid","mycobacterium tuberculosis")
from library.protein import Protein
protein_list=Protein.get_proteins(gene_list)
from library.fasta import Fasta
protein_sequence_file_path=Fasta.get_fasta(protein_list,'/home/bola/Documents/Private/BBK/project/other')
from library.blast import Blast
blast_list=Blast.do_blast(protein_sequence_file_path,90,"IPEHBUMCGVEMRF-UHFFFAOYSA-N")
print(blast_list)
if len(blast_list)==0:
    print("********Homology model next**************")
elif len(blast_list)>0:
    for blast in blast_list:
        if blast[library.definitions.DICT_BLAST_INCHIKEY_FOUND]==True:
            print("*********** direct analysis next*************")
        elif blast[library.definitions.DICT_BLAST_INCHIKEY_FOUND]==False:
            print("*********** docking next **************")


