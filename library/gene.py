# bola coker
# 17th May 2017
#
# class for gene

from Bio import Entrez
import library.definitions

class Gene:
    """Class for gene"""    

    """static method to get NCBI Ids for genes"""
    @staticmethod
    def get_genes(gene_name,organism_name):
        try:
            Entrez.email=library.definitions.EMAIL_ADDRESS
            search_string="("+gene_name+"[Gene Name]) AND "+organism_name+"[Organism]"
            handle = Entrez.esearch(db=library.definitions.GENE_DATABASE, term=search_string)
            record = Entrez.read(handle)
            handle.close()
            return record["IdList"]
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
       	    raise

   
