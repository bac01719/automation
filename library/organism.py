# bola coker
# 19th May 2017
#
# class for organism

import library.definitions
from Bio import Entrez

class Organism:
    """methods related to the organism."""

    """static method to check if organism exist"""
    @staticmethod
    def check_exist(organism_name):
        try:
            Entrez.email=library.definitions.EMAIL_ADDRESS
            handle=Entrez.esearch(db=library.definitions.ORGANISM_DATABASE, term=organism_name)
            record = Entrez.read(handle)
            if len(record["IdList"])>0:
                exist=True
            else:
                exist=False
            return exist
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise
	
    
