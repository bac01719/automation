# bola coker
# 17th May 2017
#
# class for protein
import sys
from Bio import Entrez
import library.definitions
import xml.etree.ElementTree

class Protein:
    """class for protein"""

    """static method to return list of proteins for a given gene list"""
    @staticmethod
    def get_proteins(gene_list):
        try:
            for id_number in gene_list:
                # get full record detail in xml format for each id
                id_handle = Entrez.efetch(db=library.definitions.GENE_DATABASE, id=id_number, retmode="xml")
                id_XML=id_handle.read()
                # parse xml
                root = xml.etree.ElementTree.fromstring(id_XML)
                # check if xml is for discontinued protein
                discontinued=False
                # obtain gene track status of protein
                for entrezgene in root.findall('Entrezgene'):
                    for entrezgene_track_info in entrezgene.findall("Entrezgene_track-info"):
                        for gene_track in entrezgene_track_info.find("Gene-track"):
                            if gene_track.tag == "Gene-track_status" and gene_track.get("value") == "discontinued":
                                discontinued = True
                # obtain protein name, NCBI number, and if discontinued. Add to list
                proteins=[]
                for entrezgene in root.findall('Entrezgene'):
                    for entrezgene_locus in entrezgene.findall('Entrezgene_locus'):
                        for gene_commentary_first in entrezgene_locus.findall('Gene-commentary'):
                            for gene_commentary_products in gene_commentary_first.findall('Gene-commentary_products'):
                                for gene_commentary_second in gene_commentary_products.findall('Gene-commentary'):
                                    if gene_commentary_second[0].get('value') == 'peptide':
                                       proteins.append({library.definitions.DICT_PROTEIN_NAME:gene_commentary_second[1].text,
                                                        library.definitions.DICT_PROTEIN_NCBI_ID:gene_commentary_second[2].text,
                                                        library.definitions.DICT_GENE_DISCONTINUED: discontinued})
                # return proteins list
                return proteins
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

