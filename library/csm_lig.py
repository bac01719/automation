# bola coker
# 9th July, 2017
# class to wrap csm-lig web service
#

import library.definitions as definitions
import urllib.request
import urllib.parse
import io
from library.multi_part_form import MultiPartForm
from bs4 import BeautifulSoup  # beautiful soap 4 installed with sudo python3 -m pip install beautifulsoup4
import os


class CSMLig:
    """ class to wrap csm-lig web service """

    _complex_pdb_file_path=""
    _complex_pdb_file_contents=None
    _complex_pdb_file_name=""
    _ligand_pdb_symbol=""
    _ligand_smiles=""

    def __init__(self,complex_pdb_file_path,ligand_pdb_symbol,ligand_smiles):
        global _complex_pdb_file_path  
        global _complex_pdb_file_contents
        global _complex_pdb_file_name
        global _ligand_smiles
        global _ligand_pdb_symbol
        global _binding_affinity_value
        _complex_pdb_file_path=complex_pdb_file_path
        _ligand_smiles=ligand_smiles
        _ligand_pdb_symbol=ligand_pdb_symbol
        print("\nObtaining CSM-Lig predicted binding affinity\n")
        try:
            # read in protein contents
            with open(complex_pdb_file_path, 'r') as complex_pdb_file:
                _complex_pdb_file_contents = complex_pdb_file.read()
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    """ method to web scrape csm_website """
    def get_binding_affinity(self):
        global _complex_pdb_file_path
        global _complex_pdb_file_contents
        global _complex_pdb_file_name
        global _ligand_smiles
        global _ligand_pdb_symbol
        global _binding_affinity_value
        # Create the form with simple fields
        form = MultiPartForm()
        form.add_file('pdb_file', _complex_pdb_file_path, fileHandle=io.BytesIO(_complex_pdb_file_contents.encode('utf-8')))
        form.add_field('lig_id', _ligand_pdb_symbol)
        form.add_field('smiles_str', _ligand_smiles)
        form.add_field('pred_type', 'single')
        form.add_field('pdb_zip', '')
        # Build the request, including the byte-string for the data posted
        data_bytes = bytes(form)
        try:
            mcsm_lig_results = urllib.request.Request(definitions.CSM_LIG_URL, data_bytes)
            mcsm_lig_results.add_header('Content-type', form.get_content_type())
            mcsm_lig_results.add_header('Content-length', len(data_bytes))
            mcsm_lig_response = urllib.request.urlopen(mcsm_lig_results)
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise
        try:
            # obtain predicted binding affinity
            binding_affinity_value="NOT DEFINED"
            beautiful_soup = BeautifulSoup(mcsm_lig_response, 'html.parser')
            #print(beautiful_soup.prettify())
            count = -1
            font_anchors = beautiful_soup.find_all('font')
            for font_anchor in font_anchors:
                count += 1
                if count == 0:
                    body_position = font_anchor
                    binding_affinity_value = body_position.next_element.string
        except Exception:
            print("Error with Protein: %s \nLigand Symbol: %s\nLigand SMILES\n" % \
                  (_complex_pdb_file_name,_ligand_pdb_symbol,_ligand_smiles))
        finally:
            return binding_affinity_value
