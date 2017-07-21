# bola coker
# 7th July, 2017
# class to wrap duet web service
#

import library.definitions as definitions
import urllib.request
import urllib.parse
import io
from library.multi_part_form import MultiPartForm
from bs4 import BeautifulSoup  # beautiful soap 4 installed with sudo python3 -m pip install beautifulsoup4
import os


class Duet:
    """ class to wrap duet service """

    _complex_pdb_file_path=""
    _complex_pdb_file_contents=None
    _complex_pdb_file_name=""
    _duet_results_file_path=""

    def __init__(self,data_analysis_folder,complex_pdb_file_path,mutation_file_path):
        global _complex_pdb_file_path  
        global _complex_pdb_file_contents
        global _complex_pdb_file_name
        _complex_pdb_file_path=complex_pdb_file_path
        print("\nMutation analysis using DUET website\n")
        try:
            # read in protein contents
            with open(complex_pdb_file_path, 'r') as complex_pdb_file:
                _complex_pdb_file_contents = complex_pdb_file.read()
            # check data analysis folder exist otherwise create
            if not os.path.exists(data_analysis_folder):
                os.mkdir(data_analysis_folder)
            # create DUET csv file
            duet_results_file_path = data_analysis_folder+definitions.FILE_SEPARATOR+"duet_analysis.csv"
            duet_results_file = open(duet_results_file_path, 'w')
            duet_results_file.write(
                'PDB,' + 'mCSM,' + 'SDM,' + 'DUET,' + 'WildType,' + 'Position,' + 'MutantType,' + 'Chain,' + \
                'RelativeSolventAccessibility,' + 'SecondaryStructure,' + 'SideChainHydrogenBond,' + '\n')
            with open(mutation_file_path, 'r') as mutations:
                for line in mutations:
                    mutation_line = line.strip('\n')
                    (chain, mutation) = (mutation_line.split(',')[0], mutation_line.split(',')[1])
                    _complex_pdb_file_name = os.path.basename(complex_pdb_file_path)
                    print('analysing protein:%s mutation:%s chain:%s' % (_complex_pdb_file_name, mutation, chain))
                    duet_results = self.get_duet_results(mutation, chain)
                    # write results to file
                    duet_results_file.write(
                        duet_results[0] + ',' + duet_results[1] + ',' + duet_results[2] + ',' + duet_results[3] +\
                        ',' + duet_results[4] + ',' + duet_results[5] + ',' + duet_results[6] + ',' + \
                        duet_results[7] + ',' + duet_results[8] + ',' + \
                        duet_results[9] + ',' + duet_results[10] + '\n')
            duet_results_file.close()
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    """ method to web scrape DUET website """
    def get_duet_results(self,mutation, chain):
        global _complex_pdb_file_path
        global _complex_pdb_file_contents
        global _complex_pdb_file_name
        # Create the form with simple fields
        form = MultiPartForm()
        form.add_file('wild', _complex_pdb_file_path, fileHandle=io.BytesIO(_complex_pdb_file_contents.encode('utf-8')))
        form.add_field('pdb_code', '')
        form.add_field('mutation', mutation)
        form.add_field('chain', chain)
        form.add_field('run', 'single')
        form.add_field('mutation_sys', '')
        form.add_field('chain_sys', '')
        # Build the request, including the byte-string for the data posted
        data_bytes = bytes(form)
        try:
            duet_results = urllib.request.Request(definitions.DUET_URL, data_bytes)
            duet_results.add_header('Content-type', form.get_content_type())
            duet_results.add_header('Content-length', len(data_bytes))
            duet_response = urllib.request.urlopen(duet_results)
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise
        # Obtain response
        mcsm_value = definitions.MISSING_ANALYSIS_VALUE
        sdm_value = definitions.MISSING_ANALYSIS_VALUE
        duet_value = definitions.MISSING_ANALYSIS_VALUE
        wild_type_value = definitions.MISSING_ANALYSIS_VALUE
        position_value = definitions.MISSING_ANALYSIS_VALUE
        mutant_type_value = definitions.MISSING_ANALYSIS_VALUE
        chain_value = definitions.MISSING_ANALYSIS_VALUE
        relative_solvent_accessibility_value = definitions.MISSING_ANALYSIS_VALUE
        secondary_structure_value = definitions.MISSING_ANALYSIS_VALUE
        side_chain_hydrogen_bond_value = definitions.MISSING_ANALYSIS_VALUE
        try:
            beautiful_soup = BeautifulSoup(duet_response, 'html.parser')
            # print(beautiful_soup.prettify())
            count = -1
            font_anchors = beautiful_soup.body.find_all('font')
            for font_anchor in font_anchors:
                count += 1
                if count == 0:
                    body_position = font_anchor
                    mcsm_value = body_position.next_element.string
                    body_position = body_position.next_element
                    mcsm_value = mcsm_value + body_position.next_element.string
                    body_position = body_position.next_element
                    mcsm_value = mcsm_value + body_position.next_element.next_element.string
                elif count == 1:
                    body_position = font_anchor
                    sdm_value = body_position.next_element.string
                    body_position = body_position.next_element
                    sdm_value = sdm_value + body_position.next_element.string
                    body_position = body_position.next_element
                    sdm_value = sdm_value + body_position.next_element.next_element.string
                elif count == 2:
                    body_position = font_anchor.find_next('font')
                    duet_value = body_position.next_element.string
                    body_position = body_position.next_element
                    duet_value = duet_value + body_position.next_element.string
                    body_position = body_position.next_element
                    duet_value = duet_value + body_position.next_element.next_element.string
                elif count == 3:
                    body_position = font_anchor.find_next('font').find_next('b')
                    wild_type_value = body_position.string
                    body_position = body_position.find_next('b')
                    position_value = body_position.string
                    body_position = body_position.find_next('b')
                    mutant_type_value = body_position.string
                    body_position = body_position.find_next('b')
                    chain_value = body_position.string
                    body_position = body_position.find_next('b')
                    relative_solvent_accessibility_value = body_position.string
                    body_position = body_position.find_next('b')
                    secondary_structure_value = body_position.string
                    body_position = body_position.find_next('b')
                    side_chain_hydrogen_bond_value = body_position.string
        except Exception:
            print("Error with mutation: %s chain: %s" % (mutation, chain))
        finally:
            return (_complex_pdb_file_name, mcsm_value, sdm_value, duet_value, wild_type_value, position_value,\
                    mutant_type_value, chain_value, relative_solvent_accessibility_value, secondary_structure_value,\
                    side_chain_hydrogen_bond_value)

    def get_file_path(self):
        global _duet_results_file_path
        return _duet_results_file_path

