# bola coker
# 10th July, 2017
# class to wrap mcsm_lig web service
#

import library.definitions as definitions
import library.utilities as utilities
import urllib.request
import urllib.parse
import io
from library.multi_part_form import MultiPartForm
from bs4 import BeautifulSoup  # beautiful soap 4 installed with sudo python3 -m pip install beautifulsoup4
import os
from library.csm_lig import CSMLig

class MCSMLig:
    """ class to wrap mcsm_lig service """

    _complex_pdb_file_path=""
    _complex_pdb_file_contents=None
    _complex_pdb_file_name=""
    _lig_pred_affinity=""
    _mcsm_lig_results_file_path=""

    def __init__(self,data_analysis_folder,complex_pdb_file_path,mutation_file_path,ligand_pdb_symbol,ligand_smiles):
        global _complex_pdb_file_path  
        global _complex_pdb_file_contents
        global _complex_pdb_file_name
        global _wild_type_affinity
        global _mcsm_lig_results
        global _lig_pred_affinity
        global _mcsm_lig_results_file_path
        _complex_pdb_file_path=complex_pdb_file_path
        _complex_pdb_file_name = os.path.basename(complex_pdb_file_path)
        print("\nMutation analysis using mCSM_lig website\n")
        try:
            # read in protein contents
            with open(complex_pdb_file_path, 'r') as complex_pdb_file:
                _complex_pdb_file_contents = complex_pdb_file.read()
            # check data analysis folder parent folder does not exist otherwise create
            if not os.path.exists(os.path.dirname(data_analysis_folder)):
                os.mkdir(os.path.dirname(data_analysis_folder))
            # check data analysis folder exist otherwise create
            if not os.path.exists(data_analysis_folder):
                os.mkdir(data_analysis_folder)
            # read in wild type affinity
            if definitions.AFFINITY_USE_CSMLIG==True:
                csm_lig_object=CSMLig(complex_pdb_file_path,ligand_pdb_symbol,ligand_smiles)
                _lig_pred_affinity=csm_lig_object.get_binding_affinity()
                # compute wild type as negative antilog of predicted affinity value on CSM-Lig website
                # convert to nano seconds
                #wild_type_affinity=round(pow(10,-float(_csm_lig_pred_affinity))*pow(10,9),3)
                wild_type_affinity=utilities.compute_wildtype_nanomoles(float(_lig_pred_affinity))
                affinity_column_header='CSM_Lig Affinity -log10(Kd|Ki),'
            elif definitions.AFFINITY_USE_CSMLIG==False:
                #use autodock vina free energy value instead
                vina_out_file_path=os.path.dirname(complex_pdb_file_path)+definitions.FILE_SEPARATOR+"vina.out"
                _lig_pred_affinity=utilities.get_affinity_calpermol(vina_out_file_path)
                wild_type_affinity=round(utilities.compute_binding_affinity_nanomoles(_lig_pred_affinity),5)
                affinity_column_header = 'Vina Affinity cal/mol,'
            # create mCSM_LIG csv file
            _mcsm_lig_results_file_path = data_analysis_folder+definitions.FILE_SEPARATOR+"mcsm_lig_analysis.csv"
            mcsm_lig_results_file = open(_mcsm_lig_results_file_path, 'w')
            mcsm_lig_results_file.write('PDB,'+affinity_column_header+'WildTypeAffinity nM,'+\
                                        'PredictedAffinityChange,'+'WildType,'+'Position,'+'MutantType,'+'Chain,' +\
                                        'LigandID,'+'DistanceToLigand,'+'DUETStabilityChange,'+'\n')
            with open(mutation_file_path, 'r') as mutations:
                for line in mutations:
                    mutation_line = line.strip('\n')
                    (chain, mutation) = (mutation_line.split(',')[0], mutation_line.split(',')[1])
                    print('analysing protein:%s mutation:%s chain:%s ligand_id:%s affinity:%s' %\
                          (_complex_pdb_file_name, mutation, chain, ligand_pdb_symbol, wild_type_affinity))
                    mcsm_lig_results = self.get_mcsm_lig_results(mutation, chain, ligand_pdb_symbol, wild_type_affinity)
                    # write results to file
                    mcsm_lig_results_file.write(
                        mcsm_lig_results[0] + ',' + str(mcsm_lig_results[1]) + ',' + str(mcsm_lig_results[2]) + ',' + \
                        mcsm_lig_results[3] + ',' + mcsm_lig_results[4] + ',' + mcsm_lig_results[5] + ',' + \
                        mcsm_lig_results[6] + ',' + mcsm_lig_results[7] + ',' + mcsm_lig_results[8] + ',' + \
                        mcsm_lig_results[9] + ',' + mcsm_lig_results[10] + '\n')
            mcsm_lig_results_file.close()
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    """ method to web scrape mCSM_Lig website """
    def get_mcsm_lig_results(self,mutation, chain, lig_id, wild_type_affinity):
        global _complex_pdb_file_path
        global _complex_pdb_file_contents
        global _complex_pdb_file_name
        global _lig_pred_affinity
        # Create the form with simple fields
        form = MultiPartForm()
        form.add_file('wild', _complex_pdb_file_name, fileHandle=io.BytesIO(_complex_pdb_file_contents.encode('utf-8')))
        form.add_field('pdb_code', '')
        form.add_field('mutation', mutation)
        form.add_field('chain', chain)
        form.add_field('lig_id', lig_id)
        form.add_field('affin_wt', str(wild_type_affinity))
        form.add_field('run', 'single')
        # Build the request, including the byte-string for the data posted
        data_bytes = bytes(form)
        try:
            mcsm_lig_results = urllib.request.Request(definitions.MCSM_LIG_URL, data_bytes)
            mcsm_lig_results.add_header('Content-type', form.get_content_type())
            mcsm_lig_results.add_header('Content-length', len(data_bytes))
            mcsm_lig_response = urllib.request.urlopen(mcsm_lig_results)
        except Exception as Argument:
            # error not raised as error could be due to snp
            print("An error has occurred \n%s" % Argument)
            return (_complex_pdb_file_name, _lig_pred_affinity, wild_type_affinity,\
                    definitions.SERVER_ERROR, definitions.SERVER_ERROR, definitions.SERVER_ERROR, \
                    definitions.SERVER_ERROR, definitions.SERVER_ERROR, definitions.SERVER_ERROR, \
                    definitions.SERVER_ERROR, definitions.SERVER_ERROR)
            #raise
        # Obtain response
        pred_affin_change_value = definitions.MISSING_ANALYSIS_VALUE
        wild_type_value = definitions.MISSING_ANALYSIS_VALUE
        position_value = definitions.MISSING_ANALYSIS_VALUE
        mutant_type_value = definitions.MISSING_ANALYSIS_VALUE
        chain_value = definitions.MISSING_ANALYSIS_VALUE
        ligand_id_value = definitions.MISSING_ANALYSIS_VALUE
        distance_to_ligand_value = definitions.MISSING_ANALYSIS_VALUE
        duet_stability_change_value = definitions.MISSING_ANALYSIS_VALUE
        try:
            beautiful_soup = BeautifulSoup(mcsm_lig_response, 'html.parser')
            #print(beautiful_soup.prettify())
            count = -1
            font_anchors = beautiful_soup.body.find_all('font')
            for font_anchor in font_anchors:
                count += 1
                if count == 0:
                    body_position = font_anchor
                    pred_affin_change_value = body_position.next_element.string
                    body_position = body_position.next_element
                    pred_affin_change_value = pred_affin_change_value + body_position.next_element.string
                    body_position = body_position.next_element
                    pred_affin_change_value = pred_affin_change_value + body_position.next_element.next_element.string
                elif count == 1:
                    body_position = font_anchor.find_next('b')
                    wild_type_value = body_position.string
                    body_position = body_position.find_next('b')
                    position_value = body_position.string
                    body_position = body_position.find_next('b')
                    mutant_type_value = body_position.string
                    body_position = body_position.find_next('b')
                    chain_value = body_position.string
                    body_position = body_position.find_next('b')
                    ligand_id_value = body_position.string
                    body_position = body_position.find_next('b')
                    distance_to_ligand_value = body_position.string
                    body_position = body_position.find_next('font')
                    duet_stability_change_value = body_position.string
                    duet_stability_change_value = duet_stability_change_value +\
                                                  body_position.next_element.next_element.string.lstrip()
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise
        finally:
            return (_complex_pdb_file_name, _lig_pred_affinity, wild_type_affinity, pred_affin_change_value,\
                    wild_type_value, position_value, mutant_type_value, chain_value, ligand_id_value,\
                    distance_to_ligand_value, duet_stability_change_value)

    def get_file_path(self):
        global _mcsm_lig_results_file_path
        return _mcsm_lig_results_file_path
