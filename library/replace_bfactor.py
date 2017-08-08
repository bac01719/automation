# bola Coker
# 14th July, 2017
#
# class to create pdb files with bfactor replaced with
# delta delta G values from mCSM, SDM, DUET and mCSM-Lig

import re
import library.definitions as definitions
import os

class ReplaceBFactor:
    """ class to replace pdb bfactor with delet delta G values of
        mCSM, SDM, DUET and mCSM-Lig
    """

    _complex_pdb_file_path=""
    _pdb_mcsm_file_path=""
    _pdb_sdm_file_path = ""
    _pdb_duet_file_path = ""
    _pdb_mcsm_lig_file_path = ""

    def __init__(self,complex_pdb_file_path,duet_csv_file_path,mcsm_lig_csv_file_path):
        global _complex_pdb_file_path
        global _pdb_mcsm_file_path
        global _pdb_sdm_file_path
        global _pdb_duet_file_path
        global _pdb_mcsm_lig_file_path
        _complex_pdb_file_path=complex_pdb_file_path
        print("Replace B factor in protein files\n\n")
        try:
            complex_file_name=os.path.basename(complex_pdb_file_path)
            data_analysis_folder=os.path.dirname(duet_csv_file_path)
            # create file handle for modified pdb complex with mcsm values
            _pdb_mcsm_file_path=data_analysis_folder+definitions.FILE_SEPARATOR+\
                               complex_file_name.replace(".pdb","_MCSM.pdb")
            # create file handle for modified pdb complex with sdm values
            _pdb_sdm_file_path = data_analysis_folder+definitions.FILE_SEPARATOR+\
                               complex_file_name.replace(".pdb","_SDM.pdb")
            # create file handle for modified pdb complex with duet values
            _pdb_duet_file_path = data_analysis_folder+definitions.FILE_SEPARATOR+\
                               complex_file_name.replace(".pdb","_DUET.pdb")
            # create file handle for modified pdb complex with csm lig values
            _pdb_mcsm_lig_file_path = data_analysis_folder+definitions.FILE_SEPARATOR+\
                               complex_file_name.replace(".pdb","_MCSM_LIG.pdb")
            # regular expression to read delt adelta g value and predicted affinity change from CSV files
            reg_exp = re.compile('\s?(-?\d*\.?\d*)\D*')
            # create dictionary list for mCSM, SDM and DUET values
            mcsm_dictionary={}
            sdm_dictionary={}
            duet_dictionary={}
            mcsm_lig_dictionary={}
            # iterate through DUET csv file
            with open(duet_csv_file_path,"r") as duet_csv:
                line_number=0
                for line in duet_csv:
                    line_number+=1
                    if line_number==1:
                        # skip header
                        continue
                    line_items=line.split(",")
                    if len(line_items)==11: # line has been split properly
                        duet_key = line_items[definitions.DUET_FILE_CHAIN_POS] +\
                                  "_" + line_items[definitions.DUET_FILE_WILDTYPENO_POS]
                        mcsm_matches = reg_exp.findall(line_items[definitions.MCSM_DDG_POS])
                        sdm_matches = reg_exp.findall(line_items[definitions.SDM_DDG_POS])
                        duet_matches = reg_exp.findall(line_items[definitions.DUET_DDG_POS])
                        mcsm_dictionary[duet_key]=\
                            {definitions.DICT_RESNAME: line_items[definitions.DUET_FILE_WILDTYPE_POS],\
                             definitions.DICT_CHAINID: line_items[definitions.DUET_FILE_CHAIN_POS],\
                             definitions.DICT_RESNUM: line_items[definitions.DUET_FILE_WILDTYPENO_POS],\
                             definitions.DICT_STABILITY_CHANGE: mcsm_matches[0]}
                        sdm_dictionary[duet_key] = \
                            {definitions.DICT_RESNAME: line_items[definitions.DUET_FILE_WILDTYPE_POS], \
                             definitions.DICT_CHAINID: line_items[definitions.DUET_FILE_CHAIN_POS], \
                             definitions.DICT_RESNUM: line_items[definitions.DUET_FILE_WILDTYPENO_POS], \
                             definitions.DICT_STABILITY_CHANGE: sdm_matches[0]}
                        duet_dictionary[duet_key] = \
                            {definitions.DICT_RESNAME: line_items[definitions.DUET_FILE_WILDTYPE_POS], \
                             definitions.DICT_CHAINID: line_items[definitions.DUET_FILE_CHAIN_POS], \
                             definitions.DICT_RESNUM: line_items[definitions.DUET_FILE_WILDTYPENO_POS], \
                             definitions.DICT_STABILITY_CHANGE: duet_matches[0]}
                        print("Reading mutation (%s) %s " % ("DUET CSV", duet_key))
            # create dictionary list for MCSM-LIG CSV file
            mcsm_lig_dictionary = {}
            # iterate through csm_lig csv file
            with open(mcsm_lig_csv_file_path, "r") as mcsm_lig_csv:
                line_number=0
                for line in mcsm_lig_csv:
                    line_number += 1
                    if line_number == 1:
                        # skip header
                        continue
                    line_items = line.split(",")
                    mcsm_lig_key = line_items[definitions.MCSMLIG_FILE_CHAIN_POS] + \
                               "_" + line_items[definitions.MCSMLIG_FILE_WILDTYPENO_POS]
                    mcsm_lig_matches = reg_exp.findall(line_items[definitions.MCSM_LIG_DDG_POS])
                    mcsm_lig_dictionary[mcsm_lig_key] = \
                        {definitions.DICT_RESNAME: line_items[definitions.MCSMLIG_FILE_WILDTYPE_POS], \
                         definitions.DICT_CHAINID: line_items[definitions.MCSMLIG_FILE_CHAIN_POS], \
                         definitions.DICT_RESNUM: line_items[definitions.MCSMLIG_FILE_WILDTYPENO_POS],\
                         definitions.DICT_STABILITY_CHANGE: mcsm_lig_matches[0]}
                    print("Reading mutation (%s) %s " % ("mCSM_Lig", mcsm_lig_key))
            # update pdbs
            self.update_pdb(mcsm_dictionary,_pdb_mcsm_file_path)
            self.update_pdb(sdm_dictionary, _pdb_sdm_file_path)
            self.update_pdb(duet_dictionary, _pdb_duet_file_path)
            self.update_pdb(mcsm_lig_dictionary, _pdb_mcsm_lig_file_path)
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    def update_pdb(self,csv_dictionary,new_pdb_file_path):
        global _complex_pdb_file_path
        try:
            # open new pdb_file_path
            new_pdb_file=open(new_pdb_file_path,"w")
            # open pdb complex and iterate through it
            with open(_complex_pdb_file_path,"r") as complex_pdb_file:
                for pdb_line in complex_pdb_file:
                    # select the atom records
                    if pdb_line[:definitions.PDB_ATOM_END] == "ATOM  ":
                        # get data
                        b_factor = pdb_line[definitions.PDB_BFACTOR_START-1:definitions.PDB_BFACTOR_END]
                        res_num = pdb_line[definitions.PDB_RESSEQ_START-1:definitions.PDB_RESSEQ_END]
                        res_type = pdb_line[definitions.PDB_RESNAME_START-1:definitions.PDB_RESNAME_END]
                        chain_id = pdb_line[definitions.PDB_CHAINID-1]
                        pdb_key = chain_id + "_" + res_num.strip()
                        new_pdb_line=""
                        if pdb_key in csv_dictionary.keys():
                            csv_dictionary[pdb_key][definitions.DICT_STABILITY_CHANGE]=\
                                round(float(csv_dictionary[pdb_key][definitions.DICT_STABILITY_CHANGE]),2)
                            csv_dictionary[pdb_key][definitions.DICT_STABILITY_CHANGE] = \
                                '{0:.2f}'.format(csv_dictionary[pdb_key][definitions.DICT_STABILITY_CHANGE])
                            if res_type == csv_dictionary[pdb_key][definitions.DICT_RESNAME]:
                                if "-" in csv_dictionary[pdb_key][definitions.DICT_STABILITY_CHANGE]:
                                    new_pdb_line = pdb_line[:definitions.PDB_ATOM_OCC_END] + " " +\
                                                   csv_dictionary[pdb_key][definitions.DICT_STABILITY_CHANGE]
                                else:
                                    new_pdb_line = pdb_line[:definitions.PDB_ATOM_OCC_END] + "  " + \
                                                   csv_dictionary[pdb_key][definitions.DICT_STABILITY_CHANGE]
                                print("Updating residual %s, position %s with stability change %s" % ( \
                                    csv_dictionary[pdb_key][definitions.DICT_RESNAME], res_num, \
                                    csv_dictionary[pdb_key][definitions.DICT_STABILITY_CHANGE]))
                        else:
                            new_pdb_line = pdb_line[:definitions.PDB_ATOM_OCC_END] + "  0.00"
                        new_pdb_file.write(new_pdb_line+pdb_line[definitions.PDB_BFACTOR_END:len(pdb_line)])
                    else:
                        new_pdb_file.write(pdb_line)
            new_pdb_file.close()
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    def get_file_paths(self):
        global _pdb_mcsm_file_path
        global _pdb_sdm_file_path
        global _pdb_duet_file_path
        global _pdb_mcsm_lig_file_path
        return {definitions.DICT_MCSM_FILE_PATH: _pdb_mcsm_file_path, \
                definitions.DICT_SDM_FILE_PATH: _pdb_sdm_file_path, \
                definitions.DICT_DUET_FILE_PATH: _pdb_duet_file_path, \
                definitions.DICT_MCSM_LIG_FILE_PATH: _pdb_mcsm_lig_file_path}


