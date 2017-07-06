# bola coker
# 6th June, 2016
# class do perform homology

import os
import library.definitions as definitions
from pathlib import Path
import re
from zeep import Client
from library.processfile import ProcessFile
import library.utilities as utilities
import urllib.request

class Homology:
    """ class to handle homology methods"""

    _modeller_folder=""
    _pir_file_path=""
    _modeller_config_file_path=""
    _protein_template_file_path=""
    _protein_list=[]
    _opal_response=[]

    """ initialise method of homology class """
    def __init__(self,msa_file_path):
        global _modeller_folder
        global _pir_file_path
        global _modeller_config_file_path
        global _protein_template_file_path
        global _protein_list
        try:
            # check msa file path exist
            if not os.path.exists (msa_file_path):
                raise FileNotFoundError
            # create modeller folder and contents
            msa_file_path_object=Path(msa_file_path)
            _modeller_folder=str(msa_file_path_object.parent)+definitions.FILE_SEPARATOR+"modeller"
            if not os.path.exists(_modeller_folder):
                os.mkdir(_modeller_folder)
            # create protein_list file
            _protein_template_file_path=_modeller_folder+definitions.FILE_SEPARATOR+"proteins.txt"
            protein_template_file=open(_protein_template_file_path,"w")
            # create empty protein list dictionary
            _protein_list=[]
            # create PIR alignment file from fasta alignment file
            pir_file=str(msa_file_path_object.name).replace(str(msa_file_path_object.suffix),\
                                                            definitions.MODELLER_PIR_EXTENSION)
            _pir_file_path=_modeller_folder+definitions.FILE_SEPARATOR+pir_file
            msa_file=open(msa_file_path,"r")
            pir_file=open(_pir_file_path,"w+")
            # first pass create PIR file but without asterick at end of sequence
            for line in msa_file.readlines():
                if line[0]==">":
                    line_re = re.compile(">([^|.]+)")
                    line_protein = line_re.findall(line)[0].strip()
                    if line.find("|A")>-1:  # pdb sequence for template
                        protein_template_file.write(line_protein + "\n")
                        _protein_list.append(line_protein + ".pdb")
                        line = line.replace(">", ">P1;")
                        line = line.replace("|A", "")
                        pir_file.write(line)
                        pir_file.write("structure:" + line_protein + ":FIRST:@:LAST:@::::\n")
                    else:   # target sequence
                        protein_template_file.write(line_protein+"\n")
                        line = line.replace(">", ">P1;")
                        pir_file.write(line)
                        pir_file.write("sequence:" + line_protein + ":.:.:.:.::::\n")
                else:
                    pir_file.write(line)
            # second pass for pir file to add asterick at end of sequences
            pir_file.seek(0)
            pir_file_content=pir_file.read()
            reg_object=re.compile(r"(\n>)")
            pir_file_content=reg_object.sub(r"*\n>",pir_file_content)
            reg_object=re.compile(r"(\n$)")
            pir_file_content=reg_object.sub(r"*\n",pir_file_content)
            pir_file.seek(0)
            pir_file.write(pir_file_content)
            # close files
            pir_file.close()
            msa_file.close()
            protein_template_file.close()
            # create modeller config file
            Homology.__write_modeller_config(self)
            # call opal soap web service
            Homology.__call_opal_service(self)
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise


    """ method to write xml config file for modeller opal server """
    def __write_modeller_config(self):
        global _modeller_folder
        global _modeller_config_file_path
        try:
            # get modeller licence
            modeller_licence_file=open(definitions.MODELLER_LICENCE_FILE,"r")
            modeller_licence=modeller_licence_file.readline().strip()
            _modeller_config_file_path=_modeller_folder+definitions.FILE_SEPARATOR+"ModellerScriptConfig.xml"
            modeller_config_file=open(_modeller_config_file_path,"w")
            modeller_config_file.write('<?xml version="1.0" encoding="UTF-8"?>')
            modeller_config_file.write('<modeller9v8>')
            modeller_config_file.write('<key>'+modeller_licence+'</key>')
            modeller_config_file.write('<version>2</version>')
            modeller_config_file.write('<numModel>5</numModel>')
            modeller_config_file.write('<hetAtom>0</hetAtom>')
            modeller_config_file.write('<water>0</water>')
            modeller_config_file.write('<allHydrogen>0</allHydrogen>')
            modeller_config_file.write('<veryFast>0</veryFast>')
            modeller_config_file.write('<loopInfo>0</loopInfo>')
            modeller_config_file.write('</modeller9v8>')
            modeller_config_file.close()
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    """ method to call opal webservice """
    def __call_opal_service(self):
        global _modeller_folder
        global _pir_file_path
        global _modeller_config_file_path
        global _protein_template_file_path
        global _protein_list
        global _opal_response
        try:
            opal_client = Client(definitions.OPAL_MODELLER_WSDL)
            alignment_file_content = ProcessFile.encode_file(_pir_file_path)
            python_file_content = ProcessFile.encode_file(definitions.MODELLER_AUTOLOOP_PYTHON_FILE)
            config_file_content = ProcessFile.encode_file(_modeller_config_file_path)
            protein_list_file_content = ProcessFile.encode_file(_protein_template_file_path)
            input_file=[{"name":"alignment.ali","contents":alignment_file_content},\
                       {"name":"ModellerModelling.py","contents":python_file_content},\
                       {"name":"ModellerScriptConfig.xml","contents":config_file_content},\
                       {"name":"namelist.dat","contents":protein_list_file_content}]
            # read proteins to inputFile dictionary list
            for protein in _protein_list:
                protein_file_path=str(Path(_modeller_folder).parent)+definitions.FILE_SEPARATOR+protein
                input_file.append({"name":protein,"contents":ProcessFile.encode_file(protein_file_path)})
            # call opal client
            print("\nModeller command : ModellerScriptConfig.xml\n")
            _opal_response=opal_client.service.launchJobBlocking(argList="ModellerScriptConfig.xml", \
                                                            inputFile=input_file)
            print(_opal_response)
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    """ get opal server returned error status """
    def get_status(self):
        global _opal_response
        return utilities.get_status(_opal_response)

    """ get opal server returned error message """
    def get_error(self):
        global _opal_response
        return utilities.get_error(_opal_response)

    """ save output from opal server """
    def save_output(self):
        global _opal_response
        global _modeller_folder
        try:
            # save stdOut and stdErr files
            utilities.save_std_files(_opal_response, _modeller_folder)
            # get base url
            opal_base=_opal_response['status']['baseURL']
            # get index for ok_models.dat output
            ok_models_index=-1
            counter=-1
            for name in _opal_response['jobOut']['outputFile']:
                counter+=1
                if name['name']=='ok_models.dat':
                    ok_models_index=counter
            # get modeller results ok models and save results to file
            ok_models=urllib.request.urlopen(_opal_response['jobOut']['outputFile'][ok_models_index]['url'])
            ok_models_content=ok_models.read()
            ok_models_file_path = _modeller_folder+definitions.FILE_SEPARATOR+"ok_models.dat"
            ok_models_file = open(ok_models_file_path, 'w+')
            ok_models_file.write(ok_models_content.decode("utf-8"))
            ok_models.close()
            # read ok model file into a dictionary list
            ok_models_list=[]
            ok_models_file.seek(0)
            for lines in ok_models_file.readlines()[1:]:
                line_split=lines.strip().split("\t")
                print({"pdb_name":line_split[0],"GA341":line_split[1],"zDOPE":line_split[2]})
                ok_models_list.append({"pdb_name":line_split[0],"GA341":line_split[1],"zDOPE":line_split[2]})
            ok_models_file.close()
            # download ok model's pdb files and return oath of best homology model
            best_homology_model_path=""
            zDOPE=99999.00
            for pdb in ok_models_list:
                pdb_url=urllib.request.urlopen(opal_base+"/"+pdb['pdb_name'])
                pdb_content=pdb_url.read()
                pdb_file_path=_modeller_folder+definitions.FILE_SEPARATOR+pdb["pdb_name"]
                pdb_file = open(pdb_file_path, 'w')
                pdb_file.write(pdb_content.decode("utf-8"))
                if float(pdb["zDOPE"])<float(zDOPE):
                    best_homology_model_path=pdb_file_path
                    zDOPE=pdb["zDOPE"]
                pdb_url.close()
                pdb_file.close()
            return {definitions.DICT_HOMOLOGY_BEST_PATH:best_homology_model_path,definitions.DICT_HOMOLOGY_ZDOPE:zDOPE}
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise








