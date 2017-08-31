# bola coker 28th June,2017
# prepare receptor module
# tp prepare receptor and create .pdbqt file for docking with autodock vina
#
from zeep import Client
from library.processfile import ProcessFile
import urllib.request
import library.definitions as definitions
import os
import library.utilities as utilities
from shutil import copyfile
from pathlib import Path
import time

class PrepareReceptor:

    """ class to wrap opal service prepare receptor """

    _docking_folder=""
    _opal_response=[]
    _opal_status_response=[]
    _pdbqt_file_name=""

    """ constructor for prepare_receptor opal service wrapper """
    def __init__(self,pqr_file_path,docking_folder):
        global _opal_response
        global _opal_status_response
        global _docking_folder
        global _pdbqt_file_name
        _docking_folder=docking_folder
        try:
            # initialise soap client
            opal_client = Client(definitions.OPAL_PREPRECEP_WSDL)
            # create docking folder
            if not os.path.exists(docking_folder):
                os.mkdir(docking_folder)
            # copy pqr to docking folder if not already there
            pqr_path_object=Path(pqr_file_path)
            _pdbqt_file_name=pqr_path_object.name.replace(".pqr",".pdbqt")
            docking_pqr_file_path=docking_folder+definitions.FILE_SEPARATOR+pqr_path_object.name
            if not os.path.isfile(docking_pqr_file_path):
                # copy pqr to docking folder
                copyfile(pqr_file_path,docking_pqr_file_path)
            # set up argument list for pdb2pqr
            arg_list="-r "+pqr_path_object.name+" -o "+_pdbqt_file_name+" -v " +\
                "-A "+definitions.PRERECEP_OPTION_REPAIRS+" -U "+definitions.PRERECEP_OPTION_CLEANUP +\
                " -e "+str(definitions.PRERECEPT_OPTION_NONSTDCHAIN)
            if definitions.PRERECEP_OPTION_INPUTCHGS==True:
                arg_list+=" -C"
            if definitions.PRERECEP_OPTION_PRESERVE!="":
                arg_list+=(" "+definitions.PRERECEP_OPTION_PRESERVE)
            print("\nPrepare Receptor command: prepare_receptor4.py %s\n" % arg_list)
            opal_initial_response=opal_client.service.launchJob \
                (argList=arg_list,\
                 inputFile={"name":pqr_path_object.name,"contents":ProcessFile.encode_file(pqr_file_path)})
            print(opal_initial_response)
            opal_jobid=opal_initial_response["jobID"]
            while True:
                # get opal status
                _opal_status_response=opal_client.service.queryStatus(opal_jobid)
                if utilities.opal_job_running(_opal_status_response):
                    _opal_response=opal_client.service.getOutputs(opal_jobid)
                    print(_opal_response)
                    break
                else:
                    time.sleep(definitions.OPAL_POOLING_TIME)
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    """ base class for opal clients """
    """ get opal server returned error status """
    def get_status(self):
        return _opal_status_response["code"]

    """ get opal server returned error message """
    def get_error(self):
        global _opal_response
        return utilities.get_error(_opal_response)

    """ get pdbqt results and save in file """
    def save_output(self):
        global _opal_response
        global _pdbqt_file_name
        global _docking_folder
        try:
            # save stdOut and stdErr files
            utilities.save_std_files(_opal_response, _docking_folder +\
                                     definitions.FILE_SEPARATOR + "prepare_receptor_out")
            # get pdb2pqr results
            opal_pdbqt=""
            for output in _opal_response['outputFile']:
                if output['name']==_pdbqt_file_name:
                    opal_pdbqt=urllib.request.urlopen(output['url'])
                    opal_pdbqt_contents = opal_pdbqt.read()
                    opal_pdbqt.close()
                    # save pdbqt results to file
                    pdbqt_file_path = _docking_folder + definitions.FILE_SEPARATOR + _pdbqt_file_name
                    pdbqt_file = open(pdbqt_file_path, 'w')
                    pdbqt_file.write(opal_pdbqt_contents.decode("utf-8"))
                    pdbqt_file.close()
            # return pdbqt file path
            return pdbqt_file_path
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise







