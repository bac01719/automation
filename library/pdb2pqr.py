# bola coker 28th June,2017
# pdb2pqr module
# to protonate hydrogens in protein before sending to prepare receptor webservice
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

class PDB2PQR:

    """ class to wrap opal service pdb2pqr """

    _docking_folder=""
    _opal_response=[]
    _opal_status_response=[]
    _pqr_file_name=""

    """ constructor for pdb2pqr opal service wrapper """
    def __init__(self,docking_receptor_file_path,docking_folder):
        global _opal_response
        global _opal_status_response
        global _docking_folder
        global _pqr_file_name
        _docking_folder=docking_folder
        try:
            # initialise soap client
            opal_client = Client(definitions.OPAL_PDB2PQR_WSDL)
            # create docking folder
            if not os.path.exists(docking_folder):
                os.mkdir(docking_folder)
            # copy receptor to docking folder
            path_object=Path(docking_receptor_file_path)
            pdb_file_path=docking_folder+definitions.FILE_SEPARATOR+path_object.name
            _pqr_file_name=path_object.name.replace(".pdb",".pqr")
            copyfile(docking_receptor_file_path,pdb_file_path)
            # set up argument list for pdb2pqr
            arg_list="--verbose --chain --summary --hbond --drop-water --ff=amber --ph-calc-method=propka "+\
                     path_object.name+" "+_pqr_file_name
            print("\nPDB2PQR command: pdb2pqr.py %s\n" % arg_list)
            opal_initial_response=opal_client.service.launchJob\
                (argList=arg_list,\
                 inputFile={"name":path_object.name,"contents":ProcessFile.encode_file(pdb_file_path)})
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
        return utilities.get_error(_opal_response)

    """ get pqr results and save in file """
    def save_output(self):
        global _opal_response
        global _pqr_file_name
        global _docking_folder
        try:
            # save stdOut and stdErr files
            utilities.save_std_files(_opal_response, _docking_folder+definitions.FILE_SEPARATOR+"pdb2pqr_out")
            # get pdb2pqr results
            for output in _opal_response['outputFile']:
                if output['name']==_pqr_file_name:
                    opal_pqr=urllib.request.urlopen(output['url'])
                    opal_pqr_contents = opal_pqr.read()
                    opal_pqr.close()
                    # save pqr results to file
                    pqr_file_path = _docking_folder + definitions.FILE_SEPARATOR + _pqr_file_name
                    pqr_file = open(pqr_file_path, 'w')
                    pqr_file.write(opal_pqr_contents.decode("utf-8"))
                    pqr_file.close()
                if output['name']==_pqr_file_name.replace(".pqr",".hbond"):
                    opal_hbond=urllib.request.urlopen(output['url'])
                    opal_hbond_contents = opal_hbond.read()
                    opal_hbond.close()
                    # save hbond results to file
                    hbond_file_path = _docking_folder + definitions.FILE_SEPARATOR +\
                                      _pqr_file_name.replace(".pqr",".hbond")
                    hbond_file = open(hbond_file_path, 'w')
                    hbond_file.write(opal_hbond_contents.decode("utf-8"))
                    hbond_file.close()
                if output['name']==_pqr_file_name.replace(".pqr",".propka"):
                    opal_propka=urllib.request.urlopen(output['url'])
                    opal_propka_contents = opal_propka.read()
                    opal_propka.close()
                    # save propka results to file
                    propka_file_path = _docking_folder + definitions.FILE_SEPARATOR + \
                                      _pqr_file_name.replace(".pqr", ".propka")
                    propka_file = open(propka_file_path, 'w')
                    propka_file.write(opal_propka_contents.decode("utf-8"))
                    propka_file.close()
                if output['name']==_pqr_file_name.replace(".pqr",".summary"):
                    opal_summary=urllib.request.urlopen(output['url'])
                    opal_summary_contents = opal_summary.read()
                    opal_summary.close()
                    # save summary results to file
                    summary_file_path = _docking_folder + definitions.FILE_SEPARATOR + \
                                       _pqr_file_name.replace(".pqr", ".summary")
                    summary_file = open(summary_file_path, 'w')
                    summary_file.write(opal_summary_contents.decode("utf-8"))
                    summary_file.close()
            # return pqr file path
            return pqr_file_path
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise







