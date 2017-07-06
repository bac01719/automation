# bola coker 22nd June,2017
# openbabel module
#
from zeep import Client
from library.processfile import ProcessFile
import urllib.request
import library.definitions as definitions
import os
import library.utilities as utilities

class OpenBabel():
    """ class to wrap opal service open babel """

    _docking_folder=""
    _opal_response=""

    """ constructor for openbabel opal service wrapper """
    def __init__(self,ligand_smiles,docking_folder):
        global _opal_response
        global _docking_folder
        _docking_folder=docking_folder
        try:
            # initialise soap client
            opal_client = Client(definitions.OPAL_OPENBABEL_WSDL)
            # create docking folderpath
            if not os.path.exists(docking_folder):
                os.mkdir(docking_folder)
            # store smiles in a file
            smiles_path=docking_folder+definitions.FILE_SEPARATOR+"smile.txt"
            smiles_file=open(smiles_path,"w")
            smiles_file.write(ligand_smiles+"\n")
            smiles_file.close()
            # get 3d inchi pdbqt model
            arg_options="--gen3d "
            if definitions.OBABEL_OPTION_DELHYD:
                arg_options+="-d "
            if definitions.OBABEL_OPTION_ADDHYD:
                arg_options +="-h "
            if definitions.OBABEL_OPTION_ADDHYD_PH:
                arg_options +="-p "
            if definitions.OBABEL_OPTION_CENTRE and not definitions.OBABEL_OPTION_ADDHYD_PH:
                arg_options +="-c "
            arg_list="-ifile smile.txt -iformat -ismi -ofile "+definitions.PDBQT_LIGAND_FILE+" -oformat -opdbqt "+" "+arg_options
            print("\nBabel command: babel %s\n" % arg_list)
            _opal_response=opal_client.service.launchJobBlocking \
                (argList=arg_list,\
                 inputFile={"name":"smile.txt","contents":ProcessFile.encode_file(smiles_path)})
            print(_opal_response)
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    """ base class for opal clients """
    """ get opal server returned error status """
    def get_status(self):
        global _opal_response
        return utilities.get_status(_opal_response)

    """ get opal server returned error message """
    def get_error(self):
        global _opal_response
        return utilities.get_error(_opal_response)

    """ get babel results and save in file """
    def save_ligand(self):
        global _opal_response
        global _template_file_path
        try:
            # save stdOut and stdErr files
            utilities.save_std_files(_opal_response, _docking_folder+\
                                     definitions.FILE_SEPARATOR + "smile2pdbqt_out")
            # get open babel results
            opal_ligand=""
            for output in _opal_response['jobOut']['outputFile']:
                if output['name']==definitions.PDBQT_LIGAND_FILE:
                    opal_ligand=urllib.request.urlopen(output['url'])
            opal_ligand_contents = opal_ligand.read()
            opal_ligand.close()
            # save results to file
            ligand_file_path = _docking_folder+definitions.FILE_SEPARATOR+definitions.PDBQT_LIGAND_FILE
            ligand_file = open(ligand_file_path, 'w')
            ligand_file.write(opal_ligand_contents.decode("utf-8"))
            ligand_file.close()
            return ligand_file_path
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise







