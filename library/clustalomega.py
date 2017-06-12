# bola coker
# 4th June, 2017
# class to wrap call to opal soap server
#

from zeep import Client
import library.definitions as definitions
from library.processfile import ProcessFile
import urllib.request

class ClustalOmega:
    """ class to wrap call to opal soap server - clustal omega service"""

    _opal_response=[]
    _template_file_path=""

    """ initialise class with full template file path """
    def __init__(self, fasta_file_path, template_file_path):
        global _opal_response
        global _template_file_path
        _template_file_path=template_file_path
        try:
            # initialise soap client
            opal_client=Client(definitions.OPAL_CLUSTALOMEGA_WSDL)
            fasta_file_contents=ProcessFile.encode_file(fasta_file_path)
            template_file_contents=ProcessFile.encode_file(template_file_path)
            # define arguments
            arg_list="-i templateFile.fa -o msa.fa"
            if definitions.CLUSTAL_OMEGA_DEALIGN==True:
                arg_list+=" --dealign"
            if definitions.CLUSTAL_OMEGA_FULL==True:
                arg_list+=" --full"
            if definitions.CLUSTAL_OMEGA_FULLITER==True:
                arg_list+=" --full-iter"
            if definitions.CLUSTAL_OMEGA_ITER>0:
                arg_list+=" --iter "+str(definitions.CLUSTAL_OMEGA_ITER)
            print(arg_list)
            _opal_response=opal_client.service.launchJobBlocking\
                (argList=arg_list,\
                 inputFile={"name":"templateFile.fa","contents":fasta_file_contents+template_file_contents})
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    """ get opal server returned error status """
    def get_status(self):
        global _opal_response
        return _opal_response["status"]["code"]

    """ get opal server returned error message """
    def get_error(self):
        global _opal_response
        try:
            # get opal error link and read contents
            opal_error=urllib.request.urlopen(_opal_response["jobOut"]["stdErr"])
            opal_error_contents=opal_error.read()
            opal_error.close()
            return opal_error_contents.decode("utf-8")
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    """ get clustal omega results and save in file """
    def save_msa(self):
        global _opal_response
        global _template_file_path
        try:
            # get clustal omega results
            opal_msa=urllib.request.urlopen(_opal_response['jobOut']['outputFile']['name'=='msa.fa']['url'])
            opal_msa_contents=opal_msa.read()
            opal_msa.close()
            # save results to file
            msa_file_path=_template_file_path.replace(definitions.HOMOLOGY_TEMPLATE_FILE_EXTENSION,\
                                                        definitions.MSA_FILE_EXTENSION)
            msa_file=open(msa_file_path,'w')
            msa_file.write(opal_msa_contents.decode("utf-8"))
            msa_file.close()
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise








