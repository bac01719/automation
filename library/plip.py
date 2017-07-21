# bola coker
# 16th July, 2017
#
# wrapper for PLIP website
# https://projects.biotec.tu-dresden.de/plip-web/plip/

from bs4 import BeautifulSoup
from library.multi_part_form_plip import MultiPartForm
import library.definitions as definitions
import os
import io
import urllib.request
import urllib.parse
import time

class PLIP:
    """ Easy and fast identification of noncovalent interactions between proteins and their ligands"""

    """ initialise clas with full path name of complex file """
    def __init__(self,data_analysis_complex_file_path):
        print("\nObtaining noncovalent interactions between protein and ligand\n")
        try:
            # read contents of complex file
            with open(data_analysis_complex_file_path, 'r') as complex_file_handle:
                complex_file = complex_file_handle.read()
            # Create the form with simple fields
            form = MultiPartForm()
            form.add_field('select-pdb', 'by-file')
            form.add_file('file', data_analysis_complex_file_path,\
                          fileHandle=io.BytesIO(complex_file.encode('utf-8')))
            form.add_field('enableReplaceObsolete', "1")
            form.add_field('showAdvancedOptions',"")
            form.add_field('name', "automation_job")
            form.add_field('userEmail', "")
            form.add_field('submit',"Run analysis")
            # build the request
            data_bytes = bytes(form)
            plip_request = urllib.request.Request(definitions.PLIP_URL, data_bytes)
            plip_request.add_header('Content-type', form.get_content_type())
            plip_request.add_header('Content-length', len(data_bytes))
            plip_response = urllib.request.urlopen(plip_request)
            # get url of report page and read
            plip_results_url = plip_response.geturl()
            print("PLIP results URL : %s" % plip_results_url)
            # loop results page till complete
            while True:
                print("Pooling for PLIP results page")
                # read results page
                plip_results = urllib.request.urlopen(plip_results_url)
                plip_results_contents = plip_results.read().decode("utf-8")
                if plip_results_contents.find(definitions.PLIP_RESULTS_TITLE)>0:
                    #print(plip_results_contents)
                    # get url of xml file
                    plip_xml_url=plip_results_url.replace("/result/","/download/")
                    plip_xml_url+="?filePath=outputs%2Freport.xml"
                    print("PLIP results XML URL : %s" % plip_xml_url)
                    plip_xml=urllib.request.urlopen(plip_xml_url)
                    plip_xml_content=plip_xml.read().decode("utf-8")
                    print(plip_xml.read().decode("utf-8"))
                    # store xml to file
                    xml_file_path=data_analysis_complex_file_path.replace(".pdb","_plip.xml")
                    xml_file=open(xml_file_path,'w')
                    xml_file.write(plip_xml_content)
                    xml_file.close()
                    break
                else:
                    time.sleep(definitions.PLIP_POOLING_TIME)
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise



