# bola coker
# 9th June, 2017
#
# utilities.py
# includes general functions
import urllib.request
import library.definitions as definitions
import os

""" function to return a string given a dictionary """
def return_string(dictionary_list):
    dictionary_string=""
    for key_value in dictionary_list:
        if dictionary_list[key_value]!="":
            dictionary_string+="{} ".format(key_value)+dictionary_list[key_value]
        elif dictionary_list[key_value]=="":
            dictionary_string+=key_value
        dictionary_string+=" "
    return dictionary_string

""" get opal server returned error status """
def get_status(opal_response):
    return opal_response["status"]["code"]

""" get opal server returned error message """
def get_error(opal_response):
    try:
        # get opal error link and read contents
        opal_error=urllib.request.urlopen(opal_response["jobOut"]["stdErr"])
        opal_error_contents=opal_error.read()
        opal_error.close()
        return opal_error_contents.decode("utf-8")
    except Exception as Argument:
        print("An error has occurred \n%s" % Argument)
        raise

""" store stdOut and stdErr files """
def save_std_files(opal_response,parent_folder):
    try:
        # check parent folder exist. if not create
        if not os.path.exists(parent_folder):
            os.mkdir(parent_folder)
        # read stdout
        opal_stdout = urllib.request.urlopen(opal_response['jobOut']['stdOut'])
        opal_stdout_contents = opal_stdout.read()
        opal_stdout.close()
        # save stdout results to file
        stdout_file_path = parent_folder + definitions.FILE_SEPARATOR + "stdOut.txt"
        stdout_file = open(stdout_file_path, "w")
        stdout_file.write(opal_stdout_contents.decode("utf-8"))
        stdout_file.close()
        # read stderr
        opal_stderr = urllib.request.urlopen(opal_response['jobOut']['stdErr'])
        opal_stderr_contents = opal_stderr.read()
        opal_stderr.close()
        # save stderr results to file
        stderr_file_path = parent_folder + definitions.FILE_SEPARATOR + "stdErr.txt"
        stderr_file = open(stderr_file_path, "w")
        stderr_file.write(opal_stderr_contents.decode("utf-8"))
        stderr_file.close()
    except Exception as Argument:
        print("An error has occurred \n%s" % Argument)
        raise