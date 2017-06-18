# bola coker
# 9th June, 2017
#
# utilities.py
# includes general functions
import urllib.request

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

