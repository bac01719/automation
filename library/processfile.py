# bola coker
# May 29th 2017
#
# method to act of file

import base64

""" class to process files """
class ProcessFile:

    """ method to escape a text file """
    def encode_file(file_name):
        try:
            text_file=open(file_name,"rb")
            text_file_contents=text_file.read()
            text_file.close()
            #return base64.encodebytes(text_file_contents)
            return text_file_contents
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise
