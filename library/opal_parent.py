# bola coker
# June 22nd 2017
# abstract class for opal client

import library.utilities as utilities

class Opal_Parent:

    global _opal_response

    """ base class for opal clients """
    """ get opal server returned error status """
    def get_status(self):
        return utilities.get_status(_opal_response)

    """ get opal server returned error message """
    def get_error(self,_opal_response):
        return utilities.get_error(_opal_response)

