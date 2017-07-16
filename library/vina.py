# bola coker
# 29th June, 2017
# class to perform docking
#

import library.definitions as definitions
from pathlib import Path
from zeep import Client
from library.processfile import ProcessFile
import library.utilities as utilities
import urllib.request
import time

class Vina:
    """ class to perform docking with autodock vina """

    _docking_folder=""
    _receptor_pdbqt_file_path= ""
    _vina_opal_response=[]
    _vina_opal_status_response=[]
    _ligand_path_object=""
    _receptor_ligand_pdbqt_file_path= ""
    _receptor_ligand_pdbqt_file_name=""

    """ initialisation class """
    def __init__(self, receptor_pdbqt_file_path, ligand_pdbqt_file_path, docking_folder):
        global _docking_folder
        global _receptor_pdbqt_file_path
        global _ligand_path_object
        global _vina_opal_response
        global _vina_opal_status_response
        _docking_folder=docking_folder
        _receptor_pdbqt_file_path=receptor_pdbqt_file_path
        try:
            # initialise soap client
            opal_client = Client(definitions.OPAL_VINA_WSDL)
            # get cooordinates for config file
            receptor_coords=self.get_coords()
            # write config file
            receptor_path_object=Path(receptor_pdbqt_file_path)
            _ligand_path_object=Path(ligand_pdbqt_file_path)
            config_file_path=docking_folder+definitions.FILE_SEPARATOR+"config.txt"
            config_file=open(config_file_path,"w")
            config_file.write("receptor="+receptor_path_object.name+"\n")
            config_file.write("ligand="+_ligand_path_object.name + "\n\n")
            config_file.write("center_x="+str(receptor_coords[definitions.DICT_RECEPTOR_CTRX])+"\n")
            config_file.write("center_y=" + str(receptor_coords[definitions.DICT_RECEPTOR_CTRY])+"\n")
            config_file.write("center_z=" + str(receptor_coords[definitions.DICT_RECEPTOR_CTRZ])+"\n\n")
            config_file.write("size_x=" + str(receptor_coords[definitions.DICT_RECEPTOR_SIZEX])+"\n")
            config_file.write("size_y=" + str(receptor_coords[definitions.DICT_RECEPTOR_SIZEY])+"\n")
            config_file.write("size_z=" + str(receptor_coords[definitions.DICT_RECEPTOR_SIZEZ])+"\n\n")
            config_file.write("log=vina.out\n\n")
            config_file.write("seed=2017\n")
            config_file.write("exhaustiveness="+str(definitions.VINA_OPTION_EXHAUST)+"\n")
            config_file.write("num_modes="+str(definitions.VINA_OPTION_MODES)+"\n")
            config_file.close()
            # call web service
            print("\nVina command: vina --config=config.txt\n")
            vina_opal_initial_response = opal_client.service.launchJob \
                (argList="--config=config.txt", \
                 inputFile=[{"name": "config.txt", "contents": ProcessFile.encode_file(config_file_path)}, \
                            {"name": receptor_path_object.name,"contents":\
                                ProcessFile.encode_file(receptor_pdbqt_file_path)}, \
                            {"name": _ligand_path_object.name, "contents":\
                                ProcessFile.encode_file(ligand_pdbqt_file_path)}])
            print(vina_opal_initial_response)
            vina_opal_jobid=vina_opal_initial_response["jobID"]
            while True:
                # get opal status
                _vina_opal_status_response=opal_client.service.queryStatus(vina_opal_jobid)
                if utilities.opal_job_running(_vina_opal_status_response):
                    _vina_opal_response=opal_client.service.getOutputs(vina_opal_jobid)
                    print(_vina_opal_response)
                    break
                else:
                    time.sleep(definitions.OPAL_POOLING_TIME)
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    """ method to obtain docking parameters """
    def get_coords(self):
        global _receptor_pdbqt_file_path
        max_x = -99999
        min_x = 99999
        max_y = -99999
        min_y = 99999
        max_z = -99999
        min_z = 99999
        with open(_receptor_pdbqt_file_path, "r") as receptor_pdbqt_file:
            for line in receptor_pdbqt_file:
                if line[definitions.PDB_ATOM_START-1:definitions.PDB_ATOM_END-1].strip()=="ATOM":
                    x_coord=float(line[definitions.PDB_ATOM_XSTART-1:definitions.PDB_ATOM_XEND-1])
                    y_coord = float(line[definitions.PDB_ATOM_YSTART - 1:definitions.PDB_ATOM_YEND - 1])
                    z_coord = float(line[definitions.PDB_ATOM_ZSTART - 1:definitions.PDB_ATOM_ZEND - 1])
                    if x_coord>max_x:
                        max_x=x_coord
                    elif x_coord<min_x:
                        min_x=x_coord
                    if y_coord > max_y:
                        max_y=y_coord
                    elif y_coord < min_y:
                        min_y=y_coord
                    if z_coord > max_z:
                        max_z=z_coord
                    elif z_coord < min_z:
                        min_z=z_coord
        return {definitions.DICT_RECEPTOR_CTRX:round(min_x+((max_x-min_x+1)/2),3),\
                definitions.DICT_RECEPTOR_CTRY:round(min_y+((max_y-min_y)/2),3),\
                definitions.DICT_RECEPTOR_CTRZ:round(min_z+((max_z-min_z)/2),3),\
                definitions.DICT_RECEPTOR_SIZEX:round(max_x-min_x,3),\
                definitions.DICT_RECEPTOR_SIZEY:round(max_y-min_y,3),\
                definitions.DICT_RECEPTOR_SIZEZ:round(max_z-min_z,3)}

    """ base class for opal clients """
    """ get opal server returned error status """
    def get_status(self):
        return _vina_opal_status_response["code"]

    """ get opal server returned error message """
    def get_error(self):
        global _vina_opal_response
        return utilities.get_error(_vina_opal_response)

    """ get vina results and save in file """
    def save_output(self):
        global _vina_opal_response
        global _ligand_path_object
        global _docking_folder
        global _receptor_pdbqt_file_path
        global _receptor_ligand_pdbqt_file_path
        global _receptor_ligand_pdbqt_file_name
        try:
            # save stdOut and stdErr files
            utilities.save_std_files(_vina_opal_response, _docking_folder + definitions.FILE_SEPARATOR + "vina_out")
            # get vina results
            complex_pdb_file_path=""
            opal_poses_file_name=_ligand_path_object.name.replace(".pdbqt","_out.pdbqt")
            for output in _vina_opal_response['outputFile']:
                if output['name'] == opal_poses_file_name:
                    opal_poses = urllib.request.urlopen(output['url'])
                    opal_poses_contents = opal_poses.read()
                    opal_poses.close()
                    # save results to file
                    poses_file_path = _docking_folder + definitions.FILE_SEPARATOR + opal_poses_file_name
                    poses_file = open(poses_file_path, 'w+')
                    poses_file.write(opal_poses_contents.decode("utf-8"))
                    # write complex pdbqt
                    # by copying receptor pdbqt and first pose to a new ligand_receptor pdbqt
                    _receptor_ligand_pdbqt_file_name=Path(_receptor_pdbqt_file_path).name.\
                        replace(".pdbqt","_ligand.pdbqt")
                    _receptor_ligand_pdbqt_file_path=_docking_folder+definitions.FILE_SEPARATOR+\
                                                     _receptor_ligand_pdbqt_file_name
                    receptor_ligand_pdbqt=open(_receptor_ligand_pdbqt_file_path, "w")
                    with open(_receptor_pdbqt_file_path, "r") as receptor_pdbqt:
                        for line in receptor_pdbqt:
                            receptor_ligand_pdbqt.write(line)
                    poses_file.seek(0)
                    # change ATOM in ligand.pdbqt file to HETATM
                    # assign ligand to chain A
                    while True:
                        line=poses_file.readline()
                        line=line.replace("LIG    1","LIG A  1")
                        "LIG    1"
                        if line != "MODEL 2\n":
                            if line.count("ATOM      ") == 1:
                                receptor_ligand_pdbqt.write(line.replace("ATOM      ", "HETATM    "))
                            elif line.count("ATOM     ") == 1:
                                receptor_ligand_pdbqt.write(line.replace("ATOM     ", "HETATM   "))
                            else:
                                receptor_ligand_pdbqt.write(line)
                        else:
                            break
                    # close files
                    poses_file.close()
                    receptor_ligand_pdbqt.close()
                    # create complex pdb
                    complex_pdb_file_path=self.create_pdb()
                if output['name'] == "vina.out":
                    opal_log = urllib.request.urlopen(output['url'])
                    opal_log_contents = opal_log.read()
                    opal_log.close()
                    # save results to file
                    log_file_path = _docking_folder + definitions.FILE_SEPARATOR + "vina.out"
                    log_file = open(log_file_path, 'w')
                    log_file.write(opal_log_contents.decode("utf-8"))
                    log_file.close()
            # return vina complex file
            return complex_pdb_file_path
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise

    def create_pdb(self):
        global _receptor_ligand_pdbqt_file_path
        global _receptor_ligand_pdbqt_file_name
        try:
            # initialise soap client
            opal_client = Client(definitions.OPAL_OPENBABEL_WSDL)
            complex_pdb_file_path = _receptor_ligand_pdbqt_file_path.replace(".pdbqt", ".pdb")
            complex_pdb_file_name=_receptor_ligand_pdbqt_file_name.replace(".pdbqt", ".pdb")
            arg_list = "-ifile " + _receptor_ligand_pdbqt_file_name + " -iformat -ipdbqt -ofile " + \
                       complex_pdb_file_name + " -oformat -opdb"
            print("\nBabel command: babel %s\n" % arg_list)
            babel_opal_initial_response = opal_client.service.launchJob \
                (argList=arg_list, \
                 inputFile={"name": _receptor_ligand_pdbqt_file_name, "contents":\
                     ProcessFile.encode_file(_receptor_ligand_pdbqt_file_path)})
            print(babel_opal_initial_response)
            babel_opal_jobid=babel_opal_initial_response["jobID"]
            while True:
                # get opal status
                babel_opal_status_response=opal_client.service.queryStatus(babel_opal_jobid)
                if utilities.opal_job_running(babel_opal_status_response):
                    babel_opal_response=opal_client.service.getOutputs(babel_opal_jobid)
                    print(babel_opal_response)
                    break
                else:
                    time.sleep(definitions.OPAL_POOLING_TIME)
            if babel_opal_status_response["code"]==definitions.OPAL_SUCCESS:
                # save stdOut and stdErr files
                utilities.save_std_files(babel_opal_response, _docking_folder\
                                         + definitions.FILE_SEPARATOR + "pdbqt2pdb_out")
                # get open babel results
                opal_pdb = ""
                for output in babel_opal_response['outputFile']:
                    if output['name'] == complex_pdb_file_name:
                        opal_pdb = urllib.request.urlopen(output['url'])
                opal_pdb_contents = opal_pdb.read()
                opal_pdb.close()
                # save results to file
                pdb_file = open(complex_pdb_file_path, 'w')
                pdb_file.write(opal_pdb_contents.decode("utf-8"))
                pdb_file.close()
                return complex_pdb_file_path
            else:
                print(utilities.get_error(babel_opal_response))
        except Exception as Argument:
            print("An error has occurred \n%s" % Argument)
            raise





