# bola coker
# 19/05/17
#
# definitions of constants used

# general
EMAIL_ADDRESS="bac01719.yahoo.com"
FILE_SEPARATOR='/'
HOMOLOGY_PERCENT_IDENTITY_LOWER_RANGE=30   # sequences larger than or equal to this are part of homology template
HOMOLOGY_PERCENT_IDENTITY_UPPER_RANGE=95   # sequences smaller than or equal to this are part of the homology template
HOMOLOGY_MINIMUM_RESOLUTION=2.0 # sequences with resolutions larger than this are excluded
OTHER_MINIMUM_RESOLUTION=2.0 # sequences with resolutions larger than this are not considered for docking or analysis
MINIMUM_PERCENT_IDENTITY_FOR_OTHER=95   # sequences larger than or equal to this are considered for docking or analysis
FASTA_FILE_EXTENSION=".fasta"   # fasta file extension
BLAST_FILE_EXTENSION=".blast"   # blast file extension
HOMOLOGY_TEMPLATE_FILE_EXTENSION=".template"    # homology template file extension
MSA_FILE_EXTENSION=".fa" #multisequence file alignment extension (e.g. ClustalOmega)
PDBQT_LIGAND_FILE="ligand.pdbqt"   # pdbqt file extension for ligand and receptor

# pdb file format
PDB_ATOM_START=1                # atom definition start
PDB_ATOM_END=6                  # atom definition end
PDB_ATOM_XSTART=31              # atom x coord start. format
PDB_ATOM_XEND=38                # atom x coord end. format
PDB_ATOM_YSTART=39              # atom y coord start. format
PDB_ATOM_YEND=46                # atom y coord end. format
PDB_ATOM_ZSTART=47              # atom z coord start. format
PDB_ATOM_ZEND=54                # atom z coord end. format

# blast module options
BLAST_PROGRAM='blastp'   #options are blastp and blastx
BLAST_PERCENTAGE_IDENTITY_CUTOFF=30     # sequences above or equal to this figure are returned by BLAST

# homology module
HOMOLOGY_FOLDER="homology"
HOMOLOGY_MONOMER_CHAIN="A"

# databases
ORGANISM_DATABASE="taxonomy"
GENE_DATABASE="gene"
PROTEIN_DATABASE="protein"

# dictionary keys
DICT_PROTEIN_NAME="protein_name"
DICT_PROTEIN_NCBI_ID="protein_ncbi_id"
DICT_GENE_DISCONTINUED="gene_discontinued"
DICT_BLAST_PERCENTAGE_IDENTITY="blast_percentage_identity"
DICT_BLAST_PROTEIN="blast_protein"
DICT_BLAST_SMILES_FOUND="blast_smiles_found"
DICT_BLAST_PROTEIN_RESOLUTION="blast_protein_resolution"
DICT_FASTA_SEQUENCE="fasta"
DICT_FASTA_FILE_PATH="fasta_file_path"
DICT_HOMOLOGY_ZDOPE="max_zDOPE"
DICT_HOMOLOGY_BEST_PATH="best_homology_file_path"
DICT_RECEPTOR_CTRX="receptor_center_x"
DICT_RECEPTOR_CTRY="receptor_center_y"
DICT_RECEPTOR_CTRZ="receptor_center_z"
DICT_RECEPTOR_SIZEX="receptor_size_x"
DICT_RECEPTOR_SIZEY="receptor_size_y"
DICT_RECEPTOR_SIZEZ="receptor_size_z"
DICT_FASTA_LIST="fasta_list"
DICT_FASTA_PROTEIN_FOLDER="fasta_protien_folder"

# OPAL definitions
OPAL_CLUSTALOMEGA_WSDL="/home/bola/Documents/Private/BBK/project/python/library/xsd/ClustalOmegaService.xml"
OPAL_MODELLER_WSDL="/home/bola/Documents/Private/BBK/project/python/library/xsd/Modeller9v8Service.xml"
OPAL_OPENBABEL_WSDL="/home/bola/Documents/Private/BBK/project/python/library/xsd/openbabel_2.3.1.xml"
OPAL_PDB2PQR_WSDL="/home/bola/Documents/Private/BBK/project/python/library/xsd/pdb2pqr_2.1.1.xml"
OPAL_PREPRECEP_WSDL="/home/bola/Documents/Private/BBK/project/python/library/xsd/prepare_receptor_1.5.6.xml"
OPAL_VINA_WSDL="/home/bola/Documents/Private/BBK/project/python/library/xsd/vina_1.1.2.xml"
OPAL_SUCCESS=8  # success code from opal servers

# ClustalOmega keys
CLUSTAL_OMEGA_DEALIGN=True   # Dealign input sequences
CLUSTAL_OMEGA_FULL=True # Use full distance matrix for guide-tree calculation (might be slow; mBed is default)
CLUSTAL_OMEGA_FULLITER=True # Use full distance matrix for guide-tree calculation during iteration
                            # (might be slowish; mBed is default)
CLUSTAL_OMEGA_ITER=3    # Number of (combined guide-tree/HMM) iterations (-9 for none)

# Modeller contants
MODELLER_PIR_EXTENSION=".pir"
MODELLER_LICENCE_FILE="/home/bola/Documents/Private/BBK/project/other/modeller.lic"
MODELLER_AUTOLOOP_PYTHON_FILE="/home/bola/Documents/Private/BBK/project/python/library/modeller/ModellerModelling.py"

# OpenBabel options for ligand pdbqt creation
OBABEL_OPTION_DELHYD=False  # delete hydrogens
OBABEL_OPTION_ADDHYD=False   # add hydrogen
OBABEL_OPTION_ADDHYD_PH=True   # Add Hydrogens appropriate for pH (use transforms in phmodel.txt)
OBABEL_OPTION_CENTRE=False   # Center atomic coordinates (** note ADD_HYD_PH and CENTRE are mutually exclusive)

# prepare receptor options for protein pdbqt file creation

# repairs - only one option can be selected. Options are
# bonds_hydrogens - build bonds and add hydrogens
# bonds - build a single bond from each atom with no bonds to its closest neighbour
# hydrogens - add hydrogens
# checkhydrogens - add hydrogens only if there are none already
# None - do not make any repairs
PRERECEP_OPTION_REPAIRS='checkhydrogens'

PRERECEP_OPTION_INPUTCHGS=False     # preserve all input charges i.e. do not add new charges
PRERECEP_OPTION_PRESERVE=""         # preserve all input charges on specifc atom types e.g. -p Zn, -p Fe

# clean up - multiple options can be selected. Options are
# nphs - merge charges and remove non-polar hydrogens
# lps - merge charges and remove lone pairs
# waters - remove water residues
# nonstdres - remove chains composed entirely of residues of types other than the standard 20 amino acids
# deleteAltB - remove XX@B atoms and rename XX@A atoms->XX
PRERECEP_OPTION_CLEANUP='nphs_lps_waters_nonstdres'

PRERECEPT_OPTION_NONSTDCHAIN=False   # delete every nonstd residue from any chain. Options are False and True

# autodock vina options - note blind docking is used
VINA_OPTION_EXHAUST=50    # exhaustiveness of the global search (roughly proportional to time)
VINA_OPTION_MODES=100       # maximum number of binding modes to generate







