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
DICT_BLAST_INCHIKEY_FOUND="blast_inchikey_found"
DICT_BLAST_PROTEIN_RESOLUTION="blast_protein_resolution"

# OPAL definitions
OPAL_CLUSTALOMEGA_WSDL="/home/bola/Documents/Private/BBK/project/python/library/xsd/ClustalOmegaService.xml"
OPAL_SUCCESS=8  # success code from opal servers

