################################################################################
## INPUT GENOMES
################################################################################
genetic_code_path="genome_db_test/taxonomy/new_taxdump/"

genbank_files_db_path="genome_db_test/genomes/refseq/viral/"
RefSeq_structure="True"
# genbank_files_db_path='flat_db_Alphavirus/'
# RefSeq_structure="False"
#Structure of the genbank files database : True or False
# False is a list of genbank files in a folder
# True same structure as in RefSeq:
# ── genomes
#     └── refseq
#         └── viral
#             ├── Aedes_flavivirus
#             │   └── latest_assembly_versions
#             │       └── GCF_000885715.1_ViralProj39601
#             │           └── GCF_000885715.1_ViralProj39601_genomic.gbff.gz


interpro_path="/mirror/interpro"

#Create a black list of poorly annotated cds
black_list_conflict_domain='true'
black_list_irrelevant_pattern='true'
confidence_score_treshold=4



################################################################################
## VARIABLES
################################################################################

## extraction and basic stat
tresholdSP="90"
taxon='Flaviviridae'
taxon='Alphavirus'
# taxon='Picornavirales'

# taxon='ssRNA viruses'
# taxon='Retro-transcribing viruses'
# taxon='Picornavirales'

## blast evalue start
blast_evalue='1e-5'

## filtering and mcl clustering
coverages='60'
evalues_filtering='1e-60' # 1e-50 1e-20' #'1e-140 1e-160'
inflations='1.8'
split_cluster="false"

# multiple alignment anlysis
WINDOW="25" # used to group cleavage sites

#Conflict between domain annotation and cleavage site annotation
threshold_overlap_prct=10
threshold_overlap_aa=15
ignoring_threshold_ratio=0.5
