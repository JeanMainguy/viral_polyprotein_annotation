################################################################################
## INPUT GENOMES DATA BASE
################################################################################

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

genetic_code_path="genome_db_test/taxonomy/new_taxdump/"

genbank_files_db_path="genome_db_test/genomes/refseq/viral/"

RefSeq_structure="True"

################################################################################
## VARIABLES
################################################################################

# Taxon : analysis will be performed on every genomes of the database belonging to this taxon
taxon='Flaviviridae'

## Signal Peptide length threshold in aa
tresholdSP="90"

## blast evalue start
blast_evalue='1e-5'

## filtering and mcl clustering
# Can be a list of value separated by a space
coverages='60' # minimal coverage percentage of the longest sequence between a pair of sequences
evalues_filtering='1e-60' #
inflations='1.8'

# Split cluster that have framshiffted proteins
split_cluster="false"

# multiple alignment anlysis
window="25" # used to group cleavage sites
confidence_score_treshold=4 # confidence score threshold for the group of cleavage sites

## IRRELEVANT ANNOTATION IDENTIFICATION:
# Conflict between domain annotation and cleavage site annotation
interpro_path="~/path/to_interproscan/folder/"
interpro_db="PFAM,CDD,ProDom,SMART,ProSiteProfiles"

threshold_overlap_prct=10
threshold_overlap_aa=15
ignoring_threshold_ratio=0.5

#Create a black list of poorly annotated cds
black_list_conflict_domain='true'
black_list_irrelevant_pattern='true'
