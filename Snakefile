rule create_taxonomy_file:
    input:
        ncbi_current="/mirror/ncbi/current/"
    output:
        taxonomy_file="data/taxonomy/taxonomy_virus.txt",
        taxon_id_gencode_file="data/taxonomy/taxon_id_gencode.txt",
        alternative_taxon_id="data/taxonomy/heterogeneous_taxon_id_taxonomy_virus.txt"
    shell:
        "bash scripts/create_taxonomy_file.sh {input.ncbi_current} true"

## Add as parameter taxon ...
rule extract_viral_protein:
    input:
        taxonomy_file="data/taxonomy/taxonomy_virus.txt",
    output:
        protein_db="data/viral_proteins/Viruses_protein_db.faa",
        protein_stat_file="results/stat_viral_protein/stat_proteins_Viruses.csv"
    shell:
        "bash scripts/extraction_viral_protein.sh"
