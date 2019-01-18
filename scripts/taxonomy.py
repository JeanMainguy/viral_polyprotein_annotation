import os
import logging
import gzip
import csv
from Bio import SeqIO
import sys


# logging.basicConfig(filename='refseq_genome_path.log')
def see_objet(obj):
    print('SEE OBJET ', type(obj))

    for attr in dir(obj):
        if attr[0] != '_' and attr != "seq":
            print(attr, " ", getattr(obj, attr))


def get_gb_file_from_RefSeq_db(genbank_file_db):
    # genbank_file_db = "/mirror/ncbi/current/genomes/refseq/viral/"
    # virus_name = "Ageratum_yellow_vein_China_virus"
    for virus_name in os.listdir(genbank_file_db):
        virus_path = os.path.join(genbank_file_db, virus_name)

        assembly_path = os.path.join(virus_path, 'latest_assembly_versions')
        # gb_files = []
        # Checking if the path exist..
        if not os.path.isdir(virus_path):
            logging.info('Cannot find the virus folder: '+virus_path)
            if not virus_name.endswith('.txt'):
                logging.warning('Cannot find the virus folder and not a txt file: '+virus_path)
            # raise Exception(virus_path, 'the virus folder does not exist..')
            # return []

        if not os.path.isdir(assembly_path):

            all_assembly_dir = os.path.join(virus_path, 'all_assembly_versions')
            all_assembly_content = os.listdir(all_assembly_dir)
            if len(all_assembly_content) == 1 and all_assembly_content[0] == 'suppressed':
                logging.info('No latest_assembly_versions folder for {}'.format(virus_name))
                logging.info('only suppressed dir in all assembly: {} | we can ignore this virus : {}'.format(
                    all_assembly_content, virus_name))
            else:
                logging.warning('No latest_assembly_versions folder for {}'.format(virus_name))
                logging.warning('NOT only suppressed dir in all assembly.. {} | This virus has a problem.. :'.format(
                    all_assembly_content, virus_name))
            continue

            # raise Exception(assembly_path, 'the virus folder does not have "latest_assembly_versions" directory...')
            # return []
        if len(os.listdir(assembly_path)) == 0:
            logging.warning('the path {} is empty'.format(assembly_path))

        for assembly in os.listdir(assembly_path):
            gb_file = os.path.join(assembly_path, assembly, assembly+'_genomic.gbff.gz')
            # print(gb_file)
            if os.path.exists(gb_file):
                yield gb_file
                # gb_files.append(gb_file)
            else:
                logging.warning('No gbff file have been found in the directory: '+gb_file)
                # raise Exception(assembly_path, 'No gbff file have been found in the directory')
        # return gb_files


def getTaxonomy(gb_file, error_taxon_ids):
    # polyprot = []
    # matpep = []
    # poly_di = {} # dict with poly as key and list of matpep as value

    # print('='*20)
    taxonomy_set = set()
    taxon_set = set()
    proper_open = gzip.open if gb_file.endswith('.gz') else open

    with proper_open(gb_file, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            taxonomy = record.annotations['taxonomy']
            taxonomy_set.add(";".join(taxonomy))
            organism = record.annotations['organism']
            # print(record.annotations)
            for f in record.features:
                if f.type == 'source':
                    taxon = False
                    for item in f.qualifiers['db_xref']:
                        if item.startswith('taxon'):
                            taxon = item.replace('taxon:', '')
                            taxon_set.add(taxon)
                    if not taxon:
                        logging.warning('No taxon id found for {}'.format(organism, gb_file))

                    # else: #ONLY check the first source feature | EDIT no because sometime more than one source ..
                    #     break
            # see_objet(record)
        if len(taxon_set) > 1:

            error_taxon_ids.append(';'.join(taxon_set)+'\t'+gb_file+'\n')
            taxon_sort = sorted([int(t) for t in taxon_set])
            taxon_set = {taxon_sort[0]}  # select the first taxon id from the list
            logging.warning(
                f'Taxon id is not homogenous ({taxon_sort}) in {gb_file}. Ids stored in error_taxon_ids file and the first id {taxon_sort[0]} is used')
        if len(taxonomy_set) > 1:
            sort_taxonomys = sorted(list(taxonomy_set))
            logging.warning(
                f'Taxonomy is not homogenous ({sort_taxonomys}) in {gb_file}. the first taxonomy is used: {sort_taxonomys[0]}')

        return taxonomy_set.pop(), organism, taxon_set.pop()


def createTaxonomyFile(taxonomy_file, genbank_file_db, alternative_taxon_id_file, refseq_structure):
    # """ create csv file from the genome refseq db with info on avaible genome """
    viral_taxons = {}
    error_taxon_ids = []
    # taxonomy_file = os.path.join(output_dir, 'taxonomy_virus.txt')
    nb_genome = len(os.listdir(genbank_file_db))
    last_percent = 0
    duplicated_taxon_ids = {}

    if refseq_structure:
        gb_files = get_gb_file_from_RefSeq_db(genbank_file_db)

    else:
        gb_files = (os.path.join(genbank_file_db, file) for file in os.listdir(genbank_file_db)
                    if file.endswith('.gb') or file.endswith('.gbff'))

    with open(taxonomy_file, 'w') as handle:

        for i, gb in enumerate(gb_files):
            taxonomy, organism, taxon = getTaxonomy(gb, error_taxon_ids)
            taxon = int(taxon)
            if taxon in viral_taxons:  # if the taxon id is already found in the dict then..
                duplicated_taxon_ids[taxon] = duplicated_taxon_ids.setdefault(taxon, 0) + 1
                # duplicated_taxon_ids[taxon] += 1
                # duplicated stored with underscore follow by the number of duplication
                duplicated_taxon = f'{taxon}_{duplicated_taxon_ids[taxon]}'
                logging.warning(
                    f'Duplicated taxon id: taxon id {taxon} is found in more than one genbank file.. new_taxon_id used: {duplicated_taxon}')
                taxon = duplicated_taxon
            else:
                viral_taxons[taxon] = None

            handle.write("{}\t{}\t{}\t{}\n".format(taxon, organism, taxonomy, gb))

            if (i/nb_genome)*100 > last_percent + 10 and (i/nb_genome)*100 <= 100:
                print(round((i/nb_genome)*100), "% processed")
                last_percent = (i/nb_genome)*100

    with open(alternative_taxon_id_file, 'w') as handle:
        for l in error_taxon_ids:
            handle.write(l)

    # print(round((i/nb_genome)*100), "% processed")
    return viral_taxons


def getAllRefseqFromTaxon(wanted_taxonomy, taxonomy_file, excluded_taxon=None):
    with open(taxonomy_file, 'r') as taxfl:

        for l in taxfl:
            # print(l)
            (tax_id, organism, taxonomy, genetic_code, gbff) = l.split("\t")
            # print(taxonomy)
            if wanted_taxonomy in taxonomy.split(';') and excluded_taxon not in taxonomy.split(';'):
                # print('taxonomy')
                yield {'gb_file': gbff.rstrip(), 'genetic_code': genetic_code, 'taxon_id': tax_id, 'taxonomy': taxonomy}

            if tax_id == wanted_taxonomy or organism == wanted_taxonomy:
                yield {'gb_file': gbff.rstrip(), 'genetic_code': genetic_code, 'taxon_id': tax_id, 'taxonomy': taxonomy}
                # no break here because genome with the same tax id... :-/


def getAllRefseqFromTaxonIdList(wanted_taxon_ids, taxonomy_file, excluded_taxon=None):
    with open(taxonomy_file, 'r') as taxfl:

        for l in taxfl:
            # print(l)
            (tax_id, organism, taxonomy, genetic_code, gbff) = l.split("\t")

            if tax_id in wanted_taxon_ids:
                yield {'gb_file': gbff.rstrip(), 'genetic_code': genetic_code, 'taxon_id': tax_id, 'taxonomy': taxonomy}
                # no break here because genome with the same tax id... :-/


def expectedPeptide(expected_pep_file):

    taxon_expectation = {}
    with open(expected_pep_file, "r", newline='') as fl:
        reader = csv.DictReader(fl, delimiter='\t')
        for row in reader:
            row = dict(row)
            # conversion in int
            for key in row:
                if key == 'taxon':
                    taxon = row[key]
                    continue
                row[key] = None if not row[key] else int(row[key])

            if row["peptide"] > 0:
                taxon_expectation[taxon] = row
            # print(row['first_name'], row['last_name'])
        #
        # for l in fl:
        #     line_elements = l.split("\t")
        #     expected_peptide = int(expected_peptide)
        #     polyprotein = None if not polyprotein else int(polyprotein)
        #     if expected_peptide > 0:
        #         taxon_expectation[taxon] = {"polyprotein":polyprotein, "peptide":expected_peptide}
    return taxon_expectation


def getGeneticCode(viral_taxons, geneticcode_file, tmp_output_file, output_file):
    count = 0
    genetfl = open(geneticcode_file, 'r')
    for line in genetfl:
        taxon, genetic_code = line.split('\t')
        # print(taxon, genetic_code)
        if int(taxon) in viral_taxons:
            viral_taxons[int(taxon)] = int(genetic_code.strip())
            count += 1

    print('count', count)
    print('viral taxons', len(viral_taxons))
    tmpfl = open(tmp_output_file, 'r')
    outfl = open(output_file, 'w')

    for l in tmpfl:
        line_split = l.split('\t')
        taxon = line_split[0]
        # in case of a duplicated taxon id, the raw_taxon will be written with an underscore and an int.
        # We want only the taxon id here
        if '_' in taxon:
            taxon = taxon[:taxon.find('_')]

        taxon = int(taxon)
        try:
            genetic_code = int(viral_taxons[taxon])
        # if taxon has not been found in genetic code file then the value of the dict is None and int(None) through an error
        except TypeError:
            logging.warning(
                'The genetic code of the taxon {} was not found. By default genetic code 1 is provided'.format(taxon))
            genetic_code = "1"
        newline = line_split[:-1] + [str(genetic_code), line_split[-1]]
        outfl.write('\t'.join(newline))

    outfl.close()
    tmpfl.close()
    genetfl.close()


if __name__ == '__main__':
    print('#Creation of taxonomy file')

    output_file = sys.argv[1]
    genbank_file_db = sys.argv[2]
    geneticcode_file = sys.argv[3]
    alternative_taxon_id_file = sys.argv[4]
    try:
        refseq_structure = False if sys.argv[5].upper().startswith('F') else True
    except IndexError:
        print('RefSeq structure bool info not provided.')
        print('By default genbank files are considered stored in a RefSeq structure')
        refseq_structure = True

    # print(tmp_output_file)
    path, file_name = os.path.split(output_file)
    tmp_output_file = os.path.join(path, 'tmp_file.tmp')
    logging.basicConfig(level=logging.INFO)

    viral_taxons = createTaxonomyFile(
        tmp_output_file, genbank_file_db, alternative_taxon_id_file, refseq_structure)
    getGeneticCode(viral_taxons, geneticcode_file, tmp_output_file, output_file)

    print('END of taxonomy.py')
