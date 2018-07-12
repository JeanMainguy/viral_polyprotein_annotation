import os
import logging
import re
import gzip, csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
from Bio.Alphabet import generic_protein


# logging.basicConfig(filename='refseq_genome_path.log')
def see_objet(obj):
    print('SEE OBJET ', type(obj))

    for attr in dir(obj):
        if attr[0] != '_' and attr != "seq":
            print(attr, " ", getattr(obj, attr))


def getGbffFile(virus_name, ncbi_refseq_db):
    # ncbi_refseq_db = "/mirror/ncbi/current/genomes/refseq/viral/"
    # virus_name = "Ageratum_yellow_vein_China_virus"

    virus_path = os.path.join(ncbi_refseq_db, virus_name)

    assembly_path = os.path.join(virus_path, 'latest_assembly_versions')
    gb_files = []
    # Checking if the path exist..
    if not os.path.isdir(virus_path):
        logging.info('Cannot find the virus folder: '+virus_path)
        if not virus_name.endswith('.txt'):
            logging.warning('Cannot find the virus folder and not a txt file: '+virus_path)
        # raise Exception(virus_path, 'the virus folder does not exist..')
        return []

    if not os.path.isdir(assembly_path):

        all_assembly_dir = os.path.join(virus_path, 'all_assembly_versions')
        all_assembly_content = os.listdir(all_assembly_dir)
        if len(all_assembly_content) == 1 and all_assembly_content[0] == 'suppressed':
            logging.info('No latest_assembly_versions folder for {}'.format(virus_name))
            logging.info('only suppressed dir in all assembly: {} | we can ignore this virus : {}'.format(all_assembly_content,virus_name ))
        else:
            logging.warning('No latest_assembly_versions folder for {}'.format(virus_name))
            logging.warning('NOT only suppressed dir in all assembly.. {} | This virus has a problem.. :'.format(all_assembly_content,virus_name ))


        # raise Exception(assembly_path, 'the virus folder does not have "latest_assembly_versions" directory...')
        return []
    if len(os.listdir(assembly_path)) == 0:
        logging.warning('the path {} is empty'.format(assembly_path))

    for assembly in os.listdir(assembly_path):
        gb_file = os.path.join(assembly_path, assembly, assembly+'_genomic.gbff.gz')
        # print(gb_file)
        if os.path.exists(gb_file):
            gb_files.append(gb_file)
        else:
            logging.warning('No gbff file have been found in the directory: '+gb_file)
            # raise Exception(assembly_path, 'No gbff file have been found in the directory')
    return gb_files

def getTaxonomy(gb_file, error_taxon_ids):
    # polyprot = []
    # matpep = []
    # poly_di = {} # dict with poly as key and list of matpep as value

    # print('='*20)
    taxonomy_set = set()
    taxon_set = set()
    with gzip.open(gb_file, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            taxonomy = record.annotations['taxonomy']
            taxonomy_set.add(";".join(taxonomy))
            organism = record.annotations['organism']
            # print(record.annotations)
            for f in record.features:
                if f.type == 'source':
                    taxon=False
                    for item in f.qualifiers['db_xref']:
                        if item.startswith('taxon'):
                            taxon = item.replace('taxon:', '')
                            taxon_set.add(taxon)
                    if not taxon:
                        logging.warning('No taxon id found for {}'.format(organism, gb_file))

                    # else: #ONLY check the first source feature no because sometime more than one source ..
                    #     break
            # see_objet(record)
        if len(taxon_set) > 1:

            error_taxon_ids.append(';'.join(taxon_set)+'\t'+gb_file+'\n')
            logging.warning('Taxon id is not homogenous in '+gb_file)
        if len(taxonomy_set) > 1:
            logging.warning('Taxonomy is not homogenous in '+gb_file)

        return taxonomy_set.pop(), organism, int(taxon_set.pop())

def createTaxonomyFile(taxonomy_file, ncbi_refseq_db, alternative_taxon_id_file):
    """
    create csv file from the genome refseq db with info on avaible genome
    """
    viral_taxons = {}
    error_taxon_ids = []
    # taxonomy_file = os.path.join(output_dir, 'taxonomy_virus.txt')
    nb_genome = len(os.listdir(ncbi_refseq_db))
    last_percent = 0
    with open(taxonomy_file, 'w') as handle:
        for i, virus in enumerate(os.listdir(ncbi_refseq_db)):

            gb_files = getGbffFile(virus, ncbi_refseq_db)
            # print(gb_files)
            for gb in gb_files:
                taxonomy, organism, taxon = getTaxonomy(gb, error_taxon_ids)

                handle.write( "{}\t{}\t{}\t{}\n".format(taxon, organism, taxonomy, gb))
                viral_taxons[int(taxon)] = None



            if (i/nb_genome)*100 > last_percent + 10:
                print(round((i/nb_genome)*100), "% processed")
                last_percent = (i/nb_genome)*100

                # return viral_taxons
        # print i


    with open(alternative_taxon_id_file, 'w') as handle:
        for l in error_taxon_ids:
            print(l)
            handle.write(l)

    print(round((i/nb_genome)*100), "% processed")
    return viral_taxons


def getAllRefseqFromTaxon(wanted_taxonomy, taxonomy_file="data/taxonomy/taxonomy_virus.txt", excluded_taxon = None):
    with open(taxonomy_file, 'r') as taxfl:

        for l in taxfl:
            # print(l)
            (tax_id, organism, taxonomy, genetic_code, gbff) = l.split("\t")
            # print(taxonomy)
            if wanted_taxonomy in taxonomy.split(';') and excluded_taxon not in taxonomy.split(';'):
                # print('taxonomy')
                yield {'gb_file':gbff.rstrip(), 'genetic_code': genetic_code, 'taxon_id':tax_id}

            if tax_id == wanted_taxonomy or organism == wanted_taxonomy:
                yield {'gb_file':gbff.rstrip(), 'genetic_code': genetic_code, 'taxon_id':tax_id}
                #no break here because genome with the same tax id... :-/

def expectedPeptide(expected_pep_file):

        taxon_expectation = {}
        with open(expected_pep_file, "r", newline='') as fl:
            reader = csv.DictReader(fl, delimiter='\t')
            for row in reader:
                row = dict(row)
                #conversion in int
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
    tmpfl = open(tmp_output_file)
    outfl = open(output_file, 'w')

    for l in tmpfl:
        line_split = l.split('\t')
        taxon = int(line_split[0])
        if viral_taxons[taxon]:
            genetic_code = str(viral_taxons[taxon])
        else:
            logging.warning('The genetic code of the taxon {} was not found. By default genetic code 1 is provided'.format(taxon))
            genetic_code = "1"
        newline = line_split[:-1] + [genetic_code, line_split[-1]]
        outfl.write('\t'.join(newline))

    outfl.close()
    tmpfl.close()
    genetfl.close()

if __name__ == '__main__':
    print('#Creation of taxonomy file')

    output_file= sys.argv[1]
    ncbi_refseq_db = sys.argv[2]
    geneticcode_file = sys.argv[3]
    alternative_taxon_id_file = sys.argv[4]
    # print(tmp_output_file)
    path, file_name = os.path.split(output_file)
    tmp_output_file = os.path.join(path, 'tmp_file')

    logging.basicConfig(filename='log/taxonomy_file_creation.log', level=logging.WARNING)

    viral_taxons = createTaxonomyFile( tmp_output_file, ncbi_refseq_db, alternative_taxon_id_file)
    getGeneticCode(viral_taxons, geneticcode_file,tmp_output_file, output_file)



    # print(viral_taxons)
    # taxonomy_file = "/home/user/mainguy/Documents/Data_Analysis/data/taxonomy/taxonomy_virus.txt"
    # taxon = "Nidovirales"
    # gb_file_iter = getAllRefseqFromTaxon(taxon, taxonomy_file)
    #
    # for g in gb_file_iter:
    #     print(g)


# print(len(list(gb_file_iter)))
