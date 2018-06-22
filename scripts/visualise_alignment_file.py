import taxonomy as tax
import object_analysis as obj

import os, gzip, logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, csv, re
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation
from operator import attrgetter



def visualisation(gb_file, genetic_code, gff_file, taxon_expectation, alignement_dico):

    genome = obj.Genome( gb_file)

    with gzip.open(gb_file, "rt") as handle:
    # with open("/home/user/mainguy/Documents/Data_Analysis/GCF_000885175.1_ViralMultiSegProj39867_genomic_MODIFIED.gbff", "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "genbank")):


            segment = obj.Segment(record, gb_file)
            genome.segments.append(segment)

            segment.getMatpeptidesAndPolyproteins()
            segment.checkPeptideRedundancy() #remove the redundant peptide
            segment.checkSubPeptides()
            segment.associatePepWithProt()

            segment.checkForSlippage()
            segment.identifySubProtein()

            segment.getCleavageSites()

    genome.getTaxonExpectation(taxon_expectation)
    genome.identifyExpectedElement()
    genome.getMatchObject(gff_file)
    genome.associateMatchWithPolyprotein()
    genome.visualisation(1, genetic_code)

    display_dico = {}
    display_header = []
    for i, segment in enumerate(genome.segments):
        for cds in segment.cds:
            if cds.protein_id in alignement_dico:
                starts, ends = getAnnotationBorders(cds.matchs)
                sequence = alignement_dico[cds.protein_id]
                sequence = '+'.join(sequence)

                starts_in_alignment, ends_in_alignment = getPositionInAlignment(sequence, starts, ends)
                number = '{}_{}'.format(i+1, cds.polyprotein_number)
                key = genome.expectation_node+'|'+genome.taxon_id+"|"+number
                display_dico[key] = addColorToSequence(starts_in_alignment, ends_in_alignment, sequence)

    return display_dico

def addColorToSequence(starts_in_alignment, ends_in_alignment, sequence):
    end_color = '\033[0m'
    domain_color = '\033[92m'
    cleavage_site_color = '\033[91m'

    pattern_cleavage_site = re.compile(r'([^\d])([a-z]+)')
    pattern_newline = re.compile(r'(\+)')


    for i in range(len(starts_in_alignment))[::-1]:
        start = starts_in_alignment[i]
        end = ends_in_alignment[i]
        domain_sequence = domain_color + sequence[start:end+1] + end_color
        domain_sequence = pattern_newline.sub(r"{}\1{}".format( end_color, domain_color), domain_sequence)

        sequence = sequence[:start] + domain_sequence + sequence[end+1:]

        #Deal with the new line caracther to not propagate color when the seq will be display
         # = pattern_newline.sub(r"{}\1{}".format( end_color, domain_color), sequence[start-len(domain_color):end+1+len(end_color)])
        # print([domain_sequence])
        domain_sequence = pattern_cleavage_site.sub(r"\1\2{}".format(domain_color),  sequence[start:end+1])
        sequence = sequence[:start] + domain_sequence + sequence[end+1:]


    sequence = pattern_cleavage_site.sub(r"\1{}\2{}".format(cleavage_site_color, end_color), sequence)
    # print([sequence])
    return sequence.split('+')


def getPositionInAlignment(sequence, starts, ends):
    compteur = 0
    starts_in_alignment = []
    ends_in_alignment = []

    for i, a in enumerate(sequence):
        if a not in  ['-', '+']:
            compteur += 1
            if compteur in starts:
                starts_in_alignment.append(i)
            elif compteur in ends:
                ends_in_alignment.append(i)
    # print(starts)
    # print(starts_in_alignment)
    # print(ends)
    # print(ends_in_alignment)
    return starts_in_alignment, ends_in_alignment


def getAnnotationBorders(domains):
    starts = []
    ends = []
    current_start = 0
    current_end = 0
    for d in sorted(domains, key=lambda x: x.start, reverse=False):
        if current_start <= d.start_in_prot <= current_end:
            current_end = max(d.end_in_prot, current_end)
        else:
            starts.append(d.start_in_prot)
            ends.append(d.end_in_prot)
            current_start = d.start_in_prot
            current_end = d.end_in_prot

    return starts, ends


def store_alignement_line(alignement_file):
    #Not very elegant to load all the line
    alignement_dico = {}
    taxon_ids = set()


    with open(alignement_file, "r") as handle:
        alignement_index = None

        pattern = re.compile("(\w+)\|(\d+)\|([^|]+)\|([_.\d]+)\s+(.+)")
        flag_bloc = False

        alignement_dico['file_header'] = next(handle)

        for l in handle:
            result = pattern.search(l)
            if result:
                flag_bloc = True
                protein_id = result.group(3)
                sequence = result.group(5)
                # print(sequence)
                alignement_dico.setdefault(protein_id, []).append(sequence)

                taxon_ids.add(result.group(2))
                if not alignement_index:
                    alignement_index = l.index(sequence)

            elif flag_bloc:
                flag_bloc = False
                identity = l[alignement_index:].rstrip()
                alignement_dico.setdefault('identity', []).append(identity)


        return taxon_ids, alignement_dico


def display_alignement(alignement_dico, display_dico):

    max_len = max([len(h) for h in display_dico]) + 6

    print(alignement_dico['file_header'])

    for i in range(len(alignement_dico['identity'])):
        for k, v in display_dico.items():
            print(k+' '*(max_len-len(k))+v[i])
        print(' '*(max_len)+alignement_dico['identity'][i])
        print()
if __name__ == '__main__':
    # logging.basicConfig(filename='log/genbankparser.log',level=logging.INFO)

    try:
        alignement_file = sys.argv[1]
    except IndexError:
        alignement_file = "data/alignement/Comovirus_polyprotein_1_1.0.aln"

    try:
        minimum_nb_peptide = int(sys.argv[3])
    except IndexError:
        minimum_nb_peptide = 0

    taxonomy_file="data/taxonomy/taxonomy_virus.txt"
    expected_file =  "data/taxonomy/polyprotein_expectation_by_taxon.csv"

    output_dir = '/home/user/mainguy/Documents/Data_Analysis/data/Cleavage_site_sequences'
    gff_file = '/scratch/polyproteins_interpro/interproscan_result/All_polyprotein.gff3'

    taxon_expectation = tax.expectedPeptide(expected_file)

    taxon_ids, alignement_dico = store_alignement_line(alignement_file)
    display_dico = {}
    for taxon in taxon_ids:
        # file_handle = open(os.path.join(output_dir,'cleavage_site_{}_window_{}.faa'.format(taxon, window_step_clavage_site*2)), "w")

        gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)


        for i, gb_dico in enumerate(gbff_iter):
            # print(gb_dico)

            gb_file = gb_dico['gb_file']
            genetic_code = gb_dico['genetic_code']
            # print(gb_file)

            if i%200== 0:
                # continue
                print(i)
            # print('genetic code', genetic_code)
            print(gb_file)
            display_dico_i = visualisation(gb_file, genetic_code, gff_file, taxon_expectation, alignement_dico  )
            display_dico.update(display_dico_i)

        # print(i+1, 'Genome analysed from taxon', taxon)
    display_alignement(alignement_dico, display_dico)
