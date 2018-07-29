import taxonomy as tax
import viral_genome_classes as obj
import parser_interpro_results as do
import visualisation_genome
import visualisation_protein as view_prot
import viruses_statistics as fct_stat

import os, gzip, logging, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, csv, re
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation
from operator import attrgetter

def getPositionsInAln(cds, list_obj):
    positions = []
    for obj in list_obj:
        start =cds.aln_list.index(obj.start_aa(cds)-1)
        end = cds.aln_list.index(obj.end_aa(cds)-1)
        positions.append((start,end))

    return positions


def visualisation(list_cds, aln_file, line_size):
    display_dico = {}
    for cds in list_cds:

        ##DOMAINS
        print("do")
        position_domains = getPositionsInAln(cds, cds.matchs)
        print("cs")
        position_sites = getPositionsInAln(cds, cds.cleavage_sites)
        print("cs p")
        position_predicted_sites = getPositionsInAln(cds, cds.predicted_cleavage_sites)

        sequence_list = addColorToSequence2(position_domains, position_sites, position_predicted_sites, cds.aligned_sequence)

        lines = []
        print(' line_size', line_size)
        # cut the sequence liste in piece of the line size and transform this pieces in string
        for i_line in range(0, len(sequence_list), line_size):
            line = "".join(sequence_list[i_line : i_line+line_size])
            lines.append(line)
            print(len(sequence_list[i_line : i_line+line_size]))


        for l in lines:
            print(l)
        key=cds.protein_id
        #
        display_dico[key] = lines
        # alignment dico is used only for identity and file header. not efficient at all..
        print(view_prot.visualisation_protein(cds, cds.segment, 1))

    taxon_ids, alignement_dico = store_alignement_line(aln_file, line_size)
    display_alignement(alignement_dico, display_dico)



def visualisation_old(gb_file, genetic_code, gff_file, alignement_dico, sp_treshold, taxon_id):

    genome = obj.gb_file_parser(gb_file, taxon_id, sp_treshold)

    do.getMatchObject(genome, gff_file)
    do.associateMatchWithPolyprotein(genome)
    #visualisation_genome.genomeVisualisation(genome, 1, genetic_code)

    display_dico = {}
    display_header = []
    for i, segment in enumerate(genome.segments):
        for cds in segment.cds:
            if cds.protein_id in alignement_dico:
                # print(cds, cds.matchs)
                                # print(starts, ends)
                sequence_lines = alignement_dico[cds.protein_id]
                sequence = '+'.join(sequence_lines)

                ##DOMAINS
                starts, ends = getAnnotationBorders(cds.matchs)
                starts_in_alignment, ends_in_alignment = getPositionInAlignment(sequence, starts, ends)

                sequence_clean = ''.join(sequence_lines)
                starts_in_alignment2, ends_in_alignment2 = getPositionInAlignment(sequence_clean, starts, ends)
                position_domains = [(starts_in_alignment2[i], ends_in_alignment2[i]) for i in range(len(starts_in_alignment2))]

                ##CLAVAGE SITES
                positon_cleavage_sites = []
                for s in cds.cleavage_sites:
                    positon_cleavage_sites.append(s.start_aa(cds))
                    positon_cleavage_sites.append(s.end_aa(cds))

                position_predicted_cleavage_sites = []
                cleavage_sites_in_aln, predicted_cleavage_sites_in_aln = getPositionInAlignment(sequence_clean, positon_cleavage_sites, position_predicted_cleavage_sites)
                print(cleavage_sites_in_aln)



                # lines2 = addColorToSequence2(position_domains, cleavage_sites_in_aln, position_predicted_cleavage_sites, sequence_lines)

                sequence = addColorToSequence(starts_in_alignment, ends_in_alignment, sequence)
                lines1 = sequence.split('+')

                lines = lines1


                key=genome.taxon_id+'|'+cds.protein_id
                display_dico[key] = lines

    return display_dico


def addColorToSequence2(position_domains, positon_cleavage_sites, position_predicted_cleavage_sites, aln_sequence):
    # Very unelegant to many loop for nothing but only smallvisualisation to be sure
    end_color = '\033[0m'
    domain_color = '\033[94m'
    cleavage_site_color = '\033[91m'
    predicted_site_color = '\033[93m'

    sequence_list = []

    for i, letter in enumerate(aln_sequence):

        if any(( True for start, end in positon_cleavage_sites if start <= i <= end)):
            sequence_list.append( cleavage_site_color + letter + end_color)

        elif any(( True for start, end in position_predicted_cleavage_sites if start <= i <= end)):
            sequence_list.append(predicted_site_color + letter + end_color)

        elif any(( True for start, end in position_domains if start <= i <= end)):
            sequence_list.append(domain_color + letter + end_color)

        else:
            sequence_list.append(letter)

    # print(sequence)
    return sequence_list


def addColorToSequence(starts_in_alignment, ends_in_alignment, sequence):
    end_color = '\033[0m'
    domain_color = '\033[94m'
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
    return sequence


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


def store_alignement_line(alignement_file, aln_row_len):
    #Not very elegant to load all the line
    alignement_dico = {}
    taxon_ids = set()


    with open(alignement_file, "r") as handle:
        alignement_index = None

        pattern = re.compile("(\d+)\|([^|]+)\|[\d]+\s+(.+)")
        flag_bloc = False

        file_header = next(handle)

        for l in handle:
            # print(l)
            result = pattern.search(l)
            if result:
                flag_bloc = True
                protein_id = result.group(2)
                sequence = result.group(3)
                # print(sequence)
                alignement_dico.setdefault(protein_id, []).append(sequence)

                taxon_ids.add(result.group(1))
                # print('sequence', '|'+sequence+'|')
                if not alignement_index:
                    alignement_index = l.index(sequence)


            elif flag_bloc:
                flag_bloc = False
                identity = l[alignement_index:].rstrip()
                identity += ' '*(len(sequence)-len(identity))
                alignement_dico.setdefault('identity', []).append(identity)
                # print('identity', '|'+identity+'|')

        # print(alignement_dico)
        ## Display alignement with the custom aln length per row
        alignement_dico = {id:''.join(seq) for id,seq in alignement_dico.items()}
        for id,seq in alignement_dico.items():
            seq = ''.join(seq)
            # print('len seq', id, len(seq))
            alignement_dico[id] = [seq[i:i+aln_row_len] for i in range(0, len(seq), aln_row_len)]

        alignement_dico['file_header']  = file_header
        # input()
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
    # alignement_file = sys.argv[1]
    try:
        alignement_file = sys.argv[1]
    except IndexError:
        alignement_file = "data/alignment/Viruses_evalue_1e-20coverage50_I1_8/seq_cluster40.aln"

    try:
        minimum_nb_peptide = int(sys.argv[3])
    except IndexError:
        minimum_nb_peptide = 0
    sp_treshold=90
    aln_row_len=150
    taxonomy_file="data/taxonomy/taxonomy_virus.txt"
    expected_file =  "data/taxonomy/polyprotein_expectation_by_taxon.csv"

    gff_file = 'data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'

    # taxon_expectation = tax.expectedPeptide(expected_file)

    taxon_ids, alignement_dico = store_alignement_line(alignement_file, aln_row_len)
    display_dico = {}
    print(taxon_ids)
    for taxon in taxon_ids:
        # file_handle = open(os.path.join(output_dir,'cleavage_site_{}_window_{}.faa'.format(taxon, window_step_clavage_site*2)), "w")

        gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)


        for i, gb_dico in enumerate(gbff_iter):
            # print(gb_dico)
            # print(gb_dico)
            gb_file = gb_dico['gb_file']
            genetic_code = gb_dico['genetic_code']
            taxon_id = gb_dico['taxon_id']
            # print(gb_file)
            print(taxon)
            # print('genetic code', genetic_code)
            # print(gb_file)

            # visualisation_genome.visualisation_main(gb_file, genetic_code, gff_file, 1, 0, taxon_id, sp_treshold)
            display_dico_i = visualisation_old(gb_file, genetic_code, gff_file, alignement_dico, sp_treshold,taxon_id  )
            display_dico.update(display_dico_i)


        # print(i+1, 'Genome analysed from taxon', taxon)
    display_alignement(alignement_dico, display_dico)
    # input()
