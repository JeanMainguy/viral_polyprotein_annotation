import taxonomy as tax
import viral_genome_classes as obj
import parser_interpro_results as do

import os
import logging
import sys
import re
import os.path
import collections

SCREEN_SIZE = 100


def visualisation_main(gb_file, genetic_code, gff_file, nb_line, minimum_nb_peptide, taxon_id, sp_treshold):

    genome = obj.gb_file_parser(gb_file, taxon_id, sp_treshold)

    if gff_file:
        do.getMatchObject(genome, gff_file)
        do.associateDomainsWithPolyprotein(genome)
    # for segment in genome.segments:
    #     do.getDomainOverlappingInfo(segment)

    # genome.getTaxonExpectation(taxon_expectation)
    # genome.identifyExpectedElement()
    # genome.getMatchObject(gff_file)
    # genome.associateMatchWithPolyprotein()
    # taxon_id = segment.taxon_id
    # for segment in genome.segments:
    #     if taxon_id != segment.taxon_id or  segment.taxon_id == '40687':
    #         print(gb_file)
    #         input()

    # genome.visualisation(nb_line, genetic_code)
    # if genome.numberOf("peptides") >= minimum_nb_peptide:
    genomeVisualisation(genome, nb_line, genetic_code)


def genomeVisualisation(genome, nb_line=1, genetic_code=None):
    '''
    Visualisation of the genome in a friendly way
    Example:
    Segment 1:
    ------------------------------------------------------
    ======>   ======================>               ====>
              ==========<-1====================>
    Segment 2:
    ----------------------------------------
    ======>   ======================>
    '''
    visu_string = 'GENOME of {}|{}\n'.format(genome.segments[0].organism, genome.taxon_id)
    longuest_segment_len = max([len(s) for s in genome.segments])

    conversion = longuest_segment_len/SCREEN_SIZE
    for i, segment in enumerate(genome.segments):
        visu_string += 'SEGMENT {}/{}: {} cds | {} polyprotein | {}/{} final mat_peptides\n'.format(i+1,
                                                                                                    len(genome.segments),
                                                                                                    len(segment.cds),
                                                                                                    len([
                                                                                                        cds for cds in segment.cds if cds.polyprotein]),
                                                                                                    len(segment.peptides)-len(
                                                                                                        segment.parent_peptides),
                                                                                                    len(segment.peptides))

        visu_string += f'taxonomy: {segment.record.annotations["taxonomy"]}\n'

        # print('taxonomy:{} \nnode giving the nb of peptides {}'.format(segment.record.annotations['taxonomy'], genome.expectation_node))
        # print('Expected number of mature peptide {} '.format(genome.peptide_expectation))
        # for p in segment.cds:
        #     print('Is annotation relevant? ',p.polyprotein_number, p.isAnnotationRelevant())
        visu_string += visualisation(segment, nb_line, genetic_code)
        # print(features_strg)
    print(visu_string)

# visualisation genome


def visualisation(self, nb_line, genetic_code):

    compatible_dico = collections.OrderedDict()
    compatible_dico['match'] = buildCompatibleGroup(self.matchs)[::-1]
    compatible_dico['cds'] = buildCompatibleGroup(self.cds)
    compatible_dico['pep'] = buildCompatibleGroup(self.peptides)
    compatible_dico['unannotated_region'] = buildCompatibleGroup(self.unannotated_region)

    return get_final_strings(self, nb_line, genetic_code, compatible_dico)


def getStringIndices(seq, final_seq_to_display, conversion):
    # The whole segment is visualise

    if final_seq_to_display.__class__.__name__ == 'Segment':
        return (int(seq.start/conversion), int((seq.end-1)/conversion))

    # The protein is visualise. protein reference coordinate are then returned
    if seq == final_seq_to_display:
        return (0, int((len(seq)/3)/conversion)-1)
    else:
        return (int(seq.start_aa(final_seq_to_display)/conversion), int((seq.end_aa(final_seq_to_display)-1)/conversion))


def get_final_strings(final_seq_to_display, nb_line, genetic_code, compatible_dico, size_sequence=SCREEN_SIZE):
    """
    final_seq_to_display may be a segment object or a Protein object
    the getStringIndices fct would then react accordingly
    """
    try:
        display_len = int(nb_line) * int(size_sequence)
        if final_seq_to_display.__class__.__name__ == 'Segment':
            conversion = int(len(final_seq_to_display)) / display_len
        elif final_seq_to_display.__class__.__name__ == 'Protein':
            conversion = int(len(final_seq_to_display)/3) / display_len
    except ValueError:
        display_len = int(len(final_seq_to_display) /
                          3) if nb_line == 'aa' else int(len(final_seq_to_display))
        conversion = 3 if nb_line == 'aa' else 1

    strings = []
    seq_type_dico = {
        'match': ('[', '#', ']', ' '),
        'cds': ('=', '=', '>', '-'),
        'pep': ('|', '+', '|', ' '),
        'unannotated_region': ('|', '~', '|', ' ')
    }
    color_end = '\033[0m'
    overlap_col = '\033[91m'
    # color = '$'
    # color_end= '$'
    for seq_type, compatible_group in compatible_dico.items():
        # different symbol if it protein or peptide
        left, central, right, neutral = seq_type_dico[seq_type]
        for group in compatible_group:
            string = neutral*display_len

            for seq in group:
                i_start, i_end = getStringIndices(seq, final_seq_to_display, conversion)
                string = string[:i_start] + left+central * \
                    (i_end - i_start-1) + right + string[i_end+1:]

                if seq.__class__.__name__ == 'Protein' and final_seq_to_display.__class__.__name__ == 'Segment':
                    string = format_protein_line(
                        seq, string, genetic_code, i_start, i_end, final_seq_to_display.record, nb_line)
                else:
                    name = seq.getNameToDisplay()
                    string = string[:i_start+2] + name + string[i_start+2+len(name):]
            strings.append(string)

    # REGEX
    # overlap_match re.compile(r"\[#<"g)
    # match = re.compile(r"\[#[^<]"g)
    # sub_string = re.sub(r"(\[#<[A-Za-z0-9]+#+\]?)", r"{}\1{}".format(color, color_end), sub_string)
    visu_string = ''
    unfinished_col_list = ["" for s in strings]
    for i in range(0, display_len, SCREEN_SIZE):
        visu_string += '/'*SCREEN_SIZE
        visu_string += '\n\n'
        for i_line, string in enumerate(strings):
            sub_string = string[i:i+SCREEN_SIZE]

            if unfinished_col_list[i_line]:
                sub_string = re.sub(
                    r"^([^ \[]+)", r"{}\1{}".format(unfinished_col_list[i_line], color_end), sub_string)
            sub_string = re.sub(r"(\[?#?<[A-Za-z0-9]+#*\]?)",
                                r"{}\1{}".format(overlap_col, color_end), sub_string)

            # -1 not not catch the end color character of the end of the line
            ansi_list = re.findall(r"(\x1B\[[0-?]*[ -/]*[@-~])", sub_string[:-1])

            unfinished_col_list[i_line] = '' if not ansi_list or ansi_list[-1] == color_end else ansi_list[-1]
            # if ansi_list and ansi_list[-1] != color_end:
            #     unfinished_col_dict[i_line] = ansi_list[-1]
            visu_string += f'  {sub_string}\n'

    # protein = '\n'.join(strings)

    visu_string += color_end
    # print(protein)
    # return 'protein'
    return visu_string


def format_protein_line(seq, string, genetic_code, i_start, i_end, record, nb_line):
    # print(seq.getSequenceRecord(self.organism, self.taxon_id, self.record, genetic_code).seq)
    if nb_line == 'aa':
        seqaa = seq.getSequenceAA(record, genetic_code)

        string = string[:i_start] + seqaa + string[i_start+len(seqaa):]
        return string
    if nb_line == 'nt':
        seqaa = seq.getSequenceAA(record, genetic_code)
        seqaa = ''.join([a*3 for a in seqaa])
        string = string[:i_start] + seqaa + string[i_start+len(seqaa):]
        return string

    # display positional number
    # if seq.polyprotein_number:
    #     nb = str(seq.polyprotein_number)
    #     string = string[:i_start+2] +nb+ string[i_start+2+len(nb):]
    #
    # Display protein number
    nb = str(seq.number)
    string = string[:i_start+2] + nb + string[i_start+2+len(nb):]

    if False and seq.ribosomal_slippage:
        for position, shift in seq.ribosomal_slippage.items():
            string_postion = round(position/conversion)
            string_postion = string_postion if string_postion < i_end else string_postion-1
            shift_string = '+' if shift > 0 else "-"
            shift_string += str(abs(shift))

            string = string[:string_postion-len(shift_string)+1] + \
                shift_string + string[string_postion+1:]
    return string


def buildCompatibleGroup(set_of_sequence):
    # sorted(list(set_of_sequence), key=lambda x: x.bp_obj.location.start, reverse=False)
    sequences = set_of_sequence.copy()
    # Build the set of non overlaping protein

    for i, seq in enumerate(sequences):

        for seq_next in list(sequences)[i+1:]:

            if not seq.overlap(seq_next):
                seq.non_overlapping_prot.add(seq_next)
                seq_next.non_overlapping_prot.add(seq)

    # Build list of sequences that don't overlap
    compatible_groupes = []
    while len(sequences) > 0:
        seq = sequences.pop()
        group = [seq]
        # finding the sequence that are compatible with seq and within each other
        # set of sequence that are compatible with the sequence stored in the list group
        group_compatible = seq.non_overlapping_prot
        while group_compatible:
            try:
                # extarct one seq that is compatibe with group
                compatible_seq = (group_compatible & sequences).pop()
            except KeyError:
                break
            # build new set 'group_compatible' that store compatible set
            group_compatible = group_compatible & compatible_seq.non_overlapping_prot & sequences

            group.append(compatible_seq)
            sequences.remove(compatible_seq)

        compatible_groupes.append(sorted(group, key=lambda x: x.end, reverse=True))
    return compatible_groupes


if __name__ == '__main__':
    # logging.basicConfig(filename='log/genbankparser.log',level=logging.INFO)
    try:
        nb_line = sys.argv[2]
    except IndexError:
        nb_line = 1

    try:
        minimum_nb_peptide = int(sys.argv[3])
    except IndexError:
        minimum_nb_peptide = 0

    sp_treshold = 90
    taxonomy_file = "results_db_viral_2018-10-19/genomes_index/taxonomy_virus.txt"

    gff_file = 'results_db_viral_2018-10-19/intermediate_files/interproscan_results/domains_viral_sequences.gff3'

    if not os.path.isfile(gff_file):
        logging.warning('///No interproscan result file found///')
        gff_file = None

    try:
        taxon = sys.argv[1]
    except IndexError:
        # taxon = 'ssRNA viruses'
        taxon = "Togaviridae"
        # taxon = 'Rubivirus'
        # taxon = "Marafivirus"
        # # taxon= "Togaviridae"
        # taxon='ssRNA positive-strand viruses, no DNA stage'
        # taxon='Alphavirus'
        # taxon="11036"

    try:
        excluded_taxon = sys.argv[4]
    except IndexError:
        excluded_taxon = None

    # taxon_expectation = tax.expectedPeptide(expected_file)

    # file_handle = open(os.path.join(output_dir,'cleavage_site_{}_window_{}.faa'.format(taxon, window_step_clavage_site*2)), "w")

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file, excluded_taxon)

    print('*-'*100)
    print('VISUALISATION OF THE RefSEq GENOME FROM THE TAXON {} THAT HAVE AT LEAST {} ANNOTATED PEPTIDE'.format(
        taxon, minimum_nb_peptide))
    print('*-'*100)
    for i, gb_dico in enumerate(gbff_iter):
        # print(gb_dico)

        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        taxon_id = gb_dico['taxon_id']
        print(gb_file)

        if i % 200 == 0:
            # continue
            print(i)
        # print('genetic code', genetic_code)
        # print(gb_file)
        visualisation_main(gb_file, genetic_code, gff_file, nb_line,
                           minimum_nb_peptide, taxon_id, sp_treshold)
        input()

    print(i+1, 'Genome analysed from taxon', taxon)
