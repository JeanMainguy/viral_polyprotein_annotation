SCREEN_SIZE = 200
import taxonomy as tax
import viral_genome_classes as obj
import parser_interpro_results as do

import os, gzip, logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, csv, re
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation
from operator import attrgetter
import collections


def visualisation_main(gb_file, genetic_code, gff_file, nb_line, minimum_nb_peptide, taxon_id, sp_treshold):

    genome = obj.gb_file_parser(gb_file, taxon_id, sp_treshold)

    do.getMatchObject(genome, gff_file)
    do.associateMatchWithPolyprotein(genome)

    print(genome.matchs)
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
    print('GENOME of {}|{}'.format(genome.organism, genome.taxon_id))
    longuest_segment_len = max([len(s) for s in genome.segments])

    conversion = longuest_segment_len/SCREEN_SIZE
    for i, segment in enumerate(genome.segments):
        print('SEGMENT {}/{}: {} cds | {} polyprotein | {}/{} final mat_peptides'.format(i+1,
            len(genome.segments),
            len(segment.cds),
            len(segment.polyproteins),
            len(segment.peptides)-len(segment.parent_peptides),
            len(segment.peptides)))
        print('taxonomy:{} \nnode giving the nb of peptides {}'.format(genome.taxonomy, genome.expectation_node))
        print('Expected number of mature peptide {} '.format(genome.peptide_expectation))
        # for p in segment.cds:
        #     print('Is annotation relevant? ',p.polyprotein_number, p.isAnnotationRelevant())
        visualisation(segment, nb_line, genetic_code)
        # print(features_strg)

#visualisation genome
def visualisation(self, nb_line, genetic_code):

    try:
        display_len = int(nb_line)*SCREEN_SIZE
        conversion = int(len(self.record)) / display_len
    except ValueError:
        display_len = int(len(self.record)/3) if nb_line == 'aa' else int(len(self.record))
        conversion = 3 if nb_line == 'aa' else 1



    compatible_dico = collections.OrderedDict()
    compatible_dico['match'] = buildCompatibleGroup(self, self.matchs)[::-1]
    compatible_dico['cds'] = buildCompatibleGroup(self, self.cds)
    compatible_dico['pep'] = buildCompatibleGroup(self, self.peptides)
    compatible_dico['unannotated_region'] = buildCompatibleGroup(self, self.unannotated_region)

    strings = []
    seq_type_dico = {
        'match': ('[', '#', ']' , ' '),
        'cds': ('=', '=', '>', '-') ,
        'pep':('|', '+', '|' , ' '),
        'unannotated_region':('|', '~', '|' , ' ')
        }
    color_end = '\033[0m'
    overlap_col = '\033[91m'

    # color = '$'
    # color_end= '$'
    for seq_type, compatible_group in compatible_dico.items():
        left, central, right, neutral = seq_type_dico[seq_type] # different symbol if it protein or peptide
        for group in compatible_group:
            # print('GROUPPP')
            string = neutral*display_len
            # print(string)

            for seq in group:
                i_start, i_end = getStringIndices(seq, conversion)

                string = string[:i_start] + left+central*(i_end - i_start-1)+ right + string[i_end+1:]

                if seq.__class__.__name__ ==  'Protein':

                    # print(seq.getSequenceRecord(self.organism, self.taxon_id, self.record, genetic_code).seq)
                    if  nb_line == 'aa':
                        seqaa = seq.getSequenceAA(self.record, genetic_code)

                        string = string[:i_start] + seqaa + string[i_start+len(seqaa):]
                        continue
                    if  nb_line == 'nt':
                        seqaa = seq.getSequenceAA(self.record, genetic_code)
                        seqaa = ''.join([a*3 for a in seqaa])
                        string = string[:i_start] + seqaa + string[i_start+len(seqaa):]
                        continue

                    #display positional number
                    # if seq.polyprotein_number:
                    #     nb = str(seq.polyprotein_number)
                    #     string = string[:i_start+2] +nb+ string[i_start+2+len(nb):]
                    #
                    #Display protein number
                    nb = str(seq.number)
                    string = string[:i_start+2] +nb+ string[i_start+2+len(nb):]

                    if False and seq.ribosomal_slippage:
                        for position, shift in seq.ribosomal_slippage.items():
                            string_postion = round(position/conversion)
                            string_postion = string_postion if string_postion<i_end else string_postion-1
                            shift_string = '+' if shift > 0 else "-"
                            shift_string += str(abs(shift))

                            string = string[:string_postion-len(shift_string)+1] + shift_string + string[string_postion+1:]

                else:
                    name = seq.getNameToDisplay()
                    string = string[:i_start+2] +name+ string[i_start+2+len(name):]
            # if seq_type == "match" and nb_line == 1:
                # print(string)

                # string = re.sub(r"(\[#<[A-Za-z0-9]+#+\]?)", r"{}\1{}".format(color, color_end), string)
            # print(string)
            strings.append(string)

    ##REGEX
    # overlap_match re.compile(r"\[#<"g)
    # match = re.compile(r"\[#[^<]"g)
    # sub_string = re.sub(r"(\[#<[A-Za-z0-9]+#+\]?)", r"{}\1{}".format(color, color_end), sub_string)

    unfinished_col_list = ["" for s in strings]
    for i in range(0, display_len, SCREEN_SIZE):
        print('/'*SCREEN_SIZE)
        print()
        for i_line, string in enumerate(strings):
            sub_string = string[i:i+SCREEN_SIZE]
            if unfinished_col_list[i_line]:
                sub_string = re.sub(r"^([^ \[]+)", r"{}\1{}".format(unfinished_col_list[i_line], color_end), sub_string)
            sub_string = re.sub(r"(\[?#?<[A-Za-z0-9]+#*\]?)", r"{}\1{}".format(overlap_col, color_end), sub_string)

            ansi_list = re.findall(r"(\x1B\[[0-?]*[ -/]*[@-~])" ,sub_string[:-1]) # -1 not not catch the end color character of the end of the line

            unfinished_col_list[i_line] = '' if not ansi_list or ansi_list[-1] == color_end else ansi_list[-1]
            # if ansi_list and ansi_list[-1] != color_end:
            #     unfinished_col_dict[i_line] = ansi_list[-1]
            print(sub_string)

    # protein = '\n'.join(strings)

    print(color_end)
    # print(protein)
    # return 'protein'

def getStringIndices(self, conversion):
    # print(conversion,' start', self.start/conversion, 'end', (self.end-1)/conversion)
    # print(conversion,' start', int(self.start/conversion), 'end', int((self.end-1)/conversion))
    return (int(self.start/conversion), int((self.end-1)/conversion))
    # return round(self.bp_obj.location.start/conversion), round(self.bp_obj.location.end/conversion)



def buildCompatibleGroup(self, set_of_sequence):
        sequences = set_of_sequence.copy() #sorted(list(set_of_sequence), key=lambda x: x.bp_obj.location.start, reverse=False)
        #Build the set of non overlaping protein

        for i, seq in enumerate(sequences):

            for seq_next in list(sequences)[i+1:]:

                if not seq.overlap(seq_next):
                    seq.non_overlapping_prot.add(seq_next)
                    seq_next.non_overlapping_prot.add(seq)

        #Build list of sequences that don't overlap
        compatible_groupes = []
        used_seq = []
        while len(sequences) > 0:
            seq = sequences.pop()
            group = [seq]
            # finding the sequence that are compatible with seq and within each other
            group_compatible = seq.non_overlapping_prot # set of sequence that are compatible with the sequence stored in the list group
            while group_compatible:
                try:
                    compatible_seq = (group_compatible & sequences).pop() #extarct one seq that is compatibe with group
                except KeyError:
                    break
                group_compatible = group_compatible & compatible_seq.non_overlapping_prot & sequences # build new set 'group_compatible' that store compatible set

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


    sp_treshold=90
    taxonomy_file="data/taxonomy/taxonomy_virus.txt"
    expected_file =  "data/taxonomy/polyprotein_expectation_by_taxon.csv"

    gff_file = 'data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'

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
    print('VISUALISATION OF THE RefSEq GENOME FROM THE TAXON {} THAT HAVE AT LEAST {} ANNOTATED PEPTIDE'.format(taxon,minimum_nb_peptide ))
    print('*-'*100)
    for i, gb_dico in enumerate(gbff_iter):
        # print(gb_dico)

        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        taxon_id = gb_dico['taxon_id']
        print(gb_file)

        if i%200== 0:
            # continue
            print(i)
        # print('genetic code', genetic_code)
        # print(gb_file)
        visualisation_main(gb_file, genetic_code, gff_file, nb_line, minimum_nb_peptide, taxon_id, sp_treshold)

    print(i+1, 'Genome analysed from taxon', taxon)
