#!/usr/bin/env python3

import gzip
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation

SCREEN_SIZE = 200

###########
# Functions
##########


def gb_file_parser(gb_file, taxon_id, sp_treshold):
    # Main function that manage the parsing of the gb_file

    genome = Genome(gb_file, taxon_id)

    # Manage compress and not compress gb file
    proper_open = gzip.open if gb_file.endswith('.gz') else open

    with proper_open(gb_file, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "genbank")):
            segment = Segment(record, gb_file)
            genome.segments.append(segment)

            segment.getMatpeptidesAndPolyproteins()

            if not segment.peptides:  # if no peptide annotation we don't need to do the next step of the loop
                continue

            segment.checkPeptideRedundancy()  # remove the redundant peptide
            segment.checkSubPeptides()
            segment.associatePepWithProt()

            # segment.checkForSlippage()
            segment.identifySubProtein()

            segment.getCleavageSites()
            segment.identifyAnnotatedPolyproteins(sp_treshold)

    return genome


def getGenomicLocation(sequence, start_in_seq, end_in_seq):
    """
        We know the start and end position of a mature peptide/interpro domain in the protein sequence
        and we want to know the start and end position of this peptide in the genome.
        Take care of potential frameshift in the sequence
    """

    assert len(sequence.bp_obj.location) / \
        3 >= end_in_seq, f'{len(sequence.bp_obj.location) / 3} >= {end_in_seq}'
    # print(sequence)
    # print('GET GENOMIC FROM', start_in_seq, end_in_seq)
    none_location = FeatureLocation(0, 0, ref='tmp')
    previous_seq_location = none_location
    pep_location = none_location  # to be able to add FeatureLocation from it
    for seq_location in sequence.bp_obj.location.parts:

        if start_in_seq - 1 < len(seq_location + previous_seq_location)/3:
            # petide starts before are in the current seq_location

            if len(pep_location):
                # pep_location contains already a valid Featurelocation found earlier
                # that means start peptide is before the current seq_location
                # Then the new Featurelocation starts where the seq_location starts
                start_pep_nt = seq_location.start
            else:
                # the peptide start inside the current location
                start_pep_nt = seq_location.start + \
                    (start_in_seq-1)*3 - len(previous_seq_location)

            if end_in_seq <= len(seq_location + previous_seq_location)/3:
                # The end of the peptide is inside the current location
                # We can add a Feature Location and break the loop
                # The peptide doesn't continue further
                end_in_seq = start_pep_nt + \
                    (end_in_seq - start_in_seq + 1)*3 - len(pep_location)
                pep_location += FeatureLocation(start_pep_nt,
                                                end_in_seq, strand=seq_location.strand)
                break
            else:
                # The peptide doesn't end in the current location
                # We add a intermediate Feature Location to pep_location
                pep_location += FeatureLocation(start_pep_nt,
                                                seq_location.end, strand=seq_location.strand)

        previous_seq_location += seq_location

    pep_location.parts.remove(none_location)  # clean pep location

    if len(pep_location.parts) == 1:
        pep_location = pep_location.parts[0]
    return pep_location


def getProteinLocation(peptide, cds):
    # print('cds', cds.bp_obj.location)
    # print("pep", peptide.bp_obj.location, peptide.bp_obj.location.parts[0].start,
    #       peptide.bp_obj.location.parts[-1].end, len(peptide.bp_obj.location))
    assert len(cds) % 3 == 0 and len(peptide) % 3 == 0
    start_pep_nt = peptide.bp_obj.location.parts[0].start  # peptide.start
    len_previous_parts = 0  # length of the cds from its start to the current location
    # print('start_pep_nt', start_pep_nt)
    for cds_location_part in cds.bp_obj.location.parts:
        # print("cds_location_part", cds_location_part)
        if start_pep_nt in cds_location_part:
            # print('  START IS IN LOCATION')
            len_cds_until_pep_start = len_previous_parts + \
                (start_pep_nt - cds_location_part.start)
            # print("   len_cds_until_pep_start", len_cds_until_pep_start, len_cds_until_pep_start % 3)
            if len_cds_until_pep_start % 3 == 0:  # if the cds is in frame with the peptide start
                start_aa = len_cds_until_pep_start/3 + 1
                end_aa = start_aa + len(peptide)/3 - 1
                return int(start_aa), int(end_aa)
                break
            # else:
            #     logging.warning('Problem with peptide protein location')
            #
        len_previous_parts += len(cds_location_part)


############
# Class
############


class Genome:
    """Class for storing info about genome from a gb file"""

    def __init__(self, gb_file, taxon_id):
        """Init method of genome obj"""
        self.gb_file = gb_file
        self.segments = []
        self.matchs = []
        self.taxonomy = None
        self.organism = None
        self.taxon_id = taxon_id
        # self.expectation_node = None
        # self.peptide_expectation = None
        # self.polyprotein_expectation = None
        # self.variable_polyprotein_expectation = None
        # self.expectation_info = {}
        Protein.COUNTER = 0
        Peptide.COUNTER = 0

    def __len__(self):
        return sum([len(s) for s in self.segments])

    def hasEnoughPeptide(self):
        total_final_peptides = sum(
            [len(segment.peptides)-len(segment.parent_peptides) for segment in self.segments])
        return True if self.peptide_expectation <= total_final_peptides else False


class Segment:
    """Segment consists of a distinct record from the biopython parsing of a gb file"""

    POSSIBLE_TYPE = set()

    def __init__(self, record, gb_file):
        self.record = record
        self.gb_file = gb_file
        self.source = None
        self.polyproteins = set()  # polyproteins list
        self.parental_polyproteins = []
        self.peptides = set()  # mat_peptides list
        self.matchs = set()
        self.parent_peptides = set()
        self.sub_peptides = set()
        self.cds = set()
        self.unannotated_region = set()
        self.taxon_id = None
        self.organism = record.annotations['organism']
        self.cleavage_sites = []  # list of cleavage_site object
        self.relevant_annotation = False

        self.polyprotein_expectation = None

    def __len__(self):
        return len(self.source)

    # def getSubProteins(self):
    #     return [p for p in self.cds if p.parental_prot ] # return protein that have parental prot

    # def getPolyWithFeature(self, feature):
    #     return [p for p in self.polyproteins if getattr(p, feature)]

    def identifyAnnotatedPolyproteins(self, sp_treshold):
        for cds in self.cds:
            strand = cds.bp_obj.strand
            if len(cds.peptides) == 0:
                continue
            else:
                cds.polyprotein = True
            # (self.start, self.end) if self.bp_obj.strand == 1 else (self.end, self.start)
            cds_start, cds_end = cds.realStart(), cds.realEnd()

            if 0 < len(cds.peptides) <= 2:
                # SIGNAL P Identification
                for pep in cds.peptides:
                    pep_start, pep_end = (pep.start, pep.end) if cds.bp_obj.strand == 1 else (
                        pep.end, pep.start)
                    if cds_start == pep_start and len(pep)/3 < sp_treshold:
                        # if pep.start_aa(cds) != 1:
                        #     print('LENGTH IS NOT 1')
                        #     print(pep.start_aa(cds))
                        cds.polyprotein = False
                        cds.non_polyprotein_explanation = "Signal Peptide"

                # Peptide cooveering all the length identification
                if len(cds.peptides) == 1 and (pep.end_aa(cds) == len(cds)/3-1 or pep.end_aa(cds) == len(cds)/3):
                    if pep.end_aa(cds) == len(cds)/3:
                        logging.info('Peptide annotation include stop codon in {}|{}'.format(
                            self.taxon_id, cds.protein_id))

                    # cds_start_sp_threshold = cds_start + (sp_treshold*3)*strand
                    # pep start has to be between the cds start and the cds + signal p strshold
                    # if  min(cds_start_sp_threshold,cds_start ) <= pep_start <= max(cds_start_sp_threshold,cds_start ):
                    if pep.start_aa(cds) < sp_treshold:
                        # potential protein with signal peptide where only the mature peptide is annotated
                        # it happens that the whole sequence is covered by a signle mat_peptide with the exception of the first and last condon
                        # example: 11886 Rous sarcoma virus	Viruses;Retro-transcribing viruses;Retroviridae;Orthoretrovirinae;Alpharetrovirus
                        cds.polyprotein = False
                        logging.info('single peptide annotaion covering the whole CDS in {}|{}'.format(
                            self.taxon_id, cds.protein_id))
                        if 1 <= pep.start_aa(cds) <= 2:  # if it starts at 1 or 2 in protein sequence

                            cds.non_polyprotein_explanation = "single peptide annotation covering the whole CDS"
                        else:
                            cds.non_polyprotein_explanation = "single peptide annotation covering almost the whole CDS"

            # intein/extein Identification
            if len(cds.peptides) == 2:
                # sort by start and select the first peptide no matter if the strand is -1
                # because in case of intein the peptide cover the begining and the end of the CDS
                pep, pep_middle = sorted(list(cds.peptides), key=lambda x: x.start, reverse=False)

                cds_start_area = [cds.realStart(), cds.realStart()+12*cds.bp_obj.strand]
                cds_start_area.sort()
                # to identify a protein with an intein we use the mat peptide annotation that flank the mat peptide of the intein
                # This peptide annotation covers the begining of the CDS and the end, it has then 2 parts and its length is smaller than the CDS.
                # Need to be a bit flexible here because in  the genome of 1163482 the extein strat 9 nt after the actual start..
                if cds_start_area[0] <= pep.realStart() <= cds_start_area[1] and cds.realEnd() - 3*cds.bp_obj.strand == pep.realEnd():
                    if len(pep) < len(cds)-3 and len(pep.location.parts) > 1:
                        cds.polyprotein = False
                        cds.non_polyprotein_explanation = "Intein outline: extein surounds intein"
                    # in some genome the mature peptide of the intein cover all the CDS and the intein peptide a small part
                    if len(cds)-3 >= len(pep) >= len(cds)-12 and pep_middle.start in pep.location and pep_middle.end in pep.location:
                        cds.polyprotein = False
                        cds.non_polyprotein_explanation = "Intein outline: extein includes intein"

            if len(cds.cleavage_sites) == 0 and len(cds.peptides) > 0 and cds.polyprotein:
                print(self.taxon_id)
                print(cds)
                [print(pep) for pep in cds.peptides]
                print(cds.polyprotein)
                cds.polyprotein = False
                logging.warning('No cleavage site identify  and no explanation in {} from {}'.format(
                    cds.protein_id, self.taxon_id))
                cds.non_polyprotein_explanation += "No cleavage site identify and no explanation"
                # print('NO CEAVAGE SITE AND NO EXPLANATION')
                # input()

    def getMatpeptidesAndPolyproteins(self):
        for feat in self.record.features:
            Segment.POSSIBLE_TYPE.add(feat.type)

            if feat.type == 'source':
                for item in feat.qualifiers['db_xref']:
                    if item.startswith('taxon'):
                        self.taxon_id = item.replace('taxon:', '')

                self.source = feat

            elif feat.type == "mat_peptide" or feat.type == "sig_peptide" or feat.type == "proprotein":
                self.peptides.add(Peptide(feat))
            elif feat.type == 'CDS':
                prot_obj = Protein(feat, self)
                self.cds.add(prot_obj)

    def associatePepWithProt(self):
        # Should be changed at the end to catch polyprot that are not annotated
        # Here we want to find the polyprot that have a mat peptide so we are pretty sure that thz protein is a polyprotein
        # EDIT: The protein that has some peptide is no longer consider as polyprotein
        # To be consider as polyprotein the cds undergoes the identifyExpectedPolyprotein() function
        # No Anymore the case

        for cds in self.cds:
            # print('CDS', cds)
            for pep in self.peptides:
                # print('  PEP', pep)
                if pep in cds:
                    cds.peptides.add(pep)
                    pep.polyproteins.add(cds)

            if cds.peptides:
                cds.proteinCoverage()
                self.unannotated_region.update(cds.unannotated_region)

        # SMALL check up to be sure that every peptides have been assigned to at least one protein
        for pep in self.peptides:
            if not pep.polyproteins:
                logging.warning("The peptide {} {} was not assigned to any cds in {}".format(
                    pep.number, str(pep.bp_obj.location), self.taxon_id))

    def checkSubPeptides(self):
        # Check if peptide is included in a bigger peptide
        for i, pep in enumerate(self.peptides):
            for pep_next in self.peptides:
                if pep != pep_next and pep in pep_next:
                    self.parent_peptides.add(pep_next)
                    pep_next.parent_peptide = True
                    self.sub_peptides.add(pep)

    def checkPeptideRedundancy(self):
        # Remove Peptide that have similat start and end...
        peptides = list(self.peptides)

        for i, pep in enumerate(peptides):

            for pep_next in peptides[i+1:]:

                # if pep.location.start ==  pep_next.location.start and pep.location.end ==  pep_next.location.end:
                if pep.location == pep_next.location:
                    self.peptides.remove(pep_next)
                    pep.redundant_pep.append(pep_next)

    def identifySubProtein(self):
        list_prot = list(self.cds)

        for i, prot in enumerate(list(list_prot)):
            for prot_next in list(list_prot)[i+1:]:

                if prot.isIncludedIn(prot_next):
                    prot.parental_prot.append(prot_next)
                    prot_next.sub_prot.append(prot)
                    prot.checkforAlternativeStart(prot_next)

                elif prot_next.isIncludedIn(prot):
                    prot_next.parental_prot.append(prot)
                    prot.sub_prot.append(prot_next)
                    prot_next.checkforAlternativeStart(prot)

    def getCleavageSites(self):

        # detected cleavage site located at
        # the border of the cds are ignore
        # need to give a margin: default is 3 aa at the begining and end of cds
        margin = 0  # 3  # in aa

        cleavage_site_dict = {}

        peptide_list = sorted(list(self.peptides | self.unannotated_region),
                              key=lambda x: x.start, reverse=False)

        for pep in peptide_list:

            strand = pep.bp_obj.strand

            for cds in pep.polyproteins:
                margin_end = (len(cds)-3)/3 - margin
                type = "start"
                end_aa_cs = pep.start_aa(cds)
                start_aa_cs = end_aa_cs - 1
                # print("PEP", pep.number, pep.position_prot_relative)
                # print(type, start_aa_cs, end_aa_cs)
                # input()
                if margin < start_aa_cs < margin_end:
                    cs_location = getGenomicLocation(cds, start_aa_cs, end_aa_cs)
                    tuple_location = (cs_location.start, cs_location.end,
                                      cs_location.strand)
                    if tuple_location in cleavage_site_dict:
                        site = cleavage_site_dict[tuple_location]
                        site.update(type, pep, {cds})
                    else:
                        cleavage_site = CleavageSite(cs_location, {cds}, pep, self.taxon_id, type)
                        cleavage_site_dict[tuple_location] = cleavage_site

                type = 'end'  # end of the peptide give the cs
                start_aa_cs = pep.end_aa(cds)
                end_aa_cs = start_aa_cs + 1
                # print("PEP", pep.number, pep.position_prot_relative)
                # print(type, start_aa_cs, end_aa_cs)
                # input()
                if margin < end_aa_cs < (len(cds)-3)/3 - margin:  # -3 to remove the stop codon
                    cs_location = getGenomicLocation(cds, start_aa_cs, end_aa_cs)
                    tuple_location = (cs_location.start, cs_location.end,
                                      cs_location.strand)
                    if tuple_location in cleavage_site_dict:
                        site = cleavage_site_dict[tuple_location]
                        site.update(type, pep, {cds})
                    else:
                        cleavage_site = CleavageSite(cs_location, {cds}, pep, self.taxon_id, type)
                        cleavage_site_dict[tuple_location] = cleavage_site

        self.cleavages_sites_new = list(cleavage_site_dict.values())

    def getCleavageSites_old(self):

        start_cleavage_sites = {}
        # print("-*********** GET CS ***************")
        # (self.peptides | self.unannotated_region):  # sorted(list(self.peptides | self.unannotated_region), key=lambda x: x.start, reverse=False):
        for pep in sorted(list(self.peptides | self.unannotated_region), key=lambda x: x.start, reverse=False):
            strand = pep.bp_obj.strand
            # if pep.number != 5:
            #     continue
            # print('PEPTIDE',pep.number)
            for type, border in [("start", pep.start), ('end', pep.end)]:

                start = border - 3
                end = border + 3
                # if pep.number == 4:
                # print(" ",type, border)
                # print("  CS position", start, end)
                # polyprot that are compatible with the cleavage site
                polyproteins = {
                    poly for poly in pep.polyproteins if poly.start+9 < start < poly.end-9}
                # print('  border compatible with', len(polyproteins))
                if polyproteins:
                    partial_location = False
                    previous_part_end = False
                    previous_part = False
                    overlap = False
                    cleavage_site = False
                    for poly in polyproteins:
                        # if pep.number == 4:
                        #     print(poly)
                        # print('  CHECK in poly', poly.number, poly.bp_obj.location)
                        for i, part in enumerate(poly.bp_obj.location.parts):

                            if previous_part and end in part and start not in part and start not in previous_part:  # cleavage site overlap an intron
                                logging.info('Cleavage site from {} is overlaping an intron in genome {}'.format(
                                    poly.protein_id, self.taxon_id))
                                overlap = True
                                previous_part_end = previous_part.end

                            if overlap and end in part:
                                if partial_location is not False:
                                    # print('WE are in case where border was in the previous part')
                                    # input()
                                    len_cs_previous_part = len(partial_location)
                                    len_cs_actual_part = 6 - len_cs_previous_part
                                    new_end = part.start + len_cs_actual_part
                                    location = partial_location + \
                                        FeatureLocation(part.start, new_end, strand=strand)
                                    tuple_position = (start, new_end)

                                elif previous_part_end is not False:
                                    # print('WE are in case where border is in the current part')
                                    # input()
                                    partial_location = FeatureLocation(
                                        part.start, end, strand=strand)
                                    # print('partial_location',partial_location)
                                    len_previous_part = 6 - len(partial_location)
                                    new_start = previous_part_end - len_previous_part
                                    # print("len_previous_part", len_previous_part)
                                    # print("new_start",new_start )
                                    location = FeatureLocation(
                                        new_start, previous_part_end, strand=strand) + partial_location
                                    tuple_position = (new_start, end)

                                if tuple_position in start_cleavage_sites:
                                    site = start_cleavage_sites[tuple_position]
                                    # this cleavage site has been already treated
                                    site.update(type, pep, {poly})
                                else:  # we create the CS object
                                    cleavage_site = CleavageSite(
                                        location, polyproteins, pep, self.taxon_id, type)
                                    # we store the start and end of the CS in the dico to update it and no recreate when it is nececessary
                                    start_cleavage_sites[tuple_position] = cleavage_site
                                    # print(cleavage_site)
                                    # print(location, len(location))
                                    # input()
                                    assert(len(location) == 6), 'Cleavage site has a length of {} location {}'.format(
                                        len(location), location)
                                partial_location = False
                                previous_part_end = False
                                overlap = False
                                break  # we break and check a new polyprotein

                            if start in part and end in part:  # cleavage site is included in the prot part
                                if (start, end) in start_cleavage_sites:
                                    site = start_cleavage_sites[(start, end)]
                                    # this cleavage site has been already treated
                                    site.update(type, pep, {poly})
                                else:  # creation of the cleavage site object
                                    location = FeatureLocation(start, end, strand=strand)
                                    cleavage_site = CleavageSite(
                                        location, polyproteins, pep, self.taxon_id, type)
                                    start_cleavage_sites[(start, end)] = cleavage_site
                                break  # we break and check a new polyprotein

                            if start in part and end not in part:
                                overlap = True
                                if border-1 in part:  # the cleavage site overlap two part of the protein
                                    # print("border", border, 'in part', part)
                                    partial_location = FeatureLocation(
                                        start, part.end, strand=strand)
                                    # print('modulo start prot part',part.start%3)
                                    # print('modulo start CS', start%3)
                                    # print("partial location", partial_location, len(partial_location))
                                    # # we check the next part of the protein -->
                                    continue

                                else:
                                    # print("border", border, 'not in part', part)
                                    previous_part_end = part.end
                            previous_part = part

                    # print("Cleavage site:", type, border,' start end', start, end, '\t',  cleavage_site.start_aa(poly), cleavage_site.end_aa(poly))
                    if cleavage_site:
                        self.cleavage_sites.append(cleavage_site)

    def isSegmentAnnotationRelevant(self):
        return False if any([p for p in self.polyproteins if not p.isAnnotationRelevant()]) else True


class Sequence:

    def __contains__(self, pep):
        # test if a peptide is in self. Self can be a CDS or another peptide
        # the qualifiers gene is no more used to know if a peptide belong to a polyprotein
        # Only the position of the peptide tel us if the peptide belongs to a poly
        # pep.location.parts[0].start need to use part to get the real start position and the min and max position in the sequence
        # otherwise when the sequence is circular like in Woodchuck hepatitis virus|35269

        pep_start = pep.location.parts[0].start
        pep_end = pep.location.parts[-1].end
        if not (pep.location.parts[0].start in self.bp_obj.location and pep.location.parts[-1].end-1 in self.bp_obj.location and pep.bp_obj.strand == self.bp_obj.strand):

            return False
        # print("\n","***"*3,str(self)[:8],"***"*3)
        # print("PROT location", self.bp_obj.location)
        # print("PEP location", pep.bp_obj.location)
        pep_parts = iter(pep.location.parts)
        pep_part = next(pep_parts)
        case2 = False
        shift = 0
        protein_start = self.bp_obj.location.parts[0].start
        previous_end = protein_start
        prot_len = 0
        for sub_location in self.bp_obj.location.parts:
            # shift += sub_location.start - previous_end
            #
            # print('\n==Protein Sub location==', sub_location)
            # #
            # print('Shift', shift)
            # print('pep part', pep_part)
            # print('prot part', sub_location)
            if case2:
                # print('\n-WE are in case 2', )

                if pep_part.start == sub_location.start:  # check if the 2 parts start are the same. So if the join() is similar
                    # print('--the start are the same')
                    # print("--pep parts", pep_parts)
                    try:
                        pep_part = next(pep_parts)
                    except StopIteration:
                        if pep_end in sub_location:  # No more part for the peptide and it end is included in the subprot prat
                            # print('---end is included in the prot part')
                            # print('TRUUE')
                            return True
                        else:
                            return False
                    else:
                        # print('--We go to next part')

                        continue
                else:  # The peptide in somehow not folowing the prot part
                    # print('FALSE')
                    # print('-The peptide in somehow not folowing the prot part')
                    return False
            # Case 1 : The peptide is included in one part of the protein
            if pep_start in sub_location and pep_end in sub_location and len(pep.location.parts) == 1 and pep_start % 3 == (sub_location.start - shift) % 3:
                # print("-the peptide is included in the subprot part")
                return True

            # print(" pep_part.start%3 ",  pep_part.start%3 )
            # print("(sub_location.start - shift)%3", sub_location.start, '-' ,shift,')%3', (sub_location.start - shift)%3)
            # Case 2 Peptide is overlapping this part with the next one
            if pep_start in sub_location and pep_end not in sub_location and pep_part.end == sub_location.end and pep_part.start % 3 == (sub_location.start - shift) % 3:
                try:
                    pep_part = next(pep_parts)
                except StopIteration:
                    return True  # the peptide end where the ribosomal shifting occurs
                case2 = True
                # print('-the peptide is overlapping the prot part')

            prot_len += len(sub_location)
            shift = prot_len % 3

        # SPECIAL CASE: intein found in some CDS of dsDNA (example Phicbkvirus)
        if pep.realStart() == self.realStart() and pep.realEnd() == self.realEnd() - 3*self.bp_obj.strand and len(pep) < len(self)-3 and len(pep.location.parts) > 1:
            # print("SPECIAL CASE: intein ")
            return True

        # print("False ..")
        return False

    def __str__(self):

        string = 'Sequence: from {} to {}'.format(self.start, self.end)
        # string += '\n'.join([str(p) for p in self.peptides])
        # print('ddd')
        return string

    def __len__(self):
        return len(self.bp_obj)
        # return self.end - self.start +1

    def overlap(self, seq):
        if self.start <= seq.end <= self.end or seq.start <= self.end <= seq.end:
            return True
        return False

    def realStart(self):
        # return position of the real start of the sequence independantly of the strand.
        return self.start if self.bp_obj.strand == 1 else self.end

    def realEnd(self):
        # return position the real end of the sequence independantly of the strand.
        return self.end if self.bp_obj.strand == 1 else self.start

    def getGenomicPositions(self, prot_start):
        # TODO CHANGE INTO getGenomicLocation
        # for uncovered region and match to get the genomic position and not the prot relative position
        self.start = prot_start - 1 + self.start_in_prot*3 - 2 - 1
        self.end = prot_start + self.end_in_prot*3 - 1

    def getNameToDisplay(self):
        return 'None'


class Protein(Sequence):

    COUNTER = 0

    def __init__(self, biopyth_obj, segment):
        # self.start = start
        # self.end = end

        self.bp_obj = biopyth_obj
        self.segment = segment
        self.start, self.end = (biopyth_obj.location.start, biopyth_obj.location.end)

        self.gene_type = biopyth_obj.type
        self.protein_id = "Unknown" if 'protein_id' not in biopyth_obj.qualifiers else biopyth_obj.qualifiers[
            'protein_id'][0]
        self.product = "Unknown" if 'product' not in biopyth_obj.qualifiers else biopyth_obj.qualifiers[
            'product'][0]
        self.peptides = set()
        self.predicted_mat_peptides = []  # list of obj Predicted_peptide
        self.matchs = []
        self.unannotated_len = len(biopyth_obj)/3
        self.polyprotein = False
        self.non_polyprotein_explanation = ''
        self.unannotated_region = []
        self.non_overlapping_prot = set()  # for the visualisation
        # self.ribosomal_slippage = {}
        # self.readthrough = {}
        self.parental_prot = []
        self.sub_prot = []
        # self.alternative_start = False
        self.cleavage_sites = []
        self.predicted_cleavage_sites = []
        self.black_listed_cleavage_sites = []
        self.border_sites = []  # site that are at the close border of CDS
        self.polyprotein_number = 0.0
        self.sequence = ''

        self.annotation_quality = 0

        Protein.COUNTER += 1
        self.number = Protein.COUNTER

    def __len__(self):
        return len(self.bp_obj)

    def __str__(self):

        # status += str(self.polyprotein_number)
        string = 'protein {}: {}nt | {}aa | location:{} \n'.format(
            self.number, len(self), int(len(self)/3), self.bp_obj.location)
        string += 'mature paptides {} | unannotated region {}\n'.format(
            len(self.peptides), len(self.unannotated_region))

        string += f'polyprotein outline: {self.polyprotein}  {self.non_polyprotein_explanation}\n'
        string += f'product: {self.product}\n'
        if self.parental_prot:
            string += f'is included in other cds: {"|".join([p.protein_id for p in self.parental_prot])}\n'
        if self.sub_prot:
            string += f'includs other cds: {"|".join([s.protein_id for s in self.sub_prot])}\n'

        return string

    def get_unnannotated_version(self):
        unnannotated_version = Protein(self.bp_obj, self.segment)
        unnannotated_version.annotation_version = self
        return unnannotated_version

    def getSequenceAA(self, record, genetic_code):

        if self.sequence:
            return self.sequence

        elif 'translation' in self.bp_obj.qualifiers:
            seq = self.bp_obj.qualifiers["translation"][0]
        else:
            seq = self.bp_obj.location.extract(record).seq.translate(table=genetic_code)
            seq = str(seq)
        seq = seq.upper()
        if not self.polyprotein:
            return seq
        for site in self.cleavage_sites:
            site_position = site.start_aa(self)-1  # -1 to be in base 0

            site_extraction = site.bp_obj.location.extract(record).seq.translate(table=genetic_code)
            site_seq = str(self.bp_obj.location.extract(record).seq.translate(
                table=genetic_code, to_stop=False))[site_position:site_position+2]

            if site_extraction.upper() != site_seq.upper():
                # print(site_position)

                logging.warning('{}:cleavage site ({}) is not similar in extraction from genome sequence and the protein sequence: extraction:{} and prot_seq:{}'.format(
                    self.protein_id, site.bp_obj.location, site_extraction, site_seq))
            if site_position < 0:
                logging.warning('{}:The cleavage site position in the sequence is wrong {}'.format(
                    self.protein_id, site_position))

            # Lower case cleavage sites
            if any([True for peptide in site.peptide_positions['left'] if type(peptide).__name__ != "UnannotatedRegion"]):
                seq = seq[:site_position] + seq[site_position:site_position+1].lower() + \
                    seq[site_position+1:]
            if any([True for peptide in site.peptide_positions['right'] if type(peptide).__name__ != "UnannotatedRegion"]):
                seq = seq[:site_position+1] + seq[site_position +
                                                  1:site_position+2].lower() + seq[site_position+2:]

        self.sequence = seq

        return seq

    def isAnnotationRelevant(self):
        # Simple rule to be sure that the annotation is relevant
        #  the number of uncvered region < number of peptide
        # To not take into acount this kind of annotation
        # ---==1.0============================================================================================================================================================================>------------------------
        #           |+5|         |+6+++|     |+7++++|
        #    |~8~~~~|  |~9~~~~~~~|     |~10~~|      |~11~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

        if len(self.unannotated_region) >= len([p for p in self.peptides if not p.parent_peptide]):
            return False
        else:
            return True

    def getSequenceRecord(self, organism, header, genetic_code, subPosition=False):

        record = self.segment.record
        seq = self.getSequenceAA(record, genetic_code)

        # print(seq)
        if subPosition:
            seq = seq[subPosition[0]:subPosition[1]]
            header += '|{}:{}'.format(subPosition[0], subPosition[1])

        return SeqRecord(Seq(seq, generic_protein), id=header, description="{}|{}".format(self.product, organism))
        # print("PROTTTT")
        # print(seq_to_write)
        # SeqIO.write(seq_to_write, file_handle,"fasta")

    def proteinCoverage(self):
        # Check if the mat peptide associated with the protein are covering all the prot or not
        # If not create unannotated_region to fill the gaps
        strand = self.bp_obj.strand
        self.unannotated_len = 0
        iter_parts = iter(self.bp_obj.location.parts)
        protpart = next(iter_parts)
        current_po = protpart.start
        partial_location = None
        # for pep in sorted(list(self.peptides), key=lambda x: x.start, reverse=False) :
        peptides = sorted(list(self.peptides), key=lambda x: x.start, reverse=False)
        pep_i = 0
        while pep_i < len(peptides):
            pep = peptides[pep_i]
            if current_po < pep.start:  # if pep start is after the current positon then a unnatotade region need to be created

                if protpart.end < pep.start:  # that means the next peptide start is located in the next part of the protein
                    location = FeatureLocation(current_po, protpart.end, strand=strand)
                    partial_location = partial_location + location if partial_location else location

                    protpart = next(iter_parts)
                    current_po = protpart.start
                    continue

                else:
                    # an unannotated position is created.
                    location = FeatureLocation(current_po, pep.start, strand=strand)
                    if partial_location:
                        location = partial_location + location
                        partial_location = None
                    # in case of strand - the start is the end and then we don't want to take into account the stop codon
                    if not(strand == -1 and len(location) <= 3):

                        unannotated_seq = UnannotatedRegion(location, self, strand)
                        self.unannotated_region.append(unannotated_seq)
                    current_po = pep.end

            elif current_po < pep.end:
                if protpart.end < pep.end:  # the peptide is overlapping two part of the prot
                    protpart = next(iter_parts)  # we take the next part

                current_po = pep.end
            pep_i += 1

            # else: #Another case where curent_po > pep.end but then we don't change the current po and we move to the next pep
            #     current_po =  current_po

        if current_po < protpart.end:  # we didn't not arrive at the end of the protein
            location = FeatureLocation(current_po, protpart.end, strand=strand)
            final_location = partial_location + location if partial_location else location

            for protpart in iter_parts:
                # print(protpart)
                final_location += FeatureLocation(protpart.start, protpart.end, strand=strand)

            if len(final_location) > 3:
                unannotated_seq = UnannotatedRegion(final_location, self, strand)

                self.unannotated_region.append(unannotated_seq)
            # unannotated_seq = UnannotatedRegion(current_po, int(len(self.bp_obj)/3 -1), self) #SeqFeature(FeatureLocation(current_po, int(len(self.bp_obj)/3 -1)), type="unannotated_region", qualifiers={'note':'Position given in amino acid', "start_aa":current_po-1, 'end_aa':int(len(self.bp_obj)/3 -1-1)})
            # unannotated_seq = Peptide(unannotated_seq_feature)

            # self.unannotated_len += len(self.bp_obj)/3 - current_po
        # poly.qualifiers['unannotated_region'] = unannotated_region
        # poly.qualifiers['peptide_coverage'] = len(poly)/3 - unannotated_len -1 # -1 because we dont count the stop codon
        # print(unannotated_position)

    def isIncludedIn(self, poly_next):

        next_location = poly_next.bp_obj.location
        start, end = self.bp_obj.location.start, self.bp_obj.location.end
        if start in next_location and end in next_location and start % 3 == next_location.start % 3:
            return True
        return False

    def checkforAlternativeStart(self, poly):

        if self.end == poly.end and self.start != poly.start and self.start % 3 == poly.start % 3:
            self.alternative_start = True


class Peptide(Sequence):

    COUNTER = 0

    def __init__(self, bp_obj):
        assert len(bp_obj) % 3 == 0, "length of the peptide nucleotide sequence is not a multiple of 3"
        self.bp_obj = bp_obj
        self.location = bp_obj.location
        self.start = bp_obj.location.start
        self.end = bp_obj.location.end
        self.redundant_pep = []
        self.non_overlapping_prot = set()  # for the visualisation
        self.parent_peptide = False

        self.polyproteins = set()

        self.position_prot_relative = {}

        # Domains attr related: key=domain value len of the domain on the peptide
        self.included_domains = {'fully': [], "partially": []}
        self.overlapped_by_domain = {'left': [], "right": []}

        Peptide.COUNTER += 1
        self.number = Peptide.COUNTER
        # self.genome_start = None
        # self.genome_end = None

    def __str__(self):

        string = 'Peptide {}: from {} to {} | belongs to {}\n'.format(
            self.number, self.start, self.end, [prot.number for prot in self.polyproteins])
        string += '{} nt | {}aa | location :{}  \n'.format(
            len(self), len(self)/3, self.bp_obj.location)
        string += '  Position in protein\n'
        for prot in self.polyproteins:
            string += '    protein {} from {} to {}\n'.format(
                prot.number, self.start_aa(prot), self.end_aa(prot))

        return string

    def start_aa(self, prot):
        if prot.protein_id not in self.position_prot_relative:
            self.get_position_prot_relative(prot)

        return self.position_prot_relative[prot.protein_id][0]

    def end_aa(self, prot):

        if prot.protein_id not in self.position_prot_relative:
            self.get_position_prot_relative(prot)

        return self.position_prot_relative[prot.protein_id][1]

    def get_position_prot_relative(self, prot):
        if prot.protein_id not in self.position_prot_relative:
            start_aa, end_aa = getProteinLocation(self, prot)
            self.position_prot_relative[prot.protein_id] = (start_aa, end_aa)

        # return self.position_prot_relative[prot.protein_id]

    def getProteinPosition_old(self, prot):

        len_previous_part = 0
        bp_prot = prot.bp_obj
        # Searching the peptides position protein relative
        # Due to ribosomal slippage the conversion is not trivial
        for subprotpart in bp_prot.location.parts:
            # print("subprotpart ", subprotpart)
            # print(self.start, self.start%3)
            # print(subprotpart.start, subprotpart.start%3 )
            len_from_prot_start_to_pep_start = len_previous_part + subprotpart.start + self.start
            if subprotpart.start <= self.start <= subprotpart.end and len_from_prot_start_to_pep_start % 3 == 1:

                pstart = len_previous_part + self.start-subprotpart.start + 1
                p_end = pstart + len(self.bp_obj) - 1

                if p_end <= bp_prot.location.end:
                    if self.bp_obj.strand == -1:
                        pstart, p_end = len(prot)-p_end+1, len(prot)-pstart+1
                    print(pstart, p_end)
                    print('postion aa', (pstart-1)/3+1, (p_end-2-1)/3+1)
                    input()
                    self.position_prot_relative[prot.protein_id] = (
                        int((pstart-1)/3+1), int((p_end-2-1)/3+1))
                    break
            # In some genome the cleavage site is located at the location of the ribosomal slippage
            # Then the first codon of the site is in frame of the next subprotpart but is not included in this part
            # and consequently if above don't catch it
            # I test here the end codon if it is included in the part and is in frame
            elif subprotpart.start <= self.end <= subprotpart.end:
                pstart = len_previous_part + self.start-subprotpart.start + 1
                p_end = pstart + len(self.bp_obj) - 1

                if p_end <= bp_prot.location.end:
                    self.position_prot_relative[prot.protein_id] = (
                        int((pstart-1)/3+1), int((p_end-2-1)/3+1))
                    # print(pstart, p_end)
                    # print('postion aa', (pstart-1)/3+1, (p_end-2-1)/3+1)
                    # input()
                    break

            len_previous_part += len(subprotpart)

    def getNameToDisplay(self):
        return str(self.number)


class Predicted_peptide(Peptide):
    COUNTER = 0

    def __init__(self, cds, start_in_prot, end_in_prot):

        self.start_in_prot = start_in_prot
        self.end_in_prot = end_in_prot
        self.strand = cds.bp_obj.location.strand
        # HERE ONLY SIMPLE CASE WHERE NO SHIFT IN SEQUENCE
        # NEED PARTS IN LOCATION WHEN THERE IS A SHIFFTING
        # location = partial_location + \
        #     FeatureLocation(part.start, new_end, strand=strand)
        self.location = getGenomicLocation(cds, start_in_prot, end_in_prot)

        self.start = self.location.start
        self.end = self.location.end
        self.bp_obj = SeqFeature(self.location,
                                 strand=self.location.strand,
                                 type="mat_peptide",
                                 qualifiers={"note": 'predicted mature peptide'})

        self.polyproteins = {cds}

        self.position_prot_relative = {}
        Predicted_peptide.COUNTER += 1
        self.number = Predicted_peptide.COUNTER
        # self.start_aa(cds)

        # print('real po', start_in_prot, end_in_prot)

        # print('fct', self.position_prot_relative[cds.protein_id])
        # assert self.start_aa(cds) == start_in_prot, self.end_aa(cds) == end_in_prot

    def __str__(self):
        return 'predicted peptide object'


class UnannotatedRegion(Peptide):

    def __init__(self, location, protein, strand):
        self.start = location.start
        self.end = location.end
        self.polyproteins = [protein]  # this region belongs to this protein

        self.non_overlapping_prot = set()
        self.position_prot_relative = {}

        # Domains attr related: key=domain value len of the domain on the peptide
        self.included_domains = {'fully': [], "partially": []}
        self.overlapped_by_domain = {'left': [], "right": []}

        Peptide.COUNTER += 1
        self.number = Peptide.COUNTER

        self.bp_obj = SeqFeature(location,
                                 strand=strand,
                                 type="unannotated_region",
                                 qualifiers={})


class CleavageSite(Peptide):
    COUNTER = 0

    def __init__(self, location, proteins, peptide, taxon_id, type):
        self.location = location
        self.start = location.start
        self.end = location.end
        self.proteins = proteins  # Set of protein
        self.peptides = {peptide}
        # cleavage site is created from the end of a peptide (the peptide is then on the left of the cleavage site)
        # Or it is created from the start of a peptide, the peptide is then on the right of the cleavage site.
        self.peptide_positions = {'left': set(), 'right': set()}
        side = 'right' if type == 'start' else 'left'
        self.peptide_positions[side].add(peptide)  # method .update() fills also the dict..

        self.position_prot_relative = {}
        self.taxon_id = taxon_id
        self.bp_obj = SeqFeature(location,
                                 type="cleavage site",
                                 strand=location.strand,
                                 qualifiers={})

        # NOTIFY PROTEINS THAT THEY hAVE A CLEAVAGE SITE
        for poly in proteins:
            poly.cleavage_sites.append(self)

        # Domain annotations :
        self.overlapping_domains = {}  # dict of overlapping domains with left and right distance

        CleavageSite.COUNTER += 1
        self.number = CleavageSite.COUNTER

        # alignment analysis
        self.start_in_aln = None
        self.quality = 0
        self.neighboring_sites = set()
        self.cds_of_aln = []

    def __eq__(x, y):
        return x.__key() == y.__key()

    def __str__(self):
        str_cs = f'Cleavage site {self.number}: {self.start} => {self.end} ({len(self)}nt)'
        for cds in self.proteins:
            str_cs += f'Protein{cds.number}: {self.start_aa(cds)} => {self.end_aa(cds)}\n'
        str_cs += f'And have been made from peptide {[p.number for p in self.peptides]}\n'
        return str_cs

    def peptide_composition(self):
        compo = set()
        for pep in self.peptides:
            compo.add(pep.bp_obj.type)
        return '|'.join(compo)

    def update(self, type, pep, polyproteins):

        # new proteins that haven't been notify of the cleavage site
        polyproteins = polyproteins - self.proteins

        # NOTIFY PROT
        for poly in polyproteins:
            poly.cleavage_sites.append(self)
        # cleavage site is created from the end of a peptide (the peptide is then on the left of the cleavage site)
        # Or it is created from the start of a peptide, the peptide is then on the right of the cleavage site.
        side = 'right' if type == 'start' else 'left'
        self.peptide_positions.setdefault(side, set()).add(pep)
        self.peptides.add(pep)
        self.proteins |= polyproteins

    def __key(self):
        # For the moment only taxon id is used. Normally they have different taxon id
        return (self.start, self.end, self.taxon_id)

    def __hash__(self):
        return hash(self.__key())

    # def extractSeq(self, organism, taxon, record, genetic_code, window_step):
    #     #window in aa to extract
    #     sites = []
    #     # proteins = sorted(list(self.proteins), key=lambda x: x.polyprotein_number, reverse=True)
    #
    #     for protein in proteins:
    #         # print(protein.polyprotein_number)
    #         position_prot_relative = int((self.position - protein.start)/3) #posiiton starting at 0
    #         if (self.position - protein.start)%3 != 0:
    #             logging.warning('cleavage_site seems to not point the first nt of a codon '+str(self))
    #
    #         # print(position_prot_relative)
    #         header = taxon+'|protein_{}'.format(protein.polyprotein_number)
    #         seq =  protein.getSequenceRecord('', header, record, genetic_code, (position_prot_relative-window_step,position_prot_relative+window_step +1))
    #         # print(dir(seq))
    #         seq.description = ''
    #         sites.append(seq)
    #         # seq.id = header #seq.id +'prot'+str(protein.number)
    #     for i, seq in enumerate(sites):
    #         for seq_next in sites[i+1:]:
    #             if seq_next.seq != seq.seq:
    #                 logging.warning("A same cleavage site give different sequences")
    #
    #     return seq
    #         # print(position_prot_relative/3)


class PredictedCleavageSite(CleavageSite):
    COUNTER = 0

    def __init__(self, start_in_prot, end_in_prot, cds, group, confidence_score):
        self.start_in_prot = int(start_in_prot)
        self.end_in_prot = int(end_in_prot)
        self.protein = cds  # Set of protein
        self.cleavage_site_group = group
        self.confidence_score = confidence_score
        # self.bp_obj = SeqFeature(location,
        #     type="cleavage site",
        #     strand=location.strand,
        #     qualifiers={})
        # NOTIFY PROTEINS THAT it has a predicted CLEAVAGE SITE
        cds.predicted_cleavage_sites.append(self)

        PredictedCleavageSite.COUNTER += 1
        self.number = PredictedCleavageSite.COUNTER
        # alignment analysis
        self.start_in_aln = None

        # Get positions genome relative
        # = self.get_genome_coordinates(cds)
        self.location = getGenomicLocation(cds, start_in_prot, end_in_prot)
        # assert len(self.location) % 3 == 0
        self.start = self.location.start
        self.end = self.location.end
        self.bp_obj = SeqFeature(self.location,
                                 strand=self.location.strand,
                                 type="misc_feature",
                                 qualifiers={"note": 'predicted cleavage site'})

    def get_genome_coordinates(self, cds):
        # print(cds.bp_obj.location)
        if len(cds.bp_obj.location.parts) > 1:
            print(cds.bp_obj.location)
            a = ''' # different parts in protein
            # mat peptide that overlap should also have more that one part
            # coordinate for cleavage site should take the different part intoo account
            # with the phase.
            # For now simple case... '''
            input(a)
        start = self.start_in_prot*3 + cds.start
        end = start + self.end_in_prot*3
        return start, end

    def start_aa(self, cds):
        if cds != self.protein:
            raise Error(
                'the cds of the domain and the cds given to the match object to retrieve its position are different')
        return self.start_in_prot

    def end_aa(self, cds):
        if cds != self.protein:
            raise Error(
                'the cds of the domain and the cds given to the match object to retrieve its position are different')
        return self.end_in_prot


class Match(Sequence):
    def __init__(self, taxid, seqid, method, start, end, score, Dbxref, Name, match_id, signature_desc):
        self.taxid = taxid
        self.seqid = seqid
        self.method = method
        self.start_in_prot = int(start)
        self.end_in_prot = int(end)
        self.score = score
        self.Dbxref = Dbxref
        self.name = Name
        self.match_id = match_id
        self.signature_desc = signature_desc
        self.protein = None

        self.duplicated = False

        self.overlapping = False
        self.partially_included_in = {}  # all the peptides that have a part of the domains annotation in their seq
        self.belongs_to_peptides = []
        self.right_overlaps_peptide = 0  # distance overlap on the right with the peptide
        self.left_overlaps_peptide = 0  # distance on the left
        self.right_overlap = {}
        self.left_overlap = {}

        self.overlapped_cleavage_sites = {}

        self.non_overlapping_prot = set()

    def __len__(self):
        return self.end - self.start + 1

    def __str__(self):
        string = '\n==Domain {} from {} and {}\n'.format(
            self.name, self.start_in_prot, self.end_in_prot)
        return string

    def start_aa(self, cds):
        if cds != self.protein:
            raise Error(
                'the cds of the domain and the cds given to the match object to retrieve its position are different')
        return self.start_in_prot

    def end_aa(self, cds):
        if cds != self.protein:
            raise Error(
                'the cds of the domain and the cds given to the match object to retrieve its position are different')
        return self.end_in_prot

    def isEqualTo(self, other_match):
        # If 2 match has same genomic start and end and have teh same interprot id they are duplicated
        # they would be only different in the seqid attr and potentially in the protein coordinate
        if self.start == other_match.start and self.end == other_match.end and self.Dbxref and other_match.Dbxref:
            return True

        return False

    def getNameToDisplay(self):

        if self.duplicated:
            name = 'd'*len(self.name)
        else:
            name = self.name
        return name if not self.overlapping else '<' + name

    def OverlappingDistance(self):
        return self.left_overlaps_peptide + self.right_overlaps_peptide
