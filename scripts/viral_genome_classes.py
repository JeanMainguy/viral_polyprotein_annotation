#!/usr/bin/env python3

import taxonomy as tax
import os, gzip, logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, csv, re, collections
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation
from operator import attrgetter

SCREEN_SIZE = 200


class Genome:
    def __init__(self, gb_file):
        self.gb_file = gb_file
        self.segments = []
        self.matchs = []
        self.taxonomy = None
        self.organism = None
        self.taxon_id = None
        self.expectation_node = None
        self.peptide_expectation = None
        self.polyprotein_expectation = None
        self.variable_polyprotein_expectation = None
        self.expectation_info = {}
        Protein.COUNTER = 0
        Peptide.COUNTER = 0

    def __len__(self):
        return sum([len(s) for s in self.segments])


    def extractObjectFromGbff(self, gb_file):

        with gzip.open(gb_file, "rt") as handle:
            for i, record in enumerate(SeqIO.parse(handle, "genbank")):
                segment = Segment(record)
                genome.segments.append(segment)


    def hasEnoughPeptide(self):
        total_final_peptides = sum([len(segment.peptides)-len(segment.parent_peptides) for segment in self.segments])
        return True if self.peptide_expectation <=  total_final_peptides else False

    #
    # def getTaxonExpectation(self, taxon_expectation):
    #     """
    #     We presume that the taxonomy is homogenous in the differents segments of a genome
    #     For one genome it is not the case but we ignore this case
    #     """
    #     # print('TAXON ....')
    #     if not self.segments:
    #         logging.warning( 'ERROR no segment stored in genome')
    #         return 'ERROR no segment stored in genome'
    #     segment = self.segments[0]
    #     record = segment.record
    #
    #
    #     taxonomy = record.annotations['taxonomy']
    #     self.taxonomy = taxonomy
    #     self.organism = record.annotations['organism']
    #     self.taxon_id = segment.taxon_id
    #     # print("TAXONOMY")
    #     # print(taxonomy)
    #     self.expectation_node='Not_Found'
    #     # self.peptide_expectation = 1
    #     for node in taxonomy:
    #
    #         if node in taxon_expectation:
    #             self.expectation_node = node
    #             self.peptide_expectation = taxon_expectation[node]['peptide']
    #             self.polyprotein_expectation = taxon_expectation[node]['polyprotein']
    #             self.variable_polyprotein_expectation = taxon_expectation[node]['variable_polyprotein_expectation']
    #             self.expectation_info = taxon_expectation[node]
    #     # print(len(self.segments)," == ",self.expectation_info['segment'])
    #     if not  self.expectation_info: # if the genome taxonmy has no match with any taxon where we expect peptides then we stop here
    #         return
    #     if len(self.segments) > 0 and len(self.segments) == self.expectation_info['segment']:
    #         for i in range(self.expectation_info["segment"]):
    #             # BIG approximation here! the segment in the genbank file need to be in the same order all the time
    #             # WARNING
    #
    #             nb_polyprotein_key = 'segment{}_polyprotein'.format(i+1)
    #             nb_peptide_key = 'segment{}_peptide'.format(i+1)
    #             # print(nb_polyprotein_key,  self.expectation_info[nb_peptide_key])
    #             self.segments[i].polyprotein_expectation = self.expectation_info[nb_polyprotein_key]
    #             self.segments[i].peptide_expectation = self.expectation_info[nb_peptide_key]
    #             # print( self.segments[i].polyprotein_expectation)
    #     # print('expectation node', self.expectation_node)
    #     # print(taxonomy)


    def associateMatchWithPolyprotein(self):
        for segment in self.segments:
            for poly in segment.polyproteins:
                # poly.matchs = [match for match in self.matchs if match.seqid == poly.protein_id]
                for match in self.matchs:
                    if match.seqid == poly.protein_id:
                        poly.matchs.append(match)
                        match.protein = poly
                        match.getGenomicPositions(poly.start)
                        segment.matchs.add(match)
            segment.getDomainOverlappingInfo()
            segment.identifyDuplicatedMatch()


    # def identifyExpectedElement(self):
    #     #The question is when their is more than 1 segment with more than 1 polyprotein...
    #     for segment in self.segments:
    #         segment.identifyExpectedPolyprotein(self.polyprotein_expectation, self.variable_polyprotein_expectation)
    #
    #
    # def isGenomeAnnotationRelevant(self):
    #     return False if any([s for s in self.segments if not s.isSegmentAnnotationRelevant() ]) else True


class Segment:

    POSSIBLE_TYPE = set()

    def __init__(self, record, gb_file):
        self.record = record
        self.gb_file = gb_file
        self.source = None
        self.polyproteins = set() # polyproteins list
        self.parental_polyproteins = []
        self.peptides = set() # mat_peptides list
        self.matchs = set()
        self.parent_peptides = set()
        self.sub_peptides = set()
        self.cds = set()
        self.unannotated_region = set()
        self.taxon_id = None
        self.organism = record.annotations['organism']
        self.cleavage_sites = [] # list of cleavage_site object
        self.relevant_annotation = False

        self.polyprotein_expectation = None



    def __len__(self):
        return len(self.source)

    #
    # def identifyExpectedPolyprotein(self, polyprotein_expectation, variable_polyprotein_expectation):
    #     # print("self.polyprotein_expectation", self.polyprotein_expectation)
    #     if self.polyprotein_expectation is not None:
    #         # print("self.polyprotein_expectation", self.polyprotein_expectation)
    #         polyprotein_expectation = self.polyprotein_expectation
    #         # print("polyprotein_expectation", polyprotein_expectation)
    #     if not polyprotein_expectation:
    #         # print(self.polyprotein_expectation)
    #
    #         return
    #     # if the segment has it own polyprotein_expectation then we take them as reference
    #
    #
    #     # print(polyprotein_expectation, 'Expected polyprotein in this genome')
    #     parental_prot =  [ prot for prot in self.cds if not prot.parental_prot] # if this prot is not included in another one
    #     if variable_polyprotein_expectation and len(parental_prot) > 1:
    #         parental_prot = sorted(parental_prot, key=lambda x: len(x), reverse=False)
    #         length_difference = [len(parental_prot[i]) - len(parental_prot[i-1]) for i in range(len(parental_prot))[1:]]
    #         index_boundary  = length_difference.index(max(length_difference))
    #         self.parental_polyproteins = parental_prot[index_boundary+1:]
    #
    #
    #     elif len(parental_prot) < polyprotein_expectation:
    #         logging.warning('Not enough parental protein found {}/{}'.format(len(parental_prot), polyprotein_expectation))
    #         self.parental_polyproteins = parental_prot
    #     elif len(parental_prot) == polyprotein_expectation:
    #         self.parental_polyproteins = parental_prot
    #     elif len(parental_prot) > polyprotein_expectation:# rule: take the biggest protein until the number match the expectation and consider them as polyprotein
    #         self.parental_polyproteins = sorted(parental_prot, key=lambda x: len(x), reverse=True)[:polyprotein_expectation]
    #
    #     self.parental_polyproteins = sorted(self.parental_polyproteins, key=lambda x: x.start, reverse=False)
    #     #giving a positional number to the polyprotein
    #     #if the polyprotein is included it get the positional number of the parent + a decimal number between 0.1
    #
    #
    #     # self.polyproteins = []
    #     for i, p_poly in enumerate(self.parental_polyproteins):
    #         p_poly.polyprotein_number = i+1.0
    #         self.polyproteins.add(p_poly)
    #         # print(p_poly)
    #         for j, sub_poly in enumerate(sorted(p_poly.sub_prot, key=lambda x: len(x), reverse=True)):
    #
    #             sub_poly.polyprotein_number = (i+1)+(j+1)/10
    #             self.polyproteins.add(sub_poly)
    #             # print(sub_poly)
    #             if j > 9: logging.error('polyprotein_number is not working') # subpolyprotein > 9
    #

    # def getSubProteins(self):
    #     return [p for p in self.cds if p.parental_prot ] # return protein that have parental prot


    # def getPolyWithFeature(self, feature):
    #     return [p for p in self.polyproteins if getattr(p, feature)]



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
                prot_obj = Protein(feat)
                self.cds.add(prot_obj)


    def associatePepWithProt(self, sp_treshold=90):
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
                cds.ProteinCoverage()
                self.unannotated_region.update(cds.unannotated_region)
                cds.polyprotein = True
                cds.checkForSignalP(sp_treshold)


        ##SMALL check up to be sure that every peptides have been assigned to at least one protein
        for pep in self.peptides:
            if not pep.polyproteins:
                logging.warning("The peptide {} {} was not assigned to any cds in {}".format(pep.number, str(pep.bp_obj.location), self.taxon_id))


    def checkSubPeptides(self):
        #Check if peptide is included in a bigger peptide
        for i, pep in enumerate(self.peptides):
            for pep_next in self.peptides:
                #try to determine if pep is included in pep_next
                if pep == pep_next or pep.location.strand != pep_next.location.strand:
                    continue

                if pep_next.start <= pep.start <= pep_next.end and pep_next.start <= pep.end <= pep_next.end:

                    if not(pep_next.start%3 == pep.start%3 and pep.end%3 == pep_next.end%3):

                        logging.warning('Peptide {} is included in peptide {} but does not share the same frame'.format(pep.number, pep_next.number))

                        continue

                     # and pep_next.start%3 == pep.start%3 and pep.end%3 == pep_next.end%3:
                    # print('\n', pep.bp_obj,'\nin\n', pep_next.bp_obj, '\n' )
                    self.parent_peptides.add(pep_next)

                    pep_next.parent_peptide = True
                    self.sub_peptides.add(pep)
                    # print(pep_next)
                    # print(pep)


    def checkPeptideRedundancy(self):
        #Remove Peptide that have similat start and end...
        peptides = list(self.peptides)

        for i, pep in enumerate(peptides):

            for pep_next in peptides[i+1:]:

                # if pep.location.start ==  pep_next.location.start and pep.location.end ==  pep_next.location.end:
                if  pep.location == pep_next.location:
                    self.peptides.remove(pep_next)
                    pep.redundant_pep.append(pep_next)

    # def checkForSlippage(self):
    #     """
    #     Find the protein that are included in another protein due to ribosomal_slippage are alternative start
    #     ribosomal_slippage attr is a dico with as a key the position of the end of the part of the protein and as a key the shift
    #     """
    #
    #     for poly in self.cds:
    #         # print(poly)
    #         # print(poly.bp_obj.location)
    #         if len(poly.bp_obj.location.parts) == 1:
    #             # print('PROTEIN HAS ONLY ONE PART')
    #             continue
    #         first_part = poly.bp_obj.location.parts[0]
    #         # print(first_part)
    #         for part_next in poly.bp_obj.location.parts[1:]:
    #             # print(part_next)
    #             # print('CHECK NEXT PART')
    #             shift = part_next.start - first_part.end
    #             # print('SHIFT of ',shift)
    #
    #             if abs(shift) < 3:
    #
    #                 poly.ribosomal_slippage[first_part.end] = shift
    #             elif shift == 3:
    #                 poly.readthrough[first_part.end] = shift
    #                 logging.info('Shift of 3.. Should be a readthrough but better to check {}'.format(self.gb_file))
    #             else:
    #
    #                 logging.info('shift of {}, probably a spliced gene  {}'.format(shift, self.gb_file))
    #
    #             first_part = part_next

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

    #
    # def writeAnnotatedProteins(self, file_handle, genetic_code):
    #     for polyprotein in self.polyproteins:
    #         if polyprotein.peptides: # we write every protein that have at least one peptide annotation
    #             seq_to_write = polyprotein.getSequenceRecord(self.organism, self.taxon_id, self.record, genetic_code)
    #             SeqIO.write(seq_to_write, file_handle,"fasta")


    # def writeIdentifiedPolyproteins(self, file_handle, genetic_code):
    #     for polyprotein in self.polyproteins:
    #         seq_to_write = polyprotein.getSequenceRecord(self.organism, self.taxon_id, self.record, genetic_code)
    #         SeqIO.write(seq_to_write, file_handle,"fasta")
    #

    # def writeCleavageSite(self, file_handle, genetic_code, window_step):
    #     for cleavage_site in self.cleavage_sites:
    #         site_seq = cleavage_site.extractSeq(self.organism, self.taxon_id, self.record, genetic_code, window_step) #extractSeq(self, organism, taxon, record, genetic_code, window_step)
    #         # print(site_seq)
    #         SeqIO.write(site_seq, file_handle, "fasta")


    def getCleavageSites(self):

        start_cleavage_sites = {}

        for pep in sorted(list(self.peptides | self.unannotated_region), key=lambda x: x.start, reverse=False): #(self.peptides | self.unannotated_region):  # sorted(list(self.peptides | self.unannotated_region), key=lambda x: x.start, reverse=False):
            strand = pep.bp_obj.strand

            for type, border in [("start", pep.start), ('end', pep.end)]:

                start = border - 3
                end = border + 3

                # print(type, border)
                # print("CS position", start, end)
                polyproteins = {poly for poly in pep.polyproteins if poly.start+9 < start < poly.end-9} #polyprot that are compatible with the cleavage site
                if polyproteins:
                    if start in start_cleavage_sites:
                        site = start_cleavage_sites[start]

                        site.update(type, pep, polyproteins) # this cleavage site has been already treated
                        continue
                    partial_location = False
                    previous_part_end = False
                    overlap = False
                    for poly in polyproteins:
                        for part in poly.bp_obj.location.parts:

                            if overlap and end in part:
                                if partial_location is not False:
                                    # print('WE are in case where border was in the previous part')
                                    # input()
                                    len_cs_previous_part = len(partial_location)
                                    len_cs_actual_part = 6 - len_cs_previous_part
                                    new_end = part.start + len_cs_actual_part
                                    location = partial_location + FeatureLocation(part.start, new_end, strand=strand)
                                    tuple_position = (start, new_end)

                                elif previous_part_end is not False:
                                    # print('WE are in case where border is in the current part')
                                    # input()
                                    partial_location = FeatureLocation(part.start, end, strand=strand)
                                    len_previous_part = 6 - len(partial_location)
                                    new_start = previous_part_end - len_previous_part
                                    location = FeatureLocation(new_start, previous_part_end, strand=strand) + partial_location
                                    tuple_position = (new_start, end)


                                if tuple_position in start_cleavage_sites:
                                    site = start_cleavage_sites[tuple_position]
                                    site.update(type, pep, {poly}) # this cleavage site has been already treated
                                else: # we create the CS object
                                    cleavage_site = CleavageSite(location, polyproteins, pep, self.taxon_id)
                                    start_cleavage_sites[tuple_position] =  cleavage_site # we store the start and end of the CS in the dico to update it and no recreate when it is nececessary
                                    # print(cleavage_site)
                                    # print(location, len(location))
                                    # input()
                                    assert(len(location) == 6), 'Cleavage site has a length of {} location {}'.format(len(location), location)
                                partial_location = False
                                previous_part_end = False
                                overlap = False
                                break # we break and check a new polyprotein

                            if start in part and end in part: # cleavage site is included in the prot part
                                if (start, end) in start_cleavage_sites:
                                    site = start_cleavage_sites[(start, end)]
                                    site.update(type, pep, {poly}) # this cleavage site has been already treated
                                else: # creation of the cleavage site object
                                    location = FeatureLocation(start, end, strand=strand)
                                    cleavage_site = CleavageSite(location, polyproteins, pep, self.taxon_id)
                                    start_cleavage_sites[(start, end)] =  cleavage_site
                                break # we break and check a new polyprotein

                            if start in part and end not in part:
                                overlap = True
                                if border-1 in part: # the cleavage site overlap two part of the protein
                                    # print("border", border, 'in part', part)
                                    partial_location = FeatureLocation(start, part.end, strand=strand)
                                    # print('modulo start prot part',part.start%3)
                                    # print('modulo start CS', start%3)
                                    # print("partial location", partial_location, len(partial_location))
                                    # # we check the next part of the protein -->

                                    continue

                                else:
                                    # print("border", border, 'not in part', part)
                                    previous_part_end = part.end



                    # print("Cleavage site:", type, border,' start end', start, end, '\t',  cleavage_site.start_aa(poly), cleavage_site.end_aa(poly))
                    start_cleavage_sites[start] =  cleavage_site
                    self.cleavage_sites.append(cleavage_site)

    def isSegmentAnnotationRelevant(self):
        return False if any([p for p in self.polyproteins if not p.isAnnotationRelevant() ]) else True


class Sequence:

    def __str__(self):

        string = 'Sequence: from {} to {}'.format(self.start,self.end)
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
        #for uncovered region and match to get the genomic position and not the prot relative position
        self.start = prot_start -1 + self.start_in_prot*3 -2 -1
        self.end = prot_start + self.end_in_prot*3 -1


    def getNameToDisplay(self):
        return 'None'


class Protein(Sequence):

    COUNTER=0

    def __init__(self, biopyth_obj):
        # self.start = start
        # self.end = end

        self.bp_obj = biopyth_obj
        self.start, self.end = (biopyth_obj.location.start, biopyth_obj.location.end)

        self.gene_type = biopyth_obj.type
        self.protein_id = "Unknown" if 'protein_id' not in biopyth_obj.qualifiers else biopyth_obj.qualifiers['protein_id'][0]
        self.product =  "Unknown" if 'product' not in biopyth_obj.qualifiers else biopyth_obj.qualifiers['product'][0]
        self.peptides = set()
        self.matchs = []
        self.unannotated_len = len(biopyth_obj)/3
        self.polyprotein = False
        self.non_polyprotein_explanation = ''
        self.unannotated_region = []
        self.non_overlapping_prot = set() #for the visualisation
        # self.ribosomal_slippage = {}
        # self.readthrough = {}
        self.parental_prot = []
        self.sub_prot = []
        # self.alternative_start = False
        self.cleavage_sites = []
        self.polyprotein_number = 0.0
        self.sequence = ''

        Protein.COUNTER+=1
        self.number = Protein.COUNTER


    def __len__(self):
        return len(self.bp_obj)


    def __str__(self):

        # status += str(self.polyprotein_number)
        string = 'Protein {}: {}nt | {}aa | {} \n'.format(self.number,len(self), len(self)/3, self.bp_obj.location)
        string += '  Peptide: {}\n  UnannotatedRegion {} \n'.format(len(self.peptides), len(self.unannotated_region))
        # for p in self.peptides:
        #     string += str(p) +"\n"
        # string += "Position_number {}\n".format(self.polyprotein_number)
        # string += "{}:{}\n".format(self.gene_type, ' and '.join(self.reasons))


        # string += '{} annotated peptides | {} covered by peptide\n'.format(len(self.peptides), 100 - (self.unannotated_len/(len(self)/3))*100) if self.polyprotein else ''
        # string += '{} domain annotations\n'.format(len(self.matchs)) if self.polyprotein else ''
        # string += '\n'.join([str(p) for p in sorted(list(self.peptides), key=lambda x: x.start, reverse=False)])
        # print('ddd')
        return string


    def __contains__(self, pep):
        ##test if a peptide is in polyprotein
        # the qualifiers gene is no more used to know if a peptide belong to a polyprotein
        # Only the position of the peptide tel us if the peptide belongs to a poly
        # Start < end always. The coordinate are given according the strand of the gene
        # pep = seq.bp_obj # extract biopython object of peptide objet
        # seq_location = seq if not seq.__class__.__name__ ==  'Peptide' else seq.bp_obj.location
        # print('---------------------------')
        # print('PEPTIDE: ', pep.number, pep.location)
        #
        # print('Protein: ', self.number, self.bp_obj.location)

        # print('strand',  pep.location.strand ==  self.bp_obj.location.strand)
        # print("start", pep.location.parts[0].start , pep.location.parts[0].start in self.bp_obj.location)
        # print('end',pep.location.parts[-1].end,pep.location.parts[-1].end in self.bp_obj.location )
        # pep.location.parts[0].start need to use part to get the real start position and the min and max position in the sequence
        # otherwise when the sequence is circular like in Woodchuck hepatitis virus|35269
        pep_start = pep.location.parts[0].start
        pep_end = pep.location.parts[-1].end

        if not (pep.location.parts[0].start in self.bp_obj.location and pep.location.parts[-1].end-1 in self.bp_obj.location and pep.bp_obj.strand ==  self.bp_obj.strand):
            # print("direct F")
            return False

        pep_parts = iter(pep.location.parts)
        pep_part = next(pep_parts)
        case2 = False
        shift = 0
        protein_start =  self.bp_obj.location.parts[0].start
        previous_end = protein_start
        prot_len = 0
        for sub_location in self.bp_obj.location.parts:
            # shift += sub_location.start - previous_end

            # print('\n==Protein Sub location==', sub_location)
            #
            # print('Shift', shift)
            # print('pep part', pep_part)
            # print('prot part', sub_location)
            if case2:
                # print('\n-WE are in case 2', )

                if pep_part.start == sub_location.start: # check if the 2 parts start are the same. So if the join() is similar
                    # print('--the start are the same')
                    # print("--pep parts", pep_parts)
                    try:
                        pep_part = next(pep_parts)
                    except StopIteration:
                        if pep_end in sub_location: # No more part for the peptide and it end is included in the subprot prat
                            # print('---end is included in the prot part')
                            # print('TRUUE')
                            return True
                        else:
                            return False
                    else:
                        # print('--We go to next part')

                        continue
                else: #The peptide in somehow not folowing the prot part
                    # print('-The peptide in somehow not folowing the prot part')
                    # print('FALSE')
                    return False
            # Case 1 : The peptide is included in one part of the protein
            if pep_start in sub_location and pep_end in sub_location and len(pep.location.parts) == 1 and pep_start%3 == (sub_location.start - shift)%3:
                # print("-the peptide is included in the subprot part")
                return True

            # print(" pep_part.start%3 ",  pep_part.start%3 )
            # print("(sub_location.start - shift)%3", sub_location.start, '-' ,shift,')%3', (sub_location.start - shift)%3)
            # Case 2 Peptide is overlapping this part with the next one
            if pep_start in sub_location and pep_end not in sub_location and pep_part.end == sub_location.end and pep_part.start%3 == (sub_location.start - shift)%3:
                try:
                    pep_part = next(pep_parts)
                except StopIteration:
                    return True # the peptide end where the ribosomal shifting occurs
                case2 = True
                # print('-the peptide is overlapping the prot part')

            previous_end = sub_location.end
            prot_len += len(sub_location)
            shift = prot_len%3
        # print("False ..")
        ## SPECIAL CASE: intein found in some CDS of dsDNA (example Phicbkvirus)
        if pep.realStart() == self.realStart() and pep.realEnd() == self.realEnd() -3*self.bp_obj.strand and len(pep) < len(self)-3 and len(pep.location.parts)> 1:
            return True


        return False

    def checkForSignalP(self, sp_treshold):
        cds_start, cds_end = (self.start, self.end) if self.bp_obj.strand == 1 else (self.end, self.start)

        if  0 < len(self.peptides) <= 2:
            for pep in self.peptides:
                pep_start, pep_end = (pep.start, pep.end) if self.bp_obj.strand == 1 else (pep.end, pep.start)
                if cds_start == pep_start and len(pep)/3<sp_treshold:
                    self.polyprotein = False
                    self.non_polyprotein_explanation = "Signal Peptide"

            if len(self.peptides) == 1 and pep_end == cds_end-3:
                if pep_start < cds_start + sp_treshold*3:
                    # potential protein with signal peptide where only the mature peptide is annotated
                    # it happens that the whole sequence is covered by a signle mat_peptide with the exception of the first and last condon
                    # example: 11886 Rous sarcoma virus	Viruses;Retro-transcribing viruses;Retroviridae;Orthoretrovirinae;Alpharetrovirus
                    self.polyprotein = False
                    if cds_start <= pep_start <= cds_start +3:
                        self.non_polyprotein_explanation = "single mat_peptide covering the whole CDS"
                    else:
                        self.non_polyprotein_explanation = "single mat peptide covering almost the whole CDS"


        #Check for intein
        if len(self.peptides) == 2:

            # sort by start and select the first peptide no matter if the strand is -1
            # because in case of intein the peptide cover the begining and the end of the CDS
            pep, pep_middle = sorted(list(self.peptides), key=lambda x: x.start, reverse=False)
            # to identify a protein with an intein we use the mat peptide annotation that flank the mat peptide of the intein
            # This peptide annotation covers the begining of the CDS and the end, it has then 2 parts and its length is smaller than the CDS.

            if pep.realStart() == self.realStart() and pep.realEnd() == self.realEnd() -3*self.bp_obj.strand:
                if len(pep) < len(self)-3 and len(pep.location.parts)> 1:

                    self.polyprotein = False
                    self.non_polyprotein_explanation = "Intein outline extein suround intein"
                # in some genome the mature peptide of the intein cover all the CDS and the intein peptide a small part
                if  len(pep) == len(self)-3 and  pep_middle.start in pep.location and pep_middle.end in pep.location :
                    self.polyprotein = False
                    self.non_polyprotein_explanation = "Intein outline extein include intein"



        # if self.non_polyprotein_explanation == "":
        #     print(self.non_polyprotein_explanation)
        #     print(self)
        #     print(self.protein_id)
        #     for p in self.peptides:
        #         print(p)
        #     print("non_polyprotein_explanation:", self.non_polyprotein_explanation)
        #     # input()

    def getSequenceAA(self, record, genetic_code):

        if self.sequence:
            return self.sequence

        elif 'translation' in self.bp_obj.qualifiers:
            seq = self.bp_obj.qualifiers["translation"][0]
        else:
            seq = self.bp_obj.location.extract(record).seq.translate(table=genetic_code)
            seq = str(seq)
        seq = seq.upper()

        for site in self.cleavage_sites:
            site_position = site.start_aa(self)-1 # -1 to be in base 0
            # print('SITE\n',site)
            # print('location',site.bp_obj.location)
            #
            # print('PROT\n', self)
            # print('location',self.bp_obj.location)

            site_extraction = site.bp_obj.location.extract(record).seq.translate(table=genetic_code)
            site_seq =  str(self.bp_obj.location.extract(record).seq.translate(table=genetic_code, to_stop=False))[site_position:site_position+2]

            # print(site.bp_obj.location.extract(record).seq)
            # print( str(self.bp_obj.location.extract(record).seq)[(site_position-8)*3:(site_position+8)*3])

            # print(site_seq)
            # print(site_extraction)
            if site_extraction.upper() != site_seq.upper():
                # print(site_position)

                logging.warning('{}:cleavage site ({}) is not similar in extraction from genome sequence and the protein sequence: extraction:{} and prot_seq:{}'.format(self.protein_id,site.bp_obj.location, site_extraction, site_seq))
            if site_position < 0:
                logging.warning('{}:The cleavage site position in the sequence is wrong {}'.format(self.protein_id,site_position))

            seq = seq[:site_position] + seq[site_position:site_position+2].lower() + seq[site_position+2:]

            # print(seq[site_position-8:site_position+8])
        self.sequence = seq

        return seq


    # def extractCleavageSites(self, record, genetic_code, window):
    #     cleavage_sites_seq = []
    #     if self.isAnnotationRelevant():
    #         for site in self.cleavage_sites:
    #             site_position = site.peptide.end_aa(self)-1 # -1 to be in base 0
    #             print(site.peptide.bp_obj.location.extract(record).seq.translate(table=genetic_code))
    #
    #             cleavage_sites_seq.append(seq[site_position-window:site_position+window])
    #
    #             print(site_position, seq[site_position-10:site_position+11])


    def isAnnotationRelevant(self):
        #Simple rule to be sure that the annotation is relevant
        #  the number of uncvered region < number of peptide
        #To not take into acount this kind of annotation
        # ---==1.0============================================================================================================================================================================>------------------------
        #           |+5|         |+6+++|     |+7++++|
        #    |~8~~~~|  |~9~~~~~~~|     |~10~~|      |~11~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

        if len(self.unannotated_region) >= len([p for p in self.peptides if not p.parent_peptide]):
             return False
        else:
            return True


    def getSequenceRecord(self,organism, header, record, genetic_code, subPosition=False):


        seq = self.getSequenceAA(record, genetic_code)

        # print(seq)
        if subPosition:
            seq = seq[subPosition[0]:subPosition[1]]
            header += '|{}:{}'.format(subPosition[0],subPosition[1])

        return  SeqRecord(Seq(seq,generic_protein) , id=header, description="{}|{}".format(self.product, organism))
        # print("PROTTTT")
        # print(seq_to_write)
        # SeqIO.write(seq_to_write, file_handle,"fasta")


    def ProteinCoverage(self):
        ##Check if the mat peptide associated with the protein are covering all the prot or not
        #If not create unannotated_region to fill the gaps
        unannotated_position = []
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
            if current_po < pep.start: # if pep start is after the current positon then a unnatotade region need to be created

                if protpart.end < pep.start: #that means the next peptide start is located in the next part of the protein
                    location = FeatureLocation(current_po, protpart.end, strand=strand)
                    partial_location = partial_location + location if partial_location else location

                    protpart = next(iter_parts)
                    current_po = protpart.start
                    continue

                else:
                    #an unannotated position is created.
                    location = FeatureLocation(current_po, pep.start, strand=strand)
                    if partial_location:
                        location = partial_location + location
                        partial_location = None
                    if not(strand == -1 and len(location) <= 3): # in case of strand - the start is the end and then we don't want to take into account the stop codon

                        unannotated_seq = UnannotatedRegion(location, self, strand)
                        self.unannotated_region.append(unannotated_seq)
                    current_po = pep.end


            elif current_po < pep.end:
                if protpart.end < pep.end: #the peptide is overlapping two part of the prot
                    protpart = next(iter_parts) # we take the next part

                current_po = pep.end
            pep_i += 1

            # else: #Another case where curent_po > pep.end but then we don't change the current po and we move to the next pep
            #     current_po =  current_po


        if current_po < protpart.end: # we didn't not arrive at the end of the protein
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
        if start in next_location and end in next_location and start%3 == next_location.start%3:
            return True
        return False


    def checkforAlternativeStart(self, poly):

        if self.end ==  poly.end and  self.start !=  poly.start and self.start%3 ==  poly.start%3:
             self.alternative_start = True


class Peptide(Sequence):

    COUNTER = 0
    def __init__(self, bp_obj):
        #by default the start and end is in amino acid
        self.bp_obj = bp_obj
        self.location = bp_obj.location
        self.start = bp_obj.location.start
        self.end =  bp_obj.location.end
        self.redundant_pep = []
        self.non_overlapping_prot = set() #for the visualisation
        self.parent_peptide = False

        self.polyproteins = set()

        self.position_prot_relative = {}
        Peptide.COUNTER +=1
        self.number = Peptide.COUNTER
        # self.genome_start = None
        # self.genome_end = None


    def __str__(self):

        string = 'Peptide {}: from {} to {} | belongs to {}\n'.format(self.number, self.start,self.end, [prot.number for prot in self.polyproteins])
        string += '{} nt | {}aa | location :{}  \n'.format(len(self), len(self)/3, self.bp_obj.location)
        string += '  Position in protein\n'
        for prot in self.polyproteins:
            string += '    protein {} from {} to {}\n'.format(prot.number, self.start_aa(prot), self.end_aa(prot))

        return string


    def start_aa(self, prot):
        if prot.protein_id not in self.position_prot_relative:
            self.getProteinPosition(prot)

        return self.position_prot_relative[prot.protein_id][0]


    def end_aa(self, prot):

        if prot.protein_id not in self.position_prot_relative:
            self.getProteinPosition(prot)

        return self.position_prot_relative[prot.protein_id][1]


    def get_position_prot_relative(self, prot):
        if prot.protein_id not in self.position_prot_relative:
            self.getProteinPosition(prot)

        return self.position_prot_relative[prot.protein_id]


    def getProteinPosition(self, prot):
        # print("==================getProteinPosition======")
        len_previous_part = 0
        bp_prot = prot.bp_obj
        ## Searching the peptides position protein relative
        ## Due to ribosomal slippage the conversion is not trivial
        for subprotpart in bp_prot.location.parts:
            # print("subprotpart ", subprotpart)
            # print(self.start, self.start%3)
            # print(subprotpart.start, subprotpart.start%3 )

            if subprotpart.start <=  self.start <= subprotpart.end:

                pstart = len_previous_part +  self.start-subprotpart.start +1
                p_end = pstart + len(self.bp_obj) -1

                if p_end <= bp_prot.location.end:
                    if self.bp_obj.strand == -1:
                        pstart, p_end = len(prot)-p_end+1, len(prot)-pstart+1
                    self.position_prot_relative[prot.protein_id] =  (int((pstart-1)/3+1), int((p_end-2-1)/3+1))
                    break
            # In some genome the cleavage site is located at the location of the ribosomal slippage
            # Then the first codon of the site is in frame of the next subprotpart but is not included in this part
            # and consequently if above don't catch it
            # I test here the end codon if it is included in the part and is in frame
            elif subprotpart.start <=  self.end <= subprotpart.end:
                pstart = len_previous_part +  self.start-subprotpart.start +1
                p_end = pstart + len(self.bp_obj) -1

                if p_end <= bp_prot.location.end:
                    self.position_prot_relative[prot.protein_id] =  (int((pstart-1)/3+1), int((p_end-2-1)/3+1))
                    break

            len_previous_part += len(subprotpart)


    def getNameToDisplay(self):
        return str(self.number)


class UnannotatedRegion(Peptide):

    def __init__(self, location, protein, strand):
        self.start =  location.start
        self.end   = location.end
        self.polyproteins = [protein] # this region belongs to this protein

        self.non_overlapping_prot = set()
        self.position_prot_relative = {}

        Peptide.COUNTER +=1
        self.number = Peptide.COUNTER

        self.bp_obj = SeqFeature(location,
            strand=strand,
            type="unannotated_region",
            qualifiers={})


class CleavageSite(Peptide):
    COUNTER = 0
    def __init__(self, location, proteins, peptide, taxon_id):
        self.start = location.start
        self.end = location.end
        self.proteins = proteins # Set of protein
        self.peptides = {peptide}
        self.position_prot_relative = {}
        self.taxon_id = taxon_id
        self.bp_obj = SeqFeature(location,
            type="cleavage site",
            strand=location.strand,
            qualifiers={})
        #NOTIFY PROTEINS THAT THEY hAVE A CLEAVAGE SITE
        for poly in proteins:
            poly.cleavage_sites.append(self)

        CleavageSite.COUNTER +=1
        self.number = CleavageSite.COUNTER

    def __eq__(x, y):
        return x.__key() == y.__key()

    def __str__(self):
        return 'Cleavage site {} from {} to {} ({}nt) belongs to {}, And have been made from peptide {}'.format(self.number, self.start, self.end, len(self), [p.number for p in self.proteins], [p.number for p in self.peptides] )

    def update(self, type, pep, polyproteins):

        polyproteins = polyproteins - self.proteins # new proteins that haven't been notify of the cleavage site

        for poly in polyproteins:
            poly.cleavage_sites.append(self)
        self.peptides.add(pep)
        self.proteins |= polyproteins

    def __key(self):
        return (self.start, self.taxon_id) # For the moment only taxon id is used. Normally they have different taxon id

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


class Match(Sequence):
    def __init__(self, taxid, seqid,  method, start, end, score, Dbxref, Name, matchID, signature_desc):
        self.taxid =taxid
        self.seqid=seqid
        self.method = method
        self.start_in_prot = int(start)
        self.end_in_prot = int(end)
        self.score=score
        self.Dbxref=Dbxref
        self.name =Name
        self.matchID = matchID
        self.signature_desc = signature_desc
        self.protein = None

        self.duplicated = False

        self.overlapping = False
        self.partially_included_in = {} # all the peptides that have a part of the domains annotation in their seq
        self.including_peptides = []
        self.right_overlaps_peptide = 0 # distance overlap on the right with the peptide
        self.left_overlaps_peptide = 0 # distance on the left

        self.non_overlapping_prot = set()

    def __len__(self):
        return self.end - self.start +1

    def __str__(self):
        string = '\n==Domain {} from {} and {}\n'.format(self.name, self.start_in_prot, self.end_in_prot)
        # string += 'Is included in: ' + str([p.number for p in self.including_peptides]) + '\n'
        # for p in self.including_peptides:
        #     string += '  pep:'+str(p.number) +':   '+ str(p.get_position_prot_relative( self.protein) ) + '\n'
        # string += '\nIs overlaping on the right: ' + str({p.number:dist for p, dist in self.right_overlaps_peptides.items()}) + '\n'
        # for p, dist in self.right_overlaps_peptides.items():
        #     string +='  pep:'+str(p.number) +'  overlap of '  + str(dist) +':   '+str(p.get_position_prot_relative( self.protein) ) + '\n'
        #     string += "  "+str(p)+'\n'
        # string += '\nIs overlaping on the left: ' + str({p.number:dist for p, dist in self.left_overlaps_peptides.items()}) + '\n'
        # for p, dist in self.left_overlaps_peptides.items():
        #
        #     string +='  pep:'+str(p.number) +'  overlap of '  + str(dist) +':   '+str(p.get_position_prot_relative( self.protein) ) + '\n'
        #     string += "  "+str(p)+'\n'
        return string

    def isEqualTo(self, other_match):
        # If 2 match has same genomic start and end and have teh same interprot id they are duplicated
        # they would be only different in the seqid attr and potentially in the protein coordinate
        if self.start == other_match.start and  self.end == other_match.end and self.Dbxref and other_match.Dbxref:
            return True

        return False

    def getNameToDisplay(self):

        if self.duplicated:
            name = 'd'*len(self.name)
        else:
            name = self.name
        return name if not self.overlapping else '<' + name

    # def get_csv_dico(self, match_header):
    #     dico = {}
    #     for attribute_name in match_header:
    #         if hasattr(self,attribute_name):
    #             attribute = getattr(self, attribute_name)
    #
    #             if callable(attribute):
    #                 dico[attribute_name] = attribute()
    #             else:
    #                 dico[attribute_name] = attribute
    #
    #     return dico
    def OverlappingDistance(self):
        return self.left_overlaps_peptide + self.right_overlaps_peptide
