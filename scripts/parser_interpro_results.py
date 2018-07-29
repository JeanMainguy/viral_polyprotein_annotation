import re, logging
import viral_genome_classes as obj


def getMatchObject(genome, gff_file):
    #retreive all the match that belong to the genome
    flag = False
    with open(gff_file) as fl:
        for l in fl:

            if l.startswith(">") or l.startswith("##FASTA"):
                break
            if l.startswith('##'):
                flag = True if l.startswith('##sequence-region {}|'.format(genome.taxon_id)) else False
                continue
            if flag:
                (taxseqid, method, feature_type, start, end, score, strand, phase, attributes) = l.split("\t")
                #match = {"start":int(begin), "end":int(end), "app":method, "seqid":seqid, "taxid":taxid}
                if feature_type == 'protein_match':
                    (taxid, seqid, length) = taxseqid.split("|")
                    re_result = re.search("ID=match\$([\d]+)_[\d]+_[\d]+", attributes)
                    matchID = "" if not re_result else int(re_result.group(1))

                    # signature_desc
                    re_result = re.search("signature_desc=([^;]+)", attributes)
                    signature_desc= "" if not re_result else re_result.group(1)

                    # Domain name
                    re_result = re.search("Name=([\d\w]+)", attributes)
                    Name = "" if not re_result else re_result.group(1)

                    # Dbxref
                    re_result = re.search('Dbxref="([^"]+)"', attributes)
                    Dbxref= "" if not re_result else re_result.group(1)

                    genome.matchs.append(obj.Match(taxid, seqid,  method, start, end, score, Dbxref, Name, matchID, signature_desc))

def associateMatchWithPolyprotein(genome):
    for segment in genome.segments:
        for cds in segment.cds:
            # cds.matchs = [match for match in self.matchs if match.seqid == cds.protein_id]
            for match in genome.matchs:
                if match.seqid == cds.protein_id:
                    cds.matchs.append(match)
                    match.protein = cds
                    match.getGenomicPositions(cds.start)
                    segment.matchs.add(match)
        getDomainOverlappingInfo(segment)
        identifyDuplicatedMatch(segment)

    #small checking to be sure that all match have of the genome has at least one prot
    for match in genome.matchs:
        if not match.protein:
            logging.warning(f'Domain annotation {match.name} has no protein in {genome.taxon_id}')


def getDomainOverlappingInfo(segment):

    for poly in segment.cds:
        for pep in list(poly.peptides) + poly.unannotated_region:
            if getattr(pep, 'parent_peptide', False):
                continue

            pep_start = pep.start_aa(poly)
            pep_end   = pep.end_aa(poly)

            for m in poly.matchs:

                if pep_end < m.start_in_prot  or m.end_in_prot < pep_start: # the match is not on peptide seq
                    continue

                if pep_start <= m.start_in_prot and m.end_in_prot <= pep_end: # The match fall into the peptide
                    m.fully_included_in.append(pep)
                    pep.fully_included_domains[m] = (m.start_in_prot, m.end_in_prot)
                    continue
                m.overlapping = True
                # here 3 possible scenario
                # The match overlap on the left, on the right or overlap oon the left and on the right

                # We find the length of the sequence that the domain annotation share with pep
                # to later determine which peptide share the biggest sequence with annotatuon
                max_start = max(pep_start, m.start_in_prot)
                min_end = min(pep_end, m.end_in_prot)

                m.partially_included_in[pep] = min_end - max_start +1
                pep.partially_included_domains[m] = (max_start, min_end) # position of the domain inside the peptide

        for m in poly.matchs:
            #For the match that overlap 2 or more peptide
            if not m.overlapping:
                continue
            # to which peptide the match belong to? based on the length
            max_length = 0
            for pep, length in m.partially_included_in.items():

                if length > max_length:
                    max_length = length
                    peptide = pep
            m.fully_included_in.append(peptide)
            pep_start = peptide.start_aa(poly)
            pep_end   = peptide.end_aa(poly)

            if  pep_end < m.end_in_prot: # the annotation overlap on the right
                m.right_overlaps_peptide =  m.end_in_prot - pep_end

            if  m.start_in_prot < pep_start: # it overlaps on the left
                m.left_overlaps_peptide = pep_start - m.start_in_prot


def identifyDuplicatedMatch(segment):
    # Due to ribosomal_slippage 2 proteins have been given to interpro but they share a similar part and then may share same somain
    # To not count twice this domain we identify the domain that have the same start and end in the cds and have same name

    # The structure is not that efficient because the duplicated match check the other match and it is not usefull they should be remove
    matchs = list(segment.matchs)
    for i, match in enumerate(matchs):
        for match_next in matchs[i+1:]:
            if match.isEqualTo(match_next):
                match_next.duplicated = True


def getStringIndices(segment, conversion):
    # print(conversion,' start', segment.start/conversion, 'end', (segment.end-1)/conversion)
    # print(conversion,' start', int(segment.start/conversion), 'end', int((segment.end-1)/conversion))
    return (int(segment.start/conversion), int((segment.end-1)/conversion))
    # return round(segment.bp_obj.location.start/conversion), round(segment.bp_obj.location.end/conversion)
