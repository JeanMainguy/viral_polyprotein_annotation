
def getMatchObject(self, gff_file):
    #retreive all the match that belong to the genome
    flag = False
    with open(gff_file) as fl:
        for l in fl:

            if l.startswith(">") or l.startswith("##FASTA"):
                break
            if l.startswith('##sequence-region'):
                flag = True if l.startswith('##sequence-region {}|{}'.format(self.expectation_node,self.taxon_id)) else False
                continue
            if flag:
                (taxseqid, method, feature_type, start, end, score, strand, phase, attributes) = l.split("\t")
                #match = {"start":int(begin), "end":int(end), "app":method, "seqid":seqid, "taxid":taxid}
                if feature_type == 'protein_match':
                    (expectation_node, taxid, seqid, polyprotein_number) = taxseqid.split("|")
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

                    self.matchs.append(Match(taxid, seqid,  method, start, end, score, Dbxref, Name, matchID, signature_desc))


def getDomainOverlappingInfo(self):

    for poly in self.polyproteins:
        for pep in list(poly.peptides) + poly.unannotated_region:

            if getattr(pep, 'parent_peptide', False):
                continue

            pep_start = pep.start_aa(poly)
            pep_end   = pep.end_aa(poly)

            for m in poly.matchs:

                if pep_end < m.start_in_prot  or m.end_in_prot < pep_start: # the match is not on peptide seq
                    continue

                if pep_start <= m.start_in_prot and m.end_in_prot <= pep_end: # The match fall into the peptide
                    # print(m.name)
                    m.including_peptides.append(pep)
                    continue

                m.overlapping = True
                # here 3 possible scenario
                # The match overlap on the left, on the right or overlap oon the left and on the right

                # We find the length of the sequence that the domain annotation share with pep
                # to later determine which peptide share the biggest sequence with annotatuon
                max_start = max(pep_start, m.start_in_prot)
                min_end = min(pep_end, m.end_in_prot)

                m.partially_included_in[pep] = min_end - max_start +1

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
            m.including_peptides.append(peptide)
            pep_start = peptide.start_aa(poly)
            pep_end   = peptide.end_aa(poly)

            if  pep_end < m.end_in_prot: # the annotation overlap on the right
                m.right_overlaps_peptide =  m.end_in_prot - pep_end

            if  m.start_in_prot < pep_start: # it overlaps on the left
                m.left_overlaps_peptide = pep_start - m.start_in_prot


def identifyDuplicatedMatch(self):
    # Due to ribosomal_slippage 2 proteins have been given to interpro but they share a similar part and then may share same somain
    # To not count twice this domain we identify the domain that have the same start and end in the genome and have same name

    # The structure is not that efficient because the duplicated match check the other match and it is not usefull they should be remove
    matchs = list(self.matchs)
    for i, match in enumerate(matchs):
        for match_next in matchs[i+1:]:
            if match.isEqualTo(match_next):
                match_next.duplicated = True


def getStringIndices(self, conversion):
    # print(conversion,' start', self.start/conversion, 'end', (self.end-1)/conversion)
    # print(conversion,' start', int(self.start/conversion), 'end', int((self.end-1)/conversion))
    return (int(self.start/conversion), int((self.end-1)/conversion))
    # return round(self.bp_obj.location.start/conversion), round(self.bp_obj.location.end/conversion)
