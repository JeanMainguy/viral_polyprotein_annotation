# viral_polyprotein_annotation



##### Identification of annotated Polyprotein in RefSeq genome.
Base rule: all viral CSD with a mat peptide annotation are considered to be a polyprotein.
Exception with Signal Peptide

###### Exception
Problem with Caudovirales - Phicbkvirus and Cafeteriavirus.
Weird mat peptide annotation with a mat_peptide surrounding a nested mat peptide. This caused by intein annotations found in some dsDNA viruses : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2984142/.
Example:
`693272	Cafeteria roenbergensis virus BV-PW1	Viruses;dsDNA viruses, no RNA stage;Mimiviridae;Cafeteriavirus`: `/mirror/ncbi/current/genomes/refseq/viral/Cafeteria_roenbergensis_virus/latest_assembly_versions/GCF_000889395.1_ViralProj59783/GCF_000889395.1_ViralProj59783_genomic.gbff.gz`


`1211641	Caulobacter virus Karma	Viruses;dsDNA viruses, no RNA stage;Caudovirales;Siphoviridae;Phicbkvirus`

`/mirror/ncbi/current/genomes/refseq/viral/Caulobacter_virus_Rogue/latest_assembly_versions/GCF_000900555.1_ViralProj179422/GCF_000900555.1_ViralProj179422_genomic.gbff.gz`


```
    gene            209621..212353
                     /locus_tag="CcrKarma_gp333"
                     /db_xref="GeneID:13996474"
     CDS             209621..212353
                     /locus_tag="CcrKarma_gp333"
                     /note="TerL"
                     /codon_start=1
                     /transl_table=11
                     /product="intein-containing putative terminase large
                     subunit precursor"
                     /protein_id="YP_006989713.1"
                     /db_xref="GeneID:13996474"
                     /translation="MSYYPIEDRAKARSVV[...]LDLLIREIGQFPEGAHDDQVDAMTQYLRWAKSKRTRFGARKVGSMG"
     mat_peptide     join(209621..210010,211034..212350)
                     /locus_tag="CcrKarma_gp333"
                     /product="terminase large subunit"
                     /protein_id="YP_007001268.1"
     mat_peptide     210011..211033
                     /locus_tag="CcrKarma_gp333"
                     /product="intein"
                     /protein_id="YP_007001269.1"
```

`1211640	Caulobacter phage CcrColossus	Viruses;dsDNA viruses, no RNA stage;Caudovirales;Siphoviridae	11	/mirror/ncbi/current/genomes/refseq/viral/Caulobacter_phage_CcrColossus/latest_assembly_versions/GCF_000899635.1_ViralProj179419/GCF_000899635.1_ViralProj179419_genomic.gbff.gz`


Unrelevant sig peptide and mat peptide annoatation. Sig peptide annotation doesn't start at the begining. The mat peptide could be valid. Potential annotation mistake. It looks like they forget to add a sig peptide annotation in YP_294026.1
181082	YP_294025.1	-1		1	2	2	False	False	None		264	['Viruses', 'dsDNA viruses, no RNA stage', 'Phycodnaviridae', 'Coccolithovirus']
181082	YP_294026.1	1		1	1	1	False	False	None		303	['Viruses', 'dsDNA viruses, no RNA stage', 'Phycodnaviridae', 'Coccolithovirus']

```
gene            complement(241349..241612)
                /locus_tag="EhV269A"
                /db_xref="GeneID:3654832"
CDS             complement(241349..241612)
                /locus_tag="EhV269A"
                /codon_start=1
                /product="hypothetical protein"
                /protein_id="YP_294025.1"
                /db_xref="UniProtKB/TrEMBL:Q4A2L2"
                /db_xref="GeneID:3654832"
                /translation="MALFACLRIYSGFLRETDQENNRAIPVTPPHILINTIPPFDPGP
                PNAMDVNITTTRPNTTIIHPLALNILLFKYLTIKKYHKTVSVL"
sig_peptide     complement(241406..241483)
                /locus_tag="EhV269A"
                /inference="protein motif:SignalP:2.0"
                /note="Signal peptide predicted for EhV270 by SignalP 2.0
                HMM (Signal peptide probability 0.913, signal anchor
                probability 0.086) with cleavage site probability 0.757
                between residues 25 and 26"
gene            241407..241709
                /locus_tag="EhV270"
                /db_xref="GeneID:3654833"
CDS             241407..241709
                /locus_tag="EhV270"
                /codon_start=1
                /product="hypothetical protein"
                /protein_id="YP_294026.1"
                /db_xref="UniProtKB/TrEMBL:Q4A2L1"
                /db_xref="GeneID:3654833"
                /translation="MFSASGWIIVVLGLVVVMLTSIAFGGPGSKGGIVLISMWGGVTG
                IALLFSWSVSRRNPEYILRQANSAMGRQFTSVRRGTSRIGNSVRARLGSVRRTQVP"
mat_peptide     241482..241706
                /locus_tag="EhV270"
                /product="hypothetical protein"
                /protein_id="YP_002296226.1"

```



mat_peptide annotation covering all the cds expect the first codon and the stop codon. In this case there is no cleavage site annotation and consequently this kind of annotaion should not be considered as a potential polyprotein.

example in  `11886	Rous sarcoma virus	Viruses;Retro-transcribing viruses;Retroviridae;Orthoretrovirinae;Alpharetrovirus	1	/mirror/ncbi/current/genomes/refseq/viral/Rous_sarcoma_virus/latest_assembly_versions/GCF_000855425.1_ViralProj14978/GCF_000855425.1_ViralProj14978_genomic.gbff.gz` with the protein `NP_056888.1`

```
   CDS             7129..8709
                     /gene="src"
                     /locus_tag="RSVgp4"
                     /note="p60-SRC phosphoprotein"
                     /codon_start=1
                     /product="p60 src"
                     /protein_id="NP_056888.1"
                     /db_xref="GeneID:1491925"
                     /translation="MG[...]PACVLEVAE"
     mat_peptide     7132..8706
                     /gene="src"
                     /locus_tag="RSVgp4"
                     /product="pp60 SRC"
                     /protein_id="NP_955616.1"
```

Potential signal Peptide wrongly annotated in Halorubrum pleomorphic virus 1 2 3 and 6 (1156719 1156720 1156721 1156722).
The mat peptide cover almost all the length except the first 100 aa.

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3384331/
Related haloarchaeal pleomorphic viruses contain different genome types

> The next protein encoded in the conserved cluster of genes is VP4-like major structural protein (Figure 1), forming the spike structure on the virion surface. This protein is processed during maturation. In order to predict the translation start site, the empirically determined N-termini (Table 1) together with Signal P and Tat find programs were used. Three of the manually annotated VP4-like protein precursors were predicted to have a twin-arginine signal sequence suggesting export of the protein in a folded state. Most of the VP4-like protein precursors share only ∼20% identity (Supplementary Table S8). Even though VP4-like polypeptides are much more diverse than VP3-like proteins, the secondary structures are predicted to be similar. All of them are also predicted to contain a coiled-coil region in the C-terminus of the protein just preceding the transmembrane domain that serves as a membrane anchor.

```
gene            1132..2592
                   /gene="4"
                   /locus_tag="HGPV-1_gp04"
                   /db_xref="GeneID:11948291"
   CDS             1132..2592
                   /gene="4"
                   /locus_tag="HGPV-1_gp04"
                   /codon_start=1
                   /transl_table=11
                   /product="VP4 precursor"
                   /protein_id="YP_005454298.1"
                   /db_xref="GeneID:11948291"
                   /translation="MND[...]N"
   mat_peptide     1246..2589
                   /gene="4"
                   /locus_tag="HGPV-1_gp04"
                   /product="hypothetical protein"
                   /protein_id="YP_005454316.1"

```

In Clostridium phage phi3626 a major capsid protein NP_612835.1 has a mature peptide annotation starting at 114 aa of the CDS. The threshold for SP is 90 aa then this protein is not considered by the program has a SP annotation protein.

  `190478	Clostridium phage phi3626	Viruses;dsDNA viruses, no RNA stage;Caudovirales;Siphoviridae	11	/mirror/ncbi/current/genomes/refseq/viral/Clostridium_phage_phi3626/latest_assembly_versions/GCF_000839145.1_ViralProj14166/GCF_000839145.1_ViralProj14166_genomic.gbff.gz`

according to the paper attached to the gb file the part cleaved from the protein may have a function. Then this protein may be consider has a polyprotein. ?

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC135250/

Genomic Analysis of Clostridium perfringens Bacteriophage φ3626, Which Integrates into guaA and Possibly Affects Sporulation

> Posttranslational processing of the major head protein is frequently found in bacteriophages. During capsid maturation, φ3626 Cps is processed by removal of the first 114 residues. A similar processing can be observed in many other phages which lack a scaffold protein. As described for Sfi21, φPVL (11), and HK97 (9), φ3626 Cps revealed a possible coiled-coil structure in the amino-terminal part of the protein which is removed. It has been assumed by Duda and coworkers (14) that this domain might be a functional equivalent of a scaffold, fused to the capsid protein.

Same situation in
`51369	Lactobacillus phage A2	Viruses;dsDNA viruses, no RNA stage;Caudovirales;Siphoviridae	11	/mirror/ncbi/current/genomes/refseq/viral/Lactobacillus_phage_A2/latest_assembly_versions/GCF_000848025.1_ViralProj14602/GCF_000848025.1_ViralProj14602_genomic.gbff.gz`

```
gene            4378..5580
                /locus_tag="A2p05"
                /db_xref="GeneID:951713"
CDS             4378..5580
                /locus_tag="A2p05"
                /experiment="experimental evidence, no additional details
                recorded"
                /note="ORF5"
                /citation=[8]
                /codon_start=1
                /transl_table=11
                /product="major head protein"
                /protein_id="NP_680487.1"
                /db_xref="InterPro:IPR006444"
                /db_xref="UniProtKB/TrEMBL:Q8LTC0"
                /db_xref="GeneID:951713"
                /translation="MTL[...]KA"
mat_peptide     4747..5577
                /locus_tag="A2p05"
                /product="major head protein"
                /note="ORF5"
                /citation=[8]

```

same situation in `10719	Streptomyces virus phiC31	Viruses;dsDNA viruses, no RNA stage;Caudovirales;Siphoviridae;Phic31virus	11	/mirror/ncbi/current/genomes/refseq/viral/Streptomyces_virus_phiC31/latest_assembly_versions/GCF_000848045.2_ViralProj14606/GCF_000848045.2_ViralProj14606_genomic.gbff.gz`


The virus Enterobacteria phage phi92 has one sig peptide annotation covering the first codon of the CDS and 2 mature peptides.  
The protein has a "Intramolecular triple-beta-helix chaperone domain,highly conserved tailspike protein chaperone domain" in Cterminal that is cleaved off after it helps the protein to fold correctly.
Explanation of Intramolecular triple-beta-helix chaperone domain in this paper:

Schulz, E. C., Dickmanns, A., Urlaub, H., Schmitt, A., Mühlenhoff, M., Stummeyer, K., ... & Ficner, R. (2010). Crystal structure of an intramolecular chaperone mediating triple–β-helix folding. Nature structural & molecular biology, 17(2), 210.

This protein may be considered by the program as a polyprotein even if the mechanism is not compltely similar.  

`948870	Enterobacteria phage phi92	Viruses;dsDNA viruses, no RNA stage;Caudovirales;Myoviridae	11	/mirror/ncbi/current/genomes/refseq/viral/Enterobacteria_phage_phi92/latest_assembly_versions/GCF_000914915.1_ViralProj240593/GCF_000914915.1_ViralProj240593_genomic.gbff.gz`

```
gene            79764..82520
                /gene="PHI92_gene_143"
                /locus_tag="PHI92_gp248"
                /gene_synonym="endoN92_gene"
                /db_xref="GeneID:22277995"
CDS             79764..82520
                /gene="PHI92_gene_143"
                /locus_tag="PHI92_gp248"
                /gene_synonym="endoN92_gene"
                /EC_number="3.2.1.129"
                /function="precursor protein of endosialidase 92"
                /inference="protein motif:PFAM:PF12195"
                /inference="protein motif:PFAM:PF12217"
                /inference="protein motif:PFAM:PF12218"
                /inference="protein motif:PFAM:PF12219"
                /codon_start=1
                /transl_table=11
                /product="Phi92_gp143"
                /protein_id="YP_009012475.1"
                /db_xref="GOA:I7HXG2"
                /db_xref="InterPro:IPR001724"
                /db_xref="InterPro:IPR011040"
                /db_xref="InterPro:IPR023366"
                /db_xref="InterPro:IPR024427"
                /db_xref="InterPro:IPR024428"
                /db_xref="InterPro:IPR024429"
                /db_xref="InterPro:IPR024430"
                /db_xref="InterPro:IPR030392"
                /db_xref="PDB:4HIZ"
                /db_xref="UniProtKB/TrEMBL:I7HXG2"
                /db_xref="GeneID:22277995"
                /translation="MSTT[..]SK"
sig_peptide     79764..79766
                /gene="PHI92_gene_143"
                /locus_tag="PHI92_gp249"
mat_peptide     79767..82031
                /gene="endoN92*"
                /locus_tag="PHI92_gp247"
                /product="proteolytically matured endosialidase 92*
                (endoN92*), mature tailspike
                protein,endo-alpha-2,8-sialidase"
                /function="depolymerization of poly-alpha-2,8-sialic acid
                and poly-(alpha-2,8-alpha-2,9)-sialic acid (the capsular
                polymers of Escherichia coli K1 and K92
                strains,respectively)"
                /experiment="identified by mass spectrometry (Q-TOF
                MS/MS); apparent mol. mass from SDS-PAGE: 85 kDa"
                /protein_id="YP_009019178.1"
mat_peptide     82032..82517
                /gene="endoN92-CTD"
                /locus_tag="PHI92_gp246"
                /product="endoN92-CTD"
                /function="C-terminal intramolecular chaperone domain of
                endoN92"
                /note="Intramolecular triple-beta-helix chaperone
                domain,highly conserved tailspike protein chaperone
                domain. The chaperone domain is not found in the phage
                particle"
                /protein_id="YP_009019179.1"

```

Circular genome with a gene overlapping the position 1 of the genome.
The program don't manage circular DNA.. so it doesn't properly manage to map peptide over the gene.  
It has a regular CDS with sig pep + mat pep annotation.
`46014	NP_040808.1	True	True	2	0	1	False		636	mat_peptide	717	1	['Viruses', 'dsDNA viruses, no RNA stage', 'Plasmaviridae', 'Plasmavirus']`


2 CDS with mature peptides found in the genome of a Deltaretrovirus virus. This genomes is expected to have a polyprotein and in the gb another cds annotated is found with 3 mat peptide.
Problem with an Deltaretrovirus virus. Wrong identification.  Deltaretrovirus are expected to have one polyprotein. This genomes has 3 CDS annotated with mat peptide: one which is the polyprotein with 3 peptide and 2 other one with a peptide annotaion that cover only the main exon of the gene.

`33748	Simian T-lymphotropic virus 2	Viruses;Retro-transcribing viruses;Retroviridae;Orthoretrovirinae;Deltaretrovirus	1	/mirror/ncbi/current/genomes/refseq/viral/Primate_T-lymphotropic_virus_2/latest_assembly_versions/GCF_000860005.1_ViralProj15221/GCF_000860005.1_ViralProj15221_genomic.gbff.gz`


```
gene            5153..8207
                /locus_tag="STLV2gp07"
                /db_xref="GeneID:1724972"
CDS             join(5153..5155,7167..8207)
                /locus_tag="STLV2gp07"
                /codon_start=1
                /product="tax protein"
                /protein_id="NP_056910.1"
                /db_xref="UniProtKB/TrEMBL:O70643"
                /db_xref="GeneID:1724972"
                /translation="MAHF[..]SA"
mat_peptide     7167..8204
                /gene="tax"
                /locus_tag="STLV2gp11"
                /product="tax protein (short variant)"
                /protein_id="NP_954571.1"

```

One genome of Mammarenavirus has a CDS annottated with 2 mat pep. And the two pep cover all the sequence so ideal annotation. But in expected taxon file Arenaviridae's genomes are supposed to have 0 peptides.
`11619	Junin mammarenavirus	Viruses;ssRNA viruses;ssRNA negative-strand viruses;Arenaviridae;Mammarenavirus	1	/mirror/ncbi/current/genomes/refseq/viral/Junin_mammarenavirus/latest_assembly_versions/GCF_000856545.1_ViralMultiSegProj15028/GCF_000856545.1_ViralMultiSegProj15028_genomic.gbff.gz`

In taxons where we expect polyproteins


#### Identification of miss annotation

Example of CDS with not accurate annotaion:
