#!/usr/bin/env python3

# ACMG criteria from PMID: 25741868
# What data is output from the report
# GEMINI - 
# t1.chrom, t1.start, t1.end, t1.gene, t1.impact, t.impact_severity,
# t1.cadd_scaled,t1.polyphen_score,t1.sift_score, t1.clinvar_sig,t1.max_aaf_all,t1.sub_type,t1.num_het, t1.num_hom_alt,
# t1.exon, t1.codon_change, t1.aa_change, t1.aa_length, t1.rs_ids, t1.transcript, t1.biotype, t1.vep_hgvsc, t1.vep_hgvsp,
#              gts.VARPATIENT_ID, 
#              gt_depths.VARPATIENT_ID,
#              gt_quals.VARPATIENT_ID,
#              gt_ref_depths.VARPATIENT_ID,
#              gt_alt_depths.VARPATIENT_ID, 
#              t2.lof_gcount

# - need to add t1.clinvar_causal_allele,t1.clinvar_disease_name,t1.pfam_domain,t1.is_conserved

# HPO - VARPATIENT_ID_terms.csv
# - contains HPO and OMIM terms so need to process

# gene context, mechanism of disease and previous pathogenic variant
# Timo will alter phenoparser to take the final list of genes compiled from associations with all phenotypes and pull back all the variants in those genes linked to disease
# I am going to assume that we will get back something like:
# gene, rsID, disease, aa_change, clinvar_accession
import requests
import json

# gene list can be from individual patient report or probably better to implement after phenoparser has run so as to avoid redundancy
geneList = 'IRF6, NSUN2, BRIP1, BHMT2, EIF5B, PAM, KIAA1109, NHSL1, TTC40, IGSF3, GCN1L1, DGKQ, GABRP, PLEKHG1, TGS1, RHOBTB1, TRIM51, RBBP6, MYO19, CYP4F2, GPR116, PRDM10, TUBGCP3, KBTBD2, ABCD3, MTUS2, DIDO1, CTSL, PCDH7, PPP1R10, PET112, VPS52, APBB2, NEK11, OBSCN, PTH2R, ARPP21, WDR73, KCNG2, ANKRD6, KMT2E, SERGEF, MYBPC2, ARPP21, LSG1, MTDH, AP2A1, CDC20B, ZNF555, ANKRD52, JMJD4, COL3A1, FMNL1, KLHDC1, FSCN2, UNC13A, COBL, ZNF786, PCF11, SIGLEC15, SELO, BMP5, FREM1, VPS52, TAS2R41, OR5AR1, TIPARP, AMDHD1, GPR180, ECE2, SLC6A15, C15orf61, TECPR1, FAM20B, FBXO24, GUCA1B, UBR4, TRMT10B, DSC3, FBLN1, PRKCD, SCARA5, SLC2A12, SLC6A8, HOXA4, PTPDC1, RELN, ZNF565, CEACAM19, CLASP1, LMTK3, ZSWIM4, COPS7A, ITCH, SYCE1, SNX1, SGSM1, POGK, WDR38, NAPRT1, ANKS4B, OR6Y1, PPP1R1B, ASXL3, EIF2S3, SLC8A3, C15orf39, RUSC2, PIWIL2, PI3, SLC38A6, TARS, CCDC142, PTK2, ITPR3, VASP, PPFIA2, MT-CYB, NDUFAF2, INTS2, ARC, ZNF697, HYDIN, MT-ND4, MT-ND5, HYDIN, MT-CO3, MT-ND4, HYDIN, KIR2DL4, HYDIN'
#phenotype_list = 'HP:0000202, HP:0002744, HP:0100267, HP:0100336'

for gene in geneList.split(','):
    print( gene )
    # recreate the following URL but replacing with the relevant gene
    # "https://api.omim.org/api/entry/search?search=IRF6&include=allelicVariantList&format=json&apiKey=Wqy5lssmS7uWGdpyy8H9zw"
    r_omim = requests.post('https://api.omim.org/api/entry/search?', verify=True, data={'search': gene, 'include': 'allelicVariantList', 'format':'json', 'apiKey':'Wqy5lssmS7uWGdpyy8H9zw'})
    # error if we get anything other than a 200 status code
    if ( r_omim.status_code == 200 ):
        # return type is JSON so process the text to an associative array
        parsed_json_omim = json.loads(r_omim.text)
        # allelicVariants is contained in one of the "entry" records in "entryList"
        # need to search each entrylist record (values) for allelicVariants
        for entry in parsed_json_omim['omim']['searchResponse']['entryList']:
            for value in entry.values():
                if 'allelicVariantList' in value:
                    # variants and consequences are contained within "allelicVariantList" list
                    for avar in value['allelicVariantList']:
                        # find the interesting fields we want
                        for item in ['dbSnps','name','mutations','status','clinvarAccessions','number','mimNumber']:
                            # check the field is in the allelicVariant record first
                            if item in avar['allelicVariant']:
                                itemval = avar['allelicVariant'][item]
                                print( "{} => {}".format(item,itemval) )
                                # if we are looking at the dbSNP field then head off the ensembl for the SNP data
                                if item == 'dbSnps':
                                    # this is grch37
                                    url = 'https://grch37.rest.ensembl.org/variation/human/'+itemval
                                    # again return type is JSON so bring the data back to an associative array
                                    r_ensembl = requests.get(url, data={'content-type':'application/json'}, verify=True)
                                    if ( r_ensembl.status_code == 200 ):
                                        parsed_json_ensembl = json.loads(r_ensembl.text)
                                        # interesting fields
                                        for param in ['mappings','clinical_significance','synonyms']:
                                            # again check that the field exists in the record
                                            if param in parsed_json_ensembl:
                                                # if it's in there then process the contents
                                                for data in parsed_json_ensembl[param]:
                                                    # the mappings record contains the location information we want
                                                    if param == 'mappings':
                                                        for coord in ['seq_region_name','strand','start','end','allele_string']:
                                                            print ( "\t{} => {}".format(coord,data[coord]) )
                                                    # the clinical_significance and synonyms records are arrays
                                                    else:
                                                        print ( "\t{} => {}".format(param,data) )
                                            else:
                                                print("\t{} not in this record".format(param))
                                    else:
                                        print("\tfailed with error code {}".format(r_ensembl.status_code))
                                        continue

                            else:
                                print("{} not in this record".format(item))

                else: 
                    continue
    else:
        print( "failed with error code {}".format( r_omim.status_code ) )
        continue
    break

# Amelie
#r = requests.post('https://amelie.stanford.edu/api/', verify=False, data={'genes': gene_list, 'phenotypes': phenotype_list})

# Very strong
# PVS1 null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon deletion) in a gene where LOF is a known mechanism of disease
#  Caveats:
#    •  Beware of genes where LOF is not a known disease mechanism (e.g., GFAP, MYH7)
#    •  Use caution interpreting LOF variants at the extreme 3′ end of a gene
#    •  Use caution with splice variants that are predicted to lead to exon skipping but leave the remainder of the  protein intact
#    •  Use caution in the presence of multiple transcripts
###
# We already track information pertaining to the consequence via the GEMINI framework and report this in HGVS nomenclature
#	Also report position of variant in protein and which exon, from most severely affected transcript
# We report consequence in relation to the most affected transcript not necessarily the most biologically relevant transcript
#	Maybe use expression studies or preferably just use longest transcript or CCDS transcript (GEMINI gene_detailed table)
# We don’t currently establish how intact the resultant protein would be
#	GEMINI gene_detailed table has CDS and protein lengths
# Phenoparser gives us genes linked to phenotype from OMIM and Phenolyser



# Strong
# PS1 Same amino acid change as a previously established pathogenic variant regardless of nucleotide change
#  Example: Val→Leu caused by either G>C or G>T in the same codon
#  Caveat: Beware of changes that impact splicing rather than at the amino acid/protein level
#PS2 De novo (both maternity and paternity confirmed) in a patient with the disease and no family history
#  Note: Confirmation of paternity only is insufficient. Egg donation, surrogate motherhood, errors in embryo transfer, and so on, can contribute to nonmaternity.
#PS3 Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene  product
#  Note: Functional studies that have been validated and shown to be reproducible and robust in a clinical diagnostic laboratory setting are considered the most well established.
#PS4 The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls
#  Note 1: Relative risk or OR, as obtained from case–control studies, is >5.0, and the confidence interval around the estimate of relative risk or OR does not include 1.0. See the article for detailed guidance.
#  Note 2: In instances of very rare variants where case–control studies may not reach statistical significance, the prior observation of the variant in multiple unrelated patients with the same phenotype, and its absence in controls, may be used as moderate level of evidence.
###
# We don’t currently compare our AA change with a known disease causing AA change
# 	We’d need to obtain this from OMIM (has allelicVariants field) or Phenolyser maybe via clinVar entry currently within report
#	Need to confirm transcript context (maybe via clinVar HGVS)
# We do not have parental samples



# Moderate
# PM1 Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of an enzyme) without benign variation
# PM2 Absent from controls (or at extremely low frequency if recessive) (table 6) in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium
#  Caveat: Population data for insertions/deletions may be poorly called by next-generation sequencing.
# PM3 For recessive disorders, detected in trans with a pathogenic variant
# Note: This requires testing of parents (or offspring) to determine phase.
# PM4 Protein length changes as a result of in-frame deletions/insertions in a nonrepeat region or stop-loss variants
# PM5 Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before
#  Example: Arg156His is pathogenic; now you observe Arg156Cys
#  Caveat: Beware of changes that impact splicing rather than at the amino acid/protein level.
# PM6 Assumed de novo, but without confirmation of paternity and maternity
###
# We do not currently have any functional data upon which to base decisions
# GEMINI reports allele frequency in population databases
# 	PathWest seemed to favour gNOMAD (genome aggregation database - http://gnomad.broadinstitute.org/)
#	Currently report max AF in all databases
#		Not necessarily race-matched
#	GEMINI uses old pre-compiled data
# We don’t report any protein length deviation
# 	Could in theory but may need additional work to confirm annotation in a non-repeat region
# As with PS1 we need to obtain allelicVariant/ClinVar data
# We can’t determine de novo without parental data


# Supporting
# PP1 Cosegregation with disease in multiple affected family members in a gene definitively known to cause the disease
#  Note: May be used as stronger evidence with increasing segregation data
# PP2 Missense variant in a gene that has a low rate of benign missense variation and in which missense variants
# are a common mechanism of disease
# PP3 Multiple lines of computational evidence support a deleterious effect on the gene or gene product  (conservation, evolutionary, splicing impact, etc.)
#   Caveat: Because many in silico algorithms use the same or very similar input for their predictions, each algorithm should not be counted as an independent criterion. PP3 can be used only once in any evaluation of  a variant.
# PP4 Patient’s phenotype or family history is highly specific for a disease with a single genetic etiology
# PP5 Reputable source recently reports variant as pathogenic, but the evidence is not available to the laboratory to perform an independent evaluation
###
# We do not have family data
# Could use Ensembl to query rates of benign missense variation in genes
#	Not sure where to reference for typical mechanism of a given disease
# Current pipeline plus inclusion of Denise’s work could cover PP3


# Stand-alone
#  BA1 Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium
# Strong
#  BS1 Allele frequency is greater than expected for disorder (see table 6)
#  BS2 Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked
#(hemizygous) disorder, with full penetrance expected at an early age
#  BS3 Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing   BS4 Lack of segregation in affected members of a family
#    Caveat: The presence of phenocopies for common phenotypes (i.e., cancer, epilepsy) can mimic lack of segregation among affected individuals. Also, families may have more than one pathogenic variant contributing to an autosomal dominant disorder, further confounding an apparent lack of segregation.
###
# Mode of inheritance currently only obtained from OMIM if terms are available
#	Allele frequency is reported but not automatically assessed in context of mode of inheritance
#	Population databases assumed normal but resolution makes BS2 hard
# No information available relating to functional studies
# No current family data

#Supporting
#  BP1 Missense variant in a gene for which primarily truncating variants are known to cause disease
#  BP2 Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed in cis with a pathogenic variant in any inheritance pattern
#  BP3 In-frame deletions/insertions in a repetitive region without a known function
#  BP4 Multiple lines of computational evidence suggest no impact on gene or gene product (conservation, evolutionary, splicing impact, etc.)
#    Caveat: Because many in silico algorithms use the same or very similar input for their predictions, each algorithm cannot be counted as an independent criterion. BP4 can be used only once in any evaluation of a variant.
#  BP5 Variant found in a case with an alternate molecular basis for disease
#  BP6 Reputable source recently reports variant as benign, but the evidence is not available to the laboratory to perform an independent evaluation
#  BP7 A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site AND the nucleotide is not highly conserved
###
# Variant type reported but not in context of typical consequence associated with disease
# No ability to phase to investigate cis/trans
# Could use Ensembl API to identify variant context (in repeat regions etc.)
# Current pipeline plus inclusion of Denise’s work could cover PP3
# Need to know disease first before knowing alternate molecular basis
# PathWest may report this but this is part of interpretation rather than automation
# BP7 would be ranked low by CADD

### Combining criteria
# Pathogenic
#  (i) 1 Very strong (PVS1) AND
#    (a) ≥1 Strong (PS1–PS4) OR
#    (b) ≥2 Moderate (PM1–PM6) OR
#    (c) 1 Moderate (PM1–PM6) and 1 supporting (PP1–PP5) OR
#    (d) ≥2 Supporting (PP1–PP5)
#  (ii) ≥2 Strong (PS1–PS4) OR
#  (iii) 1 Strong (PS1–PS4) AND
#      (a)≥3 Moderate (PM1–PM6) OR
#    (b)2 Moderate (PM1–PM6) AND ≥2 Supporting (PP1–PP5) OR
#    (c)1 Moderate (PM1–PM6) AND ≥4 supporting (PP1–PP5)
#Likely pathogenic
#  (i) 1 Very strong (PVS1) AND 1 moderate (PM1– PM6) OR
#  (ii) 1 Strong (PS1–PS4) AND 1–2 moderate (PM1–PM6) OR
#  (iii) 1 Strong (PS1–PS4) AND ≥2 supporting (PP1–PP5) OR
#  (iv)  ≥3 Moderate (PM1–PM6) OR
#  (v) 2 Moderate (PM1–PM6) AND ≥2 supporting
#(PP1–PP5) OR
#  (vi) 1 Moderate (PM1–PM6) AND ≥4 supporting
#(PP1–PP5)
# Benign
#  (i) 1 Stand-alone (BA1) OR
#  (ii) ≥2 Strong (BS1–BS4)
# Likely benign
#  (i) 1 Strong (BS1–BS4) and 1 supporting (BP1– BP7) OR
#  (ii) ≥2 Supporting (BP1–BP7)
# Uncertain  significance
#  (i) Other criteria shown above are not met OR
#  (ii) the criteria for benign and pathogenic are contradictory
