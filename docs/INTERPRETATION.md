# Variant interpretation

The **Variant interpretation** tab answers, for the variant that is open (or, before a comparison exists, the record selected in step 1): what does it do to transcripts and proteins, is it tied to splicing, expression, disease or a pathway, and in which tissues. It reads; it never reruns BAM counting, refits Model 1 or changes a saved result.

It has three layers, in the order they appear.

| Layer | When it runs | Network |
|---|---|---|
| Transcript and protein consequences, local sources | Automatically, a few seconds | None |
| Public databases (Open Targets, GTEx, MyVariant.info) | **Query public databases** button | Sends the variant and gene IDs only |
| SpliceAI | **Run SpliceAI** button, about 1–2 minutes | None |
| LD with GWAS Catalog variants | **LD with GWAS variants** button, about 30 seconds | None |

A public query sends the variant's chromosome, position, REF and ALT, and the Ensembl IDs of the genes it overlaps. No sample, genotype count, cohort name or result is sent. Answers are kept for the browser session and written only to the files you download.

A button pressed while the local lookup is still loading waits for it and then starts; the line under the buttons says so. If another variant is selected before that, the request is dropped and nothing is sent.

## Transcript and protein consequences

Computed in the app from the GENCODE release selected in the top bar and the dataset's reference FASTA, for every transcript within 5 kb.

- Terms are Sequence Ontology terms in Ensembl VEP's severity order; `impact` is VEP's HIGH / MODERATE / LOW / MODIFIER grouping.
- Intronic positions: `splice_donor_variant` at +1/+2, `splice_acceptor_variant` at −1/−2, `splice_donor_5th_base_variant` at +5, `splice_donor_region_variant` at +3, +4 and +6, `splice_region_variant` at +7/+8 and −3 to −8, `splice_polypyrimidine_tract_variant` at −3 to −17. Exonic bases within 3 of an intron-adjacent exon edge add `splice_region_variant`.
- Coding SNVs get the codon, amino acids, residue number and HGVS-style `c.`/`p.` labels. The stop codon counts as coding; `c.-N` and `c.*N` are 5′ and 3′ UTR positions; intronic labels are `c.N+D` or `c.N-D` from the nearer exon. chrM uses the vertebrate mitochondrial code and has no upstream/downstream rows.
- `donor_distance` and `acceptor_distance` count into the intron from its 5′ and 3′ end (1 = first or last intronic base). Negative values are exonic bases counted from the exon edge.
- `cds_start_NF` transcripts are numbered from their annotated partial CDS and are flagged. A REF allele that disagrees with the reference sequence is flagged at the top of the tab.
- Insertions, deletions and multi-base changes are classified by what the changed bases overlap (frameshift, in-frame, splice site, UTR, intron). No protein change is computed for them.

Check: 432 transcript-level calls on 266 SNVs in APOE, HBB, KRAS, PCSK9, BRCA1 and CNIH4 were compared with the VEP consequences published by Open Targets; all term sets and all 82 amino-acid changes agreed (GENCODE v48 here, Ensembl's own release there, compared on shared transcript IDs).

## Protein consequence of a splice change

A variant that changes splicing changes the protein through the transcript it produces. For every splice change tied to the variant, the tab re-splices the coding transcripts the change fits, from the reference sequence, translates them from their annotated start codon and compares the product with the annotated protein. It runs by itself after the local sources are in, and again when SpliceAI or an open result adds a change.

Where the changes come from:

- an intron that is used more with the ALT allele: a significant TOPMed sQTL intron with a positive slope, or a junction of the open result whose Model 1 slope is positive at the table's Bonferroni cut-off (*observed*);
- a splice site SpliceAI predicts to appear or to weaken (*predicted* at 0.5 or more, *weakly predicted* from 0.2). A predicted new site that is an end of an intron seen to rise with ALT supports that intron, because the RNA shows what it splices to; otherwise it is paired with the annotated site it would splice to (the next annotated acceptor downstream of a new donor, the nearest annotated donor upstream of a new acceptor);
- an essential splice site (intron positions 1, 2 and the last two) that the variant itself changes (*predicted*).

How a change is placed on a transcript: the new intron must share at least one end with an intron of that transcript and must not be one of its introns already. It is spliced out in place of the annotated introns it overlaps, so exonic bases inside it are removed and intronic bases outside it are kept. Sharing both ends is exon skipping; sharing one is an alternative donor or acceptor, which may also leave exons out. A change that would keep more than 300 nt of intron is not placed on that transcript. When a site is lost and no change shows what replaces it, the two usual outcomes are listed as *assumed outcome*: the neighbouring exon skipped, and the intron retained. When the variant's ALT allele ends up in the new mRNA it is read as ALT. When the annotated intron that the new one stands in for is itself used less with ALT (a negative TOPMed sQTL slope or Model 1 slope), `support` says so: the two together show a switch from one intron to the other, not only an addition.

What the table says, per transcript and change:

- `mrna_change`: exonic bases removed and intronic bases kept.
- `coding_effect` and `protein_change`: unchanged coding sequence; in-frame deletion, insertion or substitution (the annotated stop codon is still read in frame) with the residues concerned; premature stop or frameshift with the last unchanged residue and the number of new residues; no stop codon; start codon removed.
- `nmd`: nonsense-mediated decay is *predicted* when the new stop codon ends 50 nt or more before the last exon junction of the new mRNA, the usual rule; otherwise a shortened or changed protein may be made.
- `annotated_in`: transcripts of the selected GENCODE release that already contain the intron, with their type. An intron that ALT favours and that belongs to a transcript without a coding sequence points the same way as a predicted premature stop.

Check: for 300 MANE Select transcripts drawn at random from ten chromosomes (162 on the plus strand, 138 on the minus strand), re-splicing the annotated exons and translating gave the annotated coding sequence in all 300, with the protein length equal to the CDS length divided by three, less the stop codon. Skipping each of 832 internal coding exons in turn gave an in-frame change for all 329 whose length is a multiple of three (296 plain deletions of exactly a third of the exon's length in residues, 31 where the codon across the new junction also changes, 2 where that codon is a stop) and a frameshift, premature stop or missing stop for all 503 others.

Limits. The prediction says what the changed splice would encode, not how much of the gene's RNA is spliced that way: a strong sQTL can still move a small share of transcripts. Translation always starts at the annotated start codon; an alternative start, a second splice change on the same molecule and a new exon made of two new sites are not modelled. The 50 nt rule is a rule of thumb, and decay is rarely complete. Transcripts longer than 2.5 Mb are skipped.

## Local sources

Every path can be overridden in `config/runtime.Renviron`; a source that is missing or unreadable is listed as unavailable in the **Sources** sub-tab and skipped.

| Source | Variable | What it gives | Lookup |
|---|---|---|---|
| CADD v1.7 (GRCh38) | `REGSHINY_CADD_SNV`, `REGSHINY_CADD_INDEL` | Raw and PHRED score of the exact allele | tabix |
| Ensembl VEP 112 cache | `REGSHINY_VEP_CACHE` | rsID, ClinVar significance (2023-10), 1000 Genomes and gnomAD frequencies | tabix |
| TOPMed multi-tissue cis-QTL | `REGSHINY_TOPMED_QTL` | cis-eQTL and cis-sQTL in lung, monocytes, nasal epithelium, PBMC, T cells and whole blood | tabix |
| Cross-context sharing model | `REGSHINY_TISSUE_MODEL_EQTL`, `REGSHINY_TISSUE_MODEL_SQTL` | Effect and its certainty in each TOPMed context, from a model trained on those files | read |
| FHS whole-blood isoform-ratio sQTL | `REGSHINY_FHS_SQTL` (`discovery;replication`) | Rows with p < 1e-4, N about 2,622 and 1,093 | `look` |
| FHS whole-blood gene eQTL | `REGSHINY_FHS_EQTL` | Rows with p < 1e-7, N about 2,622 | `look` |
| GTEx v11 cis-eQTL | `REGSHINY_GTEX_EQTL`, `REGSHINY_GTEX_EGENES` | Significant variant-gene pairs in 50 tissues, and the genes tested per tissue | tabix, `look` |
| cis-pQTL summary statistics | `REGSHINY_PQTL_UKBPPP`, `REGSHINY_PQTL_FHS`, `REGSHINY_PQTL_PROTEIN_MAP` | Plasma level of a protein against the variants near its gene: UKB-PPP (2,910 proteins, about 34,000 participants) and the FHS scan (2,543 proteins); assay-to-gene table | tabix |
| GWAS Catalog associations (e114, 2025-07-21) | `REGSHINY_GWAS_CATALOG` | Published trait associations reported at the variant's position | tabix |
| ClinVar VCF (2025-07-21) | `REGSHINY_CLINVAR_VCF` | Classification, review status and conditions of the exact allele | tabix |
| ClinVar conditions per gene | `REGSHINY_GENE_DISEASE` | Pathogenic and likely pathogenic ClinVar variants of each linked gene and their conditions | `look` |
| Reactome gene-to-pathway mapping | `REGSHINY_REACTOME` | Pathways of the genes linked to the variant | `look` |
| MSigDB gene sets (2026.1) | `REGSHINY_MSIGDB` (directory) | Hallmark, canonical pathways and GO biological process of those genes | read |
| SpliceAI 1.3.1 | `REGSHINY_SPLICEAI` | Predicted splice-site gain and loss | SCC module wrapper |

`REGSHINY_TABIX` names the tabix binary; without it the app uses the one on `PATH` or the newest `htslib` module on SCC. The local sources are GRCh38 only.

**TOPMed cis-QTL.** For each of the six contexts the app reads every tested pair of the variant from the nominal all-pairs file. A pair is `significant` when its nominal p is at or below that gene's permutation threshold (tensorQTL, genes with q ≤ 0.05). `pip` is the SuSiE posterior inclusion probability when the variant sits in a credible set. `slope_alt` is per ALT allele on normalised phenotypes, so it compares across contexts but not with Model 1. `n_estimated` is derived from the minor-allele count and frequency. Whole blood also has a run down to MAF 0.1%, used when the 1% run has no row. The table keeps significant pairs, credible-set pairs, pairs with p < 0.001, the three smallest p per context, and every context's row of a phenotype that is significant somewhere.

**Shared or tissue-specific.** For every gene (eQTL) or intron (sQTL, matched by intron coordinates and gene because cluster numbers differ between runs) that is significant in at least one context, the tab reports Cochran's Q and I² of the slopes across contexts and a one-line pattern: shared, consistent effect size but not significant everywhere, context-specific, heterogeneous, or opposite directions. The contexts share donors and differ several-fold in sample size, so Q is descriptive and a small context can miss a real effect.

**Cross-context model.** The comparison above looks at each context's own significance. The model re-estimates the effect in every context from all of them together, in the manner of multivariate adaptive shrinkage (Urbut et al. 2019, *Nat Genet*). It is trained on the TOPMed files themselves by `scripts/train_tissue_model.R`, once for eQTL and once for sQTL:

- *Strong pairs*: the lead variant of about 6,000 genes or introns with q ≤ 0.05, drawn evenly from the six contexts. *Random pairs*: about 5,000 variants 20–200 kb from a lead variant, chosen without looking at their p-values. Each pair's slope and standard error is read in all six contexts.
- The correlation of estimates between contexts when there is no effect (the contexts share donors) comes from the random pairs. Sharing patterns are fixed ones (the same effect everywhere, one context only, partly correlated) plus patterns learned from the strong pairs. Their weights are fitted by EM on both sets, with a penalty that favours "no effect" when the data do not decide.
- One fifth of the phenotypes is held out. For those pairs each context is predicted from the others and compared with what was measured, against two baselines (the precision-weighted mean slope and the mean z of the other contexts). The model is fitted on the slope scale and on the z scale and the one with the smaller held-out error is kept. The numbers are stored in the model file and shown above the table.

At view time, for every gene or intron significant in at least one context, the tab shows per context the re-estimated slope (`posterior_slope`), its uncertainty and the local false sign rate `lfsr`: the probability that the effect is absent or has the other sign. `lfsr` < 0.05 is called an effect. `leading_pattern` is the kind of sharing pattern that carries most of the posterior weight and `log10_bayes_factor` compares "some effect" with "none". Two consequences: a small context can be called without passing its own threshold, because it borrows strength from the others; and a context where the pair was not tested gets a prediction, labelled as such.

Limits: the model knows six contexts from blood, airway and lung and says nothing about other tissues. It describes the association signal at this variant, which a neighbour in linkage disequilibrium can carry. Lead variants are selected for significance, so intervals on such pairs are somewhat too narrow (see the stored coverage). Retrain when the TOPMed files change:

```bash
Rscript scripts/train_tissue_model.R --type eqtl --topmed DIR --out tissue_sharing_eqtl.json --work WORKDIR
Rscript scripts/train_tissue_model.R --type sqtl --topmed DIR --out tissue_sharing_sqtl.json --work WORKDIR
```

**Colocalisation.** A variant is often associated with several molecular traits at once: the splicing of one gene, the expression of others, a protein level. Whether those are one signal or several decides how to read them. When the local sources give at least two such traits (per gene the strongest significant sQTL intron and the eQTL, each in the TOPMed context where it is strongest, up to four; and up to two proteins from the cis-pQTL scans), every pair is compared over the variants within 1 Mb of this one that both were tested on, using their full summary statistics and the method of Giambartolomei et al. (2014, *PLoS Genet*; the calculation of `coloc.abf`). It runs by itself after the local sources are in.

- `shared_cause` is the posterior probability that one variant causes both signals, `distinct_causes` that each has its own causal variant; `only_trait_1`, `only_trait_2` and `neither` complete the five hypotheses. A pair is read as shared or distinct when that probability is at least 0.8.
- `this_variant_share` is this variant's part of the shared signal and `leading_variant` the variant with the largest part.
- For each trait on its own, the status line gives this variant's posterior probability within that trait's signal (the same single-signal calculation as for the pQTL scans), which says whether the variant leads that trait.

How to read it: a shared cause ties two traits together, for example an eQTL and the level of the same gene's protein. Distinct causes between a trait this variant leads and an eQTL it does not lead say that the eQTL is probably reached through linkage disequilibrium with another causal variant; the Gene expression line then says so. Limits: at most one causal variant per trait is assumed in the region, no LD reference is used, and the priors are 1e-4 for each trait and 1e-5 for a shared variant. Several independent signals in one trait, or cohorts with different LD (TOPMed against UKB-PPP), push the result towards distinct causes, so a shared cause is the more robust of the two readings. Two introns of one gene are not compared, because intron ratios within a gene move against each other by construction. At least 50 variants must be common to both traits.

**FHS whole-blood QTL.** Isoform-ratio sQTL (transcript / gene) and gene eQTL summary statistics, effects per ALT allele. A transcript "replicates" when it appears in both the discovery and the replication set with the same sign. The gene eQTL file must be a position-sorted copy; make it once with

```bash
scripts/build_qtl_index.sh SOURCE.csv DEST.csv
```

which only reads the source.

**GTEx cis-eQTL across tissues.** The significant variant-gene pairs of all 50 GTEx v11 tissues, read from one position-indexed copy. For each gene the variant is an eQTL of, the tab gives the number of tissues where the pair is significant, the number of tissues in which GTEx tested the gene at all (`tissues_tested`, from the per-tissue gene tables), the direction of the ALT allele, the strongest tissue, the status in whole blood (significant; tested, not significant; not tested), and the number of tissues where this variant is the gene's most significant variant (`lead_variant_in`). The plot shows the slope and its 95% interval in every significant tissue. Only significant pairs are stored, so there is no slope for a tissue where the pair was tested and missed the threshold; GTEx tissues differ several-fold in donors, so a count of tissues is descriptive and a small tissue can miss an effect. Being the lead variant is a hint, not fine-mapping. The copies are made once from the release's `*.eQTLs.signif_pairs.parquet` and `*.eGenes.txt.gz` files (needs `pyarrow`):

```bash
python3 scripts/build_gtex_index.py --source GTEx_Analysis_v11_eQTL --pairs gtex_v11_eqtl_signif_pairs.b38.tsv.gz --genes gtex_v11_egenes.tsv
```

Name the pairs file `gtex_vNN_...` so the release shows in the tab.

**cis-pQTL scans.** Two scans of plasma protein levels against the variants near each protein's gene, kept with every p-value (not only the significant rows), so a variant that was tested and is not associated is told apart from one that was not tested: UKB-PPP (Olink Explore 3072, European discovery set, about 34,000 participants; Sun et al. 2023, *Nature*) and the Framingham Heart Study scan on Olink Explore HT (a project file; the model it was fitted with is not recorded with the file, so read it as supporting evidence). The table lists this variant against every nearby protein: effect per ALT allele (`beta_alt`, turned to ALT when the scan's effect allele is REF), standard error, p-value, sample size, the distance to the protein's gene and a `significance` label: experiment-wide (p ≤ 1.7e-11, UKB-PPP's published threshold of 5e-8 divided by 2,923 proteins), genome-wide (p ≤ 5e-8) or suggestive (p ≤ 1e-5).

For the proteins the variant is associated with at p ≤ 5e-8 (the four strongest per scan) the protein's whole cis window is read and fine-mapped under the assumption of one causal variant: each variant gets Wakefield's approximate Bayes factor from its effect and standard error (prior standard deviation 0.15 on the scan's standardised scale, the calculation `coloc` uses for one trait), and the factors are normalised over the window. `posterior_probability` is this variant's share, `credible_set_size` the number of variants needed to reach 95%, and `lead_variant` the variant with the largest share when it is not this one. A posterior near 1 says that, if the window has one causal variant, the data point to this one; a posterior near 0 next to a small p-value says the variant is a neighbour in linkage disequilibrium with a stronger one. The assumption fails where a window holds several independent signals, and no LD reference is used, so treat a posterior between about 0.1 and 0.9 as undecided. A variant that changes the protein's sequence can also change how the assay binds it, which shows as a cis-pQTL without a change in abundance.

The copies are made once with

```bash
scripts/build_pqtl_index.sh regenie UKB_CIS_DIR ukbppp_eur_cis_pqtl.b38.tsv.gz      # one <protein>.csv.gz per protein
scripts/build_pqtl_index.sh table FHS_CIS_PQTL.txt.gz fhs_olink_cis_pqtl.b38.tsv.gz  # one tab-separated table
python3 scripts/build_protein_map.py --ukbppp olink_protein_map_3k_v1.tsv --olink-ht olink-exht-annot.csv --out olink_protein_map.b38.tsv
```

The assay table also says where each protein is encoded. For a GWAS Catalog trait that is the blood level of one protein, the tab uses it to label the association cis (the protein's gene lies within 1 Mb of the variant) or trans (`protein_locus` in the catalog and LD tables); a protein the table does not hold stays undetermined.

**GWAS Catalog.** Every association the catalog reports at the variant's position, strongest first: trait, mapped ontology trait, p-value, effect, risk allele, study and PubMed ID. The catalog keys an association by rsID and position, not by allele; `risk_allele_is` says whether the reported risk allele is this variant's ALT, its REF, another allele or not reported. A p-value too small for a double is stored as 0 and shown as `< 1e-300`. The position-indexed copy is made once from the catalog's "all associations v1.0.2" file with

```bash
scripts/build_gwas_index.sh gwas_catalog_v1.0.2-associations_eNNN_rYYYY-MM-DD.tsv DEST.tsv.gz
```

Rows without a mapped position are left out and a row that reports several variants becomes one row each. Name the copy `..._eNNN_rYYYY-MM-DD...` so the release shows in the tab.

**LD with nearby GWAS Catalog variants.** A variant that is not in the catalog can still carry a published signal if it travels with the variant a study reported. On request the tab takes every variant the catalog reports within 500 kb (the nearest 300 at a dense locus, which is stated), reads their genotypes and this variant's from the connected dataset's own VCF with `bcftools`, and computes r² between ALT dosages over the samples that have both calls. Rows with r² ≥ 0.2 are listed with the trait, p-value and study, and `risk_allele_travels_with` says which allele of this variant the reported risk allele is found with. Only correlations and counts are kept; no genotype and no sample identifier is shown or exported. Where a position has several VCF records, the one whose ID is the catalog's rsID is used; a multi-allelic record gives one row per ALT allele. This is the LD of this cohort: a study in another ancestry can differ, and r² near 1 says the two variants cannot be told apart by association, not which one is causal.

**ClinVar.** The record of the exact REF/ALT allele from a tabix-indexed ClinVar VCF: germline classification, review status, conditions, gene and molecular consequence. Records of other alleles at the same position are counted and not shown. When the VCF is unreadable the older significance stored in the VEP cache (2023-10) is used; after a public query the per-condition records from MyVariant.info are shown instead.

**Disease genes among the linked genes.** For each linked gene (below), the number of ClinVar variants classified pathogenic or likely pathogenic, conflicting classifications left out, and the conditions they are reported for, most variants first. It answers "what happens when this gene is damaged", which is the disease context of a variant that acts through the gene. It is gene-level context: it is listed under Disease and traits as such and never changes the level. The table is made once from the ClinVar VCF with

```bash
python3 scripts/build_gene_disease_index.py --clinvar clinvar.vcf.gz --out clinvar_gene_disease_YYYY-MM-DD.tsv
```

**Linked genes and pathways.** The genes a variant can act through are the ones it overlaps, the ones whose splicing or expression is significantly associated with it in a TOPMed context or a GTEx tissue (`sQTL target`, `eQTL target`), and the ones whose plasma protein level is associated with it at p ≤ 5e-8 in a cis-pQTL scan (`cis-pQTL protein`). For each, the tab lists its Reactome pathways (all levels; `top_level` marks Reactome's top-level pathways) and its MSigDB gene sets (Hallmark; canonical pathways from KEGG, WikiPathways, PID and BioCarta; GO biological process; `genes_in_set` is the set's size, smaller being more specific). Reactome is matched by Ensembl gene ID, MSigDB by gene symbol. The Reactome copy is made once with

```bash
scripts/build_pathway_index.sh Ensembl2Reactome_All_Levels.txt DEST.tsv
```

Pathways describe the gene. They do not show that the variant changes the pathway.

**SpliceAI.** Runs the SCC module on a one-line VCF with the dataset's reference and a 500 bp window. Scores are the change in predicted splice-site probability (0 to 1); common cut-offs are 0.2 (permissive), 0.5 (recommended) and 0.8 (high precision). SpliceAI reports the exonic boundary base; `site_position` converts it to the intronic base next to the site, the coordinate used in Genotype association, and `annotated_site` says whether the selected GENCODE release has a site there. When the open result has a bin-ratio row for that site, it is shown beside the prediction. SpliceAI and its models are free for academic and not-for-profit use only.

## Public sources

| Source | What it gives |
|---|---|
| Open Targets Platform (variant) | GWAS credible sets containing the variant, with its posterior probability and the top locus-to-gene prediction; eQTL, sQTL, pQTL and single-cell eQTL credible sets; ClinVar and other curated evidence; in-silico predictors; pharmacogenomics |
| Open Targets Platform (gene) | Function, UniProt ID, Reactome pathways, Gene Ontology terms, loss-of-function constraint, most associated diseases |
| GTEx v10 | Significant single-tissue eQTL and sQTL for the variant; median expression of the overlapped genes in each tissue |
| MyVariant.info | ClinVar records with condition and review status, rsID, gnomAD frequency, REVEL |

Each source fails on its own: the others are still shown and the failure is named. Up to 2,000 GWAS and 2,000 molecular-QTL credible sets are read per variant; the count is shown when more exist. `tau` in the Tissue line is the usual tissue-specificity index of the gene's GTEx medians (0 = even, 1 = one tissue).

## Evidence summary

Seven fixed questions, each with a level, a one-sentence finding, the evidence lines behind it and the sources. A card shows its first four evidence lines; the others open on request. Above them a short reading restates the result: the best-supported molecular effect with its level, what else is supported, the disease or trait link, and which on-request steps (SpliceAI, the LD scan, the public query) have not been run for this variant. The reading adds no evidence of its own.

| Level | Meaning |
|---|---|
| Strong | Direct or fine-mapped evidence that this variant is responsible |
| Moderate | An association or a prediction, but not both |
| Limited | Weak or indirect support only |
| None found | The sources were consulted and report nothing |
| Not assessed | The sources needed were unavailable or not queried |
| Context | Describes the gene or the tissue, not evidence about the variant |

Rules:

- **Protein sequence.** Strong: a HIGH-impact term (stop gained, frameshift, stop lost, start lost, essential splice site) on a MANE Select or Ensembl canonical transcript. Moderate: the same on other transcripts only, or a missense or in-frame change. Limited: a missense change that AlphaMissense calls likely benign with CADD below 20. None found: no change to the protein sequence. The predicted consequence of a splice change (above) can raise this level, never lower it, when it changes the protein: to Strong when the change is both observed and predicted and falls on a MANE Select or Ensembl canonical transcript; to Moderate when it is observed or predicted but not both, or on other transcripts only; to Limited when it is weakly predicted or an assumed outcome. The finding then says what the variant itself does and what the splice change would encode.
- **Protein abundance (pQTL).** From the local cis-pQTL scans: Strong when the variant is associated with a protein at p ≤ 5e-8 and single-signal fine-mapping gives it a posterior probability ≥ 0.5, or when it is associated at p ≤ 1.7e-11 with the protein of a gene whose product it is itself predicted to change (a coding or essential splice-site consequence in that gene, or SpliceAI ≥ 0.5 for it). Moderate: associated at p ≤ 5e-8 without either. Limited: 5e-8 < p ≤ 1e-5 only. None found: the variant was tested against the nearby proteins and none reaches 1e-5. A variant that is in neither scan (no measured protein is encoded within the scans' window, or the variant was not tested) is Not assessed. From the GWAS Catalog: Moderate for an association at p ≤ 5e-8 with the blood level of a protein (mapped trait "level of … in blood") reported at this variant or at a variant in LD with it at r² ≥ 0.8 in this cohort; the finding says how many of those proteins are encoded near the variant (cis) or elsewhere (trans). From the public query: Strong for a cis pQTL credible set with posterior probability ≥ 0.9, Moderate for any credible set with posterior probability ≥ 0.1, Limited below that. The highest of these is shown.
- **Splicing.** *Associated* means any of: a Model 1 row of the open result (bin ratio or junction) with permutation p below 0.05 divided by the number of fitted rows of its table, sites flagged `NO_DEPTH_DROP_AT_SITE` set aside; a significant TOPMed sQTL; an FHS isoform-ratio sQTL that replicates in cis; a GTEx sQTL. *Pinned to the variant* means any of: an essential splice-site term; SpliceAI ≥ 0.5; a TOPMed sQTL SuSiE probability ≥ 0.5; an Open Targets splicing or transcript-usage credible set with probability ≥ 0.5. Strong: associated and pinned. Moderate: associated, or predicted (essential site or SpliceAI ≥ 0.5), not both. Limited: a splice-region term or SpliceAI 0.2–0.5 only.
- **Colocalisation** adds evidence lines under Protein abundance, Splicing and Gene expression for the pairs it decides, and does not change a level. One sentence is added to the Gene expression finding when an eQTL that this variant does not lead (posterior below 0.1 within that eQTL) has a different causal variant from a trait this variant does lead (posterior of at least 0.5).
- **Gene expression (eQTL).** Strong: a TOPMed or Open Targets eQTL credible set with posterior probability ≥ 0.5. Moderate: a significant eQTL in TOPMed, FHS, GTEx (the local copy or the public query) or Open Targets without that. GTEx's significant-pair files carry no fine-mapping, so a GTEx eQTL alone stops at Moderate.
- **Disease and traits.** Strong: ClinVar pathogenic or likely pathogenic with no benign assertion, or a GWAS credible set with posterior probability ≥ 0.5. Moderate: a GWAS Catalog association at p ≤ 5e-8 reported at the variant's position, a GWAS posterior probability of 0.1–0.5, ClinVar risk factor, drug response, association or conflicting assertions, or a pharmacogenomic annotation. Limited: other ClinVar records, GWAS Catalog associations above 5e-8 only, or GWAS credible sets with posterior probability below 0.1. After an LD scan, r² ≥ 0.8 with another reported variant at p ≤ 5e-8 also counts as Moderate, and r² of 0.2–0.8 as Limited. A catalog association shows that the variant or a neighbour in linkage disequilibrium is associated, so on its own it stops at Moderate; fine-mapping needs the public query, and the card says when it was not run. None found needs both the local GWAS Catalog and a ClinVar source to have been read; otherwise an empty answer reads Not assessed. A last evidence line names the linked genes in which ClinVar has pathogenic variants, with the leading condition; it is marked as gene-level context and does not enter the level.
- **Pathway** is a Context line: for each linked gene, its Reactome top-level pathways, or its most specific gene sets when Reactome has none. It works without the public query; after one, Open Targets' counts for the overlapped genes are added.
- **Tissue** is a Context line: where the QTL effect is seen (each TOPMed context's own significance, then the cross-context model's reading after "model:"; then, per gene, the GTEx tissues where the variant is an eQTL out of those that tested the gene, the direction, and the status in whole blood) and where the gene is expressed. This cohort's RNA is whole blood.

What the levels are not: a clinical classification, a causal proof, or a replacement for reading the tables. A QTL association can be carried by a neighbour in linkage disequilibrium, which is why association alone stops at Moderate. Fine-mapping probabilities depend on the reference panel and the study. Model 1 is unadjusted for covariates and relatedness.

## Exports

**Interpretation JSON** holds everything shown (consequences, the predicted splice outcomes with the splice changes they were computed from, the colocalisation of the variant's molecular traits, every local and public table, SpliceAI, the LD table, the cohort RNA evidence used, the summary). **Evidence summary CSV** and **Consequences CSV** hold those two tables. Public and local reference content keeps its own terms of use; cite the sources.
