# Evidence summary presentation

The summary has Table (default) and Map views inside the existing Evidence summary section. This is presentation metadata, not a new evidence model. Levels and the existing `question`, `level`, `finding`, `evidence`, and `sources` strings retain their original values. JSON export keys are unchanged. The summary CSV adds only `on_request_steps`.

## Step state and baseline

`interp_request_steps()` supplies button and badge labels. `interp_summary_scope()` is the single question-to-step registry:

| Question | Requested steps |
|---|---|
| Protein sequence | public, spliceai |
| Protein abundance (pQTL) | public, spliceai, ld |
| Splicing | public, spliceai |
| Gene expression (eQTL) | public |
| Disease and traits | public, ld |
| Pathway | public |
| Tissue | public |

A non-null, non-error result is After; an error stored for that variant is failed; otherwise Before. Partial source failures inside an otherwise completed public response follow the existing completion rule. The global badges remain visible after every step completes, and show the public query timestamp when present. The UI omits the old pending-step narrative sentence; the backend reading and JSON reading remain unchanged.

The earlier level is displayed only from a captured all-Before summary after local annotation, splice-outcome computation and any applicable colocalisation have settled. Cache entries are keyed by variant/annotation context and checked against their local/RNA/colocalisation inputs. A fast click that prevents this snapshot, changed inputs, or an unavailable snapshot suppresses the earlier-level note. No pre-SpliceAI outcome is inferred by filtering a later table.

## Evidence metadata and layout

`interp_evidence_tag()` tags each evidence-producing expression with a source key, step, role and optional public gene identifiers. `interp_evidence_combine()` preserves these annotations in original order. `interp_evidence_summary(structured = TRUE)` returns the unchanged summary plus parallel evidence items; its default return remains the legacy data frame. No English sentence is parsed to infer its role or provenance.

The splice-outcome source is computed here; when its calculation includes structured SpliceAI events its presentation step is SpliceAI. This describes the inputs to that computation, not proof that SpliceAI alone causes a biological effect. Gene-level pathway/expression background and overflow notices are context; not-run lines are notices.

The shared preview count remains two. Evidence-role items are preferred; if a row has no evidence, its first two context items form the preview (for example Pathway). Notices always stay in details and never become links. A three-column row places the question/level/states in about 240 px, the finding and preview together in the remaining width, and sources in about 220 px. Findings and each preview wrap normally and clamp to two lines, with no character-width cap or nowrap. Full details is at the end of the body column; opening it unclamps the finding and previews and reveals the remaining items. State/source labels and baseline notes are at least 12 px. The narrow-screen layout stacks the same fields.

The obsolete 2,600-visible-character limit is removed. The browser contract fixes the summary container at 1,260 px within a 1728×895 viewport: a long preview must expose at least 150 characters before truncation; a shorter preview must be complete; at most two previews, each at most two lines; total collapsed height at most 1,050 px.

The inline SVG uses the same source registry and metadata. Solid edges represent contributed evidence/context; dashed gray edges represent scoped steps not completed. Private/cohort sources and app computations have no external links. Public-database local copies link to the same public entity when a verified template and identifier are available, and are marked local copy. Each graph node has a title and text is bounded within its viewBox. SVG height follows its width without a max-height squeeze. Generated numeric node IDs and pure CSS :has() selectors highlight a hovered/focused node's edges and fade the others; no JavaScript was added.

## Validated link templates

Only whitelisted scalar public identifiers can enter URLs. Arbitrary caller fields (including sample, genotype and cohort fields) are ignored. Rendering makes no requests; external pages open only on a click, using `target="_blank" rel="noopener noreferrer"`.

Checked with rs4680, 22_19963748_G_A and COMT / ENSG00000093010:

| Template | Test URL | Verification | Result |
|---|---|---|---|
| `https://platform.opentargets.org/variant/{variant}` | https://platform.opentargets.org/variant/22_19963748_G_A | HTTP + entity API; user confirmed browser | 200; enabled (SPA) |
| `https://platform.opentargets.org/target/{gene_id}` | https://platform.opentargets.org/target/ENSG00000093010 | HTTP + entity API | 200; enabled (SPA) |
| `https://gtexportal.org/home/snp/{rs_id}` | https://gtexportal.org/home/snp/rs4680 | HTTP + entity API; user confirmed browser | 200; enabled for public and local-copy records (SPA) |
| `https://gtexportal.org/home/gene/{gene}` | https://gtexportal.org/home/gene/COMT | HTTP + entity API | 200; enabled (SPA) |
| `https://www.ebi.ac.uk/gwas/variants/{rs_id}` | https://www.ebi.ac.uk/gwas/variants/rs4680 | HTTP + user confirmed rs4680, position and COMT | 200; local copy |
| `https://reactome.org/content/query?q={gene}` | https://reactome.org/content/query?q=COMT | HTTP + user confirmed Search results for COMT | 200; Reactome local-copy gene entries |
| `https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}` | https://www.ncbi.nlm.nih.gov/clinvar/variation/VCV000017591 | E-utilities esummary, db=clinvar, id=17591 | 200; identity verified; page may present a bot challenge |
| `https://cadd.gs.washington.edu/snv/GRCh38-v1.7/{cadd_variant}` | https://cadd.gs.washington.edu/snv/GRCh38-v1.7/22%3A19963748_G_A | Server-rendered variant, Direct link and matching CADD API | 200; local copy (SNVs only) |

GTEx local-copy links have a title explaining that the portal shows Release V10 and may differ from the local release (such as v11).

NCBI esummary returned accession VCV000017591, dbSNP 4680, gene COMT, GRCh38 position 19963748 and canonical SPDI `NC_000022.11:19963747:G:A`. The local ClinVar parser already reads the VCF ID column into `records$variation_id`; MyVariant supplies `clinvar_variant_id`. These identifiers now accompany the evidence metadata without changing the legacy summary strings. Both paths link to ClinVar. Missing or ambiguous IDs are not guessed. A VEP-cache record without an ID remains unlinked. MyVariant-only evidence (such as REVEL) is labelled via API rather than as an unlinked external page.

MSigDB has no independently verified gene-direct template in this change and remains unlinked local copy. Aggregate pathway/ClinVar gene context without a single public entity identifier also remains unlinked. TOPMed, FHS, cohort LD/RNA, colocalisation, splice outcome and pQTL fine-mapping retain local/computed markers and no href. Every registered key/template, the notice exclusion, unsafe identifier rejection and both ClinVar ingestion paths have tests. No network request is made while rendering.

## Regression evidence

The golden JSON contains only synthetic outputs captured at baseline e6d2ed6. The contract test runs the existing backend scenarios, compares legacy strings and checks metadata/text correspondence and step scope. UI tests cover per-variant cache restoration, failures, exact/missing baselines and the actual CSV download. Browser tests use seven synthetic rows (eight evidence items each, approximately 9,000 characters), measure in a fixed 1,260 px container at 1728×895, check each preview via character rectangles respecting two-line clipping, and retain all text in the DOM. The fixture includes a context-only Pathway row and all-before/all-after scenarios. Screenshots and metrics are written under `test-artifacts/evidence-summary-1260/`; no controlled data is used.
