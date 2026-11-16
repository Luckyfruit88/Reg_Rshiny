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

The shared visible-evidence constant is two. Only evidence-role items enter that preview. Findings are clamped to two lines and long evidence previews are ellipsized. One Full details disclosure reveals the full finding, unclamps the preview evidence and exposes every remaining line, including notices and context. No evidence text is deleted. The narrow-screen layout stacks the same fields.

The inline SVG uses the same source registry and metadata. Solid edges represent contributed evidence/context; dashed gray edges represent scoped steps not completed. Local sources and app computations have no external links. Each graph node has a title and text is bounded within its viewBox.

## Validated link templates

Only whitelisted scalar public identifiers can enter URLs. Arbitrary caller fields (including sample, genotype and cohort fields) are ignored. Rendering makes no requests; external pages open only on a click, using `target="_blank" rel="noopener noreferrer"`.

Checked 2026-11-16 with rs4680, 22_19963748_G_A and COMT / ENSG00000093010:

| Template | Test URL | HTTP | SPA | Enabled |
|---|---|---:|---|---|
| `https://platform.opentargets.org/variant/{variant}` | https://platform.opentargets.org/variant/22_19963748_G_A | 200 | yes | yes |
| `https://platform.opentargets.org/target/{gene_id}` | https://platform.opentargets.org/target/ENSG00000093010 | 200 | yes | yes |
| `https://gtexportal.org/home/snp/{rs_id}` | https://gtexportal.org/home/snp/rs4680 | 200 | yes | yes |
| `https://gtexportal.org/home/gene/{gene}` | https://gtexportal.org/home/gene/COMT | 200 | yes | yes |

A SPA's HTTP 200 alone is not entity verification. The Open Targets GraphQL variant and target queries also returned the exact variant, rs4680, ENSG00000093010 and COMT. GTEx `dataset/variant?snpId=rs4680&datasetId=gtex_v10` and `reference/gene?geneId=COMT&gencodeVersion=v39` independently returned those entities. The SPA destinations still merit interactive browser review.

The tested MyVariant hg38 endpoint returned 404 and is not a link template. ClinVar search returned 200 but did not verify a specific variation page, so no ClinVar template was added. Ensembl redirected to a SPA; it is not used for locally computed GENCODE evidence. Local TOPMed, FHS, GWAS Catalog copies, ClinVar copies, CADD, cohort LD/RNA and computed results remain unlinked.

## Regression evidence

The golden JSON contains only synthetic outputs captured at baseline e6d2ed6. The contract test runs the existing backend scenarios, compares legacy strings and checks metadata/text correspondence and step scope. UI tests cover per-variant cache restoration, failures, exact/missing baselines and the actual CSV download. Browser tests use seven synthetic rows (eight evidence items each, approximately 9,000 characters), measure at 1728×895, check visible text via character rectangles respecting clipping/closed details, and retain all text in the DOM. Screenshots and metrics are written under `test-artifacts/evidence-summary/`; no controlled data is used.
