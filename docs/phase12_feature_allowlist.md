# Phase 12: Experiment 2 Biological Feature Allowlist

## Scope and enforcement

This is the positive allowlist for any future Experiment 2 feature matrix. It
does not authorize full annotation or model training. Annotation may retain a
separate provenance/status table, but only the predictor fields below may enter
the feature matrix. Every unlisted VEP/cache field is rejected by default.
Do not request or parse broad `--everything` output as a feature table.

The Ensembl release-116 cache contains ClinVar 2025-09. Therefore, filtering
must happen before an annotation record is joined to labels or supplied to a
model. The parser must assert that no prohibited field or target-derived value
is present and fail closed when VEP introduces an unrecognized feature field.

## Allowed predictors

| Feature group | Allowed values | Source / rule | Conditions before use |
|---|---|---|---|
| Consequence | Most severe consequence; unique transcript consequence ontology terms | VEP release 116 `most_severe_consequence`, transcript `consequence_terms` | Retain all terms as a set/long table; never derive from ClinVar classification |
| Gene context | Stable gene ID and gene symbol | VEP transcript consequence `gene_id`, `gene_symbol` | Treat as annotation context; no target-rate, pathogenicity-rate, or label-count encoding |
| Transcript context | Stable transcript ID, MANE Select/Plus Clinical status, exon/intron ordinal/count where returned | VEP transcript `transcript_id`, `mane_select`, `exon`, `intron` | Preserve transcript-level records; any scalar consequence/prediction policy must be prespecified |
| SIFT | Prediction and score | VEP `sift_prediction`, `sift_score` | Pin cache/component provenance; scalar only for exactly one MANE Select transcript; otherwise preserve missing/ambiguous status |
| PolyPhen-2 | Prediction and score | VEP `polyphen_prediction`, `polyphen_score` | Pin cache/component provenance; scalar only for exactly one MANE Select transcript; otherwise preserve missing/ambiguous status |
| Population frequency | Exact alternate-allele gnomAD genome/exome AF, kept separately | Release-matched cache `gnomadg` and `gnomade` namespaces | Require direct cache version/checksum proof, exact normalized chromosome-position-ref-alt match, and conflict audit; never substitute ClinVar colocated values |
| Reference context | GC fraction, reference/alternate sequence, or deterministic sequence encodings | Checksum-pinned GRCh38.p14 reference FASTA | Verify original REF; define window length, strand, indel handling, and ambiguous-base policy; only true reference-derived sequence |
| Variant representation | Normalized chromosome, 1-based position, REF, ALT; derived REF/ALT lengths and substitution/indel category | Existing GRCh38 input plus deterministic normalization | Keys may join annotations; exclude from the model unless approved as explicit non-target-derived predictors |

The release-116 cache documentation lists MANE v1.5, GENCODE 50, SIFT
6.2.1, PolyPhen-2 2.2.3, and gnomAD genomes/exomes v4.1. These are cache
documentation claims; the REST pilot did not identify all component versions
per record. Local cache fields are not approved until a one-variant and
10k-local test confirms their actual output and provenance.

## Prohibited predictors and output content

The following are prohibited, including aliases, normalized forms, derived
statistics, and values copied from colocated records:

- ClinVar `CLIN_SIG` / `clin_sig`, clinical significance, pathogenicity
  classification, assertion, submitter interpretation, and review status.
- ClinVar phenotype/disease assertions, ClinVar variation/allele IDs used as
  predictors, and any annotation field whose meaning reproduces or encodes a
  ClinVar target.
- `target`, `label`, `ClinicalSignificance`, label definitions, target counts,
  gene pathogenicity rates, target-conditioned encodings, or any statistic
  calculated using train/validation/test labels.
- ClinVar `Origin`, ClinVar curation/review metadata, and submission fields as
  predictors; their availability at prediction time is not established.
- Unrecognized VEP/cache fields, arbitrary colocated source fields, and
  unversioned or non-exact allele-frequency values.
- Fabricated, sampled, inferred, or label-conditioned values for CADD,
  conservation, gnomAD, SIFT, PolyPhen-2, or sequence context.

CADD remains deferred. Conservation is not approved until its exact GRCh38
track/build, checksum/ETag, coordinate semantics, and feasible indexed-query
method are validated.

## Required parser and join assertions

1. Store the original input row identifier and preserve all 1,702,326 source
   rows in the final left join. Unsupported alleles receive a status/reason;
   they are not silently dropped.
2. Build a normalized GRCh38 chromosome:position:REF:ALT key using the pinned
   reference; preserve the original representation and record normalization.
3. Keep provenance/status fields separate from predictors. Allow only the
   names and meanings above into the feature matrix.
4. Reject ClinVar and target-related field names case-insensitively and reject
   all unrecognized fields until reviewed.
5. Assert unique output keys or explicitly preserve multiple transcript rows
   in a separate long table; report duplicate and conflicting source matches.
6. Keep missing source annotations missing with reason codes. Never fill them
   from labels or replace absence with an unsupported biological value.
7. Preserve the established chromosome partitions: train 1–17, validation
   18–20, test 21, 22, X, Y, and MT/M as in the existing implementation.

Until the local VEP release/cache and source fields pass these checks, this
allowlist is a policy, not evidence that any full-dataset feature is ready.
