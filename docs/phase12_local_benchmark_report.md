# Phase 12: Local VEP Benchmark Report

## Status: NOT RUN — BLOCKED BEFORE INSTALLATION

No VEP executable, release-116 cache, container runtime, or required
Bio::EnsEMBL Perl modules are installed. The release-116 cache was not
downloaded, and no local annotation was run. This is an intentional stop at
the pre-download storage gate, not a failed annotation result.

### Reason the local benchmark was not attempted

- Current free disk: 52 GiB.
- Release-116 GRCh38 cache compressed archive: exactly 27,644,657,162 bytes
  (about 25.74 GiB).
- Official extracted cache footprint: not published in the inspected
  documentation.
- A 20 GiB free-space operating reserve would leave only about 6.26 GiB beyond
  the compressed archive for extracted files, runtime/dependencies,
  benchmark outputs, and temporary data. Sufficiency cannot be established
  without the uncompressed footprint.
- Docker and Podman are unavailable. Native Perl 5.34.1 is among the versions
  tested by the VEP 116 README; required `Archive::Zip` and `DBI`, and the
  `gcc`/`g++`/`make` tools, are present, but
  `Bio::EnsEMBL::Registry`, `Bio::EnsEMBL::VEP`,
  `Bio::EnsEMBL::Variation`, and `Bio::EnsEMBL::Utils::Sequence` are missing.

## REST pilot baseline for a future same-key comparison

These are the measured Phase 11 REST pilot results, not local results:

| Field | REST pilot result |
|---|---:|
| Variants selected | 10,000 |
| VEP responses | 10,000 / 10,000 |
| GRCh38 reference-allele checks | 10,000 / 10,000 |
| Most severe consequence | 100.00% |
| Consequence terms | 99.99% |
| Gene symbols | 99.99% |
| MANE Select transcript IDs | 99.60% |
| SIFT | 10.59% |
| PolyPhen-2 | 10.12% |
| gnomAD genomes AF | 26.90% |
| gnomAD exomes AF | 36.69% |
| GC context | 100.00% |
| Runtime / peak RSS | 43.15 minutes / 588 MiB |

The REST payload confirmed Ensembl release 116 before/after and GRCh38.p14
assembly metadata, but did not identify all component data versions per
record. In particular, gnomAD and individual SIFT/PolyPhen-2 versions must be
confirmed directly from the local cache before those fields can be considered
version-pinned.

## Local-vs-REST comparison

| Comparable field | REST 10k coverage | Local 10k coverage / agreement | Conclusion |
|---|---:|---:|---|
| Consequence / most severe consequence | 99.99% / 100.00% | Not measured | No local comparison |
| Gene / transcript context | Gene 99.99%; MANE 99.60% | Not measured | No local comparison |
| SIFT | 10.59% | Not measured | No local comparison |
| PolyPhen-2 | 10.12% | Not measured | No local comparison |
| gnomAD genome/exome AF | 26.90% / 36.69% | Not measured | No local comparison |
| REF/allele validation | 10,000 / 10,000 | Not measured | No local comparison |
| Installation/cache size, startup, RSS, throughput | Not applicable | Not measured | No local benchmark |

No discrepancy has been assessed. REST values must not be treated as a local
cache validation or as proof of component-version identity.

## Go/no-go

**BLOCKED — STORAGE.** Reconsider local installation only when a suitable
volume can retain the verified archive/extraction workspace, outputs, and at
least 20 GiB of free operating headroom. After installation, run one
well-characterized GRCh38 variant, then the exact Phase 11 10k key set, before
designing full annotation. Do not annotate all 1,702,326 variants, train
Experiment 2, train the CNN, or download CADD as part of this phase.
