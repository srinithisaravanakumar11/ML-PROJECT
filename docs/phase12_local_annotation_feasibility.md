# Phase 12: Local Ensembl VEP 116 Feasibility

## Decision: B) BLOCKED — STORAGE

Do not download the release-116 cache or start full annotation. No VEP/cache
installation or local benchmark was attempted. The archive alone is 27,644,657,162
bytes, while the Mac currently reports 52 GiB free. Ensembl's release
documentation does not publish an extracted-cache size, so available space
cannot be shown to cover extraction, installation, output, and a safe reserve.

The decision is conservative and specific: current storage is inadequate to
prove a safe cache install. Missing software/modules are additional setup work,
but not the primary gate. This does not claim that VEP is intrinsically
incompatible with Apple Silicon.

## Environment audit

| Item | Observed value | Evidence / implication |
|---|---|---|
| Operating system | macOS 27.0.0, Darwin 27.0.0 | `uname -a` |
| Architecture | arm64 (Apple Silicon) | Kernel reports `RELEASE_ARM64_T8132` |
| RAM | 17,179,869,184 bytes = 16 GiB | `sysctl -n hw.memsize` |
| Current free disk | 52 GiB available on the data volume | `df -h .`; rounded display |
| Python | 3.9.6 | Configured project interpreter / environment snapshot |
| Perl | 5.34.1, darwin-thread-multi | `perl -v`; VEP 116 README recommends Perl >=5.22 and lists 5.34 among tested versions |
| Docker | Not installed / not on PATH | `command -v docker` |
| Podman | Not installed / not on PATH | `command -v podman` |
| VEP executable | Not installed / not on PATH | `command -v vep` |
| Bio::EnsEMBL::Registry | Missing | Direct Perl module load failed |
| Bio::EnsEMBL::VEP | Missing | Direct Perl module load failed |
| Bio::EnsEMBL::Variation | Missing | Direct Perl module load failed |
| Bio::EnsEMBL::Utils::Sequence | Missing | Direct Perl module load failed |
| Archive::Zip / DBI | Installed | VEP 116 README lists these as required Perl libraries |
| gcc / g++ / make | Installed | VEP 116 README lists these build tools |
| DBD::mysql | Missing | Listed as optional for offline cache use; required for database access |
| Set::IntervalTree | Missing | Listed as optional; recommended for speed improvements |

No Perl dependencies, containers, or scientific resources were installed during
this audit. The Phase 11 pilot's 588 MiB peak RSS is a Python REST client's
measurement and is not a local VEP RAM benchmark.

## Release-116 cache facts

| Requirement | Verified information |
|---|---|
| VEP release | Ensembl VEP 116.0; use the release-116 source/tag and matching cache |
| Human cache archive | `homo_sapiens_vep_116_GRCh38.tar.gz` |
| Official URL | `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz` |
| Exact archive size | 27,644,657,162 bytes (`Content-Length` from official FTP HEAD; 27.64 decimal GB, about 25.74 GiB) |
| Official extracted size | Not stated in the inspected Ensembl cache documentation or archive HTTP metadata |
| Assembly | GRCh38 cache; release documentation/Phase 11 endpoint metadata identifies GRCh38.p14 |
| Cache data includes ClinVar | Yes: release-116 cache documentation lists ClinVar 2025-09 |
| gnomAD | Cache documentation lists genomes v4.1 and exomes v4.1; only accessioned/available variants are represented, not every gnomAD record |
| SIFT | 6.2.1 |
| PolyPhen-2 | 2.2.3 |
| MANE | v1.5 |
| GENCODE | 50 |
| Archive timestamp metadata | Last modified 2026-04-09; ETag `66fbffe0a-64efc55395049` |

The official documentation recommends a cache matching the VEP release. The
documented manual install is a cache download followed by extraction; the VEP
installer can fetch cache/reference files too. Neither method was run here.

Official references:

- [Ensembl release-116 cache documentation](https://jun2026.archive.ensembl.org/info/docs/tools/vep/script/vep_cache.html)
- [Ensembl release-116 download and installation documentation](https://jun2026.archive.ensembl.org/info/docs/tools/vep/script/vep_download.html)
- [Ensembl VEP 116 source README](https://github.com/Ensembl/ensembl-vep/blob/release/116/README.md)
- [Release-116 GRCh38 cache archive](https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz)

## Storage calculation and stop condition

The cache download alone would use about 25.74 GiB of the currently available
52 GiB, leaving approximately 26.26 GiB while the compressed archive is still
present. Using a conservative 20 GiB free-space reserve for the OS, project,
temporary files, and safe recovery leaves only about 6.26 GiB for the
uncompressed cache increment, VEP software/dependencies, local benchmark
inputs/outputs, and working files. The uncompressed cache size is not
officially documented, so that remainder cannot be shown to be sufficient.

Required free space during a safe download/extraction is:

`compressed archive + extracted cache + VEP/runtime files + benchmark/output/temp files + operating reserve`

The archive and extracted cache may coexist during checksum verification and
extraction. Removing the archive before validating the extracted data is not
an acceptable shortcut. Current disk capacity therefore fails the safe
pre-download gate. Reconsider only after the exact extracted footprint is
measured on a separate suitable volume or more capacity is available, with at
least 20 GiB left free after installation and benchmark workspace are staged.

No archive bytes were downloaded. No local cache size, install size, startup
time, RAM use, or one-variant VEP result is claimed. The inspected official
VEP documentation does not provide a minimum/expected RAM figure for this
offline cache workload; the 16 GiB machine must therefore be benchmarked after
the storage gate is resolved.

## Production pipeline design (not launched)

Once a suitable versioned local cache and reference are available:

1. Freeze the source dataset SHA-256, schema, row identifiers, and the current
   chromosome-based split. Exclude target and label columns from annotation
   inputs.
2. Preserve all 1,702,326 rows. Assign each row a stable source-row ID and
   status; do not silently drop unsupported or invalid alleles.
3. Normalize chromosome, 1-based coordinate, REF, and ALT against a
   checksum-pinned GRCh38.p14 FASTA. Retain original alleles and record
   normalization decisions. Use the canonical normalized key for joins.
4. Split deterministically into bounded VCF/JSON batches (initial production
   benchmark will determine batch size). Record each batch's input row range,
   input checksum, tool/config/cache/reference versions, start/end time, exit
   status, and output checksum.
5. Write each completed batch to a temporary file, validate schema and row/key
   counts, then atomically rename it into a completed-batch directory. Resume
   only when the batch checksum and all version/config identifiers match.
6. Emit progress and failure counts to structured logs. Keep unsupported,
   failed, no-match, and reference-mismatch rows with explicit status/reason
   values; retry only transient local failures.
7. Parse VEP output with [the Phase 12 feature allowlist](./phase12_feature_allowlist.md).
   Keep transcript-level records separate; avoid arbitrary transcript
   aggregation; exclude ClinVar and target-related fields.
8. Merge completed batches by stable row ID and normalized key. Assert every
   source row appears exactly once; report duplicate keys, multi-transcript
   rows, conflicting annotation matches, missingness, and REF mismatches.
9. Produce a manifest with source/software/cache/reference SHA-256 values,
   exact versions and command options, environment details, per-batch logs,
   run metrics, and final output checksums.
10. Train nothing until the local 10k benchmark, REST-vs-local discrepancy
    review, full output validation, and full-cohort missingness report pass.

## Required local benchmark before full annotation

After clearing the storage gate, first install only the matching VEP 116
runtime and GRCh38 cache on a volume with measured headroom. Run one known
GRCh38 variant as a smoke test, then annotate exactly the same 10,000 sampled
keys in `results/metrics/phase11_10k_variants.csv`. Compare local and Phase 11 REST
results for consequence, gene/transcript, MANE, SIFT, PolyPhen-2, exact
population AF, and reference/allele validation. Preserve both outputs, report
per-field agreement/missingness/conflicts, and investigate each discrepancy.
Do not compare fields whose release/provenance differs as though they were
identical. Measure install/cache disk, startup, run time, peak RSS, throughput,
and output size. Only then estimate full-cohort local runtime and storage.

The 122.41-hour / 42,559-request figure is the Phase 11 **REST** extrapolation;
it is not a local VEP runtime estimate.

## Conservation-source note

UCSC's official GRCh38 100-way conservation tracks are indexed BigWig files,
but are not small on this machine:

| Track | Assembly / scope | Official file size | Last modified |
|---|---|---:|---|
| `hg38.100way.phastCons100way.bw` | GRCh38/hg38, 100-way | 5,886,377,734 bytes (about 5.48 GiB) | 2015-05-08 |
| `hg38.100way.phyloP100way.bw` | GRCh38/hg38, 100-way | 9,870,053,206 bytes (about 9.19 GiB) | 2015-05-08 |

These could be evaluated one at a time on separate storage, but neither was
downloaded or queried for scores. Both official servers accepted a one-byte
HTTP Range request (206 response), so an indexed remote-range approach is a
candidate for a later storage-light pilot. That probe establishes HTTP range
support only; it does not validate BigWig random-access tooling, API/service
terms, release immutability, coordinate semantics, throughput, or variant
coverage. The full files remain 5.48 GiB and 9.19 GiB, respectively; no
smaller versioned conservation data artifact was identified in this audit.
Any future point-query pilot must pin checksums/build, demonstrate data
coverage, and benchmark request volume before proposing full use. CADD v1.7
remains deferred and its approximately 81 GB file was not downloaded.
