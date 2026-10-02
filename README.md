# *Aspergillus nidulans* genome: from raw reads to biosynthetic gene clusters

I put this pipeline together to learn fungal genome assembly and annotation end to end, before applying it to my own marine fungal isolates. I chose *Aspergillus nidulans* FGSC A4 on purpose: it has a well-annotated reference genome, so I could check every step against a known answer. Everything ran on a laptop (Apple Silicon MacBook, 16 GB RAM) with an external drive for the large files.

## Results at a glance

| Step | What I measured | Result | Compared with |
| --- | --- | --- | --- |
| Trimming | Read pairs kept | 12,254,337 (about 80x coverage) | 13,520,670 raw pairs |
| Assembly | Size, scaffolds of 500 bp or more | 29.52 Mb in 332 scaffolds | Reference: 29.83 Mb |
| Assembly | N50 / L50 | 274.7 kb / 33 | |
| Assembly | Reference genome covered | 97.7% | QUAST |
| Completeness | BUSCO complete (eurotiales_odb12) | 97.3% | Reference genome: 97.0% |
| Repeats | Genome soft-masked | 2.28% | |
| Genes | Gene models predicted | 9,350 | Reference: 10,453 protein-coding genes |
| Genes | Predictions matching a reference protein | 93.6% (47.5% identical) | DIAMOND |
| Gene clusters | BGC regions found | 66 | 30 of them at a scaffold edge |

## Data

- **Reads:** SRA run [SRR4236261](https://www.ebi.ac.uk/ena/browser/view/SRR4236261) (BioProject PRJNA335082), Illumina HiSeq 2500, 2 x 101 bp, 13.5 million pairs (2.73 Gb, roughly 91x for a 30 Mb genome). This is JGI resequencing of an FGSC A4 line, so a few real differences from the reference are possible.
- **Reference:** RefSeq assembly GCF_000011425.1 (ASM1142v1), with its GFF3 annotation and proteins. I only used it to evaluate my results. The assembler never saw it.

## Setup

Most Bioconda packages don't have Apple Silicon builds, so I installed Intel (osx-64) versions, which run through Rosetta. Each tool got its own conda environment to avoid dependency clashes.

| Environment | Tools |
| --- | --- |
| `ncbi` | NCBI datasets |
| `genome` | FastQC, MultiQC, fastp, BUSCO |
| `spades` | SPAdes |
| `quast` | QUAST (plus setuptools, see the troubleshooting table) |
| `repeats` | RepeatModeler, RepeatMasker |
| `augustus` | AUGUSTUS, gffread, gffcompare, DIAMOND, Biopython |

Partway through I switched from conda to mamba, after a BUSCO environment spent three hours "solving". Since then I install like this, always with a dry run first:

```bash
mamba create -p ~/opt/anaconda3/envs/repeats --platform osx-64 \
  --override-channels -c conda-forge -c bioconda \
  "repeatmodeler>=2" repeatmasker --dry-run
```

`--override-channels` keeps Anaconda's default channels out of the environment, which avoids licensing issues and mixed-channel conflicts.

## Repository layout

```
nidulans_assembly/
├── 00_raw_reads/       raw FASTQ (not in the repo)
├── 01_qc_raw/          FastQC/MultiQC on raw reads
├── 02_trimmed/         trimmed reads + fastp report
├── 03_qc_trimmed/      FastQC/MultiQC on trimmed reads
├── 04_assembly/        filtered assembly (SPAdes output on external drive)
├── 05_evaluation/      QUAST and BUSCO
├── 06_repeats/         RepeatModeler and RepeatMasker (external drive)
├── 07_annotation/      masked genome, AUGUSTUS genes, proteins, DIAMOND results
├── 08_antismash/       fungiSMASH input, results and region summary
├── reference/          NCBI reference package
└── scripts/
    ├── filter_scaffolds.py
    ├── filter_augustus_gff.py
    └── summarise_antismash.py
```

Several folders are symbolic links to the external drive, because SPAdes, BUSCO and RepeatModeler produce far more temporary data than my internal disk could hold.

## Pipeline

### 1. Download

FTP downloads kept dropping partway, so I switched to HTTPS and let curl resume in a loop until each file was complete:

```bash
until curl -C - -O https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR423/001/SRR4236261/SRR4236261_1.fastq.gz; do sleep 5; done
until curl -C - -O https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR423/001/SRR4236261/SRR4236261_2.fastq.gz; do sleep 5; done
```

I checked the file sizes against ENA (829,554,401 and 831,409,525 bytes) before going on.

The reference came from NCBI:

```bash
datasets download genome accession GCF_000011425.1 --include genome,gff3,protein --filename nidulans_ref.zip
```

### 2. Quality control and trimming

```bash
fastqc -t 2 -o 01_qc_raw 00_raw_reads/*.fastq.gz
multiqc 01_qc_raw -o 01_qc_raw
```

The raw reads were good overall (101 bp, 48% GC, nothing flagged as poor quality), with two problems:

- Every overrepresented sequence was TruSeq adapter. 4.6% of R1 reads *started* with adapter, which means adapter dimers.
- R1 failed per-base quality. After a first round of trimming it still failed, and the per-position table showed why: only position 1 was bad (median Q18; everything from position 6 on was Q38). So I added `--trim_front1 1`.

Final trimming command:

```bash
fastp \
  -i 00_raw_reads/SRR4236261_1.fastq.gz \
  -I 00_raw_reads/SRR4236261_2.fastq.gz \
  -o 02_trimmed/SRR4236261_1.trimmed.fastq.gz \
  -O 02_trimmed/SRR4236261_2.trimmed.fastq.gz \
  --trim_front1 1 \
  --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
  --length_required 50 \
  --thread 4 \
  --html 02_trimmed/fastp_report.html \
  --json 02_trimmed/fastp_report.json
```

I kept quality trimming moderate on purpose, since SPAdes corrects errors itself and over-trimming just throws away coverage.

| fastp result | Value |
| --- | --- |
| Reads kept | 24.51 M of 27.04 M (90.6%) |
| Adapter dimers removed | 1.63 M reads (6.0%) |
| Too short after trimming | 2.0% |
| Low quality / too many N | 0.95% / 0.3% |
| Q30 bases, before → after | 91.5% → 94.7% |
| GC, before → after | 48.3% → 49.2% |

After trimming, FastQC passed everything except GC content and length distribution. The length warning is expected (reads now range from 50 to 100 bp), and the GC warning is common in fungal genomes because of mitochondrial DNA and repeats.

### 3. Assembly

```bash
caffeinate -i spades.py --isolate \
  -1 02_trimmed/SRR4236261_1.trimmed.fastq.gz \
  -2 02_trimmed/SRR4236261_2.trimmed.fastq.gz \
  -t 6 -m 12 -o 04_assembly/spades
```

`--isolate` is the SPAdes mode for high-coverage data from a single cultured organism. It picked k = 21, 33 and 55 by itself. The run took about 35 minutes once it had enough disk space (the first attempt died with "No space left on device"). `caffeinate -i` just keeps the Mac awake. Note that macOS can't enforce the `-m` memory limit, so I closed everything else while it ran.

SPAdes gave 2,518 scaffolds, most of them tiny. I wrote `scripts/filter_scaffolds.py` to drop anything under 500 bp and report the before/after statistics:

```bash
python3 scripts/filter_scaffolds.py \
  04_assembly/spades/scaffolds.fasta \
  04_assembly/scaffolds_min500.fasta \
  --min-length 500
```

| | Sequences | Total length | N50 |
| --- | --- | --- | --- |
| Before | 2,518 | 29,782,932 bp | 274,685 bp |
| After (≥ 500 bp) | 332 | 29,515,703 bp | 274,685 bp |

The 2,186 removed fragments add up to only 267 kb. My script's numbers matched QUAST's exactly, which was a nice independent check on both. `scaffolds_min500.fasta` is the assembly I used from here on.

### 4. Assembly evaluation

**QUAST** against the reference:

```bash
REF=reference/nidulans_ref/ncbi_dataset/data/GCF_000011425.1
quast.py 04_assembly/spades/scaffolds.fasta -r $REF/*.fna -g $REF/genomic.gff -t 6 -o 05_evaluation/quast
```

| Metric | Value |
| --- | --- |
| Total length / reference | 29,515,703 / 29,828,291 bp |
| GC / reference | 50.11% / 50.37% |
| Largest scaffold | 1,199,232 bp |
| N50 / NGA50 | 274,685 / 270,620 bp |
| Genome fraction | 97.70% |
| Duplication ratio | 1.000 |
| Misassemblies | 20 (in 15 scaffolds) |
| Mismatches / indels per 100 kb | 8.49 / 3.16 |
| Unaligned | 366 kb (20 whole + 42 partial scaffolds) |

NGA50 barely drops from N50, so the misassemblies aren't inflating the contiguity. Some of them may be real differences in this resequenced line.

**BUSCO**, run on my assembly and, as a control, on the reference genome with the same settings:

```bash
caffeinate -i busco -i 04_assembly/scaffolds_min500.fasta -m genome -l eurotiales_odb12 -c 6 \
  -o busco_spades --out_path /Volumes/GENOMICS/busco --download_path /Volumes/GENOMICS/busco_downloads
```

| | My assembly | Reference |
| --- | --- | --- |
| Complete | 97.3% | 97.0% |
| Single-copy | 97.0% | 96.8% |
| Duplicated | 0.2% | 0.2% |
| Fragmented | 0.4% | 0.4% |
| Missing | 2.3% | 2.6% |
| Internal stop codons (E) | 3.3% | 3.3% |

I ran the reference control because of that 3.3% "E" value (BUSCO genes with internal stop codons). Since the reference shows exactly the same value, it comes from how BUSCO's Miniprot step predicts *Aspergillus* genes, not from errors in my assembly.

### 5. Repeat identification and masking

RepeatModeler builds a repeat library from the genome itself, and RepeatMasker then marks those repeats so the gene predictor doesn't treat transposon genes as fungal genes.

```bash
BuildDatabase -name nidulans_db "$PROJECT/04_assembly/scaffolds_min500.fasta"
caffeinate -i RepeatModeler -database nidulans_db -threads 6 -LTRStruct 2>&1 | tee repeatmodeler.log
```

This took 48 minutes and found 42 repeat families. The classification step failed at first because the conda RepeatMasker doesn't come with the Dfam database. I downloaded only the two Dfam 4.0 files that are actually needed (`dfam40.0.h5` and `dfam40.curated.consensus.0.h5`, about 90 MB compressed), pointed `FAMDB_DATA_DIR` in `famdb.conf` at them, and reran just the classifier:

```bash
RepeatClassifier -consensi nidulans_repeats_unclassified.fa -stockholm nidulans_repeats.stk
```

| Class | Families |
| --- | --- |
| Unknown | 26 |
| DNA/TcMar-Fot1 | 5 |
| LTR/Gypsy | 3 |
| Satellite, LINE/I | 2 each |
| DNA/TcMar-Tc1, LTR/Copia, LINE, SINE/tRNA | 1 each |

Interestingly, RepeatModeler's structural LTR search found no LTR elements, yet four families were classified as LTR retrotransposons by similarity. Complete LTR elements rarely survive a short-read assembly intact, but their fragments are still recognisable.

Soft-masking:

```bash
caffeinate -i RepeatMasker -lib "$PWD/nidulans_repeats_classified.fa" \
  -xsmall -gff -pa 1 -dir masked \
  "$PROJECT/04_assembly/scaffolds_min500.fasta" 2>&1 | tee repeatmasker.log
```

| Repeat type | % of genome |
| --- | --- |
| Unclassified | 0.88 |
| Simple repeats | 0.50 |
| DNA transposons | 0.32 |
| LTR elements | 0.31 |
| LINEs | 0.15 |
| Low complexity | 0.09 |
| Satellites | 0.03 |
| SINEs | 0.01 |
| **Total masked** | **2.28** |

2.28% is low, as expected for this compact genome. To confirm the masking really worked, I counted lowercase bases in the output:

```bash
grep -v ">" masked/scaffolds_min500.fasta.masked | tr -cd 'acgtn' | wc -c   # 672252
```

That matches RepeatMasker's "bases masked" figure exactly.

### 6. Gene prediction

AUGUSTUS ships with a trained *A. nidulans* model, so I used it directly, without any RNA-seq or protein evidence:

```bash
caffeinate -i augustus --species=aspergillus_nidulans --softmasking=1 \
  --gff3=on --uniqueGeneId=true \
  07_annotation/nidulans_masked.fasta > 07_annotation/augustus.gff3
```

This gave **9,350 genes**. My first count said 28,050, which turned out to be my mistake: `awk` splits on spaces by default, so the `# start gene` and `# end gene` comment lines also matched `$3 == "gene"`. With `awk -F'\t'` the count is correct. I also wrote `scripts/filter_augustus_gff.py` to remove low-confidence short genes, but it only removed 12, so I kept the unfiltered set.

To check the predictions, I extracted the proteins and compared them with the reference proteome:

```bash
gffread -g 07_annotation/nidulans_masked.fasta -y 07_annotation/augustus_proteins.faa -S 07_annotation/augustus.gff3

diamond makedb --in $REF/protein.faa -d 07_annotation/ref_proteins
diamond blastp -q 07_annotation/augustus_proteins.faa -d 07_annotation/ref_proteins \
  -o 07_annotation/augustus_vs_ref.tsv \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
  --max-target-seqs 1 --evalue 1e-10 --threads 6

# keep matches covering at least half of both proteins
awk -F'\t' '$4/$5 >= 0.5 && $4/$6 >= 0.5' 07_annotation/augustus_vs_ref.tsv > 07_annotation/augustus_vs_ref_cov50.tsv
```

The `-S` matters: without it gffread writes stop codons as `.`, and DIAMOND hung for hours trying to load the file.

| | Count | % |
| --- | --- | --- |
| My predictions matching a reference protein | 8,752 / 9,350 | 93.6 |
| Reference proteins recovered | 8,700 / 10,453 | 83.2 |
| Identical to the reference (100% identity, same length) | 4,438 | 47.5 of predictions |
| 95% identity or more | 7,730 | 82.7 of predictions |

Because it's the same species, anything below 100% identity is mostly a gene-structure difference (a missed exon, a different start codon), which is what you'd expect from ab initio prediction without evidence. Neighbouring genes in my assembly also matched neighbouring genes in the reference, so gene order is conserved.

### 7. Biosynthetic gene clusters

I ran the fungal version of antiSMASH on its web server ([fungiSMASH](https://fungismash.secondarymetabolites.org), antiSMASH 8.0.4) with my own gene models. Two preparation steps were needed: a clean GFF3 with gene lines, and short scaffold names in both files.

```bash
gffread 07_annotation/augustus.gff3 --keep-genes -o - \
  | sed -E 's/NODE_([0-9]+)_length_[0-9]+_cov_[0-9]+\.[0-9]+/scaffold_\1/g' \
  > 08_antismash/nidulans_genes.gff3

sed -E 's/NODE_([0-9]+)_length_[0-9]+_cov_[0-9]+\.[0-9]+/scaffold_\1/g' \
  04_assembly/scaffolds_min500.fasta > 08_antismash/nidulans_genome.fasta
```

Settings: relaxed strictness, all extra features on (KnownClusterBlast against MIBiG 4.0, ClusterBlast, SubClusterBlast, Pfam, GO terms). I then summarised the 66 region GenBank files with `scripts/summarise_antismash.py` (Biopython):

```bash
python3 scripts/summarise_antismash.py antismash_results antismash_regions.tsv
```

| Cluster type | Regions |
| --- | --- |
| T1PKS | 17 |
| Terpene | 14 |
| Hybrid (two or more types) | 11 |
| NRPS-like | 8 |
| NRPS | 7 |
| Indole | 3 |
| Betalactone | 3 |
| Terpene-precursor | 2 |
| Isocyanide | 1 |
| **Total** | **66** (30 at a scaffold edge) |

Some of the known *A. nidulans* clusters came out clearly:

| Region | Type | Closest known cluster | Similarity | Comment |
| --- | --- | --- | --- | --- |
| 4.2 | T1PKS | emodin / monodictyphenone / shamixanthone | High | |
| 8.1 | T1PKS | asperthecin | High | |
| 10.1 | NRPS + T1PKS | emericellamide A/B | High | 92.5 kb on one scaffold |
| 104.2 | Terpene | clavaric acid | High | |
| 1.1 | NRPS | N-acetyltryptophan | Medium | |
| 106.1 | T1PKS | sterigmatocystin / aflatoxin | Low | Cut at both scaffold ends; core genes in order, flanking genes missing |
| 143.1 | NRPS | penicillin | Low | All three core genes seem present; the low score comes from extra flanking genes in the MIBiG entries |
| 3.1 | T1PKS | none | | A PKS with no characterised match |

The sterigmatocystin and penicillin results were a useful lesson: a "Low" similarity score can mean a fragmented cluster or simply a broadly defined reference entry, and only the gene-by-gene view tells you which.

## Methods summary

Illumina HiSeq 2500 paired-end reads (2 x 101 bp; SRA SRR4236261) from *Aspergillus nidulans* FGSC A4 were quality-checked with FastQC and MultiQC and trimmed with fastp v1.3.6 (TruSeq adapter removal, first base of read 1 removed, 3′ trimming at a 4-base window mean of Q20, minimum length 50 bp), retaining 12,254,337 read pairs. Reads were assembled with SPAdes v4.3.0 in isolate mode, and scaffolds shorter than 500 bp were removed. Assembly quality was assessed with QUAST v5.3.0 against the reference genome ASM1142v1 (GCF_000011425.1) and with BUSCO v6.1.0 (eurotiales_odb12), using the reference genome as a control.

Repeat families were identified de novo with RepeatModeler v2.0.9 (with LTR structural analysis), classified against Dfam 4.0 curated families, and used to soft-mask the assembly with RepeatMasker v4.2.4. Genes were predicted ab initio with AUGUSTUS v3.5.0 using the *A. nidulans* parameters. Predicted proteins (gffread v0.12.9) were compared with the reference proteome using DIAMOND v2.2.6 blastp (E-value ≤ 1e-10, best hit, ≥ 50% coverage of both proteins). Biosynthetic gene clusters were detected with fungiSMASH (antiSMASH v8.0.4, relaxed strictness, KnownClusterBlast against MIBiG 4.0).

## Software versions

| Tool | Version |
| --- | --- |
| fastp | 1.3.6 |
| SPAdes | 4.3.0 |
| QUAST | 5.3.0 |
| BUSCO | 6.1.0 |
| RepeatModeler | 2.0.9 |
| RepeatMasker | 4.2.4 (RMBlast 2.17.1) |
| Dfam / FamDB | 4.0 / 3.0.0 |
| AUGUSTUS | 3.5.0 |
| gffread | 0.12.9 |
| DIAMOND | 2.2.6 |
| Biopython | 1.88 |
| antiSMASH (fungiSMASH) | 8.0.4, MIBiG 4.0 |
| FastQC, MultiQC, NCBI datasets | to be added |

## Problems I ran into

Most of the time went into practical issues rather than the analysis itself. In the order I hit them:

| Problem | Cause | What fixed it |
| --- | --- | --- |
| `PermissionError: Operation not permitted` | macOS privacy protection on the Downloads folder | Gave Terminal Full Disk Access |
| curl error 18 halfway through downloads | Dropped FTP connections | HTTPS plus `curl -C -` in an `until` loop |
| `spades.py: command not found` inside its environment | Package never installed on Apple Silicon | Recreated the environment with `CONDA_SUBDIR=osx-64` |
| SPAdes stopped with error code 66 | Internal disk full | Moved output to an external drive (reformatted from NTFS, which macOS mounts read-only, to exFAT) |
| Only 266 MB free after SPAdes | macOS swap files | `conda clean --all`, emptied Trash, restarted |
| QUAST: `No module named 'distutils'` | Python 3.13 removed distutils | Installed setuptools in the environment |
| BUSCO environment "solving" for 3 hours | conda 23.7's slow solver | Used BUSCO already in another environment; mamba from then on |
| Solver offered RepeatModeler 1.0.11 | Old version selected | Asked for `"repeatmodeler>=2"` explicitly |
| RepeatClassifier / RepeatMasker: FamDB not found | No Dfam database in the conda package | Downloaded two Dfam 4.0 partitions, set `FAMDB_DATA_DIR` |
| `grep -P` not recognised | macOS grep has no Perl regex | Used awk |
| 28,050 "genes" | awk splitting comment lines on spaces | `awk -F'\t'` for anything GFF |
| DIAMOND stuck at "Loading query sequences" | `.` as stop codon in gffread output | `gffread -S` |

What I'd tell myself at the start: check disk space before every long run, test each install with a dry run, and check every important number a second way.

## Limitations

- Short reads break the assembly at repeats. That affects 30 of the 66 cluster regions, and it means the repeat content is probably underestimated.
- Gene models come from ab initio prediction alone, so about half differ from the reference in structure.
- The 366 kb of unaligned sequence hasn't been looked at yet (mitochondrial DNA, rDNA and contamination are the obvious candidates).

## Next steps

- [ ] Add FastQC, MultiQC and NCBI datasets versions
- [ ] BLAST the unaligned scaffolds
- [ ] Plot the antiSMASH region table in R
- [ ] Functional annotation of the predicted proteins (eggNOG-mapper or InterProScan)
- [ ] Apply the pipeline to my marine isolates, with Funannotate on an HPC system and long reads where possible
