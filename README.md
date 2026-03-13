# SRA-convert-to-fastq

steps for downloading SRA and convert to fastq.gz

**Dataset 1**

Lai, H., Cheng, X., Liu, Q., Luo, W., Liu, M., Zhang, M., Miao, J., Ji, Z., Lin, G. N., Song, W., Zhang, L., Bo, J., Yang, G., Wang, J., & Gao, W. Q. (2021). Single-cell RNA sequencing reveals the epithelial cell heterogeneity and invasive subpopulation in human bladder cancer. International journal of cancer, 149(12), 2099–2115. https://doi.org/10.1002/ijc.33794

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135337 (for SRA files)
https://www.ebi.ac.uk/ena/browser/view/SRP217277 (for fastq.gz)

```bash
# create a batch script
nano run_sra_batch.sh
```
```bash
# write the batch script
#!/bin/bash
#SBATCH --job-name=GSE135337_robust
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=sra_out_%j.log
#SBATCH --error=sra_err_%j.log

# point to scratch space
cd /global/scratch/$USER/7marchGSE135337
module load sra-toolkit/3.0.9

# could use the dedicated download node or Aspera to speed up downloads for the future

# all IDs in this dataset GSE135337_SRP217277 (Lai et al., 2021),  we can also batch download with e-utilities
SRR_IDS=("SRR12539462" "SRR12539463" "SRR14615558" "SRR9897621" "SRR9897622" "SRR9897623" "SRR9897624" "SRR9897625")

# loop thru IDs
for SRR in "${SRR_IDS[@]}"; do
echo "[$(date +'%Y-%m-%d %H:%M:%S')] starting $SRR"

# download concurrently w/ extraction
echo "downloading $SRR..."
prefetch $SRR --type sra --max-size 200G --output-directory . --resume yes
# define max size otherwise default >20GB downloads blocked
    
# error handling for download
if [ ! -d "./$SRR" ]; then
echo "SRR failed to download"
exit 1
fi

# concurrently extract to FASTQ - extract one at a time
wait 

# extract and clean up
{
echo "Action: Extracting $SRR..."
fasterq-dump "./$SRR" \
--threads 8 \
--split-files \
--outdir . \
--temp .

# for future downloads add a pigz compression at this step to compress directly to fastq.gz

# verify fastq and clean up source file
if ls ${SRR}*.fastq 1> /dev/null 2>&1; then
echo "cleanup source"
rm -rf "./$SRR"
echo "cleanup source done"
else
echo "srr extraction failed"
exit 1
fi
} &
done

# ensure last extraction before script ends 
wait
echo "all samnples extracted"
```

```bash
# check lines against metadata (slow way)
wc -l SRR12539462_*.fastq
# or wc -l SRR*.fastq

# check 'spots' (SEQ) against metadata (fast way)
seqkit stats -j 8 SRR*.fastq
```
```yaml
# output from wc -l ► both fastq files should have same number of lines
1436300884 SRR12539462_1.fastq
1436300884 SRR12539462_2.fastq
2872601768 total

# output from seqkit
[mii] loading StdEnv/2023 seqkit/2.5.1 ...
processed files:  2 / 2 [======================================] ETA: 0s. done
file                 format  type     num_seqs         sum_len  min_len  avg_len  max_len
SRR12539462_1.fastq  FASTQ   DNA   359,075,221   9,335,955,746       26       26       26
SRR12539462_2.fastq  FASTQ   DNA   359,075,221  54,220,358,371      151      151      151
```

```bash
# load metadata for one SRR
vdb-dump --info SRR12539462

# load metadata for entire study (SRPxxxx)
module load edirect/20.9.20231210
module load sra-toolkit/3.0.9

# loop
esearch -db sra -query SRP217277 | efetch -format runinfo | cut -d ',' -f 1 | grep "^SRR" | while read -r SRR; do
  echo -e "\n >>> ACCESSION: $SRR <<<"
  vdb-dump --info "$SRR"
  done
```

```yaml
# output metadata for one SRR
[mii] loading StdEnv/2023 gcc/12.3 sra-toolkit/3.0.9 ...
acc    : SRR12539462
path   : https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR12539462/SRR12539462
size   : 28,997,579,445
type   : Database
platf  : SRA_PLATFORM_ILLUMINA
SEQ    : 359,075,221
SCHEMA : NCBI:align:db:alignment_sorted#1.3
TIME   : 0x000000005f4940f2 (08/28/2020 13:37)
FMT    : FASTQ
FMTVER : 2.9.1
LDR    : latf-load.2.9.1
LDRVER : 2.9.1
LDRDATE: Jun 15 2018 (6/15/2018 0:0)
```

**relation btw SEQ and lines**

* (one FASTQ entry) = 1 spot (SEQ) = 1 read pair = 4 lines
* Total spots (SEQ) x 4 = Total lines 
* 359 075 221 x 4 = 1 436 300 884

```bash
# check head and tail of fastq (8 lines = 2 reads)
head -n 8 SRR12539462_*.fastq
tail -n 8 SRR12539462_*.fastq
```
```sh
# fastq 1 contains barcodes and UMIs
# each line is followed by quality score line (Phred 33 in this case)
# first 4 lines of fastq 1
==> SRR12539462_1.fastq <==
@SRR12539462.1 1 length=26
NATGCCCAGAAGGCCTGATCGTGATT
+SRR12539462.1 1 length=26
#AAAFFF7FAJJFFAJJFJAFAJ-AA
@SRR12539462.2 2 length=26
NAAGTAGTCCAGAGGACTCTAGTACA
+SRR12539462.2 2 length=26
#AA-FFJJJJJJJJFAAJFJJFJFJJ

# fastq 2 has the actual sequence and Phred score (phred 33)
# first 4 lines fastq 2
==> SRR12539462_2.fastq <==
@SRR12539462.1 1 length=151
NAATAAAGAAGAGGCTGCAGAATATGATAAACTTTTGGCCAAGAGAATGAAGGAGGCTAAGGAGAAGCGCCAGGAACAAATTGCGAGGAGACGCAGACTTTCCTCTCTGCGAGCTTCTACTTCTAAGTCTGAATCCAGTCAGAAATAAGAT
+SRR12539462.1 1 length=151
#AAAF--F--J-FF-F<-FFFJFAFJFJJAJ--7JJF<FAJFJFF<JJAA7AF-F<7AJAFA<J<FA--A-FAJA77-<<FJF<FF-AJJ-7AAFFAF<J<--AAAA-<JAJJ7-<JF<7AJ-7--7--7FJ7-77-F<F<F<<--7FA7-
@SRR12539462.2 2 length=151
NAAGAATTTCTTCAGGTTGAATTACCTAGAAGTTTGTCACTGACTTGTGTTCCTGAACTATGACACATGAATGTGTGGGCTAAGAAATAGTTCCTCTTGATAAATAAACAATTAACAAACAAAAAAAAAAAAAAAAAAAAAAAAGAATATA
+SRR12539462.2 2 length=151
#AAFFJJJ<JJFJJJJJJFFJFJJJFFJJJJJJJJJJFJJJJF-FAAJJJJJJJFFFJJJJAAJJJAAJF7AJ7AJFF7JJFAJFFJJJJJJJF<JJJJ7---FA<FJJJA7AFAJF<<-AJJJJJJJJJJJJJJJF<F7FA<A-7-----
```

| File                     | Content                                                         | Length (10X Genomics)|
| ------------------------ | --------------------------------------------------------------- | ---------------------|
| SRR12539462_1.fastq (R1) | Cell barcode (16 bp) + UMI (10-12 bp) + poly-T/adapter (~10 bp) | ~28-91 bp |
| SRR12539462_2.fastq (R2) | Actual cDNA transcript (gene expression data)                   | 91-150 bp |



* _1.fastq has barcode/UMI
* _2.fastq actual genomic sequence and quality score (Phred)
* peek at Phred score to know if it's phred 33 or phred 64
* this one has phred 33
* https://people.duke.edu/~ccc14/duke-hts-2018/bioinformatics/quality_scores.html
* **download complete - data hygiene check in footnotes**

***

**Do the same for dataset 2: a larger cross platform dataset to test integration**

Tran, M. A., Youssef, D., Shroff, S., Chowhan, D., Beaumont, K. G., Sebra, R., Mehrazin, R., Wiklund, P., Lin, J. J., Horowitz, A., Farkas, A. M., Galsky, M. D., Sfakianos, J. P., & Bhardwaj, N. (2024). Urine scRNAseq reveals new insights into the bladder tumor immune microenvironment. The Journal of experimental medicine, 221(8), e20240045. https://doi.org/10.1084/jem.20240045

https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1112449

```bash
# create batch script
nano PRJNA1112449.sh

# use transfer node in the future to speed up this 1TB download

# write batch script 
#!/bin/bash
#SBATCH --job-name=PRJNA1112449
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=sra_out_%j.log
#SBATCH --error=sra_err_%j.log

# set up path
WORKING_DIR="/global/scratch/$USER/7marchPRJNA1112449"
mkdir -p "$WORKING_DIR"
cd "$WORKING_DIR" || exit

# load SRA and eutilities
module load sra-toolkit/3.0.9
module load edirect

# define variable 
PROJECT_ACC="PRJNA1112449"

# use eutils to grab all SRRs
SRR_IDS=($(esearch -db sra -query "$PROJECT_ACC" | efetch -format runinfo | cut -d ',' -f 1 | grep SRR))

# verify SRR numbers - we expect 168
SRR_COUNT=${#SRR_IDS[@]}
echo "total SRR found: $SRR_COUNT"

# loop through SRRs
for SRR in "${SRR_IDS[@]}"; do
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Starting $SRR"

# download
echo "downloading $SRR..."
prefetch "$SRR" --max-size 200G --output-directory . --resume yes
    
# download error handling
if [ ! -d "./$SRR" ] && [ ! -f "./$SRR.sra" ]; then
echo "error: $SRR failed to download."
exit 1
fi

# sequential extraction in background
wait

# for future downloads add a pigz compression at this step to compress directly to fastq.gz

# extract and clean up 
{
echo "extracting $SRR..."
fasterq-dump "./$SRR" \
--threads 8 \
--split-files \
--outdir . \
--temp . \
--skip-technical

# verify and clean 
if ls ${SRR}*.fastq 1> /dev/null 2>&1; then
echo "$SRR extracted. clean up source."
rm -rf "./$SRR"
else
echo "Error: $SRR extraction failed"
exit 1
fi
} & 
done

# wait for last sample before exit 
wait
echo "all samples extracted"

# run job
sbatch PRJNA1112449.sh
```

```bash
# check our process logs
cat sra_out_5415604.log

# output - all SRRs found
Fetching Run IDs for PRJNA1112449...
Total SRRs found: 168
```

```bash
# job complete - hygiene check - we expect 168 x 2 = 336 fastq
ls -1d SRR*.fastq 2>/dev/null | wc -l

# output: only 332. we're missing 2 SRR extractions = 4 fastq files
332

# quick scroll through found 2 SRRs without fastq
drwxr-x---. 2 hpc6140 hpc6140 4.0K SRR29060942
drwxr-x---. 2 hpc6140 hpc6140 4.0K SRR29061015

# check error logs
grep -E "SRR29060942|SRR29061015" sra_out_5415604.log

# error log output
grep -E "SRR29060942|SRR29061015" sra_out_5415604.log
[2026-03-08 01:01:40] Starting SRR29060942
Action: Downloading SRR29060942...
2026-03-08T06:22:42 prefetch.3.0.9: 1) Downloading 'SRR29060942.lite'...
2026-03-08T06:23:10 prefetch.3.0.9: 1) failed to download 'SRR29060942.lite': RC(rcNS,rcNoTarg,rcValidating,rcConnection,rcNotFound) 
Action: Extracting SRR29060942...
Error: SRR29060942 extraction failed
[2026-03-08 06:01:30] Starting SRR29061015
Action: Downloading SRR29061015...
2026-03-08T10:22:31 prefetch.3.0.9: 1) Downloading 'SRR29061015.lite'...
2026-03-08T10:22:59 prefetch.3.0.9: 1) failed to download 'SRR29061015.lite': RC(rcNS,rcNoTarg,rcValidating,rcConnection,rcNotFound) 
Action: Extracting SRR29061015...
Error: SRR29061015 extraction failed
```

**Troubleshoot 2 missing files**
* SRA Toolkit has a preference for SRA Lite where available
* In this case, the SRA Lite files seemed present but unavailable
* So SRA Toolkit just skipped it
* We need to add a code line to grab original SRA if SRA Lite is not available

```bash
# manually get the 2 SRRs with the fallback to original SRA
# create nano script
nano PRJNA1112449troubleshoot.sh

# job script
#!/bin/bash
#SBATCH --job-name=PRJNA1112449troubleshoot
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=10:00:00
#SBATCH --output=sra_out_%j.log
#SBATCH --error=sra_err_%j.log

# pt to the directory
WORKING_DIR="/global/scratch/$USER/7marchPRJNA1112449"
mkdir -p "$WORKING_DIR"
cd "$WORKING_DIR" || exit

module load sra-toolkit/3.0.9
module load edirect

# 2 missing SRRs needed
SRR_IDS=("SRR29061015" "SRR29060942")

for SRR in "${SRR_IDS[@]}"; do
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Starting $SRR"

echo "downloading $SRR"
prefetch "$SRR" --type sra --max-size 200G --output-directory . --resume yes
# specify SRA so we don't revert to the unavailable SRA lite version

if [ ! -d "$SRR" ] && [ ! -f "$SRR.sra" ] && [ ! -f "$SRR.sralite" ]; then
echo "Error: $SRR failed to download."
exit 1
fi

wait

{
echo "extracting $SRR"
fasterq-dump "$SRR" \
--threads 4 \
--split-files \
--outdir . \
--temp . \
--skip-technical

if ls ${SRR}*.fastq 1> /dev/null 2>&1; then
echo "$SRR extracted. cleaning up."
rm -rf "$SRR" "$SRR.sra" "$SRR.sralite"
else
echo "error: $SRR extraction failed."
fi
} &
done 

# & puts this block into the background so next download begins concurrently

wait
echo "[$(date +'%Y-%m-%d %H:%M:%S')] job complete."

# run job script
sbatch PRJNA1112449troubleshoot.sh
```
```bash
# check total SRR fastq now
ls -1d SRR*.fastq 2>/dev/null | wc -l

# output : all present
336

# double check all consecutive numbers present
f=0
for i in {29060919..29061086}
do
    compgen -G "SRR$i*" >/dev/null || { 
        echo "missing: $i"
        f=1
    }
done

if [ $f -eq 0 ]; then
    echo "all present"
fi

# output
all present
```
