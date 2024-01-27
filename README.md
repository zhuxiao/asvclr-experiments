# ASVCLR experiments 

## Prerequisites

```sh
# Get ASVCLR 
$ wget -c https://github.com/zhuxiao/asvclr/releases/download/1.3.0/asvclr_1.3.0.tar.xz
$ tar -xf asvclr_1.3.0.tar.xz
$ cd asvclr_1.3.0/
$ ./auto_gen.sh
# Or get from github
$ git clone https://github.com/zhuxiao/asvclr.git
$ tar -xf asvclr_1.3.0.tar.xz
$ cd asvclr_1.3.0/
$ ./auto_gen.sh
```

And the binary file `asvclr` will be output into the folder `bin` in this package directory.

```sh
# Get svdss pbsv sniffles svim debreak and samtools
$ conda install svdss=1.0.5 pbsv=2.9.0 sniffles=2.2 svim=1.4.2 debreak=1.0.2 samtools  
# We need ngmlr v0.2.7 to align fasta or fastq with reference
$ wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz
$ tar xvzf ngmlr-0.2.7-linux-x86_64.tar.gz
$ cd ngmlr-0.2.7/
$ mkdir build ; cd build
$ cmake ..
$ make
# We also need sratoolkit to download PacBio CCS data.
# centOS
$ wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-centos_linux64.tar.gz
$ tar -zxvf sratoolkit.3.0.10-centos_linux64.tar.gz
$ cd sratoolkit.3.0.10-centos/bin/
# ubuntu
$ wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-ubuntu64.tar.gz
$ tar -zxvf sratoolkit.3.0.10-ubuntu64.tar.gz
$ cd sratoolkit.3.0.10-ubuntu64
```

And the binary file `prefetch`、 `fastq-dump`  and `fasterq-dump` will be output into the folder `bin` in this package directory.

We used  [SV_STAT](https://github.com/zhuxiao/sv_stat) to evaluate variant calling results.

```sh
$ wget -c https://github.com/zhuxiao/sv_stat/releases/download/0.8.0/sv_stat_0.8.0.tar.xz
$ tar -xf sv_stat_0.8.0.tar.xz
$ cd sv_stat_0.8.0/
$ ./autogen.sh
```

And the binary file `sv_stat` will be output into the folder `bin` in this package directory.

### Data

The reference used in our experiment is hs37d5.

#### Download reference

```sh
# Reference
$ wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
# Extract chromosomes from 1 to 22 and X and Y
$ samtools faidx hs37d5.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > hs37d5.fa
```

## HG002

Download [HG002 PacBio CCS](https://www.ncbi.nlm.nih.gov/sra/SRX5327410) data and convert them into bam files using samtools and create indexes. For convenience, we provide a shell script and a list of accession to help you obtain the fastq file (see `doc` folder). Significantly, you need to ensure that the file and the script are in the same folder. 

```sh
$ ./autofq.sh
$ ngmlr -t 32 --rg-id na24385_pb_ccs -r hs37d5.fa -q SRR885_whole.fastq -o H
G002_pacbio_ccs.sam
$ samtools view -bSh -@ 32 HG002_pacbio_ccs.sam > HG002_pacbio_ccs.bam
$ samtools sort -@ 32 -o HG002_pacbio_ccs_sorted.bam HG002_pacbio_ccs.bam
$ samtools index -@ 32 HG002_pacbio_ccs_sorted.bam HG002_pacbio_ccs_sorted.bai
# remove fastq to save storage space
$ rm -rf HG002_pacbio_ccs.sam HG002_pacbio_ccs.bam
```

#### Call structural variation

```sh
# ASVCLR
$ asvclr all -t 32 -o out_asvclr -m 20 -M 50000 hs37d5.fa HG002_PacBio_CCS_sorted.bam
```

More  detailed usage of ASVCLR can be obtained from Github([ASVCLR](https://github.com/zhuxiao/asvclr)).

You can get variant detection results in folder `4_results` and variant detection results are reported in VCF file format in this file folder: `genome_variants.vcf`.

```sh
# SVDSS
$ SVDSS index --threads 32 --fasta hs37d5.fa --index hs37d5.bwt
$ SVDSS smooth --threads 32 --bam HG002_pacbio_ccs_sorted.bam --reference hs37d5.fa --workdir $PWD
$ SVDSS search --threads 32 --index hs37d5.bwt --bam smoothed.selective.bam --workdir $PWD
$ SVDSS assemble --threads 32 --workdir $PWD --batches 9 
$ SVDSS call --threads 32 --min-sv-length 20 --workdir $PWD --bam smoothed.selective.bam --reference hs37d5.fa
# Debreak
$ debreak --thread 32 --min_size 20 --min_support 5 --bam HG002_pacbio_ccs_sorted.bam --outpath output_debreak --rescue_large_ins --poa --ref hs37d5.fa 
# SVIM
$ svim alignment --min_sv_size 20 output_svim HG002_pacbio_ccs_sorted.bam hs37d5.fa
# pbsv
$ pbsv discover -s HG002_30X_CCS HG002_pacbio_ccs_sorted.bam HG002_pacbio_ccs_sorted.svsig.gz
$ pbsv call -j 32 --ccs hs37d5.fa HG002_pacbio_ccs_sorted.svsig.gz output_pbsv.vcf
# cuteSV
$ cuteSV -t 32 -s 2 HG002_pacbio_ccs_sorted.bam hs37d5.fa output_cutesv.vcf $PWD
# Sniffles2
$ sniffles --input HG002_pacbio_ccs_sorted.bam --vcf output_min20.vcf --reference hs37d5.fa --threads 32 --minsupport 5 --minsvlen 20
```

You can get the following four results:

* **SVDSS** : `svs_poa.vcf`
* **DeBreak** : `output_debreak.vcf`
* **pbsv** : `output_pbsv.vcf`
* **Sniffles2** : `output_sniffles.vcf`
* **SVIM** : `variants.vcf`
* **cuteSV** : `output_cutesv.vcf`

#### CMRG analysis

To compare the results of calling variation to the CMRG callset :

```sh
# Download the CMRG callset an decompression
$ wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh37/StructuralVariant/HG002_GRCh37_CMRG_SV_v1.00.vcf.gz
$ gunzip HG002_GRCh37_CMRG_SV_v1.00.vcf.gz
# Run sv_stat against the CMRG callset and SV_STAT can evaluate multiple callsets simultaneously.
$ sv_stat -o output_CMRG -T "ASVCLR;SVDSS;DeBreak;pbsv;Sniffles2;SVIM;cuteSV" genome_variants.vcf svs_poa.vcf output_debreak.vcf output_pbsv.vcf output_sniffles.vcf variants.vcf output_cutesv.vcf HG002_GRCh37_CMRG_SV_v1.00.vcf hs37d5.fa 
```

In generally, the results are saved to folders under their respective tool names in `output_CMRG`.

The results of this experiment are shown in table:

|           | TP_user | TP_benchmark |  Recall  |
| :-------: | :-----: | :----------: | :------: |
|  ASVCLR   |   180   |     187      | 0.795745 |
|   SVDSS   |   212   |     209      | 0.889326 |
|  DeBreak  |   125   |     132      | 0.561702 |
|   pbsv    |   196   |     202      | 0.856574 |
| Sniffles2 |   211   |     214      | 0.910638 |
|   SVIM    |   230   |     220      | 0.936170 |
|  cuteSV   |   183   |     194      | 0.825532 |

More detailed result information can be found in `evaluation_report.html`.

#### GIAB analysis

```sh
# Get GIAB VCF Tier 1 
$ wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
gunzip HG002_SVs_Tier1_v0.6.vcf.gz
$ sv_stat -o output_giab_Tier1 -T "ASVCLR;SVDSS;DeBreak;pbsv;Sniffles2;SVIM;cuteSV" genome_variants.vcf svs_poa.vcf output_debreak.vcf output_pbsv.vcf output_sniffles.vcf variants.vcf output_cutesv.vcf HG002_SVs_Tier1_v0.6.vcf hs37d5.fa
```

In generally, the results are saved to folders under their respective tool names in `output_giab_Tier1`.

The results of analysis are shown in table:

|           | Precision | F1-score |  Recall  |
| :-------: | :-------: | :------: | :------: |
|  ASVCLR   | 0.793446  | 0.660873 | 0.566260 |
|   SVDSS   | 0.803384  | 0.583478 | 0.458088 |
|  DeBreak  | 0.817189  | 0.674359 | 0.574029 |
|   pbsv    | 0.652436  | 0.56085  | 0.491812 |
| Sniffles2 | 0.777225  | 0.670102 | 0.588932 |
|   SVIM    | 0.560031  | 0.593852 | 0.632019 |
|  cuteSV   | 0.829295  | 0.640409 | 0.521605 |

More detailed result information can be found in `evaluation_report.html`.

## Contact

If you have problems or some suggestions, please contact: zhuxiao_hit@yeah.net without hesitation. 

---- Enjoy !!! -----
