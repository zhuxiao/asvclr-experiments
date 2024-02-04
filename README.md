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
# Get svdss pbsv sniffles debreak cuteSV and samtools
$ conda install svdss=1.0.5 pbsv=2.9.0 sniffles=2.2 debreak=1.0.2 cuteSV=2.1.0 samtools 
# Get SVIM v2.0.0
$ wget -c https://github.com/eldariont/svim/archive/refs/tags/v2.0.0.tar.gz
$ tar -zxvf svim-2.0.0.tar.gz
$ cd svim-2.0.0
$ pip install .
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
$ gunzip hs37d5.fa.gz
# Extract chromosomes from 1 to 22 and X and Y
$ samtools faidx hs37d5.fa 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y > hs37d5.chroms.fa
```

## HG002

Download [HG002 PacBio CCS](https://www.ncbi.nlm.nih.gov/sra/SRX5327410) data and convert them into bam files using samtools and create indexes. For convenience, we provide a shell script and a list of accession to help you obtain the fastq file (see `doc` folder). Significantly, you need to ensure that the list file and the shell script are in the same folder. 

```sh
$ ./prefetch_fastq.sh SRR_Acc_List.txt SRR885_whole.fastq
$ ngmlr -t 32 --rg-id na24385_pb_ccs -r hs37d5.chroms.fa -q SRR885_whole.fastq -o H
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
$ asvclr all -t 32 -o out_asvclr -m 20 -M 50000 hs37d5.chroms.fa HG002_PacBio_CCS_sorted.bam
```

More  detailed usage of ASVCLR can be obtained from Github([ASVCLR](https://github.com/zhuxiao/asvclr)).

You can get variant detection results in folder `4_results` and variant detection results are reported in VCF file format in this file folder: `genome_variants.vcf`.

```sh
# SVDSS
$ SVDSS index --threads 32 --reference hs37d5.chroms.fa --index hs37d5.chroms.fmd
$ SVDSS smooth --threads 32 --reference hs37d5.chroms.fa --bam HG002_pacbio_ccs_sorted.bam > HG002_pacbio_ccs_sorted_smoothed.bam
$ samtools index HG002_pacbio_ccs_sorted_smoothed.bam
$ SVDSS search --threads 32 --index hs37d5.chroms.fmd --bam HG002_pacbio_ccs_sorted_smoothed.bam > specifics.txt
$ SVDSS call --threads 32 --reference hs37d5.chroms.fa --bam HG002_pacbio_ccs_sorted_smoothed.bam --sfs specifics.txt > output_svdss.vcf
# Debreak
$ debreak --thread 32 --min_size 20 --bam HG002_pacbio_ccs_sorted.bam --outpath output_debreak --rescue_large_ins --poa --ref hs37d5.chroms.fa 
$ cd output_debreak && mv debreak.vcf output.debreak.vcf
# SVIM
$ svim alignment --min_sv_size 20 output_svim HG002_pacbio_ccs_sorted.bam hs37d5.chroms.fa
$ cd output_svim && mv variants.vcf output_svim.vcf
# pbsv
$ pbsv discover -s HG002_30X_CCS HG002_pacbio_ccs_sorted.bam HG002_pacbio_ccs_sorted.svsig.gz
$ pbsv call -j 32 --ccs hs37d5.chroms.fa HG002_pacbio_ccs_sorted.svsig.gz output_pbsv.vcf
# cuteSV
$ cuteSV -t 32 -s 2 HG002_pacbio_ccs_sorted.bam hs37d5.chroms.fa output_cutesv.vcf $PWD
# Sniffles2
$ sniffles --input HG002_pacbio_ccs_sorted.bam --vcf output_sniffles.vcf --reference hs37d5.chroms.fa --threads 32 --minsvlen 20
```

You can get the following four results:

* **SVDSS** : `output_svdss.vcf`
* **DeBreak** : `output_debreak.vcf`
* **pbsv** : `output_pbsv.vcf`
* **Sniffles2** : `output_sniffles.vcf`
* **SVIM** : `output_svim.vcf`
* **cuteSV** : `output_cutesv.vcf`

#### CMRG analysis

To compare the results of calling variation to the CMRG callset :

```sh
# Download the CMRG callset an decompression
$ wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh37/StructuralVariant/HG002_GRCh37_CMRG_SV_v1.00.vcf.gz
$ gunzip HG002_GRCh37_CMRG_SV_v1.00.vcf.gz
# Run sv_stat against the CMRG callset and SV_STAT can evaluate multiple callsets simultaneously.
$ sv_stat -o output_CMRG -T "ASVCLR;SVDSS;DeBreak;pbsv;Sniffles2;SVIM;cuteSV" genome_variants.vcf output_svdss.vcf output_debreak.vcf output_pbsv.vcf output_sniffles.vcf output_svim.vcf output_cutesv.vcf HG002_GRCh37_CMRG_SV_v1.00.vcf hs37d5.chroms.fa 
```

In generally, the results are saved to folders under their respective tool names in `output_CMRG`.

The results of this experiment are shown in table:

|           | TP_bench | TP_user | Recall |
| :-------: | :------: | :-----: | :----: |
|  ASVCLR   |   212    |   206   |  90.2  |
|   SVDSS   |   209    |   220   |  88.9  |
|  DeBreak  |   216    |   203   |  91.9  |
|   pbsv    |   226    |   227   |  96.1  |
| Sniffles2 |   215    |   215   |  91.4  |
|   SVIM    |   220    |   241   |  93.6  |
|  cuteSV   |   202    |   190   |  85.9  |

More detailed results information can be found in `evaluation_report.html`.

#### GIAB analysis

```sh
# Get GIAB VCF Tier 1 
$ wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
gunzip HG002_SVs_Tier1_v0.6.vcf.gz
$ sv_stat -o output_giab_Tier1 -T "ASVCLR;SVDSS;DeBreak;pbsv;Sniffles2;SVIM;cuteSV" genome_variants.vcf output_svdss.vcf output_debreak.vcf output_pbsv.vcf output_sniffles.vcf output_svim.vcf output_cutesv.vcf HG002_SVs_Tier1_v0.6.vcf hs37d5.chroms.fa
```

In generally, the results are saved to folders under their respective tool names in `output_giab_Tier1`.

The results of analysis are shown in table:

|           | Precision | Recall | F1-score | TP    | Time(min) | Peak Mem.(GiB)    |
| :-------: | :-------: | :----: | :------: | ----- | --------- | ----------------- |
|  ASVCLR   |   82.8    |  59.5  |   69.3   | 42206 | 29.28     | 17.6 (threads=32) |
|   SVDSS   |   80.0    |  46.3  |   58.7   | 36813 | 124.7     | 11.8 (threads=32) |
|  DeBreak  |   84.0    |  59.0  |   69.3   | 41177 | 18.2      | 8.6 (threads=32)  |
|   pbsv    |   80.4    |  59.4  |   68.3   | 42261 | 65.18     | 21.3(threads=32)  |
| Sniffles2 |   79.4    |  59.5  |   68.0   | 42162 | 0.8       | 2.9(threads=32)   |
|   SVIM    |   58.8    |  64.1  |   61.3   | 46416 | 20.6      | 1.1 (threads=1)   |
|  cuteSV   |   84.5    |  53.0  |   65.1   | 36862 | 0.85      | 2.1(threads=32)   |

More detailed results information can be found in `evaluation_report.html`.

## Contact

If you have problems or some suggestions, please contact: zhuxiao_hit@yeah.net without hesitation. 

---- Enjoy !!! -----
