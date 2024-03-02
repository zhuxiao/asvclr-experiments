# ASVCLR experiments

ASVCLR-experiment introduces the general process of the experiment. This includes prerequisites, using different tools for variation detection, and comparing the detection results of ASVCLR with other tools. The detailed usages and information of ASVCLR can be found at GitHub ([ASVCLR](https://github.com/zhuxiao/asvclr)).

## Prerequisites

```sh
# Get ASVCLR 
$ wget -c https://github.com/zhuxiao/asvclr/releases/download/1.4.0/asvclr_1.4.0.tar.xz
$ tar -xf asvclr_1.4.0.tar.xz
$ cd asvclr_1.4.0/
$ ./auto_gen.sh
# Or get from github 
$ git clone https://github.com/zhuxiao/asvclr.git
$ cd asvclr
$ ./auto_gen.sh
$ cd bin && ln -s asvclr /usr/local/bin
```

And the binary file `asvclr` will be output into the folder `bin` in this package directory and link it to the `bin` folder under the username folder, or to the `bin` directory under the `usr` folder that requires higher permissions.

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
$ ln -s prefetch /home/usrname/bin
$ ln -s fastq-dump /home/usrname/bin
$ ln -s fasterq-dump /home/usrname/bin
# ubuntu
$ wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-ubuntu64.tar.gz
$ tar -zxvf sratoolkit.3.0.10-ubuntu64.tar.gz
$ cd sratoolkit.3.0.10-ubuntu64/bin/
$ ln -s prefetch /usr/local/bin
$ ln -s fastq-dump /usr/local/bin
$ ln -s fasterq-dump /usr/local/bin
```

And the binary file `prefetch`、 `fastq-dump` and `fasterq-dump` will be output into the folder `bin` in this package directory,and it is necessary to link to the `bin` folder under the username folder, or to the `bin` directory under the `usr` folder that requires higher permissions.

We used  [SV_STAT](https://github.com/zhuxiao/sv_stat) to evaluate variant calling results.

```sh
$ wget -c https://github.com/zhuxiao/sv_stat/releases/download/0.9.0/sv_stat_0.9.0.tar.xz
$ tar -xf sv_stat_0.9.0.tar.xz
$ cd sv_stat_0.9.0/
$ ./autogen.sh
$ ln -s 
```

And the binary file `sv_stat` will be output into the folder `bin` in this package directory.

### Data

The reference used in our experiment is hs37d5.

#### Download reference

```sh
# Reference
$ wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
$ gunzip hs37d5.fa.gz
```

## HG002

Download [HG002 PacBio CCS](https://www.ncbi.nlm.nih.gov/sra/SRX5327410) data and convert them into bam files using samtools and create indexes. For convenience, we provide a shell script and a list of accession to help you obtain the fastq file (see `script` folder). Significantly, you need to ensure that the list file and the shell script are in the same folder. 

```sh
$ ./prefetch_fastq.sh SRR_Acc_List.txt SRR885_whole.fastq
$ ngmlr -t 32 --rg-id na24385_pb_ccs -r hs37d5.fa -q SRR885_whole.fastq -o H
G002_pacbio_ccs.sam
$ samtools view -bSh -@ 32 HG002_pacbio_ccs.sam > HG002_pacbio_ccs.bam
$ samtools sort -@ 32 -o HG002_pacbio_ccs_sorted.bam HG002_pacbio_ccs.bam
$ samtools index -@ 32 HG002_pacbio_ccs_sorted.bam HG002_pacbio_ccs_sorted.bai
# remove fastq to save storage space
$ rm -rf SRR885_whole.fastq HG002_pacbio_ccs.sam HG002_pacbio_ccs.bam
```

#### Call structural variation

```sh
# ASVCLR
$ asvclr all -t 32 -o out_asvclr -m 20 -M 50000 hs37d5.fa HG002_PacBio_CCS_sorted.bam
```

More  detailed usage of ASVCLR can be obtained from GitHub ([ASVCLR](https://github.com/zhuxiao/asvclr)).

You can get variant detection results in folder `4_results` and variant detection results are reported in VCF file format in this file folder: `genome_variants.vcf`.

```sh
# SVDSS
$ SVDSS index --threads 32 --reference hs37d5.fa --index hs37d5.chroms.fmd
$ SVDSS smooth --threads 32 --reference hs37d5.fa --bam HG002_pacbio_ccs_sorted.bam > HG002_pacbio_ccs_sorted_smoothed.bam
$ samtools index HG002_pacbio_ccs_sorted_smoothed.bam
$ SVDSS search --threads 32 --index hs37d5.fmd --bam HG002_pacbio_ccs_sorted_smoothed.bam > specifics.txt
$ SVDSS call --threads 32 --reference hs37d5.fa --bam HG002_pacbio_ccs_sorted_smoothed.bam --sfs specifics.txt > output_svdss.vcf
# DeBreak
$ debreak --thread 32 --min_size 20 --bam HG002_pacbio_ccs_sorted.bam --outpath output_debreak --rescue_large_ins --poa --ref hs37d5.fa 
$ cd output_debreak && mv debreak.vcf output_debreak.vcf
# SVIM
$ svim alignment --min_sv_size 20 output_svim HG002_pacbio_ccs_sorted.bam hs37d5.fa
$ cd output_svim && mv variants.vcf output_svim.vcf
# pbsv
$ pbsv discover -s HG002_30X_CCS HG002_pacbio_ccs_sorted.bam HG002_pacbio_ccs_sorted.svsig.gz
$ pbsv call -j 32 --ccs hs37d5.fa HG002_pacbio_ccs_sorted.svsig.gz output_pbsv.vcf
# cuteSV
$ cuteSV -t 32 -s 2 HG002_pacbio_ccs_sorted.bam hs37d5.fa output_cutesv.vcf $PWD
# Sniffles2
$ sniffles --input HG002_pacbio_ccs_sorted.bam --vcf output_sniffles.vcf --reference hs37d5.fa --threads 32 --minsvlen 20
```

You can get the following four results:

* **SVDSS** : `output_svdss.vcf`
* **DeBreak** : `output_debreak.vcf`
* **pbsv** : `output_pbsv.vcf`
* **Sniffles2** : `output_sniffles.vcf`
* **SVIM** : `output_svim.vcf`
* **cuteSV** : `output_cutesv.vcf`

For convenience, the variation detection results of all tools are included in the `results` folder. Due to the large size of `output_debreak.vcf`, we chose to use bcftools to split and upload.

#### CMRG analysis

To compare the results of calling variation to the CMRG callset :

```sh
# Download the CMRG callset an decompression
$ wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh37/StructuralVariant/HG002_GRCh37_CMRG_SV_v1.00.vcf.gz
$ gunzip HG002_GRCh37_CMRG_SV_v1.00.vcf.gz
# Run sv_stat against the CMRG callset and SV_STAT can evaluate multiple callsets simultaneously.
$ sv_stat -o output_CMRG -T "ASVCLR;SVDSS;DeBreak;pbsv;Sniffles2;SVIM;cuteSV" genome_variants.vcf output_svdss.vcf output_debreak.vcf output_pbsv.vcf output_sniffles.vcf output_svim.vcf output_cutesv.vcf HG002_GRCh37_CMRG_SV_v1.00.vcf hs37d5.fa 
```

In generally, the results are saved to folders under their respective tool names in `output_CMRG`.

The results of this experiment are shown in table:

|           | SVs_bench | TP_bench | TP_user | Recall |
| :-------- | --------- | :------- | :------ | :----- |
| ASVCLR    | 235       | 228      | 222     | 0.9702 |
| SVDSS     | 235       | 209      | 220     | 0.8893 |
| DeBreak   | 235       | 216      | 203     | 0.9191 |
| pbsv      | 235       | 226      | 227     | 0.9617 |
| Sniffles2 | 235       | 215      | 215     | 0.9148 |
| SVIM      | 235       | 220      | 241     | 0.9361 |
| cuteSV    | 235       | 202      | 190     | 0.8595 |

More detailed experimental results can be seen in the `sv_stat_reports.html` file in the `output_CMRG` folder.

#### GIAB analysis

```sh
# Get GIAB VCF Tier 1 
$ wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
gunzip HG002_SVs_Tier1_v0.6.vcf.gz
$ sv_stat -o output_giab_Tier1 -s 200 -T "ASVCLR;SVDSS;DeBreak;pbsv;Sniffles2;SVIM;cuteSV" genome_variants.vcf output_svdss.vcf output_debreak.vcf output_pbsv.vcf output_sniffles.vcf output_svim.vcf output_cutesv.vcf HG002_SVs_Tier1_v0.6.vcf hs37d5.fa
```

In generally, the results are saved to folders under their respective tool names in `output_giab_Tier1`. In addition, the evaluation results of all tools are uploaded to `output_giab_Tier1` in the `analysis` folder(The `convert` folder is not included because the files is too large and only works during the evaluation process).

The results of analysis are shown in table:

|           | SVs_bench | SVs_user | TP_bench | TP_user |  FP   |  FN   | Precision | Recall   | F1-Score | Time(min) | Peak Mem.(GiB) | Threads |
| :-------: | --------- | -------- | -------- | ------- | :---: | :---: | --------- | -------- | -------- | --------- | -------------- | ------- |
|  ASVCLR   | 74012     | 53170    | 45694    | 44180   | 8990  | 28318 | 0.617386  | 0.830920 | 0.708412 | 29.3      | 17.6           | 32      |
|   SVDSS   | 74012     | 45990    | 34689    | 37221   | 8769  | 39323 | 0.468694  | 0.809328 | 0.593616 | 124.7     | 11.8           | 32      |
|  DeBreak  | 74012     | 50040    | 44078    | 41458   | 7505  | 29934 | 0.595552  | 0.846721 | 0.699266 | 18.2      | 8.6            | 32      |
|   pbsv    | 74012     | 53271    | 44494    | 42927   | 9636  | 29518 | 0.601173  | 0.816677 | 0.692547 | 65.2      | 21.3           | 32      |
| Sniffles2 | 74012     | 54261    | 44658    | 42801   | 10284 | 29354 | 0.603389  | 0.806273 | 0.690231 | 0.8       | 2.9            | 32      |
|   SVIM    | 74012     | 117556   | 48028    | 47230   | 31612 | 25984 | 0.648922  | 0.599046 | 0.622987 | 20.6      | 1.1            | 1       |
|  cuteSV   | 74012     | 45233    | 39685    | 37172   | 6406  | 34327 | 0.536197  | 0.852999 | 0.658475 | 0.9       | 2.1            | 32      |

More detailed experimental results can be seen in the `sv_stat_reports.html` file in the `output_giab_Tier1` folder.

The results can also be observed through the bar chart generated by SV-STAT, as follows:

<div style="text-align: center;">
    <img src="images\result_classification.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="images\evaluation_result.png" alt="Benchmark results between different tools" style="display: inline-block;" width="400"/>
</div>

## Contact

If you have problems or some suggestions, please contact: zhuxiao_hit@yeah.net without hesitation. 

---- Enjoy !!! -----
