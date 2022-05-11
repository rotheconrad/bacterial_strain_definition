# bacterial_strain_definition
Contains the code and workflow for the bacterial strain definition paper with Kostas Kostantinidis

# STEP 01: Get the data
---
Project directory: /storage/home/hcoda1/9/rconrad6/p-ktk3-0/04_Recombination
project root: 04_Recombination

I chose to download data from NCBI through the FTP access.
NCBI provides current genome summary for genebank or refseq

Date: Apr 20 2022
In Direcotry : 04_Recombination/01a_Data_Download_Prep

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
```

We want to focus on just archaea and bacteria so we will ignore the files above which include everything
and work with the files listed below:
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
mv assembly_summary.txt NCBI_assembly_summary_bacteria.txt
```

### Description of the assembly summary file

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/README.txt
```
I decided to go with refseq genomes for this project.

Genomes are divided into four assembly levels (Complete, Chromosome, Scaffold, and Contig).

I wrote a Python script to parse this file.
It creates summaries counting the number of genomes at each assembly level for each species.
It also creates files to download genomes for each species at each assembly level.
The download file can be created based on the number of genomes available at each level for each species.
n=10 for complete genomes will list the ftp download only for species with at least 10 complete genomes.
The script also filters for the "latest" version of the genome ignores other genomes versions.

```bash
python ../00c_Scripts/00a_Parse_NCBI_Assembly_Summary.py -i NCBI_assembly_summary_bacteria.txt -p bacteria -n 10
```

Download genomes:
Date: Apr 20 2022
In Directory : 04_Recombination/01b_Genomes
Total Species complete genomes >= 10: 330
Total genomes: 18153

#### Make directories

```bash
mkdir bacteria bacteria/Complete
```

### Download Complete level bacteria genomes to species directories

```bash
while read p; do t=bacteria; l=Complete; n=`echo -e "$p" | cut -f1`; m=`echo -e "$p" | cut -f2`; if [ ! -d ${t}/${l}/$n ]; then mkdir ${t}/${l}/$n; fi; wget -P ${t}/${l}/$n $m; done < ../01a_Data_Download_Prep/bacteria_Complete_ftps.sh
```

### Check we got them all

```bash
while read p; do t=bacteria; l=Complete; n=`echo -e "$p" | cut -f1`; m=`echo -e "$p" | cut -f2`; x=`echo $m | rev | cut -d/ -f1 | cut -d. -f2- | rev`; if [ ! -s ${t}/${l}/${n}/$x ]; then echo $n $x "NOT COMPLETE DOWNLOADING"; wget -P ${t}/${l}/$n $m; fi; done < ../01a_Data_Download_Prep/bacteria_Complete_ftps.sh
```

### Unzip them

```bash
qsub -v f=01b_Genomes/bacteria/*/*/* 00b_PBS/01_gunzip.pbs
```

# STEP 02: Run All vs All fastANI for each species
---

April 21 2022
In Directory: 04_Recombination/01b_Genomes/bacteria

### Run all vs all for each species fastANI

```bash
mkdir 00a_log fastANI_Complete
for d in Complete/*; do n=`basename $d`; qsub -v fDir=$d,oDir=fastANI_Complete,n=$n ../../00b_PBS/02a_fastANI.pbs; done
```

#### concate files

```bash
mkdir 01a_AllvAll_fastANI
for d in fastANI_Complete/*; do n=`basename $d`; cat ${d}/* >> 01a_AllvAll_fastANI/${n}.ani; echo $d; done
cat 01a_ALLvAll_fastANI/*.ani >> fastANI_Complete_All.ani
```

#### plot
```bash
python ../../00c_Scripts/02b_fastANI_scatter_pyGAM.py -i fastANI_Complete_All.ani -s All_species -o fastANI_Complete_All_95_density_pyGAM.pdf -z True -g True
python ../../00c_Scripts/02b_fastANI_scatter_pyGAM.py -i fastANI_Complete_All.ani -s All_species -o fastANI_Complete_All_98_density_pyGAM.pdf -xmin 98 -t 0.5 -z True -g True
```

![Shared genome fraction vs ANI plot for 330 species constrained at 95% ANI.](figures/fastANI_Complete_All_95_density_pyGAM.pdf)

![Shared genome fraction vs ANI plot for 330 species constrained at 98% ANI.](figures/fastANI_Complete_All_98_density_pyGAM.pdf)

#### Range fraction count
```bash
python ../../00c_Scripts/02d_fastANI_fraction_in_range.py -i fastANI_Complete_All.ani -xmin 99.2 -xmax 99.8
```

> Genome pair counts:	(A) Genome pairs in range [99.2, 99.8]: 235527
> 	(B) Total genome pairs: 4455818
> 	(C) Genome pairs >= 95% ANI: 4344982
> 	(D) Genome pairs >= 96% ANI: 4280042
> 	(E) Genome pairs >= 97% ANI: 3637508
> 	(F) Genome pairs >= 98% ANI: 2677076
> 	(G) Genome pairs >= 99% ANI: 934658

> Fraction of A in B-G:	(A) / (B) = 0.0529
> 	(A) / (C) = 0.0542
> 	(A) / (D) = 0.0550
> 	(A) / (E) = 0.0647
> 	(A) / (F) = 0.0880
> 	(A) / (G) = 0.2520

# STEP 03 ## Dorian
---

We can add dorians work here.
