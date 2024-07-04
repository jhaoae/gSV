# gSV
gSV is a general SV detector that combines the alignment-based and assembly-based methods with maximum exact match (MEM) idea.

## Dependency
Dependencies for gSV:

* python3 （tested with version 3.10.13）
* minimap2 (tested with version 2.17)
* [wtdbg2](https://github.com/ruanjue/wtdbg2)
* [copMEM](https://github.com/wbieniec/copmem) (tested with version 2.1)
  
## Installation

### Install from source

```Linux
## Get the source code
git clone https://github.com/jhaoae/gSV.git
cd gSV

## Create a conda environment for gSV
conda env create -f ./environment.yml 

## Install from source
conda activate gSV
```

### Install copMEM
```Linux
git https://github.com/wbieniec/copmem.git
cd copmem
make all
```

## Usage
### General usage
```Linux
## Run gSV
python3 gSV.py -r /path/to/reference.fa -b /path/to/bam.sort.bam -o /path/to/output --mempath /path/to/copMEM2/
```
Required parameters
```Linux
-r REF_PATH           Absolute path to reference (index required in the directory)
-b BAM_PATH           Absolute path to bam file (index required in the directory)
-o OUT_PATH           Absolute path to output
--mempath MEMPATH     Absolute path of copMEM2
```
### Run demo data
The demo data is  `./demo/DUP_INV.demo.bam`, extracted from the simulated complex SVs (DUP+INV, ID4) used in this study. 

`./demo/ground_truth_complexSV.vcf` and `./demo/ground_truth_SV.vcf` record the information of complex SVs and the corresponding sub-SVs of each complex SV.

Assume that the absolute path is `/home/gSV`, and copMEM2 is also installed in this directory.

Download reference genome [hs37d5](https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/)
```
## Run gSV
python3 gSV.py -r /home/gSV/demo/hs37d5.fa -b /home/gSV/demo/DUP_INV.demo.bam -o /home/gSV/DUP_INV --mempath /home/gSV/copMEM2 --Complex True
```


## Output files
The output directory includes:
```
workspace/              Intermediate files
gSV.vcf                 VCF file of both simple and complex SV calls
gSV_complexSV.vcf       VCF file of complex SV calls
```
Complex SVS are recorded as sub-SVs in the gSV.vcf file. 
For example, a complex SV of INV+DUP is recorded in gSV_complexSV.vcf as
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    
1       3413020 SV1     C       <COM>   .       PASS    END=3413984;SVTYPE=COMPLEX:INV:3413020-3413984-965;DUP:3413731-3413984-253;SVLEN=964
```
The corresponding record in the gSV.vcf file is
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    
1       3413020 SV1     C       <INV>   .       PASS    END=3413984;SVTYPE=INV;SVLEN=965
1       3413731 SV2     C       <DUP>   .       PASS    END=3413984;SVTYPE=tandemDUP;SVLEN=253
```

