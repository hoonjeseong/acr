# ACR (Additional Clustering Refiner)

__Additional Clustering Refiner (ACR)__, which regroups the contigs of the metagenome-assembled genomes (MAGs) using k-means clustering algorithm, to obtain MAGs with better quality. 

----
#### __Usage__
```
Usage: acr.py -g [bin folder] -c [coverage file] -o [output]

Options:
  -h, --help            show this help message and exit
  -g GENOME, --genome=GENOME
                        genome file path
  -o OUTPUT, --output=OUTPUT
                        output folder
  -c COVERAGE, --coverage=COVERAGE
                        contig coverage file
  -e EXTENSION, --extention=EXTENSION
                        genome file extension |default: fa
  -p PREFIX, --prefix=PREFIX
                        prefix |default: refine
  -s SIZE, --size=SIZE  MAG size filter |default: 500k
  -t PROCESS, --process=PROCESS
                        number of workers | default = 1
  -b BYPASS, --bypass=BYPASS
                        bypass prodigal - hmmsearch | default = N
  -j JGI, --from_jgi_cov=JGI
                        please insert Y or N (Y=using jgi coverage file from MetaBAT2) |
                        default = N
  -m GMESEUK, --run_gmesEuk=GMESEUK
                        gene prediction with gmes for Euk (Y or N) | default =
                        Y
  --target=TARGET       refiner target (Prok or Euk) | default = Both
                        (Prokaryote and Eukaryote)
```

#### __Examples__
`python acr.py -g test_SHIPPO/bin/ -o test_SHIPPO/refine -c test_SHIPPO/cov.txt`

The example data are genome bins from [Seong, et. al's] previous ocean microbiome study.

[Seong, et. al's]:https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-022-01340-w

----
#### __Requirements__

1. Please download core gene database file to ACR path. Please check thecommands below.

```
#in the [ACR path]
wget -O data.tar.gz https://figshare.com/ndownloader/files/41282157
tar -zxvf data.tar.gz
```

2. Check the Python library

ACR requires python 3.7 or above.

Here are the python library requirements

- eukcc==0.2
- EukRep==0.6.7
- scikit-learn==0.19.2
- pandas==1.1.2
- kmeans1d==0.3.1

```
1. conda (or use mamba) create -n ACR -c bioconda python=3.7 eukcc==0.2 EukRep==0.6.7 scikit-learn==0.19.2 pandas==1.1.2
2. conda activate ACR
3. pip install kmeans1d
```

3. To run ACR, the absolute paths of prodigal and hmmsearch must be written in the program.txt file as follows:

```
e.g)
prodigal:[/usr/bin/prodigal]
hmmsearch:[/usr/bin/hmmsearch]
hmmalign:[add absolute path]
hmmpress:[add absolute path]
pplacer:[add absolute path] #https://matsen.fhcrc.org/pplacer/
guppy:[add absolute path] #https://matsen.fhcrc.org/pplacer/
EukRep:[add absolute path] #https://github.com/patrickwest/EukRep
runGMES:[add absolute path] #please install GeneMark-ES [http://exon.gatech.edu/GeneMark/license_download.cgi]
```
#### __The overall schematic of the ACR refinement approach with detailed steps__

#### __Benchmark with CAMI 1 high complexity data__

#### __Benchmark with CAMI II rhizosphere data__


