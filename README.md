# ACR (Additional Clustering Refiner)

__Additional Clustering Refiner (ACR)__, which regroups the contigs of the metagenome-assembled genomes (MAGs) using k-means clustering algorithm, to obtain MAGs with better quality. 

----
#### __Usage__
```
Usage: refine_bins.py -g [bin folder] -c [coverage file] -o [output]

Options:
  -h, --help            show this help message and exit
  -g GENOME, --genome=GENOME
                        genome file path
  -o OUTPUT, --output=OUTPUT
                        output folder
  -c COVERAGE, --converage=COVERAGE
                        contig coverage file
  -e EXTENSION, --extention=EXTENSION
                        genome file extension |default: fa
  -p PREFIX, --prefix=PREFIX
                        prefix |default: refine
  -s SIZE, --size=SIZE  MAG size filter |default: 500k
  -t THREAD, --thread=THREAD
                        number of workers | default = 1
  -b BYPASS, --bypass=BYPASS
                        bypass prodigal - hmmsearch | default = N
```

#### __Examples__
`python acr.py -g test_set/bin/ -c test_set/cov.txt -o test_result`

----
#### __Require__

ACR requires python 3.6 or above. 

To run ACR, the absolute paths of prodigal and hmmsearch must be written in the program.txt file as follows:

`prodigal:[/usr/bin/prodigal]`

- prodigal
- hmmsearch
- scikit-learn
- kmeans1d
- pandas
