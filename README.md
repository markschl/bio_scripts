# Collection of scripts for the analysis of biological sequences

## Installing

Download and make sure that `bio_scripts/py` is added to PATH.

## Obtaining sequences from GenBank

`gb_downloader.py` helps downloading large amounts of sequences from GenBank. Example:

```sh
term='("Cytochrome c oxidase subunit I" OR "Cytochrome c oxidase subunit 1" OR COi OR CO1 OR COXi OR COX1) AND eukaryotes[porgn] AND ("120"[SLEN] : "10000"[SLEN])'
email=my.email@domain.org  # required if using Entrez API
gb_downloader.py -q "$term" -e "$email" -o COI.fasta
```

It is also possible to download sequences for a list of accessions:

```sh
# first download accessions only (-t acc)
term='(internal transcribed spacer OR ITS OR ITS1 OR ITS2) AND Hygrocybe[Organism]'
email=my.email@domain.org  # required if using Entrez API
gb_downloader.py -q "$term" -e "$email" -b 2000 -t acc -o Hygrocybe_acc.txt

# (then, e.g. choose the accessions for which sequences should be downloaded and remove the rest)
gb_downloader.py -a Hygrocybe_acc.txt -e "$email" -o Hygrocybe_ITS.fasta
```

## Obtaining NCBI taxonomy

### Small amounts of data

Given, that we've downloaded all *Hygrocybe* spp. ITS sequences, we also want the corresponding taxonomy.

```sh
email=my.email@domain.org  # required if using Entrez API
# get_taxonomy assumes the first column of a tab-delimited file to contain the accessions.
# Can be changed using -c option.
get_taxonomy.py -i Hygrocybe_acc.txt -e "$email" -o Hygrocybe_tax.txt
```
The start of `Hygrocybe_tax.txt` looks as follows:

```
MW374250.1      no_rank:cellular organisms      superkingdom:Eukaryota  clade:Opisthokonta      kingdom:Fungi   subkingdom:Dikarya      phylum:Basidiomycota    subphylum:Agaricomycotina     class:Agaricomycetes     subclass:Agaricomycetidae       order:Agaricales        family:Hygrophoraceae   genus:Hygrocybe species:Hygrocybe conica        varietas:Hygrocybe conica var. conica
KF545851.1      no_rank:cellular organisms      superkingdom:Eukaryota  clade:Opisthokonta      kingdom:Fungi   subkingdom:Dikarya      phylum:Basidiomycota    subphylum:Agaricomycotina     class:Agaricomycetes     subclass:Agaricomycetidae       order:Agaricales        family:Hygrophoraceae   genus:Hygrocybe species:Hygrocybe conica        varietas:Hygrocybe conica var. conica
KF545852.1      no_rank:cellular organisms      superkingdom:Eukaryota  clade:Opisthokonta      kingdom:Fungi   subkingdom:Dikarya      phylum:Basidiomycota    subphylum:Agaricomycotina     class:Agaricomycetes     subclass:Agaricomycetidae       order:Agaricales        family:Hygrophoraceae   genus:Hygrocybe species:Hygrocybe conica        varietas:Hygrocybe conica var. conica
KF545860.1      no_rank:cellular organisms      superkingdom:Eukaryota  clade:Opisthokonta      kingdom:Fungi   subkingdom:Dikarya      phylum:Basidiomycota    subphylum:Agaricomycotina     class:Agaricomycetes     subclass:Agaricomycetidae       order:Agaricales        family:Hygrophoraceae   genus:Hygrocybe species:Hygrocybe conica        varietas:Hygrocybe conica var. conica
(...)
```

Depending on the kingdom, different taxonomic ranks may be important. In the following, `extract_taxonomy.py` will obtain the names for those ranks.

```sh
ranks=superkingdom,kingdom,phylum,class,order,family,genus,species,varietas
extract_taxonomy.py -r "$ranks" Hygrocybe_tax.txt > Hygrocybe_tax_ranks.txt
```

*Note:* The specified ranks are actually the default, so we would not need to specify the `-r` argument.

### Huge amounts of data

With the large COI query, querying the Entrez API would take too much time. But there is another option: NCBI creates database dumps of the taxonomy database and uploads them here: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy. If `ftp` is configured as source, these files are downloaded and directly parsed while downloading, looking for the accessions we need information for. This is the fastest way with thousands or even millions of accessions, but it also downloads (and discards) large amounts of unnecessary data. Also, all accessions are currently loaded into memory in order to count them. This could be optimized in the future.

Note that by supplying `-s ftp,entrez`, we use NCBI FTP as primary source, but accessions not found in these files (maybe because the accessions were added to NCBI after the weekly dumps were created) are then searched against the Entrez API to obtain taxonomic lineages.

```sh
email=my.email@domain.org  # not needed for 'ftp' source, but for 'entrez'
get_taxonomy.py -s ftp,entrez -e $email -f COI.fasta -v | extract_taxonomy.py > COI_taxonomy.txt
```

*Note:* Here, we use the FASTA file directly as input

*Note 2:* Since we don't need the intermediate taxonomy format, so we just chain the `get_taxonomy.py` and `extract_taxonomy.py`.

*Progress:* With the `-v` argument supplied to `get_taxonomy.py`, the progress is printed, which is useful with such an amount of data. It may not proceed linearly because the accessions, for which taxonomy is obtained are not evenly distributed in the taxonomy dump files.

## Other (generally useful) scripts

* `tax_consensus.py`: Tries to find reasonable taxonomic lineages for clusters of sequences with taxonomic annotation (see `tax_consensus.py -h`)
* `extract_bam_region.py`: Extracts all aligned sequences from a SAM/BAM file given a BED file with coordinates relative to the SAM/BAM references.

### Specialized scripts

* `itsx2bed.py`: Takes the output of ITSx and creates BED files and a list of recognized domains from it.
* `unique_tax.py`: Searches for duplicate names within a UTAX formatted file that are on the same taxonomic level and makes them unique (see `unique_tax.py -h`)
* `assemble_sam.py`: Assembles OTU/ASV sequences obtained from sequencing the same marker(s) with different but overlapping primer combinations. Takes SAM/BAM files that can be created by USEARCH/VSEARCH `-usearch_global` of raw (errorneous reads) against OTUs/ASVs in an USEARCH-style amplicon pipeline. Each primer combination needs to be annotated as "amplicon group" in the read IDs.
* `primer_analysis.py`: Scans an aligned SAM/BAM file for primer mismatches given primers. Useful to analyse primer mismatches in environmental amplicon data for primers that were *not* used for sequencing (thus occurring within the sequence). But the reads in SAM/BAM may also be obtained from sequencing with different primer combinations, as long as the amplicon primers were trimmed before mapping the reads.
