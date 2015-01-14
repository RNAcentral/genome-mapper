# RNAcentral Genome Mapper

## About

This repository contains the code that is used internally for
mapping [RNAcentral](http://rnacentral.org) sequences to their genomic locations
based on the underlying INSDC accessions.

**Example**:

AL133244.1 (15569:15681) corresponds to human chromosome 2 (32945339:32945451).

Uses [Ensembl Perl API](http://www.ensembl.org/info/docs/api/index.html) and [Bioperl](http://www.bioperl.org/).

## Requirements

* DBD::mysql for connecting to Ensembl
* DBD::Oracle for connecting to RNAcentral

## Installation

```
# clone this repo
git clone https://github.com/RNAcentral/genome-mapper.git

# initialise submodules
cd genome-mapper
git submodule init
git submodule update

# optional: checkout the desired Ensembl branch

# make a copy of the params script and add the connection details
cp scripts/params_template.sh scripts/params.sh
```

The connection parameters for the RNAcentral database are stored in *scripts/params.sh* (excluded from version control).

## Usage

```
source setup.sh
perl genome_mapper.pl -s Homo_sapiens
```
