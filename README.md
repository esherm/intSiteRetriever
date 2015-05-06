# Get integration sites from Database
For a given sample name gets integration site and matched random controls(MRCs).

# Interface

Functions:

* getUniqueSites
* getMRCs
* setNameExists
* getReferenceGenome

## Generation of random positions on genome

* get_reference_genome finds genome object for UCSC name
* get_random_positions generate gender-specific positions for a genome
* get_N_MRCs generate MRCs for sites that are not in DB


# Testing of the library components

Run in the R console:

```bash
library(testthat)
test_dir(".")
```
