# CATH Tools [![Build Status](https://travis-ci.org/UCLOrengoGroup/cath-tools.svg?branch=master)](https://travis-ci.org/UCLOrengoGroup/cath-tools) [![Documentation Status](https://readthedocs.org/projects/cath-tools/badge/?version=latest)](https://readthedocs.org/projects/cath-tools/?badge=latest)


Protein structure comparison tools such as SSAP, as used by the [Orengo Group](https://www.ucl.ac.uk/orengo-group "Orengo Group website") in curating [CATH](http://www.cathdb.info/ "CATH website").



Resources
---------

| | |
|:-- |:-- |
| **Latest Release**          | [**DOWNLOADS**](https://github.com/UCLOrengoGroup/cath-tools/releases/latest "The latest CATH Tools release") (including 64-bit Linux executables, chmod them to be executable) |
| **User documentation**      | at [Read the Docs](http://cath-tools.readthedocs.io/en/latest/ "The CATH Tools user documentation at Read the Docs")                                                            |
| **Related Software**        | [UCLOrengoGroup/secmake](http://github.com/UCLOrengoGroup/secmake) (makes sec files, which are required by cath-ssap)                                                           |
| **Code**                    | at [GitHub](https://github.com/UCLOrengoGroup/cath-tools "The CATH Tools GitHub repository")                                                                                    |
| **Supplementary Materials** | at [GitHub](https://github.com/UCLOrengoGroup/cath-tools-supplementary "The CATH Tools Supplementary GitHub repository")                                                        |
| **Builds**                  | at [Travis-CI](https://travis-ci.org/UCLOrengoGroup/cath-tools "The CATH Tools Travis-CI builds")                                                                               |
| **Acknowledgements**        | the [bioplib](https://github.com/ACRMGroup/bioplib "Bioplib's GitHub Homepage") library, used to execute the superpositions                                                     |



Summary of Tools
----------------

### Primary tools

 * [ `cath-resolve-hits` ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-resolve-hits ) : Collapse a list of domain matches to your query sequence(s) down to the non-overlapping subset (ie domain architecture) that maximises the sum of the hits' scores.
 * [ `cath-ssap`         ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-ssap         ) : Run a SSAP pairwise structural alignment
 * [ `cath-superpose`    ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-superpose    ) : Superpose protein structures using an existing alignment
 * `build-test`                                                                                  : Perform the cath-tools tests (which should all pass, abeit with a few warnings)

### Extras

 * `cath-assign-domains` Use an SVM model on SSAP+PRC data to form a plan for assigning the domains to CATH superfamilies/folds
 * `cath-refine-align`   Iteratively refine an existing alignment by attempting to optimise SSAP score
 * `cath-score-align`    Score an existing alignment using structural data

<!--
| |
|:-- |:-- |
| `cath-assign-domains` | Use an SVM model on SSAP+PRC data to form a plan for assigning the domains to CATH superfamilies/folds |
| `cath-refine-align`   | Iteratively refine an existing alignment by attempting to optimise SSAP score                          |
| `cath-score-align`    | Score an existing alignment using structural data                                                      |
| `check-pdb`           | Check a PDB file for some potential problems                                                           |
-->

Authors
-------

The SSAP algorithm (`cath-ssap`) was devised by Christine A Orengo and William R Taylor.

Please cite: *Protein Structure Alignment*, Taylor and Orengo, Journal of Molecular Biology 208, 1-22, PMID: 2769748. ([PubMed](https://www.ncbi.nlm.nih.gov/pubmed/2769748), [Elsevier](http://www.sciencedirect.com/science/article/pii/0022283689900843))

Since then, many people have contributed to this code, most notably:

  * [Tony E Lewis](https://github.com/tonyelewis)             (2011  - ....)
  * Oliver C Redfern                                          (~2003 - 2011)
  * James E Bray, [Ian Sillitoe](https://github.com/sillitoe) (~2000 - 2003)
  * [Andrew C R Martin](https://github.com/AndrewCRMartin)    (considerable edits around 2001)

Feedback
--------

Please tell us about your cath-tools bugs/suggestions [here](https://github.com/UCLOrengoGroup/cath-tools/issues/new).


