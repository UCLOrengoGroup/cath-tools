# Summary of changes in cath-tools releases


### [v0.15.3](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.15.3) Improve compatibility of Travis-CI Ubuntu executables with older Ubuntus

2017-09-20 &nbsp; The changes to the Travis-CI build in [v0.14.2](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.14.2) meant that the resulting executables didn't work on Ubuntu &le; 16.04 (out-of-the-box).

This release improves things so that the executables should now work on Ubuntu &ge; 14.10.

See [42](https://github.com/UCLOrengoGroup/cath-tools/issues/42) and [a2d5fa12d1eb0fb0936cf0b1e2ab32ccb48aaaac](https://github.com/UCLOrengoGroup/cath-tools/commit/a2d5fa12d1eb0fb0936cf0b1e2ab32ccb48aaaac) for more info.


### [v0.15.2](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.15.2) Teach cath-resolve-hits to parse hmmscan output

2017-09-14 &nbsp; 


### [v0.15.1](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.15.1) Deprecate cath-resolve-hits' --apply-cath-rules option

2017-09-08 &nbsp; 


### [v0.14.4](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.14.4) cath-map-clusters: a new tool to map between different ways of clustering your data

2017-08-07 &nbsp; 


### [v0.14.2](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.14.2) Add MacOS executables (via the automatic Travis-CI build)

2017-07-28 &nbsp; 


### [v0.14.1](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.14.1) Update cath-resolve-hits' options to allow multiple outputs per run; deprecate some previous options

2017-05-31 &nbsp; 


### [v0.13.1](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.13.1) Make cath-resolve-hits substantially faster on complex data-sets (like titin)

2017-04-28 &nbsp; 


### [v0.13.0](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.13.0) Apply breaking changes to cath-resolve-hits: add headers to output, add JSON format, change default trim, add HMMER evalues to output

2017-03-31 &nbsp; Apply [potentially breaking] changes to cath-resolve-hits:
 * add headers and HMMER evalues to standard output, see [example](https://github.com/UCLOrengoGroup/cath-tools/blob/master/build-test-data/resolve_hits/eg_hmmsearch_out.out)
 * add option `--json-output` to output in JSON, see [example](https://github.com/UCLOrengoGroup/cath-tools/blob/master/build-test-data/resolve_hits/eg_hmmsearch_out.json)
 * change default trim from 50/30 to 30/10


### [v0.12.26](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.26) Fix bug in cath-resolve-hits's handling of short segments

2017-03-10 &nbsp; See [31](https://github.com/UCLOrengoGroup/cath-tools/issues/31) and [cca8fb43e712e710c8eaee4553fa723e137d1c2a](https://github.com/UCLOrengoGroup/cath-tools/commit/cca8fb43e712e710c8eaee4553fa723e137d1c2a) for details.


### [v0.12.25](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.25) Make cath-ssap run directly from PDB files

2017-02-28 &nbsp; 


### [v0.12.24](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.24) Improve cath-resolve-hits HTML output (eg: full result shown in single line; input data limitable)

2017-02-02 &nbsp; 


### [v0.12.23](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.23) Store superposition colour schemes under PyMOL scenes (accessed via the function keys)

2017-01-31 &nbsp; 


### [v0.12.22](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.22) Handle DNA ATOM records in PDBs

2017-01-27 &nbsp; 


### [v0.12.20](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.20) 

2017-01-12 &nbsp; 


### [v0.12.19](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.19) 

2017-01-04 &nbsp; 


### [v0.12.18](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.18) 

2016-12-20 &nbsp; 


### [v0.12.17](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.17) 

2016-11-28 &nbsp; 


### [v0.12.16](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.16) 

2016-11-24 &nbsp; 


### [v0.12.15](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.15) 

2016-11-18 &nbsp; 


### [v0.12.14](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.14) 

2016-11-08 &nbsp; 


### [v0.12.13](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.13) 

2016-11-06 &nbsp; 


### [v0.12.12](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.12) 

2016-11-04 &nbsp; 


### [v0.12.11](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.11) 

2016-11-03 &nbsp; 


### [v0.12.10](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.10) 

2016-11-01 &nbsp; 


### [v0.12.9](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.9) 

2016-10-19 &nbsp; 


### [v0.12.8](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.8) 

2016-10-19 &nbsp; 


### [v0.12.7](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.7) 

2016-09-21 &nbsp; 


### [v0.12.6](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.12.6) 

2016-09-21 &nbsp; 

