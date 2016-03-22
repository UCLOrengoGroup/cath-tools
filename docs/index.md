# CATH Tools [![Build Status](https://travis-ci.org/UCLOrengoGroup/cath-tools.svg?branch=master)](https://travis-ci.org/UCLOrengoGroup/cath-tools) [![Documentation Status](https://readthedocs.org/projects/cath-tools/badge/?version=latest)](https://readthedocs.org/projects/cath-tools/?badge=latest)

**Recent 64-bit Linux Binaries**: [cath-ssap](https://cath-tools.s3.amazonaws.com/UCLOrengoGroup/cath-tools/125/125.1/release_build/cath-ssap "A Linux cath-ssap binary from a recent commit"), [cath-superpose](https://cath-tools.s3.amazonaws.com/UCLOrengoGroup/cath-tools/125/125.1/release_build/cath-superpose "A Linux cath-superpose binary from a recent commit"), [cath-refine-align](https://cath-tools.s3.amazonaws.com/UCLOrengoGroup/cath-tools/125/125.1/release_build/cath-refine-align "A Linux cath-refine-align binary from a recent commit") (chmod these to be executable)

**Additional Binaries**:
[secmake](http://github.com/UCLOrengoGroup/secmake) (makes the sec files that cath-ssap uses)

**User documentation**: at [Read the Docs](http://cath-tools.readthedocs.org/en/latest/ "The CATH Tools user documentation at Read the Docs")

**Code**: at [GitHub](https://github.com/UCLOrengoGroup/cath-tools "The CATH Tools GitHub respository") 

**Builds**: at [Travis-CI](https://travis-ci.org/UCLOrengoGroup/cath-tools "The CATH Tools Travis-CI builds") 

Cloning the cath-tools GitHub Repo
---

This project includes [bioplib](https://github.com/ACRMGroup/bioplib "Bioplib's GitHub homepage") as a submodule. To ensure the bioplib directory gets populated, clone with:

    git clone --recursive https://github.com/UCLOrengoGroup/cath-tools.git

...or if you've already cloned, then use:

    git submodule update --init --recursive

Quickstart: Running SSAP
========================

The binary for SSAP is called `cath-ssap` and it'll tell you its usage if you run `cath-ssap --help`.

Most external users will want to use the `--protein_source_files PDB_DSSP_SEC` option (though this may become the default in the future).

You will need to tell it where to find the PDB, wolf and sec files. We recommend you use the path options to manage your list of directories to search. This can be done with command line options, but it may be easier to add suitable environment variables to your profile:

    CATH_TOOLS_PDB_PATH  .:/global/data/directories/pdb
    CATH_TOOLS_DSSP_PATH .:/global/data/directories/dssp
    CATH_TOOLS_SEC_PATH  .:/global/data/directories/sec

Extra Bioplib Notes
-------------------

At present, there's an issue that stops us using [bioplib](https://github.com/ACRMGroup/bioplib "Bioplib's GitHub Homepage") &ge; v3.0 so the Git submodule is v2.1.2.

Download the compressed archive from [v2.1.2 on GitHub](https://github.com/ACRMGroup/bioplib/archive/V2.1.2.tar.gz)
via [the project page](http://www.bioinf.org.uk/software/bioplib/), decompress the archive and cd into the new directory.

If you attempt to use &ge; v3.0, be sure to comment out the line in the Makefile that switches on use of libxml2:

    # Include libxml2
    # Note: xml2-config is installed with libxml2
    #       Link to libxml2 with -lxml2
    COPT := $(COPT) -D XML_SUPPORT $(shell xml2-config --cflags)

Otherwise the build against the bioplib libraries will be broken by unresolved dependencies to libxml2.
