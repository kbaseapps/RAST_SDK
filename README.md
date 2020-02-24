[![Build Status](https://travis-ci.org/kbaseapps/RAST_SDK.svg?branch=master)](https://travis-ci.org/kbaseapps/RAST_SDK)

# RAST_SDK
---

This module wraps RAST using the KBase SDK.  Most of the functional aspects can be found in the genome_annotation repo.
https://github.com/kbase/genome_annotation.  This module includes reference data that is pulled from RAST FTP servers.


## Running tests

Note that tests must be run in a KBase environment with an active search instance, as
GenomeFileUtils contacts search for some operations.

At least 100G of space is required to run tests.

