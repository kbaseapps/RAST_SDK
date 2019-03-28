### OVERVIEW
This module wraps the RAST annotation pipeline for KBase.

### Version 0.1.1
__Changes__
Add the GenomeSet to the objects created when doing Multiple Assemblies or Multiple Genomes

### Version 0.1.0
__Changes__
Add support of genetic code 25 and unknown domains.

### Version 0.0.18
__Changes__
Add support for AssemblySet and a new app for 'Annotate Multiple Microbial Assemblies'

### Version 0.0.17
__Changes__
Add support for GenomeSets to 'Annotate Multiple Microbial Genomes'

### Version 0.0.16
__Changes__
- Improved the report messages to give users more complete and accurate information.
- Updated the Dockerfile to the current version
- Added tests to check all of the advanced options
- Updated the catalog page for multiple genomes. Add the delimiter and reworded the description.

### Version 0.0.15
__Changes__
- added proper citations in PLOS format

### Version 0.0.13
__Changes__
- Handle annotation of new genomes with existing sso ids and update test to cover this.
- Version in the ontology event is now the tool version.
- update source for R dependency

### Version 0.0.12
__Changes__
- Remove find close neighbors option. Users should use genome tree instead.

### Version 0.0.11
__Changes__
- Save genomes though the GenomeFileUtils
- Adapt code for new style genome and add test
- Tests use GenomeAnnotationAPI to load data

### Version 0.0.10
__Changes__
- Switch KMer FigFam service to write to /tmp for HPC support

### Version 0.0.9
__Changes__
- UI now passes resolved references rather than object names
- Use reference chains to access assemblies from genome objects

### Version 0.0.8
__Changes__
- Added phispy dependencies
