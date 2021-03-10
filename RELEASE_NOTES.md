### OVERVIEW
This module wraps the RAST annotation pipeline for KBase.

### Version 1.9.2
__Changes__
Main changes are made to accomodate saving the ontology terms with the annotation_ontology_api service. More code improvements have been made to handle newly discovered errors.

Please note: The two new apps `Annotate Genome/Assembly` and `Bulk Annotate Genomes/Assemblies` have been used by users in both release and beta versions and have proved to be working as required. App `Annotate Genome/Assembly` performs tasks that the current apps `Annotate Microbial Genome` and `Annotate Microbial Assembly` do; and app `Bulk Annotate Genomes/Assemblies` does what apps `Annotate Multiple Microbial Genomes` and `Annotate Multiple Microbial Assemblies` do. The only difference is that the new apps set a fixed set of annotation workflow stages instead of allowing the user to freely select which ones among the 16 stages provided.

Therefore, in about six(6) months, we are going to deprecate the following four(4) RAST_SDK apps:
1) Annotate Microbial Genome (replaced by Annotate Genome/Assembly)
2) Annotate Microbial Assembly (replaced by Annotate Genome/Assembly)
3) Annotate Microbial Multiple Genomes (replaced by Bulk Annotate Genomes/Assemblies)
3) Annotate Microbial Multiple Assemblies (replaced by Bulk Annotate Genomes/Assemblies)


### Version 1.9.1
__Changes__
Made code improvements according to review comments on PR #81 by AJ.  Also added error handling for bulk rasting to catch exception from gfu.save_one_genome call. Included the two apps annotate_genome_assembly and bulk_annotate_genomes_assemblies. Updated documentations for apps Annotate Metagenome Assembly and Re-annotate Metagenomes and Bulk Annotate Genomes/Assemblies

### Version 0.1.8
__Changes__
Merged all development starting from AMA in annotate_metaG_rebased branch into master; then merged changes by AJ on test file/module renaming and cpanfile addition.
Upped the version number for potential release.

### Version 0.1.7
__Changes__
Found the cause and fixed it for the issue 'RAST annotation of MG assembly (meaning make fresh gene call) fails.'
Add new app `Annotate Genome/Assembly` that combines two previous apps to handle annotation of both genome and assembly; app `Bulk Annotate Genomes/Assemblies` for multiple genomes/assemblies annotation

### Version 0.1.6
__Changes__
annotate_genome(s) now takes a taxon ID and timestamp to allow exact specification of the taxon ID
rather than submitting a scientific name. The scientific name is retrieved from the Relation
Engine using the timestamp.

### Version 0.1.5
__Changes__
Prepended the scientific name display in the UI with NCBI.

### Version 0.1.4
__Changes__
Updated the scientific name selector for annotating contigsets to query the RE API for taxon names.

### Version 0.1.3
__Changes__
Unknown

### Version 0.1.2
__Changes__
Unknown

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
