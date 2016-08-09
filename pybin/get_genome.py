#!/usr/bin/python

#Usage: python data_api_test.py ws-url shock-url handle-url token ref id

import sys
import json
import re
import doekbase.data_api
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI , GenomeAnnotationClientAPI
from doekbase.data_api.sequence.assembly.api import AssemblyAPI , AssemblyClientAPI

ga = GenomeAnnotationAPI({
	'workspace_service_url' : sys.argv[1],
	'shock_service_url' : sys.argv[2],
	'handle_service_url' : sys.argv[3] 
},token = sys.argv[4],ref = sys.argv[5]);

gto = {
	'id' : sys.argv[6],
	'scientific_name' : "Unknown species",
	'domain' : "Unknown",
	'genetic_code' : 11,
	'dna_size' : 0,
	'num_contigs' : 0,
	'contig_lengths' : [],
	'contig_ids' : [],
	'source' : "KBase",
	'source_id' : sys.argv[6],
#	'md5' : "",
	'taxonomy' : "Unknown",
	'gc_content' : 0.5,
	'complete' : 1,
	'features' : []
};

taxon = {};
success = 0;
try:
	taxon = ga.get_taxon();
	success = 1;
except Exception, e:
	success = 0
	
if success == 1:
	try:
		gto['scientific_name'] = taxon.get_scientific_name()
	except Exception, e:
		success = 0
	try:
		gto['domain'] = taxon.get_domain()
	except Exception, e:
		success = 0
	try:
		gto['genetic_code'] = taxon.get_genetic_code()
	except Exception, e:
		success = 0
	try:
		gto['taxonomy'] = ",".join(taxon.get_scientific_lineage())
	except Exception, e:
		success = 0

assemb = {};
success = 0;
try:
	assemb = ga.get_assembly();
	success = 1;
except Exception, e:
	success = 0
	
if success == 1:
	gto['contigobj'] = {
		'id' : sys.argv[6],
		'name' : sys.argv[6],
		'source' : 'KBase',
		'source_id' : sys.argv[6],
		'md5' : "",
		'type' : "Genome",
		'contigs' : []
	};

	contigdata = assemb.get_contigs();
	for contigid in contigdata.keys():
		newcontig = {
			'id' : contigdata[contigid]['contig_id'],
			'length' : contigdata[contigid]['length'],
			'md5' : contigdata[contigid]['md5'],
			'sequence' : contigdata[contigid]['sequence'],
			'genetic_code' : sys.argv[7],
			'replicon_type' : "linear",
			'replicon_geometry' : "linear",
			'name' : contigdata[contigid]['contig_id'],
			'description' : contigdata[contigid]['description'],
			'complete' : 1
		};
		if contigdata[contigid]['is_circular'] == 1:
			newcontig['replicon_type'] = "circular";
			newcontig['replicon_geometry'] = "circular";
		gto['contigobj']['contigs'].append(newcontig);
	try:
		gto['dna_size'] = taxon.get_dna_size()
	except Exception, e:
		success = 0
	try:
		gto['num_contigs'] = taxon.get_number_contigs()
	except Exception, e:
		success = 0
	try:
		gto['contig_lengths'] = taxon.get_contig_lengths()
	except Exception, e:
		success = 0
	try:
		gto['contig_ids'] = taxon.get_contig_ids()
	except Exception, e:
		success = 0
	try:
		gto['gc_content'] = assemb.get_gc_content()
	except Exception, e:
		success = 0
	try:
		extsource = assemb.get_external_source_info()
		gto['contigobj']['source'] = extsource["external_source"]
		gto['contigobj']['source_id'] = extsource["external_source_id"]
		gto['contigobj']['name'] = extsource["external_source_id"]
		gto['source'] = extsource["external_source"]
		gto['source_id'] = extsource["external_source_id"]
	except Exception, e:
		success = 0
		
features = [];
success = 0;
try:
	features = ga.get_features();
	success = 1
except Exception, e:
	success = 0

prot = ga.get_proteins();

if success == 1:
	for ftrid in features.keys():
		ftrdata = features[ftrid]
		if 'feature_type' in ftrdata.keys():
			newfeature = {'id' : ftrid,'type' : ftrdata['feature_type'],'function' : "Unknown",'location' : []}
			array = ftrid.split("_");
			protid = 'protein_'+array[1];
			if array[0] == 'CDS' and protid in prot.keys():
				newfeature['protein_translation'] = prot[protid]['protein_amino_acid_sequence']
			if 'feature_ontology_terms' in ftrdata.keys():
				newfeature['ontology_terms'] = ftrdata['feature_ontology_terms']
			if 'feature_function' in ftrdata.keys():
				newfeature['function'] = ftrdata['feature_function']
			if 'feature_dna_sequence' in ftrdata.keys():
				newfeature['dna_sequence'] = ftrdata['feature_dna_sequence']
			if 'feature_locations' in ftrdata.keys():
				for loc in ftrdata['feature_locations']:
					newfeature['location'].append([loc['contig_id'],loc['start'],loc['strand'],loc['length']])
			#if 'feature_aliases' in ftrdata.keys():
			#newfeature['protein_translation'] = ftrdata['feature_aliases']
			if 'feature_md5' in ftrdata.keys():
				if len(ftrdata['feature_md5']) > 0:
					newfeature['md5'] = ftrdata['feature_md5']
			if 'feature_dna_sequence_length' in ftrdata.keys():
				newfeature['dna_sequence_length'] = ftrdata['feature_dna_sequence_length']		
			gto['features'].append(newfeature);

#print json.dumps(prot, ensure_ascii=False)

print json.dumps(gto, ensure_ascii=False)
print "SUCCESS"