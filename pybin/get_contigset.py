#!/usr/bin/python

#Usage: python data_api_test.py ws-url shock-url handle-url token ref id

import sys
import json
import doekbase.data_api
from doekbase.data_api.sequence.assembly.api import AssemblyAPI , AssemblyClientAPI

assemb = AssemblyAPI({
	'workspace_service_url' : sys.argv[1],
	'shock_service_url' : sys.argv[2],
	'handle_service_url' : sys.argv[3] 
},token = sys.argv[4],ref = sys.argv[5]);

contigset = {
	'id' : sys.argv[6],
	'name' : assemb.get_external_source_info()["external_source_id"],
	'md5' : "",
	'source_id' : assemb.get_external_source_info()["external_source_id"],
	'source' : assemb.get_external_source_info()["external_source"],
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
		'genetic_code' : 11,
		'replicon_type' : "linear",
		'replicon_geometry' : "linear",
		'name' : contigdata[contigid]['contig_id'],
		'description' : contigdata[contigid]['description'],
		'complete' : 1
	};
	if contigdata[contigid]['is_circular'] == 1:
		newcontig['replicon_type'] = "circular";
		newcontig['replicon_geometry'] = "circular";
	contigset['contigs'].append(newcontig);

print json.dumps(contigset, ensure_ascii=False)
print "SUCCESS"