/*
The SDK version of the KBaase Genome Annotation Service.
This wraps genome_annotation which is based off of the SEED annotations.
*/

module RAST_SDK {
	/*
        A binary boolean
    */
    typedef int bool;
	
	/*
		A string representing a genome id.
	*/
	typedef string genome_id;
	
	/*
		A string representing a ContigSet id.
	*/
	typedef string contigset_id;
	
	/*
		A string representing a workspace name.
	*/
	typedef string workspace_name;
	
	
	typedef structure {
	    string workspace;
	    genome_id input_genome;
	    contigset_id input_contigset;
	    int genetic_code;
	    string domain;
	    string scientific_name;
	    string output_genome;
	    bool call_features_rRNA_SEED;
	    bool call_features_tRNA_trnascan;
	    bool call_selenoproteins;
	    bool call_pyrrolysoproteins;
	    bool call_features_repeat_region_SEED;
	    bool call_features_insertion_sequences;
	    bool call_features_strep_suis_repeat;
	    bool call_features_strep_pneumo_repeat;
	    bool call_features_crispr;
	    bool call_features_CDS_glimmer3;
	    bool call_features_CDS_prodigal;
	    bool call_features_CDS_genemark;
	    bool annotate_proteins_kmer_v2;
	    bool kmer_v1_parameters;
	    bool annotate_proteins_similarity;
	    bool resolve_overlapping_features;
	    bool find_close_neighbors;
	    bool call_features_prophage_phispy;
	    bool retain_old_anno_for_hypotheticals;
	} AnnotateGenomeParams;
	
	typedef structure {
	    workspace_name workspace;
	    string id;
	} AnnotateGenomeResults;
	
	/*
		annotate genome
		params - a param hash that includes the workspace id and options
	*/
	funcdef annotate_genome(UnspecifiedObject params) returns (AnnotateGenomeResults) authentication required;
};
