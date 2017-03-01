/*
The SDK version of the KBaase Genome Annotation Service.
This wraps genome_annotation which is based off of the SEED annotations.
*/

module RAST_SDK {
	/* A boolean - 0 for false, 1 for true.
                @range (0, 1)
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
	
        /* Input for the annotate_genome function.

        Required parameters:
        workspace - the workspace of the destination (and source if src_workspace is not provided) of the genome/contigset object.
	input_genome - the genome_id for the genome to be RAST-ed;

        Optional parameters:
        src_workspace - the workspace of the source the genome/contigset object, default to workspace.
        dest_workspace - the workspace for the RAST-ed the genome, default to workspace.
        input_contigset - a contigset, defaut to null.
        genetic_code - an int representing the genetic code of the genome;
        domain - the domain of the genome;
	scientific_name - the scientific_name of the genome;
	output_genome - the id for the RAST-ed genome, default to the input_genome;
        The following are a group of bool settings for the RAST processing, default values are set in the implementation
	call_features_rRNA_SEED,
	call_features_tRNA_trnascan,
	call_selenoproteins,
	call_pyrrolysoproteins,
	call_features_repeat_region_SEED,
	call_features_insertion_sequences,
	call_features_strep_suis_repeat,
	call_features_strep_pneumo_repeat,
	call_features_crispr,
	call_features_CDS_glimmer3,
	call_features_CDS_prodigal,
	call_features_CDS_genemark,
	annotate_proteins_kmer_v2,
	kmer_v1_parameters,
	annotate_proteins_similarity,
	resolve_overlapping_features,
	find_close_neighbors,
	call_features_prophage_phispy,
	retain_old_anno_for_hypotheticals
        */
	typedef structure {
	    string workspace;
	    genome_id input_genome;
	    contigset_id input_contigset;
	    int genetic_code;
	    string domain;
	    string scientific_name;
	    genome_id output_genome;
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
	    string report_name;
        string report_ref;
	} AnnotateGenomeResults;
	
	/*
		annotate genome
		params - a param hash that includes the workspace id and options
	*/
	/*
        funcdef annotate_genome(UnspecifiedObject params) returns (AnnotateGenomeResults) authentication required;
	*/
        funcdef annotate_genome(AnnotateGenomeParams params) returns (AnnotateGenomeResults) authentication required;
};
