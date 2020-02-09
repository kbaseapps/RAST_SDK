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
	funcdef annotate_genome(AnnotateGenomeParams params) returns (AnnotateGenomeResults) authentication required;
	
	typedef structure {
		contigset_id input_contigset;
		genome_id input_genome;
		genome_id output_genome;
		int genetic_code;
		string domain;
	    string scientific_name;
	} GenomeParams;
	
	typedef structure {
	    string workspace;
	    list<GenomeParams> input_genomes;
		int genetic_code;
		string domain;
	    string scientific_name;
		string genome_text;
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
	    bool call_features_prophage_phispy;
	    bool retain_old_anno_for_hypotheticals;
	} AnnotateGenomesParams;
	
	typedef structure {
	    workspace_name workspace;
	    string report_name;
        string report_ref;
	} AnnotateGenomesResults;
	
	/*
		annotate genomes
		params - a param hash that includes the workspace id and options
	*/
	funcdef annotate_genomes(AnnotateGenomesParams params) returns (AnnotateGenomesResults) authentication required;
	
	typedef structure {
		list<string> proteins;
	} AnnotateProteinParams;
	
	typedef structure {
	    list<list<string> > functions;
	} AnnotateProteinResults;
	
	/*
		annotate proteins - returns a list of the RAST annotations for the input protein sequences
	*/
	funcdef annotate_proteins(AnnotateProteinParams params) returns (AnnotateProteinResults);


        /*
            For RAST annotating metagenomes (borrowed and simplied from ProkkaAnnotation moduel)
        /*
            Reference to an Assembly or Genome object in the workspace
            @id ws KBaseGenomeAnnotations.Assembly
            @id ws KBaseGenomes.Genome
            @id ws KBaseMetagenomes.AnnotatedMetagenomeAssembly
        */
        typedef string data_obj_ref;

        /*
            Reference to a Annotated Metagenome Assembly object in the workspace
            @id ws KBaseMetagenomes.AnnotatedMetagenomeAssembly
        */
        typedef string metagenome_ref;

        /*
            Required parameters:
                object_ref - reference to Assembly or Genome object,
                output_workspace - output workspace name,
                output_metagenome_name - output object name,
            Optional parameters:
                run_prodigal - default to 0, set to 1 if expect to run Prodigal
        */

        typedef structure {
            data_obj_ref object_ref;
            string output_workspace;
            string output_metagenome_name;
            bool run_prodigal;
        } MetagenomeAnnotateParams;

        typedef structure {
            metagenome_ref output_metagenome_ref;
            string output_workspace;
            string report_name;
            string report_ref;
        } MetagenomeAnnotateOutput;

        funcdef annotate_metagenome(MetagenomeAnnotateParams params) 
                returns (MetagenomeAnnotateOutput output) authentication required;
};
