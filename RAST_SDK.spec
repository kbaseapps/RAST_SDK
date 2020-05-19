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
	
	/* Parameters for the annotate_genome method.

		ncbi_taxon_id - the numeric ID of the NCBI taxon to which this genome belongs. If this
			is included scientific_name is ignored.
		relation_engine_timestamp_ms - the timestamp to send to the Relation Engine when looking
			up taxon information in milliseconds since the epoch.
		scientific_name - the scientific name of the genome. Overridden by ncbi_taxon_id.

		TODO: document remainder of parameters.
	*/
	typedef structure {
		string workspace;
		genome_id input_genome;
		contigset_id input_contigset;
		int genetic_code;
		string domain;
		int ncbi_taxon_id;
		int relation_engine_timestamp_ms;
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
	
	/* Parameters for the annotate_genomes method.

		ncbi_taxon_id - the numeric ID of the NCBI taxon to which this genome belongs. If this
			is included scientific_name is ignored.
		relation_engine_timestamp_ms - the timestamp to send to the Relation Engine when looking
			up taxon information in milliseconds since the epoch.
		scientific_name - the scientific name of the genome. Overridden by ncbi_taxon_id.
		
		TODO: document remainder of parameters.
	*/
	typedef structure {
		string workspace;
		list<GenomeParams> input_genomes;
		int genetic_code;
		string domain;
		int ncbi_taxon_id;
		int relation_engine_timestamp_ms;
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
        */

        typedef structure {
            data_obj_ref object_ref;
            string output_workspace;
            string output_metagenome_name;
            bool create_report;
        } MetagenomeAnnotateParams;

        typedef structure {
            metagenome_ref output_metagenome_ref;
            string output_workspace;
            string report_name;
            string report_ref;
        } MetagenomeAnnotateOutput;

        funcdef annotate_metagenome(MetagenomeAnnotateParams params) 
                returns (MetagenomeAnnotateOutput output) authentication required;

        typedef structure {
            list<data_obj_ref> input_AMAs;
            list<data_obj_ref> input_assemblies;
            string input_text;
            string output_workspace;
            string output_AMASet_name;
            bool create_report;
        } BulkAnnotateMetagenomesParams;

        typedef structure {
            data_obj_ref output_AMASet_ref;
            string output_workspace;
        } BulkMetagenomesAnnotateOutput;

        funcdef annotate_metagenomes(BulkAnnotateMetagenomesParams params)
                returns (BulkMetagenomesAnnotateOutput output) authentication required;


        /*
            Required parameters for rast_genome_assembly:
                object_ref - reference to a Genome or Assembly object,
                output_workspace - output workspace name,
                output_genome_name - output object name

            Optional parameters for rast_genome_assembly:
		ncbi_taxon_id - the numeric ID of the NCBI taxon to which this genome belongs. If this
			        is included scientific_name is ignored.
		relation_engine_timestamp_ms - the timestamp to send to the Relation Engine when looking
			up taxon information in milliseconds since the epoch.
		scientific_name - the scientific name of the genome. Overridden by ncbi_taxon_id.
        */
        typedef structure {
            data_obj_ref object_ref;
            string output_workspace;
            int ncbi_taxon_id;
            int relation_engine_timestamp_ms;
            string scientific_name;
            string output_genome_name;
            bool create_report;
        } RastGenomeAssemblyParams;

        typedef structure {
            genome_id output_genome_ref;
            string output_workspace;
            string report_name;
            string report_ref;
        } RastGenomeAssemblyOutput;

        funcdef rast_genome_assembly(RastGenomeAssemblyParams params)
                returns (RastGenomeAssemblyOutput output) authentication required;


        /*
            For RAST annotating genomes/assemblies
 
            Reference to a set of annotated Genome and/or Assembly objects in the workspace
            @id ws KBaseSearch.GenomeSet
        */
        typedef string genomeSet_ref;

        typedef structure {
            list<data_obj_ref> input_genomes;
            list<data_obj_ref> input_assemblies;
            string input_text;
            string output_workspace;
            int ncbi_taxon_id;
            int relation_engine_timestamp_ms;
            string scientific_name;
            string output_GenomeSet_name;
            bool create_report;
        } BulkRastGenomesAssembliesParams;

        typedef structure {
            genomeSet_ref output_GenomeSet_ref;
            string output_workspace;
        } BulkRastGenomesAssembliesOutput;

        funcdef rast_genomes_assemblies(BulkRastGenomesAssembliesParams params)
                returns (BulkRastGenomesAssembliesOutput output) authentication required;
};
