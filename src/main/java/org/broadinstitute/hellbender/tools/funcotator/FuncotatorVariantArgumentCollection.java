package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.barclay.argparser.Argument;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Arguments to be be used by the {@link Funcotator} {@link org.broadinstitute.hellbender.engine.GATKTool},
 * which are specific to {@link Funcotator}.  Use this collection for small mutations (SNP, Indel)
 * Created by jonn on 9/12/18.
 */
public class FuncotatorVariantArgumentCollection extends BaseFuncotatorArgumentCollection {
    private static final long serialVersionUID = 1L;

    //-----------------------------------------------------
    // Required args:
    // (See superclass for more)

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME,
            doc = "The output file format.  Either VCF, MAF, or SEG.  Please note that MAF output for germline use case VCFs is unsupported.  SEG will generate two output files: a simple tsv and a gene list."
    )
    public FuncotatorArgumentDefinitions.OutputFormatType outputFormatType;


    //-----------------------------------------------------
    // Optional args:
    // (See superclass for more)

    @Argument(
            fullName = FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME,
            optional = true,
            doc = "Ignore/drop variants that have been filtered in the input.  These variants will not appear in the output file."
    )
    public boolean removeFilteredVariants = false;

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME,
            optional = true,
            doc = "Method of detailed transcript selection.  This will select the transcript for detailed annotation (CANONICAL, ALL, or BEST_EFFECT)."
    )
    public TranscriptSelectionMode transcriptSelectionMode = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE;

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME,
            optional = true,
            doc = "File to use as a list of transcripts (one transcript ID per line, version numbers are ignored) OR A set of transcript IDs to use for annotation to override selected transcript."
    )
    public Set<String> userTranscriptIdSet = new HashSet<>();

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME,
            optional = true,
            doc = "Annotations to include in all annotated variants if the annotation is not specified in the data sources (in the format <ANNOTATION>:<VALUE>).  This will add the specified annotation to every annotated variant if it is not already present."
    )
    public List<String> annotationDefaults = new ArrayList<>();

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME,
            optional = true,
            doc = "Override values for annotations (in the format <ANNOTATION>:<VALUE>).  Replaces existing annotations of the given name with given values."
    )
    public List<String> annotationOverrides = new ArrayList<>();

    @Argument(
            fullName = FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_NAME,
            optional = true,
            minValue = 0,
            doc = "Number of base-pairs to cache when querying variants."
    )
    public int lookaheadFeatureCachingInBp = FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_DEFAULT_VALUE;


}
