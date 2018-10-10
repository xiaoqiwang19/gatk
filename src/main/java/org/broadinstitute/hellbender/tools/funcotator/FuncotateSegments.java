package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

// TODO: What about using a FeatureWalker?  How much work is that?
// TODO: Fill in the docs and do not forget Oncotator citation.
@CommandLineProgramProperties(
        summary = " (similar functionality to Oncotator).",
        oneLineSummary = "Functional Annotator",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class FuncotateSegments extends FeatureWalker<AnnotatedInterval> {
    private static final Logger logger = LogManager.getLogger(FuncotateSegments.class);
    @Argument(
            doc = "Input segment file (tab-separated values).",
            fullName = CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME
    )
    private File segmentFile;

    @ArgumentCollection
    private final FuncotatorSegmentArgumentCollection funcotatorArgs = new FuncotatorSegmentArgumentCollection();


    private FuncotatorEngine funcotatorEngine;

    @Override
    protected boolean isAcceptableFeatureType(final Class<? extends Feature> featureType) {
        return featureType.equals(AnnotatedInterval.class);
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {

        Utils.validateArg(funcotatorArgs.transcriptSelectionMode != TranscriptSelectionMode.ALL, "Cannot funcotate segments with the ALL transcript selection mode.  Please select another mode.");

        // Get our overrides for annotations:
        final LinkedHashMap<String, String> annotationDefaultsMap = FuncotatorEngine.splitAnnotationArgsIntoMap(funcotatorArgs.annotationDefaults);
        final LinkedHashMap<String, String> annotationOverridesMap = FuncotatorEngine.splitAnnotationArgsIntoMap(funcotatorArgs.annotationOverrides);

        // Next set up our transcript list:
        final Set<String> finalUserTranscriptIdSet = FuncotatorEngine.processTranscriptList(funcotatorArgs.userTranscriptIdSet);

        // Initialize the datasources (and make sure to filter to handle only segment-enabled funcotation factories.
        // Initialize all of our data sources:
        // Sort data sources to make them process in the same order each time:
        funcotatorArgs.dataSourceDirectories.sort(Comparator.naturalOrder());
        final Map<Path, Properties> configData = DataSourceUtils.getAndValidateDataSourcesFromPaths(funcotatorArgs.referenceVersion, funcotatorArgs.dataSourceDirectories);

        // Create the data sources from the input:
        // This will also create and register the FeatureInputs (created by the Data Sources)
        // with the GATK Engine, so we do not have to plumb them in after the fact.
        //  Only take datasources that support the annotation of segments.
        final List<DataSourceFuncotationFactory> dataSourceFuncotationFactories = DataSourceUtils.createDataSourceFuncotationFactoriesForDataSources(
                configData,
                annotationOverridesMap,
                funcotatorArgs.transcriptSelectionMode,
                finalUserTranscriptIdSet,
                this,
                funcotatorArgs.lookaheadFeatureCachingInBp)
                .stream()
                .filter(ff -> ff.isSupportingSegmentFuncotation())
                .collect(Collectors.toList());

        // Initialize a funcotator engine to handle segments.
        // TODO: Create some metadata for typical columns.
        funcotatorEngine = new FuncotatorEngine(funcotatorArgs,
                getBestAvailableSequenceDictionary(),
                VcfFuncotationMetadata.create(
                        new ArrayList<>()
                ),
                dataSourceFuncotationFactories
        );

        // Create a composite output renderer.
    }

    @Override
    public void apply(final AnnotatedInterval feature, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // Convert feature into a VariantContext

        // funcotate

        // write the variant context
    }

    @Override
    public File getDrivingFeatureFile() {
        return segmentFile;
    }
}
