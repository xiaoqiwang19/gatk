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
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;

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


    @Override
    protected boolean isAcceptableFeatureType(final Class<? extends Feature> featureType) {
        return featureType.equals(AnnotatedInterval.class);
    }

    @Override
    public void onTraversalStart() {
        // Initialize a funcotator engine to handle segments

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
