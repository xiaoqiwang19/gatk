package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Arrays;
import java.util.List;

/**
 * TODO: Docs
 * TODO: Document that we use a symbolic allele for the copy neutral situation
 */
public class AnnotatedIntervalToSegmentVariantContextConverter {
    private AnnotatedIntervalToSegmentVariantContextConverter() {}

    // TODO: Get rid of magic constants.
    final static List<String> callAnnotationNames = Arrays.asList(
            "CALL","Segment_Call","Call"
    );

    public final static String COPY_NEUTRAL_ALLELE_STRING = "<COPY_NEUTRAL>";
    public final static Allele COPY_NEUTRAL_ALLELE = Allele.create(COPY_NEUTRAL_ALLELE_STRING);

    /**
     * TODO: Docs
     * TODO: Test
     * @param segment
     * @param referenceContext
     * @return
     */
    public static VariantContext convert(final AnnotatedInterval segment, final ReferenceContext referenceContext) {
        final SimpleInterval startPointInterval = new SimpleInterval(segment.getContig(), segment.getStart(), segment.getStart());
        final SimpleInterval endPointInterval = new SimpleInterval(segment.getContig(), segment.getEnd(), segment.getEnd());
        final Allele refAlleleAtStart = Allele.create(referenceContext.getBases(startPointInterval), true);
        final Allele refAlleleAtEnd = Allele.create(referenceContext.getBases(endPointInterval), true);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        final CalledCopyRatioSegment.Call call = retrieveCall(segment);

        return variantContextBuilder
                .chr(segment.getContig())
                .start(segment.getStart())
                .stop(segment.getEnd())
                .alleles(Arrays.asList(
                        Allele.create(refAlleleAtStart, true),
                        convertActualSegCallToAllele(call)
                ))
                .attribute(VCFConstants.END_KEY, segment.getEnd())
                .make();
    }

    private static CalledCopyRatioSegment.Call retrieveCall(final AnnotatedInterval segment) {
        for (final String callAnnotationName : callAnnotationNames) {
            if (segment.hasAnnotation(callAnnotationName)) {
                final CalledCopyRatioSegment.Call call = Arrays.stream(CalledCopyRatioSegment.Call.values())
                        .filter(c -> c.getOutputString().equals(segment.getAnnotationValue(callAnnotationName))).findFirst().orElse(null);
                return call;
            }
        }

        return CalledCopyRatioSegment.Call.NEUTRAL;
    }

    // TODO: Document the assumption that copy neutral is being handled in the genotyping
    private static Allele convertActualSegCallToAllele(final CalledCopyRatioSegment.Call call) {
        switch (call) {
            case DELETION:
                return Allele.create(SimpleSVType.createBracketedSymbAlleleString("DEL"), false);
            case AMPLIFICATION:
                return Allele.create(SimpleSVType.createBracketedSymbAlleleString("INS"), false);
            case NEUTRAL:
                return COPY_NEUTRAL_ALLELE;
            default:
                throw new GATKException.ShouldNeverReachHereException(call.getOutputString() + " is not represented in conversion to variant context.");
        }
    }
}
