package org.broadinstitute.hellbender.utils.variant.writers;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;

/**
 * Genome-wide VCF writer for somatic (Mutect2) output
 * Merges reference blocks based on TLOD
 */
final public class SomaticGVCFWriter extends GVCFWriter {

    /**
     * Create a new GVCF writer
     *
     * Should be a non-empty list of boundaries.  For example, suppose this variable is
     *
     * [A, B, C]
     *
     * We would partition our hom-ref sites into the following bands:
     *
     * X < A
     * A <= X < B
     * B <= X < C
     * X >= C
     *
     * @param underlyingWriter the ultimate destination of the GVCF records
     * @param lodPartitions     a list of GQ partitions, this list must be non-empty and every element must be larger than previous element
     */
    public SomaticGVCFWriter(final VariantContextWriter underlyingWriter, final List<Integer> lodPartitions) {
        super(underlyingWriter, lodPartitions, 2);
    }

    /**
     * Create {@link HomRefBlock}s which will collectively accept variants of any genotype quality
     *
     * Each individual block covers a band of genotype qualities with the splits between bands occurring at values in {@code gqPartitions}.
     * There will be {@code gqPartitions.size() +1} bands produced covering the entire possible range of genotype qualities from 0 to {@link VCFConstants#MAX_GENOTYPE_QUAL}.
     *
     * @param gqPartitions proposed GQ partitions
     * @return a list of HomRefBlocks accepting bands of genotypes qualities split at the points specified in gqPartitions
     */
    @Override
    @VisibleForTesting
    public RangeMap<Integer,Range<Integer>> parsePartitions(final List<Integer> gqPartitions) {
        Utils.nonEmpty(gqPartitions);
        Utils.containsNoNull(gqPartitions, "The list of TLOD partitions contains a null integer");
        final RangeMap<Integer, Range<Integer>> result = TreeRangeMap.create();
        int lastThreshold = 0;
        for (final Integer value : gqPartitions) {
            if (value < lastThreshold) {
                throw new IllegalArgumentException(String.format("The list of TLOD partitions is out of order. Previous value is %d but the next is %d.", lastThreshold, value));
            } else if (value == lastThreshold) {
                throw new IllegalArgumentException(String.format("The value %d appears more than once in the list of TLOD partitions.", value));
            }

            result.put(Range.closedOpen(lastThreshold, value), Range.closedOpen(lastThreshold, value));
            lastThreshold = value;
        }

        result.put(Range.closedOpen(lastThreshold, Integer.MAX_VALUE), Range.closedOpen(lastThreshold, Integer.MAX_VALUE));

        return result;
    }

    @Override
    boolean genotypeCanBeMergedInCurrentBlock(final Genotype g) {
        final TLODBlock currentTLODBlock = (TLODBlock)currentBlock;
        return currentTLODBlock != null
                && currentTLODBlock.withinBounds(Double.parseDouble(g.getExtendedAttribute(GATKVCFConstants.TUMOR_LOD_KEY).toString()));
    }

    /**
     * Helper function to create a new HomRefBlock from a variant context and current genotype
     *
     * @param vc the VariantContext at the site where want to start the band
     * @param g  the genotype of the sample from vc that should be used to initialize the block
     * @return a newly allocated and initialized block containing g already
     */
    @Override
    GVCFBlock createNewBlock(final VariantContext vc, final Genotype g) {
        // figure out the GQ limits to use based on the GQ of g
        final double lod = Double.parseDouble(g.getExtendedAttribute(GATKVCFConstants.TUMOR_LOD_KEY).toString());
        final Range<Integer> partition = gqPartitions.get((int)Math.round(lod));

        if( partition == null) {
            throw new GATKException("GQ " + g + " from " + vc + " didn't fit into any partition");
        }

        // create the block, add g to it, and return it for use
        final TLODBlock block = new TLODBlock(vc, partition.lowerEndpoint(), partition.upperEndpoint());
        block.add(vc.getStart(), g);
        return block;
    }

}