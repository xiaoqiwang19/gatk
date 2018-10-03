package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.GenotypeBuilder;

public abstract class ReferenceConfidenceResult {
    public int refDepth = 0;
    public int nonRefDepth = 0;

    /**
     * @return Get the DP (sum of AD values)
     */
    int getDP() {
        return refDepth + nonRefDepth;
    }

    /**
     * Return the AD fields. Returns a newly allocated array every time.
     */
    int[] getAD() {
        return new int[]{refDepth, nonRefDepth};
    }
}
