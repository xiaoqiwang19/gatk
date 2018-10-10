package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.barclay.argparser.Argument;

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


    //-----------------------------------------------------
    // Optional args:
    // (See superclass for more)

    @Argument(
            fullName = FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME,
            optional = true,
            doc = "Ignore/drop variants that have been filtered in the input.  These variants will not appear in the output file."
    )
    public boolean removeFilteredVariants = false;




}
