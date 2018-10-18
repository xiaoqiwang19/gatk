package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;

import java.util.Arrays;

/**
 * Implements fields for use in known locatables.
 */
public class LocatableFuncotationCreator {

    // TODO: Remove magic constants
    final static String CONTIG_FIELD_NAME = "CONTIG";
    final static String START_FIELD_NAME = "START";
    final static String END_FIELD_NAME = "END";

    // TODO: Test
    // TODO: Docs (only creates the locatable fields)
    public static Funcotation create(final Locatable locatable, final Allele altAllele, final String dataSourceName, final FuncotationMetadata metadata ) {

        return TableFuncotation.create(Arrays.asList(CONTIG_FIELD_NAME, START_FIELD_NAME, END_FIELD_NAME),
                Arrays.asList(locatable.getContig(), String.valueOf(locatable.getStart()), String.valueOf(locatable.getEnd())),
                altAllele, dataSourceName, metadata);
    }

    // TODO: Add segment variant context converter?
}
