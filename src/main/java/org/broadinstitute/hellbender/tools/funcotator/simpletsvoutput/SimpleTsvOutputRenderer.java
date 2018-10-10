package org.broadinstitute.hellbender.tools.funcotator.simpletsvoutput;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;

import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

public class SimpleTsvOutputRenderer extends OutputRenderer {

    public SimpleTsvOutputRenderer(final Path outputFilePath,
                                   final List<DataSourceFuncotationFactory> dataSources,
                                   final VCFHeader inputFileHeader,
                                   final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                                   final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                                   final Set<String> toolHeaderLines,
                                   final String referenceVersion, final Set<String> excludedOutputFields) {

    }

    @Override
    public void close() {

    }

    @Override
    public void write(final VariantContext variant, final FuncotationMap txToFuncotationMap) {

    }
}
