package org.broadinstitute.hellbender.tools.funcotator.simpletsvoutput;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.funcotator.AnnotatedIntervalToSegmentVariantContextConverter;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;

public class SimpleTsvOutputRendererUnitTest extends GATKBaseTest {
    final public String SEG_CONFIG_RESOURCE = "org/broadinstitute/hellbender/tools/funcotator/simple_oncotator_seg_file.config";

    @DataProvider
    public Object[][] provideForSimpleSegFileWriting() {
        return new Object[][] {
                {FuncotatorTestUtils.createSimpleVariantContext(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        "3", 2100000, 3200000, "T",
                        AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE.getDisplayString())
                }
        };
    }

    @Test(dataProvider = "provideForSimpleSegFileWriting")
    public void testSimpleSegFileWriting(final VariantContext segVC) throws IOException {
        final File outputFile = File.createTempFile("simpleSegFileWriting", ".seg");

        final SimpleTsvOutputRenderer renderer = new SimpleTsvOutputRenderer(outputFile.toPath(),
                new LinkedHashMap<>(),
                new LinkedHashMap<>(),
                new HashSet<>(),
                "hg19", new HashSet<>(), Paths.get(SEG_CONFIG_RESOURCE));

        final FuncotationMap funcotationMap = FuncotationMap.createNoTranscriptInfo(
                Arrays.asList(TableFuncotation.create(Collections.singletonList("Gencode_19_genes"), Collections.singletonList("GENE1,GENE2"),
                        segVC.getAlternateAllele(0), "TEST", createDummySegmentFuncotationMetadata())));
        renderer.write(segVC, funcotationMap);
        renderer.close();

        final AnnotatedIntervalCollection collection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertEquals(collection.getRecords().size(), 1);
        Assert.assertEquals(collection.getRecords().get(0).getAnnotationValue("genes"), "GENE1,GENE2");
        Assert.assertEquals(collection.getRecords().get(0).getContig(), segVC.getContig());
        Assert.assertEquals(collection.getRecords().get(0).getStart(), segVC.getStart());
        Assert.assertEquals(collection.getRecords().get(0).getEnd(), segVC.getEnd());
    }

    private static FuncotationMetadata createDummySegmentFuncotationMetadata() {
        return VcfFuncotationMetadata.create(
                Collections.singletonList(new VCFInfoHeaderLine("Gencode_19_genes",
                        1, VCFHeaderLineType.String, "The genes overlapping the segment."))
        );
    }
}
