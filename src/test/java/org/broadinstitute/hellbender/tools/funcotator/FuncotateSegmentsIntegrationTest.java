package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class FuncotateSegmentsIntegrationTest extends CommandLineProgramTest {
    // TODO: Move to a funcotator specific test dir
    private static final String TEST_SUB_DIR = toolsTestDir + "/copynumber/utils/";
    // TODO: Copy and create a specific test file for funcotating segments
    private static final String SIMPLE_TEST_FILE = TEST_SUB_DIR + "merge-annotated-regions-by-annotation-legacy-test.seg";
    private static final String REF = hg19MiniReference;
    private static final String DS_PIK3CA_DIR  = largeFileTestDir + "funcotator" + File.separator + "small_ds_pik3ca" + File.separator;

    @Test
    public void testSimple() throws IOException {
        final File outputFile = File.createTempFile("funcotatesegs_simple", ".seg");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(SIMPLE_TEST_FILE);
        arguments.add("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        arguments.add("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.add("hg19");
        arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_PIK3CA_DIR);

        runCommandLine(arguments);

    }
}
