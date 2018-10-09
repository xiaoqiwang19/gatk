package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.File;
import java.io.Serializable;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

abstract class BaseFuncotatorArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output file to which annotated variants should be written.")
    public File outputFile;

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME,
            doc = "The version of the Human Genome reference to use (e.g. hg19, hg38, etc.).  This will correspond to a sub-folder of each data source corresponding to that data source for the given reference."
    )
    public String referenceVersion;

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME,
            doc = "The path to a data source folder for Funcotator.  May be specified more than once to handle multiple data source folders."
    )
    public List<String> dataSourceDirectories;

    // Optional fields ----------------------

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.EXCLUSION_FIELDS_LONG_NAME,
            optional = true,
            doc = "Fields that should not be rendered in the final output.  Only exact name matches will be excluded."
    )
    public Set<String> excludedFields = new HashSet<>();

    @Advanced
    @Hidden
    @Argument(
            fullName = FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION,
            optional = true,
            doc = "(Advanced / DO NOT USE*) If you select this flag, Funcotator will force a conversion of variant contig names from b37 to hg19.  *This option is useful in integration tests (written by devs) only."
    )
    public boolean forceB37ToHg19ContigNameConversion = false;

}
