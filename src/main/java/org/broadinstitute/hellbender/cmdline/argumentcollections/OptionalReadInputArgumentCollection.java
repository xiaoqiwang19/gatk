package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GatkInputURI;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * An argument collection for use with tools that accept zero or more input files containing reads
 * (eg., BAM/SAM/CRAM files).
 */
public final class OptionalReadInputArgumentCollection extends ReadInputArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "BAM/SAM/CRAM file containing reads",
            optional = true,
            common = true)
    private List<GatkInputURI> readFilesNames = new ArrayList<>();

    @Override
    public List<File> getReadFiles() {
        ArrayList<File> ret = new ArrayList<>();
        for (GatkInputURI fn : readFilesNames) {
            ret.add(fn.toPath().toFile());
        }
        return ret;
    }

    @Override
    public List<Path> getReadPaths() {
        ArrayList<Path> ret = new ArrayList<>();
        for (GatkInputURI fn : readFilesNames) {
            ret.add(fn.toPath());
        }
        return ret;
    }

    @Override
    public List<String> getReadFilesNames() {
        return new ArrayList<>(readFilesNames.stream().map(f -> f.getURIString()).collect(Collectors.toList()));
    }
}
