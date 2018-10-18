package org.broadinstitute.hellbender.tools.funcotator.simpletsvoutput;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.LocatableFuncotationCreator;
import org.broadinstitute.hellbender.tools.funcotator.mafOutput.MafOutputRenderer;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.codehaus.plexus.util.StringUtils;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * This class is very versatile, but as a result, it must do some lazy loading after it receives the first write command.
 *
 * This class assumes that funcotation maps will have the same fields regardless of allele.
 */
public class SimpleTsvOutputRenderer extends OutputRenderer {

    private static final Logger logger = LogManager.getLogger(MafOutputRenderer.class);

    private TableWriter<LinkedHashMap<String, String>> writer;
    private boolean isWriterInitialized;
    /**
     * Output column names to possible field names, which could appear in the funcotations.
     */
    private LinkedHashMap<String, List<String>> columnNameToAliasesMap = new LinkedHashMap<>();


    private LinkedHashMap<String, String> columnNameToFuncotationFieldMap;
    private Set<String> excludedOutputFields;
    private Path outputFilePath;
    private final String referenceVersion;
    private final LinkedHashMap<String, String> unaccountedForDefaultAnnotations;
    private final LinkedHashMap<String, String> unaccountedForOverrideAnnotations;

    public SimpleTsvOutputRenderer(final Path outputFilePath,
                                   final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                                   final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                                   final Set<String> toolHeaderLines,
                                   final String referenceVersion, final Set<String> excludedOutputFields,
                                   final Path configPath) {
        this.excludedOutputFields = excludedOutputFields;
        this.outputFilePath = outputFilePath;
        isWriterInitialized = false;
        this.referenceVersion = referenceVersion;

        try {
            final File tmpResourceFile = Resource.getResourceContentsAsFile(configPath.toString());
            columnNameToAliasesMap = createColumnNameToAliasesMap(tmpResourceFile.toPath());
        } catch (final IOException ioe) {
            throw new GATKException.ShouldNeverReachHereException("Could not read config file: " + configPath,
                    ioe);
        }

        this.unaccountedForDefaultAnnotations = unaccountedForDefaultAnnotations;
        this.unaccountedForOverrideAnnotations = unaccountedForOverrideAnnotations;
    }

    // This method uses the fact that the keys of the input are ordered.
    private void initializeWriter(final LinkedHashMap<String, String> columnNameToFieldNameMap) {
        final TableColumnCollection columns = new TableColumnCollection(columnNameToFieldNameMap.keySet());
        try {
            writer = TableUtils.writer(outputFilePath.toFile(), columns, (map, dataLine) -> {
                map.keySet().forEach(k -> dataLine.append(map.get(k)));
            });
        } catch (final IOException ioe) {
            throw new GATKException("Could not open the simple TSV writer.", ioe);
        }
    }

    @Override
    public void close() {
        try {
            writer.close();
        } catch (final IOException ioe) {
            throw new GATKException("Could not close the simple TSV output writing.", ioe);
        }
    }

    @Override
    public void write(final VariantContext variant, final FuncotationMap txToFuncotationMap) {
        if (txToFuncotationMap.getTranscriptList().size() > 1) {
            logger.warn("More than one transcript found.  This should be able to render (grouped by transcript), but you may need to do further processing.  No user action needed.");
        }

        // Ensure that all transcript-allele combinations have the same fields inside the matching funcotations.
        if (!txToFuncotationMap.doAllTxAlleleCombinationsHaveTheSameFields()) {
            throw new GATKException.ShouldNeverReachHereException("The funcotation map cannot be written by this simple output renderer.  This is almost certainly an issue for the GATK development team.");
        }

        // This class will create a set of funcotations based on the locatable info of a variant context.
        for (final String txId : txToFuncotationMap.getTranscriptList()) {
            for (final Allele allele : txToFuncotationMap.getAlleles(txId)) {
                // TODO: Get rid of null here.  Replace with facility to make a funcotation metadata.
                txToFuncotationMap.add(txId, LocatableFuncotationCreator.create(variant, allele, "SIMPLE_TSV_OUTPUT_RENDERER", null));
            }
        }

        if (!isWriterInitialized) {
            // Developer note:  the "get(0)" is okay, since we know that all transcript-allele combinations have the same fields.
            // Note that columnNameToFuncotationFieldMap will also hold the final say for what columns will get rendered in the end.
            columnNameToFuncotationFieldMap = createColumnNameToFieldNameMap(columnNameToAliasesMap, txToFuncotationMap, txToFuncotationMap.getTranscriptList().get(0), excludedOutputFields);
            initializeWriter(columnNameToFuncotationFieldMap);
        }
        try {
            for (final String txId : txToFuncotationMap.getTranscriptList()) {
                for (final Allele allele : txToFuncotationMap.getAlleles(txId)) {
                    final LinkedHashMap<String, String> columnNameToValueMap = createColumnNameToValueMap(columnNameToFuncotationFieldMap,
                            txToFuncotationMap, txId, allele, unaccountedForDefaultAnnotations,
                            unaccountedForOverrideAnnotations, excludedOutputFields);
                    writer.writeRecord(columnNameToValueMap);
                }
            }
        } catch (final IOException ioe) {
            throw new GATKException("Could not write to the simple TSV writer.", ioe);
        }
    }


    private static LinkedHashMap<String, String> createColumnNameToFieldNameMap(final LinkedHashMap<String, List<String>> columnNameToAliasMap,
                                                                            final FuncotationMap txToFuncotationMap, final String txId, final Set<String> excludedFields) {
        final LinkedHashMap<String, String> result = columnNameToAliasMap.entrySet().stream()
                .map(e -> createAliasList(columnNameToAliasMap, e.getKey()))
                .map(candidates -> Pair.of(candidates.get(0), findFieldNameInFuncotations(candidates, txToFuncotationMap, txId)))
                .collect(Collectors.toMap(p -> p.getLeft(), p -> p.getRight(),
                        (x1, x2) -> {
                            throw new IllegalArgumentException("Should not be able to have duplicate field names."); },
                        LinkedHashMap::new ));

        // Now lets find the columns that are left over and tack those on with their base name.
        final Set<String> leftoverColumns = Sets.difference(txToFuncotationMap.getFieldNames(txId), new HashSet<>(result.values()));
        final List<String> leftoverColumnsAsList = leftoverColumns.stream().sorted().collect(Collectors.toList());
        leftoverColumnsAsList.forEach(c -> result.put(c, c));

        // Remove the excluded columns
        excludedFields.forEach(e -> result.remove(e));

        return result;
    }

    /**    // TODO: Add unaccounted defaults and overrides.
     * // TODO: Test override and defaults
     * //TODO: Update parameters
     * Simple method that will produce a linked field:value hashmap using starting fields from the funcotations in a {@link FuncotationMap}
     *  and a map of funcotation field name
     * @param columnNameToFieldMap alias map of funcotation field names to output alias.  Never {@code null}
     * @param txToFuncotationMap
     * @param txId
     * @param allele
     * @return name value pairs that have been adjusted to the output coluymn names as the field.  Never {@code null}
     */
    private static LinkedHashMap<String, String> createColumnNameToValueMap(final LinkedHashMap<String, String> columnNameToFieldMap,
                                                                      final FuncotationMap txToFuncotationMap, final String txId, final Allele allele,
                                                                            final LinkedHashMap<String, String> unaccountedDefaultAnnotations,
                                                                            final LinkedHashMap<String, String> unaccountedOverrideAnnotations, final Set<String> excludedFields) {

        final LinkedHashMap<String, String> result = new LinkedHashMap<>();
        unaccountedDefaultAnnotations.entrySet().stream()
                .filter(e -> !excludedFields.contains(e.getKey()))
                .forEach(e -> result.put(e.getKey(),e.getValue()));

        // Simply grab the funcotation field using the column name to funcotation field map.  Filter out any values.
        //  Note that this will overwrite the defaults specified in the default annotations.
        columnNameToFieldMap.entrySet().stream()
                .filter(e -> !excludedFields.contains(e.getKey()))
                .map(colEntry -> Pair.of(colEntry.getKey(), txToFuncotationMap.getFieldValue(txId, colEntry.getValue(), allele)))
                .forEach(p -> result.put(p.getLeft(), p.getRight()));

        // For override annotations, just put the values in regardless of what is there.
        unaccountedOverrideAnnotations.entrySet().stream()
                .filter(e -> !excludedFields.contains(e.getKey()))
                .forEach(e -> result.put(e.getKey(), e.getValue()));

        return result;
    }

    //TODO: Docs (candidate is ordered by priority)  This method also assumes that funcotation fields are not allowed to have null values.
    // TODO: Mention that this assumes that the funcotation map has the same fields regardless of allele.
    private static String findFieldNameInFuncotations(final List<String> candidateFieldNames, final FuncotationMap txToFuncotationMap, final String txId) {
        final Set<String> funcotationFieldNames = txToFuncotationMap.getFieldNames(txId);
        return candidateFieldNames.stream()
                .filter(funcotationFieldNames::contains)
                .findFirst().orElse("");
    }

    // Make sure that the column name itself is at the front of the alias list.
    private static List<String> createAliasList(final LinkedHashMap<String, List<String>> columnNameToAliasMap, final String columnName) {
        final List<String> result = new ArrayList<>();
        result.add(columnName);
        result.addAll(columnNameToAliasMap.get(columnName));
        return result;
    }

    //TODO: Docs
    private static LinkedHashMap<String, List<String>> createColumnNameToAliasesMap(final Path resourceFile){
        final Properties configFileContents = new Properties();
        try ( final InputStream inputStream = Files.newInputStream(resourceFile, StandardOpenOption.READ) ) {
            configFileContents.load(inputStream);
        } catch (final Exception ex) {
            throw new UserException.BadInput("Unable to read from XSV config file: " + resourceFile.toUri().toString(), ex);
        }

        final Set<String> keys = configFileContents.stringPropertyNames();
        return keys.stream().collect(Collectors.toMap(Function.identity(), k -> Arrays.asList(StringUtils.split(configFileContents.getProperty(k), ",")),
                (x1, x2) -> {
                    throw new IllegalArgumentException("Should not be able to have duplicate field names."); },
                LinkedHashMap::new ));
    }

}
