package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.SimpleCount;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.RecordCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.io.Serializable;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Model the baseline ploidy per contig for germline samples given their fragment counts.
 * These should be either HDF5 or TSV count files generated by {@link CollectFragmentCounts}.
 * A table specifying priors for the ploidy per contig is also required.
 *
 *  <p>
 *      If multiple samples are input, then the output is a model file (which can be used for subsequently determining
 *      ploidy in individual samples, see below), a table containing the global read-depth parameter for each sample,
 *      and a table containing the baseline ploidy per contig for each sample.
 *  </p>
 *
 *  <p>
 *      If a single sample and a model file are input, then only read-depth and ploidy tables for that sample
 *      are output.
 *  </p>
 *
 * TODO Mehrtash can add documentation here
 *
 * <h3>Examples</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" DetermineGermlineContigPloidy \
 *   --input normal_1.counts.hdf5 \
 *   --input normal_2.counts.hdf5 \
 *   ... \
 *   --output output_dir \
 *   --outputPrefix normal_cohort
 * </pre>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" DetermineGermlineContigPloidy \
 *   --input normal_1.counts.hdf5 \
 *   --model normal_cohort.ploidyModel.tsv \
 *   ... \
 *   --output output_dir \
 *   --outputPrefix normal_1
 * </pre>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Model the baseline ploidy per contig for germline samples given their fragment counts.",
        oneLineSummary = "Model the baseline ploidy per contig for germline samples given their fragment counts.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class DetermineGermlineContigPloidy extends CommandLineProgram {
    private enum Mode {
        COHORT, CASE
    }

    private static final String DETERMINE_PLOIDY_AND_DEPTH_PYTHON_SCRIPT_PATH =
            "src/main/python/org/broadinstitute/hellbender/tools/gcnv/bin/determine_ploidy_and_depth.py";

    public static final String CONTIG_PLOIDY_PRIORS_FILE_LONG_NAME = "contigPloidyPriors";
    public static final String CONTIG_PLOIDY_PRIORS_FILE_SHORT_NAME = "priors";

    public static final String READ_DEPTH_TABLE_FILE_SUFFIX = ".readDepth.tsv";
    public static final String CONTIG_PLOIDY_TABLE_FILE_SUFFIX = ".contigPloidy.tsv";
    public static final String PLOIDY_MODEL_FILE_SUFFIX = ".ploidyModel.tsv";

    @Argument(
            doc = "Input read-count files containing integer read counts in genomic intervals for all samples.  " +
                    "Intervals must be identical and in the same order for all samples.  " +
                    "If only a single sample is specified, a model file must also be specified.  ",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            minElements = 1
    )
    private List<File> inputReadCountFiles = new ArrayList<>();

    @Argument(
            doc = "Input file specifying contig-ploidy priors.",
            fullName = CONTIG_PLOIDY_PRIORS_FILE_LONG_NAME,
            shortName = CONTIG_PLOIDY_PRIORS_FILE_SHORT_NAME
    )
    private File contigPloidyPriorsFile;

    @Argument(
            doc = "Input ploidy-model file.  If only a single sample is specified, this model will be used.  " +
                    "If multiple samples are specified, a new model will be built and this input will be ignored.",
            fullName = CopyNumberStandardArgument.MODEL_LONG_NAME,
            shortName = CopyNumberStandardArgument.MODEL_SHORT_NAME,
            optional = true
    )
    private File modelFile = null;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME,
            shortName = CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputDir;

    @Advanced
    @ArgumentCollection
    private PloidyDeterminationArgumentCollection ploidyDeterminationArgumentCollection;

    private Mode mode;

    @Override
    protected Object doWork() {
        //TODO different paths for cohort and case modes

        setModeAndValidateArguments();

        //read in count files and output intervals and sample x coverage-per-contig metadata table to temporary files
        final File sampleCoverageMetadataFile = IOUtils.createTempFile("sample-coverage-metadata", ".tsv");
        final List<SimpleInterval> intervals = getIntervalsFromFirstReadCountFile(inputReadCountFiles);
        final File intervalsFile = IOUtils.createTempFile("intervals", ".tsv");
        new SimpleIntervalCollection(intervals).write(intervalsFile);
        writeSampleCoverageMetadata(sampleCoverageMetadataFile, intervals);

        //call python inference code
        final boolean pythonReturnCode = executeDeterminePloidyAndDepthPythonScript(sampleCoverageMetadataFile, intervalsFile);

        if (!pythonReturnCode) {
            throw new UserException("Python return code was non-zero.");
        }

        logger.info("Germline contig ploidy determination complete.");

        return "SUCCESS";
    }

    private void setModeAndValidateArguments() {
        if (modelFile == null) {
            if (inputReadCountFiles.size() > 1) {
                logger.info("Multiple samples provided, running in cohort mode...");
                mode = Mode.COHORT;
            } else {
                throw new UserException("Multiple samples must be provided if a ploidy-model file is not.");
            }
        } else {
            if (inputReadCountFiles.size() > 1) {
                logger.warn("Multiple samples and a ploidy-model file were provided; the latter will be ignored...");
                mode = Mode.COHORT;
            } else {
                logger.info("A single sample and a ploidy-model file were provided, running in case mode...");
                mode = Mode.CASE;
            }
        }

        Utils.validateArg(inputReadCountFiles.size() == new HashSet<>(inputReadCountFiles).size(),
                "List of input read-count files cannot contain duplicates.");
        inputReadCountFiles.forEach(IOUtils::canReadFile);


        Utils.nonNull(outputPrefix);
        if (!new File(outputDir).exists()) {
            throw new UserException(String.format("Output directory %s does not exist.", outputDir));
        }

        ParamUtils.isPositive(ploidyDeterminationArgumentCollection.meanBiasStandardDeviation,
                "Contig-level mean bias standard deviation must be positive.");
        ParamUtils.isPositive(ploidyDeterminationArgumentCollection.mappingErrorRate,
                "Typical mapping error rate must be positive.");
        ParamUtils.isPositive(ploidyDeterminationArgumentCollection.globalPsiScale,
                "Global contig-level unexplained variance scale must be positive.");
        ParamUtils.isPositive(ploidyDeterminationArgumentCollection.samplePsiScale,
                "Sample-specific contig-level unexplained variance scale must be positive.");

    }

    private List<SimpleInterval> getIntervalsFromFirstReadCountFile(final List<File> inputReadCountFiles) {
        final File firstReadCountFile = inputReadCountFiles.get(0);
        logger.info(String.format("Retrieving intervals from first read-count file (%s)...", firstReadCountFile));
        final SimpleCountCollection readCounts = SimpleCountCollection.read(firstReadCountFile);
        return readCounts.getIntervals();
    }

    private void writeSampleCoverageMetadata(final File sampleCoverageMetadataFile,
                                             final List<SimpleInterval> intervals) {
        logger.info("Validating and aggregating metadata from input read-count files...");
        final int numSamples = inputReadCountFiles.size();
        final List<CoverageMetadata> coverageMetadatas = new ArrayList<>(numSamples);
        final List<String> contigs = intervals.stream().map(SimpleInterval::getContig).distinct().collect(Collectors.toList());
        final ListIterator<File> inputReadCountFilesIterator = inputReadCountFiles.listIterator();
        while (inputReadCountFilesIterator.hasNext()) {
            final int sampleIndex = inputReadCountFilesIterator.nextIndex();
            final File inputReadCountFile = inputReadCountFilesIterator.next();
            logger.info(String.format("Aggregating read-count file %s (%d / %d)", inputReadCountFile, sampleIndex + 1, numSamples));
            final SimpleCountCollection readCounts = SimpleCountCollection.read(inputReadCountFile);
            Utils.validateArg(readCounts.getIntervals().equals(intervals),
                    String.format("Intervals for read-count file %s do not match those in other read-count files.", inputReadCountFile));
            //calculate coverage per contig and construct coverage metadata for each sample
            coverageMetadatas.add(new CoverageMetadata(
                    readCounts.getSampleName(),
                    readCounts.getRecords().stream()
                            .collect(Collectors.groupingBy(
                                    SimpleCount::getContig,
                                    LinkedHashMap::new,
                                    Collectors.summingInt(SimpleCount::getCount)))));
        }
        new CoverageMetadataCollection(coverageMetadatas, contigs).write(sampleCoverageMetadataFile);
    }

    private static final class CoverageMetadata {
        private final String sampleName;
        private final LinkedHashMap<String, Integer> coveragePerContig;

        private CoverageMetadata(final String sampleName,
                                 final LinkedHashMap<String, Integer> coveragePerContig) {
            this.sampleName = sampleName;
            this.coveragePerContig = coveragePerContig;
        }
    }

    private static final class CoverageMetadataCollection extends RecordCollection<CoverageMetadata> {
        private static final String SAMPLE_NAME_TABLE_COLUMN = "SAMPLE_NAME";

        private CoverageMetadataCollection(final List<CoverageMetadata> coverageMetadatas,
                                           final List<String> contigs) {
            super(coverageMetadatas,
                    new TableColumnCollection(Stream.of(Collections.singletonList(SAMPLE_NAME_TABLE_COLUMN), contigs)),
                    dataLine -> new CoverageMetadata(
                            dataLine.get(SAMPLE_NAME_TABLE_COLUMN),
                            contigs.stream().collect(Collectors.toMap(
                                    Function.identity(),
                                    dataLine::getInt,
                                    (u, v) -> {throw new GATKException.ShouldNeverReachHereException("Cannot have duplicate contigs.");},   //contigs should already be distinct
                                    LinkedHashMap::new))),
                    (coverageMetadata, dataLine) -> {
                            dataLine.append(coverageMetadata.sampleName);
                            coverageMetadata.coveragePerContig.values().forEach(dataLine::append);});
        }
    }

    private boolean executeDeterminePloidyAndDepthPythonScript(final File sampleCoverageMetadataFile,
                                                               final File intervalsFile) {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final String outputDirArg = Utils.nonEmpty(outputDir).endsWith(File.separator) ? outputDir : outputDir + File.separator;    //add trailing slash if necessary
        return executor.executeScript(
                new Resource(DETERMINE_PLOIDY_AND_DEPTH_PYTHON_SCRIPT_PATH, null),
                null,
                Arrays.asList(
                        "--interval_list " + intervalsFile.getAbsolutePath(),
                        "--sample_coverage_metadata " + sampleCoverageMetadataFile.getAbsolutePath(),
                        "--contig_ploidy_prior_table " + contigPloidyPriorsFile.getAbsolutePath(),
                        "--output_contig_ploidy " + outputDirArg + outputPrefix + CONTIG_PLOIDY_TABLE_FILE_SUFFIX,
                        "--output_read_depth " + outputDirArg + outputPrefix + READ_DEPTH_TABLE_FILE_SUFFIX,
                        "--output_ploidy_model_path " + outputDirArg + outputPrefix + PLOIDY_MODEL_FILE_SUFFIX,
                        ploidyDeterminationArgumentCollection.generatePythonArgumentString()));
    }

    private static final class PloidyDeterminationArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        @Argument(
                doc = "Contig-level mean bias standard deviation.",
                fullName = "meanBiasStandardDeviation",
                minValue = 0.,
                optional = true
        )
        private double meanBiasStandardDeviation = 0.01;

        @Argument(
                doc = "Typical mapping error rate.",
                fullName = "mappingErrorRate",
                minValue = 0.,
                optional = true
        )
        private double mappingErrorRate = 0.01;

        @Argument(
                doc = "Global contig-level unexplained variance scale.",
                fullName = "globalPsiScale",
                minValue = 0.,
                optional = true
        )
        private double globalPsiScale = 0.001;

        @Argument(
                doc = "Sample-specific contig-level unexplained variance scale.",
                fullName = "samplePsiScale",
                minValue = 0.,
                optional = true
        )
        private double samplePsiScale = 0.0001;

        private String generatePythonArgumentString() {
            return String.format(
                    "--mean_bias_sd %f --mapping_error_rate %f --psi_j_scale %f --psi_s_scale %f",
                    meanBiasStandardDeviation, mappingErrorRate, globalPsiScale, samplePsiScale);
        }
    }
}
