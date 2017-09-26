package org.broadinstitute.hellbender.tools.copynumber.legacy;

import com.google.common.collect.ImmutableSet;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation.AlleleFractionKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation.AlleleFractionSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.TSVLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.segmentation.CRAFSegmentCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Model segmented copy ratio from denoised read counts and minor-allele fraction from allelic counts.",
        oneLineSummary = "Model segmented copy ratio from denoised read counts.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class ModelSegments extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    //filename tags for output
    protected static final String FILTERED_ALLELIC_COUNTS_FILE_SUFFIX = ".filtered.ac.tsv";
    protected static final String SEGMENTS_FILE_SUFFIX = ".seg";
    protected static final String COPY_RATIO_SEGMENTS_FILE_SUFFIX = ".cr" + SEGMENTS_FILE_SUFFIX;
    protected static final String ALLELE_FRACTION_SEGMENTS_FILE_SUFFIX = ".af" + SEGMENTS_FILE_SUFFIX;
    protected static final String CRAF_SEGMENTS_FILE_SUFFIX = ".craf" + SEGMENTS_FILE_SUFFIX;
    protected static final String BEGIN_FIT_FILE_TAG = ".model-begin";
    protected static final String FINAL_FIT_FILE_TAG = ".model-final";
    protected static final String CR_PARAMETER_FILE_SUFFIX = ".cr.param";
    protected static final String AF_PARAMETER_FILE_SUFFIX = ".af.param";

    protected static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME = "maxNumSegmentsPerChromosome";
    protected static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_SHORT_NAME = "maxNumSegsPerChr";

    protected static final String MINIMUM_TOTAL_ALLELE_COUNT_LONG_NAME = "minTotalAlleleCount";
    protected static final String MINIMUM_TOTAL_ALLELE_COUNT_SHORT_NAME = "minAC";

    protected static final String KERNEL_VARIANCE_COPY_RATIO_LONG_NAME = "kernelVarianceCopyRatio";
    protected static final String KERNEL_VARIANCE_COPY_RATIO_SHORT_NAME = "kernVarCR";

    protected static final String KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME = "kernelVarianceAlleleFraction";
    protected static final String KERNEL_VARIANCE_ALLELE_FRACTION_SHORT_NAME = "kernVarAF";

    protected static final String KERNEL_APPROXIMATION_DIMENSION_LONG_NAME = "kernelApproximationDimension";
    protected static final String KERNEL_APPROXIMATION_DIMENSION_SHORT_NAME = "kernApproxDim";

    protected static final String WINDOW_SIZES_LONG_NAME = "windowSizes";
    protected static final String WINDOW_SIZES_SHORT_NAME = "winSizes";

    protected static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_LONG_NAME = "numChangepointsPenaltyFactorCopyRatio";
    protected static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_SHORT_NAME = "numChangeptsPenCR";

    protected static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_LONG_NAME = "numChangepointsPenaltyFactorAlleleFraction";
    protected static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_SHORT_NAME = "numChangeptsPenAF";

    protected static final String NUM_COPY_RATIO_INTERVALS_SMALL_SEGMENT_THRESHOLD_LONG_NAME = "numCopyRatioIntervalsSmallSegmentThreshold";
    protected static final String NUM_COPY_RATIO_INTERVALS_SMALL_SEGMENT_THRESHOLD_SHORT_NAME = "numCRSmallSegTh";

    protected static final String NUM_SAMPLES_COPY_RATIO_LONG_NAME = "numSamplesCopyRatio";
    protected static final String NUM_SAMPLES_COPY_RATIO_SHORT_NAME = "numSampCR";

    protected static final String NUM_BURN_IN_COPY_RATIO_LONG_NAME = "numBurnInCopyRatio";
    protected static final String NUM_BURN_IN_COPY_RATIO_SHORT_NAME = "numBurnCR";

    protected static final String NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME = "numSamplesAlleleFraction";
    protected static final String NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME = "numSampAF";

    protected static final String NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME = "numBurnInAlleleFraction";
    protected static final String NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME = "numBurnAF";

    protected static final String INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME = "intervalThresholdCopyRatio";
    protected static final String INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME = "simThCR";

    protected static final String INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME = "intervalThresholdAlleleFraction";
    protected static final String INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME = "simThAF";

    protected static final String MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME = "maxNumIterationsSimSeg";
    protected static final String MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_SHORT_NAME = "maxIterSim";

    protected static final String NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_LONG_NAME = "numIterationsSimSegPerFit";
    protected static final String NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_SHORT_NAME = "numIterSimPerFit";

    @Argument(
            doc = "Input file containing denoised copy-ratio profile (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_FULL_NAME,
            shortName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputDenoisedCopyRatiosFile = null;

    @Argument(
            doc = "Input file containing allelic counts (output of CollectAllelicCounts).",
            fullName = CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_FULL_NAME,
            shortName = CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputAllelicCountsFile = null;

    @Argument(
            doc = "Prefix for output files. Will also be used as the sample name in downstream plots." +
                    "(Note: if this is a file path or contains slashes (/), " +
                    "the string after the final slash will be used as the sample name in downstream plots.)",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME,
            shortName = CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Maximum number of segments allowed per chromosome.",
            fullName = MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME,
            shortName = MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_SHORT_NAME,
            minValue = 1,
            optional = true
    )
    private int maxNumSegmentsPerChromosome = 100;

    @Argument(
            doc = "Minimum total count required to include site in allele-fraction segmentation.",
            fullName = MINIMUM_TOTAL_ALLELE_COUNT_LONG_NAME,
            shortName = MINIMUM_TOTAL_ALLELE_COUNT_SHORT_NAME,
            minValue = 0,
            optional = true
    )
    private int minTotalAlleleCount = 10;

    @Argument(
            doc = "Variance of Gaussian kernel for copy-ratio segmentation.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_COPY_RATIO_LONG_NAME,
            shortName = KERNEL_VARIANCE_COPY_RATIO_SHORT_NAME,
            minValue = 0.,
            optional = true
    )
    private double kernelVarianceCopyRatio = 0.;

    @Argument(
            doc = "Variance of Gaussian kernel for allele-fraction segmentation.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME,
            shortName = KERNEL_VARIANCE_ALLELE_FRACTION_SHORT_NAME,
            minValue = 0.,
            optional = true
    )
    private double kernelVarianceAlleleFraction = 1.;

    @Argument(
            doc = "Dimension of kernel approximation.  A subsample containing this number of data points " +
                    "will be taken from the copy-ratio profile and used to construct the approximation for each chromosome.  " +
                    "If the total number of data points in a chromosome is greater " +
                    "than this number, then all data points in the chromosome will be used.  " +
                    "Time complexity scales quadratically and space complexity scales linearly with this parameter.",
            fullName = KERNEL_APPROXIMATION_DIMENSION_LONG_NAME,
            shortName = KERNEL_APPROXIMATION_DIMENSION_SHORT_NAME,
            minValue = 1,
            optional = true
    )
    private int kernelApproximationDimension = 100;

    @Argument(
            doc = "Window sizes to use for calculating local changepoint costs.  " +
                    "For each window size, the cost for each data point to be a changepoint will be calculated " +
                    "assuming that it demarcates two adjacent segments of that size.  " +
                    "Including small (large) window sizes will increase sensitivity to small (large) events.  " +
                    "Duplicate values will be ignored.",
            fullName = WINDOW_SIZES_LONG_NAME,
            shortName = WINDOW_SIZES_SHORT_NAME,
            minValue = 1,
            optional = true
    )
    private List<Integer> windowSizes = new ArrayList<>(Arrays.asList(8, 16, 32, 64, 128, 256));

    @Argument(
            doc = "Factor A for the penalty on the number of changepoints per chromosome for copy-ratio segmentation.  " +
                    "Adds a penalty of the form A * C * [1 + log (N / C)], " +
                    "where C is the number of changepoints in the chromosome, " +
                    "to the cost function for each chromosome.  " +
                    "Must be non-negative.",
            fullName = NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_LONG_NAME,
            shortName = NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_SHORT_NAME,
            minValue = 0.,
            optional = true
    )
    private double numChangepointsPenaltyFactorCopyRatio = 1.;

    @Argument(
            doc = "Factor A for the penalty on the number of changepoints per chromosome for allele-fraction segmentation.  " +
                    "Adds a penalty of the form A * C * [1 + log (N / C)], " +
                    "where C is the number of changepoints in the chromosome, " +
                    "to the cost function for each chromosome.  " +
                    "Must be non-negative.",
            fullName = NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_SHORT_NAME,
            minValue = 0.,
            optional = true
    )
    private double numChangepointsPenaltyFactorAlleleFraction = 1.;

    @Argument(
            doc = "Threshold number of copy-ratio intervals for small-segment merging. " +
                    "If a segment contains strictly less than this number of copy-ratio intervals, " +
                    "it is considered small and will be merged with an adjacent segment.",
            fullName = NUM_COPY_RATIO_INTERVALS_SMALL_SEGMENT_THRESHOLD_LONG_NAME,
            shortName = NUM_COPY_RATIO_INTERVALS_SMALL_SEGMENT_THRESHOLD_SHORT_NAME,
            optional = true
    )
    private int numCopyRatioIntervalsSmallSegmentThreshold = 3;

//    @Argument(
//            doc = "Total number of MCMC samples for copy-ratio model.",
//            fullName = NUM_SAMPLES_COPY_RATIO_LONG_NAME,
//            shortName = NUM_SAMPLES_COPY_RATIO_SHORT_NAME,
//            optional = true
//    )
//    private int numSamplesCopyRatio = 100;
//
//    @Argument(
//            doc = "Number of burn-in samples to discard for copy-ratio model.",
//            fullName = NUM_BURN_IN_COPY_RATIO_LONG_NAME,
//            shortName = NUM_BURN_IN_COPY_RATIO_SHORT_NAME,
//            optional = true
//    )
//    private int numBurnInCopyRatio = 50;
//
//    @Argument(
//            doc = "Total number of MCMC samples for allele-fraction model.",
//            fullName = NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME,
//            shortName = NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME,
//            optional = true
//    )
//    private int numSamplesAlleleFraction = 100;
//
//    @Argument(
//            doc = "Number of burn-in samples to discard for allele-fraction model.",
//            fullName = NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME,
//            shortName = NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME,
//            optional = true
//    )
//    private int numBurnInAlleleFraction = 50;
//
//    @Argument(
//            doc = "Number of 95% credible-interval widths to use for copy-ratio similar-segment merging.",
//            fullName = INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME,
//            shortName = INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME,
//            optional = true
//    )
//    private double intervalThresholdCopyRatio = 4.;
//
//    @Argument(
//            doc = "Number of 95% credible-interval widths to use for allele-fraction similar-segment merging.",
//            fullName = INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME,
//            shortName = INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME,
//            optional = true
//    )
//    private double intervalThresholdAlleleFraction = 2.;
//
//    @Argument(
//            doc = "Maximum number of iterations allowed for similar-segment merging.",
//            fullName = MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME,
//            shortName = MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_SHORT_NAME,
//            optional = true
//    )
//    private int maxNumSimilarSegmentMergingIterations = 10;
//
//    @Argument(
//            doc = "Number of similar-segment--merging iterations per MCMC model refit. " +
//                    "(Increasing this will decrease runtime, but the final number of segments may be higher. " +
//                    "Setting this to 0 will completely disable model refitting between iterations.)",
//            fullName = NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_LONG_NAME,
//            shortName = NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_SHORT_NAME,
//            optional = true
//    )
//    protected int numSimilarSegmentMergingIterationsPerFit = 0;

    //initialize data/segment variables, some of which may be optional
    private CopyRatioCollection denoisedCopyRatios = null;
    private CopyRatioSegmentCollection copyRatioSegments = null;
    private AllelicCountCollection filteredAllelicCounts = null;
    private AlleleFractionSegmentCollection alleleFractionSegments = null;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        validateArguments();

        final String originalLogLevel =
                (ctx.getLocalProperty("logLevel") != null) ? ctx.getLocalProperty("logLevel") : "INFO";
        ctx.setLogLevel("WARN");

        if (inputDenoisedCopyRatiosFile == null && inputAllelicCountsFile == null) {
            throw new UserException("Must provide at least a denoised copy-ratio file or an allelic-counts file.");
        }

        if (inputDenoisedCopyRatiosFile != null) {
            readDenoisedCopyRatios();
            performCopyRatioSegmentation();
            writeSegments(copyRatioSegments, COPY_RATIO_SEGMENTS_FILE_SUFFIX);
        }
        
        if (inputAllelicCountsFile != null) {
            readAndFilterAllelicCounts();
            performAlleleFractionSegmentation();
            writeSegments(alleleFractionSegments, ALLELE_FRACTION_SEGMENTS_FILE_SUFFIX);
        }

        logger.info("Combining available copy-ratio and allele-fraction segments...");
        final CRAFSegmentCollection crafSegments = CRAFSegmentCollection.unionAndMergeSmallSegments(
                copyRatioSegments, denoisedCopyRatios, alleleFractionSegments, filteredAllelicCounts, numCopyRatioIntervalsSmallSegmentThreshold);
        writeSegments(crafSegments, CRAF_SEGMENTS_FILE_SUFFIX);

//        //initial MCMC model fitting performed by ACNVModeller constructor
//        final ACNVModeller modeller = new ACNVModeller(segmentedGenome,
//                numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction, ctx);
//
//        //write initial segments and parameters to file
//        writeACNVModeledSegmentAndParameterFiles(modeller, BEGIN_FIT_FILE_TAG);
//
//        //similar-segment merging (segment files are output for each merge iteration)
//        performSimilarSegmentMergingStep(modeller);
//
//        //write final segments and parameters to file
//        writeACNVModeledSegmentAndParameterFiles(modeller, FINAL_FIT_FILE_TAG);
//
        ctx.setLogLevel(originalLogLevel);
        logger.info("SUCCESS: ModelSegments run complete.");
    }

    private void validateArguments() {
        Utils.nonNull(outputPrefix);
        Utils.validateArg(!(inputDenoisedCopyRatiosFile == null && inputAllelicCountsFile == null),
                "Must provide at least a denoised copy-ratio profile file or an allelic-counts file.");
        if (inputDenoisedCopyRatiosFile != null) {
            IOUtils.canReadFile(inputDenoisedCopyRatiosFile);
        }
        if (inputAllelicCountsFile != null) {
            IOUtils.canReadFile(inputAllelicCountsFile);
        }
    }

    private void readDenoisedCopyRatios() {
        logger.info(String.format("Reading denoised copy-ratio profile file (%s)...", inputDenoisedCopyRatiosFile));
        denoisedCopyRatios = new CopyRatioCollection(inputDenoisedCopyRatiosFile);
    }

    private void performCopyRatioSegmentation() {
        logger.info("Starting copy-ratio segmentation...");
        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;
        copyRatioSegments = new CopyRatioKernelSegmenter(denoisedCopyRatios)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactorCopyRatio, numChangepointsPenaltyFactorCopyRatio);
    }

    private void readAndFilterAllelicCounts() {
        logger.info(String.format("Reading allelic-counts file (%s)...", inputAllelicCountsFile));
        final AllelicCountCollection unfilteredAllelicCounts = new AllelicCountCollection(inputAllelicCountsFile);
        logger.info(String.format("Filtering allelic counts with total count less than %d...", minTotalAlleleCount));
        filteredAllelicCounts = new AllelicCountCollection(
                unfilteredAllelicCounts.getSampleName(),
                unfilteredAllelicCounts.getRecords().stream()
                        .filter(ac -> ac.getRefReadCount() + ac.getAltReadCount() >= minTotalAlleleCount)
                        .collect(Collectors.toList()));
        final File filteredAllelicCountsFile = new File(outputPrefix + FILTERED_ALLELIC_COUNTS_FILE_SUFFIX);
        filteredAllelicCounts.write(filteredAllelicCountsFile);
        logger.info(String.format("Retained %d / %d sites after filtering on total count...",
                filteredAllelicCounts.getRecords().size(), unfilteredAllelicCounts.getRecords().size()));
        logger.info(String.format("Filtered allelic counts written to %s.", filteredAllelicCountsFile));
    }

    private void performAlleleFractionSegmentation() {
        logger.info("Starting allele-fraction segmentation...");
        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;
        alleleFractionSegments = new AlleleFractionKernelSegmenter(filteredAllelicCounts)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceAlleleFraction, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactorAlleleFraction, numChangepointsPenaltyFactorAlleleFraction);
    }

    private void writeSegments(final TSVLocatableCollection<?> segments,
                               final String fileSuffix) {
        final File segmentsFile = new File(outputPrefix + CRAF_SEGMENTS_FILE_SUFFIX);
        segments.write(segmentsFile);
        logger.info(String.format("Segments written to %s.", segmentsFile));
    }
}