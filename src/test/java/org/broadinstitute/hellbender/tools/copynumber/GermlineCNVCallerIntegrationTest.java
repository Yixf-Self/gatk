package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Created by slee on 11/20/17.
 */
public final class GermlineCNVCallerIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testCohort() {
        final String path = "/home/slee/working/gatk/";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(new File(path, "TCGA-02-2483-10A-01D-1494-08.counts.hdf5"))
                .addInput(new File(path, "TCGA-02-2485-10A-01D-1494-08.counts.hdf5"))
                .addInput(new File(path, "TCGA-05-4244-10A-01D-1105-08.counts.hdf5"))
                .addInput(new File(path, "TCGA-05-4249-10A-01D-1105-08.counts.hdf5"))
                .addInput(new File(path, "TCGA-05-4250-10A-01D-1105-08.counts.hdf5"))
                .addInput(new File(path, "TCGA-05-4382-10A-01D-1265-08.counts.hdf5"))
                .addInput(new File(path, "TCGA-05-4384-10A-01D-1753-08.counts.hdf5"))
                .addInput(new File(path, "TCGA-05-4389-10A-01D-1265-08.counts.hdf5"))
                .addInput(new File(path, "TCGA-05-4390-10A-01D-1753-08.counts.hdf5"))
                .addInput(new File(path, "TCGA-05-4395-10A-01D-1265-08.counts.hdf5"))
                .addArgument("L", new File(path, "agilent.22.interval_list").getAbsolutePath())
                .addArgument(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, new File(path, "Homo_sapiens_assembly19.dict").getAbsolutePath())
                .addArgument(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME, path + "test-ploidy-calls")
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, path)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-caller")
                .addArgument("imr", IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addArgument(StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG");
        runCommandLine(argsBuilder);
    }

    @Test
    public void testCase() {
        final String path = "/home/slee/working/gatk/";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(new File(path, "TCGA-02-2483-10A-01D-1494-08.counts.hdf5"))
                .addArgument(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME, path + "test-ploidy-calls")
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, path)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-caller-case")
                .addArgument(CopyNumberStandardArgument.MODEL_LONG_NAME, new File(path, "test-caller-model").getAbsolutePath())
                .addArgument("imr", IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addArgument(StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG");
        runCommandLine(argsBuilder);
    }
}