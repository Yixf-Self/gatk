package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

/**
 * Created by slee on 11/20/17.
 */
public final class GermlineCNVCallerIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final String path = "/home/slee/working/gatk/";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument("L", new File(path, "agilent.21-22-X-Y.interval_list").getAbsolutePath())
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
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, path)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test");
        runCommandLine(argsBuilder);
    }
}