package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation;

import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.TSVLocatableCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Represents a legacy allele-fraction segmentation.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionSegmentCollection extends TSVLocatableCollection<AlleleFractionSegment> {
    enum AlleleFractionSegmentTableColumn {
        SAMPLE,
        CONTIG,
        START,
        END,
        NUM_POINTS,
        MEAN_MINOR_ALLELE_FRACTION;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, AlleleFractionSegment> ALLELE_FRACTION_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
        final String contig = dataLine.get(AlleleFractionSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(AlleleFractionSegmentTableColumn.START);
        final int end = dataLine.getInt(AlleleFractionSegmentTableColumn.END);
        final int numPoints = dataLine.getInt(AlleleFractionSegmentTableColumn.NUM_POINTS);
        final double meanMinorAlleleFraction = dataLine.getDouble(AlleleFractionSegmentTableColumn.MEAN_MINOR_ALLELE_FRACTION);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new AlleleFractionSegment(interval, numPoints, meanMinorAlleleFraction);
    };

    private static final BiConsumer<AlleleFractionSegment, DataLine> ALLELE_FRACTION_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER = (alleleFractionSegment, dataLine) ->
            dataLine.append(alleleFractionSegment.getContig())
                    .append(alleleFractionSegment.getStart())
                    .append(alleleFractionSegment.getEnd())
                    .append(alleleFractionSegment.getNumPoints())
                    .append(alleleFractionSegment.getMeanMinorAlleleFraction());

    public AlleleFractionSegmentCollection(final File inputFile) {
        super(inputFile, AlleleFractionSegmentTableColumn.COLUMNS, ALLELE_FRACTION_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, ALLELE_FRACTION_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public AlleleFractionSegmentCollection(final String sampleName,
                                            final List<AlleleFractionSegment> AlleleFractionSegments) {
        super(sampleName, AlleleFractionSegments, AlleleFractionSegmentTableColumn.COLUMNS, ALLELE_FRACTION_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, ALLELE_FRACTION_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }
}