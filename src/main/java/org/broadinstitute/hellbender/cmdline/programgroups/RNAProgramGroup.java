package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Created by gauthier on 11/28/17.
 */
public final class RNAProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return HelpConstants.DOC_CAT_RNA; }
    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_RNA_SUMMARY; }
}