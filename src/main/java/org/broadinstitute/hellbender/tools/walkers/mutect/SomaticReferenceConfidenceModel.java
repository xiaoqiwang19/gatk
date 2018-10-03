package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RefVsAnyResult;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceResult;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public class SomaticReferenceConfidenceModel extends ReferenceConfidenceModel {

    private final SampleList samples;
    private final SomaticGenotypingEngine genotypingEngine;



    /**
     * Create a new ReferenceConfidenceModel
     *
     * @param samples the list of all samples we'll be considering with this model
     * @param header the SAMFileHeader describing the read information (used for debugging)
     * @param indelInformativeDepthIndelSize the max size of indels to consider when calculating indel informative depths
     */
    SomaticReferenceConfidenceModel(final SampleList samples,
                                    final SAMFileHeader header,
                                    final int indelInformativeDepthIndelSize,
                                    final SomaticGenotypingEngine genotypingEngine){
        super(samples, header, indelInformativeDepthIndelSize, 0);
        this.samples = samples;
        this.genotypingEngine = genotypingEngine;
    }

        /**
         * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
         *
         * @param ploidy target sample ploidy.
         * @param pileup the read backed pileup containing the data we want to evaluate
         * @param refBase the reference base at this pileup position
         * @param qual the min base quality for a read in the pileup at the pileup position to be included in the calculation
         * @param hqSoftClips running average data structure (can be null) to collect information about the number of high quality soft clips
         * @return a RefVsAnyResult genotype call.
         */
    @Override
    public ReferenceConfidenceResult calcGenotypeLikelihoodsOfRefVsAny(final int ploidy,
                                                                       final ReadPileup pileup,
                                                                       final byte refBase,
                                                                       final byte qual,
                                                                       final MathUtils.RunningAverage hqSoftClips,
                                                                       final boolean readsWereRealigned) {

        final SomaticRefVsAnyResult result = new SomaticRefVsAnyResult();
        final Map<String, List<GATKRead>> perSampleReadMap = new HashMap<>();
        perSampleReadMap.put(samples.getSample(0), pileup.getReads());
        final ReadLikelihoods<Allele> readLikelihoods = new ReadLikelihoods<>(samples, new IndexedAlleleList<>(Arrays.asList(Allele.create(refBase,true), Allele.NON_REF_ALLELE)), perSampleReadMap);
        final Iterator<PileupElement> pileupIter = pileup.iterator();
        for (int i = 0; i < pileup.size(); i++) {
            final PileupElement element = pileupIter.next();
            final boolean isAlt = readsWereRealigned ? isAltAfterAssembly(element, refBase) : isAltBeforeAssembly(element, refBase);
            final double nonRefLikelihood;
            if (isAlt) {
                nonRefLikelihood = QualityUtils.qualToProbLog10(qual);
                result.nonRefDepth++;
            } else {
                nonRefLikelihood = QualityUtils.qualToErrorProbLog10(qual) + MathUtils.LOG10_ONE_THIRD;
                result.refDepth++;
            }
            readLikelihoods.sampleMatrix(0).set(0, i, nonRefLikelihood);
        }
        result.lods = genotypingEngine.somaticLog10Odds(readLikelihoods.sampleMatrix(0));
        return result;
    }

    @Override
    public void addGenotypeData(final ReferenceConfidenceResult result, final GenotypeBuilder gb) {
        gb.attribute(GATKVCFConstants.TUMOR_LOD_KEY, ((SomaticRefVsAnyResult)result).lods.get(Allele.NON_REF_ALLELE));
    }

    @Override
    public void doIndelRefConfCalc(final int ploidy, final byte[] ref, final ReadPileup pileup, final int refOffset, final ReferenceConfidenceResult homRefCalc) {  }

    }
