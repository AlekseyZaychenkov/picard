/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.AbstractLocusInfo;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang.ArrayUtils;
import picard.filter.CountingFilter;
import picard.filter.CountingPairedFilter;

import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;

/**
 * Class for collecting data on reference coverage, base qualities and excluded bases from one AbstractLocusInfo object for
 * CollectWgsMetrics.
 * <p>
 * The shared code for forming result for CollectWgsMetrics is abstracted into this class.
 * Classes that extend this collector implement their logic in addInfo() method.
 * @author Mariia_Zueva@epam.com, EPAM Systems, Inc. <www.epam.com>
 */

@XmlRootElement()
public abstract class AbstractWgsMetricsCollector<T extends AbstractRecordAndOffset> implements Serializable {

    /**
     * The source CollectWgsMetrics object
     */
    final transient CollectWgsMetrics collectWgsMetrics;
    /** Count of sites with a given depth of coverage. Includes all but quality 2 bases.
     * We draw depths from this histogram when we calculate the theoretical het sensitivity.
     */

    @XmlElementWrapper
    protected final AtomicLongArray unfilteredDepthHistogramArray;
    /** Count of bases observed with a given base quality. Includes all but quality 2 bases.
     * We draw bases from this histogram when we calculate the theoretical het sensitivity.
     */
    protected final AtomicLongArray unfilteredBaseQHistogramArray;
    /**
     * Count of sites with a given depth of coverage.
     * Excludes bases with quality below MINIMUM_BASE_QUALITY (default 20).
     */
    protected final AtomicLongArray highQualityDepthHistogramArray;
    /**
     * Number of aligned bases that were filtered out because they were of low base quality (default is < 20).
     */
    protected AtomicLong basesExcludedByBaseq = new AtomicLong();
    /**
     * Number of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.
     */
    protected AtomicLong basesExcludedByOverlap = new AtomicLong();
    /**
     * Number of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x).
     */
    protected AtomicLong basesExcludedByCapping = new AtomicLong();
    /**
     * Positions with coverage exceeding this value are treated as if they had coverage at this value
     */
    protected final int coverageCap;

    protected transient IntervalList intervals;
    /*
     * This value indicates that processing will stop after specified int the metric amount of genomic bases.
     */
    private boolean usingStopAfter;
    /**
     * The number of processed genomic bases
     */
    protected long counter = 0;



    public ArrayList getDataSet(){
        ArrayList datalist = new ArrayList();
        datalist.add(this.unfilteredDepthHistogramArray);
        datalist.add(this.unfilteredBaseQHistogramArray);
        datalist.add(this.highQualityDepthHistogramArray);
        datalist.add(this.basesExcludedByBaseq);
        datalist.add(this.basesExcludedByOverlap);
        datalist.add(this.basesExcludedByCapping);
        datalist.add(this.coverageCap);
        datalist.add(this.usingStopAfter);
        datalist.add(this.counter);

        return datalist;
    }

    /**
     * Creates a collector and initializes the inner data structures
     *
     * @param collectWgsMetrics CollectWgsMetrics, that creates this collector
     * @param coverageCap       coverage cap
     */

    





    AbstractWgsMetricsCollector(CollectWgsMetrics collectWgsMetrics, final int coverageCap, final IntervalList intervals) {
        if (coverageCap <= 0) {
            throw new IllegalArgumentException("Coverage cap must be positive.");
        }
        this.collectWgsMetrics = collectWgsMetrics;
        unfilteredDepthHistogramArray = new AtomicLongArray(coverageCap + 1);
        highQualityDepthHistogramArray = new AtomicLongArray(coverageCap + 1);
        unfilteredBaseQHistogramArray = new AtomicLongArray(Byte.MAX_VALUE);
        this.coverageCap    = coverageCap;
        this.intervals      = intervals;
        this.usingStopAfter = collectWgsMetrics.STOP_AFTER > 0;
    }

    /**
     * Accumulates the data from AbstractLocusInfo in inner structures
     * @param info {@link htsjdk.samtools.util.AbstractLocusInfo} with aligned to reference position reads
     * @param ref  {@link htsjdk.samtools.reference.ReferenceSequence}
     * @param referenceBaseN true if current the value of reference base represents a no call
     */
    public abstract void addInfo(final AbstractLocusInfo<T> info, final ReferenceSequence ref, boolean referenceBaseN);

    /**
     * Adds collected metrics and depth histogram to file
     * @param file MetricsFile for result of collector's work
     * @param dupeFilter         counting filter for duplicate reads
     * @param adapterFilter      counting filter for adapter reads
     * @param mapqFilter         counting filter for mapping quality
     * @param pairFilter         counting filter for reads without a mapped mate pair
     */
    public void addToMetricsFile(final MetricsFile<WgsMetrics, Integer> file,
            final boolean includeBQHistogram,
            final CountingFilter dupeFilter,
            final CountingFilter adapterFilter,
            final CountingFilter mapqFilter,
            final CountingPairedFilter pairFilter) {
        final WgsMetrics metrics = getMetrics(dupeFilter, adapterFilter, mapqFilter, pairFilter);

        // add them to the file
        file.addMetric(metrics);
        file.addHistogram(getHighQualityDepthHistogram());
        if (includeBQHistogram) addBaseQHistogram(file);
    }

    protected void addBaseQHistogram(final MetricsFile<WgsMetrics, Integer> file) {
        file.addHistogram(getUnfilteredBaseQHistogram());
    }

    protected Histogram<Integer> getHighQualityDepthHistogram() {
        long[] highQualityDepthHistogramNotAtomicArray = new long[highQualityDepthHistogramArray.length()];
        for(int i=0; i<highQualityDepthHistogramArray.length(); i++)
            highQualityDepthHistogramNotAtomicArray[i] = highQualityDepthHistogramArray.get(i);

        return getHistogram(highQualityDepthHistogramNotAtomicArray, "coverage", "high_quality_coverage_count");
    }

    protected Histogram<Integer> getUnfilteredDepthHistogram() {
        long[] unfilteredDepthHistogramNotAtomicArray = new long[unfilteredDepthHistogramArray.length()];
        for(int i=0; i<unfilteredDepthHistogramArray.length(); i++)
            unfilteredDepthHistogramNotAtomicArray[i] = unfilteredDepthHistogramArray.get(i);

        return getHistogram(unfilteredDepthHistogramNotAtomicArray, "coverage", "unfiltered_coverage_count");
    }

    protected Histogram<Integer> getUnfilteredBaseQHistogram() {
        long[] unfilteredBaseQHistogramNotAtomicArray = new long[unfilteredBaseQHistogramArray.length()];
        for(int i=0; i<unfilteredBaseQHistogramArray.length(); i++)
            unfilteredBaseQHistogramNotAtomicArray[i] = unfilteredBaseQHistogramArray.get(i);

        return getHistogram(unfilteredBaseQHistogramNotAtomicArray, "baseq", "unfiltered_baseq_count");
    }

    protected Histogram<Integer> getHistogram(final long[] array, final String binLabel, final String valueLabel) {
        final Histogram<Integer> histogram = new Histogram<>(binLabel, valueLabel);
        for (int i = 0; i < array.length; ++i) {
            histogram.increment(i, array[i]);
        }
        return histogram;
    }

    /**
     * Creates CollectWgsMetrics.WgsMetrics - the object holding the result of CollectWgsMetrics
     *
     * @param dupeFilter     counting filter for duplicate reads
     * @param mapqFilter     counting filter for mapping quality
     * @param pairFilter     counting filter for reads without a mapped mate pair
     * @return CollectWgsMetrics.WgsMetrics with set fields
     */
    protected WgsMetrics getMetrics(final CountingFilter dupeFilter,
            final CountingFilter adapterFilter,
            final CountingFilter mapqFilter,
            final CountingPairedFilter pairFilter) {

        return collectWgsMetrics.generateWgsMetrics(
                this.intervals,
                getHighQualityDepthHistogram(),
                getUnfilteredDepthHistogram(),
                collectWgsMetrics.getBasesExcludedBy(adapterFilter),
                collectWgsMetrics.getBasesExcludedBy(mapqFilter),
                collectWgsMetrics.getBasesExcludedBy(dupeFilter),
                collectWgsMetrics.getBasesExcludedBy(pairFilter),
                basesExcludedByBaseq.get(),
                basesExcludedByOverlap.get(),
                basesExcludedByCapping.get(),
                coverageCap,
                getUnfilteredBaseQHistogram(),
                collectWgsMetrics.SAMPLE_SIZE);
    }

    /**
     * @return true, of number of processed loci exceeded the threshold, otherwise false
     */
    boolean isTimeToStop(final long processedLoci) {
        return usingStopAfter && processedLoci > collectWgsMetrics.STOP_AFTER - 1;
    }


    public AbstractWgsMetricsCollector combineDataWith(ArrayList datalist){
        AbstractWgsMetricsCollector combinedCollector
                = new CollectWgsMetrics.WgsMetricsCollector(this.collectWgsMetrics, this.coverageCap, this.intervals);

        try {
            AtomicLongArray otherUnfilteredDepthHistogramArray
                    = (AtomicLongArray) datalist.get(0);
            AtomicLongArray otherUnfilteredBaseQHistogramArray
                    = (AtomicLongArray) datalist.get(1);
            AtomicLongArray otherHighQualityDepthHistogramArray
                    = (AtomicLongArray) datalist.get(2);

            AtomicLong otherBasesExcludedByBaseq
                    = (AtomicLong) datalist.get(3);
            AtomicLong otherBasesExcludedByOverlap
                    = (AtomicLong) datalist.get(4);
            AtomicLong basesExcludedByCapping
                    = (AtomicLong) datalist.get(5);

            long otherCounter = (long) datalist.get(6);


            for(int i=0; i<combinedCollector.unfilteredDepthHistogramArray.length(); i++){
                long a = this.unfilteredDepthHistogramArray.get(i);
                long b = otherUnfilteredDepthHistogramArray.get(i);
                combinedCollector.unfilteredDepthHistogramArray.set(i, (a+b));
            }

            for(int i=0; i<combinedCollector.unfilteredDepthHistogramArray.length(); i++){
                long a = this.unfilteredBaseQHistogramArray.get(i);
                long b = otherUnfilteredBaseQHistogramArray.get(i);
                combinedCollector.unfilteredBaseQHistogramArray.set(i, (a+b));
            }

            for(int i=0; i<combinedCollector.highQualityDepthHistogramArray.length(); i++){
                long a = this.highQualityDepthHistogramArray.get(i);
                long b = otherHighQualityDepthHistogramArray.get(i);
                combinedCollector.highQualityDepthHistogramArray.set(i, (a+b));
            }

            long a = this.basesExcludedByBaseq.get();
            long b = otherBasesExcludedByBaseq.get();
            combinedCollector.basesExcludedByBaseq.set((a+b));

            a = this.basesExcludedByOverlap.get();
            b = otherBasesExcludedByOverlap.get();
            combinedCollector.basesExcludedByOverlap.set((a+b));

            a = this.basesExcludedByCapping.get();
            b = basesExcludedByCapping.get();
            combinedCollector.basesExcludedByCapping.set((a+b));

            a = this.counter;
            b = otherCounter;
            combinedCollector.counter = (a+b);

        } catch (ClassCastException cce){
            cce.getMessage();
        }


        return combinedCollector;
    }

    /**
     * Sets the counter to the current number of processed loci. Counter, must be updated
     * from outside, since we are skipping a no call reference positions outside of the collector
     *
     * @param counter number of processed loci
     */
    public void setCounter(long counter) {
        this.counter = counter;
    }

    /**
     * Checks if reference base at given position is unknown.
     *
     * @param position to check the base
     * @param ref      reference sequence
     * @return true if reference base at position represents a no call, otherwise false
     */
    boolean isReferenceBaseN(final int position, final ReferenceSequence ref) {
        final byte base = ref.getBases()[position - 1];
        return SequenceUtil.isNoCall(base);
    }
}
