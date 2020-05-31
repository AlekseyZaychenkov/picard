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
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.*;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.barclay.argparser.Argument;
import picard.PicardException;
import picard.filter.CountingFilter;
import picard.filter.CountingPairedFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.LongStream;

/**
 * Implementation of {@link picard.analysis.WgsMetricsProcessor} that gets input data from a given iterator
 * and processes it with a help of collector
 * @author Mariia_Zueva@epam.com, EPAM Systems, Inc. <www.epam.com>
 */
public class WgsMetricsProcessorImpl<T extends AbstractRecordAndOffset> implements WgsMetricsProcessor {
    /**
     * Source of input data
     */
    private final AbstractLocusIterator<T, AbstractLocusInfo<T>> iterator;
    /**
     * Accumulates the data from iterator
     */
    private final AbstractWgsMetricsCollector<T> collector;
    /**
     * ReferenceWalker for a processed reference sequence
     */
    private final ReferenceSequenceFileWalker refWalker;
    /**
     * Logger for the progress of work
     */
    private final ProgressLogger progress;


  /*  Needed for multithreading. Amount of reads in the pack.")  */
    public static final int READS_IN_PACK = 1000;  // Multithreading


    private final Log log = Log.getInstance(WgsMetricsProcessorImpl.class);

    /**
     * @param iterator  input {@link htsjdk.samtools.util.AbstractLocusIterator}
     * @param refWalker over processed reference file
     * @param collector input {@link picard.analysis.AbstractWgsMetricsCollector}
     * @param progress  logger
     */
    public WgsMetricsProcessorImpl(AbstractLocusIterator<T, AbstractLocusInfo<T>> iterator,
            ReferenceSequenceFileWalker refWalker,
            AbstractWgsMetricsCollector<T> collector,
            ProgressLogger progress) {
        this.iterator = iterator;
        this.collector = collector;
        this.refWalker = refWalker;
        this.progress = progress;
    }

    /**
     * Method gets the data from iterator for each locus and processes it with the help of collector.
     */
    @Override
    public void processFile() {
        long counter = 0;
        List<SamLocusIterator.LocusInfo> records = new ArrayList<>(READS_IN_PACK);
        ExecutorService service = Executors.newCachedThreadPool();

        while (iterator.hasNext()) {
            AbstractLocusInfo info = iterator.next();
            ReferenceSequence ref = refWalker.get(info.getSequenceIndex());
            //System.out.println("<=========Flag 01===========>   "+info);

            /* Multithreading */
            // Check that the reference is not N
            final byte base = ref.getBases()[info.getPosition() - 1];
            if (SequenceUtil.isNoCall(base))
                continue;

            records.add((SamLocusIterator.LocusInfo)info);


            // Record progress
            progress.record(info.getSequenceName(), info.getPosition());
            if (records.size() < READS_IN_PACK && iterator.hasNext())
            {
                continue;
            }

            final List<SamLocusIterator.LocusInfo> tmpRecords = records;
            records = new ArrayList<>(READS_IN_PACK);

            // Add read pack to the collector in a separate stream
            service.submit(new Runnable() {
                @Override
                public void run() {
                    for (SamLocusIterator.LocusInfo rec : tmpRecords) {
                        boolean referenceBaseN = collector.isReferenceBaseN(info.getPosition(), ref);
                        collector.addInfo((AbstractLocusInfo) rec, ref, referenceBaseN);
                        if (referenceBaseN) {
                            continue;
                        }
                      //  progress.record(info.getSequenceName(), info.getPosition());
                    }
                }
            });
            // Perhaps stop
           // if (usingStopAfter && counter > stopAfter) break;

            // Record progress
            if (collector.isTimeToStop(++counter)) {
                break;
            }

            collector.setCounter(counter);
        }
        service.shutdown();


        try {
            service.awaitTermination(1, TimeUnit.DAYS);
        } catch (InterruptedException ex) {
            Logger.getLogger(SinglePassSamProgram.class.getName()).log(Level.SEVERE,
                    null, ex);
        }


        long unfilteredBaseQHistogramSum = 0;
        long unfilteredDepthHistogramSum = 0;
        for (int i = 0; i < collector.unfilteredBaseQHistogramArray.length(); ++i) {
            unfilteredBaseQHistogramSum += collector.unfilteredBaseQHistogramArray.get(i);
        }
        for (int i = 0; i <= collector.coverageCap; ++i) {
            unfilteredDepthHistogramSum += i*collector.unfilteredDepthHistogramArray.get(i);
        }
        if (unfilteredBaseQHistogramSum != unfilteredDepthHistogramSum) {
            throw new PicardException("updated coverage and baseQ distributions unequally");
        }


        // check that we added the same number of bases to the raw coverage histogram and the base quality histograms
        List arr1 = Arrays.asList(collector.unfilteredBaseQHistogramArray);
        final long sumBaseQ = Arrays.stream(ArrayUtils.toPrimitive((Long[])arr1.toArray())).sum();

        List arr2 = Arrays.asList(collector.unfilteredDepthHistogramArray);
        final long sumDepthHisto = LongStream.rangeClosed(0, collector.coverageCap).map(i -> (i * ArrayUtils.toPrimitive((Long[])arr2.toArray())[(int) i])).sum();
        if (sumBaseQ != sumDepthHisto) {
            log.error("Coverage and baseQ distributions contain different amount of bases!");
        }
    }

    @Override
    public void addToMetricsFile(MetricsFile<WgsMetrics, Integer> file,
            boolean includeBQHistogram,
            CountingFilter dupeFilter,
            CountingFilter adapterFilter,
            CountingFilter mapqFilter,
            CountingPairedFilter pairFilter) {
        collector.addToMetricsFile(file, includeBQHistogram, dupeFilter, adapterFilter, mapqFilter, pairFilter);
    }
}
