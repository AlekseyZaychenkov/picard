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

import com.sun.management.OperatingSystemMXBean;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.*;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.barclay.argparser.Argument;
import picard.PicardException;
import picard.filter.CountingFilter;
import picard.filter.CountingPairedFilter;
import shaded.cloud_nio.com.google.api.client.util.DateTime;

import java.io.*;
import java.lang.management.ManagementFactory;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.UnknownHostException;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;
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

    /*
     program attempts call garbage collector one time per:
     */
    private static final long GC_ATTEMPT_OF_CALL_FREQUENCY = 5*60*1000; // in milliseconds


    /*
     garbage collector was already called GC_CALLED_TIMES times
    */
    private static int GC_CALLED_TIMES = 0;




    /*
    temporal variable for calling garbage collector
    */
    private static long previousCallTime = System.currentTimeMillis();
    private static long firstCallTime = System.currentTimeMillis();



    private final Log log = Log.getInstance(WgsMetricsProcessorImpl.class);


    boolean DISTRIBUTED_COMPUTING;

    boolean IS_SERVER;

    IntervalList intervalList;

    /**
     * @param iterator  input {@link htsjdk.samtools.util.AbstractLocusIterator}
     * @param refWalker over processed reference file
     * @param collector input {@link picard.analysis.AbstractWgsMetricsCollector}
     * @param progress  logger
     */
    public WgsMetricsProcessorImpl(AbstractLocusIterator<T, AbstractLocusInfo<T>> iterator,
            ReferenceSequenceFileWalker refWalker,
            AbstractWgsMetricsCollector<T> collector,
            ProgressLogger progress,
            boolean DISTRIBUTED_COMPUTING,
            boolean IS_SERVER,
            IntervalList intervalList) {
        this.iterator = iterator;
        this.collector = collector;
        this.refWalker = refWalker;
        this.progress = progress;
        this.DISTRIBUTED_COMPUTING = DISTRIBUTED_COMPUTING;
        this.IS_SERVER = IS_SERVER;
        this.intervalList = intervalList;
    }

    /**
     * Method gets the data from iterator for each locus and processes it with the help of collector.
     */
    @Override
    public void processFile() {
        final AtomicLong counter = new AtomicLong();
        final List<SamLocusIterator.LocusInfo> records = new ArrayList<>(READS_IN_PACK);
     //   ExecutorService service = Executors.newFixedThreadPool(16);

        com.sun.management.OperatingSystemMXBean osBean = ManagementFactory.getPlatformMXBean(OperatingSystemMXBean.class);
        long pack_counter = 0;

        // if the program working in mode of DISTRIBUTED_COMPUTING
        if(DISTRIBUTED_COMPUTING==true) {
            log.info("DISTRIBUTED_COMPUTING");
            String args[] = {"", ""};

            // if the program working in server mode
            if(IS_SERVER == false) {
                String[] address = Arrays.asList("127.0.0.1", "4322").toArray(new String[0]);

                if (address.length != 2) {
                    System.err.println(
                            "Usage: java EchoClient <host name> <port number>");
                }

                String hostName = address[0];
                int portNumber = Integer.parseInt(address[1]);

                try (Socket echoSocket = new Socket(hostName, portNumber);
                     ObjectOutputStream out =
                             new ObjectOutputStream(echoSocket.getOutputStream());
                     BufferedReader in =
                             new BufferedReader(
                                     new InputStreamReader(echoSocket.getInputStream()));
                     BufferedReader stdIn =
                             new BufferedReader(
                                     new InputStreamReader(System.in))
                ) {
                    String userInput;
                    while ((userInput = stdIn.readLine()) != null) {

                        out.writeObject((Object) new Integer(0));
                        out.flush();
                        String receivedWord = in.readLine();
                        System.out.println("echo: " + receivedWord);
                        if (receivedWord.equals("close"))
                            break;
                    }
                    System.out.println("The communication finished");
                } catch (UnknownHostException e) {
                    System.err.println("Don't know about host " + hostName);
                    System.exit(1);
                } catch (IOException e) {
                    e.printStackTrace();
                    System.err.println("Couldn't get I/O for the connection to " +
                            hostName);
                    System.exit(1);
                }

            } else {


                args = Arrays.asList("4322").toArray(new String[0]);
                if (args.length != 1) {
                    System.err.println("Usage: java EchoServer <port number>");
                }

                int portNumber = Integer.parseInt(args[0]);
                try {
                    ServerSocket serverSocket = new ServerSocket(Integer.parseInt(args[0]));
                    ExecutorService executorService = Executors.newCachedThreadPool();
                    // for (int i = 0; i < 3; i++) {

                    executorService.submit(() -> {
                        try (
                                Socket clientSocket = serverSocket.accept();
                                PrintWriter out = new PrintWriter(clientSocket.getOutputStream(), true);
                                ObjectInputStream in = new ObjectInputStream(new BufferedInputStream(clientSocket.getInputStream()))
                        ) {
                            String inputLine = "none";

                            Object objectInput = null;
                            //System.out.println(in.read());

                            while ((objectInput = in.readObject()) != null) {

                                if (objectInput instanceof String) {
                                    System.out.println("Is String");
                                    inputLine = (String) objectInput;
                                } else {
                                    System.out.println("Not String");
                                }

                                if (inputLine.equals("close"))
                                    break;

                                out.println("***toFormAnswer for "+inputLine);
                            }
                            //out.println("Closing connection with server");
                        } catch (IOException e) {
                            System.out.println("Exception caught when trying to listen on port " + portNumber + " or listening for a connection");
                            System.out.println(e.getMessage());
                        } catch (ClassNotFoundException ec) {
                            System.out.println("ClassNotFoundException");
                            System.out.println(ec.getMessage());
                        }
                    });
                    //}
                    executorService.awaitTermination(5, TimeUnit.MINUTES);
                    serverSocket.close();


                } catch (IOException | InterruptedException e) {
                    System.out.println("Exception caught when trying create new Cached Thread Pool or await termination");
                    System.out.println(e.getMessage());
                }
            }


        // if the program working in NOT mode of DISTRIBUTED_COMPUTING
        } else {


            System.out.println("intervalList.toString() : "+intervalList.toString());
            List<Interval> intervalCollectionList = intervalList.getIntervals();
            int c = 0;
           /* for(Interval i : intervalCollectionList) {
                System.out.println("( " + (++c) + " )");
                System.out.println("i.toString() : " + i.toString());
                System.out.println("i.getStart() : " + i.getStart());
                System.out.println("i.length() : " + i.length());
                System.out.println("i.getContig() : " + i.getContig());
                System.out.println("i.countBases(intervalCollectionList) : " + i.countBases(intervalCollectionList));
*/




            while (iterator.hasNext()) {
                //synchronized (collector) {
                    final AbstractLocusInfo<T> info = iterator.next();
                    final ReferenceSequence ref = refWalker.get(info.getSequenceIndex());


                    final boolean referenceBaseN = collector.isReferenceBaseN(info.getPosition(), ref);

                    collector.addInfo(info, ref, referenceBaseN);

                    if (referenceBaseN) {
                        continue;
                    }

                    progress.record(info.getSequenceName(), info.getPosition());
              //  }

                    if (collector.isTimeToStop(counter.incrementAndGet())) {
                        break;
                    }
                    collector.setCounter(counter.get());

            }




               /* while (iterator.hasNext()) {


                    AbstractLocusInfo info = iterator.next();
                    pack_counter++;
                    //  long placeToStart = 10 * 100 * 1000 * 1000L + 1;
                    //   if (pack_counter >= placeToStart) {
                    // System.out.println(pack_counter +" >= "+ placeToStart+" : "+(pack_counter >= placeToStart));

                    ReferenceSequence ref = refWalker.get(info.getSequenceIndex());


                    /* Multithreading */
                    // Check that the reference is not N


                /*




                    final byte base = ref.getBases()[info.getPosition() - 1];

                    if ((SequenceUtil.isNoCall(base)))
                        continue;

                    records.add((SamLocusIterator.LocusInfo) info);

                    // Record progress
                    progress.record(info.getSequenceName(), info.getPosition());
                    if ((records.size() < READS_IN_PACK && iterator.hasNext()))
                        continue;

                    final List<SamLocusIterator.LocusInfo> tmpRecords = records;
                    records = new ArrayList<>(READS_IN_PACK);


                    // Add read pack to the collector in a separate stream
                  /*  service.submit(new Runnable() {
                        @Override
                        public void run() {*/


                   /*



                            for (SamLocusIterator.LocusInfo rec : tmpRecords) {

                                boolean referenceBaseN = collector.isReferenceBaseN(info.getPosition(), ref);
                                collector.addInfo((AbstractLocusInfo) rec, ref, referenceBaseN);
                                if (referenceBaseN) {
                                    continue;
                                }
                                //  progress.record(info.getSequenceName(), info.getPosition());
                            }
                        /*}
                    });*/

                    // Perhaps stop
                    // if (usingStopAfter && counter > stopAfter) break;


                    // Record progress

            /*


                    if (collector.isTimeToStop(counter.incrementAndGet())) {
                        break;
                    }
                    collector.setCounter(counter.get());
                    //  }



               /* if ((System.currentTimeMillis() - previousCallTime) > GC_ATTEMPT_OF_CALL_FREQUENCY) {
                     System.gc();
                     log.info("Attempt to call the garbage collector");
                     GC_CALLED_TIMES++;
                     previousCallTime = System.currentTimeMillis();
                }*/


               /*     if ((System.currentTimeMillis() - previousCallTime) > 30 * 1000) {
                        //service2 = (ThreadPoolExecutor)service;
                        double diff = (((double) System.currentTimeMillis() - (double) firstCallTime) / 60000);
                        previousCallTime = System.currentTimeMillis();
                        System.out.println("Elapsed time in WgsProcessor: " + diff + " min");
                        System.out.println("info.getStart(): " + info.getStart());
                        System.out.println("pack_counter: " + pack_counter);
*/


                   /* log.info("service2.getQueue().size() "+service2.getQueue().size());
                    log.info("service2.getActiveCount() "+service2.	getActiveCount());
                    log.info("service2.getCompletedTaskCount() "+service2.getCompletedTaskCount());
                   // log.info("service2.getKeepAliveTime(NANOSECONDS) "+service2.getKeepAliveTime(TimeUnit.NANOSECONDS));

                 //   System.out.println("ProcessCpuLoad: "+osBean.getProcessCpuLoad() * 100);
                //    System.out.println("SystemCpuLoad: "+osBean.getSystemCpuLoad()  * 100);
*/
                 //   }

             //   }


          /*  }*/
         //   service.shutdown();








          /*  try {
                service.awaitTermination(1, TimeUnit.DAYS);
            } catch (InterruptedException ex) {
                Logger.getLogger(SinglePassSamProgram.class.getName()).log(Level.SEVERE,
                        null, ex);
            }
*/

        }

/*



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
        AtomicLongArray unfilteredBaseQHistogramArray = collector.unfilteredBaseQHistogramArray;
        long sum = 0;
        for(int i=0; i<unfilteredBaseQHistogramArray.length(); i++)
            sum += unfilteredBaseQHistogramArray.get(i);
        final long sumBaseQ = sum;

        AtomicLongArray unfilteredDepthHistogramArray = collector.unfilteredDepthHistogramArray;
        long[] unfilteredDepthHistogramNotAtomicArray = new long[unfilteredDepthHistogramArray.length()];
        for(int i=0; i<unfilteredDepthHistogramArray.length(); i++)
            unfilteredDepthHistogramNotAtomicArray[i] = unfilteredDepthHistogramArray.get(i);

        final long sumDepthHisto = LongStream.rangeClosed(0, collector.coverageCap).map(i -> (i * unfilteredDepthHistogramNotAtomicArray[(int)i])).sum();
        if (sumBaseQ != sumDepthHisto) {
            log.error("Coverage and baseQ distributions contain different amount of bases!");
        }

        */




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
