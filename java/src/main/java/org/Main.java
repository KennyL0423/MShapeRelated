package org;

import org.util.DataLoader;

import java.io.IOException;
import java.util.List;
import org.util.Metric;
public class Main {

    static int seqLen = 140;
    static int clusterNum = 5;

    static String dataset = "ecg";

    static String method = "frok";


    public static void main(String[] args) throws IOException {
        String csvFile = "../datasets/time_series/" + dataset + ".csv";

        int iter = 10;
        long start = System.currentTimeMillis();
        double[] RIarray = new double[iter];

        for (int i = 0; i < iter; i++) {
            List<double[]> timeSeriesData = DataLoader.readTimeSeriesFromCSV(csvFile, seqLen);
            int[] pred;
            switch (method){
                case "wfcm":
                    wFCM wfcm = new wFCM(timeSeriesData, seqLen, clusterNum);
                    pred = wfcm.fit();
                    break;
                case "frok":
                    FrOKShape frok = new FrOKShape(timeSeriesData, seqLen, clusterNum, 0.6);
                    pred = frok.fit();
                    break;
                case "pam":
                    PAM pam = new PAM(timeSeriesData, seqLen, clusterNum);
                    pred = pam.fit();
                    break;
                default:
                    throw new IllegalArgumentException("Invalid method");
            }
            int[] truth = DataLoader.readLabelsFromCSV("../labels/" + dataset + ".csv");
            RIarray[i] = Metric.calculateRandIndex(pred, truth);
        }
        long end = System.currentTimeMillis();

        double mean = Metric.calMean(RIarray);
        double std = Metric.calStd(RIarray, mean);
        System.out.println("Rand Index: " + mean + " +- " + std);

//        DataLoader.writeLabelsIntoCSV("./res/"+ dataset + "-" + method + ".csv", pred);
        System.out.println("Time taken: " + (end - start) / iter + "ms");
    }
}