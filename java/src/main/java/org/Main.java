package org;

import org.util.DataLoader;

import java.io.IOException;
import java.util.List;
import org.util.Metric;
public class Main {

    static int seqLen = 166;
    static int clusterNum = 3;


    public static void main(String[] args) throws IOException {
        long start = System.currentTimeMillis();
        String csvFile = "/Users/suyx1999/ExpData/shape/air.csv";
        List<double[]> timeSeriesData = DataLoader.readTimeSeriesFromCSV(csvFile, seqLen);

        wFCM clustering = new wFCM(timeSeriesData, seqLen, clusterNum);
//        FrOKShape clustering = new FrOKShape(timeSeriesData, seqLen, clusterNum, 0.6);
//        PAM pam = new PAM(timeSeriesData, seqLen, clusterNum);


        int[] pred = clustering.fit();
        long end = System.currentTimeMillis();


        int[] truth = DataLoader.readLabelsFromCSV("./groundtruth/air.csv");
        double ri = Metric.calculateRandIndex(pred, truth);
        System.out.println("Rand Index: " + ri);

//        DataLoader.writeLabelsIntoCSV("./res/air-wfcm.csv", pred);
        System.out.println("Time taken: " + (end - start) + "ms");
    }
}