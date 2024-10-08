package org;

import org.util.DataLoader;

import java.io.IOException;
import java.util.List;
import org.util.Metric;

enum DATASET{
    AIR("air", 166, 3),
    ECG("ecg", 140, 5),
    WAFER("Wafer", 152, 2),
    PATTERN("TwoPatterns", 128, 4),
    FACES("FacesUCR", 131, 14),
    ADIAC("Adiac", 176, 37),
    CBF("CBF", 128, 3),
    DIATOMSIZEREDUCTION("DiatomSizeReduction", 345, 4),
    DISTALPHALANXOUTLINEAGEGROUP("DistalPhalanxOutlineAgeGroup", 80, 3),
    DISTALPHALANXOUTLINECORRECT("DistalPhalanxOutlineCorrect", 80, 2),
    EARTHQUAKES("Earthquakes", 512, 2),
    ELECTRICDEVICES("ElectricDevices", 96, 7),
    FORDA("FordA", 500, 2),
    FORDB("FordB", 500, 2),
    FREEZERREGULARTRAIN("FreezerRegularTrain", 301, 2),
    HAM("Ham", 431, 2),
    ITALYPOWERDEMAND("ItalyPowerDemand", 24, 2);


    int seqLen;
    int clusterNum;
    String fileName;
    DATASET(String fileName, int seqLen, int clusterNum){
        this.fileName = fileName;
        this.seqLen = seqLen;
        this.clusterNum = clusterNum;
    }
}
public class Main {

    static String method = "pam";
    static DATASET dataset = DATASET.ELECTRICDEVICES;


    public static void main(String[] args) throws IOException {
        String csvFile = "../datasets/time_series/" + dataset.fileName + ".csv";

        int iter = 1;
        long start = System.currentTimeMillis();
        double[] RIarray = new double[iter];

        for (int i = 0; i < iter; i++) {
            List<double[]> timeSeriesData = DataLoader.readTimeSeriesFromCSV(csvFile, dataset.seqLen);
            int[] pred;
            switch (method){
                case "wfcm":
                    wFCM wfcm = new wFCM(timeSeriesData, dataset.seqLen, dataset.clusterNum);
                    pred = wfcm.fit();
                    break;
                case "frok":
                    FrOKShape frok = new FrOKShape(timeSeriesData, dataset.seqLen, dataset.clusterNum, 0.6);
                    pred = frok.fit();
                    break;
                case "pam":
                    PAM pam = new PAM(timeSeriesData, dataset.seqLen, dataset.clusterNum);
                    pred = pam.fit();
                    break;
                default:
                    throw new IllegalArgumentException("Invalid method");
            }
            int[] truth = DataLoader.readLabelsFromCSV("../labels/" + dataset.fileName + ".csv");
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