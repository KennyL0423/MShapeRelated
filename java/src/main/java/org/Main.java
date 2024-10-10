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
    ITALYPOWERDEMAND("ItalyPowerDemand", 24, 2),
    BME("BME", 128, 3),
    COFFEE("Coffee", 286, 2),
    DODGERLOOPDAY("DodgerLoopDay", 288, 7),
    DODGERLOOPGAME("DodgerLoopGame", 288, 7),
    MEAT("Meat", 448, 3),
    GUNPOINT("GunPoint", 150, 2),
    PLANE("Plane", 144, 7),
    SMOOTHSUBSPACE("SmoothSubspace", 15, 3),
    WINE("Wine", 234, 2),
    COMPUTERS("Computers", 720, 2),
    WORMS("Worms", 900, 5),
    HAPTICS("Haptics", 1092, 5),
    POWERCONS("PowerCons", 144, 2),
    INSECTEPGREGULARTRAIN("InsectEPGRegularTrain", 601, 3),
    FUNGI("Fungi", 201, 18),
    MELBOURNEPEDESTRIAN("MelbournePedestrian", 24, 10),
    CHINATOWN("Chinatown", 24, 2);


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

    static DATASET dataset = DATASET.COMPUTERS;


    public static void main(String[] args) throws IOException {
        String csvFile = "../datasets/time_series/" + dataset.fileName + ".csv";
        int iter = 2;
        long start, end;
        double mean, std;
        double[] RIarray;

        start = System.currentTimeMillis();
        RIarray = new double[iter];
        for (int i = 0; i < iter; i++) {
            List<double[]> timeSeriesData = DataLoader.readTimeSeriesFromCSV(csvFile, dataset.seqLen);
            int[] pred;
            FrOKShape frok = new FrOKShape(timeSeriesData, dataset.seqLen, dataset.clusterNum, 0.6);
            pred = frok.fit();
            int[] truth = DataLoader.readLabelsFromCSV("../labels/" + dataset.fileName + ".csv");
            RIarray[i] = Metric.calculateRandIndex(pred, truth);
        }
        end = System.currentTimeMillis();
        mean = Metric.calMean(RIarray);
        std = Metric.calStd(RIarray, mean);
        System.out.println("FROK Rand Index: " + mean + " +- " + std);
        System.out.println("FROK Time taken: " + (end - start) / iter + "ms");

        start = System.currentTimeMillis();
        RIarray = new double[iter];
        for (int i = 0; i < iter; i++) {
            List<double[]> timeSeriesData = DataLoader.readTimeSeriesFromCSV(csvFile, dataset.seqLen);
            int[] pred;
            wFCM wfcm = new wFCM(timeSeriesData, dataset.seqLen, dataset.clusterNum);
            pred = wfcm.fit();
            int[] truth = DataLoader.readLabelsFromCSV("../labels/" + dataset.fileName + ".csv");
            RIarray[i] = Metric.calculateRandIndex(pred, truth);
        }
        end = System.currentTimeMillis();
        mean = Metric.calMean(RIarray);
        std = Metric.calStd(RIarray, mean);
        System.out.println("wfcm Rand Index: " + mean + " +- " + std);
        System.out.println("wfcm Time taken: " + (end - start) / iter + "ms");


        start = System.currentTimeMillis();
        List<double[]> timeSeriesData = DataLoader.readTimeSeriesFromCSV(csvFile, dataset.seqLen);
        int[] pred;
        PAM pam = new PAM(timeSeriesData, dataset.seqLen, dataset.clusterNum);
        pred = pam.fit();
        int[] truth = DataLoader.readLabelsFromCSV("../labels/" + dataset.fileName + ".csv");

        end = System.currentTimeMillis();
        System.out.println("pam Rand Index: " + Metric.calculateRandIndex(pred, truth));
        System.out.println("pam Time taken: " + (end - start) / iter + "ms");

    }
}