package org;

import org.util.DataLoader;
import org.util.LinearRegression;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class wFCM {
    int MAX_ITERATIONS = 100;
    double THRESHOLD = 1e-4;
    int M = 2; // Fuzzification coefficient
    int alpha = 5;
    int qmax = 20;
    int seqLen;
    int clusterNum;
    List<double[]> data;
    List<List<TIG>> tigsOfAllData;
    List<List<TIG>> prototypes;
    int[] labels;


    public static void main(String[] args) {
        long start = System.currentTimeMillis();
        String csvFile = "/Users/suyx1999/ExpData/shape/air.csv";
//        String csvFile = "/Users/suyx1999/Downloads/jinfeng.csv";

        List<double[]> timeSeriesData = DataLoader.readTimeSeriesFromCSV(csvFile, 166);
        wFCM clustering = new wFCM(timeSeriesData, 166, 3);
        int[] clusterLabels = clustering.fit();

        long end = System.currentTimeMillis();
        System.out.println("Time taken: " + (end - start) + "ms");
    }

    public wFCM(List<double[]> data, int seqLen, int cluserNum) {
        this.data = data;
        this.seqLen = seqLen;
        this.clusterNum = cluserNum;
        this.tigsOfAllData = new ArrayList<>();
        this.prototypes = new ArrayList<>();
        this.labels = new int[data.size()];
    }

    public int[] fit(){
        int optimalPartNum = findOptimalPartNum();
        System.out.println("Optimal part number: " + optimalPartNum);
        tigsOfAllData = transformTIG(optimalPartNum);
        double[][] U = initialize();
        fuzzyClustering(U);
        return labels;
    }


    private int findOptimalPartNum(){
        double minCost = Double.MAX_VALUE;
        double lambda = 4.0;
        int optimalNum = -1;
        for (int i = 2; i <= qmax; i++){
            double cost = 0.0;
            List<List<TIG>> tigMatrix = transformTIG(i);
            for (List<TIG> tigArray: tigMatrix){
                for (TIG tig: tigArray){
                    cost += tig.tigL * Math.abs(tig.tigInterval[0] - tig.tigInterval[1]);
                }
                cost += lambda * i * Math.log(seqLen);
            }
            if (cost < minCost){
                minCost = cost;
                optimalNum = i;
            }
        }
        return optimalNum;
    }

    private List<List<TIG>> transformTIG(int partNum){
        int l = data.get(0).length;
        List<List<TIG>> tigMatrix = new ArrayList<>();
        for (double[] seq: data){
            List<TIG> tigArray = new ArrayList<>();
            for (int s = 0; s < partNum; s++){
                double[] part = Arrays.copyOfRange(seq, s * l / partNum, (s + 1) * l / partNum);
                LinearRegression lr = new LinearRegression(part);
                TIG tig = new TIG(lr.getSlope(), lr.findInfoGran(), part.length);
                tigArray.add(tig);
            }
            tigMatrix.add(tigArray);
        }
        return tigMatrix;
    }

    private double[][] initialize() {
        Random random = new Random();
        for (int i = 0; i < labels.length; i++) {
            labels[i] = random.nextInt(clusterNum);
        }

        // init prototypes
        List<List<Integer>> indexesByClass = new ArrayList<>();
        for (int i = 0; i < clusterNum; i++) {
            indexesByClass.add(new ArrayList<>());
        }
        for (int i = 0; i < labels.length; i++) {
            int label = labels[i];
            indexesByClass.get(label).add(i);
        }
        for (int i = 0; i < clusterNum; i++) {
            List<Integer> indexes = indexesByClass.get(i);
            int randomIndex = random.nextInt(indexes.size());
            prototypes.add(tigsOfAllData.get(indexes.get(randomIndex)));
        }

        //init matrix U
        double[][] U = new double[tigsOfAllData.size()][clusterNum];
        for (int i = 0; i < tigsOfAllData.size(); i++) {
            for (int j = 0; j < clusterNum; j++) {
                U[i][j] = calculateUij(i, j, prototypes);
            }
        }
        return U;
    }


    // Fuzzy clustering algorithm based on DTW
    private void fuzzyClustering(double[][] U) {
        int n = tigsOfAllData.size();

        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            double[][] newU = new double[n][clusterNum];
            List<List<TIG>> newPrototypes = new ArrayList<>();

            // Step 1: update prototypes
            for (int k = 0; k < clusterNum; k++) {
                List<List<TIG>> membersInCluster = new ArrayList<>();
                List<Double> u = new ArrayList<>();
                for (int i = 0; i < n; i++) {
                    if (labels[i] == k) {
                        membersInCluster.add(tigsOfAllData.get(i));
                        u.add(U[i][k]);
                    }
                }
                newPrototypes.add(gbaryCenter(membersInCluster, u));
            }

            // Step 2: update membership matrix U
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < clusterNum; j++) {
                    newU[i][j] = calculateUij(i, j, newPrototypes);
                }
            }

            // Step 3: Update labels
            for (int i = 0; i < n; i++) {
                double maxU = Double.NEGATIVE_INFINITY;
                int maxIdx = -1;
                for (int j = 0; j < clusterNum; j++) {
                   if (newU[i][j] > maxU){
                       maxU = newU[i][j];
                       maxIdx = j;
                   }
                }
                labels[i] = maxIdx;
            }

            // Step 4: Check for convergence
            double maxChange = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < clusterNum; j++) {
                    maxChange = Math.max(maxChange, Math.abs(newU[i][j] - U[i][j]));
                }
            }

            //Step 5: Iterate U and prototypes
            U = newU;
            prototypes = newPrototypes;

            System.out.println("Iteration " + iter + ", maxChange = " + maxChange);
            if (maxChange < THRESHOLD) {
                System.out.println("Converged after " + iter + " iterations.");
                break;
            }
        }
    }


    // U_ij between x_i and c_j
    private double calculateUij(int i, int j, List<List<TIG>> centroids){
        double sum = 0.0;
        double numerator = Math.pow(TIG.DTW_Omega(tigsOfAllData.get(i), centroids.get(j)), 2) +
                alpha * Math.pow(TIG.DTW_Trend(tigsOfAllData.get(i), centroids.get(j)), 2);
        double denominator;
        for (int c = 0; c < clusterNum; c++){
            denominator = Math.pow(TIG.DTW_Omega(tigsOfAllData.get(i), centroids.get(c)), 2) +
                    alpha * Math.pow(TIG.DTW_Trend(tigsOfAllData.get(i), centroids.get(c)), 2);
            sum += Math.pow(numerator / denominator,  1.0 / (M - 1));
        }
        if (Double.isNaN(sum)){
            return 1.0;
        }
        return 1.0 / sum;

    }

    private List<TIG> gbaryCenter(List<List<TIG>> members, List<Double> u){
        int tigNum = members.get(0).size();
        double[] sumTigK = new double[tigNum];
        double[] sumTigInterval0 = new double[tigNum];
        double[] sumTigInterval1 = new double[tigNum];
        double denominator = 0.0;
        for (int i = 0; i < members.size(); i++){
            List<TIG> member = members.get(i);
            for (int j = 0; j < member.size(); j++){
                sumTigK[j] += Math.pow(u.get(i), M) * member.get(j).tigK;
                sumTigInterval0[j] += Math.pow(u.get(i), M) * member.get(j).tigInterval[0];
                sumTigInterval1[j] += Math.pow(u.get(i), M) * member.get(j).tigInterval[1];
            }
            denominator += Math.pow(u.get(i), M);
        }
        List<TIG> newPrototype = new ArrayList<>();
        for (int i = 0; i < tigNum; i++){
            sumTigK[i] /= denominator;
            sumTigInterval0[i] /= denominator;
            sumTigInterval1[i] /= denominator;
            newPrototype.add(new TIG(sumTigK[i], new double[]{sumTigInterval0[i], sumTigInterval1[i]}, members.get(0).get(i).tigL));
        }
        return newPrototype;
    }

    private static class TIG{
        public double tigK; //slope
        public double[] tigInterval; // information granularity interval
        public int tigL;

        public TIG(double k, double[] interval, int l){
            this.tigK = k;
            this.tigInterval = interval;
            this.tigL = l;
        }

        public static double DTW_Omega(List<TIG> x, List<TIG> y){
            int n = x.size();
            int m = y.size();

            double[][] dtw = new double[n+1][m+1];

            for (int i = 1; i <= n; i++)
                dtw[i][0] = Double.MAX_VALUE;
            for (int i = 1; i <= m; i++)
                dtw[0][i] = Double.MAX_VALUE;

            dtw[0][0] = 0;
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    double dG = Math.max(Math.abs(x.get(i-1).tigInterval[0] - y.get(j-1).tigInterval[0]),
                            Math.abs(x.get(i-1).tigInterval[1] - y.get(j-1).tigInterval[1]));
                    dtw[i][j] = dG + Math.min(Math.min(dtw[i - 1][j],
                                    dtw[i][j - 1]),
                            dtw[i - 1][j - 1]
                    );
                }
            }
            return dtw[n][m];
        }

        public static double DTW_Trend(List<TIG> x, List<TIG> y){
            int n = x.size();
            int m = y.size();

            double[][] dtw = new double[n+1][m+1];

            for (int i = 1; i <= n; i++)
                dtw[i][0] = Double.MAX_VALUE;
            for (int i = 1; i <= m; i++)
                dtw[0][i] = Double.MAX_VALUE;

            dtw[0][0] = 0;
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    double dT = (x.get(i-1).tigL + y.get(j-1).tigL) * 1.0 / 2 * Math.abs(x.get(i-1).tigK - y.get(j-1).tigK);
                    dtw[i][j] = dT + Math.min(Math.min(dtw[i - 1][j],
                                    dtw[i][j - 1]),
                            dtw[i - 1][j - 1]
                    );
                }
            }
            return dtw[n][m];
        }


    }
}
