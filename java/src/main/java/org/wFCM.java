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
    int seqLen;
    int clusterNum;
    List<double[]> data;
    List<List<TIG>> tigsOfAllData;


    public static void main(String[] args) {
        long start = System.currentTimeMillis();
        String csvFile = "/Users/suyx1999/ExpData/shape/air.csv";
//        String csvFile = "/Users/suyx1999/Downloads/jinfeng.csv";

        List<double[]> timeSeriesData = DataLoader.readTimeSeriesFromCSV(csvFile, 166);
        wFCM clustering = new wFCM(timeSeriesData, 166, 3, 100);
        int[] clusterLabels = clustering.fit();

        long end = System.currentTimeMillis();
        System.out.println("Time taken: " + (end - start) + "ms");
    }

    public wFCM(List<double[]> data, int seqLen, int cluserNum, int max_iter) {
        this.data = data;
        this.seqLen = seqLen;
        this.clusterNum = cluserNum;
        this.MAX_ITERATIONS = max_iter;
        this.tigsOfAllData = new ArrayList<>();
    }

    public int[] fit(){
        int n = data.size();
        transformTIG();
        double[][] U = initializeMembershipMatrix(n, clusterNum);
        fuzzyClustering(data, clusterNum, U);
        int[] labels = new int[n];
        for (int i = 0; i < n; i++) {
            labels[i] = findMaxIndex(U[i]);
        }
        return labels;
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

        public static double DTW(List<TIG> x, List<TIG> y){
            int n = x.size();
            int m = y.size();

            // 创建一个 n+1 * m+1 的矩阵用于存储中间计算结果
            double[][] dtw = new double[n+1][m+1];

            // 初始化矩阵的第一个元素为无穷大
            for (int i = 1; i <= n; i++)
                dtw[i][0] = Double.MAX_VALUE;
            for (int i = 1; i <= m; i++)
                dtw[0][i] = Double.MAX_VALUE;

            dtw[0][0] = 0;
            // 计算 DTW 距离
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    double dG = Math.max(Math.abs(x.get(i-1).tigInterval[0] - y.get(j-1).tigInterval[0]),
                            Math.abs(x.get(i-1).tigInterval[1] - y.get(j-1).tigInterval[1]));
                    double dt = (x.get(i-1).tigL + y.get(j-1).tigL) * 1.0 / 2 * Math.abs(x.get(i-1).tigK - y.get(j-1).tigK);
                    dtw[i][j] = dG + dt + Math.min(Math.min(dtw[i - 1][j],
                                    dtw[i][j - 1]),
                            dtw[i - 1][j - 1]
                    );
                }
            }
            // 返回最终的 DTW 距离
            return dtw[n][m];
        }


    }

    private void transformTIG(){
        int l = data.get(0).length;
        int partNum = 4; // parameter q in the paper
        for (double[] seq: data){
            List<TIG> tigArray = new ArrayList<>();
            // divide seq into 4 parts with equal length and calculate the slope of each part
            for (int s = 0; s < partNum; s++){
                double[] part = Arrays.copyOfRange(seq, s * l / partNum, (s + 1) * l / partNum);
                LinearRegression lr = new LinearRegression(part);
                TIG tig = new TIG(lr.getSlope(), lr.findInfoGran(), part.length);
                tigArray.add(tig);
            }
            tigsOfAllData.add(tigArray);
        }
    }



    // Function to initialize the membership matrix randomly
    private double[][] initializeMembershipMatrix(int numSeries, int numClusters) {
        double[][] U = new double[numSeries][numClusters];
        Random random = new Random();
        for (int i = 0; i < numSeries; i++) {
            double sum = 0.0;
            for (int j = 0; j < numClusters; j++) {
                U[i][j] = random.nextDouble();
                sum += U[i][j];
            }
            // Normalize the membership values for each time series
            for (int j = 0; j < numClusters; j++) {
                U[i][j] /= sum;
            }
        }
        return U;
    }

    // Fuzzy clustering algorithm based on DTW
    private void fuzzyClustering(List<double[]> timeSeriesData, int numClusters, double[][] U) {
        int n = timeSeriesData.size();
        int l = timeSeriesData.get(0).length;

        // Cluster prototypes (initially set to random time series)
        double[][] clusterPrototypes = initializeClusterPrototypes(timeSeriesData, numClusters);

        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            // Step 1: Update the membership matrix U
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < numClusters; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < numClusters; k++) {
                        double dtw1 = DTW(timeSeriesData.get(i), clusterPrototypes[j]);
                        double dtw2 = DTW(timeSeriesData.get(i), clusterPrototypes[k]);
                        if (Math.abs(dtw1 - dtw2) < 1e-6) {
                            sum += 1.0;
                        } else {
                            sum += Math.pow(dtw1 / dtw2, 2 / (M - 1));
                        }
                    }
                    U[i][j] = 1.0 / sum;
                }
            }

            // Step 2: Update the cluster prototypes
            for (int j = 0; j < numClusters; j++) {
                for (int t = 0; t < l; t++) {
                    double numerator = 0.0;
                    double denominator = 0.0;
                    for (int i = 0; i < n; i++) {
                        double u_ij_m = Math.pow(U[i][j], M);
                        numerator += u_ij_m * timeSeriesData.get(i)[t];
                        denominator += u_ij_m;
                    }
                    clusterPrototypes[j][t] = numerator / denominator;
                }
            }

            // Step 3: Check for convergence
            double maxChange = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < numClusters; j++) {
                    maxChange = Math.max(maxChange, Math.abs(U[i][j] - clusterPrototypes[j][i]));
                }
            }
            if (maxChange < THRESHOLD) {
                System.out.println("Converged after " + iter + " iterations.");
                break;
            }
        }
    }

    // Function to calculate DTW distance between two time series
    private double DTW(double[] series1, double[] series2) {
        int n = series1.length;
        int m = series2.length;
        double[][] dtw = new double[n + 1][m + 1];

        // Initialize matrix with infinity
        for (double[] row : dtw) {
            Arrays.fill(row, Double.POSITIVE_INFINITY);
        }
        dtw[0][0] = 0.0;

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                double cost = Math.abs(series1[i - 1] - series2[j - 1]);
                dtw[i][j] = cost + Math.min(Math.min(dtw[i - 1][j], dtw[i][j - 1]), dtw[i - 1][j - 1]);
            }
        }

        return dtw[n][m];
    }

    // Function to find the index of the maximum value in an array
    private static int findMaxIndex(double[] array) {
        int maxIndex = 0;
        for (int i = 1; i < array.length; i++) {
            if (array[i] > array[maxIndex]) {
                maxIndex = i;
            }
        }
        return maxIndex;
    }

    // Function to initialize cluster prototypes by selecting random time series
    private static double[][] initializeClusterPrototypes(List<double[]> timeSeriesData, int numClusters) {
        double[][] prototypes = new double[numClusters][timeSeriesData.get(0).length];
        Random random = new Random();
        for (int i = 0; i < numClusters; i++) {
            int randomIndex = random.nextInt(timeSeriesData.size());
            System.arraycopy(timeSeriesData.get(randomIndex), 0, prototypes[i], 0, timeSeriesData.get(0).length);
        }
        return prototypes;
    }
}
