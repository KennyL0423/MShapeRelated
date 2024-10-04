package org.example;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class FrOKShape {
    private static final int MAX_ITERATIONS = 100;
    static int seqLen = 166;
    static int clusterNum = 3;

    public static void main(String[] args) {

        long start = System.currentTimeMillis();
//        String csvFile = "/Users/suyx1999/Downloads/jinfeng.csv";
        String csvFile = "/Users/suyx1999/ExpData/shape/air.csv";
        List<double[]> timeSeriesData = readTimeSeriesFromCSV(csvFile);
        int numSeries = timeSeriesData.size();

        // init labels
        int[] labels = new int[numSeries];
        int[] clusterLabels = FrOKShapeClustering(timeSeriesData, clusterNum);

        long end = System.currentTimeMillis();
        System.out.println("Time taken: " + (end - start) + "ms");

//        for (int i = 0; i < clusterLabels.length; i++) {
//            System.out.println("Time series " + i + " is assigned to cluster " + clusterLabels[i]);
//        }
    }

    private static List<double[]> readTimeSeriesFromCSV(String filePath) {
        List<double[]> timeSeriesData = new ArrayList<>();
        double[] seq = new double[seqLen];

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            boolean isFirstLine = true;

            int cnt = 0;
            while ((line = br.readLine()) != null) {
                if (isFirstLine) {
                    isFirstLine = false;
                    continue;
                }
                String[] values = line.split(",");

                if (values.length < 2) seq[cnt] = 0.0;
                else seq[cnt] = Double.parseDouble(values[1]);
                cnt += 1;
                if (cnt == seqLen) {
                    timeSeriesData.add(seq.clone());
                    seq = new double[seqLen];
                    cnt = 0;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return timeSeriesData;
    }

    // FrOKShape Clustering Method
    public static int[] FrOKShapeClustering(List<double[]> X, int K) {
        int n = X.size();
        int l = X.get(0).length;
        int[] labels = new int[n]; // Cluster labels

        Random random = new Random();
        for (int i = 0; i < n; i++) {
            labels[i] = random.nextInt(K);
        }

        // Cluster centers
        double[][] kCenter = new double[K][l];

        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            // Store previous labels
            int[] previousLabels = Arrays.copyOf(labels, labels.length);

            // Update cluster centers using DeterminClusteringCenter (DBA-based)
            for (int k = 0; k < K; k++) {
                kCenter[k] = DeterminClusteringCenter(X, labels, k);
            }

            // Update labels for each time series based on FrOSBD distance
            for (int i = 0; i < n; i++) {
                double minDist = Double.MAX_VALUE;
                int bestCluster = -1;
                for (int k = 0; k < K; k++) {
                    double dist = FrOSBD(X.get(i), kCenter[k]);
                    if (dist < minDist) {
                        minDist = dist;
                        bestCluster = k;
                    }
                }
                labels[i] = bestCluster;
            }
            if (Arrays.equals(previousLabels, labels)) {
                System.out.println("Converged after " + iter + " iterations.");
                break;
            }
        }
        return labels;
    }

    // Calculate the Determining Clustering Center using DBA
    private static double[] DeterminClusteringCenter(List<double[]> X, int[] labels, int k) {
        int n = X.size();
        int l = X.get(0).length;
        double[] center = new double[l];

        int count = 0;
        for (int i = 0; i < n; i++) {
            if (labels[i] == k) {
                for (int j = 0; j < l; j++) {
                    center[j] += X.get(i)[j];
                }
                count++;
            }
        }

        if (count > 0) {
            for (int j = 0; j < l; j++) {
                center[j] /= count;
            }
        }

        return center;
    }

    // Fractional Order Shape-Based Distance (FrOSBD)
    private static double FrOSBD(double[] x, double[] y) {
        // Assuming p is the fractional order (can be parameterized)
        double maxFrONCC = 0;
        int l = x.length;

        for (int shift = -l + 1; shift < l; shift++) {
            double FrONCC = fractionalCorrelation(x, y, shift);
            maxFrONCC = Math.max(maxFrONCC, FrONCC);
        }

        return 1.0 - maxFrONCC;
    }

    // Fractional Correlation Calculation (FrOCC)
    private static double fractionalCorrelation(double[] x, double[] y, int shift) {
        double sum = 0.0;

        for (int i = 0; i < x.length; i++) {
            int j = (i + shift) % y.length;
            if (j < 0) j += y.length;
            sum += x[i] * y[j];
        }

        return sum / (norm(x) * norm(y));
    }

    // Function to calculate the norm of a time series
    private static double norm(double[] series) {
        double sum = 0.0;
        for (double v : series) {
            sum += v * v;
        }
        return Math.sqrt(sum);
    }

    private static double norm(List<Double> series) {
        double sum = 0.0;
        for (double v : series) {
            sum += v * v;
        }
        return Math.sqrt(sum);
    }
}

