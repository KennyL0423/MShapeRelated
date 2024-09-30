package org.example;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class wFCM {
    private static final int MAX_ITERATIONS = 100;
    private static final double THRESHOLD = 1e-4;
    private static final int M = 2; // Fuzzification coefficient

    public static void main(String[] args) {
        String csvFile = "path_to_your_file.csv";

        // Read time series data from the CSV file
        double[][] timeSeriesData = readTimeSeriesFromCSV(csvFile);

        int numClusters = 2; // Number of clusters
        int numSeries = timeSeriesData.length; // Number of time series

        // Initialize membership matrix U
        double[][] U = initializeMembershipMatrix(numSeries, numClusters);

        // Iterate and update U and cluster prototypes
        fuzzyClustering(timeSeriesData, numClusters, U);

        // Final cluster assignment (for each time series, the highest membership value determines its cluster)
        for (int i = 0; i < numSeries; i++) {
            int assignedCluster = findMaxIndex(U[i]);
            System.out.println("Time series " + i + " belongs to cluster " + assignedCluster);
        }
    }
    private static double[][] readTimeSeriesFromCSV(String filePath) {
        List<List<Double>> timeSeriesList = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            boolean isFirstLine = true;
            while ((line = br.readLine()) != null) {
                if (isFirstLine) {
                    isFirstLine = false; // Skip the header
                    continue;
                }
                String[] values = line.split(",");
                if (values.length > 1) {
                    double value = Double.parseDouble(values[1]);
                    if (timeSeriesList.isEmpty()) {
                        timeSeriesList.add(new ArrayList<>());
                    }
                    timeSeriesList.get(0).add(value);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        double[][] timeSeriesData = new double[timeSeriesList.size()][];
        for (int i = 0; i < timeSeriesList.size(); i++) {
            List<Double> series = timeSeriesList.get(i);
            timeSeriesData[i] = series.stream().mapToDouble(Double::doubleValue).toArray();
        }
        return timeSeriesData;
    }

    // Function to initialize the membership matrix randomly
    private static double[][] initializeMembershipMatrix(int numSeries, int numClusters) {
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
    private static void fuzzyClustering(double[][] timeSeriesData, int numClusters, double[][] U) {
        int numSeries = timeSeriesData.length;
        int length = timeSeriesData[0].length;

        // Cluster prototypes (initially set to random time series)
        double[][] clusterPrototypes = initializeClusterPrototypes(timeSeriesData, numClusters);

        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            // Step 1: Update the membership matrix U
            for (int i = 0; i < numSeries; i++) {
                for (int j = 0; j < numClusters; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < numClusters; k++) {
                        double ratio = DTW(timeSeriesData[i], clusterPrototypes[j]) /
                                DTW(timeSeriesData[i], clusterPrototypes[k]);
                        sum += Math.pow(ratio, 2 / (M - 1));
                    }
                    U[i][j] = 1.0 / sum;
                }
            }

            // Step 2: Update the cluster prototypes
            for (int j = 0; j < numClusters; j++) {
                for (int t = 0; t < length; t++) {
                    double numerator = 0.0;
                    double denominator = 0.0;
                    for (int i = 0; i < numSeries; i++) {
                        double u_ij_m = Math.pow(U[i][j], M);
                        numerator += u_ij_m * timeSeriesData[i][t];
                        denominator += u_ij_m;
                    }
                    clusterPrototypes[j][t] = numerator / denominator;
                }
            }

            // Step 3: Check for convergence
            double maxChange = 0.0;
            for (int i = 0; i < numSeries; i++) {
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
    private static double DTW(double[] series1, double[] series2) {
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
    private static double[][] initializeClusterPrototypes(double[][] timeSeriesData, int numClusters) {
        double[][] prototypes = new double[numClusters][timeSeriesData[0].length];
        Random random = new Random();
        for (int i = 0; i < numClusters; i++) {
            int randomIndex = random.nextInt(timeSeriesData.length);
            System.arraycopy(timeSeriesData[randomIndex], 0, prototypes[i], 0, timeSeriesData[0].length);
        }
        return prototypes;
    }
}
