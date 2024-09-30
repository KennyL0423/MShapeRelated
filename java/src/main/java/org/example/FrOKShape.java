package org.example;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

public class FrOKShape {
    private static final int MAX_ITERATIONS = 100;

    public static void main(String[] args) {

        String csvFile = "jinfeng.csv";
        double[][] timeSeriesData = readTimeSeriesFromCSV(csvFile);
        int numClusters = 2; // cluster num
        int numSeries = timeSeriesData.length;

        // init labels
        int[] labels = new int[numSeries];
        int[] clusterLabels = FrOKShapeClustering(timeSeriesData, numClusters);

        for (int i = 0; i < clusterLabels.length; i++) {
            System.out.println("Time series " + i + " is assigned to cluster " + clusterLabels[i]);
        }
    }

    private static double[][] readTimeSeriesFromCSV(String filePath) {
        double[][] timeSeriesData = null;

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            boolean isFirstLine = true;

            int seriesCount = 0;
            while ((line = br.readLine()) != null) {
                if (isFirstLine) {
                    isFirstLine = false;
                    continue;
                }
                String[] values = line.split(",");
                double value = Double.parseDouble(values[1]);
                if (timeSeriesData == null) {
                    timeSeriesData = new double[1000][1];
                }
                timeSeriesData[seriesCount][0] = value;
                seriesCount++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return timeSeriesData;
    }

    // FrOKShape Clustering Method
    public static int[] FrOKShapeClustering(double[][] X, int K) {
        int numSeries = X.length;
        int[] labels = new int[numSeries]; // Cluster labels

        Random random = new Random();
        for (int i = 0; i < numSeries; i++) {
            labels[i] = random.nextInt(K);
        }

        // Cluster centers
        double[][] kCenter = new double[K][X[0].length];

        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            // Store previous labels
            int[] previousLabels = Arrays.copyOf(labels, labels.length);

            // Update cluster centers using DeterminClusteringCenter (DBA-based)
            for (int k = 0; k < K; k++) {
                kCenter[k] = DeterminClusteringCenter(X, labels, k);
            }

            // Update labels for each time series based on FrOSBD distance
            for (int i = 0; i < numSeries; i++) {
                double minDist = Double.MAX_VALUE;
                int bestCluster = -1;
                for (int k = 0; k < K; k++) {
                    double dist = FrOSBD(X[i], kCenter[k]);
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
    private static double[] DeterminClusteringCenter(double[][] X, int[] labels, int k) {
        int numSeries = X.length;
        int numPoints = X[0].length;
        double[] center = new double[numPoints];

        int count = 0;
        for (int i = 0; i < numSeries; i++) {
            if (labels[i] == k) {
                for (int j = 0; j < numPoints; j++) {
                    center[j] += X[i][j];
                }
                count++;
            }
        }

        if (count > 0) {
            for (int j = 0; j < numPoints; j++) {
                center[j] /= count;
            }
        }

        return center;
    }

    // Fractional Order Shape-Based Distance (FrOSBD)
    private static double FrOSBD(double[] x, double[] y) {
        // Assuming p is the fractional order (can be parameterized)
        double maxFrONCC = 0;

        for (int shift = -x.length + 1; shift < x.length; shift++) {
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
}

