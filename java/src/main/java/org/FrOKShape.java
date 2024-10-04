package org;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.util.Complex;
import org.util.DataLoader;
import org.util.FFT;

public class FrOKShape {
    int MAX_ITERATIONS;
    int seqLen;
    int clusterNum;
    double alpha;
    List<double[]> data;

    public static void main(String[] args) {
        long start = System.currentTimeMillis();
//        String csvFile = "/Users/suyx1999/Downloads/jinfeng.csv";
        String csvFile = "/Users/suyx1999/ExpData/shape/air.csv";
        List<double[]> timeSeriesData = DataLoader.readTimeSeriesFromCSV(csvFile, 166);

        FrOKShape clustering = new FrOKShape(timeSeriesData,  166, 3, 0.6, 100);
        int[] clusterLabels = clustering.fit();

        long end = System.currentTimeMillis();
        System.out.println("Time taken: " + (end - start) + "ms");

//        for (int i = 0; i < clusterLabels.length; i++) {
//            System.out.println("Time series " + i + " is assigned to cluster " + clusterLabels[i]);
//        }
    }

    public FrOKShape(List<double[]> data, int seqLen, int cluserNum, double alpha, int max_iter) {
        this.data = data;
        this.seqLen = seqLen;
        this.clusterNum = cluserNum;
        this.alpha = alpha;
        this.MAX_ITERATIONS = max_iter;
    }

    public int[] fit(){
        return FrOKShapeClustering(data, clusterNum);
    }


    // FrOKShape Clustering Method
    public int[] FrOKShapeClustering(List<double[]> X, int K) {
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
            System.out.println("Iteration " + iter + " completed.");
        }
        return labels;
    }

    // Calculate the Determining Clustering Center using DBA
    private double[] DeterminClusteringCenter(List<double[]> X, int[] labels, int k) {
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
    private double FrOSBD(double[] x, double[] y) {
        double den = FFT.norm(x) * FFT.norm(y);
        if (den < 1e-9) den = Double.MAX_VALUE;
        int l = x.length; // l
        int fft_size = (int) Math.pow(2, Integer.toBinaryString(2 * l - 1).length());

        Complex[] X_alpha = FFT.frft(x, alpha);
        Complex[] Y_alpha = FFT.frft(y, alpha);

        double[] cc =
                FFT.ifft(Complex.multiply(X_alpha, Complex.conjugate(Y_alpha)));
        double[] ncc = new double[fft_size - 1];
        for (int i = 0; i < fft_size - 1; i++)
            if (i < l - 1) ncc[i] = cc[cc.length - l + 1 + i] / den;
            else ncc[i] = cc[i - l + 1] / den;
        double maxFroncc = Arrays.stream(ncc).max().getAsDouble();
        return 1.0 - maxFroncc;
    }

}

