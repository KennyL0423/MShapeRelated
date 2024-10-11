package org;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class PAM {
    List<double[]> data;
    int seqLen;
    int clusterNum;
    int[] labels;

    public PAM(List<double[]> data, int seqLen, int clusterNum){
        this.data = data;
        this.seqLen = seqLen;
        this.clusterNum = clusterNum;
        this.labels = new int[data.size()];
    }

    public int[] fit(){
        pamClustering(data, clusterNum);
        return labels;
    }

    private void pamClustering(List<double[]> sequences, int k) {
        int n = sequences.size();
        System.out.println(n);
        List<Integer> medoids = initializeMedoids(data.size(), k);
        boolean changed = true;
        int maxIter = 100, iter = 0;

        while (changed && iter < maxIter) {
            int[] clusters = assignClusters(medoids);
            List<Integer> newMedoids = updateMedoids(clusters, medoids);
            if (newMedoids.equals(medoids)) {
                changed = false;
                medoids = newMedoids;
            } else {
                medoids = newMedoids;
            }
            iter += 1;
            System.out.println(iter);
        }


        labels = assignClusters(medoids);


    }

    private List<Integer> initializeMedoids(int numPoints, int k) {
        Random rand = new Random();
        List<Integer> medoids = new ArrayList<>();
        while (medoids.size() < k) {
            int medoid = rand.nextInt(numPoints);
            if (!medoids.contains(medoid)) {
                medoids.add(medoid);
            }
        }
        return medoids;
    }

    // 将点分配到最近的 Medoid
    private int[] assignClusters(List<Integer> medoids) {
        int[] clusters = new int[data.size()];
        for (int i = 0; i < data.size(); i++) {
            double minDist = Double.MAX_VALUE;
            int closestMedoid = -1;
            for (int medoid : medoids) {
                double dist = computeDistance(data.get(i), data.get(medoid));
                if (dist < minDist) {
                    minDist = dist;
                    closestMedoid = medoid;
                }
            }
            clusters[i] = closestMedoid;
        }
        return clusters;
    }

    private double totalCost(int[] clusters) {
        double cost = 0.0;
        for (int i = 0; i < data.size(); i++) {
            cost += computeDistance(data.get(i), data.get(clusters[i]));
        }
        return cost;
    }

    private List<Integer> updateMedoids(int[] clusters, List<Integer> medoids) {
        double bestCost = totalCost(clusters);
        for (int i = 0; i < medoids.size(); i++) {
            int currentMedoid = medoids.get(i);
            for (int j = 0; j < data.size(); j++) {
                if (!medoids.contains(j)) {
                    medoids.set(i, j);
                    int[] newClusters = assignClusters(medoids);
                    double newCost = totalCost(newClusters);
                    if (newCost < bestCost) {
                        bestCost = newCost;
                        currentMedoid = j;
                    }
                    medoids.set(i, currentMedoid);  // 恢复当前 medoid
                }
            }
        }
        return medoids;
    }

    // You need to implement this function to compute the distance between two sequences
    private double computeDistance(double[] seq1, double[] seq2) {
//        return msm(seq1, seq2, 1);
        return dtw(seq1, seq2);
    }


    public double msm(double[] a, double[] b, double c) {
        int m = a.length;
        double[][] CM = new  double[m + 1][m + 1];

        // Initialize base cases
        CM[1][1] = Math.abs(a[0] - b[0]);

        // Fill the cost matrix
        for (int i = 2; i <= m; i++) {
            CM[i][1] = CM[i - 1][1] + C(a[i - 1], a[i - 2], b[0], c);
        }

        for (int j = 2; j <= m; j++) {
            CM[1][j] = CM[1][j - 1] + C(b[j - 1], a[0], b[j - 2], c);
        }

        for (int i = 2; i <= m; i++) {
            for (int j = 2; j <= m; j++) {
                double move = CM[i - 1][j - 1] + Math.abs(a[i - 1] - b[j - 1]);
                double split = CM[i - 1][j] + C(a[i - 1], a[i - 2], b[j - 1], c);
                double merge = CM[i][j - 1] + C(b[j - 1], b[j - 2], a[i - 1], c);
                CM[i][j] = Math.min(Math.min(move, split), merge);
            }
        }

        return CM[m][m];
    }

    public double C(double x, double y, double z, double c) {
        if (y <= x && x <= z || y >= x && x >= z) {
            return c;
        } else {
            return c + Math.min(Math.abs(x - y), Math.abs(x - z));
        }
    }

    public double dtw(double[] s, double[] t) {
        int n = s.length;
        int m = t.length;
        double[][] dtw = new double[n + 1][m + 1];

        // Initialize the first row and column
        for (int i = 0; i <= n; i++) {
            dtw[i][0] = Double.POSITIVE_INFINITY;
        }
        for (int j = 0; j <= m; j++) {
            dtw[0][j] = Double.POSITIVE_INFINITY;
        }
        dtw[0][0] = 0;

        // Compute the DTW distance
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                double cost = Math.abs(s[i - 1] - t[j - 1]);
                dtw[i][j] = cost + Math.min(Math.min(dtw[i - 1][j], dtw[i][j - 1]), dtw[i - 1][j - 1]);
            }
        }

        return dtw[n][m];
    }
}