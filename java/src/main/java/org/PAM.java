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

    public void fit(){
        pamClustering(data, clusterNum, labels);
    }

    private void pamClustering(List<double[]> sequences, int k, int[] assignment) {
        int n = sequences.size();
        System.out.println(n);
        ArrayList<Integer> medoidIndexes = selectInitialMedoids(sequences, n, k);

        boolean changed = true;
        int maxIter = 6, iter = 0;
        while (changed && iter < maxIter) {
            changed = false;
            assignSeq(sequences, medoidIndexes, assignment);

            ArrayList<Integer> newMedoidIndices = new ArrayList<>();
            for (int i = 0; i < k; i++) {
                int newMedoidIndex = findNewMedoid(sequences, medoidIndexes.get(i), assignment, i);
                newMedoidIndices.add(newMedoidIndex);
                if (newMedoidIndex != medoidIndexes.get(i)) {
                    changed = true;
                }
            }

            if (changed) {
                medoidIndexes = newMedoidIndices;
            }
            iter += 1;
            System.out.println(iter);
        }
        int [] indexes = new int[medoidIndexes.size()];
        for (int i = 0; i < medoidIndexes.size(); i++) {
            indexes[i] = medoidIndexes.get(i);
        }
        this.labels = indexes;
    }

    private ArrayList<Integer> selectInitialMedoids(List<double[]> sequences, int n, int k) {
        ArrayList<Integer> indexes = new ArrayList<>();
        Random rand = new Random();
        int firstMedoid = rand.nextInt(n);
        indexes.add(firstMedoid);

        for (int i = 1; i < k; i++) {
            int maxDistance = 0;
            int maxDistanceIndex = 0;
            for (int j = 0; j < n; j++) {
                double minDistance = Double.MAX_VALUE;
                for (int index : indexes) {
                    double distance = computeDistance(sequences.get(j), sequences.get(index));
                    if (distance < minDistance) {
                        minDistance = distance;
                    }
                }
                if (minDistance > maxDistance) {
                    maxDistance = (int) minDistance;
                    maxDistanceIndex = j;
                }
            }
            indexes.add(maxDistanceIndex);
        }
        return indexes;
    }

    private void assignSeq(List<double[]> sequences, ArrayList<Integer> indexes, int[] assignment) {
        for (int i = 0; i < sequences.size(); i++) {
            double minDistance = Double.MAX_VALUE;
            int closestMedoidIndex = -1;
            for (int j = 0; j < indexes.size(); j++) {
                double distance = computeDistance(sequences.get(i), sequences.get(indexes.get(j)));
                if (distance < minDistance) {
                    minDistance = distance;
                    closestMedoidIndex = j;
                }
            }
            assignment[i] = closestMedoidIndex;
        }
    }

    private int findNewMedoid(List<double[]> sequences, int oldMedoidIndex, int[] assignment, int clusterIndex) {
        double minTotalDistance = Double.MAX_VALUE;
        int newMedoidIndex = oldMedoidIndex;
        for (int i = 0; i < sequences.size(); i++) {
            if (assignment[i] == clusterIndex) {
                double totalDistance = 0;
                for (int j = 0; j < sequences.size(); j++) {
                    if (assignment[j] == clusterIndex) {
                        totalDistance += computeDistance(sequences.get(i), sequences.get(j));
                    }
                }
                if (totalDistance < minTotalDistance) {
                    minTotalDistance = totalDistance;
                    newMedoidIndex = i;
                }
            }
        }
        return newMedoidIndex;
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