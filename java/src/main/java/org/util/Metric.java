package org.util;

public class Metric {
    public static double calculateRandIndex(int[] clustering1, int[] clustering2) {
        if (clustering1.length != clustering2.length) {
            throw new IllegalArgumentException();
        }

        int n = clustering1.length;
        int a = 0;
        int b = 0;

        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                boolean sameClusterInClustering1 = clustering1[i] == clustering1[j];
                boolean sameClusterInClustering2 = clustering2[i] == clustering2[j];

                if (sameClusterInClustering1 && sameClusterInClustering2) {
                    a++;
                } else if (!sameClusterInClustering1 && !sameClusterInClustering2) {
                    b++;
                }
            }
        }

        int totalPairs = n * (n - 1) / 2;
        return (double) (a + b) / totalPairs;
    }
}
