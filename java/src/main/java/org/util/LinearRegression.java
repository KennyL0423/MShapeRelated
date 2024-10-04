package org.util;

import java.util.Arrays;

public class LinearRegression {

    double slope;
    double intercept;
    double[] residuals;

    double[] seq;



    public LinearRegression(double[] y) {
        seq = y;
        int n = y.length;
        double sumX = 0;
        double sumY = 0;
        double sumXY = 0;
        double sumX2 = 0;

        // x is from 1 to n
        for (int i = 0; i < n; i++) {
            double x = i + 1;
            sumX += x;
            sumY += y[i];
            sumXY += x * y[i];
            sumX2 += x * x;
        }

        // Calculate the slope (a) and intercept (b)
        double a = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
        double b = (sumY - a * sumX) / n;

        // Calculate residuals
        double[] res = new double[n];
        for (int i = 0; i < n; i++) {
            double x = i + 1;
            double predictedY = a * x + b;  // y = ax + b
            res[i] = y[i] - predictedY;  // Residual = actual y - predicted y
        }

        slope = a;
        intercept = b;
        residuals = res;
    }

    public double[] findInfoGran(){
        double a, b;
        double meanv = 0.0;
        for (int i = 0; i < residuals.length; i++) {
            meanv += residuals[i];
        }
        meanv /= residuals.length;
        double va = Double.MIN_VALUE;
        double vb = Double.MIN_VALUE;

        double[] sortedRes = residuals.clone();
        Arrays.sort(sortedRes);
        int meanIdx = -1;
        for (int i = 0; i < sortedRes.length - 1; i++) {
            if (sortedRes[i] <= meanv && meanv <= sortedRes[i + 1]){
                meanIdx = i;
                break;
            }
        }
        int nminus = meanIdx + 1;
        int nplus = sortedRes.length - meanIdx - 1;
        int rightIdx = meanIdx + 1;
        int leftIdx = meanIdx;
        int leftCard = 0, rightCard = 0;
        while (true){
            double tmp_va = leftCard * 1.0 / nminus * (1- (meanv - sortedRes[leftIdx])/ (sortedRes[rightIdx] - sortedRes[leftIdx]));
            double tmp_vb = rightCard * 1.0 / nplus * (1- (sortedRes[rightIdx] - meanv)/ (sortedRes[rightIdx] - sortedRes[leftIdx]));
            if (tmp_va > va){
                va = tmp_va;
                leftCard++;
                leftIdx--;
            }
            if (tmp_vb > vb){
                vb = tmp_vb;
                rightCard++;
                rightIdx++;
            }
            if (tmp_va <= va && tmp_vb <= vb){
                break;
            }
        }
        return new double[]{va, vb};
    }

    public double getSlope() {
        return slope;
    }

    public double getIntercept() {
        return intercept;
    }

    public double[] getResiduals() {
        return residuals;
    }
}
