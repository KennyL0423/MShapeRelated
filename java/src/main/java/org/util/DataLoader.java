package org.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class DataLoader {
    public static List<double[]> readTimeSeriesFromCSV(String filePath, int seqLen) {
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
}
