package org.util;

import java.io.*;
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

    public static void writeLabelsIntoCSV(String filePath, int[] labels) throws IOException{
        File file = new File(filePath);
        File father = file.getParentFile();
        if (!father.exists()) father.mkdirs();

        BufferedWriter bufferedWriter =  new BufferedWriter(new FileWriter(file, false));
        for (int label : labels) {
            bufferedWriter.write(label + "\n");
        }
        bufferedWriter.flush();
        bufferedWriter.close();
    }

    public static int[] readLabelsFromCSV(String filePath) {
        List<Integer> labels = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                labels.add(Integer.parseInt(line));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return labels.stream().mapToInt(i -> i).toArray();
    }

}
