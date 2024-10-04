package org.util;

import java.util.Arrays;

public class FFT {
  public static Complex[] frft(double[] x, double alpha) {
    int n = x.length;
    int fft_size = (int) Math.pow(2, Integer.toBinaryString(2 * n - 1).length());
    Complex[] fftResult = fft(x, fft_size);

    // Modify FFT coefficients by fractional power of alpha
    for (int k = 0; k < fft_size; k++) {
      double angle = alpha * Math.PI * k * k / fft_size;
      Complex fractionalExponent = new Complex(Math.cos(angle), Math.sin(angle));
      fftResult[k] = fftResult[k].multiply(fractionalExponent);
    }
    return fftResult;
  }

  public static Complex[] fft(double[] a, int length) {
    Complex[] b = new Complex[a.length];
    for (int i = 0; i < a.length; i++) b[i] = new Complex(a[i], 0);
    return fft(b, length);
  }

  // FFT with padding or croping
  public static Complex[] fft(Complex[] a, int length) {
    if (length < a.length) return fft(Arrays.copyOf(a, length));
    if (length > a.length) {
      Complex[] b = new Complex[length];
      for (int i = 0; i < a.length; i++) b[i] = a[i];
      for (int i = a.length; i < length; i++) b[i] = new Complex(0, 0);
      return fft(b);
    }
    return fft(a);
  }

  private static Complex[] fft(Complex[] a) {
    int n = a.length;

    // Base case
    if (n == 1) return new Complex[] {a[0]};

    // Split the array into even and odd parts
    Complex[] even = new Complex[n / 2];
    Complex[] odd = new Complex[n / 2];
    for (int i = 0; i < n / 2; i++) {
      even[i] = a[i * 2];
      odd[i] = a[i * 2 + 1];
    }

    // Recursive FFT
    Complex[] evenResult = fft(even);
    Complex[] oddResult = fft(odd);

    // Combine
    Complex[] result = new Complex[n];
    for (int i = 0; i < n / 2; i++) {
      double angle = -2 * i * Math.PI / n;
      Complex t = new Complex(Math.cos(angle), Math.sin(angle)).multiply(oddResult[i]);
      result[i] = evenResult[i].add(t);
      result[i + n / 2] = evenResult[i].subtract(t);
    }

    return result;
  }

  public static double[] ifft(Complex[] a) {
    int n = a.length;

    // Take conjugate of the input
    for (int i = 0; i < n; i++) {
      a[i] = a[i].conjugate();
    }

    // Compute FFT using conjugated input
    Complex[] tmp = fft(a);
    double[] result = new double[n];

    // Take conjugate of the output and scale
    for (int i = 0; i < n; i++) {
      result[i] = tmp[i].conjugate().scale(1.0 / n).getReal();
    }

    return result;
  }

  public static double[] ncc(double[] x1, double[] x2) {
    double den = norm(x1) * norm(x2);
    if (den < 1e-9) den = Double.MAX_VALUE;
    int x_len = x1.length; // l
    int fft_size = (int) Math.pow(2, Integer.toBinaryString(2 * x_len - 1).length());
    double[] cc =
            FFT.ifft(Complex.multiply(FFT.fft(x1, fft_size), Complex.conjugate(FFT.fft(x2, fft_size))));
    double[] ncc = new double[fft_size - 1];
    // [-(x_len-1):] + [:x_len]
    for (int i = 0; i < fft_size - 1; i++)
      if (i < x_len - 1) ncc[i] = cc[cc.length - x_len + 1 + i] / den;
      else ncc[i] = cc[i - x_len + 1] / den;
    return ncc;
  }

  public static double norm(double[] x) {
    double res = 0.0;
    for (double v : x) {
      res += v * v;
    }
    return Math.sqrt(res);
  }





  public static void main(String[] args) {
    Complex[] a = {new Complex(1, 0), new Complex(2, 0), new Complex(3, 0), new Complex(4, 0)};

    System.out.println("Input: " + Arrays.toString(a));

    Complex[] fftResult = fft(new double[] {1, 2, 3, 4}, 4);
    System.out.println("FFT: " + Arrays.toString(fftResult));

    double[] ifftResult = ifft(fftResult);
    System.out.println("IFFT: " + Arrays.toString(ifftResult));
  }
}
