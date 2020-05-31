#include <iostream>
#include <time.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <complex>
//#include "matplotlibcpp.h"

#define _USE_MATH_DEFINES
#define TIMEPOINTS 256
#define NUM 14
#define FREQU 2000
#define MAX 256

using namespace std;

double generatorSimple (int, float*, float*, float*);
void getSimpleSignal(complex<double>*, float*, float*, float*); //get pseudo random signal based on NUM=14 harmonic signals

int log2(int N);
int check(int n);
int reverse(int N, int n);
void ordina(complex<double>* f1, int N);
void transform(complex<double>* f, int N);
void FFT(complex<double>* f, int N, double d);

int main()
{
    int n = 256;      // Number of elements
    double d = 2000;   // Value of w
    complex<double> vec[n];

    float amplitudes[NUM] = {0};
    float phathe[NUM] = {0};
    float freq[NUM] = {0};

    srand (time(NULL));
    //generator of parametersof NUM harmonics
    for (int i = 0; i < NUM; i++){
      amplitudes[i] = rand();
      phathe[i] = rand();
      freq[i] = rand() % 2000 + 1;
    }

    getSimpleSignal(&vec[0], &amplitudes[0], &phathe[0], &freq[0]);

    FFT(vec, n, d);

    cout << "Printing the FFT of the array specified:" << endl;
    for(int j = 0; j < n; j++)
        cout << vec[j] << endl;

    return 0;
}

double generatorSimple (int t, float* amplitudes, float* ph, float* fr){
    double result = 0;
    for(int i = 0; i < NUM; i++)
        result += *(amplitudes + i) * sin(*(fr + i) * t + *(ph + i) * M_PI);
    return result;
}

void getSimpleSignal(complex<double>* pValueX, float* amplitudes, float* ph, float* fr){
  for (int i = 0; i < TIMEPOINTS; i++)
    *(pValueX + i) = generatorSimple(i, amplitudes, ph, fr);
}

// Function to calculate the log2(.) of int numbers
int log2(int N){
    int k = N, i = 0;
    while(k) {
        k >>= 1;
        i++;
    }
    return i - 1;
}

// Checking if the number of element is a power of 2
int check(int n){
    return n > 0 && (n & (n - 1)) == 0;
}

// Calculating reverse number
int reverse(int N, int n){
    int j, p = 0;

    for(j = 1; j <= log2(N); j++) {
        if(n & (1 << (log2(N) - j)))
            p |= 1 << (j - 1);
    }
    return p;
}

// Using the reverse order in the array
void ordina(complex<double>* f1, int N){
    complex<double> f2[MAX];

    for(int i = 0; i < N; i++)
        f2[i] = f1[reverse(N, i)];

    for(int j = 0; j < N; j++)
        f1[j] = f2[j];
}

// Calculating the fransformation
void transform(complex<double>* f, int N){
    ordina(f, N);    // 1.Reverse order
    complex<double> *W;
    int n = 1;
    int a = 0;

    W = (complex<double> *)malloc(N / 2 * sizeof(complex<double>));
    W[1] = polar(1., -2. * M_PI / N); // 2. Calculate W(n), it aslo can be: polar(1., 2. * M_PI / N)
    W[0] = 1;

    for(int i = 2; i < N / 2; i++)
        W[i] = pow(W[1], i);

    a = N / 2;

    for(int j = 0; j < log2(N); j++) {
        for(int i = 0; i < N; i++) {
            if(!(i & n)) {
                complex<double> temp = f[i];
                complex<double> Temp = W[(i * a) % (n * a)] * f[i + n];
                f[i] = temp + Temp;
                f[i + n] = temp - Temp;
            }
        }
        n *= 2;
        a = a / 2;
    }
}

// FFT function
void FFT(complex<double>* f, int N, double d){
    transform(f, N);
    for(int i = 0; i < N; i++)
        f[i] *= d; // Multiplying by step
}
