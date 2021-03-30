/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FOURIER_H
#define DY4_FOURIER_H

// add headers as needed
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

#include <vector>
#include <math.h>

//template<typename T>
//std::vector<T> arange(T start, T stop, T step = 1){
//    std::vector<T> values;
//    for (T value = start; value < stop; value += step)
//        values.push_back(value);
//    return values;
//}

// declaration of a function prototypes
void DFT(const std::vector<float> &,
         std::vector<std::complex<float>> &);

// you should add your own IDFT
// time-permitting you can build your own function for FFT

void computeVectorMagnitude(const std::vector<std::complex<float>> &,
                            std::vector<float> &);

// provide the prototype to estimate PSD
// ...
void estimatePSD(const std::vector<float> &, int, float, std::vector<float> &, std::vector<float> &);

#endif // DY4_FOURIER_H
