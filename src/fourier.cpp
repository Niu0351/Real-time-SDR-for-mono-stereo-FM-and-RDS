/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "iofunc.h"
#include "fourier.h"
#include <vector>
#include <math.h>

template<typename T>
std::vector<T> arange(T start, T stop, T step = 1){
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

// just DFT function (no FFT yet)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
    Xf.resize(x.size(), static_cast<std::complex<float>>(0, 0));
    for (auto m = 0; m < Xf.size(); m++) {
        for (auto k = 0; k < x.size(); k++) {
            std::complex<float> expval(0, -2*PI*(k*m) / x.size());
            Xf[m] += x[k] * std::exp(expval);
        }
    }
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
    // only the positive frequencies
    Xmag.resize(Xf.size(), static_cast<float>(0));
    for (auto i = 0; i < Xf.size(); i++) {
        Xmag[i] = std::abs(Xf[i])/Xf.size();
    }
}

// add your own code to estimate the PSD
void estimatePSD(const std::vector<float> &samples, int no_FFT, float Fs, std::vector<float> &freq, std::vector<float> &psd_est)
{
    int freq_bins = no_FFT;
    float df = Fs/freq_bins;
    freq = arange<float> (0, Fs/2, df);

    std::vector<float> hann;
    hann.resize(freq_bins, 0.0);
    for (auto i=0; i<int(hann.size()); i += 1){
        hann[i] = pow(sin(i*PI/freq_bins),2);
    }

    int no_segments = int(floor(samples.size()/float(freq_bins)));

    std::vector<float> psd_list;
    psd_list.resize(no_segments*freq_bins, 0.0);
    int list_index = 0;

    for (auto k=0; k < no_segments; k += 1){

        int no_win = 0;
        std::vector<float> windowed_samples;
        windowed_samples.resize(freq_bins, 0.0);
        std::vector<std::complex<float>> Xf;
        Xf.resize(no_FFT, static_cast<std::complex<float>>(0, 0));
        std::vector<std::complex<float>> Xf_pos;
        Xf_pos.resize(int(Xf.size()/2), static_cast<std::complex<float>>(0, 0));
        std::vector<float> psd_seg;
        psd_seg.resize(Xf_pos.size(), 0.0);

        for (auto m=k*freq_bins; m < (k+1)*freq_bins; m +=1){
            windowed_samples[no_win] = samples[m] * hann[no_win];
            no_win += 1;
        }

        DFT(windowed_samples, Xf);

        for (auto m=0; m < int(no_FFT/2); m += 1){
            Xf_pos[m] = Xf[m];
            psd_seg[m] = 2*((1/(Fs*freq_bins/2)) * pow(abs(Xf_pos[m]),2));
        }

        for (auto m=0; m < int(psd_seg.size()); m += 1) {
            psd_seg[m] = 10 * (log10(psd_seg[m]));
            psd_list[list_index] = psd_seg[m];
            list_index += 1;
        }

    }

    psd_est.resize(int(freq_bins/2), 0.0);

    for (auto m=0; m < int(freq_bins/2); m += 1){
        for (auto n=0; n < no_segments; n += 1){
            psd_est[m] += psd_list[m + n*int(freq_bins/2)];
        }
        psd_est[m] /= no_segments;

    }

}
