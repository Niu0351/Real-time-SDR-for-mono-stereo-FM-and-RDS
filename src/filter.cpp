/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>
#include <chrono>

#define PI 3.14159265358979323846

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
    // allocate memory for the impulse response
    h.resize(num_taps, 0.0);

    // the rest of the code in this function is to be completed by you
    // based on your understanding and the Python code from the first lab
    float norm_cutoff = Fc/(Fs/2);

    for (auto i = 0; i < num_taps; i += 1){
        if (i == ((num_taps-1) / 2)){
            h[i] = norm_cutoff;
        }
        else{
            h[i] = norm_cutoff * ((sin(PI*norm_cutoff*(i-(num_taps-1)/2))) / (PI*norm_cutoff*(i-(num_taps-1)/2)));
        }
        h[i] *= pow((sin(i*PI/num_taps)), 2);
    }
}

void impulseResponseRootRaiseCosine(float Fs, int N_taps, std::vector<float> &impulseResponseRRC){

    //duration for each symbol
    float T_symbol = 1/2375.0;

    //roll-off factor, 0 < beta < 1
    float beta = 0.90;

    for (auto i=0;i<N_taps;i++){
        float t = (i-N_taps/2)/Fs;
        //ignore the 1/T_symbol scale factor
        if (t == 0.0) {
            impulseResponseRRC[i] = 1.0 + beta*(4/PI - 1);
        }
        else if (t == (-T_symbol/(4*beta)) || t == T_symbol/(4*beta)){
            impulseResponseRRC[i] = (beta/sqrt(2)) * (((1+2/PI)*sin(PI/(4*beta))) + ((1-2/PI) * (cos(PI/(4*beta)))));
        }
        else {
            impulseResponseRRC[i] = (sin(PI*t*(1-beta)/T_symbol) + (4*beta*(t/T_symbol)*cos(PI*t*(1+beta)/T_symbol)))/(PI*t*(1-(4*beta*t/T_symbol) * (4*beta*t/T_symbol))/T_symbol);
        }
        //std::cerr << i << " " << impulseResponseRRC[i] << std::endl;
    }

}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
    // allocate memory for the output (filtered) data
    y.resize(x.size()+h.size()-1, 0.0);

    // the rest of the code in this function is to be completed by you
    // based on your understanding and the Python code from the first lab
    for (auto n = 0;n < y.size(); n += 1){
        int n1 = n;
        for (auto m = 0; m < x.size(); m += 1){

            if (n1>=0 && n1<h.size()){
                y[n] += x[m] * h[n1];
            }
            n1 -= 1;
        }
    }
}

void blockConv(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
    int cycle = 0;
    int y_count = 0;

    for (auto n = 0; (n+10) < block_size; n += 10){
        int buf = 0;
        y[y_count] = 0;

        for (auto k = 0; k < h.size(); k += 1){
            if ((cycle+buf) < (h.size()-1)){
                y[y_count] += buf_wave[cycle+buf] * h[h.size()-k-1];
                buf += 1;
            }
            else{
                y[y_count] += x[n+k-h.size()+1] * h[h.size()-k-1];
            }
        }
        cycle += 10;
        y_count++;
    }
}

void blockConv_I(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
    int cycle = 0;
    int y_count = 0;

    for (auto n = 0; n < block_size; n += 20){
        int buf = 0;
        y[y_count] = 0;
        int x_count = 0;

        for (auto k = 0; k < h.size(); k += 1){
            if ((cycle+buf) < (h.size()-1)){
                y[y_count] += buf_wave[cycle+buf] * h[h.size()-k-1];
                buf += 1;
            }
            else{
                if (n < 202){
                    y[y_count] += x[n+k-h.size()+1-y_count*10+x_count] * h[h.size()-k-1];
                }
                else{
                    y[y_count] += x[n+x_count-h.size()+1-h.size()+1+x_count] * h[h.size()-k-1];
                }
                x_count++;
            }
        }
        cycle += 10;
        //std::cerr << y_count << " " << y[y_count] << std::endl;
        y_count++;
    }
}

// 1024*5*10/2
void blockConv_Q(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
    int cycle = 0;
    int y_count = 0;
    int base = 0;

    for (auto n = 1; n < block_size; n += 20){
        int buf = 0;
        y[y_count] = 0;
        int x_count = 0;

        for (auto k = 0; k < h.size(); k += 1){
            if ((cycle+buf) < (h.size()-1)){
                y[y_count] += buf_wave[cycle+buf] * h[h.size()-k-1];
                buf += 1;
            }
            else{
                if (n < 203){
                    y[y_count] += x[n+k-h.size()+1-y_count*10+x_count] * h[h.size()-k-1];
                }
                else {
                    y[y_count] += x[n+x_count-h.size()+1-h.size()+1+x_count] * h[h.size()-k-1];
                }
                x_count++;
            }
        }
        cycle += 10;
        y_count++;
    }
}

void blockConv1(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
    int cycle = 0;
    int y_count = 0;

    for (auto n = 0; (n+5) < block_size; n += 5){
        int buf = 0;

        for (auto k = 0; k < h.size(); k += 1){
            if ((cycle+buf) < (h.size()-1)){
                y[y_count] += buf_wave[cycle+buf] * h[h.size()-k-1];
                buf += 1;
            }
            else{
                y[y_count] += x[n+k-h.size()+1] * h[h.size()-k-1];
            }
        }
        cycle += 5;
        y_count++;
    }

}

void effBlockConv(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
    int cycle = 0;
    int y_count = 0;
    int ini = 0;
    int jump_count = 0;
    //as n >= 2500, only for those k values, would do compuations
    int jump_ini_index[24] = {19, 14, 9, 4, 23, 18, 13, 8, 3, 22, 17, 12, 7, 2, 21, 16, 11, 6, 1, 20, 15, 10, 5, 0};

    for (auto n = 0; (n+125) < block_size; n += 125){
        int buf = 0;
        int jump = 24;

        if (n >= 2500){
            ini = jump_ini_index[jump_count];
            if (jump_count < 23){
                jump_count++;
            }
            else{
                jump_count = 0;
            }
        }

        //only for position 0 24 48 72... buf_wave has values
        //as n = 2500, buf_wave values would not be used at all
        for (auto k = ini; k < h.size(); k += jump){
            if ((cycle+buf) < (h.size()-1)){
                y[y_count] += buf_wave[cycle+buf] * h[h.size()-k-1];
                buf += 24;
            }
            else{
                if (n < 2500) jump = 1;
                else jump = 24;
                if ((n+k-h.size()+1)%24 == 0 && (n+k-h.size()+1) != 0){
                    y[y_count] += x[n+k-h.size()+1] * h[h.size()-k-1];
                }
            }
        }
        y[y_count] *= 24;
        cycle += 125;
        y_count++;
    }
}

void stereoBPConv(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
    int cycle = 0;

    for (auto n = 0; n < block_size; n += 1){
        int buf = 0;

        for (auto k = 0; k < h.size(); k += 1){
            if ((cycle+buf) < (h.size()-1)) {
                y[n] += buf_wave[cycle + buf] * h[h.size() - k - 1];
                buf += 1;
            }
            else{
                y[n] += x[n+k-h.size()+1] * h[h.size()-k-1];
            }
        }
        //y[n] *= 19;
        cycle += 1;
    }
}

void stereoBPConvEff(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
    int cycle = 0;
    int ini = 0;
    int jump_count = 0;
    //as n >= 1896, only for those k values, would do compuations
    int jump_ini_index[4] = {7, 31, 55, 79};

    for (auto n = 0; n < block_size; n += 24){
        int buf = 0;
        int jump = 24;

        if (n >= 1896){
            ini = jump_ini_index[jump_count];
            if (jump_count < 4){
                jump_count++;
            }
            else{
                jump_count = 1;
            }
        }

        for (auto k = ini; k < h.size(); k += jump){
            //would finish buffer part as n = 1872, as n = 1896, there would be no buf_wave values participate in computations
            if ((cycle+buf) < (h.size()-1)){
                y[n] += buf_wave[cycle+buf] * h[h.size()-k-1];
                buf += 24;
            }
            else{
                if (n < 2000) jump = 1;
                else jump = 24;
                // as buf_wave finished, x has non-zero values as k = 7 31 55 79
                if ((n+k-h.size()+1)%24 == 0 && (n+k-h.size()+1) != 0) {
                    y[n] += x[n + k - h.size() + 1] * h[h.size() - k - 1];
                }
            }
        }
        cycle += 1;
    }
}

void rdsRecBPFConv(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
    int cycle = 0;

    for (auto n = 0; n < block_size; n += 1){
        int buf = 0;

        for (auto k = 0; k < h.size(); k += 1){
            if ((cycle+buf) < (h.size()-1)) {
                //squaring
                y[n] += buf_wave[cycle + buf] * buf_wave[cycle + buf] * h[h.size() - k - 1];
                buf += 1;
            }
            else{
                y[n] += x[n+k-h.size()+1] * x[n+k-h.size()+1] * h[h.size()-k-1];
            }
        }
        cycle += 1;
    }
}

void rdsLPFConv_3k(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
    int cycle = 0;

    for (auto n = 0; n < block_size; n += 1){
        int buf = 0;

        for (auto k = 0; k < h.size(); k += 1){
            if ((cycle+buf) < (h.size()-1)) {
                y[n*19] += buf_wave[cycle + buf] * h[h.size() - k - 1];
                buf += 1;
            }
            else{
                y[n*19] += x[n+k-h.size()+1] * h[h.size()-k-1];
            }
        }
        cycle += 1;
    }
}

void rds_3k_conv(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size){

    int cycle = 0;
    int y_count = 0;
    int ini = 0;
    int buf_jump_count = 0;
    int jump_count = 0;
    //as n < 1920, used buf_wave values
    int jump_buf_index[19] = {0, 15, 11, 7, 3, 18, 14, 10, 6, 2, 17, 13, 9, 5, 1, 16, 12, 8, 4};
    //as n >= 1920
    int jump_ini_index[19] = {17, 13, 9, 5, 1, 16, 12, 8, 4, 0, 15, 11, 7, 3, 18, 14, 10, 6, 2};

    for (auto n = 0; n < block_size; n += 80){
        int buf = 0;
        int jump = 19;

        if (n < 1920){
            ini = jump_buf_index[buf_jump_count];
            if (buf_jump_count < 18){
                buf_jump_count++;
            }
            else {
                buf_jump_count = 0;
            }
        }
        else {
            ini = jump_ini_index[jump_count];
            if (jump_count < 18){
                jump_count++;
            }
            else {
                jump_count = 0;
            }
        }

        for (auto k = ini; k < h.size(); k += jump){
            if ((cycle+buf) < (h.size()-1)){
                y[y_count] += buf_wave[cycle+buf] * h[h.size()-k-1];
                buf += 19;
            }
            else{
                if (n < 1920) jump = 1;
                else jump = 19;
                if ((n+k-h.size()+1)%19 == 0 && (n+k-h.size()+1) != 0) {
                    y[y_count] += x[n + k - h.size() + 1] * h[h.size() - k - 1];
                }
            }
        }
        y[y_count] *= 19;
        cycle += 80;
        y_count++;
    }
}

void rdsRRCConv(std::vector<float> &y, std::vector<float> &buf_wave, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
    int cycle = 0;

    for (auto n = 0; n < block_size; n += 1){
        int buf = 0;

        for (auto k = 0; k < h.size(); k += 1){
            if ((cycle+buf) < (h.size()-1)){
                y[n] += buf_wave[cycle+buf] * h[h.size()-k-1];
                buf += 1;
            }
            else{
                y[n] += x[n + k - h.size() + 1] * h[h.size() - k - 1];
            }
        }
        cycle += 1;
        y[n] *= 19;
    }
}