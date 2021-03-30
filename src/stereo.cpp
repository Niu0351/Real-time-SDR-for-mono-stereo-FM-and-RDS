//
// Created by kaka on 3/11/21.
//

#include "dy4.h"
#include "stereo.h"
#include "filter.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846

float fb_rec = 18.5e3;
float fe_rec = 19.5e3;
float fb_extr = 22e3;
float fe_extr = 54e3;

float stereo_Fs_mode0 = 2.4e5;
float stereo_Fs_mode1 = 2.5e5;
int stereo_bp_taps = 101;

float ncoScale = 2.0;
float phaseAdjust = 0.0;
float normBandwidth = 0.01;

void fmPLL(std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, float &feedbackI, float &feedbackQ, float &integrator, float &phaseEst, float &triArg, float &trigOffset, float &ncoOut_buf_I, float &ncoOut_buf_Q){

    std::vector<float> ncoOut_I(pllIn.size()+1, 0.0);
    std::vector<float> ncoOut_Q(pllIn.size()+1, 0.0);

    // scale factors for proportional/integrator terms
    float Cp = 2.666;
    float Ci = 3.555;

    //gain for the proportional term
    float Kp = normBandwidth * Cp;
    //gain for the integrator term
    float Ki = normBandwidth * normBandwidth * Ci;

    //initialize internal state
    ncoOut_I[0] = ncoOut_buf_I;
    ncoOut_Q[0] = ncoOut_buf_Q;

    for (auto i=0;i<pllIn.size();i+=1) {
        //phase detector
        float errorI = pllIn[i] * feedbackI; // complex conjugate of the
        float errorQ = pllIn[i] * (-1) * feedbackQ; //feedback complex exponential

        //four-quadrant arctangent discriminator for phase error detection
        float errorD = atan2(errorQ,errorI);

        //loop filter
        integrator += Ki * errorD;
        //update phase estimate
        phaseEst += (Kp * errorD) + integrator;
        //internal oscillator
        triArg = 2 * PI * (freq/Fs) * (trigOffset+i+1) + phaseEst;
        feedbackI = cos(triArg);
        feedbackQ = sin(triArg);
        ncoOut_I[i+1] = cos(triArg * ncoScale + phaseAdjust);
        ncoOut_Q[i+1] = sin(triArg * ncoScale + phaseAdjust);
    }
    trigOffset += pllIn.size();
    ncoOut_buf_I = ncoOut_I[ncoOut_I.size()-1];
    ncoOut_buf_Q = ncoOut_Q[ncoOut_Q.size()-1];

}

void bandRec(std::vector<float> &h_rec, int mode){

    float Fs = 0.0;
    if (mode == 0){
        Fs = stereo_Fs_mode0;
    }
    else{
        Fs = stereo_Fs_mode1;
    }

    float normCenter = (fe_rec+fb_rec)/Fs;
    float normPass = (fe_rec-fb_rec)*2/Fs;

    for (auto i = 0;i < stereo_bp_taps;i += 1){
        if (i == (stereo_bp_taps-1)/2){
            h_rec[i] = normPass;
        }
        else{
            h_rec[i] = normPass * ((sin(PI*(normPass/2)*(i-(stereo_bp_taps-1)/2)))/(PI*(normPass/2)*(i-(stereo_bp_taps-1)/2)));
        }

        h_rec[i] *= cos(PI*i*normCenter);
        h_rec[i] *= (sin(i*PI/stereo_bp_taps))*(sin(i*PI/stereo_bp_taps));
    }
}

void bandExtr(std::vector<float> &h_extr, int mode){

    float Fs = 0.0;
    if (mode == 0){
        Fs = stereo_Fs_mode0;
    }
    else{
        Fs = stereo_Fs_mode1;
    }

    float normCenter = (fe_extr+fb_extr)/Fs;
    float normPass = (fe_extr-fb_extr)*2/Fs;

    for (auto i=0;i<stereo_bp_taps;i++){
        if (i == (stereo_bp_taps-1)/2){
            h_extr[i] = normPass;
        }
        else{
            h_extr[i] = normPass * ((sin(PI*(normPass/2)*(i-(stereo_bp_taps-1)/2)))/(PI*(normPass/2)*(i-(stereo_bp_taps-1)/2)));
        }

        h_extr[i] *= cos(PI*i*normCenter);
        h_extr[i] *= (sin(i*PI/stereo_bp_taps))*(sin(i*PI/stereo_bp_taps));
    }
}

void mixer(std::vector<float> &stereo_samples, std::vector<float> &filter_demod_extr, std::vector<float> &ncoOut){
    for (auto i=0; i<stereo_samples.size();i+=24){
        stereo_samples[i] = filter_demod_extr[i] * ncoOut[i];
    }
}

//same digital filtering as mono path
void audio_samples_stereo(std::vector<float> &fm_demod, std::vector<float> &audio_data, std::vector<float> &audio_coeff, std::vector<float> &buf_wave_audio, int mode){

    //2560
    if(mode == 1){
        // 2560 * 24 = 61440
        effBlockConv(audio_data, buf_wave_audio, fm_demod, audio_coeff, fm_demod.size());
    }
    else if (mode == 0){
        // 2560
        blockConv1(audio_data, buf_wave_audio, fm_demod, audio_coeff, fm_demod.size());
    }

    int buf_count = 0;
    for (auto j = (fm_demod.size() - audio_coeff.size() + 1); j < fm_demod.size(); j++) {
        buf_wave_audio[buf_count] = fm_demod[j];
        buf_count += 1;
    }
}

//input as fm_demod, buf_wave should be defined as it is called
void stereoExtrRec(std::vector<float> &audio_data, std::vector<float> &audio_coeff, std::vector<float> &stereo_samples, std::vector<float> &fm_demod, std::vector<float> &bp_buf_wave_rec, std::vector<float> &bp_buf_wave_extr, std::vector<float> &bp_buf_wave_audio, int block_size, float freq, int mode, int block_id, float &feedbackI, float &feedbackQ, float &integrator, float &phaseEst, float &triArg, float &trigOffset, float &ncoOut_buf_I, float &ncoOut_buf_Q) {

    std::vector<float> filter_demod_rec(fm_demod.size(), 0.0);
    std::vector<float> filter_demod_extr(fm_demod.size(), 0.0);
    std::vector<float> rec_coeff(stereo_bp_taps, 0.0);
    std::vector<float> extr_coeff(stereo_bp_taps, 0.0);
    std::vector<float> ncoOut(fm_demod.size(), 0.0);

    bandRec(rec_coeff, mode);
    bandExtr(extr_coeff, mode);

    stereoBPConv(filter_demod_extr, bp_buf_wave_rec, fm_demod, rec_coeff, block_size);
    stereoBPConv(filter_demod_rec, bp_buf_wave_extr, fm_demod, extr_coeff, block_size);

    int buf_count = 0;
    //int jump = 24;
    //pos 8 32...
    for (auto i = (fm_demod.size() - rec_coeff.size() + 1); i < fm_demod.size(); i+=1) {
        bp_buf_wave_rec[buf_count] = fm_demod[i];
        bp_buf_wave_extr[buf_count] = fm_demod[i];
        buf_count += 1;
    }

    if (mode == 0) {
        fmPLL(filter_demod_rec, freq, stereo_Fs_mode0, ncoScale, phaseAdjust, normBandwidth, feedbackI, feedbackQ, integrator, phaseEst, triArg, trigOffset,ncoOut_buf_I, ncoOut_buf_Q);
    }
    else {
        fmPLL(filter_demod_rec, freq, stereo_Fs_mode1, ncoScale, phaseAdjust, normBandwidth, feedbackI, feedbackQ, integrator, phaseEst, triArg, trigOffset,ncoOut_buf_I, ncoOut_buf_Q);
    }

    mixer(stereo_samples, filter_demod_extr, ncoOut);

    audio_samples_stereo(stereo_samples, audio_data, audio_coeff, bp_buf_wave_audio, mode);

}

//to ensure that mono would finish first, need a flag in main and then call combiner in main
void combiner(std::vector<float> &stereo_data, std::vector<float> &mono_data, std::vector<float> &left_channel, std::vector<float> &right_channel){
    for (int i=0;i<stereo_data.size();i++){
        left_channel[i] = (stereo_data[i]+mono_data[i])/2;
        right_channel[i] = (stereo_data[i]-mono_data[i])/2;
    }
}