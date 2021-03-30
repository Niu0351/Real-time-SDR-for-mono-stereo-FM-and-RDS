//
// Created by kaka on 3/19/21.
//

#include "dy4.h"
#include "rds.h"
#include "filter.h"
#include "logfunc.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

float rds_fb_extr = 54e3;
float rds_fe_extr = 60e3;
float rds_fb_rec = 113.5e3;
float rds_fe_rec = 114.5e3;

float rds_Fs_mode0 = 2.4e5;
float rds_Fs_mode1 = 2.5e5;
int rds_bp_taps = 101;

//in order to get 57 kHz as ncoOutput
float rds_ncoScale = 0.5;
float rds_phaseAdjust = 0;
float rds_normBandwidth = 0.005;

void rds_fmPLL(std::vector<float> &ncoOut_I, std::vector<float> &ncoOut_Q, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, float &feedbackI, float &feedbackQ, float &integrator, float &phaseEst, float &triArg, float &trigOffset, float &ncoOut_buf_I, float &ncoOut_buf_Q){

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

void bandpassRDS(std::vector<float> &h_extr, int mode){

    float Fs = 2.4e5;

    float normCenter = (rds_fe_extr+rds_fb_extr)/Fs;
    float normPass = (rds_fe_extr-rds_fb_extr)*2/Fs;

    for (auto i=0;i<rds_bp_taps;i++){
        if (i == (rds_bp_taps-1)/2){
            h_extr[i] = normPass;
        }
        else{
            h_extr[i] = normPass * (sin(PI*(normPass/2)*(i-(rds_bp_taps-1)/2)))/(PI*(normPass/2)*(i-(rds_bp_taps-1)/2));
        }

        h_extr[i] *= cos(PI*i*normCenter);
        h_extr[i] *= sin(i*PI/rds_bp_taps)*sin(i*PI/rds_bp_taps);
    }
}

void bpRDSRec(std::vector<float> &h_rec, int mode){

    float Fs = 2.4e5;

    float normCenter = (rds_fe_rec+rds_fb_rec)/Fs;
    float normPass = (rds_fe_rec-rds_fb_rec)*2/Fs;

    for (auto i = 0;i < rds_bp_taps;i += 1){
        if (i == (rds_bp_taps-1)/2){
            h_rec[i] = normPass;
        }
        else{
            h_rec[i] = normPass * (sin(PI*(normPass/2)*(i-(rds_bp_taps-1)/2)))/(PI*(normPass/2)*(i-(rds_bp_taps-1)/2));
        }

        h_rec[i] *= cos(PI*i*normCenter);
        h_rec[i] *= sin(i*PI/rds_bp_taps)*sin(i*PI/rds_bp_taps);
    }
}

void rds_mixer(std::vector<float> &rds_samples_I, std::vector<float> &rds_samples_Q, std::vector<float> &filter_demod_extr, std::vector<float> &ncoOut_I, std::vector<float> &ncoOut_Q){
    for (auto i=0; i<filter_demod_extr.size();i+=1){
        rds_samples_I[i] = filter_demod_extr[i] * ncoOut_I[i];
        rds_samples_Q[i] = filter_demod_extr[i] * ncoOut_Q[i];
    }
}

template<typename T>
std::vector<T> arange(T start, T stop, T step = 1){
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

void rdsDataRecovery(std::vector<float> &rrc_outputs_I_3blocks, std::vector<float> &rrc_outputs_Q_3blocks, std::vector<float> &rrc_rec_I, std::vector<float> &rrc_rec_Q, int block_id, int &max_index){

    if (block_id == 1){

        int max_index_first = 0;
        int find_index = 0;
        for (auto i=0;i<24;i++){
            int flag = 0;

            for (auto j=0;j<rrc_rec_I.size()/3;j+=1){
                if (j%2 == 1 && ((rrc_outputs_I_3blocks[24*j+i] > 0 && rrc_outputs_I_3blocks[24*(j-1)+i] < 0) || (rrc_outputs_I_3blocks[24*j+i] < 0 && rrc_outputs_I_3blocks[24*(j-1)+i] > 0))) {
                    flag = 1;
                }
                else if (j%2 == 1 && ((rrc_outputs_I_3blocks[24*j+i] < 0 && rrc_outputs_I_3blocks[24*(j-1)+i] < 0) || (rrc_outputs_I_3blocks[24*j+i] > 0 && rrc_outputs_I_3blocks[24*(j-1)+i] > 0))){
                    flag = 0;
                    break;
                }

                if ((j==rrc_rec_I.size()-1) && flag == 1) {
                    find_index = i;
                }
            }

            if (abs(rrc_outputs_I_3blocks[max_index_first]) < abs(rrc_outputs_I_3blocks[find_index])) max_index_first = find_index;
        }
        max_index = max_index_first;
    }

    for (auto j=0;j<rrc_rec_I.size();j+=1){
        rrc_rec_I[j] = rrc_outputs_I_3blocks[24*j+max_index];
        rrc_rec_Q[j] = rrc_outputs_Q_3blocks[24*j+max_index];
    }

}

//put in IF data as fm_demod
void rds(std::vector<float> &rrc_outputs_I_3blocks, std::vector<float> &rrc_outputs_Q_3blocks, std::vector<float> &extr_coeff, std::vector<float> &fm_demod_rdsRec_coeff, std::vector<float> &rds_3k, std::vector<float> &rds_rationalLPF_coeff, int block_id, std::vector<float> &rds_samples_I, std::vector<float> &rds_samples_Q, std::vector<float> &fm_demod, std::vector<float> &rds_buf, std::vector<float> &rds_carrier_buf_I, std::vector<float> &rds_carrier_buf_Q, std::vector<float> &rds_carrier_buf, std::vector<float> &rds_up19_buf_I, std::vector<float> &rds_up19_buf_Q, std::vector<float> &rrc_buf_I, std::vector<float> &rrc_buf_Q, int block_size, float freq, int mode, float &feedbackI, float &feedbackQ, float &integrator, float &phaseEst, float &triArg, float &trigOffset, float &ncoOut_buf_I, float &ncoOut_buf_Q, int &max_index){

    std::vector<float> filter_demod_extr(fm_demod.size(), 0.0);
    std::vector<float> filter_demod_rec(fm_demod.size(), 0.0);
    std::vector<float> squaring_results(fm_demod.size(), 0.0);
    std::vector<float> ncoOut_I(fm_demod.size()+1, 0.0);
    std::vector<float> ncoOut_Q(fm_demod.size()+1, 0.0);
    float Fs = 2.4e5;

    //RDS Extraction Channel
    stereoBPConv(filter_demod_extr, rds_buf, fm_demod, extr_coeff, fm_demod.size());//same as stereo
    int buf_count = 0;
    for (auto i = (fm_demod.size() - extr_coeff.size() + 1); i < fm_demod.size(); i+=1) {
        rds_buf[buf_count] = fm_demod[i];
        buf_count += 1;
    }

    //RDS Carrier Recovery
    //squaring non-linear is combined into conv process
    rdsRecBPFConv(filter_demod_rec, rds_carrier_buf, filter_demod_extr, fm_demod_rdsRec_coeff, squaring_results.size());//same as stereo
    buf_count = 0;
    for (auto i = (filter_demod_extr.size() - fm_demod_rdsRec_coeff.size() + 1); i < filter_demod_extr.size(); i+=1) {
        rds_carrier_buf[buf_count] = filter_demod_extr[i];
        buf_count += 1;
    }

    //pll
    rds_fmPLL(ncoOut_I, ncoOut_Q, filter_demod_rec, freq, Fs, rds_ncoScale, rds_phaseAdjust, rds_normBandwidth, feedbackI, feedbackQ, integrator, phaseEst, triArg, trigOffset, ncoOut_buf_I, ncoOut_buf_Q);

    //mixer
    rds_mixer(rds_samples_I, rds_samples_Q, filter_demod_extr, ncoOut_I, ncoOut_Q);

    //3k LPF combined with upsampled by 19
    std::vector<float> rds_3k_results_I_up19(rds_samples_I.size()*19, 0.0);
    std::vector<float> rds_3k_results_Q_up19(rds_samples_Q.size()*19, 0.0);
    rdsLPFConv_3k(rds_3k_results_I_up19, rds_carrier_buf_I, rds_samples_I, rds_3k, rds_samples_I.size());
    rdsLPFConv_3k(rds_3k_results_Q_up19, rds_carrier_buf_Q, rds_samples_Q, rds_3k, rds_samples_I.size());

    buf_count = 0;
    for (auto i=2470;i<rds_samples_I.size();i+=19){
        rds_carrier_buf_I[buf_count*19] = rds_samples_I[i];
        rds_carrier_buf_Q[buf_count*19] = rds_samples_Q[i];
        buf_count++;
    }

    //16k LPF combined with downsampled by 80
    //non-zero values in filter_rds_up19 = 608
    std::vector<float> filter_rds_down80_I(rds_3k_results_I_up19.size()/80,0.0);
    std::vector<float> filter_rds_down80_Q(rds_3k_results_Q_up19.size()/80,0.0);

    rds_3k_conv(filter_rds_down80_I, rds_up19_buf_I, rds_3k_results_I_up19, rds_rationalLPF_coeff, rds_3k_results_I_up19.size());
    rds_3k_conv(filter_rds_down80_Q, rds_up19_buf_Q, rds_3k_results_Q_up19, rds_rationalLPF_coeff, rds_3k_results_Q_up19.size());
    buf_count = 0;
    for (auto i=(rds_3k_results_I_up19.size()-rds_rationalLPF_coeff.size());i<rds_3k_results_I_up19.size();i++){
        rds_up19_buf_I[buf_count] = rds_3k_results_I_up19[i];
        rds_up19_buf_Q[buf_count] = rds_3k_results_Q_up19[i];
        buf_count++;
    }

    //RRC
    std::vector<float> rrc_coeff(rds_bp_taps, 0.0);
    impulseResponseRootRaiseCosine(57e3, rds_bp_taps, rrc_coeff);

    std::vector<float> rrc_outputs_I(filter_rds_down80_I.size(), 0.0);
    std::vector<float> rrc_outputs_Q(filter_rds_down80_Q.size(), 0.0);
    stereoBPConv(rrc_outputs_I, rrc_buf_I, filter_rds_down80_I, rrc_coeff, filter_rds_down80_I.size());
    stereoBPConv(rrc_outputs_Q, rrc_buf_Q, filter_rds_down80_Q, rrc_coeff, filter_rds_down80_I.size());
    buf_count = 0;
    for (auto i=(filter_rds_down80_I.size()-rrc_coeff.size());i<filter_rds_down80_I.size();i++){
        rrc_buf_I[buf_count] = filter_rds_down80_I[i];
        rrc_buf_Q[buf_count] = filter_rds_down80_Q[i];
        buf_count++;
    }

    if (block_id != 0){
        for (auto i=0;i<rrc_outputs_I.size();i++){
            rrc_outputs_I_3blocks[(block_id-1)%3*rrc_outputs_I.size()+i] = rrc_outputs_I[i];
            rrc_outputs_Q_3blocks[(block_id-1)%3*rrc_outputs_Q.size()+i] = rrc_outputs_Q[i];
        }
    }

    //recovery
    //store 3 blocks of data, and starts from block 0
    std::vector<float> rrc_rec_I(1824/24, 0.0);//608*3/24
    std::vector<float> rrc_rec_Q(1824/24, 0.0);

    if (block_id%3 == 0 && block_id != 0){
        rdsDataRecovery(rrc_outputs_I_3blocks, rrc_outputs_Q_3blocks, rrc_rec_I, rrc_rec_Q, block_id, max_index);
    }

    float df = 1;
    std::vector<float> fre = arange<float> (0, rrc_outputs_I.size(), df);

    if (block_id == 24){
       logVector("rrc_output", rrc_rec_I, rrc_rec_Q); // log only positive freq
//        logVector("rrc_output1", fre, rrc_outputs_I);
        std::cerr << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png\n";
    }

}