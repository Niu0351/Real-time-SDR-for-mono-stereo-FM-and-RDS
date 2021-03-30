/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/
#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "stereo.h"
#include "rds.h"
#include <iostream>

#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>

#define PI 3.14159265358979323846

float rf_Fs_mode0 = 2.4e6;
float rf_Fs_mode1 = 2.5e6;
float rf_Fc = 100e3;
int rf_taps = 101;
int rf_decim = 10;

float audio_Fs = 48e3;
float audio_Fc = 16e3;
int audio_taps = 101;
int audio_decim = 5;
int audio_decim_1 = 125;
int mode = 0;

void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){

    std::vector<char> raw_data(num_samples, 0);
    std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));

    for (unsigned int k=0;k<num_samples;k++){
        block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
    }
}

void demodArctan(std::vector<float> &fm_demod, const std::vector<float> &I, const std::vector<float> &Q, int block_count, float &buf_I, float &buf_Q){

    for (int k=0;k<I.size();k++){
        if (k!=0){
            fm_demod[k] = (I[k] * (Q[k]-Q[k-1]) - Q[k] * (I[k]-I[k-1]))/(pow(I[k],2) + pow(Q[k],2));
            if (k == int(I.size()-1)){
                buf_I = I[k];
                buf_Q = Q[k];
            }
        }
        else if ((block_count == 0) && (k == 0)){
            fm_demod[k] = 0.0;
        }
        else if ((block_count != 0) && (k==0)){
            fm_demod[k] = (I[k] * (Q[k]-buf_Q) - Q[k] * (I[k]-buf_I))/(pow(I[k],2) + pow(Q[k],2));
        }

        //unwrap
        if (k>0){
            if ((fm_demod[k] - fm_demod[k-1]) > PI){
                fm_demod[k] = PI - (fm_demod[k] - fm_demod[k-1]);
            }
        }
    }
}


void demodArctan1(std::vector<float> &fm_demod, std::vector<float> &up_fm_demod, const std::vector<float> &I, const std::vector<float> &Q, int block_count, float &buf_I, float &buf_Q){
    int count = 0;
    for (int k=0;k<I.size();k++){
        if (k!=0){
            fm_demod[k] = (I[k] * (Q[k]-Q[k-1]) - Q[k] * (I[k]-I[k-1]))/(pow(I[k],2) + pow(Q[k],2));
            if (k == int(I.size()-1)){
                buf_I = I[k];
                buf_Q = Q[k];
            }
        }
        else if ((block_count == 0) && (k == 0)){
            fm_demod[k] = 0.0;
        }
        else if ((block_count != 0) && (k==0)){
            fm_demod[k] = (I[k] * (Q[k]-buf_Q) - Q[k] * (I[k]-buf_I))/(pow(I[k],2) + pow(Q[k],2));
        }

        //unwrap
        if (k>0){
            if ((fm_demod[k] - fm_demod[k-1]) > PI){
                fm_demod[k] = PI - (fm_demod[k] - fm_demod[k-1]);
            }
        }
        up_fm_demod[count] = fm_demod[k];
        count += 24;
    }
}

void demod(std::vector<float> &i_raw, std::vector<float> &q_raw, std::vector<float> &i_ds, std::vector<float> &q_ds, std::vector<float> &iq_data, std::vector<float> &fm_demod, std::vector<float> &fm_demod_1, std::vector<float> &rf_coeff, std::vector<float> &buf_wave_I, std::vector<float> &buf_wave_Q, float buf_I, float buf_Q, int mode, int block_count){

    blockConv_I(i_ds, buf_wave_I, iq_data, rf_coeff, iq_data.size());
    blockConv_Q(q_ds, buf_wave_Q, iq_data, rf_coeff, iq_data.size());
    int buf_count = 0;

    if(mode == 1){
        for (auto j = 49842; j < 50000; j += 2) {
            buf_wave_I[buf_count] = iq_data[j];
            buf_wave_Q[buf_count] = iq_data[j+1];
            buf_count += 1;
        }
        demodArctan1(fm_demod, fm_demod_1, i_ds, q_ds, block_count, buf_I, buf_Q);
    }
    else{
        for (auto j = 51042; j < 51200; j += 2) {
            buf_wave_I[buf_count] = iq_data[j];
            buf_wave_Q[buf_count] = iq_data[j+1];
            buf_count += 1;
        }
        demodArctan(fm_demod, i_ds, q_ds, block_count, buf_I, buf_Q);
    }

}

void audio_samples(std::vector<float> &fm_demod, std::vector<float> &audio_data, std::vector<float> &audio_coeff, std::vector<float> &buf_wave_audio, int mode){

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

struct ToAudio{
    std::vector<float> fm_demod;
    std::vector<float> fm_demod_1;
    std::vector<float> audio_coeff_mode;
};


void thread_demod(std::queue<ToAudio> &my_queue, std::queue<std::vector<float>> &rds_queue, std::mutex &my_mutex, std::mutex &rds_mutex, std::condition_variable &my_cvar, std::condition_variable &rds_cvar){
    int block_id = 0;
    //auto start_time = std::chrono::high_resolution_clock::now();

    //51200
    // 32k - 1M
    int BLOCK_SIZE = 0;
    int fm_demod_size = 0;
    int fm_demod_size_1 = 0;
    if (mode == 0){
        BLOCK_SIZE = 51200; //1024 * rf_decim * audio_decim;
        fm_demod_size = 2560; //(block_data.size()/2)/rf_decim
    }
    else if (mode ==1){
        BLOCK_SIZE = 50000; //1000 * rf_decim * audio_decim;
        fm_demod_size = 2500;
        fm_demod_size_1 = 24*2500;
    }

    //buf_wave
    std::vector<float> buf_wave_I(100, 0.0); //(rf_taps-1, 0.0)
    std::vector<float> buf_wave_Q(100, 0.0); //(rf_taps-1, 0.0)

    std::vector<float> block_data(BLOCK_SIZE, 0.0);
    std::vector<float> rf_coeff(101, 0.0); //(rf_taps, 0.0)
    std::vector<float> audio_coeff_mode0(101, 0.0); //(audio_taps, 0.0)
    std::vector<float> audio_coeff_mode1(24*audio_taps, 0.0); //(24*audio_taps, 0.0)

    if (mode ==0){
        // int(((block_data.size()/2)/rf_decim)/audio_decim);
        impulseResponseLPF(rf_Fs_mode0, rf_Fc, rf_taps, rf_coeff);
        impulseResponseLPF(24e4, audio_Fc, audio_taps, audio_coeff_mode0);
    }
    else if (mode ==1){
        // int((((block_data.size()/2)/rf_decim))*24/audio_decim_1);
        impulseResponseLPF(rf_Fs_mode1, rf_Fc, rf_taps, rf_coeff);
        impulseResponseLPF(6e6, audio_Fc, 24*audio_taps, audio_coeff_mode1);
    }

    std::vector<float> i_raw(int(BLOCK_SIZE/2), 0.0);
    std::vector<float> q_raw(int(BLOCK_SIZE/2), 0.0);

    // 25600/10 = 2560
    std::vector<float> i_ds(int(i_raw.size()/rf_decim), 0.0);
    std::vector<float> q_ds(int(q_raw.size()/rf_decim), 0.0);

    float buf_I = 0.0;
    float buf_Q = 0.0;

    ToAudio collection;
    while(true){
        readStdinBlockData(BLOCK_SIZE, block_id, block_data);
        if ((std::cin.rdstate()) != 0) {
            std::cerr << "End of input stream reached" << std::endl;
            //auto stop_time = std::chrono::high_resolution_clock::now();
            //std::chrono::duration<double, std::milli> run_time = stop_time - start_time;
            //std::cerr << "effBlockConv ran for " << run_time.count() << " milliseconds" << std::endl;
            exit(1);
        }

        //std::cerr << "Read block " << block_id << std::endl;

        // demod and filter for I Q data
        std::vector<float> fm_demod(fm_demod_size, 0.0);
        std::vector<float> fm_demod_1(fm_demod_size_1, 0.0);

        demod(i_raw, q_raw, i_ds, q_ds, block_data, fm_demod, fm_demod_1, rf_coeff, buf_wave_I, buf_wave_Q, buf_I, buf_Q, mode, block_id);

        std::unique_lock<std::mutex> my_lock(my_mutex);
        //std::unique_lock<std::mutex> rds_lock(rds_mutex);

        collection.fm_demod = fm_demod;
        collection.fm_demod_1 = fm_demod_1;
        if(mode == 0){
            collection.audio_coeff_mode = audio_coeff_mode0;
        }
        else{
            collection.audio_coeff_mode = audio_coeff_mode1;
        }

        // if(rds_queue.size() == QUEUE_BLOCKS){
        //   rds_cvar.wait(rds_lock);
        // }
        if(my_queue.size() == QUEUE_BLOCKS){
            my_cvar.wait(my_lock);
        }
        // if(rds_queue.size() == QUEUE_BLOCKS){
        //     rds_cvar.wait(rds_lock);
        // }

        my_queue.push(collection);
        if(mode == 0){
          rds_queue.push(fm_demod);
        }

        block_id++;

        my_lock.unlock();
        my_cvar.notify_one();
        //rds_lock.unlock();
        if(mode == 0){
            rds_cvar.notify_one();
        }
    }
}

void thread_audio(std::queue<ToAudio> &my_queue, std::mutex &my_mutex, std::condition_variable &my_cvar){
    std::vector<float> buf_wave_audio(2423, 0.0); //(audio_taps*24-1, 0.0)
    std::vector<float> bp_buf_wave_rec(100, 0.0); //(rf_taps-1, 0.0)
    std::vector<float> bp_buf_wave_extr(100, 0.0); //(rf_taps-1, 0.0)
    std::vector<float> bp_buf_wave_audio(2423, 0.0); //(audio_taps*24-1, 0.0)
    std::vector<short int> audio_out(0, 0);

    int count = 0;
    float feedbackI = 1.0;
    float feedbackQ = 0.0;
    float integrator = 0.0;
    float phaseEst = 0.0;
    float triArg = 0.0;
    float trigOffset = 0.0;
    float ncoOut_buf_I = 0.0;
    float ncoOut_buf_Q = 0.0;

    int audio_out_size = 0;
    int audio_size = 0;

    int block_id = 0;
    while(true){

        std::unique_lock<std::mutex>my_lock(my_mutex);
        if(my_queue.empty()){
            my_cvar.wait(my_lock);
        }

        ToAudio collection = my_queue.front();
        my_queue.pop();

        if ((std::cin.rdstate()) != 0) {    //not sure
            std::cerr << "End of input stream reached" << std::endl;
            exit(1);
        }

        if(mode == 0) audio_size = 512;
        else audio_size = 480;

        audio_out_size += audio_size*2; //audio_size
        //audio_out_size += audio_size; //audio_size
        audio_out.resize(audio_out_size);
        std::vector<float> audio_data(audio_size, 0.0);

        float feedbackI = 1.0;
        float feedbackQ = 0.0;

        std::vector<float> fm_demod = collection.fm_demod;
        std::vector<float> fm_demod_1 = collection.fm_demod_1;
        std::vector<float> audio_coeff_mode = collection.audio_coeff_mode;

        if(mode == 0){
            //audio_size = 512;
            audio_samples(fm_demod, audio_data, audio_coeff_mode, buf_wave_audio, mode);
        }
        else if (mode == 1){
            //audio_size = 480;
            audio_samples(fm_demod_1, audio_data, audio_coeff_mode, buf_wave_audio, mode);
        }

        ///////////////////////// STEREO PATH ///////////////////////
        std::vector<float> audio_data_stereo(audio_size, 0.0);
        std::vector<float> stereo_samples(audio_size, 0.0);
        std::vector<float> left_channel(audio_size, 0.0);
        std::vector<float> right_channel(audio_size, 0.0);

        stereoExtrRec(audio_data_stereo, audio_coeff_mode, stereo_samples, fm_demod, bp_buf_wave_rec, bp_buf_wave_extr, bp_buf_wave_audio, fm_demod.size(), 19e3, mode, block_id, feedbackI, feedbackQ, integrator, phaseEst, triArg, trigOffset, ncoOut_buf_I, ncoOut_buf_Q);

        combiner(audio_data_stereo, audio_data, left_channel, right_channel);

        for (unsigned int k = 0; k < left_channel.size(); k++) {
            if (std::isnan(left_channel[k])) audio_out[count] = 0;
            else audio_out[count] = static_cast<short int>(left_channel[k] * 16384);
            count++;

            if (std::isnan(right_channel[k])) audio_out[count] = 0;
            else audio_out[count] = static_cast<short int>(right_channel[k] * 16384);
            count++;
        }

        fwrite(&audio_out[block_id*left_channel.size()*2], sizeof(short int), left_channel.size()*2, stdout);

        block_id++;

        my_lock.unlock();
        my_cvar.notify_one();
    }
}

void thread_rds(std::queue<std::vector<float>> &rds_queue, std::mutex &my_mutex, std::condition_variable &rds_cvar){
    //3
    int max_index = 0;
    int buf_mach = 0;
    float feedbackI = 1.0;
    float feedbackQ = 0.0;
    int block_id = 0;

    std::vector<float> rds_buf(100, 0.0);
    std::vector<float> rds_carrier_buf(100, 0.0);
    std::vector<float> rds_up19_buf_I(1918, 0.0);//80*19-1
    std::vector<float> rds_up19_buf_Q(1918, 0.0);//80*19-1
    std::vector<float> rrc_buf_I(100, 0.0);
    std::vector<float> rrc_buf_Q(100, 0.0);
    std::vector<float> rds_carrier_buf_I(100, 0.0);
    std::vector<float> rds_carrier_buf_Q(100, 0.0);
    float integrator = 0.0;
    float phaseEst = 0.0;
    float triArg = 0.0;
    float trigOffset = 0.0;
    float ncoOut_buf_I = 0.0;
    float ncoOut_buf_Q = 0.0;

    std::vector<float> extr_coeff(101, 0.0);
    std::vector<float> fm_demod_rdsRec_coeff(101, 0.0);
    std::vector<float> rds_3k(101, 0.0);
    std::vector<float> rds_rationalLPF_coeff(19*101, 0.0);
    std::vector<float> rrc_coeff(101, 0.0);
    bandpassRDS(extr_coeff, 0);
    bpRDSRec(fm_demod_rdsRec_coeff, 0);
    impulseResponseLPF(2.4e5, 3e3, 101, rds_3k);
    impulseResponseLPF(2.4e5*19, 16e3, 19*101, rds_rationalLPF_coeff);
    impulseResponseRootRaiseCosine(57e3, 101, rrc_coeff);


    std::vector<int> sync_I(38, 0.0);
    std::vector<int> sync_Q(38, 0.0);

    std::vector<float> rrc_outputs_I_3blocks(1824, 0.0);//608*3
    std::vector<float> rrc_outputs_Q_3blocks(1824, 0.0);

    while(true){

      std::unique_lock<std::mutex>my_lock(my_mutex);
      if(rds_queue.empty()){
          rds_cvar.wait(my_lock);
      }

      std::vector<float> fm_demod = rds_queue.front();
      rds_queue.pop();

      if ((std::cin.rdstate()) != 0) {    //not sure
          std::cerr << "End of input stream reached" << std::endl;
          exit(1);
      }

      std::vector<float> rds_samples_I(fm_demod.size(), 0.0);
      std::vector<float> rds_samples_Q(fm_demod.size(), 0.0);

      if (block_id % 3 == 0 && block_id != 0){
        rrc_outputs_I_3blocks.resize(1824, 0.0);
        rrc_outputs_Q_3blocks.resize(1824, 0.0);
      }

      rds(sync_I, sync_Q, rrc_outputs_I_3blocks, rrc_outputs_Q_3blocks, rrc_coeff, extr_coeff, fm_demod_rdsRec_coeff,
            rds_3k, rds_rationalLPF_coeff, block_id, rds_samples_I, rds_samples_Q,
            fm_demod, rds_buf, rds_carrier_buf_I, rds_carrier_buf_I, rds_carrier_buf,
            rds_up19_buf_I, rds_up19_buf_Q, rrc_buf_I, rrc_buf_Q, fm_demod.size(),
            114e3, mode, feedbackI, feedbackQ, integrator, phaseEst, triArg,
            trigOffset, ncoOut_buf_I, ncoOut_buf_Q, max_index, buf_mach);

      block_id++;

      //rds_lock.unlock();
      rds_cvar.notify_one();
    }

}

int main(int argc, char* argv[])
{
    //auto start_time = std::chrono::high_resolution_clock::now();
    if (argc < 2){
        std::cerr << "Operating in default mode 0" << std::endl;
    }
    else if (argc == 2){
        mode = atoi(argv[1]);

        if (mode != 1){
            std::cerr << "Wrong mode " << mode << std::endl;
            exit(1);
        }
        else {
            std::cerr << "Usage: " << argv[0] << std::endl;
            std::cerr << "or " << std::endl;
            std::cerr << "Usge: " << argv[0] << " 1" << std::endl;
            //exit(1);
        }
    }

    std::queue<ToAudio> my_queue;
    std::queue<std::vector<float>> rds_queue;
    std::mutex my_mutex;
    std::mutex rds_mutex;
    std::condition_variable my_cvar;
    std::condition_variable rds_cvar;

    std::thread front_thread = std::thread(thread_demod, std::ref(my_queue), std::ref(rds_queue), std::ref(my_mutex), std::ref(rds_mutex), std::ref(my_cvar), std::ref(rds_cvar));
    std::thread audio_thread = std::thread(thread_audio, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar));
    std::thread rds_thread = std::thread(thread_rds, std::ref(rds_queue), std::ref(rds_mutex), std::ref(rds_cvar));

    front_thread.join();
    audio_thread.join();
    if(mode == 0){
      rds_thread.join();
    }
    return 0;
}

//std::vector<float> rds_samples(19*fm_demod.size(), 0.0);
//std::vector<float> rds_buf(79, 0.0);
//std::vector<float> rds_carrier_buf(79, 0.0);
//std::vector<float> rds_up19_buf(1519, 0.0);//80*19-1
//std::vector<float> rrc_buf(79, 0.0);
//float feedbackI = 1.0;
//float feedbackQ = 0.0;
//float integrator = 0.0;
//float phaseEst = 0.0;
//float triArg = 0.0;
//float trigOffset = 0.0;
//float ncoOut_buf = 0.0;
//rds(rds_samples, fm_demod, rds_buf, rds_carrier_buf, rds_up19_buf, rrc_buf, fm_demod.size(), 57e3, mode, feedbackI, feedbackQ, integrator, phaseEst, triArg, trigOffset, ncoOut_buf);
