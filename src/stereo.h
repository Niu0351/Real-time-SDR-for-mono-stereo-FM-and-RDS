//
// Created by kaka on 3/11/21.
//

#ifndef PROJECT_GROUP46_THURSDAY_STEREO_H
#define PROJECT_GROUP46_THURSDAY_STEREO_H
//#ifndef DY4_FILTER_H
//#define DY4_FILTER_H

#include <iostream>
#include <vector>

void fmPLL(std::vector<float> &, std::vector<float> &, float, float, float, float, float, float &, float &, float &, float &, float &, float &, float &, float &);
void bandRec(std::vector<float> &, int);
void bandExtr(std::vector<float> &, int);
void mixer(std::vector<float> &, std::vector<float> &, std::vector<float> &);
void audio_samples_stereo(std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, int);
void stereoExtrRec(std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, int, float, int, int, float &, float &, float &, float &, float &, float &, float &, float &);
void combiner(std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &);

#endif //PROJECT_GROUP46_THURSDAY_STEREO_H
