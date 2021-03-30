/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void impulseResponseRootRaiseCosine(float, int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void blockConv(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void blockConv_I(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void blockConv_Q(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void effBlockConv(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void blockConv1(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void stereoBPConv(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void stereoBPConvEff(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void rdsRecBPFConv(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void rdsLPFConv_3k(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void rds_3k_conv(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void rdsRRCConv(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);

#endif // DY4_FILTER_H
