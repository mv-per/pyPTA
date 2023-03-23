
#ifndef EOS_CPP
#define EOS_CPP

#include "./pr77.h"
#include "./eos_helper.h"

typedef std::function<mono_eos(double, double)> call_mono_eos;

typedef std::function<mix_eos(std::vector<double>, double, double)> call_mix_eos;

#endif