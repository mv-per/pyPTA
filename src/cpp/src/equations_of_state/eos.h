
#ifndef EOS_CPP
#define EOS_CPP

#include "./eos_helper.h"
#include "./pr77.h"
#include "./pr77_peneloux.h"
#include "./srk.h"
#include "./srk_peneloux.h"


typedef std::function<mono_eos(double, double)> call_mono_eos;

typedef std::function<mix_eos(std::vector<double>, double, double)> call_mix_eos;

#endif