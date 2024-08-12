#pragma once
#include <Cell.h>

extern void computeAllDens(std::vector<iisphparticle>& var_PartC);
extern void computeAllDenstwoD(std::vector<iisphparticle>& var_PartC);
extern void computeAllDensSSPH(std::vector<iisphparticle>& var_PartC);
extern void computeAllDensErr(std::vector<iisphparticle>& var_PartC);
extern float makeAlldenserrAvg(std::vector<iisphparticle>& var_PartC);