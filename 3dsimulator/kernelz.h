#pragma once
#include <Cell.h>
extern void makeAllKernelAndKernelDer(std::vector<iisphparticle>& var_PartC);
extern void makeAllKernelAndKernelDerTwoD(std::vector<iisphparticle>& var_PartC);
extern float makesinglekernel(glm::vec3& posi, glm::vec3& posj);
extern float makesinglekernel2D(glm::vec3& posi, glm::vec3& posj);
extern glm::vec3 makesinglekernelder(glm::vec3& posi, glm::vec3& posj);