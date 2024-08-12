#pragma once
#include <Cell.h>


extern void makecfltrue(std::vector<iisphparticle>& var_PartC);
extern bool checkcfl(float maxvel);
extern std::vector<glm::vec3> parseObjFile(const std::string& filePath);
