#pragma once
#include <Cell.h>


extern void makeboundmass(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void makepartmass(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void makepartmassTwoD(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void makeboundmassTwoD(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
