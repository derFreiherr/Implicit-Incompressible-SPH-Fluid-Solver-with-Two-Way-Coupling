#pragma once
#include <Cell.h>


extern void ssphAlgo(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap);
extern void ssphAlgotwoD(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap);

extern void var_iisph(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap);
extern void twoDiisph(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap);