#pragma once
#include <Cell.h>
#include<omp.h>
//uniform grid
extern void clearuniformgrid();
extern void insertparticlesinuniformgrid(std::vector<iisphparticle>& var_PartC);
extern unsigned int uniformgridhash(float x, float y, float z);
extern void findAllNeighbourscompact2D(std::vector<iisphparticle>& var_PartC);
extern void findAllNeighbourscompact3D(std::vector<iisphparticle>& var_PartC);
//spatial hashing
extern unsigned int hashFunction3D(float x, float y, float z);
extern unsigned int hashFunction2D(float x, float y);
extern void insertAllParticlesIntoHashmap(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void insertAllParticlesIntoHashmap2D(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void findAllNeighbours(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void findAllNeighbours2D(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);



