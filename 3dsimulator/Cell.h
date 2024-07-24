#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <chrono>
#include<omp.h>
#include <random>
#include "Eigen/Dense"
#include <iostream>
#include <fstream>
using namespace Eigen;
extern int var_spezialboundpart;
extern float cellsize;
extern int hashsize;
extern const int p1;
extern const int p2;
extern const int p3;
extern float maxvel;
extern float cfl;
extern bool cfl_cond;
extern float cfl_max;
extern float cfl_min;
extern float deltaTmax;
extern float deltaTmin;
extern float deltaT;
extern const float h;
extern const float p0p0;
extern const float p0;
extern float currentdens;
extern float neighSearchTime;
extern float  drawingTime; 
extern float othercomputationTime;
extern float kerneltime ;
extern float densitytime ;
extern float nonpresatime ;
extern float predveltime ;
extern float computeallsftime ;
extern float computeallafftime ;
extern float looptime ;
extern  int var_MaxParticles;
extern  int var_fluidpart;
extern  int var_nx;
extern  int var_ny;
extern  int var_nz;
extern  int var_oneway;
extern float gravity;
extern float denserrold;
extern float visc_fluid;
extern float visc_boundary ;
extern bool clampp0;
extern float max_singledens;
extern int num_rot_part;
extern float pi_F;
extern float baseRotationSpeed;
extern float angle;
extern glm::vec3 oldpos;
extern int animationneighbourscount;
extern float clampfac;
extern bool absinterrupt;
extern int compbord;
extern int totalcomp;
extern bool fixeddt; 
extern bool ignoreincomplete;
extern float surfacetension;
extern bool highlightunderpop ;
extern bool highlightextremefast;
extern bool paussimul;
extern int numofp0high;
extern bool ignorep0tolow;
extern float animationtime; 
extern bool exportanimation;
extern float upperviualbord;
extern float lowervisualbord ;
extern float presup;
extern bool clamptolow;
extern bool colorvel;
extern bool colordens;
extern bool colorpres;
extern float gammadens;
extern float gammapres;
extern float gammabound;
extern float gammapart;
extern float jitterfac;
extern const float alphaTwoD;
extern std::vector<int> usemefordens;
extern std::vector<float> maxvels;
extern bool highlightusedforcomp;
extern bool highlightdenserr;
extern float gridbreite;
extern float gridhöhe;
extern float allrigidmass;
extern glm::vec3 posofcenterofmass;
extern glm::vec3 velofcenterofmass ;
extern glm::mat3x3 rotMat;
extern glm::vec3 xCM;
extern glm::vec3 vCM;
extern glm::mat3 A;
extern glm::vec3 L;
extern glm::mat3 I_inv;
extern glm::mat3 inertiaTensor;
extern glm::vec3 omegarigidbody;
extern glm::mat3 inertiaTensorInverse;


struct Cell {
	std::vector<int> particles;
};
class iisphparticle;
struct Particle;
extern unsigned int hashFunction3D(float x, float y, float z);
extern void makecfltrue(std::vector<iisphparticle>& var_PartC);
extern bool checkcfl(float maxvel);
extern void var_init(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);

extern void insertAllParticlesIntoHashmap(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void findAllNeighbours(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void makeAllKernelAndKernelDer(std::vector<iisphparticle>& var_PartC);
extern void computeAllDens(std::vector<iisphparticle>& var_PartC);
extern void computeAllDensSSPH(std::vector<iisphparticle>& var_PartC);
extern void MakeAllNonpresA(std::vector<iisphparticle>& var_PartC); 
extern void MakeAllNonpresAwithdens(std::vector<iisphparticle>& var_PartC);
extern void PredictAllVel(std::vector<iisphparticle>& var_PartC);
extern void computeAllSF(std::vector<iisphparticle>& var_PartC);
extern void makeAllAff1(std::vector<iisphparticle>& PartC);
extern void makeAllAfffast(std::vector<iisphparticle>& PartC);
extern void makeAllAfffastwithdens(std::vector<iisphparticle>& PartC);
extern void computeAllDensErr(std::vector<iisphparticle>& var_PartC);
extern float makeAlldenserrAvg(std::vector<iisphparticle>& var_PartC);
extern void makeAllPresA(std::vector<iisphparticle>& var_PartC);
extern void makeAllPresAwithdens(std::vector<iisphparticle>& var_PartC);
extern void makeAllAPupdatePres(std::vector<iisphparticle>& var_PartC);
extern void secondloop(std::vector<iisphparticle>& var_PartC);
extern void makeAllVandP(std::vector<iisphparticle>& var_PartC);
extern void var_teslavalve(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void var_teslavalveclosed(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void var_real_teslavalve(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void ssphAlgo(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap);
extern void computeAllPres(std::vector<iisphparticle>& var_PartC);
extern void makeAllA(std::vector<iisphparticle>& var_PartC);
extern void watercolumn(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void watercolumnsmall(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void moving_boundary(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void dambreaktest(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void smalldambreaktest(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void makeboundmass(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void makepartmass(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void MakeAllNonpresAtwoD(std::vector<iisphparticle>& var_PartC);
extern void makeAllKernelAndKernelDerTwoD(std::vector<iisphparticle>& var_PartC);
extern void watercolumnTwoD(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void makeAllPresAtwoD(std::vector<iisphparticle>& var_PartC);
extern void makepartmassTwoD(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void makeboundmassTwoD(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void DambreaktestTwoD(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void makeAllPresAwithbound(std::vector<iisphparticle>& var_PartC);
extern void makeBoundPres(std::vector<iisphparticle>& var_PartC);
extern float makesinglekernel(glm::vec3& posi, glm::vec3& posj);
extern glm::vec3 makesinglekernelder(glm::vec3& posi, glm::vec3& posj);
extern void makeAllPresAwithboundtwoD(std::vector<iisphparticle>& var_PartC);
extern void makeAllAtwoD(std::vector<iisphparticle>& var_PartC);
void ssphAlgotwoD(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap);
extern unsigned int hashFunction2D(float x, float y);
extern void findAllNeighbours2D(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void insertAllParticlesIntoHashmap2D(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void makeAllAfffparallel(std::vector<iisphparticle>& PartC);
extern void findAllNeighbourscompact2D(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void makeAllAfffparallelfast(std::vector<iisphparticle>& PartC);
extern void calculatecenterofmass(std::vector<iisphparticle>& PartC);
extern void initrigidbodies(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap);
extern void updaterigidbody(std::vector<iisphparticle>& PartC);
extern glm::mat3 skewSymmetricMatrix(const glm::vec3& v);
extern void makefloatingpartmass2d(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void makefloatingpartmass(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern std::vector<glm::vec3> parseObjFile(const std::string& filePath);
extern float makesinglekernel2D(glm::vec3& posi, glm::vec3& posj);
extern void makeBoundPres2D(std::vector<iisphparticle>& var_PartC);