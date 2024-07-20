#include <vector>
#include <tuple>
#include <unordered_map>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include <chrono>
#include <Cell.h>
#include "Eigen/Dense"
using namespace Eigen;

extern const float p0;
extern const float h;
extern const int MaxParticles;
extern float viscosity;
extern float deltaT;
extern float cloloroffset;
extern const float alpha;
extern float gammafloat;
extern float omega;
extern float searchRadius;
extern float denistyerrormax;
extern int currentitermax;
extern int currentiter;
extern std::vector<float> densitysnew;
extern std::vector<float> densitys;
extern std::vector<float> currenttimes;
extern bool makesingle;
extern int boundarya;
extern int obstaclea;
extern float maxavgdensdeviation;
extern bool singlewall;
extern int watercolheight;
class iisphparticle
{
public:
	Matrix3d inertiaTensor;
	Matrix3d inertiaTensorInverse;
	Vector3d position;
	Quaterniond orientation;
	Vector3d linearMomentum;
	Vector3d angularMomentum;
	Vector3d velocity;
	Vector3d angularVelocity;
	Matrix3d rotationMatrix;
	double mass;
	bool isfloatingboundary = false;


	bool denstolow = false;
	bool computeme = true;
	bool drawme = false;
	
	glm::vec3 pos = glm::vec3(0, 0, 0);
	glm::vec3 vel = glm::vec3(0, 0, 0);
	glm::vec3 acc = glm::vec3(0, 0, 0);
	glm::vec3 nonpresA = glm::vec3(0, 0, 0);
	glm::vec3 predictedVel = glm::vec3(0, 0, 0);
	glm::vec3 presA = glm::vec3(0, 0, 0);
	unsigned char r = 44, g = 2, b = 25, a = 10; // Color
	float size = h, density = p0, pressure = 0, predictedDens = 0, densityerror = 0, pressureiter = 0, sf = 0, Aff = 0, AP = 0;
	float m = (p0 * h * h *h);
	//float cameradistance; // *Squared* distance to the camera. if dead : -1.0f
	std::vector<float>Kernel;
	std::vector<glm::vec3> KernelDer;
	std::vector<std::tuple<int, float, glm::vec3>> IdNdistNsub;
	std::vector<std::pair<int, glm::vec3>> IdKernelder;
	std::vector<std::tuple<int, glm::vec3, glm::vec3>> IdNSubKernelder;
	bool isboundary = 1;
	bool ismovingboundary = false;
	int index;
	void resetvalues();
	void makeKernel();
	void makeKernelDer();
	void computeDens();
	void predictVel();
	void computeDensErr();
	void updatePres();
	void makeV();
	void makeP();
	/*
	bool operator<(const iisphparticle& that) const {
		// Sort in reverse order : far particles drawn first.
		return this->cameradistance > that.cameradistance;
	}
	*/
};
extern void var_iisph(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap);
extern void twoDiisph(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap);


