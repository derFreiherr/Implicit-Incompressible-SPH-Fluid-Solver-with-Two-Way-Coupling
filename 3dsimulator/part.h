#include <vector>
#include <tuple>
#include <unordered_map>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include <chrono>
#include <Cell.h>
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
extern float k;
struct Particle {
	glm::vec3 pos, vel, acc, nonpresA, predictedVel, presA;
	unsigned char r = 44, g = 2, b = 25, a = 10; // Color
	float size = h, density, pressure = p0, predictedDens = 0, densityerror = 0, pressureiter = 0, sf = 0, Aff = 0, AP = 0;
	float m = (p0 * h * h * h);
	float cameradistance; // *Squared* distance to the camera. if dead : -1.0f
	std::vector<float>Kernel;
	std::vector<glm::vec3> KernelDer;
	std::vector<std::tuple<int, float, glm::vec3>> IdNdistNsub;
	std::vector<std::pair<int, glm::vec3>> IdKernelder;
	std::vector<std::tuple<int, glm::vec3, glm::vec3>> IdNSubKernelder;
	bool isboundary = 1;
	int index;
	void makeKernel();
	void makeKernelDer();
	void computeDens();
	void makeDens();
	void makePres();
	void makeA(Particle ParticlesContainer[]);
	void makeNonpresA(Particle ParticlesContainer[]);
	void computedivofvelchange(Particle PartC[]);
	void updatePres();
	void computeDensErr();
	void makePresA(Particle PartC[]);
	void predictVel();
	void makeV();
	void makeP();
	void makeVIISPH();
	bool operator<(const Particle& that) const {
		// Sort in reverse order : far particles drawn first.
		return this->cameradistance > that.cameradistance;
	}
};
