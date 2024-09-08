#include <helperz.h>
#include <iisphparticle.h>
#include <part.h>

void makecfltrue(std::vector<iisphparticle>& ParticlesContainer) {
	maxvel = 0;

	for (int i = 0; i < var_MaxParticles; i++) {
		if (glm::length(ParticlesContainer[i].vel) > maxvel && (ParticlesContainer[i].pos.y < (upperviualbord-h) && ParticlesContainer[i].pos.y > (lowervisualbord+h))) {
			maxvel = glm::length(ParticlesContainer[i].vel);
		}
	}
	if (fixeddt == false) {
		if (maxvel > 0) {
			deltaT = h * cfl_max / maxvel;
		}
		deltaT = glm::min(deltaT, deltaTmax);
	}
	cfl = maxvel * deltaT / h;
}


bool checkcfl(float maxvel) {
	cfl = maxvel * deltaT / h;
	if (cfl <= cfl_max) {
		cfl_cond = true;
	}
	else {
		cfl_cond = false;
	}
	return cfl_cond;
}


std::vector<glm::vec3> parseObjFile(const std::string& filePath) {
	std::vector<glm::vec3> vertices;
	std::ifstream file(filePath);
	if (!file.is_open()) {
		std::cerr << "Failed to open file: " << filePath << std::endl;
		return vertices;
	}

	std::string line;
	while (std::getline(file, line)) {
		if (line.substr(0, 2) == "v ") {
			std::istringstream s(line.substr(2));
			glm::vec3 vertex;
			s >> vertex.x; s >> vertex.y; s >> vertex.z;
			vertices.push_back(vertex);
		}
	}
	file.close();
	return vertices;
}
