#include "iisphparticle.h"
#include<omp.h>
#include <iostream>
#include <string>
void iisphparticle::resetvalues() {
	vel = glm::vec3(0,0,0), acc = glm::vec3(0, 0, 0), nonpresA = glm::vec3(0, 0, 0), predictedVel = glm::vec3(0, 0, 0), presA = glm::vec3(0, 0, 0);
	r = 44, g = 2, b = 25; // Color
	size = h, density = p0, pressure = 0, predictedDens = 0, densityerror = 0, pressureiter = 0, sf = 0, Aff = 0, AP = 0;
	m = (p0 * h * h * h );
	Kernel.clear();
	KernelDer.clear();
	IdNdistNsub.clear();
	IdKernelder.clear();
	IdNSubKernelder.clear();
}
void iisphparticle::makeKernel()
{
	Kernel.clear();
	for (const auto& neig : IdNdistNsub) {
		float d = std::get<1>(neig) / h;
		float t1 = std::max((1 - d), 0.f);
		float t2 = std::max((2 - d), 0.f);
		Kernel.push_back(alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
	}
}

void iisphparticle::makeKernelDer()
{
	IdNSubKernelder.clear();
	for (const auto& neig : IdNdistNsub) {
		if (std::get<1>(neig) == 0.) {
			IdNSubKernelder.push_back(std::make_tuple(std::get<0>(neig), std::get<2>(neig), glm::vec3(0, 0, 0)));
		}
		else {
			float d = std::get<1>(neig) / h;
			float t1 = std::max((1 - d), 0.f);
			float t2 = std::max((2 - d), 0.f);
			IdNSubKernelder.push_back(std::make_tuple(std::get<0>(neig), std::get<2>(neig), alpha * std::get<2>(neig) / (std::get<1>(neig)) * (-3 * t2 * t2 + 12 * t1 * t1)));
		}
	}
}

void iisphparticle::computeDens()
{
	float dens = 0;
	for (float kern : Kernel) {
		dens += kern * m;
	}
	if (isboundary == false) {
		density = std::max(dens, 0.f);
	}
	else
	{
		density = p0;
	}
}

void iisphparticle::predictVel()
{
	predictedVel = vel + (deltaT * nonpresA);
}

void iisphparticle::computeDensErr()
{
	densityerror = AP - sf;
}


void iisphparticle::updatePres()
{
	if (Aff != 0) {
		pressureiter = glm::max(pressureiter + (omega * (sf - AP) / Aff), 0.f);
	}
	else {
		pressureiter = glm::max(pressureiter, 0.f);
	}
}

void iisphparticle::makeV()
{
	vel += deltaT * (presA + nonpresA);
	float normv = glm::length(vel);
	r = 30;
	g = glm::min(30 + (cloloroffset * normv), 143.f);
	b = 148;
}

void iisphparticle::makeP()
{
	pos += vel * deltaT;
}







