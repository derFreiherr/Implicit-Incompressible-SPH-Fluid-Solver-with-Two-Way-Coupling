#include <part.h>
void Particle::makeKernel()
{
	Kernel.clear();
	for (std::tuple<int, float, glm::vec3>neig : IdNdistNsub) {
		float d = std::get<1>(neig) / h;
		float t1 = std::max((1 - d), 0.f);
		float t2 = std::max((2 - d), 0.f);
		Kernel.push_back(alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
	}
}
void Particle::makeKernelDer()
{
	IdNSubKernelder.clear();
	for (std::tuple<int, float, glm::vec3>neig : IdNdistNsub) {
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
void Particle::computeDens()
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
void Particle::makeDens()
{
	float dens = 0;
	for (float kern : Kernel) {
		dens += kern * m;
	}
	if (isboundary == false) {
		density = glm::max(dens, p0);
	}
	else
	{
		density = p0;
	}
}
void Particle::makePres()
{
	pressure = k * ((density / p0) - 1);
}
void Particle::makeA(Particle ParticlesContainer[])
{
	glm::vec3 PresA(0.f, 0.f, 0.f);
	glm::vec3 ViscA(0.f, 0.f, 0.f);
	for (std::tuple<int, glm::vec3, glm::vec3>neig : IdNSubKernelder) {
		if (ParticlesContainer[std::get<0>(neig)].isboundary) {
			PresA -= m * (pressure / (density * density) + pressure / (p0 * p0)) * std::get<2>(neig);
			ViscA += (m / p0) * (vel * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
		}
		else {
			PresA -= m * (pressure / (density * density) + ParticlesContainer[std::get<0>(neig)].pressure / (ParticlesContainer[std::get<0>(neig)].density * ParticlesContainer[std::get<0>(neig)].density)) * std::get<2>(neig);
			ViscA += (m / ParticlesContainer[std::get<0>(neig)].density) * ((vel - ParticlesContainer[std::get<0>(neig)].vel) * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
		}

	}
	ViscA = 2 * viscosity * ViscA + glm::vec3(0, -9.81, 0);
	acc = PresA + ViscA;
}
void Particle::makeNonpresA(Particle ParticlesContainer[])
{
	glm::vec3 ViscAf(0.f, 0.f, 0.f);
	glm::vec3 ViscAb(0.f, 0.f, 0.f);
	for (std::tuple<int, glm::vec3, glm::vec3>neig : IdNSubKernelder) {
		if (ParticlesContainer[std::get<0>(neig)].isboundary) {
			ViscAb += (ParticlesContainer[std::get<0>(neig)].m / p0) * (vel * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
		}
		else {
			ViscAf += (ParticlesContainer[std::get<0>(neig)].m / /*ParticlesContainer[std::get<0>(neig)].density*/p0) * ((vel - ParticlesContainer[std::get<0>(neig)].vel) * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
		}
	}
	nonpresA = 2 * viscosity * (ViscAb + ViscAf) + glm::vec3(0, -9.81, 0);
}
void Particle::computedivofvelchange(Particle PartC[]) {
	AP = 0;
	float APF = 0;
	float APB = 0;
	for (std::tuple<int, glm::vec3, glm::vec3>neig : IdNSubKernelder) {
		if (index != std::get<0>(neig)) {
			if (PartC[std::get<0>(neig)].isboundary) {
				APB += PartC[std::get<0>(neig)].m * glm::dot(presA, std::get<2>(neig));
			}
			else {
				APF += PartC[std::get<0>(neig)].m * glm::dot((presA - PartC[std::get<0>(neig)].presA), std::get<2>(neig));
			}
		}
	}
	AP = deltaT * deltaT * ((APB)+APF);
}
void Particle::makePresA(Particle PartC[]) {
	presA = glm::vec3(0, 0, 0);
	glm::vec3 PresAf(0.f, 0.f, 0.f);
	glm::vec3 PresAb(0.f, 0.f, 0.f);
	for (std::tuple<int, glm::vec3, glm::vec3>neig : IdNSubKernelder) {
		if (index != std::get<0>(neig)) {
			if (PartC[std::get<0>(neig)].isboundary) {
				PresAb -= PartC[std::get<0>(neig)].m * (2 * pressureiter / (p0 * p0)) * std::get<2>(neig);
			}
			else {
				PresAf -= PartC[std::get<0>(neig)].m * (pressureiter / (p0 * p0) + PartC[std::get<0>(neig)].pressureiter / (p0 * p0)) * std::get<2>(neig);
			}
		}
	}
	presA = (PresAf + (gammafloat * PresAb));
}
void Particle::updatePres() {
	if (Aff != 0) {
		pressureiter = glm::max(pressureiter + (omega * (sf - AP) / Aff), 0.f);
	}
	else {
		pressureiter = glm::max(pressureiter, 0.f);
	}
}
void Particle::makeV()
{
	vel += deltaT * acc;
	float normv = glm::length(vel);
	r = 30;
	g = glm::min(30 + (cloloroffset * normv), 143.f);
	b = 148;
}
void Particle::makeVIISPH() {
	if (glm::length(presA) > 100 * glm::length(nonpresA)) {
		presA = presA * 1.f / 1000.f;
	}
	vel += deltaT * (presA + nonpresA);
	float normv = glm::length(vel);
	r = 30;
	g = glm::min(30 + (cloloroffset * normv), 143.f);
	b = 148;
}
void Particle::predictVel()
{
	predictedVel = vel + (deltaT * nonpresA);
}
void Particle::computeDensErr() {
	densityerror = AP - sf;
}
void Particle::makeP() {
	pos += vel * deltaT;
}

