#include <Cell.h>
#include <iisphparticle.h>
#include <part.h>
unsigned int hashFunction3D(float x, float y, float z)
{
	int ix = static_cast<int>(std::floor(x / cellsize));
	int iy = static_cast<int>(std::floor(y / cellsize));
	int iz = static_cast<int>(std::floor(z / cellsize));
	return ((ix * p1) ^ (iy * p2) ^ (iz * p3)) % hashsize;
}

unsigned int hashFunction2D(float x, float y)
{
	int ix = static_cast<int>(std::floor(x / cellsize));
	int iy = static_cast<int>(std::floor(y / cellsize));
	return ((ix * p1) ^ (iy * p2)) % hashsize;
}

unsigned int uniformgridhash(float x, float y, float z) {
	float xi = (2+x) / cellsize;
	float yi = (100+y) / cellsize;
	float zi = (100+z) / cellsize;
	return (xi + yi * gridbreite + zi * gridhöhe);
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

void makecfltrue(std::vector<iisphparticle>& ParticlesContainer) {
	maxvel = 0;

	for (int i = 0; i < var_MaxParticles; i++) {
		if (glm::length(ParticlesContainer[i].vel) > maxvel && ParticlesContainer[i].pos.y < upperviualbord-1 && ParticlesContainer[i].pos.y > lowervisualbord+1) {
			maxvel = glm::length(ParticlesContainer[i].vel);
		}
	}
	if (fixeddt == false) {
		if (maxvel > 0) {
			deltaT = h * cfl_max / maxvel;
		}
		deltaT = glm::min(deltaT, deltaTmax);
	}
	/*
	if (exportanimation) {
		if (animationtime + deltaT > 0.039998) {
			float dtold = deltaT;
			deltaT = 0.040002 - animationtime;
			if (deltaT < 0.002) {
				animationtime = 0.04;
				deltaT = dtold;
			}
		}
	}
	*/
	cfl = maxvel * deltaT / h;
}
void ssphAlgo(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap) {
	auto start_o = std::chrono::high_resolution_clock::now();
	makecfltrue(PartC);
	//deltaT = glm::min(deltaT, deltaTmax);
	auto start_n = std::chrono::high_resolution_clock::now();
	hashmap.clear();
	insertAllParticlesIntoHashmap(PartC, hashmap);
	//findneighbours
	findAllNeighbours(PartC, hashmap);
	auto end_n = std::chrono::high_resolution_clock::now();
	//makeKernel and Kernelder
	makeAllKernelAndKernelDer(PartC);
	auto end_kern = std::chrono::high_resolution_clock::now();
	//compute density
	computeAllDensSSPH(PartC);
	auto end_dens = std::chrono::high_resolution_clock::now();
	//compute all pres
	computeAllPres(PartC);
	//make acc
	makeAllA(PartC);
	//update velocity and Position
	makeAllVandP(PartC);
	updaterigidbody(PartC);
	auto end_o = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_o - start_o);
	othercomputationTime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_n - start_n);
	neighSearchTime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_kern - end_n);
	kerneltime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_dens - end_kern);
	densitytime = duration.count();
	densitys.push_back(0.f);
	densitysnew.push_back(denserrold);
}
void ssphAlgotwoD(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap) {
	auto start_o = std::chrono::high_resolution_clock::now();
	makecfltrue(PartC);
	//deltaT = glm::min(deltaT, deltaTmax);
	auto start_n = std::chrono::high_resolution_clock::now();
	hashmap.clear();
	insertAllParticlesIntoHashmap2D(PartC, hashmap);
	//findneighbours
	findAllNeighbours2D(PartC, hashmap);
	auto end_n = std::chrono::high_resolution_clock::now();
	//makeKernel and Kernelder
	makeAllKernelAndKernelDerTwoD(PartC);
	auto end_kern = std::chrono::high_resolution_clock::now();
	//compute density
	computeAllDensSSPH(PartC);
	auto end_dens = std::chrono::high_resolution_clock::now();
	//compute all pres
	computeAllPres(PartC);
	//make acc
	makeAllA(PartC);
	//update velocity and Position
	makeAllVandP(PartC);
	auto end_o = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_o - start_o);
	othercomputationTime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_n - start_n);
	neighSearchTime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_kern - end_n);
	kerneltime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_dens - end_kern);
	densitytime = duration.count();
	densitys.push_back(0.f);
	densitysnew.push_back(denserrold);
}
void computeAllPres(std::vector<iisphparticle>&  var_PartC) {
#pragma omp parallel for
for (int i = 0; i < var_PartC.size(); ++i) {
	iisphparticle& Part = var_PartC[i];
		if (var_PartC[i].isboundary == false || var_PartC[i].isfloatingboundary) {
			var_PartC[i].pressure = k * ((var_PartC[i].density / p0) - 1);
		}
	}
}

void makeAllA(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
for (int i = 0; i < var_PartC.size(); ++i) {
	iisphparticle& Part = var_PartC[i];
		if (Part.isboundary == false || Part.isfloatingboundary) {
			Part.presA = glm::vec3(0, 0, 0);
			glm::vec3 PresAf(0.f, 0.f, 0.f);
			glm::vec3 PresAb(0.f, 0.f, 0.f);
			glm::vec3 ViscAf(0.f, 0.f, 0.f);
			glm::vec3 ViscAb(0.f, 0.f, 0.f);
			for (std::tuple<int, glm::vec3, glm::vec3>& neig : Part.IdNSubKernelder) {
				if (Part.index != std::get<0>(neig)) {
					if (var_PartC[std::get<0>(neig)].isboundary && !var_PartC[std::get<0>(neig)].ismovingboundary) {
						PresAb -= var_PartC[std::get<0>(neig)].m * (2 * Part.pressure / (p0p0)) * std::get<2>(neig);
						ViscAb += (var_PartC[std::get<0>(neig)].m / p0) * (Part.vel * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
					}
					else {
						PresAf -= var_PartC[std::get<0>(neig)].m * (Part.pressure / (p0p0) + var_PartC[std::get<0>(neig)].pressure / (p0p0)) * std::get<2>(neig);
						ViscAf += (var_PartC[std::get<0>(neig)].m / var_PartC[std::get<0>(neig)].density) * ((Part.vel - var_PartC[std::get<0>(neig)].vel) * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
					}
				}
			}
			Part.presA = (PresAf + (gammapres * PresAb));
			Part.nonpresA = 2.f * (visc_boundary * ViscAb + visc_fluid * ViscAf) + gravity * glm::vec3(0, 1, 0);
		}
	}
}

void makeAllAtwoD(std::vector<iisphparticle>& var_PartC) {
//#pragma omp parallel for
	for (iisphparticle& Part : var_PartC) {
		if (Part.isboundary == false) {
			Part.presA = glm::vec3(0, 0, 0);
			glm::vec3 PresAf(0.f, 0.f, 0.f);
			glm::vec3 PresAb(0.f, 0.f, 0.f);
			glm::vec3 ViscAf(0.f, 0.f, 0.f);
			glm::vec3 ViscAb(0.f, 0.f, 0.f);
			for (std::tuple<int, glm::vec3, glm::vec3>& neig : Part.IdNSubKernelder) {
				if (Part.index != std::get<0>(neig)) {
					if (var_PartC[std::get<0>(neig)].isboundary && !var_PartC[std::get<0>(neig)].ismovingboundary) {
						PresAb -= var_PartC[std::get<0>(neig)].m * (2 * Part.pressure / (p0p0)) * std::get<2>(neig);
						ViscAb += (var_PartC[std::get<0>(neig)].m / p0) * (Part.vel * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
					}
					else {
						PresAf -= var_PartC[std::get<0>(neig)].m * (Part.pressure / (p0p0) + var_PartC[std::get<0>(neig)].pressure / (p0p0)) * std::get<2>(neig);
						ViscAf += (var_PartC[std::get<0>(neig)].m / var_PartC[std::get<0>(neig)].density) * ((Part.vel - var_PartC[std::get<0>(neig)].vel) * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
					}
				}
			}
			Part.presA = (PresAf + (gammapres * PresAb)) * glm::vec3(1, 1, 0);
			Part.nonpresA = ( (visc_boundary * ViscAb + visc_fluid * ViscAf) + gravity * glm::vec3(0, 1, 0))* glm::vec3(1, 1, 0);
		}
	}
}
void PredictAllVel(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isboundary == false || Part.isfloatingboundary) {
			Part.predictedVel = Part.vel + (deltaT * Part.nonpresA);
		}
	}
}
void MakeAllNonpresA(std::vector<iisphparticle>& var_PartC ) {
#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if ((Part.isboundary == false || Part.isfloatingboundary)) {
			glm::vec3 ViscAf(0.f, 0.f, 0.f);
			glm::vec3 ViscAb(0.f, 0.f, 0.f);
			glm::vec3 surfacetens(0.f, 0.f, 0.f);
			if (Part.isfloatingboundary) {
				for (const auto& neig : Part.IdNSubKernelder) {
					if (var_PartC[std::get<0>(neig)].isboundary && !var_PartC[std::get<0>(neig)].ismovingboundary) {
						ViscAb += (var_PartC[std::get<0>(neig)].m / p0) * (Part.vel * makesinglekernel(Part.pos, var_PartC[std::get<0>(neig)].pos));
					}
					else {
						ViscAf += (var_PartC[std::get<0>(neig)].m / p0) * ((Part.vel - var_PartC[std::get<0>(neig)].vel) * makesinglekernel(Part.pos, var_PartC[std::get<0>(neig)].pos));
					}
				}
				Part.nonpresA = (1.f * ViscAb + 1.f * ViscAf) + gravity * glm::vec3(0, 1, 0);
			}
			else {
				for (const auto& neig : Part.IdNSubKernelder) {
					if (var_PartC[std::get<0>(neig)].isboundary && !var_PartC[std::get<0>(neig)].ismovingboundary) {
						ViscAb += (var_PartC[std::get<0>(neig)].m / p0) * (Part.vel * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
					}
					else {
						ViscAf += (var_PartC[std::get<0>(neig)].m / /*ParticlesContainer[std::get<0>(neig)].density*/p0) * ((Part.vel - var_PartC[std::get<0>(neig)].vel) * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
					}
				}
				for (const auto& neig : Part.IdNdistNsub) {
					if (var_PartC[std::get<0>(neig)].isboundary == false && var_PartC[std::get<0>(neig)].isfloatingboundary == false) {
						float kern = makesinglekernel(Part.pos, var_PartC[std::get<0>(neig)].pos);
						surfacetens += var_PartC[std::get<0>(neig)].m * std::get<2>(neig) * kern;
					}
				}
				Part.nonpresA = (visc_boundary * ViscAb + visc_fluid * ViscAf) + gravity * glm::vec3(0, 1, 0) - ((surfacetension / Part.m) * surfacetens);
			}
			
			
			
			
			
		}
	}
}
void MakeAllNonpresAtwoD(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < (var_fluidpart+var_spezialboundpart); ++i) {
		iisphparticle& Part = var_PartC[i];
		glm::vec3 ViscAf(0.f, 0.f, 0.f);
		glm::vec3 ViscAb(0.f, 0.f, 0.f);
		glm::vec3 surfacetens(0.f, 0.f, 0.f);
		for (const auto& neig : Part.IdNSubKernelder) {
			if (var_PartC[std::get<0>(neig)].isboundary && !var_PartC[std::get<0>(neig)].ismovingboundary && !!var_PartC[std::get<0>(neig)].isfloatingboundary) {
				ViscAb += (var_PartC[std::get<0>(neig)].m / p0) * (Part.vel * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
			}
			else {
				ViscAf += (var_PartC[std::get<0>(neig)].m / /*ParticlesContainer[std::get<0>(neig)].density*/p0) * ((Part.vel - var_PartC[std::get<0>(neig)].vel) * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
			}
		}
		for (const auto& neig : Part.IdNdistNsub) {
			if (var_PartC[std::get<0>(neig)].isboundary == false && Part.isfloatingboundary == false && var_PartC[std::get<0>(neig)].isfloatingboundary == false) {
				float kern = makesinglekernel2D(Part.pos, var_PartC[std::get<0>(neig)].pos);
				surfacetens += var_PartC[std::get<0>(neig)].m * std::get<2>(neig) * kern;
			}
		}
		Part.nonpresA = (visc_boundary * ViscAb + visc_fluid * ViscAf) + gravity * glm::vec3(0, 1, 0) - ((surfacetension / Part.m) * surfacetens) * glm::vec3(1.f, 1.f, 0.f);
	}
}
void MakeAllNonpresAwithdens(std::vector<iisphparticle>& var_PartC) {
//#pragma omp parallel for
	for (iisphparticle& Part : var_PartC) {
		if (Part.isboundary == false) {
			glm::vec3 ViscAf(0.f, 0.f, 0.f);
			glm::vec3 ViscAb(0.f, 0.f, 0.f);
			glm::vec3 surfacetens(0.f, 0.f, 0.f);
			for (const auto& neig : Part.IdNSubKernelder) {
				if (var_PartC[std::get<0>(neig)].isboundary && !var_PartC[std::get<0>(neig)].ismovingboundary) {
					ViscAb += (var_PartC[std::get<0>(neig)].m / p0) * (Part.vel * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
				}
				else {
					ViscAf += (var_PartC[std::get<0>(neig)].m / var_PartC[std::get<0>(neig)].density) * ((Part.vel - var_PartC[std::get<0>(neig)].vel) * std::get<1>(neig) / (std::get<1>(neig) * std::get<1>(neig) + 0.01f * h * h)) * std::get<2>(neig);
				}
			}
			for (const auto& neig : Part.IdNdistNsub) {
				if (var_PartC[std::get<0>(neig)].isboundary == false && Part.computeme == false) {
					float d = std::get<1>(neig) / h;
					float t1 = std::max((1 - d), 0.f);
					float t2 = std::max((2 - d), 0.f);
					surfacetens += var_PartC[std::get<0>(neig)].m * std::get<1>(neig) * (alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				}
			}
			Part.nonpresA = 2.f * (visc_boundary * ViscAb + visc_fluid * ViscAf) + gravity * glm::vec3(0, 1, 0) - ((surfacetension / Part.m) * surfacetens);
		}
	}
}
void computeAllDens(std::vector<iisphparticle>& var_PartC) {
	denserrold = 0;
	numofp0high = 0;
	usemefordens.clear();
	//float deabugdens = 0;
//#pragma omp parallel for reduction(+:denserrold, numofp0high)
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		Part.density = 0;
		if (Part.isboundary == false || Part.isfloatingboundary) {
			Part.denstolow = true;
			float dens = 0;
			/*
			for (float& kern : Part.Kernel) {
				dens += kern * Part.m;
			}
			*/
			for (const auto& neig : Part.IdNdistNsub) {
				dens += var_PartC[std::get<0>(neig)].m * makesinglekernel(Part.pos, var_PartC[std::get<0>(neig)].pos);
			}
			Part.density = dens;
			if (dens >= p0 && Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord &&Part.isboundary == false && Part.isfloatingboundary == false) {
				Part.denstolow = false;
				numofp0high++;
				if (absinterrupt) {
					denserrold += 100*(std::abs(Part.density - p0)/p0);
					usemefordens.push_back(Part.index);
					//deabugdens = std::max(deabugdens, Part.density);

				}
				else {
					denserrold += 100*(Part.density - p0)/p0;
					usemefordens.push_back(Part.index);
				}
			}

			/*
			if (clampp0 == true) {
				Part.density = std::max(dens, p0 - p0 * denistyerrormax / 100 * clampfac);
				Part.denstolow = false;
				numofp0high++;
				if (absinterrupt) {
					denserrold += std::abs(Part.density - p0);
				}
				else {
					denserrold += Part.density - p0;
				}
			}
			if (Part.denstolow == true && clamptolow == true) {
				Part.density = std::max(dens, p0 - p0 * denistyerrormax / 100 * clampfac);
			}
			*/
		}
		if (Part.isboundary == true && Part.isfloatingboundary == false) {
			Part.density = p0;
			Part.denstolow = true;
		}
		
	}
	//std::cout << deabugdens << "   "<< usemefordens.size()<< std::endl;
	denserrold = (denserrold / usemefordens.size());
	//denserrold = (denserrold / var_fluidpart);
}
void computeAllDenstwoD(std::vector<iisphparticle>& var_PartC) {
	denserrold = 0;
	numofp0high = 0;
	usemefordens.clear();
	//float deabugdens = 0;
//#pragma omp parallel for reduction(+:denserrold, numofp0high)
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		Part.density = 0;
		if (Part.isboundary == false || Part.isfloatingboundary) {
			Part.denstolow = true;
			float dens = 0;
			/*
			for (float& kern : Part.Kernel) {
				dens += kern * Part.m;
			}
			*/
			for (const auto& neig : Part.IdNdistNsub) {
				dens += var_PartC[std::get<0>(neig)].m * makesinglekernel2D(Part.pos, var_PartC[std::get<0>(neig)].pos);
			}
			Part.density = dens;
			if (dens >= p0 && Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord && Part.isboundary == false && Part.isfloatingboundary == false) {
				Part.denstolow = false;
				numofp0high++;
				if (absinterrupt) {
					denserrold += 100 * (std::abs(Part.density - p0) / p0);
					usemefordens.push_back(Part.index);
					//deabugdens = std::max(deabugdens, Part.density);

				}
				else {
					denserrold += 100 * (Part.density - p0) / p0;
					usemefordens.push_back(Part.index);
				}
			}

			/*
			if (clampp0 == true) {
				Part.density = std::max(dens, p0 - p0 * denistyerrormax / 100 * clampfac);
				Part.denstolow = false;
				numofp0high++;
				if (absinterrupt) {
					denserrold += std::abs(Part.density - p0);
				}
				else {
					denserrold += Part.density - p0;
				}
			}
			if (Part.denstolow == true && clamptolow == true) {
				Part.density = std::max(dens, p0 - p0 * denistyerrormax / 100 * clampfac);
			}
			*/
		}
		if (Part.isboundary == true && Part.isfloatingboundary == false) {
			Part.density = p0;
			Part.denstolow = true;
		}

	}
	//std::cout << deabugdens << "   "<< usemefordens.size()<< std::endl;
	denserrold = (denserrold / usemefordens.size());
	//denserrold = (denserrold / var_fluidpart);
}
void computeAllDensSSPH(std::vector<iisphparticle>& var_PartC) {
	denserrold = 0;
	numofp0high = 0;
//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isboundary == false || Part.isfloatingboundary) {
			float dens = 0;
			for (float& kern : Part.Kernel) {
				dens += kern * Part.m;
			}
			dens = std::max(dens, p0);
			Part.density = dens;
			Part.denstolow = false;
			numofp0high++;
			if (absinterrupt) {
				denserrold += std::abs(Part.density - p0);
			}
			else {
				denserrold += Part.density - p0;
			}
		}
		if (Part.isboundary == true && Part.isfloatingboundary == false) {
			Part.density = p0;
		}
	}
	denserrold = 100 * (denserrold / numofp0high) / p0;
}
void makeAllKernelAndKernelDer(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isboundary == false || Part.isboundary == true) {
			Part.IdNSubKernelder.clear();
			for (const auto& neig : Part.IdNdistNsub) {
				if (std::get<1>(neig) == 0.) {
					Part.IdNSubKernelder.push_back(std::make_tuple(std::get<0>(neig), std::get<2>(neig), glm::vec3(0, 0, 0)));
				}
				else {
					float d = std::get<1>(neig) / h;
					float t1 = std::max((1 - d), 0.f);
					float t2 = std::max((2 - d), 0.f);
					Part.IdNSubKernelder.push_back(std::make_tuple(std::get<0>(neig), std::get<2>(neig), 0.999221087f* alpha * std::get<2>(neig) / (d*h*h) * (-3 * t2 * t2 + 12 * t1 * t1)));
				}
			}
			Part.Kernel.clear();
//#pragma omp parallel for
			for (const auto& neig : Part.IdNdistNsub) {
				float d = std::get<1>(neig) / h;
				float t1 = std::max((1 - d), 0.f);
				float t2 = std::max((2 - d), 0.f);
				if (var_PartC[std::get<0>(neig)].isboundary) {
					Part.Kernel.push_back(0.999221087*gammadens * alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				}
				else {
					Part.Kernel.push_back(0.999221087*alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				}
			}
		}
		
	}
}
void makeAllKernelAndKernelDerTwoD(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isboundary == false || Part.isboundary == true) {
			Part.IdNSubKernelder.clear();
			for (const auto& neig : Part.IdNdistNsub) {
				if (std::get<1>(neig) == 0.) {
					Part.IdNSubKernelder.push_back(std::make_tuple(std::get<0>(neig), std::get<2>(neig) * glm::vec3(1.f, 1.f, 0), glm::vec3(0, 0, 0)));
				}
				else {
					float d = std::get<1>(neig) / h;
					float t1 = std::max((1 - d), 0.f);
					float t2 = std::max((2 - d), 0.f);
					Part.IdNSubKernelder.push_back(std::make_tuple(std::get<0>(neig), std::get<2>(neig) * glm::vec3(1.f, 1.f, 0), 0.99911389780816045467516953288827f * alphaTwoD * std::get<2>(neig) / (d * h * h) * (-3 * t2 * t2 + 12 * t1 * t1)));
				}
			}
			Part.Kernel.clear();
			//#pragma omp parallel for
			for (const auto& neig : Part.IdNdistNsub) {
				float d = std::get<1>(neig) / h;
				float t1 = std::max((1 - d), 0.f);
				float t2 = std::max((2 - d), 0.f);
				if (var_PartC[std::get<0>(neig)].isboundary) {
					Part.Kernel.push_back(0.9991389780816045467516953288827 * gammadens * alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				}
				else {
					Part.Kernel.push_back(0.9991389780816045467516953288827 * alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				}

			}
		}
	}
}
/*
void findAllNeighbours(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	totalcomp = 0;
//#pragma omp parallel for schedule(auto)
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if (Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord) {
			Part.IdNdistNsub.clear();
			Part.drawme = false;
			Part.computeme = false;
			Part.Aff = 0;
			std::unordered_set<int> uniqueIds; // Set to track unique IDs
			for (int x = -1; x <= 1; x++) {
				for (int y = -1; y <= 1; y++) {
					for (int z = -1; z <= 1; z++) {
						int hash = hashFunction3D(Part.pos.x + x * searchRadius, Part.pos.y + y * searchRadius, Part.pos.z + z * searchRadius);
						if (hashmap.find(hash) != hashmap.end()) {
							for (int idx : hashmap[hash].particles) {
								if (uniqueIds.find(idx) == uniqueIds.end()) { // Check if ID is already added
									float distX = var_PartC[idx].pos.x - Part.pos.x;
									float distY = var_PartC[idx].pos.y - Part.pos.y;
									float distZ = var_PartC[idx].pos.z - Part.pos.z;
									float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
									if (distance < searchRadius) {
										Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
										uniqueIds.insert(idx); // Mark this ID as added
										if (var_PartC[idx].isboundary && Part.isboundary == false) {
											Part.drawme = true;
										}
									}
								}
							}
						}
					}
				}
			}
			if (Part.IdNdistNsub.size() < animationneighbourscount && Part.isboundary == false) {
				Part.drawme = true;
			}
			if (Part.IdNdistNsub.size() > compbord && Part.isboundary == false) {
				Part.computeme = true;
				totalcomp++;
			}
		}
	}
}
*/
void findAllNeighbours(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	totalcomp = 0;
#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		Part.IdNdistNsub.clear();
		Part.drawme = false;
		Part.computeme = false;
		Part.Aff = 0;
		std::unordered_set<int> uniqueIds; // Set to track unique IDs

		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				for (int z = -1; z <= 1; z++) {
					int hash = hashFunction3D(Part.pos.x + x * searchRadius, Part.pos.y + y * searchRadius, Part.pos.z + z * searchRadius);
					if (hashmap.find(hash) != hashmap.end()) {
						for (int idx : hashmap[hash].particles) {
							if (uniqueIds.find(idx) == uniqueIds.end()) { // Check if ID is already added
								/*
								float distX = var_PartC[idx].pos.x - Part.pos.x;
								float distY = var_PartC[idx].pos.y - Part.pos.y;
								float distZ = var_PartC[idx].pos.z - Part.pos.z;
								float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
								*/
								float distance = glm::distance(var_PartC[idx].pos, Part.pos);
								if (distance <= searchRadius) {
									Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
									uniqueIds.insert(idx); // Mark this ID as added
									//if (var_PartC[idx].isboundary && Part.isboundary == false) {
									//	Part.drawme = true;
									//}

								}
							}
						}
					}
				}
			}
		}
		if (Part.IdNdistNsub.size() < animationneighbourscount) {
			Part.drawme = true;
		}
		if (Part.IdNdistNsub.size() > compbord) {
			Part.computeme = true;
			totalcomp++;
		}
	}
}




void findAllNeighbours2D(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	totalcomp = 0;
#pragma omp parallel for
	for (int i = 0; i < (var_fluidpart + var_spezialboundpart); ++i) {
		iisphparticle& Part = var_PartC[i];
		Part.IdNdistNsub.clear();
		Part.drawme = false;
		Part.computeme = false;
		Part.Aff = 0;
		std::unordered_set<int> uniqueIds; // Set to track unique IDs
		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				int hash = hashFunction2D(Part.pos.x + x * searchRadius, Part.pos.y + y * searchRadius);
				if (hashmap.find(hash) != hashmap.end()) {
					for (int idx : hashmap[hash].particles) {
						if (uniqueIds.find(idx) == uniqueIds.end()) { // Check if ID is already added
							float distX = var_PartC[idx].pos.x - Part.pos.x;
							float distY = var_PartC[idx].pos.y - Part.pos.y;
							float distance = sqrt(distX * distX + distY * distY);
							if (distance <= searchRadius) {

								Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
								uniqueIds.insert(idx); // Mark this ID as added
								//if (var_PartC[idx].isboundary && Part.isboundary == false) {
								//	Part.drawme = true;
								//}

							}
						}
					}
				}
			}
		}
		if (Part.IdNdistNsub.size() < animationneighbourscount) {
			Part.drawme = true;
		}
		if (Part.IdNdistNsub.size() > compbord) {
			Part.computeme = true;
			totalcomp++;
		}
		
	}
}

void findAllNeighbourscompact2D(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	totalcomp = 0;
#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		Part.IdNdistNsub.clear();
		Part.drawme = false;
		Part.computeme = false;
		Part.Aff = 0;
		std::unordered_set<int> uniqueIds; // Set to track unique IDs
		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				int hash = uniformgridhash(Part.pos.x + x * searchRadius, Part.pos.y + y * searchRadius, Part.pos.z + 0 * searchRadius);
				if (hashmap.find(hash) != hashmap.end()) {
					for (int idx : hashmap[hash].particles) {
						if (uniqueIds.find(idx) == uniqueIds.end()) { // Check if ID is already added
							float distX = var_PartC[idx].pos.x - Part.pos.x;
							float distY = var_PartC[idx].pos.y - Part.pos.y;
							float distance = sqrt(distX * distX + distY * distY);
							if (distance < searchRadius) {
								Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
								uniqueIds.insert(idx); // Mark this ID as added
								if (var_PartC[idx].isboundary && Part.isboundary == false) {
									Part.drawme = true;
								}
							}
						}
					}
				}
			}
		}
		if (Part.IdNdistNsub.size() < animationneighbourscount && Part.isboundary == false) {
			Part.drawme = true;
		}
		if (Part.IdNdistNsub.size() > compbord && Part.isboundary == false) {
			Part.computeme = true;
			totalcomp++;
		}
	}
}



void insertAllParticlesIntoHashmap(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
#pragma omp parallel for
	for (int i = 0; i < var_MaxParticles; ++i) {
		if (var_PartC[i].pos.y < upperviualbord && var_PartC[i].pos.y > lowervisualbord) {
			int hash = hashFunction3D(var_PartC[i].pos.x, var_PartC[i].pos.y, var_PartC[i].pos.z);
#pragma omp critical
			{
				hashmap[hash].particles.push_back(i);
			}
		}
	}
}
void insertAllParticlesIntoHashmap2D(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
#pragma omp parallel for
	for (int i = 0; i < var_MaxParticles; ++i) {
		if (var_PartC[i].pos.y < upperviualbord && var_PartC[i].pos.y > lowervisualbord) {
			int hash = hashFunction2D(var_PartC[i].pos.x, var_PartC[i].pos.y);
#pragma omp critical
			{
				hashmap[hash].particles.push_back(i);
			}
		}
	}
}
void computeAllSF(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < var_MaxParticles; i++) {
		if (var_PartC[i].isboundary == false) {
			float tmpf = 0;
			float tmpb = 0.f;
			var_PartC[i].sf = 0.f;
			for (std::tuple<int, glm::vec3, glm::vec3>&neig : var_PartC[i].IdNSubKernelder) {
				if (var_PartC[std::get<0>(neig)].isboundary && var_PartC[std::get<0>(neig)].index != var_PartC[i].index) {
					tmpb += var_PartC[std::get<0>(neig)].m * glm::dot(var_PartC[i].predictedVel - var_PartC[std::get<0>(neig)].vel, std::get<2>(neig));
				}
				if (var_PartC[std::get<0>(neig)].index != var_PartC[i].index && var_PartC[std::get<0>(neig)].isboundary == false) {
					tmpf += var_PartC[std::get<0>(neig)].m * glm::dot((var_PartC[i].predictedVel - var_PartC[std::get<0>(neig)].predictedVel), std::get<2>(neig));
				}
			}
			var_PartC[i].sf = p0 - var_PartC[i].density - (deltaT * tmpf) - (deltaT * tmpb);
		}
	}
}
void makeAllAff1(std::vector<iisphparticle>& PartC) {
	//big sum over all fluidpart i
//#pragma omp parallel for
	for (int i = 0; i < var_MaxParticles; i++) {
		if (PartC[i].isboundary == false) {
			PartC[i].Aff = 0.f;
			float if_jf1 = 0;
			float if_jf2 = 0;
			float if_jb = 0;
			// sum over fluid neighbours jf from if exept if
			for (std::tuple<int, glm::vec3, glm::vec3>&jf : PartC[i].IdNSubKernelder) {
				if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != PartC[i].index) {
					glm::vec3  jjf = glm::vec3(0.f, 0.f, 0.f);
					glm::vec3  jjb = glm::vec3(0.f, 0.f, 0.f);
					// sum over all neigbours jj from if exept if
					for (std::tuple<int, glm::vec3, glm::vec3>&jj : PartC[i].IdNSubKernelder) {
						// fluidneighbours jjf
						if (PartC[std::get<0>(jj)].isboundary == false && PartC[std::get<0>(jj)].index != PartC[i].index) {
							//sum( jjf_m /p0*p0 * dWdt_jf_jjf)
							jjf += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
						// boundaryneighbours jjb
						if (PartC[std::get<0>(jj)].isboundary) {
							//sum( jjb_m /p0*p0 * dWdt_jf_jjb)
							jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
					}
					if_jf1 += PartC[std::get<0>(jf)].m * glm::dot((-jjf - 2 * gammapres * jjb), std::get<2>(jf));
				}
			}
			// sum over boundary neighbours jb from if
			for (std::tuple<int, glm::vec3, glm::vec3>&jb : PartC[i].IdNSubKernelder) {
				if (PartC[std::get<0>(jb)].isboundary) {
					glm::vec3  jjf = glm::vec3(0.f, 0.f, 0.f);
					glm::vec3  jjb = glm::vec3(0.f, 0.f, 0.f);
					// sum over all neigbours jj from if
					for (std::tuple<int, glm::vec3, glm::vec3>&jj : PartC[i].IdNSubKernelder) {
						// fluidneighbours jjf
						if (PartC[std::get<0>(jj)].isboundary == false) {
							//sum( jjf_m /p0*p0 * dWdt_jb_jjf)
							jjf += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
						// boundaryneighbours jjb
						if (PartC[std::get<0>(jj)].isboundary) {
							//sum( jjb_m /p0*p0 * dWdt_jb_jjb)
							jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
					}
					if_jb += glm::dot(PartC[std::get<0>(jb)].m * (-jjf - 2 * gammapres * jjb), std::get<2>(jb));
				}
			}
			// sum over fluid neighbours jf from if
			for (std::tuple<int, glm::vec3, glm::vec3>&jf : PartC[i].IdNSubKernelder) {
				if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != PartC[i].index) {
					glm::vec3 Kerndelder_jf_if = glm::vec3(0, 0, 0);
					// sum over all neigbours jj from jf
					for (std::tuple<int, glm::vec3, glm::vec3>&tmp : PartC[std::get<0>(jf)].IdNSubKernelder) {
						// if jj == if
						if (PartC[std::get<0>(tmp)].index == PartC[i].index) {
							Kerndelder_jf_if = std::get<2>(tmp);
							// sum ( jf_m * (if_m/(p0*p0) * dWdt_jf_if)* dWdt_if_jf
							if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((PartC[i].m / (p0p0) * Kerndelder_jf_if), std::get<2>(jf));
						}
					}
				}
			}
			PartC[i].Aff = (deltaT * deltaT) * (if_jf1 + if_jf2 + if_jb);
		}
	}
}
void makeAllAfffast(std::vector<iisphparticle>& PartC) {
	glm::vec3 Kerndelder_jf_if = glm::vec3(0, 0, 0);
	glm::vec3  jjf = glm::vec3(0.f, 0.f, 0.f);
	glm::vec3  jjb = glm::vec3(0.f, 0.f, 0.f);
	//big sum over all fluidpart i
//#pragma omp parallel for
	for (int i = 0; i < var_MaxParticles; i++) {
		if (PartC[i].isboundary == false) {
			PartC[i].Aff = 0.f;
			float if_jf1 = 0;
			float if_jf2 = 0;
			float if_jb = 0;
			auto& neighbors = PartC[i].IdNSubKernelder;
			// sum over neighbours jf from if exept if
			for (std::tuple<int, glm::vec3, glm::vec3>& jf : neighbors) {
				if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != PartC[i].index) {
					jjf = glm::vec3(0.f, 0.f, 0.f);
					jjb = glm::vec3(0.f, 0.f, 0.f);
					// sum over all neigbours jj from if exept if
					for (std::tuple<int, glm::vec3, glm::vec3>& jj : neighbors) {
						// fluidneighbours jjf
						if (PartC[std::get<0>(jj)].isboundary == false && PartC[std::get<0>(jj)].index != PartC[i].index) {
							//sum( jjf_m /p0*p0 * dWdt_jf_jjf)
							jjf += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
						// boundaryneighbours jjb
						if (PartC[std::get<0>(jj)].isboundary) {
							//sum( jjb_m /p0*p0 * dWdt_jf_jjb)
							jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
					}
					if_jf1 += PartC[std::get<0>(jf)].m * glm::dot((-jjf - 2 * gammapres * jjb), std::get<2>(jf));
				}
				if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != PartC[i].index) {
					Kerndelder_jf_if = glm::vec3(0, 0, 0);
					// sum over all neigbours jj from jf
					for (std::tuple<int, glm::vec3, glm::vec3>& tmp : PartC[std::get<0>(jf)].IdNSubKernelder) {
						// if jj == if
						if (PartC[std::get<0>(tmp)].index == PartC[i].index) {
							Kerndelder_jf_if = std::get<2>(tmp);
							// sum ( jf_m * (if_m/(p0*p0) * dWdt_jf_if)* dWdt_if_jf
							if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((PartC[i].m / (p0p0) * Kerndelder_jf_if), std::get<2>(jf));
						}
					}
				}
				if (PartC[std::get<0>(jf)].isboundary) {
					jjf = glm::vec3(0.f, 0.f, 0.f);
					jjb = glm::vec3(0.f, 0.f, 0.f);
					// sum over all neigbours jj from if
					for (std::tuple<int, glm::vec3, glm::vec3>& jj : neighbors) {
						// fluidneighbours jjf
						if (PartC[std::get<0>(jj)].isboundary == false) {
							//sum( jjf_m /p0*p0 * dWdt_jb_jjf)
							jjf += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
						// boundaryneighbours jjb
						if (PartC[std::get<0>(jj)].isboundary) {
							//sum( jjb_m /p0*p0 * dWdt_jb_jjb)
							jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
					}
					if_jb += glm::dot(PartC[std::get<0>(jf)].m * (-jjf - 2 * gammapres * jjb), std::get<2>(jf));
				}
			}
			PartC[i].Aff = (deltaT * deltaT) * (if_jf1 + if_jf2 + if_jb);
		}
	}
}
void makeAllAfffparallel(std::vector<iisphparticle>& PartC) {
	glm::vec3 jjf;
	glm::vec3 jjb;

#pragma omp parallel for private(jjf, jjb) 
	for (int i = 0; i < (var_fluidpart); ++i) {
		iisphparticle& Part = PartC[i];
		Part.Aff = 0.f;
		float if_jf1 = 0;
		float if_jf2 = 0;
		float if_jb = 0;
		auto& neighbors = Part.IdNSubKernelder;

		for (std::tuple<int, glm::vec3, glm::vec3>& jf : neighbors) {
			if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != Part.index) {
				jjf = glm::vec3(0.f, 0.f, 0.f);
				jjb = glm::vec3(0.f, 0.f, 0.f);

				for (std::tuple<int, glm::vec3, glm::vec3>& jj : neighbors) {
					if (PartC[std::get<0>(jj)].isboundary == false && PartC[std::get<0>(jj)].index != Part.index) {
						jjf += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
					}
					if (PartC[std::get<0>(jj)].isboundary) {
						jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
					}
				}
				if_jf1 += PartC[std::get<0>(jf)].m * glm::dot((-jjf - 2 * gammapres * jjb), std::get<2>(jf));
				if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((Part.m / (p0p0) * -std::get<2>(jf)), std::get<2>(jf));
			}
			if (PartC[std::get<0>(jf)].isboundary == true && PartC[std::get<0>(jf)].index != Part.index) {
				jjf = glm::vec3(0.f, 0.f, 0.f);
				jjb = glm::vec3(0.f, 0.f, 0.f);

				for (std::tuple<int, glm::vec3, glm::vec3>& jj : neighbors) {
					if (PartC[std::get<0>(jj)].isboundary == false) {
						jjf += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
					}
					if (PartC[std::get<0>(jj)].isboundary) {
						jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
					}
				}
				if_jb += glm::dot(PartC[std::get<0>(jf)].m * (-jjf - 2 * gammapres * jjb), std::get<2>(jf));
			}
		}
		Part.Aff = (deltaT * deltaT) * (if_jf1 + if_jf2 + if_jb);

	}

#pragma omp parallel for private(jjf, jjb) 
	for (int i = var_fluidpart; i < (var_fluidpart +var_spezialboundpart); ++i) {
		iisphparticle& Part = PartC[i];
		Part.Aff = 0.f;
		float if_jf1 = 0;
		float if_jf2 = 0;
		float if_jb = 0;
		auto& neighbors = Part.IdNSubKernelder;

		for (std::tuple<int, glm::vec3, glm::vec3>& jf : neighbors) {
			if (PartC[std::get<0>(jf)].isboundary == false) {
				//std::cout << Part.precompAff << "firstcase" << std::endl;
				jjf = glm::vec3(0.f, 0.f, 0.f);
				jjb = glm::vec3(0.f, 0.f, 0.f);
				for (std::tuple<int, glm::vec3, glm::vec3>& jj : neighbors) {
					if (PartC[std::get<0>(jj)].isboundary == false && PartC[std::get<0>(jj)].index != Part.index && PartC[std::get<0>(jj)].isfloatingboundary == false) {
						jjf += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
					}
					if (PartC[std::get<0>(jj)].isboundary) {
						jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
					}
				}
				if_jf1 += PartC[std::get<0>(jf)].m * glm::dot((-jjf - 2 * gammapres * jjb), std::get<2>(jf));
				if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((Part.m / (p0p0) * -std::get<2>(jf)), std::get<2>(jf));
			}
			if (PartC[std::get<0>(jf)].isboundary == true) {
				//std::cout << Part.precompAff << "secondtcase" << std::endl;
				jjf = glm::vec3(0.f, 0.f, 0.f);
				jjb = glm::vec3(0.f, 0.f, 0.f);

				for (std::tuple<int, glm::vec3, glm::vec3>& jj : neighbors) {
					if (PartC[std::get<0>(jj)].isboundary == false && PartC[std::get<0>(jj)].isfloatingboundary == false) {
						jjf += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
					}
					if (PartC[std::get<0>(jj)].isboundary) {
						jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
					}
				}
				if_jb += glm::dot(PartC[std::get<0>(jf)].m * (-jjf - 2 * gammapres * jjb), std::get<2>(jf));
			}
		}
		Part.Aff = (deltaT * deltaT) * (if_jf1 + if_jf2 +if_jb);
	}
}

void makeAllAfffparallelfast(std::vector<iisphparticle>& PartC) {
	glm::vec3 Kerndelder_jf_if;
	glm::vec3 jjf;
	glm::vec3 jjb;

#pragma omp parallel for private(Kerndelder_jf_if, jjf, jjb) 
	for (int i = 0; i < var_MaxParticles; ++i) {
		iisphparticle& Part = PartC[i];
		if (Part.isboundary == false) {
			Part.Aff = 0.f;
			float if_jf1 = 0;
			float if_jf2 = 0;
			float if_jb = 0;
			auto& neighbors = Part.IdNSubKernelder;

			for (std::tuple<int, glm::vec3, glm::vec3>& jf : neighbors) {
				if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != Part.index) {
					jjf = glm::vec3(0.f, 0.f, 0.f);
					jjb = glm::vec3(0.f, 0.f, 0.f);

					for (std::tuple<int, glm::vec3, glm::vec3>& jj : neighbors) {
						if (PartC[std::get<0>(jj)].isboundary == false && PartC[std::get<0>(jj)].index != Part.index) {
							jjf += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
						if (PartC[std::get<0>(jj)].isboundary) {
							jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
					}
					if_jf1 += PartC[std::get<0>(jf)].m * glm::dot((-jjf - 2 * gammapres * jjb), std::get<2>(jf));
				}
				if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != Part.index) {
					Kerndelder_jf_if = glm::vec3(0, 0, 0);

					for (std::tuple<int, glm::vec3, glm::vec3>& tmp : PartC[std::get<0>(jf)].IdNSubKernelder) {
						if (PartC[std::get<0>(tmp)].index == Part.index) {
							Kerndelder_jf_if = std::get<2>(tmp);
							if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((Part.m / (p0p0) * Kerndelder_jf_if), std::get<2>(jf));
						}
					}
				}
				if (PartC[std::get<0>(jf)].isboundary) {
					jjf = glm::vec3(0.f, 0.f, 0.f);
					jjb = glm::vec3(0.f, 0.f, 0.f);

					for (std::tuple<int, glm::vec3, glm::vec3>& jj : neighbors) {
						if (PartC[std::get<0>(jj)].isboundary == false) {
							jjf += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
						if (PartC[std::get<0>(jj)].isboundary) {
							jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
					}
					if_jb += glm::dot(PartC[std::get<0>(jf)].m * (-jjf - 2 * gammapres * jjb), std::get<2>(jf));
				}
			}
			Part.Aff = (deltaT * deltaT) * (if_jf1 + if_jf2 + if_jb);
		}
	}
}
void makeAllAfffastwithdens(std::vector<iisphparticle>& PartC) {
	glm::vec3 Kerndelder_jf_if = glm::vec3(0, 0, 0);
	glm::vec3  jjf = glm::vec3(0.f, 0.f, 0.f);
	glm::vec3  jjb = glm::vec3(0.f, 0.f, 0.f);
	//big sum over all fluidpart i
//#pragma omp parallel for
	for (int i = 0; i < var_MaxParticles; i++) {
		if (PartC[i].isboundary == false) {
			PartC[i].Aff = 0.f;
			float if_jf1 = 0;
			float if_jf2 = 0;
			float if_jb = 0;
			auto& neighbors = PartC[i].IdNSubKernelder;
			// sum over neighbours jf from if exept if
			for (std::tuple<int, glm::vec3, glm::vec3>& jf : neighbors) {
				if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != PartC[i].index) {
					jjf = glm::vec3(0.f, 0.f, 0.f);
					jjb = glm::vec3(0.f, 0.f, 0.f);
					// sum over all neigbours jj from if exept if
					for (std::tuple<int, glm::vec3, glm::vec3>& jj : neighbors) {
						// fluidneighbours jjf
						if (PartC[std::get<0>(jj)].isboundary == false && PartC[std::get<0>(jj)].index != PartC[i].index) {
							//sum( jjf_m /p0*p0 * dWdt_jf_jjf)
							jjf += (PartC[std::get<0>(jj)].m / (PartC[std::get<0>(jj)].density * PartC[std::get<0>(jj)].density)) * std::get<2>(jj);
						}
						// boundaryneighbours jjb
						if (PartC[std::get<0>(jj)].isboundary) {
							//sum( jjb_m /p0*p0 * dWdt_jf_jjb)
							jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
					}
					if_jf1 += PartC[std::get<0>(jf)].m * glm::dot((-jjf - 2 * gammapres * jjb), std::get<2>(jf));
				}
				if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != PartC[i].index) {
					Kerndelder_jf_if = glm::vec3(0, 0, 0);
					// sum over all neigbours jj from jf
					if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((PartC[i].m / (PartC[i].density * PartC[i].density) *-std::get<2>(jf)), std::get<2>(jf));
					//for (std::tuple<int, glm::vec3, glm::vec3>& tmp : PartC[std::get<0>(jf)].IdNSubKernelder) {
						// if jj == if
						//if (PartC[std::get<0>(tmp)].index == PartC[i].index) {
							//Kerndelder_jf_if = std::get<2>(tmp);
							// sum ( jf_m * (if_m/(p0*p0) * dWdt_jf_if)* dWdt_if_jf
							//if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((PartC[i].m / (PartC[i].density * PartC[i].density) * Kerndelder_jf_if), std::get<2>(jf));
						//}
					//}
				}
				if (PartC[std::get<0>(jf)].isboundary) {
					jjf = glm::vec3(0.f, 0.f, 0.f);
					jjb = glm::vec3(0.f, 0.f, 0.f);
					// sum over all neigbours jj from if
					for (std::tuple<int, glm::vec3, glm::vec3>& jj : neighbors) {
						// fluidneighbours jjf
						if (PartC[std::get<0>(jj)].isboundary == false) {
							//sum( jjf_m /p0*p0 * dWdt_jb_jjf)
							jjf += (PartC[std::get<0>(jj)].m / (PartC[std::get<0>(jj)].density * PartC[std::get<0>(jj)].density)) * std::get<2>(jj);
						}
						// boundaryneighbours jjb
						if (PartC[std::get<0>(jj)].isboundary) {
							//sum( jjb_m /p0*p0 * dWdt_jb_jjb)
							jjb += (PartC[std::get<0>(jj)].m / (p0p0)) * std::get<2>(jj);
						}
					}
					if_jb += glm::dot(PartC[std::get<0>(jf)].m * (-jjf - 2 * gammapres * jjb), std::get<2>(jf));
				}
			}
			PartC[i].Aff = (deltaT * deltaT) * (if_jf1 + if_jf2 + if_jb);
		}
	}
}
void computeAllDensErr(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < var_MaxParticles; i++) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isboundary == false) {
			Part.densityerror = (Part.AP - Part.sf)/p0;
		}
	}
}

float makeAlldenserrAvg(std::vector<iisphparticle>& var_PartC) {
	float densityerroraverage = 0.f;
	max_singledens = 0;
	if (ignorep0tolow) {
#pragma omp parallel for reduction(+:densityerroraverage)
		for (int i = 0; i < usemefordens.size(); ++i) {
			iisphparticle& Part = var_PartC[usemefordens[i]];
			if (absinterrupt) {
				densityerroraverage += 100 * std::abs(Part.densityerror);
				max_singledens = std::max(max_singledens, std::abs(Part.densityerror));
			}
			else {
				densityerroraverage += 100 * Part.densityerror;
				max_singledens = std::max(max_singledens, std::abs(Part.densityerror));
			}
		}
		max_singledens = 100 * max_singledens;
		return densityerroraverage / usemefordens.size();
	}
	/*
	if (ignoreincomplete) {
#pragma omp parallel for reduction(+:densityerroraverage)
		for (int i = 0; i < var_MaxParticles; ++i) {
			iisphparticle& Part = var_PartC[i];
			if (Part.isboundary == false && Part.computeme == true) {
				if (absinterrupt) {
					densityerroraverage += std::abs(Part.densityerror);
				}
				else {
					densityerroraverage += Part.densityerror;
				}
			}
		}
		return densityerroraverage / totalcomp;
	}
	
	if (ignorep0tolow) {
//#pragma omp parallel for reduction(+:densityerroraverage)
		for (int i = 0; i < var_MaxParticles; ++i) {
			iisphparticle& Part = var_PartC[i];
			//if (Part.isboundary == false && Part.denstolow == false) {
			if (Part.density >= p0 /* && Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord  && Part.isboundary == false) {
				if (absinterrupt) {
					densityerroraverage += std::abs(Part.densityerror);
				}
				else {
					densityerroraverage += Part.densityerror;
				}
			}
		}
		return (100* (densityerroraverage / numofp0high)/p0);
	}
	
	if (ignoreincomplete == false && ignorep0tolow == false) {
#pragma omp parallel for reduction(+:densityerroraverage)
		for (int i = 0; i < var_MaxParticles; ++i) {
			iisphparticle& Part = var_PartC[i];
			if (Part.isboundary == false) {
				if (absinterrupt) {
					densityerroraverage += std::abs(Part.densityerror);
				}
				else {
					densityerroraverage += Part.densityerror;
				}
			}
		}
		return densityerroraverage / var_fluidpart;
	}
	*/
	return 0.f; // default return in case all conditions are false
}
void makeAllPresA(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < (var_fluidpart + var_spezialboundpart); i++) {
		iisphparticle& Part = var_PartC[i];
		Part.presA = glm::vec3(0, 0, 0);
		glm::vec3 PresAf(0.f, 0.f, 0.f);
		glm::vec3 PresAb(0.f, 0.f, 0.f);
		for (std::tuple<int, glm::vec3, glm::vec3>& neig : Part.IdNSubKernelder) {
			if (Part.index != std::get<0>(neig)) {
				if (var_PartC[std::get<0>(neig)].isboundary) {
					PresAb -= var_PartC[std::get<0>(neig)].m * (2 * Part.pressureiter / (p0p0)) * std::get<2>(neig);
				}
				else {
					PresAf -= var_PartC[std::get<0>(neig)].m * (Part.pressureiter / (p0p0)+var_PartC[std::get<0>(neig)].pressureiter / (p0p0)) * std::get<2>(neig);
				}
			}
		}
		Part.presA = (PresAf + (gammapres * PresAb));
	}
}
void makeAllPresAwithbound(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < (var_fluidpart + var_spezialboundpart); i++) {
		iisphparticle& Part = var_PartC[i];
		Part.presA = glm::vec3(0, 0, 0);
		glm::vec3 PresAf(0.f, 0.f, 0.f);
		glm::vec3 PresAb(0.f, 0.f, 0.f);
		for (std::tuple<int, glm::vec3, glm::vec3>& neig : Part.IdNSubKernelder) {
			if (Part.index != std::get<0>(neig)) {
				if (var_PartC[std::get<0>(neig)].isboundary) {
					PresAb -= var_PartC[std::get<0>(neig)].m * (Part.pressureiter / (p0p0)+var_PartC[std::get<0>(neig)].pressureiter / (p0p0)) * std::get<2>(neig);
				}
				else {
					PresAf -= var_PartC[std::get<0>(neig)].m * (Part.pressureiter / (p0p0)+var_PartC[std::get<0>(neig)].pressureiter / (p0p0)) * std::get<2>(neig);
				}
			}
		}
		Part.presA = (PresAf + (gammapres * PresAb));
	}
}
void makeBoundPres(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = var_fluidpart + var_spezialboundpart; i < var_MaxParticles; i++) {
		iisphparticle& Part = var_PartC[i];
		float firstsum = 0;
		float secondsum = 0;
		float thirdsum = 0;
		for (auto& neigh : Part.IdNdistNsub) {
			if (var_PartC[std::get<0>(neigh)].isboundary == false) {
				firstsum += var_PartC[std::get<0>(neigh)].pressureiter * makesinglekernel(Part.pos, var_PartC[std::get<0>(neigh)].pos);
				secondsum += var_PartC[std::get<0>(neigh)].density * std::get<1>(neigh) * makesinglekernel(Part.pos, var_PartC[std::get<0>(neigh)].pos);
				thirdsum += makesinglekernel(Part.pos, var_PartC[std::get<0>(neigh)].pos);
			}
		}
		Part.pressureiter = (firstsum + gravity * secondsum) / thirdsum;
	}
}
void makeBoundPres2D(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = var_fluidpart + var_spezialboundpart; i < var_MaxParticles; i++) {
		iisphparticle& Part = var_PartC[i];
		float firstsum = 0;
		float secondsum = 0;
		float thirdsum = 0;
		for (auto& neigh : Part.IdNdistNsub) {
			if (var_PartC[std::get<0>(neigh)].isboundary == false) {
				firstsum += var_PartC[std::get<0>(neigh)].pressureiter * makesinglekernel2D(Part.pos, var_PartC[std::get<0>(neigh)].pos);
				secondsum += var_PartC[std::get<0>(neigh)].density * std::get<1>(neigh) * makesinglekernel2D(Part.pos, var_PartC[std::get<0>(neigh)].pos);
				thirdsum += makesinglekernel(Part.pos, var_PartC[std::get<0>(neigh)].pos);
			}
		}
		Part.pressureiter = (firstsum + gravity * secondsum) / thirdsum;
	}
}

float makesinglekernel(glm::vec3& posi, glm::vec3& posj) {
	float distX = posi.x - posj.x;
	float distY = posi.y - posj.y;
	float distZ = posi.z - posj.z;
	float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
	float d = distance / h;
	float t1 = std::max((1 - d), 0.f);
	float t2 = std::max((2 - d), 0.f);
	return alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)) * 1.0000280157848;
}
float makesinglekernel2D(glm::vec3& posi, glm::vec3& posj) {
	float distX = posi.x - posj.x;
	float distY = posi.y - posj.y;
	float distance = sqrt(distX * distX + distY * distY);
	float d = distance / h;
	float t1 = std::max((1 - d), 0.f);
	float t2 = std::max((2 - d), 0.f);
	return alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)) * 1.0000280157848;
}


glm::vec3 makesinglekernelder(glm::vec3& posi, glm::vec3& posj) {
	float distX = posi.x - posj.x;
	float distY = posi.y - posj.y;
	float distZ = posi.z - posj.z;
	float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
	float d = distance / h;
	float t1 = std::max((1 - d), 0.f);
	float t2 = std::max((2 - d), 0.f);
	return alphaTwoD * (posi - posj) / (d*h*h)* (-3 * t2 * t2 + 12 * t1 * t1);
}


void makeAllPresAtwoD(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < (var_fluidpart + var_spezialboundpart); i++) {
		iisphparticle& Part = var_PartC[i];
		Part.presA = glm::vec3(0, 0, 0);
		glm::vec3 PresAf(0.f, 0.f, 0.f);
		glm::vec3 PresAb(0.f, 0.f, 0.f);
		for (std::tuple<int, glm::vec3, glm::vec3>& neig : Part.IdNSubKernelder) {
			if (Part.index != std::get<0>(neig)) {
				if (var_PartC[std::get<0>(neig)].isboundary) {
					PresAb -= var_PartC[std::get<0>(neig)].m * (2 * Part.pressureiter / (p0p0)) * std::get<2>(neig);
				}
				else {
					PresAf -= var_PartC[std::get<0>(neig)].m * (Part.pressureiter / (p0p0)+var_PartC[std::get<0>(neig)].pressureiter / (p0p0)) * std::get<2>(neig);
				}
			}
		}
		Part.presA = (PresAf + (gammapres * PresAb)) * glm::vec3(1.f, 1.f, 0);
	}
}
void makeAllPresAwithboundtwoD(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < (var_fluidpart + var_spezialboundpart); i++) {
		iisphparticle& Part = var_PartC[i];
		Part.presA = glm::vec3(0, 0, 0);
		glm::vec3 PresAf(0.f, 0.f, 0.f);
		glm::vec3 PresAb(0.f, 0.f, 0.f);
		for (std::tuple<int, glm::vec3, glm::vec3>& neig : Part.IdNSubKernelder) {
			if (Part.index != std::get<0>(neig)) {
				if (var_PartC[std::get<0>(neig)].isboundary) {
					PresAb -= var_PartC[std::get<0>(neig)].m * (Part.pressureiter / (p0p0)+var_PartC[std::get<0>(neig)].pressureiter / (p0p0)) * std::get<2>(neig);
				}
				else {
					PresAf -= var_PartC[std::get<0>(neig)].m * (Part.pressureiter / (p0p0)+var_PartC[std::get<0>(neig)].pressureiter / (p0p0)) * std::get<2>(neig);
				}
			}
		}
		Part.presA = (PresAf + (gammapres * PresAb));
	}
}
void makeAllPresAwithdens(std::vector<iisphparticle>& var_PartC) {
//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isboundary == false) {
			Part.presA = glm::vec3(0, 0, 0);
			glm::vec3 PresAf(0.f, 0.f, 0.f);
			glm::vec3 PresAb(0.f, 0.f, 0.f);
			for (std::tuple<int, glm::vec3, glm::vec3>& neig : Part.IdNSubKernelder) {
				if (Part.index != std::get<0>(neig)) {
					if (var_PartC[std::get<0>(neig)].isboundary) {
						PresAb -= var_PartC[std::get<0>(neig)].m * (2 * Part.pressureiter / (Part.density * Part.density)) * std::get<2>(neig);
					}
					else {
						PresAf -= var_PartC[std::get<0>(neig)].m * (Part.pressureiter / (Part.density * Part.density) + var_PartC[std::get<0>(neig)].pressureiter / (var_PartC[std::get<0>(neig)].density * var_PartC[std::get<0>(neig)].density)) * std::get<2>(neig);
					}
				}
			}
			Part.presA = (PresAf + (gammapres * PresAb));
		}
	}
}
void makeAllAPupdatePres(std::vector<iisphparticle>& var_PartC) {
	float tmpap = 0;
//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isboundary == false) {
			Part.AP = 0;
			float APF = 0;
			float APB = 0;
			for (std::tuple<int, glm::vec3, glm::vec3>&neig : Part.IdNSubKernelder) {
				if (Part.index != std::get<0>(neig)) {
					if (var_PartC[std::get<0>(neig)].isboundary) {
						APB += var_PartC[std::get<0>(neig)].m * glm::dot(Part.presA, std::get<2>(neig));
					}
					else {
						APF += var_PartC[std::get<0>(neig)].m * glm::dot((Part.presA - var_PartC[std::get<0>(neig)].presA), std::get<2>(neig));
					}
				}
			}
			Part.AP = deltaT * deltaT * ((APB)+APF);
			tmpap = deltaT * deltaT * ((APB)+APF);
			if (Part.Aff != 0) {
				Part.pressureiter = glm::max(Part.pressureiter + (omega * (Part.sf - tmpap) / Part.Aff), 0.f);
			}
			else {
				Part.pressureiter = glm::max(Part.pressureiter, 0.f);
			}
		}
	}
}
/*
void secondloop(std::vector<iisphparticle>& var_PartC) {
	float tmpap = 0;
//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isboundary == false) {
			Part.AP = 0;
			float APF = 0;
			float APB = 0;
			for (std::tuple<int, glm::vec3, glm::vec3>& neig : Part.IdNSubKernelder) {
				if (Part.index != std::get<0>(neig)) {
					if (var_PartC[std::get<0>(neig)].isboundary) {
						APB += var_PartC[std::get<0>(neig)].m * glm::dot(Part.presA, std::get<2>(neig));
					}
					else {
						APF += var_PartC[std::get<0>(neig)].m * glm::dot((Part.presA - var_PartC[std::get<0>(neig)].presA), std::get<2>(neig));
					}
				}
			}
			Part.AP = deltaT * deltaT * ((APB)+APF);
			tmpap = deltaT * deltaT * ((APB)+APF);
			Part.densityerror = tmpap - Part.sf;
			if (Part.Aff != 0) {
				Part.pressureiter = glm::max(Part.pressureiter + (omega * (Part.sf - tmpap) / Part.Aff), 0.f);
			}
			else {
				Part.pressureiter = glm::max(Part.pressureiter, 0.f);
			}
		}
	}
}
*/
void secondloop(std::vector<iisphparticle>& var_PartC) {
	
#pragma omp parallel for
	for (int i = 0; i < (var_fluidpart+var_spezialboundpart); i++) {
		iisphparticle& Part = var_PartC[i];
		Part.AP = 0;
		float APF = 0;
		float APB = 0;
		const auto& neighbors = Part.IdNSubKernelder;
		for (const auto& neig : neighbors) {
			int neig_idx = std::get<0>(neig);
			const glm::vec3& kernel_deriv = std::get<2>(neig);
			if (Part.index != neig_idx) {
				if (var_PartC[neig_idx].isboundary) {
					APB += var_PartC[neig_idx].m * glm::dot(Part.presA - var_PartC[neig_idx].presA, kernel_deriv);
				}
				else {
					APF += var_PartC[neig_idx].m * glm::dot((Part.presA - var_PartC[neig_idx].presA), kernel_deriv);
				}
			}
		}
		float deltaT2 = deltaT * deltaT;
		float APB_APF = APB + APF;
		Part.AP = deltaT2 * APB_APF;
		float tmpap = deltaT2 * APB_APF;
		Part.densityerror = (Part.AP - Part.sf) / p0;

		//if (Part.Aff >= 0.000000001 || Part.Aff <= -0.000000001) 
		if (Part.Aff != 0)
		{
			Part.pressureiter = glm::max(Part.pressureiter + (omega * (Part.sf - tmpap) / Part.Aff), 0.f);
			//Part.pressureiter = Part.pressureiter + (omega * (Part.sf - tmpap) / Part.Aff);
		}
		else {
			Part.pressureiter = glm::max(Part.pressureiter, 0.f);
			//Part.pressureiter = Part.pressureiter ;
		}
	}
}
void makeAllVandP(std::vector<iisphparticle>& var_PartC) {
	glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.00f), baseRotationSpeed * deltaT, glm::vec3(0, 1, 0));

#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		auto& Part = var_PartC[i];
		if (Part.isboundary == false && Part.isfloatingboundary == false) {
			if (Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord) {
				Part.vel += deltaT * (Part.presA + Part.nonpresA);
				Part.r = 30;
				Part.g = 143;
				Part.b = 148;
				Part.pos += Part.vel * deltaT;
				if (colorvel) {
					Part.g = glm::min(30 + (cloloroffset * glm::length(Part.vel)), 143.f);
				}
				if (colordens) {
					Part.r = 255;
					Part.g = 200 - glm::min(10 * cloloroffset * glm::length(100 * std::abs(Part.density - p0) / p0), 199.f);
					Part.b = 200 - glm::min(10 * cloloroffset * glm::length(100 * std::abs(Part.density - p0) / p0), 199.f);
				}
				if (colorpres) {
					Part.r = 255;
					Part.g = 200 - glm::min(0.05f * cloloroffset * glm::length(Part.pressureiter / p0), 199.f);
					Part.b = 200 - glm::min(0.05f * cloloroffset * glm::length(Part.pressureiter / p0), 199.f);
				}
			}

			if (glm::length(Part.vel) > 40.f && highlightextremefast) {
				Part.g = 0;
				Part.r = 255;
			}
			if (!Part.computeme && highlightunderpop) {
				Part.g = 255;
				Part.r = 0;
				Part.b = 0;
			}
			if (highlightusedforcomp) {
				if (Part.denstolow == false) {
					Part.g = 255;
					Part.r = 0;
					Part.b = 0;
				}
			}
			if (highlightdenserr) {
				if (std::abs(std::abs(100 * (Part.density - p0) / p0) - std::abs(100*Part.densityerror)) >= denistyerrormax) {
					Part.g = 0;
					Part.r = 255;
					Part.b = 0;
				}
			}
		}
		if (Part.ismovingboundary == true) {
			glm::vec3 centeredPos = Part.pos - glm::vec3((1) * h / 2.0f, (1) * h / 2.0f, 0.0f);
			glm::vec4 rotatedPosition = rotationMatrix * glm::vec4(centeredPos, 1.0f);
			oldpos = Part.pos;
			Part.pos = glm::vec3(rotatedPosition) + glm::vec3((1) * h / 2.0f, (1) * h / 2.0f, 0.0f);
			Part.vel = (Part.pos - oldpos) / deltaT;
		}
	}
}
void makeboundmass(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	hashmap.clear();
	insertAllParticlesIntoHashmap(var_PartC, hashmap);
//#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isboundary == true && Part.isfloatingboundary ==false) {
			Part.IdNdistNsub.clear();
			std::unordered_set<int> uniqueIds; // Set to track unique IDs
			for (int x = -1; x <= 1; x++) {
				for (int y = -1; y <= 1; y++) {
					for (int z = -1; z <= 1; z++) {
						int hash = hashFunction3D(Part.pos.x + x * searchRadius, Part.pos.y + y * searchRadius, Part.pos.z + z * searchRadius);
						if (hashmap.find(hash) != hashmap.end()) {
							for (int idx : hashmap[hash].particles) {
								if (uniqueIds.find(idx) == uniqueIds.end()) { // Check if ID is already added
									float distX = var_PartC[idx].pos.x - Part.pos.x;
									float distY = var_PartC[idx].pos.y - Part.pos.y;
									float distZ = var_PartC[idx].pos.z - Part.pos.z;
									float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
									if (distance < searchRadius && var_PartC[idx].isboundary) {
										Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
										uniqueIds.insert(idx); // Mark this ID as added
									}
								}
							}
						}
					}
				}
			}
		}
	}
//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isboundary == true && Part.isfloatingboundary == false) {
			float kernsum = 0;
			Part.Kernel.clear();
//#pragma omp parallel for
			for (const auto& neig : Part.IdNdistNsub) {
				float d = std::get<1>(neig) / h;
				float t1 = std::max((1 - d), 0.f);
				float t2 = std::max((2 - d), 0.f);
				Part.Kernel.push_back(alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				kernsum += alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1));
			}
			Part.m = gammabound * p0 / std::max( kernsum, 0.0000000000000001f);
		}

	}
}


void makepartmass(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	hashmap.clear();
	insertAllParticlesIntoHashmap(var_PartC, hashmap);
//#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isboundary == false || Part.isfloatingboundary == false) {
			Part.IdNdistNsub.clear();
			std::unordered_set<int> uniqueIds; // Set to track unique IDs
			for (int x = -1; x <= 1; x++) {
				for (int y = -1; y <= 1; y++) {
					for (int z = -1; z <= 1; z++) {
						int hash = hashFunction3D(Part.pos.x + x * searchRadius, Part.pos.y + y * searchRadius, Part.pos.z + z * searchRadius);
						if (hashmap.find(hash) != hashmap.end()) {
							for (int idx : hashmap[hash].particles) {
								if (uniqueIds.find(idx) == uniqueIds.end()) { // Check if ID is already added
									float distX = var_PartC[idx].pos.x - Part.pos.x;
									float distY = var_PartC[idx].pos.y - Part.pos.y;
									float distZ = var_PartC[idx].pos.z - Part.pos.z;
									float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
									if (distance < searchRadius) {
										Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
										uniqueIds.insert(idx); // Mark this ID as added
									}
								}
							}
						}
					}
				}
			}
		}
	}
//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isboundary == false && Part.isfloatingboundary == false) {
			float kernsum = 0;
			Part.Kernel.clear();
//#pragma omp parallel for
			for (const auto& neig : Part.IdNdistNsub) {
				float d = std::get<1>(neig) / h;
				float t1 = std::max((1 - d), 0.f);
				float t2 = std::max((2 - d), 0.f);
				Part.Kernel.push_back(alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)) );
				kernsum += alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)) ;
			}
			Part.m = p0 / std::max(gammapart * kernsum, 0.0000000000000001f);
		}

	}
}
void makeboundmassTwoD(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	hashmap.clear();
	insertAllParticlesIntoHashmap(var_PartC, hashmap);
//#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isboundary == true && Part.isfloatingboundary == false) {
			Part.IdNdistNsub.clear();
			std::unordered_set<int> uniqueIds; // Set to track unique IDs
			for (int x = -1; x <= 1; x++) {
				for (int y = -1; y <= 1; y++) {
					for (int z = -1; z <= 1; z++) {
						int hash = hashFunction3D(Part.pos.x + x * searchRadius, Part.pos.y + y * searchRadius, Part.pos.z + z * searchRadius);
						if (hashmap.find(hash) != hashmap.end()) {
							for (int idx : hashmap[hash].particles) {
								if (uniqueIds.find(idx) == uniqueIds.end()) { // Check if ID is already added
									float distX = var_PartC[idx].pos.x - Part.pos.x;
									float distY = var_PartC[idx].pos.y - Part.pos.y;
									float distZ = var_PartC[idx].pos.z - Part.pos.z;
									float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
									if (distance < searchRadius && var_PartC[idx].isboundary) {
										Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
										uniqueIds.insert(idx); // Mark this ID as added
									}
								}
							}
						}
					}
				}
			}
		}
	}
//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isboundary == true && Part.isfloatingboundary == false) {
			float kernsum = 0;
			Part.Kernel.clear();
//#pragma omp parallel for
			for (const auto& neig : Part.IdNdistNsub) {
				float d = std::get<1>(neig) / h;
				float t1 = std::max((1 - d), 0.f);
				float t2 = std::max((2 - d), 0.f);
				Part.Kernel.push_back(alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)) );
				kernsum += alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)) ;
			}
			Part.m = gammabound * p0 / std::max( kernsum, 0.0000000000000001f);
		}

	}
}


void makepartmassTwoD(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	hashmap.clear();
	insertAllParticlesIntoHashmap(var_PartC, hashmap);
//#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isboundary == false) {
			Part.IdNdistNsub.clear();
			std::unordered_set<int> uniqueIds; // Set to track unique IDs
			for (int x = -1; x <= 1; x++) {
				for (int y = -1; y <= 1; y++) {
					for (int z = -1; z <= 1; z++) {
						int hash = hashFunction3D(Part.pos.x + x * searchRadius, Part.pos.y + y * searchRadius, Part.pos.z + z * searchRadius);
						if (hashmap.find(hash) != hashmap.end()) {
							for (int idx : hashmap[hash].particles) {
								if (uniqueIds.find(idx) == uniqueIds.end()) { // Check if ID is already added
									float distX = var_PartC[idx].pos.x - Part.pos.x;
									float distY = var_PartC[idx].pos.y - Part.pos.y;
									float distZ = var_PartC[idx].pos.z - Part.pos.z;
									float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
									if (distance < searchRadius) {
										Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
										uniqueIds.insert(idx); // Mark this ID as added
									}
								}
							}
						}
					}
				}
			}
		}
	}
//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isboundary == false && Part.isfloatingboundary == false) {
			float kernsum = 0;
			Part.Kernel.clear();
//#pragma omp parallel for
			for (const auto& neig : Part.IdNdistNsub) {
				float d = std::get<1>(neig) / h;
				float t1 = std::max((1 - d), 0.f);
				float t2 = std::max((2 - d), 0.f);
				Part.Kernel.push_back(alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)) );
				kernsum += alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)) ;
			}
			Part.m = p0 / std::max(gammapart * kernsum, 0.0000000000000001f);
		}

	}
}

void makefloatingpartmass(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	hashmap.clear();
	insertAllParticlesIntoHashmap(var_PartC, hashmap);
	//#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isfloatingboundary == true) {
			Part.IdNdistNsub.clear();
			std::unordered_set<int> uniqueIds; // Set to track unique IDs
			for (int x = -1; x <= 1; x++) {
				for (int y = -1; y <= 1; y++) {
					for (int z = -1; z <= 1; z++) {
						int hash = hashFunction3D(Part.pos.x + x * searchRadius, Part.pos.y + y * searchRadius, Part.pos.z + z * searchRadius);
						if (hashmap.find(hash) != hashmap.end()) {
							for (int idx : hashmap[hash].particles) {
								if (uniqueIds.find(idx) == uniqueIds.end()) { // Check if ID is already added
									float distX = var_PartC[idx].pos.x - Part.pos.x;
									float distY = var_PartC[idx].pos.y - Part.pos.y;
									float distZ = var_PartC[idx].pos.z - Part.pos.z;
									float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
									if (distance < searchRadius && var_PartC[idx].isfloatingboundary) {
										Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
										uniqueIds.insert(idx); // Mark this ID as added
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isfloatingboundary == true) {
			float kernsum = 0;
			Part.Kernel.clear();
			//#pragma omp parallel for
			for (const auto& neig : Part.IdNdistNsub) {
				float d = std::get<1>(neig) / h;
				float t1 = std::max((1 - d), 0.f);
				float t2 = std::max((2 - d), 0.f);
				Part.Kernel.push_back(alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				kernsum += alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1));
			}
			Part.m = 0.7 * p0 / std::max(kernsum, 0.0000000000000001f);
		}

	}
}
void makefloatingpartmass2d(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	hashmap.clear();
	insertAllParticlesIntoHashmap(var_PartC, hashmap);
	//#pragma omp parallel for
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isfloatingboundary == true) {
			Part.IdNdistNsub.clear();
			std::unordered_set<int> uniqueIds; // Set to track unique IDs
			for (int x = -1; x <= 1; x++) {
				for (int y = -1; y <= 1; y++) {
					for (int z = -1; z <= 1; z++) {
						int hash = hashFunction3D(Part.pos.x + x * searchRadius, Part.pos.y + y * searchRadius, Part.pos.z + z * searchRadius);
						if (hashmap.find(hash) != hashmap.end()) {
							for (int idx : hashmap[hash].particles) {
								if (uniqueIds.find(idx) == uniqueIds.end()) { // Check if ID is already added
									float distX = var_PartC[idx].pos.x - Part.pos.x;
									float distY = var_PartC[idx].pos.y - Part.pos.y;
									float distZ = var_PartC[idx].pos.z - Part.pos.z;
									float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
									if (distance < searchRadius && var_PartC[idx].isfloatingboundary) {
										Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
										uniqueIds.insert(idx); // Mark this ID as added
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isfloatingboundary == true) {
			float kernsum = 0;
			Part.Kernel.clear();
			//#pragma omp parallel for
			for (const auto& neig : Part.IdNdistNsub) {
				float d = std::get<1>(neig) / h;
				float t1 = std::max((1 - d), 0.f);
				float t2 = std::max((2 - d), 0.f);
				Part.Kernel.push_back(alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				kernsum += alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1));
			}
			Part.m = 0.7 * p0 / std::max(kernsum, 0.0000000000000001f);
		}

	}
}


void initrigidbodies(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap) {
	posofcenterofmass = glm::vec3(0.f, 0.f, 0.f);
	velofcenterofmass = glm::vec3(0.f, 0.f, 0.f);
	rotMat = glm::mat3(1.0f); // Identity matrix
	xCM = glm::vec3(0.f, 0.f, 0.f);
	vCM = glm::vec3(0.f, 0.f, 0.f);
	//A = glm::mat3(1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f); // Identity matrix
	A = glm::mat3(1.f);
	//A[0][0] = 1.f;
	//A[1][1] = 1.f;
	//A[2][2] = 1.f;
	//A = glm::mat3(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f);
	L = glm::vec3(0.f, 0.f, 0.f);
	I_inv = glm::mat3(1.0f); // Identity matrix
	inertiaTensor = glm::mat3(0.0f);
	omegarigidbody = glm::vec3(0.f, 0.f, 0.f);
	allrigidmass = 0.0f;

	// Compute center of mass
	for (int i = 0; i < PartC.size(); ++i) {
		iisphparticle& Part = PartC[i];
		if (Part.isfloatingboundary) {
			Part.m = Part.m*0.9;
			allrigidmass += Part.m;
			xCM += Part.m * Part.pos;
			//std::cout << "Particle " << i << " mass: " << Part.m << ", position: (" << Part.pos.x << ", " << Part.pos.y << ", " << Part.pos.z << ")" << std::endl;
		}
	}

	if (allrigidmass > 0) {
		xCM /= allrigidmass;
	}
	//std::cout << "Center of mass: (" << xCM.x << ", " << xCM.y << ", " << xCM.z << "), Total mass: " << allrigidmass << std::endl;

	// Compute inertia tensor
	for (int i = 0; i < PartC.size(); ++i) {
		iisphparticle& Part = PartC[i];
		if (Part.isfloatingboundary) {
			glm::vec3 r = Part.pos - xCM;
			Part.relpos = r;
			float mass = Part.m;
			glm::mat3 rr = glm::outerProduct(r, r);
			inertiaTensor -= mass * rr;
		}
	}
	//std::cout << "Inertia Tensor: \n"
	//	<< inertiaTensor[0][0] << " " << inertiaTensor[0][1] << " " << inertiaTensor[0][2] << "\n"
	//	<< inertiaTensor[1][0] << " " << inertiaTensor[1][1] << " " << inertiaTensor[1][2] << "\n"
	//	<< inertiaTensor[2][0] << " " << inertiaTensor[2][1] << " " << inertiaTensor[2][2] << "\n";

	if (glm::determinant(inertiaTensor) != 0) {
		inertiaTensorInverse = glm::inverse(inertiaTensor);
		I_inv = inertiaTensorInverse;
		//std::cout << "Inertia Tensor Inverse: \n"
		//	<< I_inv[0][0] << " " << I_inv[0][1] << " " << I_inv[0][2] << "\n"
		//	<< I_inv[1][0] << " " << I_inv[1][1] << " " << I_inv[1][2] << "\n"
		//	<< I_inv[2][0] << " " << I_inv[2][1] << " " << I_inv[2][2] << "\n";
	}

}

void updaterigidbody(std::vector<iisphparticle>& PartC) {
	glm::vec3 torguefac = torque;
	torque= glm::vec3(0.f, 0.f, 0.f);
	glm::vec3 linforce(0.f, 0.f, 0.f);
	// Parallelisierte Schleife mit OpenMP
//#pragma omp parallel for private(torque, linforce)
	for (int i = var_fluidpart; i < PartC.size(); ++i) {
		iisphparticle& Part = PartC[i];
		if (Part.isfloatingboundary) {
			glm::vec3 force = ((Part.presA) + Part.nonpresA) * Part.m; // Annahme, dass Kräfte Beschleunigungen sind
			glm::vec3 r = Part.pos - xCM; 
			torque += glm::cross(r, force);;
			linforce += force;
		}
	}

	
	//torque = torque*0.01f;
	
	vCM += deltaT * linforce / allrigidmass;
	xCM += deltaT * vCM;
	
	// Orientierungsmatrix aktualisieren
	glm::mat3 skewOmega = skewSymmetricMatrix(omegarigidbody);
	A += deltaT * skewOmega * A;
	A = A;
	// Sicherstellen, dass die Orientierungsmatrix orthonormal bleibt using Gram-Schmidt process
	glm::vec3 col1 = A[0];
	glm::vec3 col2 = A[1];
	glm::vec3 col3 = A[2];

	col1 = glm::normalize(A[0]);

	col2 = A[1] - glm::dot(col1, A[1]) * col1;
	col2 = glm::normalize(col2);

	col3 = A[2] - glm::dot(col1, A[2]) * col1 - glm::dot(col2, A[2]) * col2;
	col3 = glm::normalize(col3);
	
	A = glm::mat3(col1, col2, col3);
	L += deltaT * torque;

	if (glm::determinant(A) != 0) {
		glm::mat3 A_T = glm::transpose(A);
		glm::mat3 I_inv = A * inertiaTensorInverse * A_T;
		omegarigidbody = I_inv * L;
	}
	else {
		// Handle the case where A is not invertible
		omegarigidbody = glm::vec3(0.f); // or some other appropriate action
	}
	// Partikelpositionen und -geschwindigkeiten aktualisieren
//#pragma omp parallel for
	for (int i = var_fluidpart; i < PartC.size(); ++i) {
		iisphparticle& Part = PartC[i];
		if (Part.isfloatingboundary) {
			glm::vec3 relpostmp = A * (Part.relpos);
			Part.pos = xCM + relpostmp; // Position mit linearer Geschwindigkeit aktualisieren
			Part.vel = vCM + glm::cross(omegarigidbody, relpostmp);
		}
	}
}

// Helper function to create a skew-symmetric matrix from a vector
glm::mat3 skewSymmetricMatrix(const glm::vec3& v) {
	return glm::mat3(
		0, -v.z, v.y,
		v.z, 0, -v.x,
		-v.y, v.x, 0
	);
}

glm::mat3 skewSymmetricMatrix2d(const glm::vec3& v) {
	return glm::mat3(
		0, -v.z, v.y,
		v.z, 0, -v.x,
		-v.y, v.x, 0
	);
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