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


void var_iisph(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap) {
	if (makesingle) {
		auto start_o = std::chrono::high_resolution_clock::now();
		makecfltrue(PartC);

		auto start_n = std::chrono::high_resolution_clock::now();
		//hashmap.clear();
		//insertAllParticlesIntoHashmap(PartC, hashmap);
		//findneighbours
		//findAllNeighbours(PartC, hashmap);
		clearuniformgrid();
		insertparticlesinuniformgrid(PartC);
		findAllNeighbourscompact3D(PartC);
		auto end_n = std::chrono::high_resolution_clock::now();
		//makeKernel and Kernelder
		makeAllKernelAndKernelDer(PartC);
		auto end_kern = std::chrono::high_resolution_clock::now();
		//compute density
		computeAllDens(PartC);
		auto end_dens = std::chrono::high_resolution_clock::now();
		//nonpresAcceleartion
		MakeAllNonpresA(PartC);
		auto end_nonpresa = std::chrono::high_resolution_clock::now();
		PredictAllVel(PartC);
		auto end_predvel = std::chrono::high_resolution_clock::now();
		//compute sf
		computeAllSF(PartC);
		auto end_sf = std::chrono::high_resolution_clock::now();
		//compute diagonal element
		makeAllAfffparallel(PartC);
		auto end_aff = std::chrono::high_resolution_clock::now();
		// set p0f
#pragma omp parallel for
		for (int i = 0; i < var_MaxParticles; i++) {
			if (PartC[i].isboundary == false) {
				//PartC[i].pressureiter = 0;
				if (PartC[i].Aff != 0) {
					PartC[i].pressureiter = glm::max((omega * PartC[i].sf / PartC[i].Aff), 0.f);
				}
				if (PartC[i].Aff == 0) {
					PartC[i].pressureiter = 0;
				}
			}
		}
		computeAllDensErr(PartC);
		int l = 1;
		while (makeAlldenserrAvg(PartC) > denistyerrormax || l < 2/* || 100 * max_singledens / p0 > 50*/) {
			//first for
			makeAllPresA(PartC);
			//second for
			secondloop(PartC);
			l++;
			if (l > currentitermax) {
				//paussimul = true;
				break;
			}
		}
		auto end_loop = std::chrono::high_resolution_clock::now();
		currentiter = l;

		//update velocity and Position
		makeAllVandP(PartC);
		updaterigidbody(PartC);
		auto end_o = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_o - end_loop);
		othercomputationTime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_n - start_n);
		neighSearchTime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_kern - end_n);
		kerneltime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_dens - end_kern);
		densitytime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_nonpresa - end_dens);
		nonpresatime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_predvel - end_nonpresa);
		predveltime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_sf - end_predvel);
		computeallsftime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_aff - end_sf);
		computeallafftime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_loop - end_aff);
		looptime = duration.count();
		densitysnew.push_back(makeAlldenserrAvg(PartC));
		densitys.push_back(denserrold);
	}
	else
	{
		auto start_o = std::chrono::high_resolution_clock::now();
		makecfltrue(PartC);

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
		computeAllDens(PartC);
		auto end_dens = std::chrono::high_resolution_clock::now();
		//nonpresAcceleartion
		MakeAllNonpresA(PartC);
		auto end_nonpresa = std::chrono::high_resolution_clock::now();
		PredictAllVel(PartC);
		auto end_predvel = std::chrono::high_resolution_clock::now();
		//compute sf
		computeAllSF(PartC);
		auto end_sf = std::chrono::high_resolution_clock::now();
		//compute diagonal element
		makeAllAfffparallel(PartC);
		auto end_aff = std::chrono::high_resolution_clock::now();
		// set p0f
#pragma omp parallel for
		for (int i = 0; i < var_MaxParticles; i++) {
			
			if (PartC[i].isboundary == false) {
				//PartC[i].pressureiter = 0;
				if (PartC[i].Aff != 0) {
					PartC[i].pressureiter = glm::max((omega * PartC[i].sf / PartC[i].Aff), 0.f);
				}
				if (PartC[i].Aff == 0) {
					PartC[i].pressureiter = 0;
				}
			}
		}
		computeAllDensErr(PartC);
		int l = 1;
		while (makeAlldenserrAvg(PartC) > denistyerrormax || l < 2/* || 100 * max_singledens / 
			
			
			> 50*/) {
			makeBoundPres(PartC);
			//first for
			makeAllPresAwithbound(PartC);
			//second for
			secondloop(PartC);
			l++;
			if (l > currentitermax) {
				//paussimul = true;
				break;
			}
		}
		auto end_loop = std::chrono::high_resolution_clock::now();
		currentiter = l;

		//update velocity and Position
		makeAllVandP(PartC);
		updaterigidbody(PartC);
		auto end_o = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_o - end_loop);
		othercomputationTime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_n - start_n);
		neighSearchTime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_kern - end_n);
		kerneltime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_dens - end_kern);
		densitytime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_nonpresa - end_dens);
		nonpresatime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_predvel - end_nonpresa);
		predveltime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_sf - end_predvel);
		computeallsftime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_aff - end_sf);
		computeallafftime = duration.count();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(end_loop - end_aff);
		looptime = duration.count();
		densitysnew.push_back(makeAlldenserrAvg(PartC));
		densitys.push_back(denserrold);
	}
	
}


void twoDiisph(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap) {
	auto start_o = std::chrono::high_resolution_clock::now();
	makecfltrue(PartC);

	auto start_n = std::chrono::high_resolution_clock::now();
	hashmap.clear();
	insertAllParticlesIntoHashmap2D(PartC, hashmap);
	//findneighbours
	findAllNeighbours2D(PartC, hashmap);
	//clearuniformgrid();
	//insertparticlesinuniformgrid(PartC);
	//findAllNeighbourscompact2D(PartC);
	auto end_n = std::chrono::high_resolution_clock::now();
	//makeKernel and Kernelder
	makeAllKernelAndKernelDerTwoD(PartC);
	auto end_kern = std::chrono::high_resolution_clock::now();
	//compute density
	computeAllDenstwoD(PartC);
	auto end_dens = std::chrono::high_resolution_clock::now();
	//nonpresAcceleartion
	MakeAllNonpresAtwoD(PartC);
	auto end_nonpresa = std::chrono::high_resolution_clock::now();
	PredictAllVel(PartC);
	auto end_predvel = std::chrono::high_resolution_clock::now();
	//compute sf
	computeAllSF(PartC);
	auto end_sf = std::chrono::high_resolution_clock::now();
	//compute diagonal element
	makeAllAfffparallel(PartC);
	auto end_aff = std::chrono::high_resolution_clock::now();
	// set p0f
#pragma omp parallel for
	for (int i = 0; i < var_MaxParticles; i++) {
		if (PartC[i].isboundary == false) {
			//PartC[i].pressureiter = 0;
			if (PartC[i].Aff != 0) {
				PartC[i].pressureiter = glm::max((omega * PartC[i].sf / PartC[i].Aff), 0.f);
			}
			if (PartC[i].Aff == 0) {
				PartC[i].pressureiter = 0;
			}
		}
	}
	computeAllDensErr(PartC);
	int l = 1;
	while (makeAlldenserrAvg(PartC) > denistyerrormax || l < 2 ||max_singledens > 100) {
		if (makesingle == false) {
			makeBoundPres2D(PartC);
			makeAllPresAwithboundtwoD(PartC);
		}
		else {
			makeAllPresAtwoD(PartC);
		}
		//second for
		secondloop(PartC);
		l++;
		if (l > currentitermax) {
			//paussimul = true;
			break;
		}
	}
	auto end_loop = std::chrono::high_resolution_clock::now();
	currentiter = l;
		//update velocity and Position
	makeAllVandP(PartC);
	updaterigidbody2d(PartC);
	auto end_o = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_o - end_loop);
	othercomputationTime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_n - start_n);
	neighSearchTime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_kern - end_n);
	kerneltime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_dens - end_kern);
	densitytime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_nonpresa - end_dens);
	nonpresatime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_predvel - end_nonpresa);
	predveltime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_sf - end_predvel);
	computeallsftime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_aff - end_sf);
	computeallafftime = duration.count();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end_loop - end_aff);
	looptime = duration.count();
	densitysnew.push_back(makeAlldenserrAvg(PartC));
	densitys.push_back(denserrold);
}








