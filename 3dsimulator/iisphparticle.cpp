#include "iisphparticle.h"
#include<omp.h>
#include <iostream>
#include <string>
void iisphparticle::resetvalues() {
	vel = glm::vec3(0,0,0), acc = glm::vec3(0, 0, 0), nonpresA = glm::vec3(0, 0, 0), predictedVel = glm::vec3(0, 0, 0), presA = glm::vec3(0, 0, 0);
	r = 44, g = 2, b = 25; // Color
	size = h, density = p0, pressure = 0, predictedDens = 0, densityerror = 0, pressureiter = 0, sf = 0, Aff = 0, AP = 0;
	m = (p0 * h * h);
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
		makeAllAfffast(PartC);
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
		auto end_o = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_o - start_o);
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
		while (makeAlldenserrAvg(PartC) > denistyerrormax || l < 2/* || 100 * max_singledens / p0 > 50*/) {
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
		auto end_o = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_o - start_o);
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
	auto end_n = std::chrono::high_resolution_clock::now();
	//makeKernel and Kernelder
	makeAllKernelAndKernelDerTwoD(PartC);
	auto end_kern = std::chrono::high_resolution_clock::now();
	//compute density
	computeAllDens(PartC);
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
			makeBoundPres(PartC);
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
	auto end_o = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_o - start_o);
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





void var_init(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	ParticlesContainer.resize(0);
	var_fluidpart = var_nx * var_ny * var_nz;
	if (singlewall) {
		gammafloat = 0.7;
		var_MaxParticles = var_nx * var_ny * var_nz + (var_oneway * var_oneway * 6);
	}
	else {
		gammafloat = 1;
		var_MaxParticles = var_nx * var_ny * var_nz + (var_oneway * var_oneway * 12);
	}
	
	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	for (int x = 0; x < (var_nx ); x++) {
		for (int y = 0; y < (var_ny ); y++) {
			for (int z = 0; z < (var_nz ); z++) {
				////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 2) * (h*1), (y + 1) * (h*1), (z + 1) * (h*1));
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	//border
	for (int ii = 0; ii < var_oneway; ii++) {
		for (int j = 0; j < var_oneway; j++) {
			//1upper
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, var_oneway * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//3right
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((var_oneway-1) * h, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, 0, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//5front
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, (var_oneway - 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(0, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//7back2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (j - 1) * h, (0 - 1));
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//8right2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((var_oneway - 2) * h, (ii - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//9lower2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (0 - 1), (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//10front2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (j - 1) * h, (var_oneway - 2) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//11 left2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((0 - 1), (ii - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//12upper2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (var_oneway - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}

	/*
	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].density = p0;
		ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
		ParticlesContainer[i].nonpresA = glm::vec3(0, 0, 0);
		ParticlesContainer[i].presA = glm::vec3(0, 0, 0);
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
	*/
}

void var_teslavalveclosed(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	ParticlesContainer.resize(0);
	gammafloat = 1.f;
	var_fluidpart = 8 * 9 * 195;
	// fluid + startfloor + (start walls) + big floor + big walls + obstacles
	var_MaxParticles = 8 * 9 * 195 + (11 * 11 * 2 + 7 * 11 * 201) + (23 * 150 * 4 + 4 * 150 * 6 + 12 * 6 + 23 * 6) + (12 * 16 * 8 + 12 * 8 * 8) + (11 * 11 * 2 + 7 * 11 * 201 - 11 * 6);
	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	for (int x = -1; x < (8 - 1); x++) {
		for (int y = -1; y < (195 - 1); y++) {
			for (int z = -1; z < (9 - 1); z++) {
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 2) * h, (y + 2) * h, (z + 3) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	// back long obstacle
	float goback = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		goback -= 0.5;

	}
	// front long obstacle
	float gofront = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		gofront += 0.5;

	}
	// back short obstacle
	goback = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		goback -= 0.5;

	}
	// front long obstacle
	gofront = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		gofront += 0.5;

	}
	//start container floor
	float golower = -1;
	for (int ii = 0; ii < 11; ii++) {
		for (int j = 0; j < 11; j++) {
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 159) * h, golower * h, (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 159) * h, (golower * h - 1), (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
		golower += 0.1;
	}
	//startcontainer walls right
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 201; j++) {
			//2back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((160) * h, (j + 5) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, ((11)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((9 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, 1 * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((10 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, ((12)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}


	//big container floor and ceiling
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 23; j++) {
			//lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-2) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+4) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//lower 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-3) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+5) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	// big front and back container walls
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 6; j++) {
			//back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-4) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (17) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (18) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	/*
	// big left
	for (int ii = 0; ii < 23; ii++) {
		for (int j = 0; j < 6; j++) {
			//left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((ii)-5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	*/
	//big right
	for (int ii = 0; ii < 12; ii++) {
		for (int j = 0; j < 6; j++) {
			if (ii > 5) {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			else {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}

		}
	}
	//big left
	for (int ii = 0; ii < 12; ii++) {
		for (int j = 0; j < 6; j++) {
			if (ii > 5) {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			else {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}

		}
	}
	//start container floor
	golower = 0;
	for (int ii = 0; ii < 11; ii++) {
		for (int j = 0; j < 11; j++) {
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, golower * h, (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (golower * h - 1), (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
		golower -= 0.1;
	}
	//startcontainer walls left
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 201; j++) {
			//2back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(0, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, ((11)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(9 * h, (j + 5) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, 1 * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(-1 * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, ((12)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}


	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		if (ParticlesContainer[i].isboundary) {
			//ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
}

void var_teslavalve(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	ParticlesContainer.resize(0);
	gammafloat = 1.f;
	var_fluidpart = 8 * 9 * 195;
	// fluid + startfloor + (start walls) + big floor + big walls + obstacles
	var_MaxParticles = 8 * 9 *195 + (11 * 11 * 2 + 7 * 11 * 201) + (23 * 150 * 4 + 4 * 150 * 6 + 12 * 6 + 23 * 6) + (12 * 16 * 8 + 12 * 8 * 8) + (11 * 11 * 2 + 7 * 11 * 201 - 11*6) ;
	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	for (int x = -1; x < (8 - 1); x++) {
		for (int y = -1; y < (195 - 1); y++) {
			for (int z = -1; z < (9 - 1); z++) {
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 162) * h, (y + 2) * h, (z + 3) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	// back long obstacle
	float goback = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		goback -= 0.5;

	}
	// front long obstacle
	float gofront = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		gofront += 0.5;

	}
	// back short obstacle
	goback = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		goback -= 0.5;

	}
	// front long obstacle
	gofront = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		gofront += 0.5;

	}
	//start container floor
	float golower = -1;
	for (int ii = 0; ii < 11; ii++) {
		for (int j = 0; j < 11; j++) {
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii+159) * h, golower * h, (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii+159) * h, (golower * h - 1), (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
		golower += 0.1;
	}
	//startcontainer walls right
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 201; j++) {
			//2back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii+160) * h, (j - 1) * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((160)*h, (j +5) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii+160) * h, (j - 1) * h, ((11)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((9+160) * h, (j -1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii+160) * h, (j - 1) * h, 1 * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((10 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii+160) * h, (j - 1) * h, ((12)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}


	//big container floor and ceiling
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 23; j++) {
			//lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-2) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+4) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//lower 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-3) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+5) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	// big front and back container walls
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 6; j++) {
			//back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-4) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (17) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (18) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	/*
	// big left 
	for (int ii = 0; ii < 23; ii++) {
		for (int j = 0; j < 6; j++) {
			//left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((ii)-5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	*/
	//big right
	for (int ii = 0; ii < 12; ii++) {
		for (int j = 0; j < 6; j++) {
			if (ii > 5) {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			else {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}

		}
	}
	//big left
	for (int ii = 0; ii < 12; ii++) {
		for (int j = 0; j < 6; j++) {
			if (ii > 5) {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			else {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}

		}
	}
	//start container floor
	golower = 0;
	for (int ii = 0; ii < 11; ii++) {
		for (int j = 0; j < 11; j++) {
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, golower * h, (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (golower * h - 1), (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
		golower -= 0.1;
	}
	//startcontainer walls left
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 201; j++) {
			//2back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(0, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, ((11)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(9 * h, (j + 5) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, 1 * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(-1 * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, ((12)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}


	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
}
void var_real_teslavalve(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	ParticlesContainer.resize(0);
	gammafloat = 1.f;
	var_fluidpart = 8 * 9 * 65 + 5*7*140;
	// fluid + startfloor + (start walls) + big floor + big walls + obstacles
	var_MaxParticles = 27106+ 8 * 9 * 65 + 5 * 7 * 140;//8 * 9 * 25 + (11 * 11 * 2 + 7 * 11 * 101) + (23 * 150 * 4 + 4 * 150 * 6 + 12 * 6 + 23 * 6) + (12 * 16 * 8 + 12 * 8 * 8)+10000;
	
	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	//standing water
	for (int x = -1; x < (8 - 1); x++) {
		for (int y = -1; y < (65 - 1); y++) {
			for (int z = -1; z < (9 - 1); z++) {
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 162) * h, (y + 2) * h, (z + 3) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	//water on the ground
	for (int x = 5; x < (145 ); x++) {
		for (int y = -3; y < 2; y++) {
			for (int z = 0; z < 7; z++) {
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 10) * h, (y + 2) * h, (z + 3) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	// back long obstacle
	float goback = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		goback -= 0.5;

	}
	// front long obstacle
	float gofront = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		gofront += 0.5;

	}
	// back short obstacle
	goback = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		goback -= 0.5;

	}
	// front long obstacle
	gofront = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		gofront += 0.5;

	}
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 6; j++) {
			if ((ii < 56) || (ii > 85 && ii < 116)) {
				//front
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (12) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;

				//front2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (11) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			if (ii < 26 || (ii > 56 && ii < 86) || ii > 115) {
				//back
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//back2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (0) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;

			}
			//outer
			if ((ii > 67 && ii < 98) || (ii > 126)) {
				//front
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (17) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;

				//front2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (18) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			if (((ii < 68 && ii >36) || (ii > 98 && ii < 130 )) ) {
				//back
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//back2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-6) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;

			}
		}
	}

	//start container floor
	float golower = -1;
	for (int ii = 0; ii < 11; ii++) {
		for (int j = 0; j < 11; j++) {
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 159) * h, golower * h, (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 159) * h, (golower * h - 1), (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
		golower += 0.1;
	}
	//startcontainer walls right
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 101; j++) {
			//2back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((160) * h, (j + 5) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, ((11)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((9 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, 1 * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((10 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, ((12)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}


	//big container floor and ceiling
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 23; j++) {
			//lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-2) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+4) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//lower 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-3) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+5) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	// 
	
	// big left 
	for (int ii = 0; ii < 23; ii++) {
		for (int j = 0; j < 6; j++) {
			//left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((ii)-5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	//big left
	for (int ii = 0; ii < 12; ii++) {
		for (int j = 0; j < 6; j++) {
			if (ii > 5) {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			else {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}

		}
	}
	std::cout << i << std::endl;
	
	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
}


void watercolumn(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	// Clear the container
	ParticlesContainer.resize(0);

	// Define boundary and fluid parameters
	int var_boundarypart = 29 * 29 + 25 * 2 * (watercolheight + 100) + 27 * 2 * (watercolheight + 100);
	gammafloat = 0.7;
	if (!singlewall) {
		//gammafloat = 1;
		var_boundarypart = 29 * 29 * 2 + 25 * 2 * (watercolheight + 100) + 27 * 2 * (watercolheight + 100) + 27 * 2 * (watercolheight + 100) + 29 * 2 * (watercolheight + 100);
	}

	// Parameters for the fluid particles
	var_nx = 25;
	var_nz = 25;
	var_fluidpart = var_nx * watercolheight * var_nz;
	var_MaxParticles = var_fluidpart + var_boundarypart;
	hashsize = var_MaxParticles;

	// Resize the container to fit all particles
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	std::srand(std::time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-jitterfac, jitterfac);
	// Create fluid particles
	for (int x = 0; x < var_nx; x++) {
		for (int y = 0; y < watercolheight; y++) {
			for (int z = 0; z < var_nz; z++) {
				double rando = dis(gen);
				double randomnum = jitterfac * static_cast<double>(std::rand()) / (RAND_MAX);
				ParticlesContainer[i].pos = glm::vec3((x + rando) * h, (y + rando) * h, (z + rando) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = false;
				ParticlesContainer[i].index = i;
				i++;
			}
		}
	}
	std::cout << "fluid: " << i << std::endl;
	float boden = 0;
	// Create bottom boundary particles
	for (int x = -2; x <= 26; x++) {
		for (int y = -2; y <= -1; y++) {
			for (int z = -2; z <= 26; z++) {
				if (y > -2 || !singlewall) {
					ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (z)*h);
					ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
					ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
					ParticlesContainer[i].a = boundarya;
					ParticlesContainer[i].isboundary = true;
					ParticlesContainer[i].index = i;
					i++;
					boden++;
				}
			}
		}
	}
	std::cout << "boden: " << boden << std::endl;
	// Create front and back single
	float frontandback = 0;
	for (int x = 0; x <= 24; x++) {
		for (int y = 0; y < watercolheight + 100; y++) {
			ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (-1) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			frontandback++;
			i++;
			ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (25) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			frontandback++;
			i++;
		}
	}
	std::cout << " vorder und rück wände single : " << frontandback << std::endl;
	// Create front and back double
	frontandback = 0;
	if (!singlewall) {
		for (int x = -1; x < 26; x++) {
			for (int y = 0; y < watercolheight + 100; y++) {
				ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (-2) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				frontandback++;
				ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (26) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				frontandback++;
			}
		}
		std::cout << " vorder und rück wände double : " << frontandback << std::endl;
	}
	int leftandright = 0;
	// Create left and right single
	for (int z = -1; z <= 25; z++) {
		for (int y = 0; y < watercolheight + 100; y++) {
			ParticlesContainer[i].pos = glm::vec3((-1) * h, y * h, (z)*h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
			ParticlesContainer[i].pos = glm::vec3((25) * h, y * h, (z)*h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
		}
	}
	std::cout << " left und right wände single : " << leftandright << std::endl;
	// Create left and right double
	leftandright = 0;
	if (!singlewall) {
		for (int z = -2; z <= 26; z++) {
			for (int y = 0; y < watercolheight + 100; y++) {
				ParticlesContainer[i].pos = glm::vec3((-2) * h, y * h, (z)*h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				leftandright++;
				ParticlesContainer[i].pos = glm::vec3((26) * h, y * h, (z)*h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				leftandright++;
			}
		}
		std::cout << " left und right wände double : " << leftandright << std::endl;
	}

	// Finalize particle properties
	for (int j = 0; j < var_MaxParticles; j++) {
		ParticlesContainer[j].ismovingboundary = false;
		ParticlesContainer[j].resetvalues();
		if (ParticlesContainer[j].isboundary) {
			ParticlesContainer[j].a = boundarya;
		}
		else {
			ParticlesContainer[j].a = 250;
		}
	}
}

void watercolumnsmall(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	// Clear the container
	ParticlesContainer.resize(0);

	// Define boundary and fluid parameters
	int var_boundarypart = 19 * 19  + 15 * 2 * (watercolheight + 20)  + 17 * 2 * (watercolheight + 20);
	gammafloat = 0.7;
	if (!singlewall) {
		//gammafloat = 1;
		var_boundarypart = 19 * 19 *2 + 15 * 2 * (watercolheight + 20) + 17 * 2 * (watercolheight + 20) + 17 * 2 * (watercolheight + 20) + 19 * 2 * (watercolheight + 20);
	}

	// Parameters for the fluid particles
	var_nx = 15;
	var_nz = 15;
	var_fluidpart = var_nx * watercolheight * var_nz;
	var_MaxParticles = var_fluidpart + var_boundarypart;
	hashsize = var_MaxParticles;

	// Resize the container to fit all particles
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	std::srand(std::time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-jitterfac, jitterfac);
	// Create fluid particles
	for (int x = 0; x < var_nx; x++) {
		for (int y = 0; y < watercolheight; y++) {
			for (int z = 0; z < var_nz; z++) {
				double rando = dis(gen);
				double randomnum = jitterfac * static_cast<double>(std::rand()) / (RAND_MAX );
				ParticlesContainer[i].pos = glm::vec3((x+rando) * h, (y + rando) * h, (z + rando) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = false;
				ParticlesContainer[i].index = i;
				i++;
			}
		}
	}
	std::cout << "fluid: " << i << std::endl;
	float boden = 0;
	// Create bottom boundary particles
	for (int x = -2; x <= 16; x++) {
		for (int y = -2; y <= -1; y++) {
			for (int z = -2; z <= 16; z++) {
				if (y > -2 || !singlewall) {
					ParticlesContainer[i].pos = glm::vec3((x) * h, y * h, (z ) * h);
					ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
					ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
					ParticlesContainer[i].a = boundarya;
					ParticlesContainer[i].isboundary = true;
					ParticlesContainer[i].index = i;
					i++;
					boden++;
				}
			}
		}
	}
	std::cout << "boden: " << boden << std::endl;
	// Create front and back single
	float frontandback = 0;
	for (int x = 0; x <= 14; x++) {
		for (int y = 0; y < watercolheight + 20; y++) {
			ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (-1)*h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			frontandback++;
			i++;
			ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (15) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			frontandback++;
			i++;
		}
	}
	std::cout << " vorder und rück wände single : " << frontandback << std::endl;
	// Create front and back double
	frontandback = 0;
	if (!singlewall) {
		for (int x = -1; x < 16; x++) {
			for (int y = 0; y < watercolheight + 20; y++) {
				ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (-2) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				frontandback++;
				ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (16) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				frontandback++;
			}
		}
		std::cout << " vorder und rück wände double : " << frontandback << std::endl;
	}
	int leftandright = 0;
	// Create left and right single
	for (int z = -1; z <= 15; z++) {
		for (int y = 0; y < watercolheight + 20; y++) {
			ParticlesContainer[i].pos = glm::vec3((-1)*h, y * h, (z) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
			ParticlesContainer[i].pos = glm::vec3((15)*h, y * h, (z) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
		}
	}
	std::cout << " left und right wände single : " << leftandright << std::endl;
	// Create left and right double
	leftandright = 0;
	if (!singlewall) {
		for (int z = -2; z <= 16; z++) {
			for (int y = 0; y < watercolheight + 20; y++) {
				ParticlesContainer[i].pos = glm::vec3((-2)*h, y * h, (z) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				leftandright++;
				ParticlesContainer[i].pos = glm::vec3((16)*h, y * h, (z) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				leftandright++;
			}
		}
		std::cout << " left und right wände double : " << leftandright << std::endl;
	}

	// Finalize particle properties
	for (int j = 0; j < var_MaxParticles; j++) {
		ParticlesContainer[j].ismovingboundary = false;
		ParticlesContainer[j].resetvalues();
		if (ParticlesContainer[j].isboundary) {
			ParticlesContainer[j].a = boundarya;
		}
		else {
			ParticlesContainer[j].a = 250;
		}
	}
}

void watercolumnTwoD(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	// Clear the container
	ParticlesContainer.resize(0);

	// Define boundary and fluid parameters
	int var_boundarypart = 19 + 2 * (watercolheight + 100) ;
	gammafloat = 0.7;
	if (!singlewall) {
		//gammafloat = 1;
		var_boundarypart = var_boundarypart = 2* 19 + 2* 2 * (watercolheight + 100);
	}

	// Parameters for the fluid particles
	var_nx = 15;
	var_nz = 1;
	var_fluidpart = var_nx * watercolheight * var_nz;
	var_MaxParticles = var_fluidpart + var_boundarypart;
	hashsize = var_MaxParticles;

	// Resize the container to fit all particles
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	std::srand(std::time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-jitterfac, jitterfac);
	// Create fluid particles
	for (int x = 0; x < var_nx; x++) {
		for (int y = 0; y < watercolheight; y++) {
			double rando = dis(gen);
			double randomnum = jitterfac * static_cast<double>(std::rand()) / (RAND_MAX);
			ParticlesContainer[i].pos = glm::vec3((x + rando) * h, (y + rando) * h,0* h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = 250;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].index = i;
			i++;
		}
	}
	std::cout << "fluid: " << i << std::endl;
	float boden = 0;
	// Create bottom boundary particles
	for (int x = -2; x <= 16; x++) {
		for (int y = -2; y <= -1; y++) {
			
			if (y > -2 || !singlewall) {
				ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (0)*h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				boden++;
				
			}
		}
	}
	std::cout << "boden: " << boden << std::endl;
	int leftandright = 0;
	// Create left and right single

	for (int y = 0; y < watercolheight + 100; y++) {
		ParticlesContainer[i].pos = glm::vec3((-1) * h, y * h, (0)*h);
		ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
		ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
		ParticlesContainer[i].a = boundarya;
		ParticlesContainer[i].isboundary = true;
		ParticlesContainer[i].index = i;
		i++;
		leftandright++;
		ParticlesContainer[i].pos = glm::vec3((15) * h, y * h, (0)*h);
		ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
		ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
		ParticlesContainer[i].a = boundarya;
		ParticlesContainer[i].isboundary = true;
		ParticlesContainer[i].index = i;
		i++;
		leftandright++;
	}

	std::cout << " left und right wände single : " << leftandright << std::endl;
	// Create left and right double
	leftandright = 0;
	if (!singlewall) {
		for (int y = 0; y < watercolheight + 100; y++) {
			ParticlesContainer[i].pos = glm::vec3((-2) * h, y * h, (0) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
			ParticlesContainer[i].pos = glm::vec3((16) * h, y * h, (0) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
		}
		std::cout << " left und right wände double : " << leftandright << std::endl;
	}

	// Finalize particle properties
	for (int j = 0; j < var_MaxParticles; j++) {
		ParticlesContainer[j].ismovingboundary = false;
		ParticlesContainer[j].resetvalues();
		if (ParticlesContainer[j].isboundary) {
			ParticlesContainer[j].a = boundarya;
		}
		else {
			ParticlesContainer[j].a = 250;
		}
	}
}

void DambreaktestTwoD(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	// Clear the container
	ParticlesContainer.resize(0);

	// Define boundary and fluid parameters
	int var_boundarypart = 150 + 2 * (100);
	if (!singlewall) {
		//gammafloat = 1;
		var_boundarypart = var_boundarypart = 2 * 150 + 2 * 2 * (100);
	}

	// Parameters for the fluid particles
	var_nz = 1;
	var_fluidpart = var_nx * var_ny * 1;
	var_MaxParticles = var_fluidpart + var_boundarypart;
	hashsize = var_MaxParticles;

	// Resize the container to fit all particles
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	std::srand(std::time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-jitterfac, jitterfac);
	// Create fluid particles
	for (int x = 0; x < var_nx; x++) {
		for (int y = 0; y < var_ny; y++) {
			double rando = dis(gen);
			ParticlesContainer[i].pos = glm::vec3((x + rando) * h, (y + rando) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = 250;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].index = i;
			i++;
		}
	}
	std::cout << "fluid: " << i << std::endl;
	float boden = 0;
	// Create bottom boundary particles
	for (int x = -2; x <= 148; x++) {
		for (int y = -2; y <= -1; y++) {

			if (y > -2 || !singlewall) {
				ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (0) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				boden++;

			}
		}
	}
	std::cout << "boden: " << boden << std::endl;
	int leftandright = 0;
	// Create left and right single

	for (int y = 0; y <  100; y++) {
		ParticlesContainer[i].pos = glm::vec3((-1) * h, y * h, (0) * h);
		ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
		ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
		ParticlesContainer[i].a = boundarya;
		ParticlesContainer[i].isboundary = true;
		ParticlesContainer[i].index = i;
		i++;
		leftandright++;
		ParticlesContainer[i].pos = glm::vec3((147) * h, y * h, (0) * h);
		ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
		ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
		ParticlesContainer[i].a = boundarya;
		ParticlesContainer[i].isboundary = true;
		ParticlesContainer[i].index = i;
		i++;
		leftandright++;
	}

	std::cout << " left und right wände single : " << leftandright << std::endl;
	// Create left and right double
	leftandright = 0;
	if (!singlewall) {
		for (int y = 0; y <  100; y++) {
			ParticlesContainer[i].pos = glm::vec3((-2) * h, y * h, (0) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
			ParticlesContainer[i].pos = glm::vec3((148) * h, y * h, (0) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
		}
		std::cout << " left und right wände double : " << leftandright << std::endl;
	}

	// Finalize particle properties
	for (int j = 0; j < var_MaxParticles; j++) {
		ParticlesContainer[j].ismovingboundary = false;
		ParticlesContainer[j].resetvalues();
		if (ParticlesContainer[j].isboundary) {
			ParticlesContainer[j].a = boundarya;
		}
		else {
			ParticlesContainer[j].a = 250;
		}
	}
}




void moving_boundary(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	ParticlesContainer.resize(0);
	var_nx = var_oneway - 2;
	var_nz = var_oneway - 2;
	var_ny = 15;
	var_fluidpart = var_nx * var_ny * var_nz - std::min(num_rot_part, var_ny * num_rot_part/8) ;
	if (singlewall) {
		gammafloat = 0.7;
		var_MaxParticles = var_fluidpart + (var_oneway * var_oneway * 6) + num_rot_part;
	}
	else {
		gammafloat = 1;
		var_MaxParticles = var_nx * var_ny * var_nz + (var_oneway * var_oneway * 12)+ num_rot_part;
	}

	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	for (int x = -num_rot_part / 8; x < (var_nx- num_rot_part / 8); x++) {
		for (int y = 0; y < (var_ny); y++) {
			for (int z = -num_rot_part / 8; z < (var_nz- num_rot_part / 8); z++) {
				if (!(((x+1) >= 0 && (x+1) < num_rot_part / 8) && ((y+1) >= 2 && (y+1) < 10) && (z+1) == 0)) {
					////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
					ParticlesContainer[i].pos = glm::vec3((x + 1) * (h * 1), (y + 1) * (h * 1), (z + 1) * (h * 1));
					ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
					ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
					ParticlesContainer[i].a = 250;
					ParticlesContainer[i].isboundary = 0;
					ParticlesContainer[i].index = i;
					i += 1;
				}
				
			}
		}
	}
	std::cout << "numfluid" << i << std::endl;
	std::cout << "num oneway " << var_oneway << std::endl;
	std::cout << "numfluid diff " << i- ((var_oneway- 2)* (var_oneway-2) * 10 - num_rot_part) << std::endl;

	int nfluid = i;
	//moving stuff
	for (int ii = 0; ii < num_rot_part/8; ++ii) {
		for (int jj = 0; jj < 8; ++jj){
			glm::vec3 velocity(deltaT/2, deltaT / 2, 0.0f); // Geschwindigkeit in die korrekte Richtung
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].vel = velocity;
			ParticlesContainer[i].pos = glm::vec3((0 + ii) * h, (2 + jj) * h, 0 * h);
			ParticlesContainer[i].ismovingboundary = true;
			i += 1;
		}
		
	}
	std::cout << "nummoving" << i- nfluid << std::endl;
	int nmbound = i- nfluid;
	//border
	for (int ii = 0; ii < var_oneway; ii++) {
		for (int j = 0; j < var_oneway; j++) {
			//1upper
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii- num_rot_part / 8) * h, var_oneway * h, (j- num_rot_part / 8) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii- num_rot_part / 8) * h, j * h,- num_rot_part / 8 *h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//3right
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((var_oneway - 1- num_rot_part / 8) * h, ii * h, (j- num_rot_part / 8) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii- num_rot_part / 8) * h, 0, (j- num_rot_part / 8) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//5front
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii- num_rot_part / 8) * h, j * h, (var_oneway - 1 - num_rot_part / 8) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((0- num_rot_part / 8)*h, ii * h, (j- num_rot_part / 8) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//7back2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1- num_rot_part / 8) * h, (j - 1) * h, (-num_rot_part / 8 - 1));
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//8right2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((var_oneway - 2- num_rot_part / 8) * h, (ii - 1) * h, (j - 1 - num_rot_part / 8) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//9lower2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1- num_rot_part / 8) * h, (0 - 1), (j - 1- num_rot_part / 8) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//10front2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1- num_rot_part / 8) * h, (j - 1) * h, (var_oneway - 0- num_rot_part / 8) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//11 left2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((0 - 1- num_rot_part / 8), (ii - 1) * h, (j - 1- num_rot_part / 8) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//12upper2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1- num_rot_part / 8) * h, (var_oneway - 1) * h, (j - 1- num_rot_part / 8) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	std::cout << "numbound" << i - nfluid- nmbound << std::endl;
	int nbound = i - nfluid- nmbound;
	int totalp = i;
	std::cout << "numdiff" << var_MaxParticles - nfluid - nmbound - nbound << std::endl;

	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
		if (ParticlesContainer[i].ismovingboundary) {
			ParticlesContainer[i].a = 50;
		}
	}
}
void dambreaktest(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	ParticlesContainer.resize(0);
	var_nx = 54;
	var_ny = 30;
	var_nz = 68;
	var_fluidpart = var_nx * var_ny * var_nz;
	if (singlewall) {
		gammafloat = 0.7;
		var_MaxParticles = var_nx * var_ny * var_nz + (180*70 + 2* 180*60 + 2* 70*60);
	}
	else {
		gammafloat = 1;
		var_MaxParticles = var_nx * var_ny * var_nz + 2*(180 * 70 + 2 * 180 * 60 + 2 * 70 * 60);
	}

	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	for (int x = 0; x < (var_nx); x++) {
		for (int y = 0; y < (var_ny); y++) {
			for (int z = 0; z < (var_nz); z++) {
				////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 1) * (h * 1), (y + 1) * (h * 1), (z + 1) * (h * 1));
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	//border back and front
	for (int ii = 0; ii < 180; ii++) {
		for (int j = 0; j < 60; j++) {
			//2back
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//5front
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, (70 - 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//7back2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (j - 1) * h, (0 - 1));
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//10front2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (j - 1) * h, (70 -0) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	//border left and right 
	for (int ii = 0; ii < 60; ii++) {
		for (int j = 0; j < 70; j++) {
			//3right
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((180 - 1) * h, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(0, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//8right2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((180 - 2) * h, (ii - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//11 left2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((0 - 1), (ii - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	// border bottom
	for (int ii = 0; ii < 180; ii++) {
		for (int j = 0; j < 70; j++) {

			//4lower
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, 0, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//9lower2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii ) * h, (0 - 1)*h, (j ) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
}
void smalldambreaktest(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	ParticlesContainer.resize(0);
	var_nx = 54;
	var_ny = 30;
	var_nz = 33;
	var_fluidpart = var_nx * var_ny * var_nz;
	if (singlewall) {
		gammafloat = 1;
		var_MaxParticles = var_nx * var_ny * var_nz + (180 * 35 + 2 * 180 * 60 + 2 * 35 * 60)+ (2*10*20+10*20+2*10*10);
	}
	else {
		gammafloat = 1;
		var_MaxParticles = var_nx * var_ny * var_nz + 2 * (180 * 35 + 2 * 180 * 60 + 2 * 35 * 60) + (2 * 10 * 20 + 10 * 20 + 2 * 10 * 10);
	}

	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	for (int x = 0; x < (var_nx); x++) {
		for (int y = 0; y < (var_ny); y++) {
			for (int z = 0; z < (var_nz); z++) {
				////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 2) * (h * 1), (y + 1) * (h * 1), (z + 1) * (h * 1));
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	//border back and front
	for (int ii = 0; ii < 180; ii++) {
		for (int j = 0; j < 60; j++) {
			//2back
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//5front
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, (35 - 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//7back2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (j - 1) * h, (0 - 1));
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//10front2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (j - 1) * h, (35 - 0) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	//border left and right 
	for (int ii = 0; ii < 60; ii++) {
		for (int j = 0; j < 35; j++) {
			//3right
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((180 - 1) * h, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(1, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//8right2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((180 - 2) * h, (ii - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//11 left2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((0 ), (ii - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	// border bottom
	for (int ii = 0; ii < 180; ii++) {
		for (int j = 0; j < 35; j++) {

			//4lower
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, 0, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//9lower2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii)*h, (0 - 1) * h, (j)*h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	// box
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 20; j++) {
			//upper
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii+120) * h, 10+1*h, (j+5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//8right2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((130 - 1) * h, (ii+1 ) * h, (j +5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//11 left2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((120)*h, (ii+1 ) * h, (j +5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (j < 10) {
				//2back
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii+120) * h, (j+1) * h, 5*h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//5front
				////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii+120) * h, (j+1) * h, (24) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
}