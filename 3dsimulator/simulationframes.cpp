#include <simulationframes.h>
#include <iisphparticle.h>
#include <part.h>


void ssphAlgo(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap) {
	auto start_o = std::chrono::high_resolution_clock::now();
	makecfltrue(PartC);
	//deltaT = glm::min(deltaT, deltaTmax);
	auto start_n = std::chrono::high_resolution_clock::now();
	if (uniformgridneighsearch) {
		clearuniformgrid();
		insertparticlesinuniformgrid(PartC);
		findAllNeighbourscompact3D(PartC);
	}
	else {
		hashmap.clear();
		insertAllParticlesIntoHashmap(PartC, hashmap);
		findAllNeighbours(PartC, hashmap);
	}
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
	if (uniformgridneighsearch) {
		clearuniformgrid();
		insertparticlesinuniformgrid(PartC);
		findAllNeighbourscompact2D(PartC);
	}
	else {
		hashmap.clear();
		insertAllParticlesIntoHashmap2D(PartC, hashmap);
		findAllNeighbours2D(PartC, hashmap);
	}
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



void var_iisph(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap) {
	auto start_o = std::chrono::high_resolution_clock::now();
	makecfltrue(PartC);
	auto start_n = std::chrono::high_resolution_clock::now();
	if (uniformgridneighsearch) {
		clearuniformgrid();
		insertparticlesinuniformgrid(PartC);
		findAllNeighbourscompact3D(PartC);
	}
	else {
		hashmap.clear();
		insertAllParticlesIntoHashmap(PartC, hashmap);
		findAllNeighbours(PartC, hashmap);
	}
	
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
	convergence.clear();
	convergence.push_back(makeAlldenserrAvg(PartC));
	int l = 1;
	while (makeAlldenserrAvg(PartC) > denistyerrormax || l < currentitermin/* || 100 * max_singledens / p0 > 50*/) {
		//first for
		if (pressuredextrapolatebound) {
			makeBoundPres(PartC);
			makeAllPresAwithbound(PartC);
		}
		else {
			makeAllPresA(PartC);
		}
		//second for
		secondloop(PartC);
		l++;
		if (l > currentitermax) {
			//paussimul = true;
			break;
		}
		convergence.push_back(makeAlldenserrAvg(PartC));
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


void twoDiisph(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap) {
	auto start_o = std::chrono::high_resolution_clock::now();
	makecfltrue(PartC);

	auto start_n = std::chrono::high_resolution_clock::now();
	if (uniformgridneighsearch) {
		clearuniformgrid();
		insertparticlesinuniformgrid(PartC);
		findAllNeighbourscompact2D(PartC);
	}
	else {
		hashmap.clear();
		insertAllParticlesIntoHashmap2D(PartC, hashmap);
		findAllNeighbours2D(PartC, hashmap);
	}	
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
	convergence.clear();
	convergence.push_back(makeAlldenserrAvg(PartC));
	int l = 1;
	while (makeAlldenserrAvg(PartC) > denistyerrormax || l < currentitermin || max_singledens > 100) {
		if (pressuredextrapolatebound) {
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
		convergence.push_back(makeAlldenserrAvg(PartC));
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


