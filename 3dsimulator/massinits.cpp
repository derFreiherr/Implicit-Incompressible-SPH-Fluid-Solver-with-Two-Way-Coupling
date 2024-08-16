#include <massinits.h>
#include <iisphparticle.h>
#include <part.h>

void makeboundmass(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
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
				Part.Kernel.push_back(alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				kernsum += alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1));
			}
			Part.m = gammabound * p0 / std::max(kernsum, 0.0000000000000001f);
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
				Part.Kernel.push_back(alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				kernsum += alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)) * 0.9991389780816045467516953288829;
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
				Part.Kernel.push_back(alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				kernsum += alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1));
			}
			Part.m = gammabound * p0 / std::max(kernsum, 0.0000000000000001f);
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
				Part.Kernel.push_back(alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				kernsum += alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1))* 0.99913884;
			}// 0.9991389780816045467516953288829
			Part.m = p0 / std::max(gammapart * kernsum, 0.0000000000000001f);
		}

	}
}
