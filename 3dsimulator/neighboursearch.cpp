#include <neighboursearch.h>
#include <iisphparticle.h>
#include <part.h>
#include<omp.h>
//uniform grid
void clearuniformgrid() {
#pragma omp parallel for
	for (int i = 0; i < uniformgidvec1D.size(); i++) {
		uniformgidvec1D[i].clear();
	}
}

void insertparticlesinuniformgrid(std::vector<iisphparticle>& var_PartC) {
	for (int i = 0; i < var_MaxParticles; i++) {
		iisphparticle& Part = var_PartC[i];
		uniformgidvec1D[uniformgridhash(Part.pos.x, Part.pos.y, Part.pos.z)].push_back(Part.index);
	}
}
unsigned int uniformgridhash(float x, float y, float z) {
	int xi = static_cast<int>(std::floor((x-(overallminxpos-searchRadius)) / cellsize));
	int yi = static_cast<int>(std::floor((y - (overallminypos-searchRadius)) / cellsize));
	int zi = static_cast<int>(std::floor((z - (overallminzpos-searchRadius)) / cellsize));
	return static_cast<int>(xi + yi * gridbreite + zi * gridhöhe*gridbreite);
}

void findAllNeighbourscompact2D(std::vector<iisphparticle>& var_PartC) {
	totalcomp = 0;
#pragma omp parallel for
	for (int i = 0; i < (var_MaxParticles); ++i) {
		iisphparticle& Part = var_PartC[i];
		Part.IdNdistNsub.clear();
		Part.drawme = false;
		Part.computeme = false;
		Part.Aff = 0;
		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				int cellindex = uniformgridhash(Part.pos.x + (x * searchRadius), Part.pos.y + (y * searchRadius), 0.f);
				for (int idx : uniformgidvec1D[cellindex]) {
					float distX = var_PartC[idx].pos.x - Part.pos.x;
					float distY = var_PartC[idx].pos.y - Part.pos.y;
					float distance = sqrt(distX * distX + distY * distY);
					if (distance <= searchRadius) {
						Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));

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
void findAllNeighbourscompact3D(std::vector<iisphparticle>& var_PartC) {
	totalcomp = 0;
#pragma omp parallel for
	for (int i = 0; i < (var_MaxParticles); ++i) {
		iisphparticle& Part = var_PartC[i];
		Part.IdNdistNsub.clear();
		Part.drawme = false;
		Part.computeme = false;
		Part.Aff = 0;
		//std::unordered_set<int> uniqueIds;
		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				for (int z = -1; z <= 1; z++) {
					int cellindex = uniformgridhash(Part.pos.x + (x * searchRadius), Part.pos.y + (y * searchRadius), Part.pos.z + (z * searchRadius));
					//std::cout << cellindex << "   size" << uniformgidvec1D.size() << std::endl;
					for (int idx : uniformgidvec1D[cellindex]) {
						//if (uniqueIds.find(idx) == uniqueIds.end()) {
							float distX = var_PartC[idx].pos.x - Part.pos.x;
							float distY = var_PartC[idx].pos.y - Part.pos.y;
							float distZ = var_PartC[idx].pos.z - Part.pos.z;
							float distance = sqrt(distX * distX + distY * distY + distZ *distZ);
							if (distance <= searchRadius) {
								Part.IdNdistNsub.push_back(std::make_tuple(idx, distance, Part.pos - var_PartC[idx].pos));
								//uniqueIds.insert(idx);
							}
						//}
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

//spatial hashing
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
		if (Part.IdNdistNsub.size() < animationneighbourscount && (Part.isboundary == false || Part.ismovingboundary == true)) {
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
		if (Part.IdNdistNsub.size() < animationneighbourscount) {
			Part.drawme = true;
		}
		if (Part.IdNdistNsub.size() > compbord) {
			Part.computeme = true;
			totalcomp++;
		}

	}
}

