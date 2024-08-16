#include <rigidbodies.h>
#include <iisphparticle.h>
#include <part.h>

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
			//Part.m = Part.m*massfac;
			allrigidmass += Part.m * massfac;
			xCM += Part.m * massfac * Part.pos;
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
			float mass = Part.m *massfac;
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
	torque = glm::vec3(0.f, 0.f, 0.f);
	glm::vec3 linforce(0.f, 0.f, 0.f);
	// Parallelisierte Schleife mit OpenMP
//#pragma omp parallel for private(torque, linforce)
	for (int i = var_fluidpart; i < PartC.size(); ++i) {
		iisphparticle& Part = PartC[i];
		if (Part.isfloatingboundary) {
			glm::vec3 force = ((Part.presA) + Part.nonpresA * massfac) * Part.m ; // Annahme, dass Kräfte Beschleunigungen sind
			glm::vec3 r = Part.pos - xCM;
			torque += glm::cross(r, force);;
			linforce += force;
		}
	}


	//torque = torque * 0.9f;

	vCM += (deltaT * linforce / allrigidmass);
	xCM += (deltaT * vCM);

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
			Part.pos = (xCM + relpostmp); // Position mit linearer Geschwindigkeit aktualisieren
			Part.vel = (vCM + glm::cross(omegarigidbody, relpostmp));
		}
	}
}
void updaterigidbody2d(std::vector<iisphparticle>& PartC) {
	glm::vec3 torguefac = torque;
	torque = glm::vec3(0.f, 0.f, 0.f);
	glm::vec3 linforce(0.f, 0.f, 0.f);
	// Parallelisierte Schleife mit OpenMP
//#pragma omp parallel for private(torque, linforce)
	for (int i = var_fluidpart; i < PartC.size(); ++i) {
		iisphparticle& Part = PartC[i];
		if (Part.isfloatingboundary) {
			glm::vec3 force = ((Part.presA) + Part.nonpresA * massfac) * Part.m ; // Annahme, dass Kräfte Beschleunigungen sind
			glm::vec3 r = Part.pos - xCM;
			torque += glm::cross(r, force);;
			linforce += force;
		}
	}


	torque = torque * 0.1f;
	//torque = -torque;
	vCM += glm::vec3(1, 1, 0) * (deltaT * linforce / (allrigidmass));
	xCM += glm::vec3(1, 1, 0) * (deltaT * vCM);

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
			Part.pos = glm::vec3(1, 1, 0) * (xCM + relpostmp); // Position mit linearer Geschwindigkeit aktualisieren
			Part.vel = glm::vec3(1, 1, 0) * (vCM + glm::cross(omegarigidbody, relpostmp));
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
			Part.m = 0.9 * p0 / std::max(kernsum, 0.0000000000000001f);
		}

	}
}
