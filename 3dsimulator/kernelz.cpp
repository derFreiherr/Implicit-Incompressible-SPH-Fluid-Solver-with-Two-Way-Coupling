#include <kernelz.h>
#include <iisphparticle.h>
#include <part.h>


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
					Part.IdNSubKernelder.push_back(std::make_tuple(std::get<0>(neig), std::get<2>(neig), 0.999221087f * alpha * std::get<2>(neig) / (d * h * h) * (-3 * t2 * t2 + 12 * t1 * t1)));
				}
			}
			Part.Kernel.clear();
			//#pragma omp parallel for
			for (const auto& neig : Part.IdNdistNsub) {
				float d = std::get<1>(neig) / h;
				float t1 = std::max((1 - d), 0.f);
				float t2 = std::max((2 - d), 0.f);
				if (var_PartC[std::get<0>(neig)].isboundary) {
					Part.Kernel.push_back(0.999221087 * gammadens * alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
				}
				else {
					Part.Kernel.push_back(0.999221087 * alpha * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)));
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
	return alphaTwoD * ((t2 * t2 * t2) - (4 * t1 * t1 * t1)) * 0.99911389780816045467516953288827f;
}


glm::vec3 makesinglekernelder(glm::vec3& posi, glm::vec3& posj) {
	float distX = posi.x - posj.x;
	float distY = posi.y - posj.y;
	float distZ = posi.z - posj.z;
	float distance = sqrt(distX * distX + distY * distY + distZ * distZ);
	float d = distance / h;
	float t1 = std::max((1 - d), 0.f);
	float t2 = std::max((2 - d), 0.f);
	return alphaTwoD * (posi - posj) / (d * h * h) * (-3 * t2 * t2 + 12 * t1 * t1);
}