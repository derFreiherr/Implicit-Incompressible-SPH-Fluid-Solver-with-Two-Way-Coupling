#include<boundarypres.h>
#include <iisphparticle.h>
#include <part.h>

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
		Part.pressureiter = (firstsum + secondsum) / thirdsum;
	}
}
void makeBoundPres2D(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 
		
		+ var_spezialboundpart; i < var_MaxParticles; i++) {
		iisphparticle& Part = var_PartC[i];
		float firstsum = 0;
		float secondsum = 0;
		float thirdsum = 0;
		for (auto& neigh : Part.IdNdistNsub) {
			if (var_PartC[std::get<0>(neigh)].isboundary == false) {
				firstsum += var_PartC[std::get<0>(neigh)].pressureiter * makesinglekernel2D(Part.pos, var_PartC[std::get<0>(neigh)].pos);
				secondsum += var_PartC[std::get<0>(neigh)].density * std::get<1>(neigh) * makesinglekernel2D(Part.pos, var_PartC[std::get<0>(neigh)].pos);
				thirdsum += makesinglekernel2D(Part.pos, var_PartC[std::get<0>(neigh)].pos);
			}
		}
		Part.pressureiter = (firstsum + secondsum) / thirdsum;
	}
}

