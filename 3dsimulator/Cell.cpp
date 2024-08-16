#include <Cell.h>
#include <iisphparticle.h>
#include <part.h>





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
		if (Part.isboundary == false  || Part.isfloatingboundary) {
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
				Part.nonpresA = (.01f * ViscAb + .01f * ViscAf) + gravity * glm::vec3(0, 1, 0);
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
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		if ((Part.isboundary == false || Part.isfloatingboundary)) {
			glm::vec3 ViscAf(0.f, 0.f, 0.f);
			glm::vec3 ViscAb(0.f, 0.f, 0.f);
			glm::vec3 surfacetens(0.f, 0.f, 0.f);
			if (Part.isfloatingboundary) {
				for (const auto& neig : Part.IdNSubKernelder) {
					if (var_PartC[std::get<0>(neig)].isboundary && !var_PartC[std::get<0>(neig)].ismovingboundary) {
						ViscAb += (var_PartC[std::get<0>(neig)].m / p0) * (Part.vel * makesinglekernel2D(Part.pos, var_PartC[std::get<0>(neig)].pos));
					}
					else {
						ViscAf += (var_PartC[std::get<0>(neig)].m / p0) * ((Part.vel - var_PartC[std::get<0>(neig)].vel) * makesinglekernel(Part.pos, var_PartC[std::get<0>(neig)].pos));
					}
				}
				Part.nonpresA = ((10.f * ViscAb + 1.f * ViscAf) + gravity * glm::vec3(0, 1, 0))* glm::vec3(1, 1, 0);
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
				Part.nonpresA = ((visc_boundary * ViscAb + visc_fluid * ViscAf) + gravity * glm::vec3(0, 1, 0) - ((surfacetension / Part.m) * surfacetens))*glm::vec3(1,1,0);
			}
		}
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
			var_PartC[i].sf = p0  - var_PartC[i].density - (deltaT * tmpf) - (deltaT * tmpb);
		}
	}
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
				APB += var_PartC[neig_idx].m * glm::dot(Part.presA - var_PartC[neig_idx].presA, kernel_deriv);
				/*
				if (var_PartC[neig_idx].isboundary) {
					APB += var_PartC[neig_idx].m * glm::dot(Part.presA - var_PartC[neig_idx].presA, kernel_deriv);
				}
				else {
					APF += var_PartC[neig_idx].m * glm::dot((Part.presA - var_PartC[neig_idx].presA), kernel_deriv);
				}
				*/
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
			if (Part.pos.y > upperviualbord || Part.pos.y < lowervisualbord) {
				Part.pos.y = upperviualbord + 1;
				Part.pos.x = 0;
				Part.pos.z = 0;
				Part.vel = glm::vec3(0, 0, 0);
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


