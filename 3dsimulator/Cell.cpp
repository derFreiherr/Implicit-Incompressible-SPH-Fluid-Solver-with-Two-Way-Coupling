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
				Part.nonpresA = ((1.f * ViscAb + 1.f * ViscAf) + gravity * glm::vec3(0, 1, 0))* glm::vec3(1, 1, 0);
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


