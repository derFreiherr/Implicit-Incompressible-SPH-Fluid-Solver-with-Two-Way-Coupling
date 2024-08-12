#include <Affz.h>
#include <iisphparticle.h>
#include <part.h>

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
			for (std::tuple<int, glm::vec3, glm::vec3>& jf : PartC[i].IdNSubKernelder) {
				if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != PartC[i].index) {
					glm::vec3  jjf = glm::vec3(0.f, 0.f, 0.f);
					glm::vec3  jjb = glm::vec3(0.f, 0.f, 0.f);
					// sum over all neigbours jj from if exept if
					for (std::tuple<int, glm::vec3, glm::vec3>& jj : PartC[i].IdNSubKernelder) {
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
			for (std::tuple<int, glm::vec3, glm::vec3>& jb : PartC[i].IdNSubKernelder) {
				if (PartC[std::get<0>(jb)].isboundary) {
					glm::vec3  jjf = glm::vec3(0.f, 0.f, 0.f);
					glm::vec3  jjb = glm::vec3(0.f, 0.f, 0.f);
					// sum over all neigbours jj from if
					for (std::tuple<int, glm::vec3, glm::vec3>& jj : PartC[i].IdNSubKernelder) {
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
			for (std::tuple<int, glm::vec3, glm::vec3>& jf : PartC[i].IdNSubKernelder) {
				if (PartC[std::get<0>(jf)].isboundary == false && PartC[std::get<0>(jf)].index != PartC[i].index) {
					glm::vec3 Kerndelder_jf_if = glm::vec3(0, 0, 0);
					// sum over all neigbours jj from jf
					for (std::tuple<int, glm::vec3, glm::vec3>& tmp : PartC[std::get<0>(jf)].IdNSubKernelder) {
						// if jj == if
						if (PartC[std::get<0>(tmp)].index == PartC[i].index) {
							Kerndelder_jf_if = std::get<2>(tmp);
							// sum ( jf_m * (if_m/(p0*p0) * dWdt_jf_if)* dWdt_if_jf
							if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((PartC[i].m / (p0p0)*Kerndelder_jf_if), std::get<2>(jf));
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
							if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((PartC[i].m / (p0p0)*Kerndelder_jf_if), std::get<2>(jf));
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
					if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((PartC[i].m / (PartC[i].density * PartC[i].density) * -std::get<2>(jf)), std::get<2>(jf));
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
	for (int i = var_fluidpart; i < (var_fluidpart + var_spezialboundpart); ++i) {
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
		Part.Aff = (deltaT * deltaT) * (if_jf1 + if_jf2 + if_jb);
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
							if_jf2 += PartC[std::get<0>(jf)].m * glm::dot((Part.m / (p0p0)*Kerndelder_jf_if), std::get<2>(jf));
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
