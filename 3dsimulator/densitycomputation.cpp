#include <densitycomputation.h>
#include <iisphparticle.h>
#include <part.h>

void computeAllDens(std::vector<iisphparticle>& var_PartC) {
	denserrold = 0;
	numofp0high = 0;
	usemefordens.clear();
	if (usewholefluidforcalc == false) {
		for (int i = 0; i < var_PartC.size(); ++i) {
			iisphparticle& Part = var_PartC[i];
			Part.density = 0;
			if (Part.isboundary == false /* || Part.ismovingboundary == true*/) {
				Part.denstolow = true;
				float dens = 0;
				for (float& kern : Part.Kernel) {
					dens += kern * Part.m;
				}
				Part.density = dens;


				if (dens >= p0 && Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord && Part.isfloatingboundary == false && Part.ismovingboundary == false) {
					Part.denstolow = false;
					numofp0high++;
					if (absinterrupt) {
						denserrold += 100 * (std::abs(Part.density - p0) / p0);
						usemefordens.push_back(Part.index);
					}
					else {
						denserrold += 100 * (Part.density - p0) / p0;
						usemefordens.push_back(Part.index);
					}
				}
			}
			else {
				Part.density = p0;
				Part.denstolow = true;
			}
		}
	}
	else {
		for (int i = 0; i < var_PartC.size(); ++i) {
			iisphparticle& Part = var_PartC[i];
			Part.density = 0;
			if (Part.isboundary == false /* || Part.ismovingboundary == true*/) {
				Part.denstolow = true;
				float dens = 0;
				for (float& kern : Part.Kernel) {
					dens += kern * Part.m;
				}
				Part.density = dens;


				if (Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord && Part.isfloatingboundary == false && Part.ismovingboundary == false) {
					Part.denstolow = false;
					if (Part.density >= p0) {
						numofp0high++;
					}
					if (absinterrupt) {
						denserrold += 100 * (std::abs(Part.density - p0) / p0);
						usemefordens.push_back(Part.index);
					}
					else {
						denserrold += 100 * (Part.density - p0) / p0;
						usemefordens.push_back(Part.index);
					}
				}
			}
			else {
				Part.density = p0;
				Part.denstolow = true;
			}
		}
	}
	if (usewholefluidtodevide) {
		denserrold = (denserrold /var_fluidpart);
	}
	else {
		denserrold = (denserrold / usemefordens.size());
	}
}
void computeAllDenstwoD(std::vector<iisphparticle>& var_PartC) {
	denserrold = 0;
	numofp0high = 0;
	usemefordens.clear();
	if (usewholefluidforcalc == false) {
		for (int i = 0; i < var_PartC.size(); ++i) {
			iisphparticle& Part = var_PartC[i];
			Part.density = 0;
			if (Part.isboundary == false) {
				Part.denstolow = true;
				float dens = 0;

				for (float& kern : Part.Kernel) {
					dens += kern * Part.m;
				}
				Part.density = dens;
				if (dens >= p0 && Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord&& Part.isfloatingboundary == false) {
					Part.denstolow = false;
					numofp0high++;
					if (absinterrupt) {
						denserrold += 100 * (std::abs(Part.density - p0) / p0);
						usemefordens.push_back(Part.index);
					}
					else {
						denserrold += 100 * (Part.density - p0) / p0;
						usemefordens.push_back(Part.index);
					}
				}
			}
			else {
				Part.density = p0;
				Part.denstolow = true;
			}

		}
	}
	else {
		for (int i = 0; i < var_PartC.size(); ++i) {
			iisphparticle& Part = var_PartC[i];
			Part.density = 0;
			if (Part.isboundary == false) {
				Part.denstolow = true;
				float dens = 0;

				for (float& kern : Part.Kernel) {
					dens += kern * Part.m;
				}
				Part.density = dens;
				if (Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord && Part.isfloatingboundary == false) {
					Part.denstolow = false;
					if (Part.density >= p0) {
						numofp0high++;
					}
					if (absinterrupt) {
						denserrold += 100 * (std::abs(Part.density - p0) / p0);
						usemefordens.push_back(Part.index);
					}
					else {
						denserrold += 100 * (Part.density - p0) / p0;
						usemefordens.push_back(Part.index);
					}
				}
			}
			else {
				Part.density = p0;
				Part.denstolow = true;
			}

		}
	}
	if (usewholefluidtodevide) {
		denserrold = (denserrold / var_fluidpart);
	}
	else {
		denserrold = (denserrold / usemefordens.size());
	}
}

void computeAllDensSSPH(std::vector<iisphparticle>& var_PartC) {
	denserrold = 0;
	numofp0high = 0;
	//#pragma omp parallel for
	for (auto& Part : var_PartC) {
		if (Part.isboundary == false || Part.isfloatingboundary) {
			float dens = 0;
			for (float& kern : Part.Kernel) {
				dens += kern * Part.m;
			}
			dens = std::max(dens, p0);
			Part.density = dens;
			Part.denstolow = false;
			numofp0high++;
			if (absinterrupt) {
				denserrold += std::abs(Part.density - p0);
			}
			else {
				denserrold += Part.density - p0;
			}
		}
		if (Part.isboundary == true && Part.isfloatingboundary == false) {
			Part.density = p0;
		}
	}
	denserrold = 100 * (denserrold / numofp0high) / p0;
}
void computeAllDensErr(std::vector<iisphparticle>& var_PartC) {
#pragma omp parallel for
	for (int i = 0; i < var_MaxParticles; i++) {
		iisphparticle& Part = var_PartC[i];
		if (Part.isboundary == false ) {
			Part.densityerror = (Part.AP - Part.sf) / p0;
		}
	}
}


float makeAlldenserrAvg(std::vector<iisphparticle>& var_PartC) {
	float densityerroraverage = 0.f;
	max_singledens = 0;
#pragma omp parallel for reduction(+:densityerroraverage)
	for (int i = 0; i < usemefordens.size(); ++i) {
		iisphparticle& Part = var_PartC[usemefordens[i]];
		if (absinterrupt) {
			densityerroraverage += 100 * std::abs(Part.densityerror);
			max_singledens = std::max(max_singledens, std::abs(Part.densityerror));
		}
		else {
			densityerroraverage += 100 * Part.densityerror;
			max_singledens = std::max(max_singledens, std::abs(Part.densityerror));
		}
	}
	max_singledens = 100 * max_singledens;
	if (usewholefluidtodevide) {
		return densityerroraverage /var_fluidpart;
	}
	else {
		return densityerroraverage / usemefordens.size();
	}
	return 0;
}

