#include <densitycomputation.h>
#include <iisphparticle.h>
#include <part.h>

void computeAllDens(std::vector<iisphparticle>& var_PartC) {
	denserrold = 0;
	numofp0high = 0;
	usemefordens.clear();
	//float deabugdens = 0;
//#pragma omp parallel for private(denserrold, numofp0high, usemefordens)
	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		Part.density = 0;
		if (Part.isboundary == false /* || Part.ismovingboundary == true*/) {
			Part.denstolow = true;
			float dens = 0;
			
			for (float& kern : Part.Kernel) {
				dens += kern * Part.m;
			}
			/*
			for (const auto& neig : Part.IdNdistNsub) {
				dens += var_PartC[std::get<0>(neig)].m * makesinglekernel(Part.pos, var_PartC[std::get<0>(neig)].pos);
			}
			*/
			Part.density = dens;
			if (dens >= p0 && Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord && Part.isfloatingboundary == false && Part.ismovingboundary==false) {
				Part.denstolow = false;
				numofp0high++;
				if (absinterrupt) {
					denserrold += 100 * (std::abs(Part.density - p0) / p0);
					usemefordens.push_back(Part.index);
					//deabugdens = std::max(deabugdens, Part.density);

				}
				else {
					denserrold += 100 * (Part.density - p0) / p0;
					usemefordens.push_back(Part.index);
				}
			}

			/*
			if (clampp0 == true) {
				Part.density = std::max(dens, p0 - p0 * denistyerrormax / 100 * clampfac);
				Part.denstolow = false;
				numofp0high++;
				if (absinterrupt) {
					denserrold += std::abs(Part.density - p0);
				}
				else {
					denserrold += Part.density - p0;
				}
			}
			if (Part.denstolow == true && clamptolow == true) {
				Part.density = std::max(dens, p0 - p0 * denistyerrormax / 100 * clampfac);
			}
			*/
		}
		if (Part.isboundary == true/* && Part.ismovingboundary == false*/) {
			Part.density = p0;
			Part.denstolow = true;
		}

	}
	//std::cout << deabugdens << "   "<< usemefordens.size()<< std::endl;
	denserrold = (denserrold / usemefordens.size());
	//denserrold = (denserrold / var_fluidpart);
}
void computeAllDenstwoD(std::vector<iisphparticle>& var_PartC) {
	denserrold = 0;
	numofp0high = 0;
	usemefordens.clear();
	//float deabugdens = 0;

	for (int i = 0; i < var_PartC.size(); ++i) {
		iisphparticle& Part = var_PartC[i];
		Part.density = 0;
		if (Part.isboundary == false) {
			Part.denstolow = true;
			float dens = 0;
			
			for (float& kern : Part.Kernel) {
				dens += kern * Part.m;
			}
			/*
			for (const auto& neig : Part.IdNdistNsub) {
				dens += var_PartC[std::get<0>(neig)].m * makesinglekernel2D(Part.pos, var_PartC[std::get<0>(neig)].pos);
			}
			*/
			Part.density = dens;
			if (dens >= p0 && Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord&& Part.isfloatingboundary == false) {
				Part.denstolow = false;
				numofp0high++;
				if (absinterrupt) {
					denserrold += 100 * (std::abs(Part.density - p0) / p0);
					usemefordens.push_back(Part.index);
					//deabugdens = std::max(deabugdens, Part.density);

				}
				else {
					denserrold += 100 * (Part.density - p0) / p0;
					usemefordens.push_back(Part.index);
				}
			}

			/*
			if (clampp0 == true) {
				Part.density = std::max(dens, p0 - p0 * denistyerrormax / 100 * clampfac);
				Part.denstolow = false;
				numofp0high++;
				if (absinterrupt) {
					denserrold += std::abs(Part.density - p0);
				}
				else {
					denserrold += Part.density - p0;
				}
			}
			if (Part.denstolow == true && clamptolow == true) {
				Part.density = std::max(dens, p0 - p0 * denistyerrormax / 100 * clampfac);
			}
			*/
		}
		if (Part.isboundary == true) {
			Part.density = p0;
			Part.denstolow = true;
		}

	}
	//std::cout << deabugdens << "   "<< usemefordens.size()<< std::endl;
	denserrold = (denserrold / usemefordens.size());
	//denserrold = (denserrold / var_fluidpart);
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
	if (ignorep0tolow) {
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
		return densityerroraverage / usemefordens.size();
		//return densityerroraverage / var_fluidpart;
	}
	/*
	if (ignoreincomplete) {
#pragma omp parallel for reduction(+:densityerroraverage)
		for (int i = 0; i < var_MaxParticles; ++i) {
			iisphparticle& Part = var_PartC[i];
			if (Part.isboundary == false && Part.computeme == true) {
				if (absinterrupt) {
					densityerroraverage += std::abs(Part.densityerror);
				}
				else {
					densityerroraverage += Part.densityerror;
				}
			}
		}
		return densityerroraverage / totalcomp;
	}

	if (ignorep0tolow) {
//#pragma omp parallel for reduction(+:densityerroraverage)
		for (int i = 0; i < var_MaxParticles; ++i) {
			iisphparticle& Part = var_PartC[i];
			//if (Part.isboundary == false && Part.denstolow == false) {
			if (Part.density >= p0 /* && Part.pos.y < upperviualbord && Part.pos.y > lowervisualbord  && Part.isboundary == false) {
				if (absinterrupt) {
					densityerroraverage += std::abs(Part.densityerror);
				}
				else {
					densityerroraverage += Part.densityerror;
				}
			}
		}
		return (100* (densityerroraverage / numofp0high)/p0);
	}

	if (ignoreincomplete == false && ignorep0tolow == false) {
#pragma omp parallel for reduction(+:densityerroraverage)
		for (int i = 0; i < var_MaxParticles; ++i) {
			iisphparticle& Part = var_PartC[i];
			if (Part.isboundary == false) {
				if (absinterrupt) {
					densityerroraverage += std::abs(Part.densityerror);
				}
				else {
					densityerroraverage += Part.densityerror;
				}
			}
		}
		return densityerroraverage / var_fluidpart;
	}
	*/
	return 0.f; // default return in case all conditions are false
}

