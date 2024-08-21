#include <scenarios.h>
#include <iisphparticle.h>
#include <part.h>


void var_init(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos = 10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	ParticlesContainer.resize(0);
	var_fluidpart = var_nx * var_ny * var_nz;
	var_spezialboundpart = 0;
	if (singlewall) {
		gammafloat = 0.7;
		var_MaxParticles = var_fluidpart+ (var_oneway * var_oneway *2 + (var_oneway+2)* (var_oneway + 2)*2+ (var_oneway + 1) * (var_oneway ) * 2);
	}
	else {
		gammafloat = 1;
		var_MaxParticles = var_nx * var_ny * var_nz + (var_oneway * var_oneway * 12);
	}

	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	std::srand(std::time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-jitterfac, jitterfac);
	for (int x = 1; x < (var_nx+1); x++) {
		for (int y = 1; y < (var_ny+1); y++) {
			for (int z = 1; z < (var_nz+1); z++) {
				////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				double rando = dis(gen);
				ParticlesContainer[i].pos = glm::vec3((x + rando) * h, (y + rando) * h, (z + rando) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].isfloatingboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	std::cout << i << std::endl;
	int fluid = i;
	//border up and low
	for (int ii = 0; ii < var_oneway; ii++) {
		for (int j = 0; j < var_oneway; j++) {
			//1upper
			ParticlesContainer[i].pos = glm::vec3(ii * h, var_oneway * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			ParticlesContainer[i].pos = glm::vec3(ii * h, 0, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//9lower2
				ParticlesContainer[i].pos = glm::vec3(ii * h, - 1*h, j * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//12upper2
				ParticlesContainer[i].pos = glm::vec3(ii * h, (var_oneway +1) * h, j * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	//border left and right
	for (int ii = -1; ii < (var_oneway+1); ii++) {
		for (int j = -1; j < (var_oneway+1); j++) {
			//3right
			ParticlesContainer[i].pos = glm::vec3((var_oneway ) * h, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			ParticlesContainer[i].pos = glm::vec3(-1*h, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//3right
				ParticlesContainer[i].pos = glm::vec3((var_oneway + 1) * h, ii * h, j * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//6 left
				ParticlesContainer[i].pos = glm::vec3(-2*h, ii * h, j * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	//border front and back
	for (int ii = 0; ii < var_oneway; ii++) {
		for (int j = 0; j < var_oneway+1; j++) {
			//2back
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, -1*h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//5front
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, (var_oneway) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//2back
				ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, -2 * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//5front
				ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, (var_oneway+1)*h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	std::cout << i - fluid << "border" << std::endl;
	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		overallminxpos = std::min(overallminxpos, ParticlesContainer[i].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[i].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[i].pos.z);
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}

}

void var_teslavalveclosed(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos = 10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	ParticlesContainer.resize(0);
	var_spezialboundpart = 0;
	gammafloat = 1.f;
	var_fluidpart = 8 * 9 * 245;
	// fluid + startfloor + (start walls) + big floor + big walls + obstacles
	var_MaxParticles = 8 * 9 * 245 + 48392+ 72+7000+192;
	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	for (int x = -1; x < (8 - 1); x++) {
		for (int y = -1; y < (245 - 1); y++) {
			for (int z = -1; z < (9 - 1); z++) {
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 2) * h, (y + 2) * h, (z + 3) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	std::cout << i << std::endl;
	// back long obstacle
	float goback = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;

		}
		goback -= 0.5;

	}
	// front long obstacle
	float gofront = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;

		}
		gofront += 0.5;

	}
	// back short obstacle
	goback = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//zero1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 20) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//zero2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 21) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;

		}
		goback -= 0.5;

	}
	// front short obstacle
	gofront = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//zeero1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 20) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//teero2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 21) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			ParticlesContainer[i].isobstacle = true;
			i += 1;

		}
		gofront += 0.5;

	}
	//start container floor
	float golower = -1;
	for (int ii = 0; ii < 11; ii++) {
		for (int j = 0; j < 11; j++) {
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 159) * h, golower * h, (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 159) * h, (golower - 1) * h, (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
		golower += 0.1;
	}
	//startcontainer walls right
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 251; j++) {
			//2back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((160) * h, (j + 5) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, ((11)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((9 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, 1 * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((10 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, ((12)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}


	//big container floor and ceiling
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 23; j++) {
			//lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-2) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+4) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//lower 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-3) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+5) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	// big front and back container walls
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 6; j++) {
			//back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-4) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (17) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (18) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	/*
	// big left
	for (int ii = 0; ii < 23; ii++) {
		for (int j = 0; j < 6; j++) {
			//left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((ii)-5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	*/
	//big right
	for (int ii = 0; ii < 12; ii++) {
		for (int j = 0; j < 6; j++) {
			if (ii > 5) {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			else {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}

		}
	}

	//big left
	for (int ii = 0; ii < 12; ii++) {
		for (int j = 0; j < 6; j++) {
			if (ii > 5) {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(9 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			else {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				ParticlesContainer[i].pos = glm::vec3(9 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}

		}
	}
	//start container floor
	golower = 0;
	for (int ii = 0; ii < 9; ii++) {
		for (int j = 0; j < 9; j++) {
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 1) * h, (golower-0)*h, (j + 2) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 2) * h, (golower -1) * h, (j + 2) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
		golower -= 0.1;
	}
	//startcontainer walls left
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 251; j++) {
			//2back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(0, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, ((11)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(9 * h, (j + 5) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, 1 * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(-1 * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, ((12)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}


	std::cout << i << std::endl;

	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		overallminxpos = std::min(overallminxpos, ParticlesContainer[i].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[i].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[i].pos.z);
		if (ParticlesContainer[i].isboundary) {
			//ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
}

void var_teslavalve(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos = 10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	ParticlesContainer.resize(0);
	var_spezialboundpart = 0;
	gammafloat = 1.f;
	var_fluidpart = 8 * 9 * 245;
	// fluid + startfloor + (start walls) + big floor + big walls + obstacles
	var_MaxParticles = 8 * 9 * 245 +48392+72 + 7000 + 192;
	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	for (int x = -1; x < (8 - 1); x++) {
		for (int y = -1; y < (245 - 1); y++) {
			for (int z = -1; z < (9 - 1); z++) {
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 162) * h, (y + 2) * h, (z + 3) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	// back long obstacle
	float goback = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		goback -= 0.5;

	}
	// front long obstacle
	float gofront = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		gofront += 0.5;

	}
	// back short obstacle
	goback = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//zero1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 20) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//zero2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 21) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (3 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		goback -= 0.5;

	}
	// front short obstacle
	gofront = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//zeero1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 20) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//teero2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 21) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (9 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		gofront += 0.5;

	}
	//start container floor
	float golower = -1;
	for (int ii = 0; ii < 11; ii++) {
		for (int j = 0; j < 11; j++) {
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 159) * h, golower * h, (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 159) * h, (golower - 1) * h, (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
		golower += 0.1;
	}
	//startcontainer walls right
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 251; j++) {
			//2back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((160) * h, (j + 5) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, ((11)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((9 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, 1 * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((10 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, ((12)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}


	//big container floor and ceiling
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 23; j++) {
			//lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-2) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+4) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//lower 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-3) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+5) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	// big front and back container walls
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 6; j++) {
			//back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-4) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (17) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (18) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	/*
	// big left
	for (int ii = 0; ii < 23; ii++) {
		for (int j = 0; j < 6; j++) {
			//left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((ii)-5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	*/
	//big right
	for (int ii = 0; ii < 12; ii++) {
		for (int j = 0; j < 6; j++) {
			if (ii > 5) {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			else {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}

		}
	}

	//big left
	for (int ii = 0; ii < 12; ii++) {
		for (int j = 0; j < 6; j++) {
			if (ii > 5) {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(9 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			else {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				ParticlesContainer[i].pos = glm::vec3(9 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}

		}
	}
	//start container floor
	golower = 0;
	for (int ii = 0; ii < 9; ii++) {
		for (int j = 0; j < 9; j++) {
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 1) * h, (golower - 0) * h, (j + 2) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 2) * h, (golower - 1) * h, (j + 2) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
		golower -= 0.1;
	}
	//startcontainer walls left
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 251; j++) {
			//2back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(0, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, ((11)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(9 * h, (j + 5) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, 1 * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(-1 * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, (j - 1) * h, ((12)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}




	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		overallminxpos = std::min(overallminxpos, ParticlesContainer[i].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[i].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[i].pos.z);
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
}


void var_real_teslavalve(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos = 10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	ParticlesContainer.resize(0);
	var_spezialboundpart = 0;
	gammafloat = 1.f;
	var_fluidpart = 8 * 9 * 65 + 5 * 7 * 140;
	// fluid + startfloor + (start walls) + big floor + big walls + obstacles
	var_MaxParticles = 27106 + 8 * 9 * 65 + 5 * 7 * 140;//8 * 9 * 25 + (11 * 11 * 2 + 7 * 11 * 101) + (23 * 150 * 4 + 4 * 150 * 6 + 12 * 6 + 23 * 6) + (12 * 16 * 8 + 12 * 8 * 8)+10000;

	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	//standing water
	for (int x = -1; x < (8 - 1); x++) {
		for (int y = -1; y < (65 - 1); y++) {
			for (int z = -1; z < (9 - 1); z++) {
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 162) * h, (y + 2) * h, (z + 3) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	//water on the ground
	for (int x = 5; x < (145); x++) {
		for (int y = -3; y < 2; y++) {
			for (int z = 0; z < 7; z++) {
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 10) * h, (y + 2) * h, (z + 3) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	// back long obstacle
	float goback = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		goback -= 0.5;

	}
	// front long obstacle
	float gofront = 0;
	for (int ii = 0; ii < 16; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 35) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 36) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 65) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 66) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 95) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 96) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 125) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 126) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		gofront += 0.5;

	}
	// back short obstacle
	goback = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (1 + goback) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		goback -= 0.5;

	}
	// front long obstacle
	gofront = 0;
	for (int ii = 0; ii < 8; ii++) {
		for (int j = 0; j < 6; j++) {
			//first1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 50) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//first2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 51) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 80) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//second2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 81) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 110) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//third2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 111) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth1
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 141) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;
			//forth2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 142) * h, (j - 1) * h, (11 + gofront) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].a = obstaclea;
			i += 1;

		}
		gofront += 0.5;

	}
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 6; j++) {
			if ((ii < 56) || (ii > 85 && ii < 116)) {
				//front
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (12) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;

				//front2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (11) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			if (ii < 26 || (ii > 56 && ii < 86) || ii > 115) {
				//back
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//back2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (0) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;

			}
			//outer
			if ((ii > 67 && ii < 98) || (ii > 126)) {
				//front
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (17) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;

				//front2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (18) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			if (((ii < 68 && ii >36) || (ii > 98 && ii < 130))) {
				//back
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//back2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 11) * h, (j - 2) * h, (-6) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya; // boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;

			}
		}
	}

	//start container floor
	float golower = -1;
	for (int ii = 0; ii < 11; ii++) {
		for (int j = 0; j < 11; j++) {
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 159) * h, golower * h, (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 159) * h, (golower * h - 1), (j + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
		golower += 0.1;
	}
	//startcontainer walls right
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 101; j++) {
			//2back
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((160) * h, (j + 5) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, ((11)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((9 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//2back2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, 1 * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 right
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((10 + 160) * h, (j - 1) * h, (ii + 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			// 5 front2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 160) * h, (j - 1) * h, ((12)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}


	//big container floor and ceiling
	for (int ii = 0; ii < 150; ii++) {
		for (int j = 0; j < 23; j++) {
			//lower
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-2) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+4) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//lower 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (-3) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//upper 2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 10) * h, (+5) * h, (j - 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	// 

	// big left 
	for (int ii = 0; ii < 23; ii++) {
		for (int j = 0; j < 6; j++) {
			//left
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(10 * h, (j - 2) * h, ((ii)-5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			i += 1;
		}
	}
	//big left
	for (int ii = 0; ii < 12; ii++) {
		for (int j = 0; j < 6; j++) {
			if (ii > 5) {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((23 - ii)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
			else {
				//3 left
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3(159 * h, (j - 2) * h, ((ii)-5) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}

		}
	}
	std::cout << i << std::endl;

	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		overallminxpos = std::min(overallminxpos, ParticlesContainer[i].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[i].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[i].pos.z);
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
}

void watercolumn(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos = 10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	// Clear the container
	ParticlesContainer.resize(0);

	// Define boundary and fluid parameters
	int var_boundarypart = 29 * 29 + 25 * 2 * (watercolheight + 10) + 27 * 2 * (watercolheight + 10);
	gammafloat = 0.7;
	if (!singlewall) {
		//gammafloat = 1;
		var_boundarypart = 29 * 29 * 2 + 25 * 2 * (watercolheight + 10) + 27 * 2 * (watercolheight + 10) + 27 * 2 * (watercolheight + 100) + 29 * 2 * (watercolheight + 100);
	}

	// Parameters for the fluid particles
	var_nx = 25;
	var_nz = 25;
	var_fluidpart = var_nx * watercolheight * var_nz;
	var_spezialboundpart = 0;
	if (addfloating) {
		var_spezialboundpart = 1302;
	}
	
	var_MaxParticles = var_fluidpart + var_boundarypart + var_spezialboundpart;
	hashsize = var_MaxParticles;

	// Resize the container to fit all particles
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	std::srand(std::time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-jitterfac, jitterfac);
	// Create fluid particles
	for (int x = 0; x < var_nx; x++) {
		for (int y = 0; y < watercolheight; y++) {
			for (int z = 0; z < var_nz; z++) {
				double rando = dis(gen);
				double randomnum = jitterfac * static_cast<double>(std::rand()) / (RAND_MAX);
				ParticlesContainer[i].pos = glm::vec3((x + rando) * h, (y + rando) * h, (z + rando) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = false;
				ParticlesContainer[i].index = i;
				i++;
			}
		}
	}
	std::cout << "fluid: " << i << std::endl;
	if (addfloating) {
		float j = i;
		std::vector<glm::vec3> duckVertices = parseObjFile("../exportdata/duck2.obj");

		if (duckVertices.empty()) {
			std::cerr << "No vertices found in the file." << std::endl;
		}

		for (const auto& vertex : duckVertices) {
			ParticlesContainer[i].pos = (vertex + glm::vec3(10.f, watercolheight + 10.125f, 10.f) / 20.f) * h * 20.f; // Skalierung
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].isfloatingboundary = true;
			ParticlesContainer[i].index = i;
			i++;
		}
		std::cout << "object: " << i - j << std::endl;
	}
	
	float boden = 0;
	// Create bottom boundary particles
	for (int x = -2; x <= 26; x++) {
		for (int y = -2; y <= -1; y++) {
			for (int z = -2; z <= 26; z++) {
				if (y > -2 || !singlewall) {
					ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (z)*h);
					ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
					ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
					ParticlesContainer[i].a = boundarya;
					ParticlesContainer[i].isboundary = true;
					ParticlesContainer[i].index = i;
					i++;
					boden++;
				}
			}
		}
	}
	std::cout << "boden: " << boden << std::endl;
	// Create front and back single
	float frontandback = 0;
	for (int x = 0; x <= 24; x++) {
		for (int y = 0; y < watercolheight + 10; y++) {
			ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (-1) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			frontandback++;
			i++;
			ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (25) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			frontandback++;
			i++;
		}
	}
	std::cout << " vorder und rück wände single : " << frontandback << std::endl;
	// Create front and back double
	frontandback = 0;
	if (!singlewall) {
		for (int x = -1; x < 26; x++) {
			for (int y = 0; y < watercolheight + 10; y++) {
				ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (-2) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				frontandback++;
				ParticlesContainer[i].pos = glm::vec3((x)*h, y * h, (26) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				frontandback++;
			}
		}
		std::cout << " vorder und rück wände double : " << frontandback << std::endl;
	}
	int leftandright = 0;
	// Create left and right single
	for (int z = -1; z <= 25; z++) {
		for (int y = 0; y < watercolheight + 10; y++) {
			ParticlesContainer[i].pos = glm::vec3((-1) * h, y * h, (z)*h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
			ParticlesContainer[i].pos = glm::vec3((25) * h, y * h, (z)*h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
		}
	}
	std::cout << " left und right wände single : " << leftandright << std::endl;
	// Create left and right double
	leftandright = 0;
	if (!singlewall) {
		for (int z = -2; z <= 26; z++) {
			for (int y = 0; y < watercolheight + 10; y++) {
				ParticlesContainer[i].pos = glm::vec3((-2) * h, y * h, (z)*h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				leftandright++;
				ParticlesContainer[i].pos = glm::vec3((26) * h, y * h, (z)*h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				leftandright++;
			}
		}
		std::cout << " left und right wände double : " << leftandright << std::endl;
	}
	std::cout << i << std::endl;
	// Finalize particle properties
	for (int j = 0; j < var_MaxParticles; j++) {
		ParticlesContainer[j].ismovingboundary = false;
		ParticlesContainer[j].resetvalues();
		overallminxpos = std::min(overallminxpos, ParticlesContainer[j].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[j].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[j].pos.z);
		if (ParticlesContainer[j].isboundary) {
			ParticlesContainer[j].a = boundarya;
		}
		else {
			ParticlesContainer[j].a = 250;
		}
	}
}
void watercolumnsmall(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos = 10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	// Clear the container
	ParticlesContainer.resize(0);

	// Define boundary and fluid parameters
	int var_boundarypart =  19 * 19 + 15 * 2 * (watercolheight + 20) + 17 * 2 * (watercolheight + 20);
	gammafloat = 0.7;
	if (!singlewall) {
		//gammafloat = 1;
		var_boundarypart =  19 * 19 * 2 + 15 * 2 * (watercolheight + 20) + 17 * 2 * (watercolheight + 20) + 17 * 2 * (watercolheight + 20) + 19 * 2 * (watercolheight + 20);
	}

	// Parameters for the fluid particles
	var_nx = 15;
	var_nz = 15;
	var_fluidpart = var_nx * watercolheight * var_nz;
	var_spezialboundpart = 0;
	if (addfloating) {
		var_spezialboundpart += 125;
		var_boundarypart += 10;
	}
	var_MaxParticles = var_fluidpart + var_boundarypart + var_spezialboundpart;
	hashsize = var_MaxParticles;

	// Resize the container to fit all particles
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	std::srand(std::time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-jitterfac, jitterfac);
	const double offset = 10.0;

	// Create fluid particles
	for (int x = 0; x < var_nx; x++) {
		for (int y = 0; y < watercolheight; y++) {
			for (int z = 0; z < var_nz; z++) {
				double rando = dis(gen);
				ParticlesContainer[i].pos = glm::vec3((x + rando + offset) * h, (y + rando + offset) * h, (z + rando + offset) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = false;
				ParticlesContainer[i].index = i;
				i++;
			}
		}
	}
	if (addfloating) {
		for (float xi = 0; xi < 5; xi += 1) {
			for (float yi = 0; yi < 5; yi += 1) {
				for (float zi = 0; zi < 5; zi += 1) {
					ParticlesContainer[i].pos = glm::vec3((xi + 5 + offset) * h, (10 + watercolheight + yi + offset) * h, (zi + 5 + offset) * h);
					ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
					ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
					ParticlesContainer[i].a = boundarya;
					ParticlesContainer[i].isboundary = false;
					ParticlesContainer[i].isfloatingboundary = true;
					ParticlesContainer[i].index = i;
					i++;
				}
			}
		}
		std::cout << "Number of rigidboundary particles: " << i << std::endl;
		for (float xi = 0; xi < 5; xi += 1) {
			for (float yi = 0; yi < 2; yi += 1) {
				ParticlesContainer[i].pos = glm::vec3((xi + 5 + offset) * h, (9 + watercolheight + offset) * h, (yi + 5 + offset) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].isfloatingboundary = false;
				ParticlesContainer[i].index = i;
				i++;
			}
		}
	}
	

	float boden = 0;
	// Create bottom boundary particles
	for (int x = -2; x <= 16; x++) {
		for (int y = -2; y <= -1; y++) {
			for (int z = -2; z <= 16; z++) {
				if (y > -2 || !singlewall) {
					ParticlesContainer[i].pos = glm::vec3((x + offset) * h, (y + offset) * h, (z + offset) * h);
					ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
					ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
					ParticlesContainer[i].a = boundarya;
					ParticlesContainer[i].isboundary = true;
					ParticlesContainer[i].index = i;
					i++;
					boden++;
				}
			}
		}
	}
	std::cout << "boden: " << boden << std::endl;

	// Create front and back single
	float frontandback = 0;
	for (int x = 0; x <= 14; x++) {
		for (int y = 0; y < watercolheight + 20; y++) {
			ParticlesContainer[i].pos = glm::vec3((x + offset) * h, (y + offset) * h, (-1 + offset) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			frontandback++;
			i++;
			ParticlesContainer[i].pos = glm::vec3((x + offset) * h, (y + offset) * h, (15 + offset) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			frontandback++;
			i++;
		}
	}
	std::cout << "vorder und rück wände single: " << frontandback << std::endl;

	// Create front and back double
	frontandback = 0;
	if (!singlewall) {
		for (int x = -1; x < 16; x++) {
			for (int y = 0; y < watercolheight + 20; y++) {
				ParticlesContainer[i].pos = glm::vec3((x + offset) * h, (y + offset) * h, (-2 + offset) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				frontandback++;
				ParticlesContainer[i].pos = glm::vec3((x + offset) * h, (y + offset) * h, (16 + offset) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				frontandback++;
			}
		}
		std::cout << "vorder und rück wände double: " << frontandback << std::endl;
	}

	int leftandright = 0;
	// Create left and right single
	for (int z = -1; z <= 15; z++) {
		for (int y = 0; y < watercolheight + 20; y++) {
			ParticlesContainer[i].pos = glm::vec3((-1 + offset) * h, (y + offset) * h, (z + offset) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
			ParticlesContainer[i].pos = glm::vec3((15 + offset) * h, (y + offset) * h, (z + offset) * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
		}
	}
	std::cout << "left und right wände single: " << leftandright << std::endl;

	// Create left and right double
	leftandright = 0;
	if (!singlewall) {
		for (int z = -2; z <= 16; z++) {
			for (int y = 0; y < watercolheight + 20; y++) {
				ParticlesContainer[i].pos = glm::vec3((-2 + offset) * h, (y + offset) * h, (z + offset) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				leftandright++;
				ParticlesContainer[i].pos = glm::vec3((16 + offset) * h, (y + offset) * h, (z + offset) * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				leftandright++;
			}
		}
		std::cout << "left und right wände double: " << leftandright << std::endl;
	}

	// Finalize particle properties
	for (int j = 0; j < var_MaxParticles; j++) {
		ParticlesContainer[j].ismovingboundary = false;
		ParticlesContainer[j].resetvalues();
		overallminxpos = std::min(overallminxpos, ParticlesContainer[j].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[j].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[j].pos.z);
		if (ParticlesContainer[j].isboundary) {
			ParticlesContainer[j].a = boundarya;
		}
		else {
			ParticlesContainer[j].a = 250;
		}
		if (ParticlesContainer[j].isfloatingboundary) {
			ParticlesContainer[j].a = 200;
			ParticlesContainer[j].r = 200;
			ParticlesContainer[j].b = 0;
			ParticlesContainer[j].g = 0;
		}
	}
}
void moving_boundary(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos = 10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	ParticlesContainer.resize(0);
	var_spezialboundpart = 0/*num_rot_part*/;
	var_nx = var_oneway - 2;
	var_nz = var_oneway - 2;
	var_ny = 20;
	var_fluidpart = var_nx * var_ny * var_nz - std::min(num_rot_part, (var_ny-1) * num_rot_part / 8);
	if (singlewall) {
		gammafloat = 0.7;
		var_MaxParticles = var_fluidpart +num_rot_part+ (var_oneway * var_oneway * 2 + (var_oneway + 2) * (var_oneway + 2) * 2 + (var_oneway + 1) * (var_oneway) * 2);
	}
	else {
		gammafloat = 1;
		var_MaxParticles = var_nx * var_ny * var_nz + (var_oneway * var_oneway * 12) + num_rot_part;
	}

	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	std::srand(std::time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-jitterfac, jitterfac);
	for (int x = -num_rot_part / 8; x < (var_nx - (num_rot_part / 8)); x++) {
		for (int y = 0; y < (var_ny); y++) {
			for (int z = -num_rot_part / 8; z < (var_nz - (num_rot_part / 8)); z++) {
				if (!(((x + 1) >= 0 && (x + 1) < num_rot_part / 8) && ((y + 1) >= 2 && (y + 1) < 10) && (z + 1) == 0)) {
					////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
					double rando = dis(gen);
					ParticlesContainer[i].pos = glm::vec3((x + 1+rando) * (h * 1), (y + 1 + rando) * (h * 1), (z + 1 + rando) * (h * 1));
					ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
					ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
					ParticlesContainer[i].a = 250;
					ParticlesContainer[i].isboundary = 0;
					ParticlesContainer[i].index = i;
					i += 1;
				}

			}
		}
	}
	std::cout << "numfluid" << i << std::endl;
	std::cout << "numfluid diff " << i - ((var_oneway - 2) * (var_oneway - 2) * 10 - num_rot_part) << std::endl;

	int nfluid = i;
	//moving stuff
	for (int ii = 0; ii < num_rot_part / 8; ++ii) {
		for (int jj = 0; jj < 8; ++jj) {
			glm::vec3 velocity(deltaT / 2, deltaT / 2, 0.0f); // Geschwindigkeit in die korrekte Richtung
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].vel = velocity;
			ParticlesContainer[i].pos = glm::vec3((0 + ii) * h, (2 + jj) * h, 0 * h);
			ParticlesContainer[i].ismovingboundary = true;
			ParticlesContainer[i].index = i;
			i += 1;
		}

	}
	std::cout << "nummoving" << i - nfluid << std::endl;
	int nmbound = i - nfluid;
	//border up and low
	for (int ii = 0; ii < var_oneway; ii++) {
		for (int j = 0; j < var_oneway; j++) {
			//1upper
			ParticlesContainer[i].pos = glm::vec3((ii - (num_rot_part / 8)) * h, var_oneway * h, (j - (num_rot_part / 8)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//4lower
			ParticlesContainer[i].pos = glm::vec3((ii - (num_rot_part / 8)) * h, 0, (j - (num_rot_part / 8)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//9lower2
				ParticlesContainer[i].pos = glm::vec3((ii - (num_rot_part / 8)) * h, -1 * h, (j - (num_rot_part / 8)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//12upper2
				ParticlesContainer[i].pos = glm::vec3((ii - (num_rot_part / 8)) * h, (var_oneway + 1) * h, (j - (num_rot_part / 8)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	//border left and right
	for (int ii = -1; ii < (var_oneway + 1); ii++) {
		for (int j = -1; j < (var_oneway + 1); j++) {
			//3right
			ParticlesContainer[i].pos = glm::vec3((var_oneway - (num_rot_part / 8))*h, ii * h, (j - (num_rot_part / 8)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			ParticlesContainer[i].pos = glm::vec3(( - 1 - (num_rot_part / 8)) * h, ii * h, (j - (num_rot_part / 8)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//3right
				ParticlesContainer[i].pos = glm::vec3((var_oneway + 1 - (num_rot_part / 8)) * h, ii * h, (j - (num_rot_part / 8)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//6 left
				ParticlesContainer[i].pos = glm::vec3(( - 2 - (num_rot_part / 8)) * h, ii * h, (j - (num_rot_part / 8)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	//border front and back
	for (int ii = 0; ii < var_oneway; ii++) {
		for (int j = 0; j < var_oneway + 1; j++) {
			//2back
			ParticlesContainer[i].pos = glm::vec3((ii - (num_rot_part / 8)) * h, j * h, ( - 1 - (num_rot_part / 8)) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//5front
			ParticlesContainer[i].pos = glm::vec3((ii - (num_rot_part / 8)) * h, j * h, (var_oneway - (num_rot_part / 8))*h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//2back
				ParticlesContainer[i].pos = glm::vec3((ii - (num_rot_part / 8)) * h, j * h, ( - 2 - (num_rot_part / 8)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//5front
				ParticlesContainer[i].pos = glm::vec3((ii - (num_rot_part / 8)) * h, j * h, (var_oneway + 1 - (num_rot_part / 8)) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	std::cout << "numbound" << i - nfluid - nmbound << std::endl;
	int nbound = i - nfluid - nmbound;
	int totalp = i;
	std::cout << "numdiff" << var_MaxParticles - nfluid - nmbound - nbound << std::endl;

	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		overallminxpos = std::min(overallminxpos, ParticlesContainer[i].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[i].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[i].pos.z);
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
		if (ParticlesContainer[i].ismovingboundary) {
			ParticlesContainer[i].a = obstaclea;
		}
	}
}
void dambreaktest(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos = 10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	ParticlesContainer.resize(0);
	var_spezialboundpart = 0;
	var_nx = 54;
	var_ny = 30;
	var_nz = 68;
	var_fluidpart = var_nx * var_ny * var_nz;
	if (singlewall) {
		gammafloat = 0.7;
		var_MaxParticles = var_nx * var_ny * var_nz + (180 * 70 + 2 * 180 * 60 + 2 * 70 * 60);
	}
	else {
		gammafloat = 1;
		var_MaxParticles = var_nx * var_ny * var_nz + 2 * (180 * 70 + 2 * 180 * 60 + 2 * 70 * 60);
	}

	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	for (int x = 0; x < (var_nx); x++) {
		for (int y = 0; y < (var_ny); y++) {
			for (int z = 0; z < (var_nz); z++) {
				////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 1) * (h * 1), (y + 1) * (h * 1), (z + 1) * (h * 1));
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	//border back and front
	for (int ii = 0; ii < 180; ii++) {
		for (int j = 0; j < 60; j++) {
			//2back
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//5front
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, (70 - 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//7back2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (j - 1) * h, (0 - 1)*h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//10front2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (j - 1) * h, (70 - 0) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	//border left and right 
	for (int ii = 0; ii < 60; ii++) {
		for (int j = 0; j < 70; j++) {
			//3right
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((180 - 1) * h, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(0, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//8right2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((180 - 2) * h, (ii - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//11 left2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((0 - 1)*h, (ii - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	// border bottom
	for (int ii = 0; ii < 180; ii++) {
		for (int j = 0; j < 70; j++) {

			//4lower
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, 0, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//9lower2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii)*h, (0 - 1) * h, (j)*h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		overallminxpos = std::min(overallminxpos, ParticlesContainer[i].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[i].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[i].pos.z);
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
}

void smalldambreaktest(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos = 10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	ParticlesContainer.resize(0);
	var_spezialboundpart = 0;
	var_nx = 54;
	var_ny = 30;
	var_nz = 33;
	var_fluidpart = var_nx * var_ny * var_nz;
	if (singlewall) {
		gammafloat = 1;
		var_MaxParticles = var_nx * var_ny * var_nz + (180 * 35 + 2 * 180 * 60 + 2 * 35 * 60) + (2 * 10 * 20 + 10 * 20 + 2 * 10 * 10);
	}
	else {
		gammafloat = 1;
		var_MaxParticles = var_nx * var_ny * var_nz + 2 * (180 * 35 + 2 * 180 * 60 + 2 * 35 * 60) + (2 * 10 * 20 + 10 * 20 + 2 * 10 * 10);
	}

	hashsize = var_MaxParticles;
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	for (int x = 0; x < (var_nx); x++) {
		for (int y = 0; y < (var_ny); y++) {
			for (int z = 0; z < (var_nz); z++) {
				////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((x + 2) * (h * 1), (y + 1) * (h * 1), (z + 1) * (h * 1));
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = 250;
				ParticlesContainer[i].isboundary = 0;
				ParticlesContainer[i].index = i;
				i += 1;
			}
		}
	}
	//border back and front
	for (int ii = 0; ii < 180; ii++) {
		for (int j = 0; j < 60; j++) {
			//2back
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, 0);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//5front
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, j * h, (35 - 1) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//7back2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (j - 1) * h, (0 - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//10front2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii - 1) * h, (j - 1) * h, (35 - 0) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	//border left and right 
	for (int ii = 0; ii < 60; ii++) {
		for (int j = 0; j < 35; j++) {
			//3right
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((180 - 1) * h, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//6 left
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(1 * h, ii * h, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//8right2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((180 - 2) * h, (ii - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//11 left2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((0), (ii - 1) * h, (j - 1) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	// border bottom
	for (int ii = 0; ii < 180; ii++) {
		for (int j = 0; j < 35; j++) {

			//4lower
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3(ii * h, 0, j * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (singlewall == false) {
				//9lower2
				//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii)*h, (0 - 1) * h, (j)*h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	// box
	for (int ii = 0; ii < 10; ii++) {
		for (int j = 0; j < 20; j++) {
			//upper
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((ii + 120) * h, (10 + 1) * h, (j + 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//8right2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((130 - 1) * h, (ii + 1) * h, (j + 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			//11 left2
			//ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
			ParticlesContainer[i].pos = glm::vec3((120) * h, (ii + 1) * h, (j + 5) * h);
			ParticlesContainer[i].index = i;
			ParticlesContainer[i].isboundary = true;
			i += 1;
			if (j < 10) {
				//2back
			////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 120) * h, (j + 1) * h, 5 * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
				//5front
				////ParticlesContainer[i].cameradistance = glm::length2(ParticlesContainer[i].pos - CameraPosition);
				ParticlesContainer[i].pos = glm::vec3((ii + 120) * h, (j + 1) * h, (24) * h);
				ParticlesContainer[i].index = i;
				ParticlesContainer[i].isboundary = true;
				i += 1;
			}
		}
	}
	for (int i = 0; i < var_MaxParticles; i++) {
		ParticlesContainer[i].resetvalues();
		ParticlesContainer[i].ismovingboundary = false;
		ParticlesContainer[i].isfloatingboundary = false;
		overallminxpos = std::min(overallminxpos, ParticlesContainer[i].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[i].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[i].pos.z);
		if (ParticlesContainer[i].isboundary) {
			ParticlesContainer[i].a = boundarya;
		}
		else
		{
			ParticlesContainer[i].a = 250;
		}
	}
}

//_________________________________________________________________________________________________________________________________________________________2d

void watercolumnTwoD(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos = 10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	// Clear the container
	ParticlesContainer.resize(0);

	// Define boundary and fluid parameters
	int var_boundarypart = 19 + 2 * (watercolheight + 100);
	gammafloat = 0.7;
	if (!singlewall) {
		//gammafloat = 1;
		var_boundarypart = 2 * 19 + 2 * 2 * (watercolheight + 100) + 40;
	}

	// Parameters for the fluid particles
	var_spezialboundpart = 0;
	if (addfloating) {
		var_spezialboundpart = 36;
	}

	var_nx = 15;
	var_nz = 1;
	var_fluidpart = var_nx * watercolheight * var_nz;
	var_MaxParticles = var_fluidpart + var_boundarypart + var_spezialboundpart;
	hashsize = var_MaxParticles;

	// Resize the container to fit all particles
	ParticlesContainer.resize(var_MaxParticles);
	for (int l = 0; l < ParticlesContainer.size(); l++) {
		ParticlesContainer[l].isfloatingboundary = false;
	}
	int i = 0;
	std::srand(std::time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-jitterfac, jitterfac);
	const double offset = 5.0;

	// Create fluid particles
	for (int x = 0; x < var_nx; x++) {
		for (int y = 0; y < watercolheight; y++) {
			double rando = dis(gen);
			ParticlesContainer[i].pos = glm::vec3((x + rando + offset) * h, (y + rando + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = 250;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].index = i;
			i++;
		}
	}
	if (addfloating) {
		for (float xi = 0; xi < 4.5; xi = xi + 0.5) {
			ParticlesContainer[i].pos = glm::vec3((xi + 2 + offset) * h, (watercolheight + 1 + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].isfloatingboundary = true;
			ParticlesContainer[i].index = i;
			i++;

			ParticlesContainer[i].pos = glm::vec3((0 + 2 + offset) * h, (watercolheight + 1 + xi + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].isfloatingboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			ParticlesContainer[i].pos = glm::vec3((4 + 2 + offset) * h, (watercolheight + 1 + xi + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].isfloatingboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			ParticlesContainer[i].pos = glm::vec3((xi + 2 + offset) * h, (watercolheight + 1 + 4 + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].isfloatingboundary = true;
			ParticlesContainer[i].index = i;
			i++;
		}
	}

	std::cout << "fluid: " << i << std::endl;
	float boden = 0;
	// Create bottom boundary particles
	for (int x = -2; x <= 16; x++) {
		for (int y = -2; y <= -1; y++) {
			if (y > -2 || !singlewall) {
				ParticlesContainer[i].pos = glm::vec3((x + offset) * h, (y + offset) * h, 0 * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				boden++;
			}
		}
	}
	std::cout << "boden: " << boden << std::endl;
	int leftandright = 0;
	// Create left and right single

	for (int y = 0; y < watercolheight + 100; y++) {
		ParticlesContainer[i].pos = glm::vec3((-1 + offset) * h, (y + offset) * h, 0 * h);
		ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
		ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
		ParticlesContainer[i].a = boundarya;
		ParticlesContainer[i].isboundary = true;
		ParticlesContainer[i].index = i;
		i++;
		leftandright++;
		ParticlesContainer[i].pos = glm::vec3((15 + offset) * h, (y + offset) * h, 0 * h);
		ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
		ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
		ParticlesContainer[i].a = boundarya;
		ParticlesContainer[i].isboundary = true;
		ParticlesContainer[i].index = i;
		i++;
		leftandright++;
	}

	std::cout << " left und right wände single : " << leftandright << std::endl;
	// Create left and right double
	leftandright = 0;
	if (!singlewall) {
		for (int y = 0; y < watercolheight + 100; y++) {
			ParticlesContainer[i].pos = glm::vec3((-2 + offset) * h, (y + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
			ParticlesContainer[i].pos = glm::vec3((16 + offset) * h, (y + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
		}
		std::cout << " left und right wände double : " << leftandright << std::endl;
	}

	// Finalize particle properties
	for (int j = 0; j < var_MaxParticles; j++) {
		ParticlesContainer[j].ismovingboundary = false;
		ParticlesContainer[j].resetvalues();
		ParticlesContainer[j].m = (p0 * h * h);
		overallminxpos = std::min(overallminxpos, ParticlesContainer[j].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[j].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[j].pos.z);
		if (ParticlesContainer[j].isboundary) {
			ParticlesContainer[j].a = boundarya;
		}
		else {
			ParticlesContainer[j].a = 250;
		}
		if (ParticlesContainer[j].isfloatingboundary) {
			ParticlesContainer[j].a = 200;


			ParticlesContainer[j].r = 200;
			ParticlesContainer[j].b = 0;
			ParticlesContainer[j].g = 0;
		}
	}
}
void DambreaktestTwoD(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer) {
	overallminxpos =10000;
	overallminypos = 10000;
	overallminzpos = 10000;
	// Clear the container
	ParticlesContainer.resize(0);

	// Define boundary and fluid parameters
	int var_boundarypart = 2 * 150 + 2 * (100);
	if (!singlewall) {
		//gammafloat = 1;
		var_boundarypart = 2 * 2 * 150 + 2 * 2 * (100);
	}

	// Parameters for the fluid particles
	var_spezialboundpart = 0;
	if (addfloating) {
		var_spezialboundpart = 36;
	}

	var_nz = 1;
	var_fluidpart = var_nx * var_ny * 1;
	var_MaxParticles = var_fluidpart + var_boundarypart + var_spezialboundpart;
	hashsize = var_MaxParticles;

	// Resize the container to fit all particles
	ParticlesContainer.resize(var_MaxParticles);
	int i = 0;
	std::srand(std::time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-jitterfac, jitterfac);
	const double offset = 10.0;

	// Create fluid particles
	for (int x = 0; x < var_nx; x++) {
		for (int y = 0; y < var_ny; y++) {
			double rando = dis(gen);
			ParticlesContainer[i].pos = glm::vec3((x + rando + offset) * h, (y + rando + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = 250;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].index = i;
			i++;
		}
	}
	if (addfloating) {
		for (float xi = 0; xi < 4.5; xi = xi + 0.5) {
			ParticlesContainer[i].pos = glm::vec3((xi + 2 + var_nx + offset) * h, (1 + var_ny + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].isfloatingboundary = true;
			ParticlesContainer[i].index = i;
			i++;

			ParticlesContainer[i].pos = glm::vec3((0 + 2 + var_nx + offset) * h, (1 + xi + var_ny + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].isfloatingboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			ParticlesContainer[i].pos = glm::vec3((4 + 2 + var_nx + offset) * h, (1 + xi + var_ny + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].isfloatingboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			ParticlesContainer[i].pos = glm::vec3((xi + 2 + var_nx + offset) * h, (1 + 4 + var_ny + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = false;
			ParticlesContainer[i].isfloatingboundary = true;
			ParticlesContainer[i].index = i;
			i++;
		}
	}

	std::cout << "fluid: " << i << std::endl;
	float boden = 0;
	// Create bottom and top boundary particles
	for (int x = -2; x <= 148; x++) {
		for (int y = -2; y <= -1; y++) {
			if (y > -2 || !singlewall) {
				ParticlesContainer[i].pos = glm::vec3((x + offset) * h, (y + offset) * h, 0 * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				boden++;
				ParticlesContainer[i].pos = glm::vec3((x + offset) * h, (100 + y + offset) * h, 0 * h);
				ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
				ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
				ParticlesContainer[i].a = boundarya;
				ParticlesContainer[i].isboundary = true;
				ParticlesContainer[i].index = i;
				i++;
				boden++;
			}
		}
	}
	std::cout << "boden: " << boden << std::endl;
	int leftandright = 0;
	// Create left and right single

	for (int y = 0; y < 100; y++) {
		ParticlesContainer[i].pos = glm::vec3((-1 + offset) * h, (y + offset) * h, 0 * h);
		ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
		ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
		ParticlesContainer[i].a = boundarya;
		ParticlesContainer[i].isboundary = true;
		ParticlesContainer[i].index = i;
		i++;
		leftandright++;
		ParticlesContainer[i].pos = glm::vec3((147 + offset) * h, (y + offset) * h, 0 * h);
		ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
		ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
		ParticlesContainer[i].a = boundarya;
		ParticlesContainer[i].isboundary = true;
		ParticlesContainer[i].index = i;
		i++;
		leftandright++;
	}

	std::cout << " left und right wände single : " << leftandright << std::endl;
	// Create left and right double
	leftandright = 0;
	if (!singlewall) {
		for (int y = 0; y < 100; y++) {
			ParticlesContainer[i].pos = glm::vec3((-2 + offset) * h, (y + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
			ParticlesContainer[i].pos = glm::vec3((148 + offset) * h, (y + offset) * h, 0 * h);
			ParticlesContainer[i].vel = glm::vec3(0, 0, 0);
			ParticlesContainer[i].acc = glm::vec3(0, 0, 0);
			ParticlesContainer[i].a = boundarya;
			ParticlesContainer[i].isboundary = true;
			ParticlesContainer[i].index = i;
			i++;
			leftandright++;
		}
		std::cout << " left und right wände double : " << leftandright << std::endl;
	}

	// Finalize particle properties
	for (int j = 0; j < var_MaxParticles; j++) {
		ParticlesContainer[j].ismovingboundary = false;
		ParticlesContainer[j].resetvalues();
		ParticlesContainer[j].m = (p0 * h * h);
		overallminxpos = std::min(overallminxpos, ParticlesContainer[j].pos.x);
		overallminypos = std::min(overallminypos, ParticlesContainer[j].pos.y);
		overallminzpos = std::min(overallminzpos, ParticlesContainer[j].pos.z);
		if (ParticlesContainer[j].isboundary) {
			ParticlesContainer[j].a = boundarya;
		}
		else {
			ParticlesContainer[j].a = 250;
		}
		if (ParticlesContainer[j].isfloatingboundary) {
			ParticlesContainer[j].a = 200;
			ParticlesContainer[j].r = 200;
			ParticlesContainer[j].b = 0;
			ParticlesContainer[j].g = 0;
		}
	}
}
