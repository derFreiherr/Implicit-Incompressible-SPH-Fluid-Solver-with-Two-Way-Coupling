#pragma once
#include <Cell.h>
//3d
extern void var_init(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void var_teslavalveclosed(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void var_teslavalve(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void var_real_teslavalve(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void watercolumn(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void watercolumnsmall(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void moving_boundary(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void dambreaktest(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void smalldambreaktest(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
//2d
extern void watercolumnTwoD(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);
extern void DambreaktestTwoD(glm::vec3 CameraPosition, std::vector<iisphparticle>& ParticlesContainer);