#pragma once
#include <Cell.h>
extern glm::mat3 skewSymmetricMatrix(const glm::vec3& v);
extern glm::mat3 skewSymmetricMatrix2d(const glm::vec3& v);
extern void initrigidbodies(std::vector<iisphparticle>& PartC, std::unordered_map<int, Cell>& hashmap);
extern void updaterigidbody(std::vector<iisphparticle>& PartC);
extern void updaterigidbody2d(std::vector<iisphparticle>& PartC);
extern void makefloatingpartmass(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);
extern void makefloatingpartmass2d(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap);