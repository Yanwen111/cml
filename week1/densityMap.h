#pragma once

#include <vector>

#include <glm/glm.hpp>

// Class that stores the density readings
// and other related info
class DensityMap {
public:
	int dim;
	int radius;

	// 3D array that stores the data
	std::vector<std::vector<std::vector<float>>> cells;

	// Constructor
	DensityMap(int dim);

	// Adds a line of data between p1 and p2
	void addLine(glm::vec3 p1, glm::vec3 p2, std::vector<float> vals);

	// Overwrites everything with zeroes
	void clear();

	// Returns the vertices in a form useful to OpenGL
	std::vector<float> getVertices();
};

// Not being used right now, but maybe in the future
// to get smoother shading
float pointLineDistance(glm::vec3 a, glm::vec3 b, glm::vec3 v);