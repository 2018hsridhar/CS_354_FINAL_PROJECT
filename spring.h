#ifndef SPRING_H
#define SPRING_H 

#include <vector>
#include <stdlib.h>
#include <string>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_transform.hpp>

struct DIST
{
	int id;
	float dist;
};

struct Mass
{
	int m_id;
	float m;	
	glm::vec3 old_pos;
	glm::vec3 curr_pos;
	glm::vec3 vel;
	glm::vec3 force;
	std::vector<int> neighbors;
	std::vector<int> springs;

	/*! Constructs a mass with a position, set of neighboring masses, and set of springs. */
	Mass(int new_id, glm::vec3 _pos);

	/*! Methods for mass struct */
	void zero_out_forces(); 						// reset all force values to 0
	void applyForce(glm::vec3(force));
	void updateVel(float time_diff);
	void updatePos(float time_diff);
	int getID();

};

class Spring
{
	int s_id;
	Mass *A;
	Mass *B;
	float kCoeff = 1.0; // find some number online, or just test our yourself
	float equil_len;
	float disp; // (-) displacement = compressed, (+) displacement = stretched

public:
	/*! Constructs a spring, with 2 masses. */
	Spring(int s_id, Mass *a, Mass *b);
	Mass* getMassB();
	/*! Methods for spring class */
	float calc_SpringForce();

};

#endif
