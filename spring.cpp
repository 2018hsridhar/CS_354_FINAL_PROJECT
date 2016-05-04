#include "spring.h"
#include <fstream>
#include <iostream>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>

Mass::Mass(int new_id, glm::vec3 _pos)
{
	m_id = new_id;
	curr_pos = _pos;
	old_pos = curr_pos;
}

void Mass::zero_out_forces()
{
	force = glm::vec3(0,0,0);
}

void Mass::applyForce(glm::vec3(_newForce))
{
	force += _newForce;
}

/*! check out equations for this again */
void Mass::updateVel(float time_diff)
{
	vel += force * (1.0f / m) * time_diff;
}

void Mass::updatePos(float time_diff)
{
	glm::vec3 acc = force * (1.0f / m);
	old_pos = curr_pos;
	curr_pos = curr_pos + (vel * time_diff);
}

int Mass::getID()
{
	return m_id;
}

Spring::Spring(int new_id, Mass *_A, Mass *_B)
{
	A = _A;
	B = _B;
	s_id = new_id;	
	disp = 0; 
	equil_len = glm::distance(_A->curr_pos, _B->curr_pos);
	
}

Mass* Spring::getMassB()
{
	return B;
}

// this returns the magnitude of the spring force!!!
float Spring::calc_SpringForce()
{
	return kCoeff * disp;
}



