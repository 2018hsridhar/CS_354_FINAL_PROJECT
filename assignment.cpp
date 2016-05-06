#include <limits>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>
#include <spring.h>

// OpenGL library includes
#include <GL/glew.h>
#include <GLFW/glfw3.h>

float curTime = 0.0;
float timeStep = 0.0009f;  // analogy to 60 fps
float floor_coeff = 7.5;
float floor_cutoff = -1.60;
glm::vec3 gravity = glm::vec3(0,-0.245,0);
//glm::vec3 gravity = glm::vec3(0,-9.8,0);
bool max_id_floorHit = false;

int max_x_mass_id = 0; //surface vertic
int max_y_mass_id = 0; //
int max_z_mass_id = 0; //
std::string file_name = "obj/sphere.obj";
std::string file_name_base = "obj/sphere_tiny.obj";

bool force_applied = false;

std::vector<Spring> mass_springs;

std::ostream& operator<<(std::ostream& os, const glm::vec2& v) {
  os << glm::to_string(v);
  return os;
}

std::ostream& operator<<(std::ostream& os, const glm::vec3& v) {
  os << glm::to_string(v);
  return os;
}

std::ostream& operator<<(std::ostream& os, const glm::vec4& v) {
  os << glm::to_string(v);
  return os;
}

int window_width = 800, window_height = 600;
const std::string window_title = "Menger";
int max_id = 0; // for vertex that initial force is applied to

// VBO and VAO descriptors.

// We have these VBOs available for each VAO.
// note :: add plane VBOs here
enum {
  kVertexBuffer,
  kIndexBuffer,
  // plane_kVertexBuffer,
  // plane_kIndexBuffer,
  // small_sphere_kVertexBuffer,
  // small_sphere_kIndexBuffer,
  kNumVbos,
};

// These are our VAOs.
enum {
  kMengerVao,
  kPlaneVao,
  kSmallVao,
  kNumVaos
};

GLuint array_objects[kNumVaos];  // This will store the VAO descriptors.
GLuint buffer_objects[kNumVaos][kNumVbos];  // These will store VBO descriptors.

float last_x = 0.0f, last_y = 0.0f, current_x = 0.0f, current_y = 0.0f;
bool drag_state = false;
int current_button = -1;
float camera_distance = 3.0;
float pan_speed = 0.1f;
float roll_speed = 0.1f;
float rotation_speed = 0.05f;
float zoom_speed = 0.1f;
bool fps_mode = false;

////////// Useful globals ////////////////////
glm::vec3 eye (0,0,camera_distance);
glm::vec3 look (0,0,-1);
glm::vec3 up (0,1,0);
int L = 0; // global variable (helps with easy manipulation)
float aspect = 0.0f;
std::vector<glm::vec4> menger_vertices;
std::vector<glm::uvec3> menger_faces;

///////////////////////////////////////////////////////
/*! for generating the small sphers to help visualize deformations */
std::vector<glm::vec4> base_sphere_vertices;
std::vector<glm::uvec3> base_sphere_faces;
glm::vec4 base_center = glm::vec4(0,0.005,0,1); // homogenous coordinates ( world matrix system ) - needs to be small sphere center

std::vector<glm::vec4> small_sphere_vertices;
std::vector<glm::uvec3> small_sphere_faces;

///////////////////////////////////////////////////////

std::vector<Mass*> masses;
std::vector<int> visited_Masses;
glm::vec3 init_Force = glm::vec3(0,5,0);

std::vector<glm::vec4> plane_vertices;
std::vector<glm::uvec3> plane_faces;
////////////////////////////////////////////

const char* vertex_shader =
    "#version 330 core\n"
    "uniform vec4 light_position;"
    "in vec4 vertex_position;"
    "out vec4 vs_light_direction;"
    "void main() {"
    "gl_Position = vertex_position;"
    "vs_light_direction = light_position - gl_Position;"
    "}";

// already have ccess to data here
const char* geometry_shader =
    "#version 330 core\n"
    "layout (triangles) in;"
    "layout (triangle_strip, max_vertices = 3) out;"
    "uniform mat4 projection;"
    "uniform mat4 view;"
    "in vec4 vs_light_direction[];"
    "out vec4 normal;"
    "out vec4 light_direction;"
    "out vec4 world_position;"
    "void main() {"
    "int n = 0;"
    "vec3 a = gl_in[0].gl_Position.xyz;"
    "vec3 b = gl_in[1].gl_Position.xyz;"
    "vec3 c = gl_in[2].gl_Position.xyz;"
    "vec3 u = normalize(b - a);"
    "vec3 v = normalize(c - a);"
    "normal = normalize(vec4(normalize(cross(u, v)), 0.0));"
    "for (n = 0; n < gl_in.length(); n++) {"
    "light_direction = normalize(vs_light_direction[n]);"
    "world_position = gl_in[n].gl_Position;"					// passed to the plane fragment shader
    "gl_Position = projection * view * gl_in[n].gl_Position;" // at this point gl_Position data is changed
    "EmitVertex();"
    "}"
    "EndPrimitive();"
    "}";

const char* fragment_shader =
    "#version 330 core\n"
    "in vec4 normal;"
    "in vec4 light_direction;"
    "out vec4 fragment_color;"
    "void main() {"
    "vec4 color = normal;"
    "float dot_nl = dot(normalize(light_direction), normalize(normal));"
    "dot_nl = clamp(dot_nl, 0.0, 1.0);"
    "fragment_color = clamp(dot_nl * color, 0.0, 1.0);"
    "}";

const char* small_spheres_fragment_shader =
    "#version 330 core\n"
    "in vec4 normal;"
    "in vec4 light_direction;"
    "out vec4 small_spheres_fragment_color;"
    "void main() {"
    "small_spheres_fragment_color = vec4(0,1,0,1);"
    "}";


// how to pass unteransformed world coordinates from
// geoemtry shader to plane_fragment shader

// make suyre to use the right coordinates
// use basic color (i.e.. solid white ) to check the plane
// multi sampling?
const char* plane_fragment_shader =
    "#version 330 core\n"
    "in vec4 world_position;"
    "in vec4 light_direction;"
    "in vec4 normal;"
    "out vec4 plane_fragment_color;"
    "float checkSize = 1;"
    "void main() {"
    "float fModResult = mod(floor(checkSize * world_position.x) + floor(checkSize * world_position.z), 2.0);"
    "if(fModResult < 1.0) {"
	"plane_fragment_color = vec4(1,1,1,1);"
    "} else { "
	"plane_fragment_color = vec4(0,0,0,1);"
    "}"
    "float dot_nl = dot(normalize(light_direction), normalize(normal));"
    "dot_nl = clamp(dot_nl, 0.0, 1.0);"
    "plane_fragment_color = clamp(dot_nl * plane_fragment_color, 0.0, 1.0);"
    "}";

const char* OpenGlErrorToString(GLenum error) {
  switch (error) {
    case GL_NO_ERROR:
      return "GL_NO_ERROR";
      break;
    case GL_INVALID_ENUM:
      return "GL_INVALID_ENUM";
      break;
    case GL_INVALID_VALUE:
      return "GL_INVALID_VALUE";
      break;
    case GL_INVALID_OPERATION:
      return "GL_INVALID_OPERATION";
      break;
    case GL_OUT_OF_MEMORY:
      return "GL_OUT_OF_MEMORY";
      break;
    default:
      return "Unknown Error";
      break;
  }
  return "Unicorns Exist";
}

#define CHECK_SUCCESS(x) \
  if (!(x)) {            \
    glfwTerminate();     \
    exit(EXIT_FAILURE);  \
  }

#define CHECK_GL_SHADER_ERROR(id)                                           \
  {                                                                         \
    GLint status = 0;                                                       \
    GLint length = 0;                                                       \
    glGetShaderiv(id, GL_COMPILE_STATUS, &status);                          \
    glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);                         \
    if (!status) {                                                          \
      std::string log(length, 0);                                           \
      glGetShaderInfoLog(id, length, nullptr, &log[0]);                     \
      std::cerr << "Line :" << __LINE__ << " OpenGL Shader Error: Log = \n" \
                << &log[0];                                                 \
      glfwTerminate();                                                      \
      exit(EXIT_FAILURE);                                                   \
    }                                                                       \
  }

#define CHECK_GL_PROGRAM_ERROR(id)                                           \
  {                                                                          \
    GLint status = 0;                                                        \
    GLint length = 0;                                                        \
    glGetProgramiv(id, GL_LINK_STATUS, &status);                             \
    glGetProgramiv(id, GL_INFO_LOG_LENGTH, &length);                         \
    if (!status) {                                                           \
      std::string log(length, 0);                                            \
      glGetProgramInfoLog(id, length, nullptr, &log[0]);                     \
      std::cerr << "Line :" << __LINE__ << " OpenGL Program Error: Log = \n" \
                << &log[0];                                                  \
      glfwTerminate();                                                       \
      exit(EXIT_FAILURE);                                                    \
    }                                                                        \
  }

#define CHECK_GL_ERROR(statement)                                             \
  {                                                                           \
    { statement; }                                                            \
    GLenum error = GL_NO_ERROR;                                               \
    if ((error = glGetError()) != GL_NO_ERROR) {                              \
      std::cerr << "Line :" << __LINE__ << " OpenGL Error: code  = " << error \
                << " description =  " << OpenGlErrorToString(error);          \
      glfwTerminate();                                                        \
      exit(EXIT_FAILURE);                                                     \
    }                                                                         \
  }

// nbote :: some of these were supposed to be negative
glm::mat4 Perspective(float fovy, float aspect, float near, float far) {
	float elem_1_1 = 1.0 / (aspect * tan(0.5f * fovy));
	float elem_2_2 = 1.0 / tan(0.5f * fovy);
	float elem_3_3 = -1 * (near+far)/(far-near);
	float elem_3_4 = -1 * (2*near*far)/(far-near);
 
  return (glm::transpose(glm::mat4(glm::vec4(elem_1_1, 0, 0, 0),
                   glm::vec4(0, elem_2_2, 0,0),
                   glm::vec4(0,0,elem_3_3, elem_3_4),
                   glm::vec4(0,0,-1,0))));
}

// eye, center, up
/* NOTE :: GO TO VOUGA AND GET THIS DEBUGGED !!!! */
/* note :: Vouga's slide was not fully accurate */
glm::mat4 LookAt(const glm::vec3& eye, const glm::vec3& at, const glm::vec3& up) {

	glm::vec3 look = glm::normalize(at - eye);
	glm::vec3 tan  = glm::cross(look,up); 

	return glm::transpose(glm::mat4(
		glm::vec4(tan.x, tan.y, tan.z, -1 * glm::dot(tan, eye)),
		glm::vec4(up.x, up.y, up.z, -1 * glm::dot(up, eye)),
		glm::vec4(-1 * look.x,-1 * look.y, -1 * look.z, glm::dot(look, eye)),
        	glm::vec4(0.000000, 0.000000, 0.000000, 1.000000)
	));

}


/* Function to draw the plane */
void drawPlane()
{
	// plane above - normal above (wouldn't see anything)
	// try to debug using a finite square
	plane_vertices.push_back(glm::vec4(0,-2,0,1)); // point at the center ( need w = 1 to define plane at y = -2) 
	plane_vertices.push_back(glm::vec4(100,-2,100,1)); // four other vertices for triangels covering quadrants 
	plane_vertices.push_back(glm::vec4(-100,-2,100,1));
	plane_vertices.push_back(glm::vec4(-100,-2,-100,1));
	plane_vertices.push_back(glm::vec4(100,-2,-100,1));

	// may need to debug these faces
    plane_faces.push_back(glm::uvec3(0,3,2));
    plane_faces.push_back(glm::uvec3(0,4,3));
    plane_faces.push_back(glm::uvec3(0,1,4));
    plane_faces.push_back(glm::uvec3(0,2,1));

}

void CreatePlane(){
	drawPlane();
}

	void LoadObj(const std::string& file, std::vector<glm::vec4>& vertices, std::vector<glm::uvec3>& indices) {
		std::ifstream in(file);
		int i = 0, j = 0;
		glm::vec4 vertex = glm::vec4(0.0, 0.0, 0.0, 1.0);
		glm::uvec3 face_indices = glm::uvec3(0, 0, 0);
		while (in.good()) {
			char c = in.get();
			switch (c) {
			case 'v':
				in >> vertex[0] >> vertex[1] >> vertex[2];
				vertices.push_back(vertex);

          ////////////////////////////
        	// Masses and vertices have a one-to-one mapping
			// if we are working with small spheres, do not push their masses
			if(file.compare(file_name) == 0)
				masses.push_back(new Mass(i,glm::vec3(vertices[i])));
          ////////////////////////////
            i++;
				break;
			case 'f':
				in >> face_indices[0] >> face_indices[1] >> face_indices[2];
				face_indices -= 1;
				indices.push_back(face_indices);
				break;
			default:
				break;
			}
		}
		in.close();
	}

	//////////////////////////////////////
	std::vector<int>  getNeighbors(int x)
	{
		
		// 1. iterate over all vertices
		std::vector<int> closest_vertices;
		std::vector<DIST> distances;
		for(int i = 0; i< masses.size(); i++)
		{
			if(i != x) {
				float distance = glm::distance(masses[x]->curr_pos, masses[i]->curr_pos);
				struct DIST new_dist_obj;
				new_dist_obj.dist = distance;
				new_dist_obj.id = masses[i]->m_id;
				distances.push_back(new_dist_obj);
		}
		}

		// 2. calculate distances
		std::sort(distances.begin(),distances.end(), 
		[](DIST x, DIST y)->bool
			{
				return x.dist < y.dist; 
			}
		);

			//3. get 6 closest vertices
		for(int j = 0; j < 6; j++) {
		closest_vertices.push_back(distances[j].id);
		}
			return closest_vertices;
	}

	// need to keep a set of visited masses, based on their ids
	// timestep t = 0.017
	// make sure this is not incorrect

	/* PROBABLY BUGGY */
	void calc_NetForces(int mass_id)
	{
			if(curTime <= 0.010 && force_applied==false)
		{
			if(!max_id_floorHit) {
				masses[max_id]->applyForce(init_Force);
			} else {
				std::cout << "max id floor hit is now true" << std::endl;
			}
			force_applied = true;
		}
				

			//calculate spring forces, add them
			std::vector<int> my_springs = masses[mass_id]->springs;
			for(int i = 0; i < my_springs.size(); i++) 
			{
	
				// get spring masses and vectors for masses
				Mass *self = mass_springs[my_springs[i]].getMassA();
				Mass *other = mass_springs[my_springs[i]].getMassB();

				glm::vec3 v_1 = self->old_pos - other->old_pos;
				glm::vec3 v_2 = self->curr_pos - other->curr_pos;

				// calculate how far spring has been displaced 
				float disp = glm::length(v_2) - glm::length(v_1);
				mass_springs[i].setDisplacement(disp);
				glm::vec3 spring_force = glm::vec3(1,1,1);
				// std::cout << "disp = " << disp << std::endl;

				// calculate spring force, based on displacement (q_{i+1}) 
				if(disp >= 0) 
				{
					glm::vec3 force_dir = glm::normalize(other->curr_pos - self->curr_pos);
					spring_force = force_dir * mass_springs[i].calc_SpringForce();
				}
				else 
				{
					glm::vec3 force_dir = glm::normalize(self->curr_pos - other->curr_pos);
					spring_force = force_dir * mass_springs[i].calc_SpringForce();
				}

				glm::vec3 dampening = mass_springs[i].calc_Dampening(masses[mass_id]->vel);
				//std::cout << "spring force = " << to_string(spring_force) << std::endl;
				//std::cout << "damepning force = " << to_string(dampening) << std::endl;

				// apply spring forces ( to A and B)
				self->applyForce(spring_force);
				//self->applyForce(dampening);
				other->applyForce(-1.0f * spring_force);

			}

			

		}

		void Load_SpringSystem()
		{

			// vertex list has a direct relationship to the map list ( one-to-one)
				for(int i = 0; i < masses.size(); i++) {
					masses[i]->neighbors = getNeighbors(i);
				}
		
			// at frame ( time t = 0), apply one force in (x,y,z) direction to 
			// surface element mass whose (x,y,z) postion is largest positively 
			float max_x =  std::numeric_limits<float>::min(); //surface vertic
			float max_y =  std::numeric_limits<float>::min(); //
			float max_z =  std::numeric_limits<float>::min(); //

			for(int i = 0; i < masses.size(); i++) {
				if(masses[i]->curr_pos.x > max_x){ 
					max_x_mass_id = i;
				}
				if(masses[i]->curr_pos.y > max_y){ 
					max_y_mass_id = i;
				}
				if(masses[i]->curr_pos.z > max_z){ 
					max_z_mass_id = i;
				}
			}
				max_id = max_y_mass_id;
		}
		
	void load_smallSpheres()
	{
			for(int i = 0; i < menger_vertices.size(); i++)
			{
				// calculate vector, from base_center, to menger_vertex 
				glm::vec4 offset_vector = menger_vertices[i] - base_center;
				//std::cout << to_string(offset_vector) << std::endl;
				for(int j = 0; j < base_sphere_vertices.size();j++) 
				{
					glm::vec4 shifted_vertex = 	base_sphere_vertices[j] + offset_vector;
					small_sphere_vertices.push_back(shifted_vertex);
				}
			}

			//for(int i = 0; i < menger_faces.size();i++) 
			for(int i = 0; i < menger_vertices.size();i++) 
			{
				// becuase our faces are of form (v_i,v_j,v_k) 
				// sphere 1 :: (v_i,v_j,v_k)
				// sphere 2 :: (v_i + 62,v_j + 62,v_k + 62)
				int face_offset = 62 * i; 
				for(int j = 0; j < base_sphere_faces.size();j++)  // base_sphere_faces has 120 faces 
				{
					glm::vec3 init_faceValues = base_sphere_faces[j];
					glm::vec3 new_faceValues = init_faceValues + glm::vec3(face_offset,face_offset,face_offset);
					small_sphere_faces.push_back(new_faceValues);
				}
			}
	}

	void setupSpring()
	{
			int index =0;
			for(int x =0; x<masses.size(); x++)
			{
				std::vector<int> spring;
				for(int y = 0; y<masses[x]->neighbors.size(); ++y)
				{  
					int num = masses[x]->neighbors[y];
					Spring s(index, masses[x], masses[num]);
					spring.push_back(index);
					mass_springs.push_back(s);
					index++;
				}
				masses[x]->springs = spring;
			}

			// for(int x =0; x<masses.size(); x++)
			// {
			//     std::cout<<masses[x]->m_id<<std::endl;
				
			//     for(int y = 0; y<masses[x]->springs.size(); ++y)
			//     {
			//      std::cout<<"A    "<<(mass_springs[x*6+y].getMassA())->m_id<<std::endl;           
			//        std::cout<<"B  "<<(mass_springs[x*6+y].getMassB())->m_id<<std::endl; 
			//     }
			//     std::cout<<"\n"<<std::endl;
			// }
		}



		void ErrorCallback(int error, const char* description) {
		std::cerr << "GLFW Error: " << description << "\n";
		}

		void KeyCallback(GLFWwindow* window, int key, int scancode, int action,
				int mods) {

		if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
			glfwSetWindowShouldClose(window, GL_TRUE);
		else if (key == GLFW_KEY_W && action != GLFW_RELEASE) {

				glm::vec3 center_old = eye + (camera_distance * look);

				camera_distance -= zoom_speed;
				glm::vec3 center = eye + (camera_distance * look);

				look = glm::normalize(center - eye);
				eye += (center_old - center);

	} else if (key == GLFW_KEY_S && action != GLFW_RELEASE) {

				glm::vec3 center_old = eye + (camera_distance * look);

				camera_distance += zoom_speed;
				glm::vec3 center = eye + (camera_distance * look);

				look = glm::normalize(center - eye);
				eye += (center_old - center);

		} else if (key == GLFW_KEY_A && action != GLFW_RELEASE) {

				glm::vec3 tangent = glm::cross(look,up);

				glm::vec3 center = eye + (camera_distance * look);
				glm::vec3 new_center = center - (pan_speed * tangent);			

				glm::vec3 new_eye = eye + (new_center - center);
				eye = new_eye;


		} else if (key == GLFW_KEY_D && action != GLFW_RELEASE) {

				glm::vec3 tangent = glm::cross(look,up);

				glm::vec3 center = eye + (camera_distance * look);
				glm::vec3 new_center = center + (pan_speed * tangent);			

				glm::vec3 new_eye = eye + (new_center - center);
				eye = new_eye;

		} else if (key == GLFW_KEY_LEFT && action != GLFW_RELEASE) {

			glm::vec3 up_rotated = glm::rotate(up, roll_speed, look);
			up = up_rotated;

		} else if (key == GLFW_KEY_RIGHT && action != GLFW_RELEASE) {

			glm::vec3 up_rotated = glm::rotate(up, -1 * roll_speed, look);
			up = up_rotated;

		} else if (key == GLFW_KEY_DOWN && action != GLFW_RELEASE) {

			glm::vec3 tangent = glm::cross(look,up);

			glm::vec3 center = eye + (camera_distance * look);
			glm::vec3 new_center = center - (pan_speed * up);			

			glm::vec3 new_eye = eye + (new_center - center);
			eye = new_eye;

		} else if (key == GLFW_KEY_UP && action != GLFW_RELEASE) {

				glm::vec3 tangent = glm::cross(look,up);

				glm::vec3 center = eye + (camera_distance * look);
				glm::vec3 new_center = center + (pan_speed * up);			

				glm::vec3 new_eye = eye + (new_center - center);
				eye = new_eye;

		} else if (key == GLFW_KEY_C && action != GLFW_RELEASE)
			fps_mode = !fps_mode;
	}

	void MousePosCallback(GLFWwindow* window, double mouse_x, double mouse_y) {
	last_x = current_x;
	last_y = current_y;
	current_x = mouse_x;
	current_y = mouse_y;
	float delta_x = current_x - last_x;
	float delta_y = current_y - last_y;
	if (sqrt(delta_x * delta_x + delta_y * delta_y) < 1e-15) return;
	if (drag_state && current_button == GLFW_MOUSE_BUTTON_LEFT) {
		glm::vec3 mouse_direction = glm::normalize(glm::vec3(delta_x, delta_y, 0.0f));
		glm::vec3 tan = glm::cross(look,up);
		glm::vec3 md_world3 = (mouse_direction.x * tan) - (mouse_direction.y * up); // really it is just a change of bases
		glm::vec3 rotation_axis = glm::normalize(glm::cross(look, md_world3));	

		if(!fps_mode){ // orbital mode

			// note :: your eye stays in the same spot
			// but your look and up vectors ROTATE
			// note :: orbital is ROTATE AROUND OBJECT, not REVOLVE THE EYE
			glm::vec3 up_rotation = glm::rotate(up, rotation_speed, rotation_axis);	
			glm::vec3 look_rotation = glm::rotate(look, rotation_speed, rotation_axis);	
			glm::vec3 eye_rotation = glm::rotate(eye, rotation_speed, rotation_axis);	

			eye = eye_rotation;
			up = up_rotation;
			look = look_rotation;
		} else { // fps mode ( note :: if look changes up changes too)

			glm::vec3 look_rotation = glm::rotate(look, rotation_speed, rotation_axis);	
			look = look_rotation;

			glm::vec3 up_rotation = glm::rotate(up, rotation_speed, rotation_axis);	
			up = up_rotation;
		}

	} else if (drag_state && current_button == GLFW_MOUSE_BUTTON_RIGHT) {
			if(delta_y > 0) {

				glm::vec3 center_old = eye + (camera_distance * look);

				camera_distance -= zoom_speed;
				glm::vec3 center = eye + (camera_distance * look);

				look = glm::normalize(center - eye);
				eye += (center_old - center);
			} else {
				glm::vec3 center_old = eye + (camera_distance * look);

				camera_distance += zoom_speed;
				glm::vec3 center = eye + (camera_distance * look);

				look = glm::normalize(center - eye);
				eye += (center_old - center);
			}
	}
	}

	void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
	drag_state = (action == GLFW_PRESS);
	current_button = button;
	}

	int main(int argc, char* argv[]) {

	// SETTING UP OpenGL Context
	if (!glfwInit()) exit(EXIT_FAILURE);

	glfwSetErrorCallback(ErrorCallback);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SAMPLES, 4);

	GLFWwindow* window = glfwCreateWindow(window_width, window_height,
											&window_title[0], nullptr, nullptr);
	CHECK_SUCCESS(window != nullptr);
	glfwMakeContextCurrent(window);
	glewExperimental = GL_TRUE;
	CHECK_SUCCESS(glewInit() == GLEW_OK);
	glGetError();  // clear GLEW's error for it

	glfwSetKeyCallback(window, KeyCallback);
	glfwSetCursorPosCallback(window, MousePosCallback);
	glfwSetMouseButtonCallback(window, MouseButtonCallback);
	glfwSwapInterval(1);
	const GLubyte* renderer = glGetString(GL_RENDERER);  // get renderer string
	const GLubyte* version = glGetString(GL_VERSION);    // version as a string
	std::cout << "Renderer: " << renderer << "\n";
	std::cout << "OpenGL version supported:" << version << "\n";



		// load up vertices of mesh, and of the small spheres
	LoadObj(file_name, menger_vertices, menger_faces);
	LoadObj(file_name_base, base_sphere_vertices,base_sphere_faces);

	// Create spring systems, and small spheres
	Load_SpringSystem();
	setupSpring();
	load_smallSpheres();
	
	// Plane
	CreatePlane();
	std::cout << "Loaded plane and vertices geometries" << std::endl; 

	// Setup our VAO array
	CHECK_GL_ERROR(glGenVertexArrays(kNumVaos, array_objects));

	/*****************************************
	 *****************************************
	 *****************************************
	 *****************************************
	 *****************************************
	 */

	// Setup the menger array object.
	// Switch to the kMenger VAO.
	CHECK_GL_ERROR(glBindVertexArray(array_objects[kMengerVao]));

	// Generate buffer objects for kMengerVao
	CHECK_GL_ERROR(glGenBuffers(kNumVbos, &buffer_objects[kMengerVao][0]));

	// Setup vertex data for kMenger VBOs
	CHECK_GL_ERROR(
		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kVertexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
								sizeof(float) * menger_vertices.size() * 4,
								&menger_vertices[0], GL_STATIC_DRAW));
	CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
	CHECK_GL_ERROR(glEnableVertexAttribArray(0));

	// Setup element array buffer. (kMenger faces data )
	CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
								buffer_objects[kMengerVao][kIndexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
								sizeof(uint32_t) * menger_faces.size() * 3,
								&menger_faces[0], GL_STATIC_DRAW));

	// Let's create our Menger SHADER program.
	GLuint program_id = 0;
	CHECK_GL_ERROR(program_id = glCreateProgram());

	// Setup vertex shader.
	GLuint vertex_shader_id = 0;
	const char* vertex_source_pointer = vertex_shader;
	CHECK_GL_ERROR(vertex_shader_id = glCreateShader(GL_VERTEX_SHADER));
	CHECK_GL_ERROR(
		glShaderSource(vertex_shader_id, 1, &vertex_source_pointer, nullptr));
	glCompileShader(vertex_shader_id);
	CHECK_GL_SHADER_ERROR(vertex_shader_id);

	// Setup geometry shader.
	GLuint geometry_shader_id = 0;
	const char* geometry_source_pointer = geometry_shader;
	CHECK_GL_ERROR(geometry_shader_id = glCreateShader(GL_GEOMETRY_SHADER));
	CHECK_GL_ERROR(
		glShaderSource(geometry_shader_id, 1, &geometry_source_pointer, nullptr));
	glCompileShader(geometry_shader_id);
	CHECK_GL_SHADER_ERROR(geometry_shader_id);

	// Setup fragment shader.
	GLuint fragment_shader_id = 0;
	const char* fragment_source_pointer = fragment_shader;
	CHECK_GL_ERROR(fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
	CHECK_GL_ERROR(
		glShaderSource(fragment_shader_id, 1, &fragment_source_pointer, nullptr));
	glCompileShader(fragment_shader_id);
	CHECK_GL_SHADER_ERROR(fragment_shader_id);

		// ATTACH SHADERS
	CHECK_GL_ERROR(glAttachShader(program_id, vertex_shader_id));
	CHECK_GL_ERROR(glAttachShader(program_id, fragment_shader_id));
	CHECK_GL_ERROR(glAttachShader(program_id, geometry_shader_id));

	// Bind attributes. ( part of linking step )
	CHECK_GL_ERROR(glBindAttribLocation(program_id, 0, "vertex_position"));
	CHECK_GL_ERROR(glBindFragDataLocation(program_id, 0, "fragment_color"));
	glLinkProgram(program_id);
	CHECK_GL_PROGRAM_ERROR(program_id);

	// Get the uniform locations.
	GLint projection_matrix_location = 0;
	CHECK_GL_ERROR(projection_matrix_location =
						glGetUniformLocation(program_id, "projection"));
	GLint view_matrix_location = 0;
	CHECK_GL_ERROR(view_matrix_location =
						glGetUniformLocation(program_id, "view"));
	GLint light_position_location = 0;
	CHECK_GL_ERROR(light_position_location =
						glGetUniformLocation(program_id, "light_position"));

	///////////////////////////////////////////////
	///////////////////////////////////////////////
	///////////////////////////////////////////////
	//small sphere set up  

	CHECK_GL_ERROR(glBindVertexArray(array_objects[kSmallVao]));
	// Generate buffer objects for kPlaneVao 
	CHECK_GL_ERROR(glGenBuffers(kNumVbos, &buffer_objects[kSmallVao][0]));

	// Let's create our plane SHADER program.
	GLuint small_sphere_program_id = 2;
	CHECK_GL_ERROR(small_sphere_program_id = glCreateProgram());

	// Setup vertex data for kPlane VBOs
	CHECK_GL_ERROR(
		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kSmallVao][kVertexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
								sizeof(float) * small_sphere_vertices.size() * 4,
								&small_sphere_vertices[0], GL_STATIC_DRAW));
	CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
	CHECK_GL_ERROR(glEnableVertexAttribArray(0));

	// Setup element array buffer. (kPlane faces data )
	CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
								buffer_objects[kSmallVao][kIndexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
								sizeof(uint32_t) * small_sphere_faces.size() * 3,
								&small_sphere_faces[0], GL_STATIC_DRAW));

	// Setup fragment shader.
	GLuint small_sphere_fragment_shader_id = 2;
	const char* small_sphere_fragment_source_pointer = small_spheres_fragment_shader; 
	CHECK_GL_ERROR(small_sphere_fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
	CHECK_GL_ERROR(
		glShaderSource(small_sphere_fragment_shader_id, 1, &small_sphere_fragment_source_pointer, nullptr));
	glCompileShader(small_sphere_fragment_shader_id);
	CHECK_GL_SHADER_ERROR(small_sphere_fragment_shader_id);

	// ATTACH SHADERS
	CHECK_GL_ERROR(glAttachShader(small_sphere_program_id, vertex_shader_id));
	CHECK_GL_ERROR(glAttachShader(small_sphere_program_id, geometry_shader_id));
	CHECK_GL_ERROR(glAttachShader(small_sphere_program_id, small_sphere_fragment_shader_id));

	// Bind attributes. ( linking step )
	CHECK_GL_ERROR(glBindAttribLocation(small_sphere_program_id, 0, "vertex_position"));
	CHECK_GL_ERROR(glBindFragDataLocation(small_sphere_program_id, 0, "plane_fragment_color")); 
		
	glLinkProgram(small_sphere_program_id);
	CHECK_GL_PROGRAM_ERROR(small_sphere_program_id);

	// Get the uniform locations. [ not sure if this also needs to be copied too ]
	GLint small_sphere_projection_matrix_location = 2;
	CHECK_GL_ERROR(small_sphere_projection_matrix_location =
						glGetUniformLocation(small_sphere_program_id, "projection"));
	GLint small_sphere_view_matrix_location = 2;
	CHECK_GL_ERROR(small_sphere_view_matrix_location =
						glGetUniformLocation(small_sphere_program_id, "view"));
	GLint small_sphere_light_position_location = 2;
	CHECK_GL_ERROR(small_sphere_light_position_location =
						glGetUniformLocation(small_sphere_program_id, "light_position"));
	

	///////////////////////////////////////////////
	///////////////////////////////////////////////
	///////////////////////////////////////////////
	// Setup the plane array object.
	// Switch to the kPlaneVAO.
	CHECK_GL_ERROR(glBindVertexArray(array_objects[kPlaneVao]));
	// Generate buffer objects for kPlaneVao 
	CHECK_GL_ERROR(glGenBuffers(kNumVbos, &buffer_objects[kPlaneVao][0]));

	// Let's create our plane SHADER program.
	GLuint plane_program_id = 1;
	CHECK_GL_ERROR(plane_program_id = glCreateProgram());

	// Setup vertex data for kPlane VBOs
	CHECK_GL_ERROR(
		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kPlaneVao][kVertexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
								sizeof(float) * plane_vertices.size() * 4,
								&plane_vertices[0], GL_STATIC_DRAW));
	CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
	CHECK_GL_ERROR(glEnableVertexAttribArray(0));

	// Setup element array buffer. (kPlane faces data )
	CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
								buffer_objects[kPlaneVao][kIndexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
								sizeof(uint32_t) * plane_faces.size() * 3,
								&plane_faces[0], GL_STATIC_DRAW));

	// Setup fragment shader.
	GLuint plane_fragment_shader_id = 1;
	const char* plane_fragment_source_pointer = plane_fragment_shader; 
	CHECK_GL_ERROR(plane_fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
	CHECK_GL_ERROR(
		glShaderSource(plane_fragment_shader_id, 1, &plane_fragment_source_pointer, nullptr));
	glCompileShader(plane_fragment_shader_id);
	CHECK_GL_SHADER_ERROR(plane_fragment_shader_id);

	// ATTACH SHADERS
	CHECK_GL_ERROR(glAttachShader(plane_program_id, vertex_shader_id));
	CHECK_GL_ERROR(glAttachShader(plane_program_id, geometry_shader_id));
	CHECK_GL_ERROR(glAttachShader(plane_program_id, plane_fragment_shader_id));

	// Bind attributes. ( linking step )
	CHECK_GL_ERROR(glBindAttribLocation(plane_program_id, 0, "vertex_position"));
	CHECK_GL_ERROR(glBindFragDataLocation(plane_program_id, 0, "plane_fragment_color")); 
		
	glLinkProgram(plane_program_id);
	CHECK_GL_PROGRAM_ERROR(plane_program_id);

	// Get the uniform locations. [ not sure if this also needs to be copied too ]
	GLint plane_projection_matrix_location = 1;
	CHECK_GL_ERROR(plane_projection_matrix_location =
						glGetUniformLocation(plane_program_id, "projection"));
	GLint plane_view_matrix_location = 1;
	CHECK_GL_ERROR(plane_view_matrix_location =
						glGetUniformLocation(plane_program_id, "view"));
	GLint plane_light_position_location = 1;
	CHECK_GL_ERROR(plane_light_position_location =
						glGetUniformLocation(plane_program_id, "light_position"));

	///////////////////// SET UP SMALL SPHERES SHADER PROGRAMS ////////////////////

	// Setup the small spheres array object.
	// Switch to the kSmallVao.
	// CHECK_GL_ERROR(glBindVertexArray(array_objects[kSmallVao]));
	// Generate buffer objects for kPlaneVao 
	// CHECK_GL_ERROR(glGenBuffers(kNumVbos, &buffer_objects[kSmallVao][0]));

	// // Let's create our plane SHADER program.
	// GLuint small_spheres_program_id = 1;
	// CHECK_GL_ERROR(small_spheres_program_id = glCreateProgram());

	// //send data to VBO for small spheres 	

	// CHECK_GL_ERROR(
	// 	glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kSmallVao][kVertexBuffer]));
	// CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
	// 			sizeof(float) * small_sphere_vertices.size() * 4,
	// 			&small_sphere_vertices[0], GL_STATIC_DRAW));
	// std::cout << "Here One, lin [910], before glEnableVertexAttribArray(0) for small_sphere_VAO" << std::endl;
	// CHECK_GL_ERROR(glEnableVertexAttribArray(0)); // -- why does this cause a SEG FAULT in the program later??

	// // Setup element array buffer. (kMenger faces data )
	// CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
	// 							buffer_objects[kSmallVao][kIndexBuffer]));
	// CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
	// 							sizeof(uint32_t) * small_sphere_faces.size() * 3,
	// 							&small_sphere_faces[0], GL_STATIC_DRAW));

	
	// //Setup fragment shader.
	// GLuint small_spheres_fragment_shader_id = 1;
	// const char* small_spheres_fragment_source_pointer = small_spheres_fragment_shader; 
	// CHECK_GL_ERROR(small_spheres_fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
	// CHECK_GL_ERROR(
	// 	glShaderSource(small_spheres_fragment_shader_id, 1, &small_spheres_fragment_source_pointer, nullptr));
	// glCompileShader(small_spheres_fragment_shader_id);
	// CHECK_GL_SHADER_ERROR(small_spheres_fragment_shader_id);

	// // ATTACH SHADERS
	// CHECK_GL_ERROR(glAttachShader(small_spheres_program_id, vertex_shader_id)); // going to change vertex shader!! or fragment!!
	// CHECK_GL_ERROR(glAttachShader(small_spheres_program_id, geometry_shader_id));
	// CHECK_GL_ERROR(glAttachShader(small_spheres_program_id, small_spheres_fragment_shader_id));

	// // Bind attributes. ( linking step )
	// CHECK_GL_ERROR(glBindAttribLocation(small_spheres_program_id, 0, "vertex_position"));
	// CHECK_GL_ERROR(glBindFragDataLocation(small_spheres_program_id, 0, "small_shaders_fragment_color")); 
		
	// glLinkProgram(small_spheres_program_id);
	// CHECK_GL_PROGRAM_ERROR(small_spheres_program_id);

	// // Get the uniform locations. [ not sure if this also needs to be copied too ]
	// GLint small_spheres_projection_matrix_location = 1;
	// CHECK_GL_ERROR(small_spheres_projection_matrix_location =
	// 					glGetUniformLocation(small_spheres_program_id, "projection"));
	// GLint small_spheres_view_matrix_location = 1;
	// CHECK_GL_ERROR(small_spheres_view_matrix_location =
	// 					glGetUniformLocation(small_spheres_program_id, "view"));
	// GLint small_spheres_light_position_location = 1;
	// CHECK_GL_ERROR(small_spheres_light_position_location =
	// 					glGetUniformLocation(small_spheres_program_id, "light_position"));
	
		//////////////////////////////////////////////
		//////////////////////////////////////////////
		//////////////////////////////////////////////

	glfwSwapInterval(1);
	glm::vec4 light_position = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);

	// rendering loop - swich between VAO and VBOs

	while (!glfwWindowShouldClose(window)) {

		// Setup some basic window stuff.
		glfwGetFramebufferSize(window, &window_width, &window_height);
		glViewport(0, 0, window_width, window_height);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_MULTISAMPLE);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDepthFunc(GL_LESS);



		//////////////////////////////////////////////
		//////////////////////////////////////////////
		//////////////////////////////////////////////

		glm::vec3 center = glm::vec3(eye.x + camera_distance * look.x, eye.y + camera_distance * look.y, eye.z + camera_distance * look.z);
		glm::mat4 view_matrix = LookAt(eye, center, up); 
		aspect = static_cast<float>(window_width) / window_height;
		glm::mat4 projection_matrix = Perspective(45.0f, aspect, 0.0001f, 1000.0f);

		for(int k = 0; k < 20; k++)
		{

			///////////////////////////////////////////////////
			///////////////////////////////////////////////////
			/* SIMULATION CODE 
			* [1] set all NET forces to 0
			* 	[1.1] set positive force to largest surface element = glm::length();
			* [2] calculate q_i for each m_i
			* [3] calculate and apply all external forces to each mass
			* [4] calculate v_i for each m_i
			* [5] reset menger vertices to new, physically simulated, vertices
			*/

			// [1] set all NET forces to 0
			// [2] calculate q_i for each m_i
				for(int j = 0; j < masses.size(); j++) {
					masses[j]->zero_out_forces();
					//std::cout << "Mass [ " << j << "] has position = " << to_string(masses[j]->curr_pos) << std::endl;
					masses[j]->updatePos(timeStep*0.6f);
					//std::cout << "Mass [ " << j << "] update POS has position = " << to_string(masses[j]->curr_pos) << std::endl;
				}

			// [3] calculate and apply ALL external forces to each mass
			for(int j = 0; j < masses.size(); j++) {
				calc_NetForces(masses[j]->m_id);
			}

			if(masses[max_id]->curr_pos.y < floor_cutoff)
				max_id_floorHit = true;
			for(int j = 0; j < masses.size(); j++) {
				masses[j]->applyForce(gravity);
				if(masses[j]->curr_pos.y < floor_cutoff)
				{
					glm::vec3 mass_floor_dist = masses[j]->curr_pos - glm::vec3(0,floor_cutoff,0) ;
					float floor_dist = glm::length(mass_floor_dist);
					masses[j]->applyForce(floor_coeff * glm::vec3(0,floor_dist,0));
				}
			}

			// [4] calculate v_i for each m_i
			for(int j = 0; j < masses.size(); j++) {
				masses[j]->updateVel(timeStep * 0.6f);
				//std::cout << "Mass [ " << j << "] has velocity = " << masses[j]->vel << std::endl;
			}

		// [5] reset menger vertices to new, physically simulated, vertices
			menger_vertices.clear();
			for(int j = 0; j < masses.size(); j++) {
				menger_vertices.push_back(glm::vec4(masses[j]->curr_pos,1));
			} 

			// std::vector<glm::vec4> small_sphere_vertices;
			// std::vector<glm::uvec3> small_sphere_faces;
			// base_sphere_vertices = the small sphere that helps start the process
			// base_sphere_faces);
			// base_center ( vector of the center ) = the center of said small sphere
			// kSmall vao

			// given list of vertices, calculate distances and then set offset of small sphere
			small_sphere_vertices.clear();
			small_sphere_faces.clear();
			load_smallSpheres();
		}

		/*! INSERT OPENGL CODE HERE TO REPASS VERTICES */
		//each bind buffer call is separete 
  		// CHECK_GL_ERROR(glBindVertexArray(array_objects[kMengerVao]));

  		// CHECK_GL_ERROR(
   	//    		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kVertexBuffer]));
	   //  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
				// 	sizeof(float) * menger_vertices.size() * 4,
				// 	&menger_vertices[0], GL_STATIC_DRAW));

	   //  CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
				// 					sizeof(uint32_t) * menger_faces.size() * 3,
				// 					&menger_faces[0], GL_STATIC_DRAW));

		//bind to VAO for small spheres
		CHECK_GL_ERROR(glBindVertexArray(array_objects[kSmallVao]));

  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kSmallVao][kVertexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
					sizeof(float) * small_sphere_vertices.size() * 4,
					&small_sphere_vertices[0], GL_STATIC_DRAW));
	
		// Setup element array buffer. (kMenger faces data )
		CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
									buffer_objects[kSmallVao][kIndexBuffer]));
		CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
									sizeof(uint32_t) * small_sphere_faces.size() * 3,
									&small_sphere_faces[0], GL_STATIC_DRAW));

	
	///////////////////////////////////////////////////
	///////////////////////////////////////////////////

    //Switch to our Menger VAO.
    // CHECK_GL_ERROR(glBindVertexArray(array_objects[kMengerVao]));

    // // Use our program.
    // CHECK_GL_ERROR(glUseProgram(program_id));

    // // Pass uniforms in.
    // CHECK_GL_ERROR(glUniformMatrix4fv(projection_matrix_location, 1, GL_FALSE,
    //                                   &projection_matrix[0][0]));
    // CHECK_GL_ERROR(glUniformMatrix4fv(view_matrix_location, 1, GL_FALSE,
    //                                   &view_matrix[0][0]));
    // CHECK_GL_ERROR(
    //     glUniform4fv(light_position_location, 1, &light_position[0]));

    // // Draw our triangles.
    // CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES, menger_faces.size() * 3,
    //                               GL_UNSIGNED_INT, 0));

    /////////////////////////////////////////////////// 
    /////////////////////////////////////////////////// 
    // Switch to our plane Menger VAO.
    CHECK_GL_ERROR(glBindVertexArray(array_objects[kPlaneVao]));

    // Use our program.
    CHECK_GL_ERROR(glUseProgram(plane_program_id));

    // Pass uniforms in for the plane.  
    CHECK_GL_ERROR(glUniformMatrix4fv(plane_projection_matrix_location, 1, GL_FALSE,
                                      &projection_matrix[0][0]));
    CHECK_GL_ERROR(glUniformMatrix4fv(plane_view_matrix_location, 1, GL_FALSE,
                                      &view_matrix[0][0]));
    CHECK_GL_ERROR(
        glUniform4fv(plane_light_position_location, 1, &light_position[0]));

    // Draw plane's triangles.
    CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES, plane_faces.size() * 3,
                                  GL_UNSIGNED_INT, 0));

    /////////////////////////////////////////////////// 
    /////////////////////////////////////////////////// 
	// Switch to our small spheres VAO.
    CHECK_GL_ERROR(glBindVertexArray(array_objects[kSmallVao]));

    // Use our program.
    CHECK_GL_ERROR(glUseProgram(small_sphere_program_id));

    // Pass uniforms in.
    CHECK_GL_ERROR(glUniformMatrix4fv(small_sphere_projection_matrix_location, 1, GL_FALSE,
                                      &projection_matrix[0][0]));
    CHECK_GL_ERROR(glUniformMatrix4fv(small_sphere_view_matrix_location, 1, GL_FALSE,
                                      &view_matrix[0][0]));
    CHECK_GL_ERROR(
        glUniform4fv(small_sphere_light_position_location, 2, &light_position[0]));

    // Draw our triangles. ( THIS IS SEG FAULTING!!!)
    CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES, small_sphere_faces.size() * 3,
                                  GL_UNSIGNED_INT, 0));
	
	//////////////////////////////////////////////////// 
	// update the current time , based on time step
	curTime += timeStep;
 
    // Poll and swap.
    glfwPollEvents();
    glfwSwapBuffers(window);
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
