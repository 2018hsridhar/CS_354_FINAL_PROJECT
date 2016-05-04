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

// Menger info.
const int kMinLevel = 0;
const int kMaxLevel = 4;
const glm::vec3 kMengerMinBounds = glm::vec3(-0.5f);
const glm::vec3 kMengerMaxBounds = glm::vec3(0.5f);

// VBO and VAO descriptors.

// We have these VBOs available for each VAO.
// note :: add plane VBOs here
enum {
  kVertexBuffer,
  kIndexBuffer,
  plane_kVertexBuffer,
  plane_kIndexBuffer,
  kNumVbos,
};

// These are our VAOs.
enum {
  kMengerVao,
  kPlaneVao,
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

std::vector<Mass*> masses;

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

/* Function drawCube(vertices,faces, minx-z, maxx-z) */
void drawCube(std::vector<glm::vec4>& vertices,
		    std::vector<glm::uvec3>& faces,
		    float minx, float miny, float minz,
		    float maxx, float maxy, float maxz ) {

    /* ADD VERTICES TO SPONGE LISTS */
    vertices.push_back(glm::vec4(minx, miny, maxz,1.0f)); // A = 0 
    vertices.push_back(glm::vec4(maxx, miny, maxz,1.0f)); // B = 1 

    vertices.push_back(glm::vec4(minx, maxy, maxz,1.0f)); // C = 2 
    vertices.push_back(glm::vec4(maxx, maxy, maxz,1.0f)); // D = 3 

    vertices.push_back(glm::vec4(minx, miny, minz,1.0f)); // E = 4 
    vertices.push_back(glm::vec4(maxx, miny, minz,1.0f)); // F = 5 

    vertices.push_back(glm::vec4(minx, maxy, minz,1.0f)); // G = 6 
    vertices.push_back(glm::vec4(maxx, maxy, minz,1.0f)); // H = 7
    /*********************************/

    /* ADD FACES TO SPONGE LISTS */
    /* NOTE :: 2 triangle per each face (8 faces ) */
	int offset = menger_vertices.size() - 8;

    // front facing face
    faces.push_back(glm::uvec3(offset + 0, offset + 1, offset + 2));  // ABC
    faces.push_back(glm::uvec3(offset + 1, offset + 3, offset + 2));  // BDC
    
	// top facing face
    faces.push_back(glm::uvec3(offset + 2, offset + 3, offset + 6)); // CDG
    faces.push_back(glm::uvec3(offset + 6, offset + 3, offset + 7)); //  GDH
    
    // bottom facing face  
    faces.push_back(glm::uvec3(offset + 0, offset + 1, offset + 4));  // ABE
    faces.push_back(glm::uvec3(offset + 4,offset + 1, offset + 5));  //FBE TESTED SET {BEF, BFE, EFB, }

    // left facing face
    faces.push_back(glm::uvec3(offset + 0, offset + 6, offset + 2));  // ACG  TESTED SET { ACG, AGC }
    faces.push_back(glm::uvec3(offset + 0, offset +4, offset + 6));  // AEG TESTED SET { AEG} , 

    // right facing face 
    faces.push_back(glm::uvec3(offset+ 1, offset + 7, offset + 3)); //  BHD TESTED SET { BEH, BHE }
    faces.push_back(glm::uvec3(offset + 5, offset + 7,offset + 1 )); //  BHF 
   
    // rear facing face 
    faces.push_back(glm::uvec3(offset + 4, offset + 5, offset + 6)); //  EGF TESTED SET  { EGF,EFG } 
    faces.push_back(glm::uvec3(offset + 6, offset + 5, offset + 7)); //  GFH TESTED SET  {  }  

    /*********************************/

    }

    // RECURSIVE METHOD TO GENERATE ALL CUBE COORDINATES 
    // note :: test for method in test_spon_gen.cc
	// note :: you need to clear vertices and faces list
    void generateCubes(int level,
			    std::vector<glm::vec4>& vertices,
		    std::vector<glm::uvec3>& faces,
			    float minx, float miny, float minz,
			    float maxx, float maxy, float maxz ) {
	    int x, y, z;
	    if(level == 0){ // base case
		    drawCube(vertices,faces, 
			    minx, miny, minz, 
			    maxx, maxy, maxz);
	    } else {
		    // calculate intervals needed for cube 
		    float x_int = std::fabs(maxx - minx) / 3.0;
		    float y_int = std::fabs(maxy - miny) / 3.0;
		    float z_int = std::fabs(maxz - minz) / 3.0;

		    // Iterate over each cube's 27 sub cubes , generating 20 valid coords
			int pairs_offset = 0;
		    for(y = 0; y < 3; y++){
			    for(z = 0; z < 3; z++) {
				    for(x = 0; x < 3; x++){
						if(y == 0 || y == 2) {
							if(x != 1 || z != 1) { // negation of x == 1 && z == 1
								float new_x_min = minx + (x_int * x);
								float new_x_max = minx + (x_int * (x+1));
								float new_y_min = miny + (y_int * y);
								float new_y_max = miny + (y_int * (y+1));
								float new_z_min = minz + (z_int * z);
								float new_z_max = minz + (z_int * (z+1));

								generateCubes(level - 1,
										vertices, faces,
									new_x_min, new_y_min, new_z_min,
									new_x_max, new_y_max, new_z_max);
									pairs_offset++;	
							}
					    } 
						else if (y == 1) {
								if(x != 1 && z != 1) {
									float new_x_min = minx + (x_int * x);
									float new_x_max = minx + (x_int * (x+1));
									float new_y_min = miny + (y_int * y);
									float new_y_max = miny + (y_int * (y+1));
									float new_z_min = minz + (z_int * z);
									float new_z_max = minz + (z_int * (z+1));
								
								generateCubes(level - 1,
									vertices, faces,
									new_x_min, new_y_min, new_z_min,
									new_x_max, new_y_max, new_z_max);
						    	}	
					    }
				    }
			    }
		    }
	    }
    }

    void CreateMenger(std::vector<glm::vec4>& vertices,
		    std::vector<glm::uvec3>& faces) {
	    //std::cout << "Creating a Menger Sponge ... yeah right!\n";

		// you clear out your vertices only when changing levels (encapsulated here for better design)
		vertices.clear();
		faces.clear();

	    /* better idea --- see if the merger sponge coords set can be derived */
	    /* and then draw that merger sponge */
	    /* idea of a cube object to work with (helps abstract functionality of code) ? */

		// not sure if drawing the place should actually be here //
	    generateCubes(L,vertices, faces,
		    kMengerMinBounds.x, kMengerMinBounds.y, kMengerMinBounds.z,
		    kMengerMaxBounds.x, kMengerMaxBounds.y, kMengerMaxBounds.z);
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
          masses.push_back(new Mass(i,glm::vec3(menger_vertices[i])));
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

	std::vector<int>  getNeighbors(int x)
	{
		
		// 1. iterate over all vertices
    std::vector<int> closest_vertices;
    std::vector<DIST> distances;
    for(int i = 0; i< masses.size(); i++)
    {
       if(i != x) {
          float distance = glm::distance(masses[x]->pos, masses[i]->pos);
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

	 void Load_SpringSystem()
	 {
   //std::vector<glm::vec4> menger_vertices;
   //std::vector<glm::uvec3> menger_faces;

	  // vertex list has a direct relationship to the map list ( one-to-one)
	 	 for(int i = 0; i < masses.size(); i++)
	 	 {
       		std::vector<int> six_neighbors = getNeighbors(i); 
       		masses[i]->neighbors = six_neighbors;
	 	 }

     for(int x =0; x<masses.size(); x++)
     {
        std::cout<<masses[x]->m_id<<std::endl;
        //std::cout<<masses[x]->pos.x <<"   "<< masses[x]->pos.y<<"   "<< masses[x]->pos.z<<std::endl;
        for(int y = 0; y<masses[x]->neighbors.size(); ++y)
        {           
            std::cout<<masses[x]->neighbors[y];
        }
        std::cout<<"end neiboers \n"<<std::endl;
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

        for(int x =0; x<masses.size(); x++)
        {
            std::cout<<masses[x]->m_id<<std::endl;
            std::cout<<masses[x]->pos.x <<"   "<< masses[x]->pos.y<<"   "<< masses[x]->pos.z<<std::endl;
            for(int y = 0; y<masses[x]->springs.size(); ++y)
            {           
               std::cout<<(mass_springs[x*6+y].getMassB())->m_id;
            }
            std::cout<<"\n"<<std::endl;
        }
    }



    void ErrorCallback(int error, const char* description) {
    std::cerr << "GLFW Error: " << description << "\n";
    }

    void KeyCallback(GLFWwindow* window, int key, int scancode, int action,
		    int mods) {

    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
	glfwSetWindowShouldClose(window, GL_TRUE);
    else if (key == GLFW_KEY_W && action != GLFW_RELEASE) {

		// note :: the differences are what get's manipulated and how (center vs eye)

	//	if(!fps_mode) { // orbital mode

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

		// note :: ther DOES NOT EXIST A DISTINCTION BETWEEN 
		// FPS VS ORBITAL MODE WHEN MANIPULATING (CENTER|EYE)
		// SINCE THE CENTER IS NEVER TRACKED AND LOOK IS A VECTOR CONST IN DIRECTion

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

	// why do the keyboard presses no longer work though?
    else if (key == GLFW_KEY_0 && action != GLFW_RELEASE) {

	    L = 0;
		CreateMenger(menger_vertices, menger_faces);

  		CHECK_GL_ERROR(glBindVertexArray(array_objects[kMengerVao]));

		// each bind buffer call is separete 
  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kVertexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
					sizeof(float) * menger_vertices.size() * 4,
					&menger_vertices[0], GL_STATIC_DRAW));

  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kIndexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
					sizeof(uint32_t) * menger_faces.size() * 3,
					&menger_faces[0], GL_STATIC_DRAW));

    } else if (key == GLFW_KEY_1 && action != GLFW_RELEASE) {

		L = 1;
		CreateMenger(menger_vertices, menger_faces);

  		CHECK_GL_ERROR(glBindVertexArray(array_objects[kMengerVao]));

		// each bind buffer call is separete 
  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kVertexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
					sizeof(float) * menger_vertices.size() * 4,
					&menger_vertices[0], GL_STATIC_DRAW));

  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kIndexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
					sizeof(uint32_t) * menger_faces.size() * 3,
					&menger_faces[0], GL_STATIC_DRAW));


    } else if (key == GLFW_KEY_2 && action != GLFW_RELEASE) {

		L = 2;
		CreateMenger(menger_vertices, menger_faces);

  		CHECK_GL_ERROR(glBindVertexArray(array_objects[kMengerVao]));

		// each bind buffer call is separete 
  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kVertexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
					sizeof(float) * menger_vertices.size() * 4,
					&menger_vertices[0], GL_STATIC_DRAW));

  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kIndexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
					sizeof(uint32_t) * menger_faces.size() * 3,
					&menger_faces[0], GL_STATIC_DRAW));

    } else if (key == GLFW_KEY_3 && action != GLFW_RELEASE) {

		L = 3;
		CreateMenger(menger_vertices, menger_faces);

  		CHECK_GL_ERROR(glBindVertexArray(array_objects[kMengerVao]));

		// each bind buffer call is separete 
  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kVertexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
					sizeof(float) * menger_vertices.size() * 4,
					&menger_vertices[0], GL_STATIC_DRAW));

  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kIndexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
					sizeof(uint32_t) * menger_faces.size() * 3,
					&menger_faces[0], GL_STATIC_DRAW));

    } else if (key == GLFW_KEY_4 && action != GLFW_RELEASE) {

		L = 4;
		CreateMenger(menger_vertices, menger_faces);

  		CHECK_GL_ERROR(glBindVertexArray(array_objects[kMengerVao]));

		// each bind buffer call is separete 
  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kVertexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
					sizeof(float) * menger_vertices.size() * 4,
					&menger_vertices[0], GL_STATIC_DRAW));

  		CHECK_GL_ERROR(
   	   		glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kMengerVao][kIndexBuffer]));
	    CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
					sizeof(uint32_t) * menger_faces.size() * 3,
					&menger_faces[0], GL_STATIC_DRAW));
	}		
}

// I'm not sure what the significance of screen coordinates are here
	// i.e. why can't this just be in normal camera coordinates?
// what is orbital and fps mode?
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
		// note :: do not use radians!
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

	// but how to do this per each frame
	// in fact how do I know that zoom_speed works by frame?
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



	//Load geometries to render

  //CreateMenger(menger_vertices, menger_faces);
  std::string file_name = "obj/sphere.obj";
  LoadObj(file_name, menger_vertices, menger_faces);
  Load_SpringSystem();
  setupSpring();
  CreatePlane();
  std::cout << "Loaded plane and vertices geometries" << std::endl; 

  // Setup our VAO array
  CHECK_GL_ERROR(glGenVertexArrays(kNumVaos, array_objects));

///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////

// Setup the menger array object.
  // Switch to the kMenger VAO.
  // always have to switch and bind right buffer
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
  // Setup the plane array object.
  // Switch to the kPlaneVAO.
  CHECK_GL_ERROR(glBindVertexArray(array_objects[kPlaneVao]));
  // Generate buffer objects for kMengerVao
  CHECK_GL_ERROR(glGenBuffers(kNumVbos, &buffer_objects[kPlaneVao][0]));

  // Let's create our plane SHADER program.
  GLuint plane_program_id = 1;
  CHECK_GL_ERROR(plane_program_id = glCreateProgram());

  // Setup vertex data for kMenger VBOs
// why is this causing an issue?  (keyboard controls not working at this point )
  CHECK_GL_ERROR(
      glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kPlaneVao][plane_kVertexBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(float) * plane_vertices.size() * 4,
                              &plane_vertices[0], GL_STATIC_DRAW));
  CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
  CHECK_GL_ERROR(glEnableVertexAttribArray(0));

  // Setup element array buffer. (kMenger faces data )
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
  //CHECK_GL_ERROR(glAttachShader(plane_program_id, fragment_shader_id));
  CHECK_GL_ERROR(glAttachShader(plane_program_id, plane_fragment_shader_id));

  // Bind attributes. ( linking step )
  CHECK_GL_ERROR(glBindAttribLocation(plane_program_id, 0, "vertex_position"));
  CHECK_GL_ERROR(glBindFragDataLocation(plane_program_id, 0, "plane_fragment_color")); // need 0 (numerical pos for gpu) (things don't match!)
  //CHECK_GL_ERROR(glBindFragDataLocation(plane_program_id, 0, "fragment_color")); // need 0 (numerical pos for gpu) (things don't match!)

	// there is an issue with plane shader program itself at the moment
	
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

    // Switch to our Menger VAO.
    CHECK_GL_ERROR(glBindVertexArray(array_objects[kMengerVao]));

    // Use our program.
    CHECK_GL_ERROR(glUseProgram(program_id));

    // Pass uniforms in.
    CHECK_GL_ERROR(glUniformMatrix4fv(projection_matrix_location, 1, GL_FALSE,
                                      &projection_matrix[0][0]));
    CHECK_GL_ERROR(glUniformMatrix4fv(view_matrix_location, 1, GL_FALSE,
                                      &view_matrix[0][0]));
    CHECK_GL_ERROR(
        glUniform4fv(light_position_location, 1, &light_position[0]));

    // Draw our triangles.
    //CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES, menger_faces.size() * 3,
    CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES, menger_faces.size() * 3,
                                  GL_UNSIGNED_INT, 0));

    /////////////////////////////////////////////////// 
    /////////////////////////////////////////////////// 
    /////////////////////////////////////////////////// 
    // Switch to our plane Menger VAO.
    CHECK_GL_ERROR(glBindVertexArray(array_objects[kPlaneVao]));

    // Use our program.
    CHECK_GL_ERROR(glUseProgram(plane_program_id));

    // Pass uniforms in for the plane.  (maybe should be plane_(projection/view_matrix)_location and plane_light_position
    CHECK_GL_ERROR(glUniformMatrix4fv(plane_projection_matrix_location, 1, GL_FALSE,
                                      &projection_matrix[0][0]));
    CHECK_GL_ERROR(glUniformMatrix4fv(plane_view_matrix_location, 1, GL_FALSE,
                                      &view_matrix[0][0]));
    CHECK_GL_ERROR(
        glUniform4fv(plane_light_position_location, 1, &light_position[0]));
    //std::cout << light_position << std::endl;

    // Draw plane's triangles.
    CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES, plane_faces.size() * 3,
                                  GL_UNSIGNED_INT, 0));

	// Go back to the old VAO ( menger vertices and faces )
	// not needed to switch back 
  	//CHECK_GL_ERROR(glBindVertexArray(array_objects[kMengerVao]));
 
    // Poll and swap.
    glfwPollEvents();
    glfwSwapBuffers(window);
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
