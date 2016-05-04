#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/spline.hpp>
#include <glm/gtx/string_cast.hpp>
#include <json_spirit/json_spirit_reader_template.h>

const std::string window_title = "Scene Graph";
int window_width = 800, window_height = 600;
float aspect = static_cast<float>(window_width) / window_height;
float current_normal_length = 0.005f;
int current_key_frame = 0;
bool pause_enabled = true;
int max_key_frame = 0;

template <typename T>
struct Allocator {
  Allocator() : counts(), pointers() {}
  ~Allocator() {
    for (int i = 0; i < pointers.size(); ++i) {
      if (counts[i] > 1)
        delete[] pointers[i];
      else
        delete pointers[i];
    }
  }
  T* Allocate(int count) {
    if (count > 1)
      pointers.push_back(new T[count]);
    else
      pointers.push_back(new T);
    counts.push_back(count);
    return pointers.back();
  }
  std::vector<T*> pointers;
  std::vector<int> counts;
};

template <typename T>
T* Alloc(int count) {
  static Allocator<T> allocator;
  return allocator.Allocate(count);
}

const char* solid_vertex_shader =
    "#version 330 core\n"
    "in mat4 model;"
    "in vec3 diffuse_color;"
    "in int shade_mode;"
    "in int render_mode;"
    "in vec4 vertex_position;"
    "in vec4 vertex_normal;"
    "out vec3 vs_color;"
    "out vec4 vs_normal;"
    "out int vs_shade_mode;"
    "out int vs_render_mode;"
    "out vec4 vs_position;"
    "void main() {"
    "gl_Position = model * vertex_position;"
    "vs_color = diffuse_color;"
    "vs_normal = transpose(inverse(model)) * vertex_normal;"
    "vs_shade_mode = shade_mode;"
    "vs_render_mode = render_mode;"
    "vs_position = gl_Position;"
    "}";

const char* solid_geometry_shader =
    "#version 330 core\n"
    "layout (triangles) in;"
    "layout (triangle_strip, max_vertices = 3) out;"
    "uniform mat4 projection;"
    "uniform mat4 view;"
    "in vec3 vs_color[];"
    "in vec4 vs_normal[];"
    "in vec4 vs_position[];"
    "in int vs_shade_mode[];"
    "in int vs_render_mode[];"
    "out vec4 face_normal;"
    "out vec4 normal;"
    "out vec4 world_position;"
    "flat out int shading_mode;"
    "out vec3 color;"
    "void main() {"
    "if (vs_render_mode[0] != 0) {"
    " return;"
    "}"
    "int n = 0;"
    "vec3 a = gl_in[0].gl_Position.xyz;"
    "vec3 b = gl_in[1].gl_Position.xyz;"
    "vec3 c = gl_in[2].gl_Position.xyz;"
    "vec3 u = normalize(b - a);"
    "vec3 v = normalize(c - a);"
    "face_normal = normalize(vec4(normalize(cross(u, v)), 0.0));"
    "for (n = 0; n < gl_in.length(); n++) {"
    "world_position = vs_position[n];"
    "color = vs_color[n];"
    "normal = vs_normal[n];"
    "shading_mode = vs_shade_mode[n];"
    "gl_Position = projection * view * gl_in[n].gl_Position;"
    "EmitVertex();"
    "}"
    "EndPrimitive();"
    "}";

const char* solid_fragment_shader =
    "#version 330 core\n"
    "const int kMaxLights = 8;"
    "uniform mat4 lights[kMaxLights];"
    "uniform int num_lights;"
    "in vec4 normal;"
    "in vec4 face_normal;"
    "in vec3 color;"
    "in vec4 world_position;"
    "flat in int shading_mode;"
    "out vec4 fragment_color;"
    "void main() {"
    "vec3 tot = vec3(0);"
    "vec4 n = face_normal;"
    "float d;"
    "float dot_nl = 0.0;"
    "float atten = 1.0;"
    "fragment_color = vec4(color, 1.0);"
    "if (shading_mode != 0) {"
    "if (shading_mode == 2) {"
    "  n = normal;"
    "}"
    "float light_type;"
    "vec3 light_color;"
    "vec3 light_direction;"
    "vec3 light_position;"
    "vec3 a;"
    "mat4 light_matrix;"
    "for (int i = 0; i < kMaxLights; ++i) {"
    "  light_matrix = lights[i];"
    "  light_color = light_matrix[0].xyz;"
    "  light_position = light_matrix[1].xyz;"
    "  light_direction = light_matrix[2].xyz;"
    "  a = light_matrix[3].xyz;"
    "  light_type = light_matrix[3].w;"
    "  if (light_type >= 0) {"
    "    if (light_type == 0) {"
    "      light_direction = light_position - world_position.xyz;"
    "      d = length(light_direction);"
    "      atten = 1.0 / (a.x + d * a.y + d * d * a.z);"
    "    } else {"
    "      atten = 1.0;"
    "    }"
    "    dot_nl = dot(normalize(light_direction), normalize(n.xyz));"
    "    dot_nl = clamp(dot_nl, 0.0, 1.0);"
    "    tot +=  atten * dot_nl * light_color;"
    "  }"
    "}"
    "fragment_color = clamp(vec4(tot * color, 1.0), 0, 1);"
    "}"
    "}";

const char* wire_vertex_shader =
    "#version 330 core\n"
    "in mat4 model;"
    "in vec3 diffuse_color;"
    "in int render_mode;"
    "in vec4 vertex_position;"
    "out vec3 vs_color;"
    "out int vs_render_mode;"
    "void main() {"
    "gl_Position = model * vertex_position;"
    "vs_color = diffuse_color;"
    "vs_render_mode = render_mode;"
    "}";

const char* wire_geometry_shader =
    "#version 330 core\n"
    "layout (triangles) in;"
    "layout (line_strip, max_vertices = 3) out;"
    "uniform mat4 projection;"
    "uniform mat4 view;"
    "in vec3 vs_color[];"
    "in vec4 vs_position[];"
    "in int vs_render_mode[];"
    "out vec3 color;"
    "void main() {"
    "if (vs_render_mode[0] != 1) {"
    " return;"
    "}"
    "for (int n = 0; n < gl_in.length(); n++) {"
    "color = vs_color[n];"
    "gl_Position = projection * view * gl_in[n].gl_Position;"
    "EmitVertex();"
    "}"
    "EndPrimitive();"
    "}";

const char* wire_fragment_shader =
    "#version 330 core\n"
    "in vec3 color;"
    "out vec4 fragment_color;"
    "void main() {"
    " fragment_color = vec4(color, 1.0);"
    "}";

const char* normal_vertex_shader =
    "#version 330 core\n"
    "in mat4 model;"
    "in vec3 diffuse_color;"
    "in int normal_display_mode;"
    "in vec4 vertex_position;"
    "in int render_mode;"
    "in int shade_mode;"
    "in vec4 vertex_normal;"
    "out vec4 vs_vertex_normal;"
    "out vec3 vs_color;"
    "out int vs_normal_display_mode;"
    "void main() {"
    "gl_Position = model * vertex_position;"
    "vs_color = diffuse_color;"
    "vs_normal_display_mode = normal_display_mode;"
    "vs_vertex_normal = transpose(inverse(model)) * vertex_normal;"
    "}";

const char* normal_geometry_shader =
    "#version 330 core\n"
    "layout (triangles) in;"
    "layout (line_strip, max_vertices = 8) out;"
    "uniform mat4 projection;"
    "uniform mat4 view;"
    "uniform float normal_length;"
    "in vec3 vs_color[];"
    "in vec4 vs_position[];"
    "in int vs_normal_display_mode[];"
    "in vec4 vs_vertex_normal[];"
    "out vec3 color;"
    "void main() {"
    "vec3 a = gl_in[0].gl_Position.xyz;"
    "vec3 b = gl_in[1].gl_Position.xyz;"
    "vec3 c = gl_in[2].gl_Position.xyz;"
    "vec3 u = normalize(b - a);"
    "vec3 v = normalize(c - a);"
    "vec3 face_normal = normalize(cross(u, v));"
    "vec3 center = (a + b + c) / 3.0;"
    "for (int n = 0; n < 3; ++n) {"
    "if (vs_normal_display_mode[n] == 1 || vs_normal_display_mode[n] == 3) {"
    "  gl_Position = gl_in[n].gl_Position;"
    "  gl_Position = projection * view * gl_Position;"
    "  color = vec3(0, 1, 0);"
    "  EmitVertex();"
    "  gl_Position = gl_in[n].gl_Position + normal_length * "
    "vs_vertex_normal[n];"
    "  gl_Position = projection * view * gl_Position;"
    "  color = vec3(0, 1, 0);"
    "  EmitVertex();"
    "  EndPrimitive();"
    "}"
    "}"
    "if (vs_normal_display_mode[0] == 2 || vs_normal_display_mode[0] == 3) {"
    " gl_Position = vec4(center, 1.0);"
    " gl_Position = projection * view * gl_Position;"
    " color = vec3(0, 0, 1);"
    " EmitVertex();"
    " gl_Position = vec4(center, 1.0) + normal_length * vec4(face_normal,  "
    "0.0);"
    " gl_Position = projection * view * gl_Position;"
    " color = vec3(0, 0, 1);"
    " EmitVertex();"
    " EndPrimitive();"
    "}"
    "}";

const char* normal_fragment_shader =
    "#version 330 core\n"
    "in vec3 color;"
    "out vec4 fragment_color;"
    "void main() {"
    " fragment_color = vec4(color, 1.0);"
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

namespace glm {
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

std::ostream& operator<<(std::ostream& os, const glm::mat4& v) {
  os << glm::to_string(v);
  return os;
}

std::ostream& operator<<(std::ostream& os, const glm::mat3& v) {
  os << glm::to_string(v);
  return os;
}
}  // namespace glm

struct TransformFrame {
  glm::vec3 translate;
  glm::vec4 rotate;
  glm::vec3 scale;
  int t;
};

struct TransformFrames {
  TransformFrame* frames;
  int num_frames;
};

struct MaterialFrame {
  glm::vec3 diffuse_color;
  int t;
};

struct MaterialFrames {
  MaterialFrame* frames;
  int num_frames;
};

struct LightFrame {
  glm::vec3 color;
  int t;
};

struct LightFrames {
  LightFrame* frames;
  int num_frames;
};

struct BoundingBox {
  BoundingBox()
      : min(glm::vec3(-std::numeric_limits<float>::max())),
        max(glm::vec3(std::numeric_limits<float>::max())) {}
  glm::vec3 min;
  glm::vec3 max;
};

struct Mesh {
  Mesh()
      : id(0),
        num_instances(0),
        vertices(),
        faces(),
        vertex_normals(),
        face_normals(),
        model_matrices(),
        bounds() {}
  ~Mesh() {}
  int id;
  int num_instances;
  std::vector<glm::vec4> vertices;
  std::vector<glm::uvec3> faces;
  std::vector<glm::vec4> vertex_normals;
  std::vector<glm::vec4> face_normals;
  // These are instanced: one value per instance of the mesh.
  std::vector<glm::mat4> model_matrices;
  std::vector<glm::vec3> diffuse_colors;
  std::vector<int> shade_modes;
  std::vector<int> render_modes;
  std::vector<int> normal_display_modes;
  BoundingBox bounds;
};

void LoadObj(const std::string& file, std::vector<glm::vec4>& vertices,
             std::vector<glm::uvec3>& indices) {
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

std::ostream& operator<<(std::ostream& os, const TransformFrame& frame) {
  os << "[" << frame.t << "] translate = " << frame.translate
     << " rotate = " << frame.rotate << " scale = " << frame.scale << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const TransformFrames& frames) {
  os << "frames = ";
  if (frames.frames == nullptr || frames.num_frames <= 0)
    os << " null \n";
  else
    os << "\n";
  for (int i = 0; i < frames.num_frames; ++i) os << frames.frames[i] << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const MaterialFrame& frame) {
  os << "[" << frame.t << "]  diffuse_color= " << frame.diffuse_color << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const MaterialFrames& frames) {
  os << "frames = ";
  if (frames.frames == nullptr || frames.num_frames <= 0)
    os << " null \n";
  else
    os << "\n";
  for (int i = 0; i < frames.num_frames; ++i) os << frames.frames[i] << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const LightFrame& frame) {
  os << "[" << frame.t << "]  color= " << frame.color << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const LightFrames& frames) {
  os << "frames = ";
  if (frames.frames == nullptr || frames.num_frames <= 0)
    os << " null \n";
  else
    os << "\n";
  for (int i = 0; i < frames.num_frames; ++i) os << frames.frames[i] << "\n";
  return os;
}

namespace scene_graph {

struct SceneNode {
  enum NodeType {
    kAttribute,
    kCamera,
    kGeometry,
    kLight,
    kObject,
    kTransform,
    kNumObjectTypes
  };

  // For attribute nodes.
  struct Attribute {
    Attribute()
        : has_shade_mode(false),
          shade_mode(kShadeFlat),
          has_render_mode(false),
          render_mode(kRenderModeSolid),
          has_interpolate_mode(false),
          interpolate_mode(kInterpolateLinear),
          has_normal_display_mode(false),
          normal_display_mode(kNormalDisplayNone),
          has_frames(false),
          frames() {}
    // Lighting mode attribute.
    enum ShadeMode {
      kShadeNone,
      kShadeFlat,
      kShadePhong,
      kNumShadeModes
    };

    // Model rendering attribute.
    enum RenderMode {
      kRenderModeSolid,
      kRenderModeWireFrame,
      kNumRenderModes
    };

    // Keyframe interpolation attribute.
    enum InterpolateMode {
      kInterpolateLinear,
      kInterpolateSpline,
      kNumInterpolateModes
    };

    // Normal display attribute.
    enum NormalDisplayMode {
      kNormalDisplayNone,
      kNormalDisplayFace,
      kNormalDisplayVertex,
      kNormalDisplayAll,
      kNumNormalDisplayModes
    };

    static const char* kShadeModeStringValues[kNumShadeModes];
    static const char* kRenderModeStringValues[kNumRenderModes];
    static const char* kInterpolateModeStringValues[kNumInterpolateModes];
    static const char* kNormalDisplayModeStringValues[kNumNormalDisplayModes];

    bool has_shade_mode;
    ShadeMode shade_mode;
    bool has_render_mode;
    RenderMode render_mode;
    bool has_interpolate_mode;
    InterpolateMode interpolate_mode;
    bool has_normal_display_mode;
    NormalDisplayMode normal_display_mode;

    // Material attributes.
    bool has_frames;
    MaterialFrames frames;
  };

  struct Camera {
    enum ProjectionType {
      kProjectionOrthographic,
      kProjectionPerspective,
      kNumProjectionTypes
    };
    static const char* kProjectionTypeStringValues[kNumProjectionTypes];
    ProjectionType projection;
    float fov;     // for perspective
    float height;  // for orthographic
    float near;
    float far;
    int camera_id;
  };

  struct Geometry {
    const char* mesh;
    int mesh_id;
    int instance_id;
  };

  struct Light {
    enum LightType {
      kLightPoint,
      kLightDirectional,
      kNumLightTypes
    };
    static const char* kLightTypeStringValues[kNumLightTypes];
    LightType light_type;
    LightFrames frames;
    glm::vec3 attenuation;  // 0 - constant, 1 - linear, 2 - quadratic
    glm::vec3 position;     // recomputed
    glm::vec3 direction;    // recomputed
    glm::vec3 color;        // recomputed
    int light_id;           // in range [0,  kMaxLights - 1]
  };

  struct Transform {
    // Keyframed rotations.
    TransformFrames frames;
  };

  SceneNode()
      : type(kObject),
        children(nullptr),
        num_children(0),
        parent(nullptr),
        model(glm::mat4(1.0f)),
        attribute_info() {}
  SceneNode(const SceneNode& node)
      : type(node.type),
        children(node.children),
        num_children(node.num_children),
        parent(node.parent),
        model(node.model),
        attribute_info(node.attribute_info) {
    if (node.type == kAttribute)
      attribute = node.attribute;
    else if (node.type == kCamera)
      camera = node.camera;
    else if (node.type == kGeometry)
      geometry = node.geometry;
    else if (node.type == kLight)
      light = node.light;
    else if (node.type == kTransform)
      transform = node.transform;
  }

  ~SceneNode() {}

  int type;
  SceneNode* children;
  int num_children;
  SceneNode* parent;
  glm::mat4 model;
  Attribute attribute_info;

  union {
    Attribute attribute;
    Camera camera;
    Geometry geometry;
    Light light;
    Transform transform;
  };
};

const char* SceneNode::Attribute::kShadeModeStringValues
    [SceneNode::Attribute::ShadeMode::kNumShadeModes] = {"none", "flat",
                                                         "phong"};
const char* SceneNode::Attribute::kRenderModeStringValues
    [SceneNode::Attribute::RenderMode::kNumRenderModes] = {"solid", "wire"};

const char* SceneNode::Attribute::kInterpolateModeStringValues
    [SceneNode::Attribute::InterpolateMode::kNumInterpolateModes] = {"linear",
                                                                     "spline"};

const char* SceneNode::Attribute::kNormalDisplayModeStringValues
    [SceneNode::Attribute::NormalDisplayMode::kNumNormalDisplayModes] = {
        "none", "vertex", "face", "both"};

const char* SceneNode::Light::kLightTypeStringValues
    [SceneNode::Light::LightType::kNumLightTypes] = {"point", "directional"};

const char* SceneNode::Camera::kProjectionTypeStringValues
    [SceneNode::Camera::kNumProjectionTypes] = {"orthographic", "perspective"};

struct SceneGraph {
  // Buffer types, attributes first, then uniforms.
  enum BufferTypes {
    kVertexBuffer,             // attribute
    kNormalBuffer,             // attribute
    kDiffuseColorBuffer,       // attribute
    kShadeModeBuffer,          // attribute
    kRenderModeBuffer,         // attribute
    kNormalDisplayModeBuffer,  // attribute
    kModelMatrixBuffer,        // attribute
    kIndexBuffer,              // non-attribute
    kNumBufferTypes
  };
  SceneGraph()
      : root(nullptr),
        meshes(),
        mesh_nodes(),
        light_nodes(),
        camera_nodes(kMaxCameras, nullptr),
        mesh_lookup(),
        array_objects(),
        buffer_objects(),
        light_info(kMaxLights,
                   glm::mat4(glm::vec4(0.0f), glm::vec4(0.0f), glm::vec4(0.0f),
                             glm::vec4(0.0f, 0.0f, 0.0f, -1.0f))),
        current_camera(-1),
        num_lights(0) {}
  ~SceneGraph() {}
  SceneNode* root;
  std::vector<Mesh*> meshes;
  std::vector<SceneNode*> mesh_nodes;
  std::vector<SceneNode*> light_nodes;
  std::vector<SceneNode*> camera_nodes;
  std::unordered_map<std::string, Mesh*> mesh_lookup;
  std::vector<GLuint> array_objects;
  std::vector<GLuint> buffer_objects;
  std::vector<glm::mat4> light_info;
  int current_camera;
  int num_lights;
  const static int kMaxCameras;
  const static int kMaxLights;
};
const int SceneGraph::kMaxCameras = 10;
const int SceneGraph::kMaxLights = 8;

struct ProgramInfo {
  GLuint program_id;
  GLuint vertex_shader_id;
  GLuint fragment_shader_id;
  GLuint geometry_shader_id;
  const char* vertex_shader_source;
  const char* fragment_shader_source;
  const char* geometry_shader_source;
  GLint projection_matrix_location;
  GLint view_matrix_location;
  GLint lights_location;
  GLint normal_length_location;
};

void InitProgramGL(const char* vertex_shader, const char* fragment_shader,
                   const char* geometry_shader, ProgramInfo* program_info) {
  // Shader setup.

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

  // Program Setup

  // Let's create triangle our program.
  GLuint program_id = 0;
  CHECK_GL_ERROR(program_id = glCreateProgram());
  CHECK_GL_ERROR(glAttachShader(program_id, vertex_shader_id));
  CHECK_GL_ERROR(glAttachShader(program_id, fragment_shader_id));
  CHECK_GL_ERROR(glAttachShader(program_id, geometry_shader_id));

  // Bind attributes.
  CHECK_GL_ERROR(glBindAttribLocation(program_id, SceneGraph::kVertexBuffer,
                                      "vertex_position"));
  CHECK_GL_ERROR(glBindAttribLocation(program_id, SceneGraph::kNormalBuffer,
                                      "vertex_normal"));
  CHECK_GL_ERROR(glBindAttribLocation(
      program_id, SceneGraph::kDiffuseColorBuffer, "diffuse_color"));
  CHECK_GL_ERROR(glBindAttribLocation(program_id, SceneGraph::kShadeModeBuffer,
                                      "shade_mode"));
  CHECK_GL_ERROR(glBindAttribLocation(program_id, SceneGraph::kRenderModeBuffer,
                                      "render_mode"));
  CHECK_GL_ERROR(glBindAttribLocation(
      program_id, SceneGraph::kNormalDisplayModeBuffer, "normal_display_mode"));
  CHECK_GL_ERROR(glBindAttribLocation(program_id,
                                      SceneGraph::kModelMatrixBuffer, "model"));
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
  GLint lights_location = 0;
  CHECK_GL_ERROR(lights_location = glGetUniformLocation(program_id, "lights"));

  GLint normal_length_location = 0;
  CHECK_GL_ERROR(normal_length_location =
                     glGetUniformLocation(program_id, "normal_length"));

  // Save the data about this program instance.
  program_info->program_id = program_id;
  program_info->vertex_shader_id = vertex_shader_id;
  program_info->fragment_shader_id = fragment_shader_id;
  program_info->geometry_shader_id = geometry_shader_id;
  program_info->vertex_shader_source = vertex_source_pointer;
  program_info->fragment_shader_source = fragment_source_pointer;
  program_info->geometry_shader_source = geometry_source_pointer;
  program_info->projection_matrix_location = projection_matrix_location;
  program_info->view_matrix_location = view_matrix_location;
  program_info->lights_location = lights_location;
  program_info->normal_length_location = normal_length_location;
}

void InitMeshGL(GLuint* vao, GLuint* buffer_objects, Mesh* mesh) {
  // Setup our VAO.
  CHECK_GL_ERROR(glGenVertexArrays(1, vao));

  // Switch to the VAO.
  CHECK_GL_ERROR(glBindVertexArray(*vao));

  // Generate buffer objects
  CHECK_GL_ERROR(glGenBuffers(SceneGraph::kNumBufferTypes, buffer_objects));

  // Setup vertex data in a VBO.
  CHECK_GL_ERROR(
      glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[SceneGraph::kVertexBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(glm::vec4) * mesh->vertices.size(),
                              &mesh->vertices[0], GL_STATIC_DRAW));
  CHECK_GL_ERROR(glVertexAttribPointer(SceneGraph::kVertexBuffer, 4, GL_FLOAT,
                                       GL_FALSE, 0, 0));
  CHECK_GL_ERROR(glEnableVertexAttribArray(SceneGraph::kVertexBuffer));

  // Setup normal data in a VBO.
  CHECK_GL_ERROR(
      glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[SceneGraph::kNormalBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(glm::vec4) * mesh->vertices.size(),
                              &mesh->vertex_normals[0], GL_STATIC_DRAW));
  CHECK_GL_ERROR(glVertexAttribPointer(SceneGraph::kNormalBuffer, 4, GL_FLOAT,
                                       GL_FALSE, 0, 0));
  CHECK_GL_ERROR(glEnableVertexAttribArray(SceneGraph::kNormalBuffer));

  // Setup color data in a VBO.
  // NOTE: These are instanced, to support multiple instances of each mesh.
  CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER,
                              buffer_objects[SceneGraph::kDiffuseColorBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(glm::vec3) * mesh->diffuse_colors.size(),
                              &mesh->diffuse_colors[0], GL_STATIC_DRAW));
  // Enable this attribute.
  CHECK_GL_ERROR(glEnableVertexAttribArray(SceneGraph::kDiffuseColorBuffer));
  // We need one attribute.
  CHECK_GL_ERROR(glVertexAttribPointer(SceneGraph::kDiffuseColorBuffer, 3,
                                       GL_FLOAT, GL_FALSE, sizeof(glm::vec3),
                                       0));
  // Instance this color.
  CHECK_GL_ERROR(glVertexAttribDivisor(SceneGraph::kDiffuseColorBuffer, 1));

  // Setup shade mode data in a VBO.
  // NOTE: These are instanced, to support multiple instances of each mesh.
  CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER,
                              buffer_objects[SceneGraph::kShadeModeBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(int) * mesh->shade_modes.size(),
                              &mesh->shade_modes[0], GL_STATIC_DRAW));
  // Enable this attribute.
  CHECK_GL_ERROR(glEnableVertexAttribArray(SceneGraph::kShadeModeBuffer));
  // We need one attribute.
  CHECK_GL_ERROR(glVertexAttribIPointer(SceneGraph::kShadeModeBuffer, 1, GL_INT,
                                        sizeof(int), 0));
  // Instance this smooth shading option.
  CHECK_GL_ERROR(glVertexAttribDivisor(SceneGraph::kShadeModeBuffer, 1));

  // Setup render mode data in a VBO.
  // NOTE: These are instanced, to support multiple instances of each mesh.
  CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER,
                              buffer_objects[SceneGraph::kRenderModeBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(int) * mesh->render_modes.size(),
                              &mesh->render_modes[0], GL_STATIC_DRAW));

  // Enable this attribute.
  CHECK_GL_ERROR(glEnableVertexAttribArray(SceneGraph::kRenderModeBuffer));
  // We need one attribute.
  CHECK_GL_ERROR(glVertexAttribIPointer(SceneGraph::kRenderModeBuffer, 1,
                                        GL_INT, sizeof(int), 0));
  // Instance this smooth shading option.
  CHECK_GL_ERROR(glVertexAttribDivisor(SceneGraph::kRenderModeBuffer, 1));

  // Setup normal display mode data in a VBO.
  // NOTE: These are instanced, to support multiple instances of each mesh.
  CHECK_GL_ERROR(glBindBuffer(
      GL_ARRAY_BUFFER, buffer_objects[SceneGraph::kNormalDisplayModeBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(int) * mesh->normal_display_modes.size(),
                              &mesh->normal_display_modes[0], GL_STATIC_DRAW));

  // Enable this attribute.
  CHECK_GL_ERROR(
      glEnableVertexAttribArray(SceneGraph::kNormalDisplayModeBuffer));
  // We need one attribute.
  CHECK_GL_ERROR(glVertexAttribIPointer(SceneGraph::kNormalDisplayModeBuffer, 1,
                                        GL_INT, sizeof(int), 0));
  // Instance this smooth shading option.
  CHECK_GL_ERROR(
      glVertexAttribDivisor(SceneGraph::kNormalDisplayModeBuffer, 1));

  // Setup matrix data in a VBO.
  // NOTE: These are instanced, to support multiple instances of each mesh.
  CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER,
                              buffer_objects[SceneGraph::kModelMatrixBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(glm::mat4) * mesh->model_matrices.size(),
                              &mesh->model_matrices[0], GL_STATIC_DRAW));
  // Loop over each column of the matrix.
  for (int i = 0; i < 4; ++i) {
    // Enable each attribute.
    CHECK_GL_ERROR(
        glEnableVertexAttribArray(SceneGraph::kModelMatrixBuffer + i));
    // We need one attribute per column.
    CHECK_GL_ERROR(glVertexAttribPointer(SceneGraph::kModelMatrixBuffer + i, 4,
                                         GL_FLOAT, GL_FALSE, sizeof(glm::mat4),
                                         (void*)(sizeof(glm::vec4) * i)));
    // Instance this matrix.
    CHECK_GL_ERROR(
        glVertexAttribDivisor(SceneGraph::kModelMatrixBuffer + i, 1));
  }

  // Setup element array buffer.
  CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
                              buffer_objects[SceneGraph::kIndexBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                              sizeof(glm::uvec3) * mesh->faces.size(),
                              &mesh->faces[0], GL_STATIC_DRAW));
}

void InitSceneGraphGL(SceneGraph* scene_graph) {
  int num_meshes = scene_graph->meshes.size();

  // We will just do one VAO per mesh.
  scene_graph->array_objects.resize(num_meshes);

  // Each VAO will have four VBOs: vertices, indices,  normals, and model
  // matrices (instanced).
  scene_graph->buffer_objects.resize(num_meshes * SceneGraph::kNumBufferTypes);

  for (int i = 0; i < scene_graph->meshes.size(); ++i) {
    InitMeshGL(&scene_graph->array_objects[i],
               &scene_graph->buffer_objects[i * SceneGraph::kNumBufferTypes],
               scene_graph->meshes[i]);
  }
}

glm::vec4 deBoorsTransform(int t, int i, int k, TransformFrame* u) {
  if (k == 0)
    return u[i].rotate;
  else {
    glm::vec4 vPrevRot = deBoorsTransform(t, i, k - 1, u);
    glm::vec3 prevAxis = glm::normalize(glm::vec3(vPrevRot));
    float prevAngle = glm::radians(vPrevRot.w);
    glm::quat qPrevRot = angleAxis(prevAngle, prevAxis);
    glm::vec4 vThisRot = deBoorsTransform(t, i + 1, k - 1, u);
    glm::vec3 thisAxis = glm::normalize(glm::vec3(vThisRot));
    float thisAngle = glm::radians(vThisRot.w);
    glm::quat qThisRot = angleAxis(thisAngle, thisAxis);
    float ratioTicks =
        ((float)(u[i + 2].t - t)) / (u[i + 2].t - u[i - 2 + k].t);
    glm::quat qSlerpRot = glm::slerp(qThisRot, qPrevRot, ratioTicks);
    /* glm::slerp occasionally delivers NaNs. lerp if this happens */
    if ((isnan(glm::axis(qSlerpRot).x)) || (isnan(glm::axis(qSlerpRot).y)) ||
        (isnan(glm::axis(qSlerpRot).z)) || (isnan(glm::angle(qSlerpRot))))
      qSlerpRot = glm::lerp(qThisRot, qPrevRot, ratioTicks);
    return glm::vec4(glm::axis(qSlerpRot), glm::degrees(glm::angle(qSlerpRot)));
  }
}

glm::vec3 deBoorsMaterial(int t, int i, int k, MaterialFrame* u) {
  if (k == 0) {
    return u[i].diffuse_color;
  } else {
    glm::vec3 below = deBoorsMaterial(t, i, k - 1, u);
    glm::vec3 above = deBoorsMaterial(t, i + 1, k - 1, u);
    float alpha = ((float)(u[i + 2].t - t)) / (u[i + 2].t - u[i - 2 + k].t);
    float beta = ((float)(t - u[i - 2 + k].t)) / (u[i + 2].t - u[i - 2 + k].t);
    return (alpha * below + (1.f - alpha) * above);
  }
}

glm::vec3 deBoorsLight(int t, int i, int k, LightFrame* u) {
  if (k == 0) {
    return u[i].color;
  } else {
    glm::vec3 below = deBoorsLight(t, i, k - 1, u);
    glm::vec3 above = deBoorsLight(t, i + 1, k - 1, u);
    float alpha = ((float)(u[i + 2].t - t)) / (u[i + 2].t - u[i - 2 + k].t);
    float beta = ((float)(t - u[i - 2 + k].t)) / (u[i + 2].t - u[i - 2 + k].t);
    return (alpha * below + (1.f - alpha) * above);
  }
}

glm::vec3 deBoorsTransform(int t, int i, int k, TransformFrame* u, bool b) {
  if (k == 0) {
    return ((b) ? u[i].scale : u[i].translate);
  } else {
    glm::vec3 below = deBoorsTransform(t, i, k - 1, u, b);
    glm::vec3 above = deBoorsTransform(t, i + 1, k - 1, u, b);
    float alpha = ((float)(u[i + 2].t - t)) / (u[i + 2].t - u[i - 2 + k].t);
    float beta = ((float)(t - u[i - 2 + k].t)) / (u[i + 2].t - u[i - 2 + k].t);
    return (alpha * below + (1.f - alpha) * above);
  }
}

void InterpolateLight(int this_frame, int interpolate_mode,
                      const SceneNode* light_node, glm::vec3* color) {
  if (light_node->light.frames.frames == nullptr ||
      light_node->light.frames.num_frames == 0) {
    *color = glm::vec3(1.0f);
    return;
  }
  LightFrame* keyFrames = light_node->light.frames.frames;
  int numKeyFrames = light_node->light.frames.num_frames;
  int nextKeyFrameIdx = 0;
  while ((nextKeyFrameIdx < numKeyFrames) &&
         (keyFrames[nextKeyFrameIdx].t <= this_frame)) {
    ++nextKeyFrameIdx;
  }
  if (numKeyFrames == nextKeyFrameIdx) {
    --nextKeyFrameIdx;
  }
  if ((nextKeyFrameIdx <= 0)) {
    *color = glm::vec3(keyFrames[0].color);
    return;
  } else if (this_frame > keyFrames[numKeyFrames - 1].t) {
    *color = glm::vec3(keyFrames[numKeyFrames - 1].color);
  } else {
    if (interpolate_mode == light_node->attribute.kInterpolateLinear) {
      float ratioTicks =
          ((float)(this_frame - keyFrames[nextKeyFrameIdx - 1].t) /
           ((keyFrames[nextKeyFrameIdx].t - keyFrames[nextKeyFrameIdx - 1].t)));
      *color = glm::mix(keyFrames[nextKeyFrameIdx - 1].color,
                        keyFrames[nextKeyFrameIdx].color, ratioTicks);

    } else if (interpolate_mode == light_node->attribute.kInterpolateSpline) {
      LightFrame* paddedkeyFrames = new LightFrame[numKeyFrames + 4];
      paddedkeyFrames[0] = keyFrames[0];
      paddedkeyFrames[1] = keyFrames[0];
      for (int idx = 0; idx < numKeyFrames; idx++)
        paddedkeyFrames[idx + 2] = keyFrames[idx];
      paddedkeyFrames[numKeyFrames + 2] = keyFrames[numKeyFrames - 1];
      paddedkeyFrames[numKeyFrames + 3] = keyFrames[numKeyFrames - 1];
      *color = deBoorsLight(this_frame, nextKeyFrameIdx, 3, paddedkeyFrames);
      delete[] paddedkeyFrames;
    }
  }
}

void InterpolateMaterial(int this_frame, int interpolate_mode,
                         const SceneNode* material_node,
                         glm::vec3* diffuse_color) {
  if (!material_node->attribute_info.has_frames ||
      material_node->attribute_info.frames.frames == nullptr ||
      material_node->attribute_info.frames.num_frames == 0) {
    *diffuse_color = glm::vec3(0.5f);
    return;
  }
  MaterialFrame* keyFrames = material_node->attribute_info.frames.frames;
  int numKeyFrames = material_node->attribute_info.frames.num_frames;
  int nextKeyFrameIdx = 0;
  while ((nextKeyFrameIdx < numKeyFrames) &&
         (keyFrames[nextKeyFrameIdx].t <= this_frame)) {
    ++nextKeyFrameIdx;
  }
  if (numKeyFrames == nextKeyFrameIdx) {
    --nextKeyFrameIdx;
  }
  if ((nextKeyFrameIdx <= 0)) {
    *diffuse_color = glm::vec3(keyFrames[0].diffuse_color);
    return;
  } else if (this_frame > keyFrames[numKeyFrames - 1].t) {
    *diffuse_color = glm::vec3(keyFrames[numKeyFrames - 1].diffuse_color);
  } else {
    if (interpolate_mode == material_node->attribute.kInterpolateLinear) {
      float ratioTicks =
          ((float)(this_frame - keyFrames[nextKeyFrameIdx - 1].t) /
           ((keyFrames[nextKeyFrameIdx].t - keyFrames[nextKeyFrameIdx - 1].t)));
      *diffuse_color =
          glm::mix(keyFrames[nextKeyFrameIdx - 1].diffuse_color,
                   keyFrames[nextKeyFrameIdx].diffuse_color, ratioTicks);

    } else if (interpolate_mode ==
               material_node->attribute.kInterpolateSpline) {
      MaterialFrame* paddedkeyFrames = new MaterialFrame[numKeyFrames + 4];
      paddedkeyFrames[0] = keyFrames[0];
      paddedkeyFrames[1] = keyFrames[0];
      for (int idx = 0; idx < numKeyFrames; idx++)
        paddedkeyFrames[idx + 2] = keyFrames[idx];
      paddedkeyFrames[numKeyFrames + 2] = keyFrames[numKeyFrames - 1];
      paddedkeyFrames[numKeyFrames + 3] = keyFrames[numKeyFrames - 1];
      *diffuse_color =
          deBoorsMaterial(this_frame, nextKeyFrameIdx, 3, paddedkeyFrames);
      delete[] paddedkeyFrames;
    }
  }
}

void InterpolateTransform(int this_frame, int interpolate_mode,
                          const SceneNode* transform_node, glm::vec4* rotate,
                          glm::vec3* scale, glm::vec3* translate) {
  if (transform_node->transform.frames.frames == nullptr ||
      transform_node->transform.frames.num_frames == 0) {
    *scale = glm::vec3(1.0f);
    *rotate = glm::vec4(0.0f, 1.0f, 0.0f, 0.0f);
    *translate = glm::vec3(1.0f);
    return;
  }
  TransformFrame* keyFrames = transform_node->transform.frames.frames;
  int numKeyFrames = transform_node->transform.frames.num_frames;
  int nextKeyFrameIdx = 0;
  while ((nextKeyFrameIdx < numKeyFrames) &&
         (keyFrames[nextKeyFrameIdx].t <= this_frame)) {
    ++nextKeyFrameIdx;
  }
  if (numKeyFrames == nextKeyFrameIdx) {
    --nextKeyFrameIdx;
  }
  if ((nextKeyFrameIdx <= 0)) {
    *rotate = glm::vec4(keyFrames[0].rotate);
    *scale = glm::vec3(keyFrames[0].scale);
    *translate = glm::vec3(keyFrames[0].translate);
    return;
  } else if (this_frame > keyFrames[numKeyFrames - 1].t) {
    *rotate = glm::vec4(keyFrames[numKeyFrames - 1].rotate);
    *scale = glm::vec3(keyFrames[numKeyFrames - 1].scale);
    *translate = glm::vec3(keyFrames[numKeyFrames - 1].translate);
  } else {
    if (interpolate_mode == transform_node->attribute.kInterpolateLinear) {
      float ratioTicks =
          ((float)(this_frame - keyFrames[nextKeyFrameIdx - 1].t) /
           ((keyFrames[nextKeyFrameIdx].t - keyFrames[nextKeyFrameIdx - 1].t)));

      *scale = glm::mix(keyFrames[nextKeyFrameIdx - 1].scale,
                        keyFrames[nextKeyFrameIdx].scale, ratioTicks);
      *translate = glm::mix(keyFrames[nextKeyFrameIdx - 1].translate,
                            keyFrames[nextKeyFrameIdx].translate, ratioTicks);

      // previous rotation -> quaternion
      glm::vec3 prevAxis =
          glm::normalize(glm::vec3(keyFrames[nextKeyFrameIdx - 1].rotate));
      float prevAngle = glm::radians(keyFrames[nextKeyFrameIdx - 1].rotate.w);
      glm::quat qPrevRot = angleAxis(prevAngle, prevAxis);

      // next rotation -> quaternion
      glm::vec3 nextAxis =
          glm::normalize(glm::vec3(keyFrames[nextKeyFrameIdx].rotate));
      float nextAngle = glm::radians(keyFrames[nextKeyFrameIdx].rotate.w);
      glm::quat qNextRot = angleAxis(nextAngle, nextAxis);

      // LERP
      glm::quat qInterpolatedRot = glm::lerp(qPrevRot, qNextRot, ratioTicks);
      *rotate = glm::vec4(glm::axis(qInterpolatedRot),
                          glm::degrees(glm::angle(qInterpolatedRot)));

    } else if (interpolate_mode ==
               transform_node->attribute.kInterpolateSpline) {
      TransformFrame* paddedkeyFrames = new TransformFrame[numKeyFrames + 4];
      paddedkeyFrames[0] = keyFrames[0];
      paddedkeyFrames[1] = keyFrames[0];
      for (int idx = 0; idx < numKeyFrames; idx++)
        paddedkeyFrames[idx + 2] = keyFrames[idx];
      paddedkeyFrames[numKeyFrames + 2] = keyFrames[numKeyFrames - 1];
      paddedkeyFrames[numKeyFrames + 3] = keyFrames[numKeyFrames - 1];
      *scale = deBoorsTransform(this_frame, nextKeyFrameIdx, 3, paddedkeyFrames,
                                true);
      *translate = deBoorsTransform(this_frame, nextKeyFrameIdx, 3,
                                    paddedkeyFrames, false);
      *rotate =
          deBoorsTransform(this_frame, nextKeyFrameIdx, 3, paddedkeyFrames);
      delete[] paddedkeyFrames;
    }
  }
}

void CombineAttributes(const SceneNode::Attribute& inherited,
                       const SceneNode::Attribute& overrides,
                       SceneNode::Attribute* destination) {
  *destination = inherited;
  if (overrides.has_shade_mode) {
    destination->has_shade_mode = true;
    destination->shade_mode = overrides.shade_mode;
  }
  if (overrides.has_render_mode) {
    destination->has_render_mode = true;
    destination->render_mode = overrides.render_mode;
  }
  if (overrides.has_interpolate_mode) {
    destination->has_interpolate_mode = true;
    destination->interpolate_mode = overrides.interpolate_mode;
  }
  if (overrides.has_normal_display_mode) {
    destination->has_normal_display_mode = true;
    destination->normal_display_mode = overrides.normal_display_mode;
  }
  // NOTE: need to delete carefully.
  if (overrides.has_frames) {
    destination->has_frames = true;
    destination->frames = overrides.frames;
  }
}

void UpdateSceneNode(int key_frame, const glm::mat4& model,
                     const SceneNode::Attribute& attribute_info,
                     SceneGraph* scene_graph, SceneNode* node) {
  // std::cout << "UpdateSceneNode\n";
  node->model = model;
  node->attribute_info = attribute_info;
  // std::cout << "attribute_info = " << attribute_info.frames << "\n";
  switch (node->type) {
    case SceneNode::kAttribute: {
      // std::cout << "attribute\n";
      CombineAttributes(attribute_info, node->attribute, &node->attribute_info);
    } break;
    case SceneNode::kCamera:
      // std::cout << "camera\n";
      break;
    case SceneNode::kGeometry: {
      int mesh_id = node->geometry.mesh_id;
      int instance_id = node->geometry.instance_id;
      // std::cout << "geometry: mesh_id = " << mesh_id
      //<< " instance_id = " << instance_id
      //<< " shade_mode = " << node->attribute_info.shade_mode << "\n";
      Mesh* mesh = scene_graph->meshes[mesh_id];
      mesh->model_matrices[instance_id] = node->model;
      MaterialFrame material_frame;
      // std::cout << "attribute_info = " << attribute_info.frames
      //<< " has_frames = " << attribute_info.has_frames << "\n";
      InterpolateMaterial(key_frame, attribute_info.interpolate_mode, node,
                          &material_frame.diffuse_color);
      mesh->diffuse_colors[instance_id] = material_frame.diffuse_color;
      mesh->shade_modes[instance_id] = node->attribute_info.shade_mode;
      mesh->render_modes[instance_id] = node->attribute_info.render_mode;
      mesh->normal_display_modes[instance_id] =
          node->attribute_info.normal_display_mode;
    } break;
    case SceneNode::kLight: {
      // std::cout << "light\n";
      InterpolateLight(key_frame, attribute_info.interpolate_mode, node,
                       &node->light.color);
      if (node->light.light_type == SceneNode::Light::kLightPoint) {
        node->light.position = glm::vec3(glm::column(model, 3));
        // std::cout << "position = " << node->light.position << "\n";
      } else
        node->light.direction = glm::vec3(glm::column(model, 2));
      float type = static_cast<float>(node->light.light_type);
      scene_graph->light_info[node->light.light_id] =
          glm::mat4(glm::vec4(node->light.color, 0.0f),
                    glm::vec4(node->light.position, 1.0f),
                    glm::vec4(node->light.direction, 0.0f),
                    glm::vec4(node->light.attenuation, type));
    } break;
    case SceneNode::kTransform: {
      // std::cout << "transform\n";
      glm::vec4 rotate;
      glm::vec3 scale;
      glm::vec3 translate;
      InterpolateTransform(key_frame, attribute_info.interpolate_mode, node,
                           &rotate, &scale, &translate);
      // std::cout << "rotate = " << rotate << "\n";
      // std::cout << "scale = " << scale << "\n";
      // std::cout << "translate = " << translate << "\n";
      // std::cout << "angle = " << rotate[3] << " axis = " <<
      // glm::vec3(rotate)
      //<< "\n";

      glm::mat4 transform =
          glm::translate(translate) *
          glm::rotate(glm::radians(rotate[3]), glm::vec3(rotate)) *
          glm::scale(scale);
      node->model = node->model * transform;
    } break;
    default:
      break;
  }

  // std::cout << "UpdateSceneNode: num_children = " << node->num_children <<
  //"\n";
  // std::cout << "UpdateSceneNode: node->model = " << node->model << "\n";
  // std::cout << "UpdateSceneNode:: shade_mode = "
  //<< node->attribute_info.shade_mode << "\n";

  for (int i = 0; i < node->num_children; ++i)
    UpdateSceneNode(key_frame, node->model, node->attribute_info, scene_graph,
                    &node->children[i]);
}

void UpdateSceneGraph(int key_frame, SceneGraph* scene_graph) {
  glm::mat4 identity = glm::mat4(1.0f);
  SceneNode::Attribute default_attribute_info;
  UpdateSceneNode(key_frame, identity, default_attribute_info, scene_graph,
                  scene_graph->root);
}

void RenderSceneGraph(const std::vector<ProgramInfo*>& programs,
                      const SceneGraph* scene_graph) {
  // Get our camera.
  const SceneNode* camera_node =
      scene_graph->camera_nodes[scene_graph->current_camera];
  glm::mat4 camera_matrix = camera_node->model;
  glm::mat3 R_T = glm::transpose(glm::mat3(camera_matrix));
  glm::vec3 eye = glm::vec3(camera_matrix[3]);
  glm::mat4 view_matrix = glm::mat4(R_T);
  view_matrix[3] = glm::vec4(-R_T * eye, 1.0f);
  glm::mat4 projection_matrix = glm::mat4(1.0f);
  if (camera_node->camera.projection ==
      SceneNode::Camera::kProjectionPerspective) {
    projection_matrix =
        glm::perspective(camera_node->camera.fov, aspect,
                         camera_node->camera.near, camera_node->camera.far);
    //    std::cout << "Rendering with perspective\n";
  } else {
    float near = camera_node->camera.near;
    float far = camera_node->camera.far;
    float height = camera_node->camera.height;
    float width = aspect * height;
    projection_matrix = glm::ortho(-width / 2.0f, width / 2.0f, -height / 2.0f,
                                   height / 2.0f, near, far);
  }
  //  std::cout << "view = " << view_matrix << "\n";
  //  std::cout << "projection = " << projection_matrix << "\n";
  for (int program_index = 0; program_index < programs.size();
       ++program_index) {
    const ProgramInfo* program = programs[program_index];
    GLuint program_id = program->program_id;
    GLint projection_matrix_location = program->projection_matrix_location;
    GLint view_matrix_location = program->view_matrix_location;
    GLint lights_location = program->lights_location;
    GLint normal_length_location = program->normal_length_location;
    CHECK_GL_ERROR(glUseProgram(program_id));
    // Use our program.
    for (int mesh_index = 0; mesh_index < scene_graph->meshes.size();
         ++mesh_index) {
      GLuint vao = scene_graph->array_objects[mesh_index];
      const Mesh* mesh = scene_graph->meshes[mesh_index];

      // std::cout << "Binding to vao = " << vao << " mesh_index = " <<
      // mesh_index
      //<< " program_id = " << program_id << "\n";
      // Switch to our object VAO.
      CHECK_GL_ERROR(glBindVertexArray(vao));

      // Copy instanced data.

      // Update diffuse colors.
      CHECK_GL_ERROR(glBindBuffer(
          GL_ARRAY_BUFFER,
          scene_graph->buffer_objects[mesh_index * SceneGraph::kNumBufferTypes +
                                      SceneGraph::kDiffuseColorBuffer]));
      CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                                  mesh->num_instances * sizeof(glm::vec3),
                                  &mesh->diffuse_colors[0], GL_STATIC_DRAW));

      // Update shade modes.
      CHECK_GL_ERROR(glBindBuffer(
          GL_ARRAY_BUFFER,
          scene_graph->buffer_objects[mesh_index * SceneGraph::kNumBufferTypes +
                                      SceneGraph::kShadeModeBuffer]));
      CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                                  mesh->num_instances * sizeof(int),
                                  &mesh->shade_modes[0], GL_STATIC_DRAW));

      // Update render modes.
      CHECK_GL_ERROR(glBindBuffer(
          GL_ARRAY_BUFFER,
          scene_graph->buffer_objects[mesh_index * SceneGraph::kNumBufferTypes +
                                      SceneGraph::kRenderModeBuffer]));
      CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                                  mesh->num_instances * sizeof(int),
                                  &mesh->render_modes[0], GL_STATIC_DRAW));

      // Update normal render modes.
      CHECK_GL_ERROR(glBindBuffer(
          GL_ARRAY_BUFFER,
          scene_graph->buffer_objects[mesh_index * SceneGraph::kNumBufferTypes +
                                      SceneGraph::kNormalDisplayModeBuffer]));
      CHECK_GL_ERROR(
          glBufferData(GL_ARRAY_BUFFER, mesh->num_instances * sizeof(int),
                       &mesh->normal_display_modes[0], GL_STATIC_DRAW));
      // for (int i = 0; i < mesh->num_instances; ++i) {
      // std::cout << "shade_mode " << i << " = " << mesh->shade_modes[i]
      //<< "\n";
      //}

      // Update model view matrices.
      CHECK_GL_ERROR(glBindBuffer(
          GL_ARRAY_BUFFER,
          scene_graph->buffer_objects[mesh_index * SceneGraph::kNumBufferTypes +
                                      SceneGraph::kModelMatrixBuffer]));
      CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                                  mesh->num_instances * sizeof(glm::mat4),
                                  &mesh->model_matrices[0], GL_STATIC_DRAW));

      // std::cout << "Passing " << mesh->model_matrices.size()
      //<< " model matrices "
      //<< " and we have " << mesh->num_instances
      //<< " instances.  These should match.\n";

      // for (int i = 0; i < mesh->model_matrices.size(); ++i) {
      // std::cout << "mesh->model_matrices[" << i
      //<< "] = " << mesh->model_matrices[i] << "\n";
      //}
      // Pass uniforms in.
      CHECK_GL_ERROR(glUniformMatrix4fv(projection_matrix_location, 1, GL_FALSE,
                                        &projection_matrix[0][0]));
      CHECK_GL_ERROR(glUniformMatrix4fv(view_matrix_location, 1, GL_FALSE,
                                        &view_matrix[0][0]));
      CHECK_GL_ERROR(glUniformMatrix4fv(lights_location, SceneGraph::kMaxLights,
                                        GL_FALSE,
                                        &scene_graph->light_info[0][0][0]));
      CHECK_GL_ERROR(
          glUniform1f(normal_length_location, current_normal_length));

      // Draw our triangles.
      CHECK_GL_ERROR(
          glDrawElementsInstanced(GL_TRIANGLES, mesh->faces.size() * 3,
                                  GL_UNSIGNED_INT, 0, mesh->num_instances));
    }
  }
}

template <typename T>
inline void ParseValue(const json_spirit::Value* value, T* result) {
  *result = static_cast<T>(static_cast<std::string>(value->get_str()));
}

template <>
inline void ParseValue<glm::vec3>(const json_spirit::Value* value,
                                  glm::vec3* result) {
  std::stringstream ss(value->get_str());
  char dummy;
  float a, b, c;
  ss >> dummy >> a >> dummy >> b >> dummy >> c >> dummy;
  *result = glm::vec3(a, b, c);
}

template <>
inline void ParseValue<glm::vec4>(const json_spirit::Value* value,
                                  glm::vec4* result) {
  std::stringstream ss(value->get_str());
  char dummy;
  float a, b, c, d;
  ss >> dummy >> a >> dummy >> b >> dummy >> c >> dummy >> d >> dummy;
  *result = glm::vec4(a, b, c, d);
}

template <>
inline void ParseValue<float>(const json_spirit::Value* value, float* result) {
  *result = std::stof(value->get_str());
}

template <>
inline void ParseValue<int>(const json_spirit::Value* value, int* result) {
  *result = std::stoi(value->get_str());
}

void LookupString(const std::string& value, const char** strings, int count,
                  int* result) {
  for (int i = 0; i < count; ++i)
    if (value == std::string(strings[i])) *result = i;
}

template <>
inline void ParseValue<SceneNode::Attribute::RenderMode>(
    const json_spirit::Value* value, SceneNode::Attribute::RenderMode* result) {
  std::string string_value = value->get_str();
  LookupString(string_value, SceneNode::Attribute::kRenderModeStringValues,
               SceneNode::Attribute::kNumRenderModes,
               reinterpret_cast<int*>(result));
}

template <>
inline void ParseValue<SceneNode::Attribute::ShadeMode>(
    const json_spirit::Value* value, SceneNode::Attribute::ShadeMode* result) {
  std::string string_value = value->get_str();
  LookupString(string_value, SceneNode::Attribute::kShadeModeStringValues,
               SceneNode::Attribute::kNumShadeModes,
               reinterpret_cast<int*>(result));
}

template <>
inline void ParseValue<SceneNode::Attribute::InterpolateMode>(
    const json_spirit::Value* value,
    SceneNode::Attribute::InterpolateMode* result) {
  std::string string_value = value->get_str();
  LookupString(string_value, SceneNode::Attribute::kInterpolateModeStringValues,
               SceneNode::Attribute::kNumInterpolateModes,
               reinterpret_cast<int*>(result));
}

template <>
inline void ParseValue<SceneNode::Attribute::NormalDisplayMode>(
    const json_spirit::Value* value,
    SceneNode::Attribute::NormalDisplayMode* result) {
  std::string string_value = value->get_str();
  LookupString(string_value,
               SceneNode::Attribute::kNormalDisplayModeStringValues,
               SceneNode::Attribute::kNumNormalDisplayModes,
               reinterpret_cast<int*>(result));
}

template <>
inline void ParseValue<SceneNode::Light::LightType>(
    const json_spirit::Value* value, SceneNode::Light::LightType* result) {
  std::string string_value = value->get_str();
  LookupString(string_value, SceneNode::Light::kLightTypeStringValues,
               SceneNode::Light::kNumLightTypes,
               reinterpret_cast<int*>(result));
}

template <>
inline void ParseValue<SceneNode::Camera::ProjectionType>(
    const json_spirit::Value* value,
    SceneNode::Camera::ProjectionType* result) {
  std::string string_value = value->get_str();
  LookupString(string_value, SceneNode::Camera::kProjectionTypeStringValues,
               SceneNode::Camera::kNumProjectionTypes,
               reinterpret_cast<int*>(result));
}

template <>
inline void ParseValue<TransformFrames>(const json_spirit::Value* value,
                                        TransformFrames* result) {
  const json_spirit::Array& frame_array = value->get_array();
  //  TransformFrame* frames = new TransformFrame[frame_array.size()];
  TransformFrame* frames = Alloc<TransformFrame>(frame_array.size());
  result->frames = frames;
  result->num_frames = frame_array.size();
  for (uint32_t i = 0; i < frame_array.size(); ++i) {
    TransformFrame* frame = &frames[i];
    frame->translate = glm::vec3(0.0f);
    frame->rotate = glm::vec4(0.0f, 1.0f, 0.0f, 0.0f);
    frame->scale = glm::vec3(1.0f);
    frame->t = -1;
    const json_spirit::Object& object = frame_array[i].get_obj();
    std::unordered_map<std::string, const json_spirit::Value*> values;
    for (json_spirit::Object::size_type i = 0; i != object.size(); ++i) {
      const json_spirit::Pair& pair = object[i];
      const std::string& name = pair.name_;
      const json_spirit::Value& value = pair.value_;
      values[name] = &value;
    }
    if (values.find("translate") != values.end()) {
      // std::cout << "translate = " << values.at("translate")->get_str() <<
      // "\n";
      ParseValue<glm::vec3>(values.at("translate"), &frame->translate);
    }
    if (values.find("rotate") != values.end()) {
      ParseValue<glm::vec4>(values.at("rotate"), &frame->rotate);
    }
    if (values.find("scale") != values.end()) {
      ParseValue<glm::vec3>(values.at("scale"), &frame->scale);
    }
    if (values.find("t") != values.end()) {
      frame->t = std::stoi(values.at("t")->get_str());
    }
    if (frame->t == -1) {
      if (i == 0)
        frame->t = 0;
      else
        frame->t = frames[i - 1].t + 1;
    }
    max_key_frame = std::max(frame->t, max_key_frame);
  }
}

template <>
inline void ParseValue<MaterialFrames>(const json_spirit::Value* value,
                                       MaterialFrames* result) {
  const json_spirit::Array& frame_array = value->get_array();
  // MaterialFrame* frames = new MaterialFrame[frame_array.size()];
  MaterialFrame* frames = Alloc<MaterialFrame>(frame_array.size());
  result->frames = frames;
  result->num_frames = frame_array.size();
  for (uint32_t i = 0; i < frame_array.size(); ++i) {
    MaterialFrame* frame = &frames[i];
    frame->diffuse_color = glm::vec3(0.5f);
    frame->t = -1;
    const json_spirit::Object& object = frame_array[i].get_obj();
    std::unordered_map<std::string, const json_spirit::Value*> values;
    for (json_spirit::Object::size_type i = 0; i != object.size(); ++i) {
      const json_spirit::Pair& pair = object[i];
      const std::string& name = pair.name_;
      const json_spirit::Value& value = pair.value_;
      values[name] = &value;
    }
    if (values.find("diffuse_color") != values.end()) {
      // std::cout << "diffuse_color = " <<
      // values.at("diffuse_color")->get_str()
      //<< "\n";
      ParseValue<glm::vec3>(values.at("diffuse_color"), &frame->diffuse_color);
    }
    if (values.find("t") != values.end()) {
      frame->t = std::stoi(values.at("t")->get_str());
    }
    if (frame->t == -1) {
      if (i == 0)
        frame->t = 0;
      else
        frame->t = frames[i - 1].t + 1;
    }
    max_key_frame = std::max(frame->t, max_key_frame);
  }
}

template <>
inline void ParseValue<LightFrames>(const json_spirit::Value* value,
                                    LightFrames* result) {
  const json_spirit::Array& frame_array = value->get_array();
  // LightFrame* frames = new LightFrame[frame_array.size()];
  LightFrame* frames = Alloc<LightFrame>(frame_array.size());
  result->frames = frames;
  result->num_frames = frame_array.size();
  for (uint32_t i = 0; i < frame_array.size(); ++i) {
    LightFrame* frame = &frames[i];
    frame->color = glm::vec3(0.5f);
    frame->t = -1;
    const json_spirit::Object& object = frame_array[i].get_obj();
    std::unordered_map<std::string, const json_spirit::Value*> values;
    for (json_spirit::Object::size_type i = 0; i != object.size(); ++i) {
      const json_spirit::Pair& pair = object[i];
      const std::string& name = pair.name_;
      const json_spirit::Value& value = pair.value_;
      values[name] = &value;
    }
    if (values.find("color") != values.end()) {
      // std::cout << "color = " << values.at("color")->get_str() << "\n";
      ParseValue<glm::vec3>(values.at("color"), &frame->color);
    }
    if (values.find("t") != values.end()) {
      frame->t = std::stoi(values.at("t")->get_str());
    }
    if (frame->t == -1) {
      if (i == 0)
        frame->t = 0;
      else
        frame->t = frames[i - 1].t + 1;
    }
    max_key_frame = std::max(frame->t, max_key_frame);
  }
}

std::ostream& operator<<(std::ostream& os,
                         const SceneNode::Attribute::ShadeMode& mode) {
  std::string value = SceneNode::Attribute::kShadeModeStringValues[mode];
  os << value;
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const SceneNode::Attribute::RenderMode& mode) {
  std::string value = SceneNode::Attribute::kRenderModeStringValues[mode];
  os << value;
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const SceneNode::Attribute::InterpolateMode& mode) {
  std::string value = SceneNode::Attribute::kInterpolateModeStringValues[mode];
  os << value;
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const SceneNode::Attribute::NormalDisplayMode& mode) {
  std::string value =
      SceneNode::Attribute::kNormalDisplayModeStringValues[mode];
  os << value;
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const SceneNode::Light::LightType& mode) {
  std::string value = SceneNode::Light::kLightTypeStringValues[mode];
  os << value;
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const SceneNode::Camera::ProjectionType& mode) {
  std::string value = SceneNode::Camera::kProjectionTypeStringValues[mode];
  os << value;
  return os;
}

template <typename T>
inline bool LoadValue(
    const std::string& key,
    const std::unordered_map<std::string, const json_spirit::Value*>& values,
    T default_value, T* result) {
  *result = default_value;
  bool loaded = false;
  if (values.find(key) != values.end()) {
    ParseValue<T>(values.at(key), result);
    loaded = true;
  }
  // std::cout << "LoadValue: " << key << " -> " << *result << "\n";
  return loaded;
}

void LoadSceneNode(const json_spirit::Object& object, SceneNode* parent,
                   SceneNode* node, SceneGraph* scene_graph) {
  node->parent = parent;
  node->num_children = 0;
  node->children = nullptr;
  std::unordered_map<std::string, const json_spirit::Value*> values;
  for (json_spirit::Object::size_type i = 0; i != object.size(); ++i) {
    const json_spirit::Pair& pair = object[i];
    const std::string& name = pair.name_;
    const json_spirit::Value& value = pair.value_;
    values[name] = &value;
  }
  std::string type;
  if (values.find("type") != values.end()) {
    type = values["type"]->get_str();
    // std::cout << "Got a type: " << type << "\n";
    if (type == "attribute") {
      node->type = SceneNode::kAttribute;
      node->attribute.has_shade_mode =
          LoadValue<SceneNode::Attribute::ShadeMode>(
              "shade_mode", values, SceneNode::Attribute::kShadeFlat,
              &node->attribute.shade_mode);

      node->attribute.has_render_mode =
          LoadValue<SceneNode::Attribute::RenderMode>(
              "render_mode", values, SceneNode::Attribute::kRenderModeSolid,
              &node->attribute.render_mode);

      node->attribute.has_interpolate_mode =
          LoadValue<SceneNode::Attribute::InterpolateMode>(
              "interpolate_mode", values,
              SceneNode::Attribute::kInterpolateLinear,
              &node->attribute.interpolate_mode);

      node->attribute.has_normal_display_mode =
          LoadValue<SceneNode::Attribute::NormalDisplayMode>(
              "normal_display_mode", values,
              SceneNode::Attribute::kNormalDisplayNone,
              &node->attribute.normal_display_mode);

      MaterialFrames default_frames = {nullptr, 0};
      node->attribute.has_frames = LoadValue<MaterialFrames>(
          "frames", values, default_frames, &node->attribute.frames);
      // std::cout << "node->attribute.frames = " << node->attribute.frames
      //<< "\n";
    } else if (type == "camera") {
      node->type = SceneNode::kCamera;
      LoadValue<SceneNode::Camera::ProjectionType>(
          "projection", values,
          SceneNode::Camera::ProjectionType::kProjectionPerspective,
          &node->camera.projection);
      LoadValue<float>("fov", values, 45.0f, &node->camera.fov);
      LoadValue<float>("height", values, 1.0f, &node->camera.height);
      LoadValue<float>("near", values, 0.001f, &node->camera.near);
      LoadValue<float>("far", values, 0.001f, &node->camera.far);
      bool has_id = LoadValue<int>("id", values, -1, &node->camera.camera_id);
      if (node->camera.camera_id != -1) {
        // std::cout << " node->camera.camera_id = " << node->camera.camera_id
        //<< "\n";
        // std::cout << " camera_nodes.size() = "
        //<< scene_graph->camera_nodes.size() << "\n";
        scene_graph->camera_nodes[node->camera.camera_id] = node;
      } else {
        std::cerr << "No camera id was given!\n";
        exit(0);
      }
      if (scene_graph->current_camera == -1)
        scene_graph->current_camera = node->camera.camera_id;
    } else if (type == "geometry") {
      node->type = SceneNode::kGeometry;
      node->geometry.mesh = nullptr;
      node->geometry.mesh_id = -1;
      node->geometry.instance_id = -1;
      std::string file_name = "";
      LoadValue<std::string>("mesh", values, "", &file_name);
      Mesh* mesh = nullptr;
      if (scene_graph->mesh_lookup.find(file_name) ==
          scene_graph->mesh_lookup.end()) {
        // mesh = new Mesh;
        mesh = Alloc<Mesh>(1);
        // std::cout << "allocated mesh = " << mesh << "\n";
        LoadObj(file_name, mesh->vertices, mesh->faces);
        mesh->id = scene_graph->meshes.size();
        mesh->num_instances = 0;
        // std::cout << "Loaded mesh '" << file_name << "' with  "
        //<< mesh->vertices.size() << " vertices and "
        //<< mesh->faces.size() << " faces.\n";
        mesh->bounds.min = glm::vec3(std::numeric_limits<float>::max());
        mesh->bounds.max = glm::vec3(-std::numeric_limits<float>::max());
        for (int i = 0; i < mesh->vertices.size(); ++i) {
          mesh->bounds.min =
              glm::min(glm::vec3(mesh->vertices[i]), mesh->bounds.min);
          mesh->bounds.max =
              glm::max(glm::vec3(mesh->vertices[i]), mesh->bounds.max);
        }
        mesh->face_normals.resize(mesh->faces.size(), glm::vec4(0.0f));
        mesh->vertex_normals.resize(mesh->vertices.size(), glm::vec4(0.0f));
        for (int i = 0; i < mesh->faces.size(); ++i) {
          glm::uvec3 ijk = mesh->faces[i];
          glm::vec3 a = glm::vec3(mesh->vertices[ijk[0]]);
          glm::vec3 b = glm::vec3(mesh->vertices[ijk[1]]);
          glm::vec3 c = glm::vec3(mesh->vertices[ijk[2]]);
          glm::vec3 u = b - a;
          glm::vec3 v = c - a;
          glm::vec3 n = glm::cross(u, v);
          mesh->face_normals[i] = glm::vec4(n, 0.0f);
          mesh->vertex_normals[ijk[0]] += glm::vec4(n, 0.0f);
          mesh->vertex_normals[ijk[1]] += glm::vec4(n, 0.0f);
          mesh->vertex_normals[ijk[2]] += glm::vec4(n, 0.0f);
        }
        for (int i = 0; i < mesh->vertex_normals.size(); ++i)
          mesh->vertex_normals[i] = glm::normalize(mesh->vertex_normals[i]);
        // std::cout << "'" << file_name
        //<< "':mesh->bounds.min = " << mesh->bounds.min << "\n";
        // std::cout << "'" << file_name
        //<< "':mesh->bounds.max = " << mesh->bounds.max << "\n";
        scene_graph->meshes.push_back(mesh);
        scene_graph->mesh_lookup[file_name] = mesh;
      } else {
        // std::cout << "Reusing instance of '" << file_name << "'\n";
      }
      mesh = scene_graph->mesh_lookup.at(file_name);
      ++mesh->num_instances;
      mesh->model_matrices.push_back(glm::mat4(1.0f));
      mesh->diffuse_colors.push_back(glm::vec3(1.0f));
      mesh->shade_modes.push_back(SceneNode::Attribute::ShadeMode::kShadeFlat);
      mesh->render_modes.push_back(
          SceneNode::Attribute::RenderMode::kRenderModeSolid);
      mesh->normal_display_modes.push_back(
          SceneNode::Attribute::NormalDisplayMode::kNormalDisplayNone);
      std::string* mesh_name = Alloc<std::string>(1);
      *mesh_name = file_name;
      node->geometry.mesh = mesh_name->c_str();
      node->geometry.mesh_id = mesh->id;
      node->geometry.instance_id = mesh->num_instances - 1;
      scene_graph->mesh_nodes.push_back(node);
      // std::cout << "mesh = " << mesh << "\n";
      // std::cout << "file_name = " << file_name << "\n";
      // std::cout << "node->geometry.mesh_id = " << node->geometry.mesh_id
      //<< "\n";
      // std::cout << "node->geometry.mesh = " << node->geometry.mesh << "\n";
    } else if (type == "light") {
      if (scene_graph->num_lights >= SceneGraph::kMaxLights) return;
      ++scene_graph->num_lights;
      node->type = SceneNode::kLight;
      LoadValue<SceneNode::Light::LightType>("light_type", values,
                                             SceneNode::Light::kLightPoint,
                                             &node->light.light_type);
      LoadValue<glm::vec3>("attenuation", values, glm::vec3(1.0f, 0.0f, 0.0f),
                           &node->light.attenuation);
      LightFrames default_frames = {nullptr, 0};
      LoadValue<LightFrames>("frames", values, default_frames,
                             &node->light.frames);
      node->light.light_id = scene_graph->num_lights;
      scene_graph->light_nodes.push_back(node);
      scene_graph->light_info.push_back(glm::mat4(0.0f));
      // std::cout << "node->light.frames = " << node->light.frames << "\n";
    } else if (type == "object") {
      node->type = SceneNode::kObject;
      // Do nothing here.
    } else if (type == "transform") {
      node->type = SceneNode::kTransform;
      TransformFrames default_frames = {nullptr, 0};
      LoadValue<TransformFrames>("frames", values, default_frames,
                                 &node->transform.frames);
      // std::cout << "node->transform.frames = " << node->transform.frames
      //<< "\n";
    } else {
      // std::cerr << "Invalid type field found:" << type << "!\n";
      return;
    }
    if (values.find("children") != values.end()) {
      const json_spirit::Array& child_array =
          values.at("children")->get_array();
      node->num_children = child_array.size();
      //      node->children = new SceneNode[node->num_children];
      node->children = Alloc<SceneNode>(node->num_children);
      for (int i = 0; i < node->num_children; ++i)
        LoadSceneNode(child_array[i].get_obj(), node, &node->children[i],
                      scene_graph);
    }
  } else {
    std::cerr << "No type field found!\n";
    exit(0);
  }
}

void LoadSceneGraph(const std::string& file_name, SceneGraph* scene_graph) {
  //  SceneNode* root = new SceneNode;
  SceneNode* root = Alloc<SceneNode>(1);
  root->type = SceneNode::kObject;
  std::ifstream is(file_name);
  json_spirit::Value value;
  json_spirit::read_stream(is, value);
  const json_spirit::Array& node_array = value.get_array();
  root->parent = nullptr;
  root->num_children = node_array.size();
    root->children = new SceneNode[root->num_children];
  root->children = Alloc<SceneNode>(root->num_children);
   std::cout << "Loading " << root->num_children << " children.\n";
  for (uint32_t i = 0; i < node_array.size(); ++i)
    LoadSceneNode(node_array[i].get_obj(), root, &root->children[i],
                  scene_graph);
  scene_graph->root = root;
}

}  // namespace scene_graph

float last_x = 0.0f, last_y = 0.0f, current_x = 0.0f, current_y = 0.0f;
bool drag_state = false;
int current_button = -1;

void ErrorCallback(int error, const char* description) {
  std::cerr << "GLFW Error: " << description << "\n";
}

scene_graph::SceneGraph scene_graph_instance;
void KeyCallback(GLFWwindow* window, int key, int scancode, int action,
                 int mods) {
  int last_camera = scene_graph_instance.current_camera;
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GL_TRUE);
  else if (key == GLFW_KEY_0 && action == GLFW_PRESS)
    scene_graph_instance.current_camera = 0;
  else if (key == GLFW_KEY_1 && action == GLFW_PRESS)
    scene_graph_instance.current_camera = 1;
  else if (key == GLFW_KEY_2 && action == GLFW_PRESS)
    scene_graph_instance.current_camera = 2;
  else if (key == GLFW_KEY_3 && action == GLFW_PRESS)
    scene_graph_instance.current_camera = 3;
  else if (key == GLFW_KEY_4 && action == GLFW_PRESS)
    scene_graph_instance.current_camera = 4;
  else if (key == GLFW_KEY_5 && action == GLFW_PRESS)
    scene_graph_instance.current_camera = 5;
  else if (key == GLFW_KEY_6 && action == GLFW_PRESS)
    scene_graph_instance.current_camera = 6;
  else if (key == GLFW_KEY_7 && action == GLFW_PRESS)
    scene_graph_instance.current_camera = 7;
  else if (key == GLFW_KEY_8 && action == GLFW_PRESS)
    scene_graph_instance.current_camera = 8;
  else if (key == GLFW_KEY_9 && action == GLFW_PRESS)
    scene_graph_instance.current_camera = 9;
  else if (key == GLFW_KEY_LEFT &&
           (action == GLFW_PRESS || action == GLFW_REPEAT))
    current_key_frame -= 4;
  else if (key == GLFW_KEY_RIGHT &&
           (action == GLFW_PRESS || action == GLFW_REPEAT))
    current_key_frame += 4;
  else if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
    pause_enabled = !pause_enabled;
  else if (key == GLFW_KEY_UP && action == GLFW_PRESS)
    current_normal_length *= 2.0f;
  else if (key == GLFW_KEY_DOWN && action == GLFW_PRESS)
    current_normal_length /= 2.0f;
  else if (key == GLFW_KEY_LEFT_BRACKET && action == GLFW_PRESS) {
    do {
      scene_graph_instance.current_camera =
          (scene_graph::SceneGraph::kMaxCameras +
           scene_graph_instance.current_camera - 1) %
          scene_graph::SceneGraph::kMaxCameras;
    } while (scene_graph_instance.camera_nodes
                 [scene_graph_instance.current_camera] == nullptr);
  } else if (key == GLFW_KEY_RIGHT_BRACKET && action == GLFW_PRESS) {
    do {
      scene_graph_instance.current_camera =
          (scene_graph::SceneGraph::kMaxCameras +
           scene_graph_instance.current_camera + 1) %
          scene_graph::SceneGraph::kMaxCameras;
    } while (scene_graph_instance.camera_nodes
                 [scene_graph_instance.current_camera] == nullptr);
  }

  if (current_normal_length < 0.005f) current_normal_length = 0.005f;
  if (current_key_frame < 0) current_key_frame = 0;
  if (current_key_frame > max_key_frame + 3)
    current_key_frame = max_key_frame + 3;

  if (scene_graph_instance.camera_nodes[scene_graph_instance.current_camera] ==
      nullptr)
    scene_graph_instance.current_camera = last_camera;
}

void MousePosCallback(GLFWwindow* window, double mouse_x, double mouse_y) {
  last_x = current_x;
  last_y = current_y;
  current_x = mouse_x;
  current_y = mouse_y;
  float delta_x = current_x - last_x;
  float delta_y = current_y - last_y;
  if (sqrt(delta_x * delta_x + delta_y * delta_y) < 1e-15) return;
}

void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
  drag_state = (action == GLFW_PRESS);
  current_button = button;
}

int main(int argc, char* argv[]) {
  std::string file_name(argv[1]);
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

  // Load our scene.
  scene_graph::LoadSceneGraph(file_name, &scene_graph_instance);

  // OK, we have a scene and we have loaded all of the meshes,
  // now setup all of our VAOs, VBOs, and all that jazz.
  // We are not compiling and linking programs...yet.
  scene_graph::InitSceneGraphGL(&scene_graph_instance);

  // Setup programs here.
  std::vector<scene_graph::ProgramInfo*> programs;

  // Compile and link the program that renders solid shaded meshes.
  scene_graph::ProgramInfo solid_program_info;
  scene_graph::InitProgramGL(solid_vertex_shader, solid_fragment_shader,
                             solid_geometry_shader, &solid_program_info);
  programs.push_back(&solid_program_info);

  scene_graph::ProgramInfo wire_program_info;
  scene_graph::InitProgramGL(wire_vertex_shader, wire_fragment_shader,
                             wire_geometry_shader, &wire_program_info);
  programs.push_back(&wire_program_info);

  scene_graph::ProgramInfo normal_program_info;
  scene_graph::InitProgramGL(normal_vertex_shader, normal_fragment_shader,
                             normal_geometry_shader, &normal_program_info);
  programs.push_back(&normal_program_info);

  glfwSwapInterval(1);
  float aspect = 0.0f;
  current_key_frame = 0;
  while (!glfwWindowShouldClose(window)) {
    // Setup some basic window stuff.
    glfwGetFramebufferSize(window, &window_width, &window_height);
    glViewport(0, 0, window_width, window_height);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDepthFunc(GL_LESS);

    aspect = static_cast<float>(window_width) / window_height;

    // Update the scene graph.
    scene_graph::UpdateSceneGraph(current_key_frame, &scene_graph_instance);

    // Render!
    scene_graph::RenderSceneGraph(programs, &scene_graph_instance);
    if (!pause_enabled) ++current_key_frame;
    if (current_key_frame > max_key_frame + 3) {
      current_key_frame = max_key_frame + 3;
      pause_enabled = true;
    }
    // Poll and swap.
    glfwPollEvents();
    glfwSwapBuffers(window);
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
