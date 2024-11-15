#include "raylib.h"
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace graph {
struct Vec3 {
  float x, y, z;
  Vec3() {}
  Vec3(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }
  friend std::ostream &operator<<(std::ostream &out, const graph::Vec3 &V3);
  Vec3 &operator+(const Vec3 &rhs) {
    this->x += rhs.x;
    this->y += rhs.y;
    this->z += rhs.z;
    return *this;
  }
  Vec3 &operator-(const Vec3 &rhs) {
    this->x -= rhs.x;
    this->y -= rhs.y;
    this->z -= rhs.z;
    return *this;
  }
};

struct Viewpoint {
  Vec3 coords;
  float rho, theta, phi;
  Viewpoint() {}
  Viewpoint(Vec3 pos) {
    coords = pos;
    this->updateSphere();
  }

  void updateSphere() {
    this->rho = std::sqrt(std::pow(coords.x, 2) + std::pow(coords.y, 2) +
                          std::pow(coords.z, 2));
    if (this->rho < 10 * FLT_EPSILON) {
      this->rho = 0;
      this->theta = 0;
      this->phi = 0;
      return;
    }
    if (coords.x > 0) {
      this->theta = std::atan(coords.y / coords.x);
    } else if (coords.x < 0) {
      this->theta = M_PI + std::atan(coords.y / coords.x);
    } else if (coords.x < 10 * FLT_EPSILON && coords.y >= 0) {
      this->theta = M_PI / 2;
    } else if (coords.x < 10 * FLT_EPSILON && coords.y < 0) {
      this->theta = 3 * M_PI / 2;
    };

    this->phi = std::acos(coords.z / this->rho);
  }

  void updateCartesian() {
    coords = Vec3(rho * std::sin(phi) * std::cos(theta),
                  rho * std::sin(phi) * std::sin(theta), rho * std::cos(phi));
  }
  void move(const Vec3 &direction) {
    coords = coords + direction;
    this->updateSphere();
  }
  void move(const float &rho, const float &theta, const float &phi) {
    this->rho += rho;
    this->theta += theta;
    this->phi += phi;
    this->updateCartesian();
  }
  friend std::ostream &operator<<(std::ostream &out, const Viewpoint &Camera);
};

struct Vertex {
  Vec3 world_coords;
  Vec3 view_coords;
  Vertex() {}
  Vertex(Vec3 coords) { this->world_coords = coords; }
  Vertex(float x, float y, float z) { this->world_coords = Vec3(x, y, z); }

  void world_to_view(const Viewpoint &Camera) {
    float x_v = -world_coords.x * std::sin(Camera.theta) +
                world_coords.y * std::cos(Camera.theta);
    float y_v =
        -world_coords.x * std::cos(Camera.phi) * std::cos(Camera.theta) -
        world_coords.y * std::cos(Camera.phi) * std::sin(Camera.theta) +
        world_coords.z * std::sin(Camera.phi);
    float z_v = Camera.rho -
                world_coords.x * std::sin(Camera.phi) * std::cos(Camera.theta) -
                world_coords.y * std::sin(Camera.phi) * std::sin(Camera.theta) -
                world_coords.z * std::cos(Camera.phi);
    this->view_coords = Vec3(x_v, y_v, z_v);
  }
  friend std::ostream &operator<<(std::ostream &out, const Vertex &Ver);
};

struct Edge {
  int startIdx, endIdx;
  friend std::ostream &operator<<(std::ostream &out, const Edge &Edge);
  Edge() {}
  Edge(int start, int end) {
    this->startIdx = start;
    this->endIdx = end;
  }
};

class WFModel {
public:
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  WFModel() {}
  WFModel(std::vector<Vertex> vertices, std::vector<Edge> edges) {
    this->vertices = vertices;
    this->edges = edges;
  };

  void computeSphereCoords(const Viewpoint &Camera) {
    // std::cout << "Inside compute -> " << Camera << "^ \n";
    for (auto &Ver : this->vertices) {
      Ver.world_to_view(Camera);
    }
  }

  friend std::ostream &operator<<(std::ostream &out, const WFModel &Figure);
  friend std::ifstream &operator>>(std::ifstream &fin, const WFModel &Figure);
};

class Scene {
public:
  Viewpoint *Camera;
  WFModel *Figure;
  Scene(Viewpoint &View, WFModel &Figure) {
    this->Camera = &View;
    this->Figure = &Figure;
    Figure.computeSphereCoords(*Camera);
  }
  void change_view(const Vec3 &direction) {
    this->Camera->move(direction);
    Figure->computeSphereCoords(*Camera);
  }

  void change_view(const float &rho, const float &theta, const float &phi) {
    this->Camera->move(rho, theta, phi);
    Figure->computeSphereCoords(*Camera);
  }

  void draw(const int &xShift, const int &yShift) {
    float d = 300.;
    // float xShift = 500.;
    // float yShift = 300.;
    // draw vertices
    for (auto &vertex : Figure->vertices) {
      DrawCircle(d * vertex.view_coords.x / vertex.view_coords.z + xShift,
                 -d * vertex.view_coords.y / vertex.view_coords.z + yShift, 5,
                 BLUE);
    }
    // draw edges
    for (auto &edge : Figure->edges) {
      DrawLine(d * Figure->vertices[edge.startIdx].view_coords.x /
                       Figure->vertices[edge.startIdx].view_coords.z +
                   xShift,
               -d * Figure->vertices[edge.startIdx].view_coords.y /
                       Figure->vertices[edge.startIdx].view_coords.z +
                   yShift,
               d * Figure->vertices[edge.endIdx].view_coords.x /
                       Figure->vertices[edge.endIdx].view_coords.z +
                   xShift,
               -d * Figure->vertices[edge.endIdx].view_coords.y /
                       Figure->vertices[edge.endIdx].view_coords.z +
                   yShift,
               BLUE);
    }
  }
  friend std::ostream &operator<<(std::ostream &out, const Scene &Scene);
};

std::ostream &operator<<(std::ostream &out, const Vec3 &V3) {
  std::cout << "(" << V3.x << "," << V3.y << "," << V3.z << ")";
  return out;
}

std::ostream &operator<<(std::ostream &out, const Vertex &Ver) {
  std::cout << "World: " << Ver.world_coords << " Camera: " << Ver.view_coords
            << "\n";
  return out;
}

std::ostream &operator<<(std::ostream &out, const Edge &Edge) {
  std::cout << "Start: " << Edge.startIdx << " End: " << Edge.endIdx;
  return out;
}

std::ostream &operator<<(std::ostream &out, const Viewpoint &Camera) {
  std::cout << "\tCamera position: (" << Camera.rho << ", " << Camera.theta
            << ", " << Camera.phi << ")\n";
  return out;
}

std::ostream &operator<<(std::ostream &out, const WFModel &Figure) {
  std::cout << "Number of vertices = " << Figure.vertices.size() << "\n";
  for (int i = 0; i < Figure.vertices.size(); i++) {
    std::cout << "Vertex #" << i << " " << Figure.vertices[i];
  }
  return out;
}

std::ifstream &operator>>(std::ifstream &f, WFModel &Figure) {
  int nVertices, nEdges;
  std::string buffer;

  std::getline(f, buffer);
  nVertices = std::stoi(buffer);

  std::vector<Vertex> vertices(nVertices);

  for (int i = 0; i < nVertices; i++) {
    std::getline(f, buffer);
    std::stringstream line(buffer);
    std::vector<std::string> cell(3);
    for (int j = 0; j < 3; j++) {
      std::getline(line, cell[j], ',');
    }
    vertices[i] = Vertex(
        Vec3(std::stof(cell[0]), std::stof(cell[1]), std::stof(cell[2])));
  }

  Figure.vertices = vertices;

  std::getline(f, buffer);
  nEdges = std::stoi(buffer);

  std::vector<Edge> edges(nEdges);
  for (int i = 0; i < nEdges; i++) {
    std::getline(f, buffer);
    std::stringstream line(buffer);
    std::vector<std::string> cell(2);
    for (int j = 0; j < 2; j++) {
      std::getline(line, cell[j], ' ');
    }
    // -1 shift to align with vector index
    edges[i] = Edge(std::stoi(cell[0]) - 1, std::stoi(cell[1]) - 1);
  }
  Figure.edges = edges;
  return f;
}
std::ostream &operator<<(std::ostream &out, const Scene &Scene) {
  std::cout << "Scene:\n";
  std::cout << *Scene.Camera;
  std::cout << *Scene.Figure;
  return out;
}
} // namespace graph

int main(int argc, char *argv[]) {

  graph::WFModel Model;

  std::string filename;
  if (argc == 2) {
    filename = argv[1];
  } else {
    filename = "cube.dat";
  }
  std::ifstream f(filename);

  f >> Model;
  f.close();

  // std::vector<graph::Vertex> vertices = {
  //     graph::Vertex(1., 1., -1.),   graph::Vertex(-1., 1, -1),
  //     graph::Vertex(-1., -1., -1.), graph::Vertex(1., -1., -1.),
  //     graph::Vertex(1., 1., 1.),    graph::Vertex(-1., 1., 1.),
  //     graph::Vertex(-1., -1., 1.),  graph::Vertex(1., -1., 1.),
  // };
  // std::vector<graph::Edge> edges = {
  //     graph::Edge(0, 1), graph::Edge(1, 2), graph::Edge(2, 3), graph::Edge(3,
  //     0), graph::Edge(4, 5), graph::Edge(5, 6), graph::Edge(6, 7),
  //     graph::Edge(7, 4), graph::Edge(0, 4), graph::Edge(1, 5), graph::Edge(2,
  //     6), graph::Edge(3, 7),
  // };

  graph::Vec3 pos1(3., 3., 4.);
  graph::Viewpoint Camera1(pos1);

  graph::Vec3 pos2(10., 10., 10.);
  graph::Viewpoint Camera2(pos2);
  graph::Scene Scene(Camera2, Model);

  int width = 900, height = 700;
  InitWindow(width, height, "Screen projection");
  SetTargetFPS(60);

  while (!WindowShouldClose()) {
    BeginDrawing();
    {
      ClearBackground(BLACK);
      // Scene.change_view(graph::Vec3(-0.01, -0.05, 0.01));
      Scene.change_view(0, 0.01, 0);
      Scene.draw(width / 2, height / 2);
    }
    EndDrawing();
  }
  return 0;
}
