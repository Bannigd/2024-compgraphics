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
  Vec3 worldCoords;
  Vec3 viewCoords;
  Vertex() {}
  Vertex(Vec3 coords) { this->worldCoords = coords; }
  Vertex(float x, float y, float z) { this->worldCoords = Vec3(x, y, z); }

  void world2View(const Viewpoint &Camera) {
    float x_v = -worldCoords.x * std::sin(Camera.theta) +
                worldCoords.y * std::cos(Camera.theta);
    float y_v = -worldCoords.x * std::cos(Camera.phi) * std::cos(Camera.theta) -
                worldCoords.y * std::cos(Camera.phi) * std::sin(Camera.theta) +
                worldCoords.z * std::sin(Camera.phi);
    float z_v = Camera.rho -
                worldCoords.x * std::sin(Camera.phi) * std::cos(Camera.theta) -
                worldCoords.y * std::sin(Camera.phi) * std::sin(Camera.theta) -
                worldCoords.z * std::cos(Camera.phi);
    this->viewCoords = Vec3(x_v, y_v, z_v);
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

class Triangle {
public:
  int idxA, idxB, idxC;
  // Vertex a;
  // Vertex b;
  // Vertex c;
  bool onScreen; // true if triangle is before viewpoint
  Triangle(int a, int b, int c) {
    this->idxA = a - 1;
    this->idxB = b - 1;
    this->idxC = c - 1;
  }
  friend std::ostream &operator<<(std::ostream &out, const Triangle &Triangle);
};

class WFModel {
public:
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  std::vector<Triangle> triangles;
  WFModel() {}
  WFModel(std::vector<Vertex> vertices, std::vector<Edge> edges,
          std::vector<Triangle> triangles) {
    this->vertices = vertices;
    this->edges = edges;
    this->triangles = triangles;
  };

  void computeSphereCoords(const Viewpoint &Camera) {
    for (auto &Ver : this->vertices) {
      Ver.world2View(Camera);
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
  void changeView(const Vec3 &direction) {
    this->Camera->move(direction);
    Figure->computeSphereCoords(*Camera);
  }

  void changeView(const float &rho, const float &theta, const float &phi) {
    this->Camera->move(rho, theta, phi);
    Figure->computeSphereCoords(*Camera);
  }

  void draw(const int &xShift, const int &yShift) {
    float d = 300.;

    // draw vertices
    for (auto &vertex : Figure->vertices) {
      DrawCircle(d * vertex.viewCoords.x / vertex.viewCoords.z + xShift,
                 -d * vertex.viewCoords.y / vertex.viewCoords.z + yShift, 5,
                 BLUE);
    }

    // draw edges
    // for (auto &edge : Figure->edges) {
    //   DrawLine(d * Figure->vertices[edge.startIdx].viewCoords.x /
    //                    Figure->vertices[edge.startIdx].viewCoords.z +
    //                xShift,
    //            -d * Figure->vertices[edge.startIdx].viewCoords.y /
    //                    Figure->vertices[edge.startIdx].viewCoords.z +
    //                yShift,
    //            d * Figure->vertices[edge.endIdx].viewCoords.x /
    //                    Figure->vertices[edge.endIdx].viewCoords.z +
    //                xShift,
    //            -d * Figure->vertices[edge.endIdx].viewCoords.y /
    //                    Figure->vertices[edge.endIdx].viewCoords.z +
    //                yShift,
    //            BLUE);
    // }
    // draw triangles
    for (auto &trig : Figure->triangles) {
      Vertex A = Figure->vertices[trig.idxA];
      Vertex B = Figure->vertices[trig.idxB];
      Vertex C = Figure->vertices[trig.idxC];

      DrawLine(d * A.viewCoords.x / A.viewCoords.z + xShift,
               -d * A.viewCoords.y / A.viewCoords.z + yShift,
               d * B.viewCoords.x / B.viewCoords.z + xShift,
               -d * B.viewCoords.y / B.viewCoords.z + yShift, BLUE);
      DrawLine(d * A.viewCoords.x / A.viewCoords.z + xShift,
               -d * A.viewCoords.y / A.viewCoords.z + yShift,
               d * C.viewCoords.x / C.viewCoords.z + xShift,
               -d * C.viewCoords.y / C.viewCoords.z + yShift, BLUE);
      DrawLine(d * B.viewCoords.x / B.viewCoords.z + xShift,
               -d * B.viewCoords.y / B.viewCoords.z + yShift,
               d * C.viewCoords.x / C.viewCoords.z + xShift,
               -d * C.viewCoords.y / C.viewCoords.z + yShift, BLUE);
    }
  }
  friend std::ostream &operator<<(std::ostream &out, const Scene &Scene);
};

inline std::ostream &operator<<(std::ostream &out, const Vec3 &V3) {
  std::cout << "(" << V3.x << "," << V3.y << "," << V3.z << ")";
  return out;
}

inline std::ostream &operator<<(std::ostream &out, const Vertex &Ver) {
  std::cout << "World: " << Ver.worldCoords << " Camera: " << Ver.viewCoords
            << "\n";
  return out;
}

inline std::ostream &operator<<(std::ostream &out, const Edge &Edge) {
  std::cout << "Start: " << Edge.startIdx << " End: " << Edge.endIdx;
  return out;
}

inline std::ostream &operator<<(std::ostream &out, const Triangle &Triangle) {
  std::cout << "A: " << Triangle.idxA << " B: " << Triangle.idxB
            << " C: " << Triangle.idxC << '\n';
  return out;
}

inline std::ostream &operator<<(std::ostream &out, const Viewpoint &Camera) {
  std::cout << "\tCamera position: (" << Camera.rho << ", " << Camera.theta
            << ", " << Camera.phi << ")\n";
  return out;
}

inline std::ostream &operator<<(std::ostream &out, const WFModel &Figure) {
  std::cout << "Number of vertices = " << Figure.vertices.size() << "\n";
  for (int i = 0; i < Figure.vertices.size(); i++) {
    std::cout << "Vertex #" << i << " " << Figure.vertices[i];
  }
  return out;
}

inline std::ifstream &operator>>(std::ifstream &f, WFModel &Figure) {
  int nVertices, nEdges;
  std::string buffer;

  // read verices
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

  // read edges
  std::getline(f, buffer);
  nEdges = std::stoi(buffer);
  std::vector<Edge> edges(nEdges);

  for (int i = 0; i < nEdges; i++) {
    std::getline(f, buffer);
    std::stringstream line(buffer);
    std::vector<std::string> cell(2);
    for (int j = 0; j < 2; j++) {
      std::getline(line, cell[j], ',');
    }
    // -1 shift to align with vector index
    edges[i] = Edge(std::stoi(cell[0]) - 1, std::stoi(cell[1]) - 1);
  }
  Figure.edges = edges;

  // read triangles
  std::vector<Triangle> triangles;
  while (std::getline(f, buffer)) {
    // std::getline(f, buffer);
    std::stringstream line(buffer);
    std::vector<std::string> cell(3);
    for (int j = 0; j < 3; j++) {
      std::getline(line, cell[j], ',');
    }
    triangles.push_back(
        Triangle(std::stoi(cell[0]), std::stoi(cell[1]), std::stoi(cell[2])));
  }
  Figure.triangles = triangles;
  return f;
}
inline std::ostream &operator<<(std::ostream &out, const Scene &Scene) {
  std::cout << "Scene:\n";
  std::cout << *Scene.Camera;
  std::cout << *Scene.Figure;
  return out;
}
} // namespace graph
