#include "raylib.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define MAX_PREC 1e-6

namespace graph {
struct Vec3 {
  double x, y, z;
  Vec3() {}
  Vec3(double x, double y, double z) {
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
  double rho, theta, phi;
  Viewpoint() {}
  Viewpoint(Vec3 pos) {
    coords = pos;
    this->updateSphere();
  }

  void updateSphere() {
    this->rho = std::sqrt(std::pow(coords.x, 2) + std::pow(coords.y, 2) +
                          std::pow(coords.z, 2));
    if (this->rho < 10 * MAX_PREC) {
      this->rho = 0;
      this->theta = 0;
      this->phi = 0;
      return;
    }
    if (coords.x > 0) {
      this->theta = std::atan(coords.y / coords.x);
    } else if (coords.x < 0) {
      this->theta = M_PI + std::atan(coords.y / coords.x);
    } else if (coords.x < MAX_PREC && coords.y >= 0) {
      this->theta = M_PI / 2;
    } else if (coords.x < MAX_PREC && coords.y < 0) {
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
  void move(const double &rho, const double &theta, const double &phi) {
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
  Vertex(double x, double y, double z) { this->worldCoords = Vec3(x, y, z); }

  void world2View(const Viewpoint &Camera) {
    double x_v = -worldCoords.x * std::sin(Camera.theta) +
                 worldCoords.y * std::cos(Camera.theta);
    double y_v =
        -worldCoords.x * std::cos(Camera.phi) * std::cos(Camera.theta) -
        worldCoords.y * std::cos(Camera.phi) * std::sin(Camera.theta) +
        worldCoords.z * std::sin(Camera.phi);
    double z_v = Camera.rho -
                 worldCoords.x * std::sin(Camera.phi) * std::cos(Camera.theta) -
                 worldCoords.y * std::sin(Camera.phi) * std::sin(Camera.theta) -
                 worldCoords.z * std::cos(Camera.phi);
    this->viewCoords = Vec3(x_v, y_v, z_v);
  }
  friend std::ostream &operator<<(std::ostream &out, const Vertex &Ver);
  void draw(int i, const int &d, const int &xShift, const int &yShift) {
    DrawCircle(d * viewCoords.x / viewCoords.z + xShift,
               -d * viewCoords.y / viewCoords.z + yShift, 5, RED);
    DrawText(std::to_string(i).c_str(),
             d * viewCoords.x / viewCoords.z + xShift + 5,
             -d * viewCoords.y / viewCoords.z + yShift + 5, 15, WHITE);
  }
};

struct Edge {
  Vertex &start, &end;
  friend std::ostream &operator<<(std::ostream &out, const Edge &Edge);
  Edge(Vertex &st, Vertex &en) : start(st), end(en) {
    this->start = st;
    this->end = en;
  }
  Edge& operator=(const Edge& e) {
    this->start = e.start;
    this->end = e.end;
    return *this;
  }
  void draw(const int &d, const int &xShift, const int &yShift) {
    DrawLine(d * start.viewCoords.x / start.viewCoords.z + xShift,
             -d * start.viewCoords.y / start.viewCoords.z + yShift,
             d * end.viewCoords.x / end.viewCoords.z + xShift,
             -d * end.viewCoords.y / end.viewCoords.z + yShift, BLUE);
  }
};

class Triangle {
public:
  Vertex &A;
  Vertex &B;
  Vertex &C;
  Vec3 n;
  double a, b, c, h;
  bool onScreen; // true if triangle is before viewpoint
  Triangle(Vertex &A, Vertex &B, Vertex &C) : A(A), B(B), C(C) {
    this->A = A;
    this->B = B;
    this->C = C;
  }

  Triangle& operator=(const Triangle& t) {
    this->A = t.A;
    this->B = t.B;
    this->C = t.C;
    return *this;
  }

  void computePlane() {
    double xA = A.viewCoords.x;
    double yA = A.viewCoords.y;
    double zA = A.viewCoords.z;
    double xB = B.viewCoords.x;
    double yB = B.viewCoords.y;
    double zB = B.viewCoords.z;
    double xC = C.viewCoords.x;
    double yC = C.viewCoords.y;
    double zC = C.viewCoords.z;
    this->a = yA * (zB - zC) - yB * (zA - zC) + yC * (zA - zB);
    this->b = -(xA * (zB - zC) - xB * (zA - zC) + xC * (zA - zB));
    this->c = xA * (yB - yC) - xB * (yA - yC) + xC * (yA - yB);
    this->h = xA * (yB * zC - yC * zB) - xB * (yA * zC - yC * zA) +
              xC * (yA * zB - yB * zA);
    double r = std::sqrt(a * a + b * b + c * c);
    this->a = this->a / r;
    this->b = this->b / r;
    this->c = this->c / r;
    this->h = this->h / r;

    this->n = Vec3(this->a, this->b, this->c);
  }

  friend std::ostream &operator<<(std::ostream &out, const Triangle &Triangle);
  void draw(const int &d, const int &xShift, const int &yShift) {
    DrawLine(d * A.viewCoords.x / A.viewCoords.z + xShift,
             -d * A.viewCoords.y / A.viewCoords.z + yShift,
             d * B.viewCoords.x / B.viewCoords.z + xShift,
             -d * B.viewCoords.y / B.viewCoords.z + yShift, GREEN);
    DrawLine(d * A.viewCoords.x / A.viewCoords.z + xShift,
             -d * A.viewCoords.y / A.viewCoords.z + yShift,
             d * C.viewCoords.x / C.viewCoords.z + xShift,
             -d * C.viewCoords.y / C.viewCoords.z + yShift, GREEN);
    DrawLine(d * C.viewCoords.x / C.viewCoords.z + xShift,
             -d * C.viewCoords.y / C.viewCoords.z + yShift,
             d * B.viewCoords.x / B.viewCoords.z + xShift,
             -d * B.viewCoords.y / B.viewCoords.z + yShift, GREEN);
  }
};

class WFModel {
public:
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  std::vector<Triangle> triangles;
  WFModel() {}
  WFModel(std::vector<Vertex> &vertices, std::vector<Edge> &edges,
          std::vector<Triangle> &triangles) {
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
  int d = 300;
  int xShift, yShift;
  Scene(Viewpoint &View, WFModel &Figure, int xShift, int yShift) {
    this->Camera = &View;
    this->Figure = &Figure;
    this->xShift = xShift;
    this->yShift = yShift;
    Figure.computeSphereCoords(*Camera);
  }
  void changeView(const Vec3 &direction) {
    this->Camera->move(direction);
    Figure->computeSphereCoords(*Camera);
  }

  void changeView(const double &rho, const double &theta, const double &phi) {
    this->Camera->move(rho, theta, phi);
    Figure->computeSphereCoords(*Camera);
  }

  bool testVisibility(const Edge &edge) {
    bool visible = true;
    for (Triangle &trig : this->Figure->triangles) {
      trig.computePlane();
      if (trig.h < -MAX_PREC) {
        // std::cout << "skipped trig with h<0\n";
        continue;
      }
      Vec3 P = edge.start.viewCoords;
      Vec3 Q = edge.end.viewCoords;
      Vec3 A = trig.A.viewCoords;
      Vec3 B = trig.B.viewCoords;
      Vec3 C = trig.C.viewCoords;

      // TEST 1
      double eps1 = MAX_PREC + MAX_PREC * trig.h;
      double hP = (trig.a * P.x + trig.b * P.y + trig.c * P.z);
      double hQ = (trig.a * Q.x + trig.b * Q.y + trig.c * Q.z);

      if (hP <= (trig.h + eps1) && hQ <= (trig.h + eps1)) {
        continue;
      }

      // TEST 2
      double K1 = P.y * Q.z - Q.y * P.z;
      double K2 = Q.x * P.z - P.x * Q.z;
      double K3 = P.x * Q.y - Q.x * P.y;
      double dA = K1 * A.x + K2 * A.y + K3 * A.z;
      double dB = K1 * B.x + K2 * B.y + K3 * B.z;
      double dC = K1 * C.x + K2 * C.y + K3 * C.z;
      int eA = dA > MAX_PREC ? 1 : dA < -MAX_PREC ? -1 : 0;
      int eB = dB > MAX_PREC ? 1 : dB < -MAX_PREC ? -1 : 0;
      int eC = dC > MAX_PREC ? 1 : dC < -MAX_PREC ? -1 : 0;

      if (abs(eA + eB + eC) >= 2) {
        continue;
      }
      // TEST 3

      bool Poutside = false, Qoutside = false, outside = false;

      double lamI = 1., lamJ = 0.;

      for (int i = 0; i < 3; i++) {
        double lambda, mu;
        double C1 = A.y * B.z - B.y * A.z;
        double C2 = A.z * B.x - B.z * A.x;
        double C3 = A.x * B.y - B.x * A.y;
        double Cpos = C1 * C.x + C2 * C.y + C3 * C.z;
        double Ppos = C1 * P.x + C2 * P.y + C3 * P.z;
        double Qpos = C1 * Q.x + C2 * Q.y + C3 * Q.z;
        bool Pbeyond, Qbeyond;
        if (Cpos > MAX_PREC) {
          Pbeyond = Ppos < -MAX_PREC;
          Qbeyond = Qpos < -MAX_PREC;
          // one point is opposite to C and another is opposite or on edge
          outside =
              (Pbeyond && Qpos < MAX_PREC) || (Qbeyond && Ppos < MAX_PREC);
        } else if (Cpos < -MAX_PREC) {
          Pbeyond = Ppos > MAX_PREC;
          Qbeyond = Qpos > MAX_PREC;
          outside =
              (Pbeyond && Qpos > -MAX_PREC) || (Qbeyond && Ppos > -MAX_PREC);
        } else {
          // if point C is on plane EAB, then PQ is visible
          outside = true;
        }
        if (outside) {
          break;
        }

        Poutside = Poutside || Pbeyond;
        Qoutside = Qoutside || Qbeyond;
        lambda =
            std::abs(Qpos - Ppos) <= MAX_PREC ? 1e7 : -Ppos / (Qpos - Ppos);
        mu = std::abs(dB - dA) <= MAX_PREC ? 1e7 : -dA / (dB - dA);

        if (mu >= -MAX_PREC && mu <= 1 + MAX_PREC && lambda >= -MAX_PREC &&
            lambda <= 1 + MAX_PREC) {
          if (lambda < lamI)
            lamI = lambda;
          if (lambda > lamJ)
            lamJ = lambda;
        }
        Vec3 temp1 = A;
        A = B;
        B = C;
        C = temp1;
        double temp2 = dA;
        dA = dB;
        dB = dC;
        dC = temp2;
      }

      if (outside) {
        continue;
      }

      // TEST 4
      if (!(Poutside || Qoutside)) { // P and Q both inside
        visible = false;
        break;
      }

      // TEST 5: PQ interescts pyramid EABC before plane ABC

      double r1 = Q.x - P.x;
      double r2 = Q.y - P.y;
      double r3 = Q.z - P.z;

      Vec3 vI = Vec3(P.x + lamI * r1, P.y + lamI * r2, P.z + lamI * r3);
      Vec3 vJ = Vec3(P.x + lamJ * r1, P.y + lamJ * r2, P.z + lamJ * r3);
      if (trig.a * vI.x + trig.b * vI.y + trig.c * vI.z < trig.h - eps1) {
        continue;
      }
      if (trig.a * vJ.x + trig.b * vJ.y + trig.c * vJ.z < trig.h - eps1) {
        continue;
      }

      // TEST 6

      if (Poutside) {
        Vertex I;
        I.viewCoords = vI;
        Edge IP = Edge(I, edge.start);
        if (this->testVisibility(IP)) {
          IP.draw(this->d, this->xShift, this->yShift);
        }
      }

      if (Qoutside) {
        Vertex J;
        J.viewCoords = vJ;
        Edge JP = Edge(J, edge.end);
        if (this->testVisibility(JP)) {
          JP.draw(this->d, this->xShift, this->yShift);
        }
      }
      // all checks failed so triangle blocks view of this edge => we dont draw
      // it.
      visible = false;
      break;
    }
    return visible;
  }
  void draw() {

    // int i = 1;
    // for (auto &vertex : Figure->vertices) {
    //   vertex.draw(i, this->d, this->xShift, this->yShift);
    //   i++;
    // }

    // for (auto &trig : Figure->triangles) {
    //   trig.draw(d, xShift, yShift);
    // }

    for (auto &edge : Figure->edges) {
      if (this->testVisibility(edge)) {
        edge.draw(this->d, this->xShift, this->yShift);
      }
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
  std::cout << "Start: " << Edge.start << " End: " << Edge.end;
  return out;
}

inline std::ostream &operator<<(std::ostream &out, const Triangle &Triangle) {
  std::cout << "A: " << Triangle.A << " B: " << Triangle.B
            << " C: " << Triangle.C << '\n';
  return out;
}

inline std::ostream &operator<<(std::ostream &out, const Viewpoint &Camera) {
  std::cout << "\tCamera position: (" << Camera.rho << ", " << Camera.theta
            << ", " << Camera.phi << ")\n";
  return out;
}

inline std::ostream &operator<<(std::ostream &out, const WFModel &Figure) {
  std::cout << "Number of vertices = " << Figure.vertices.size() << "\n";
  for (size_t i = 0; i < Figure.vertices.size(); i++) {
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
  std::vector<Edge> edges;

  for (int i = 0; i < nEdges; i++) {
    std::getline(f, buffer);
    std::stringstream line(buffer);
    std::vector<std::string> cell(2);
    for (int j = 0; j < 2; j++) {
      std::getline(line, cell[j], ',');
    }
    // -1 shift to align with vector index
    int startIdx = std::stoi(cell[0]) - 1;
    int endIdx = std::stoi(cell[1]) - 1;
    edges.push_back(Edge(Figure.vertices[startIdx], Figure.vertices[endIdx]));
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
    int idx1 = std::stoi(cell[0]) - 1;
    int idx2 = std::stoi(cell[1]) - 1;
    int idx3 = std::stoi(cell[2]) - 1;
    triangles.push_back(Triangle(Figure.vertices[idx1], Figure.vertices[idx2],
                                 Figure.vertices[idx3]));
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
