#include "graph.h"
#include "raylib.h"
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

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

  graph::Vec3 pos1(3., 3., 4.);
  graph::Viewpoint Camera1(pos1);

  graph::Vec3 pos2(10., 10., 10.);
  graph::Viewpoint Camera2(pos2);
  graph::Scene Scene(Camera1, Model);

  int width = 900, height = 700;
  InitWindow(width, height, "Screen projection");
  SetTargetFPS(60);

  while (!WindowShouldClose()) {
    BeginDrawing();
    {
      ClearBackground(BLACK);
      // Scene.change_view(graph::Vec3(-0.01, -0.05, 0.01));
      Scene.changeView(0, 0.01, 0);
      Scene.draw(width / 2, height / 2);
    }
    EndDrawing();
  }
  return 0;
}
