#include "graph.h"
#include "raylib.h"
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  std::cout << std::fixed;
  std::cout << std::setprecision(16);

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

  graph::Vec3 pos1(-4, .33, 2.);
  graph::Viewpoint Camera1(pos1);

  int width = 900, height = 700;
  graph::Vec3 pos2(8., 11., 8.);
  graph::Viewpoint Camera2(pos2);

  graph::Vec3 pos3(280., 210., 280.);
  graph::Viewpoint Camera3(pos3);

  graph::Scene Scene(Camera2, Model, width / 2, height / 2);

  InitWindow(width, height, "Screen projection");
  SetTargetFPS(60);

  while (!WindowShouldClose()) {
    BeginDrawing();
    {
      ClearBackground(BLACK);
      if (!IsKeyUp(KEY_SPACE)) {
        Scene.changeView(0, 0.01, 0.);
      }
      if (!IsKeyUp(KEY_W)) {
        Scene.changeView(-1, 0., 0.);
      }
      if (!IsKeyUp(KEY_S)) {
        Scene.changeView(1, 0., 0.);
      }
      if (!IsKeyUp(KEY_Q)) {
        Scene.changeView(graph::Vec3(0., 0.5, 1));
      }
      if (!IsKeyUp(KEY_A)) {
        Scene.changeView(graph::Vec3(0., -0.5, -1));
      }
      Scene.draw();
    }
    EndDrawing();
  }
  return 0;
}
