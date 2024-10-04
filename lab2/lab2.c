#include <stdio.h>

#include <raylib.h>
#include<stdlib.h>

typedef struct Vector2_int {
    int x;
    int y;
} Vector2_int;

void my_drawline(Vector2_int p1, Vector2_int p2, Color color);
void my_drawcircle(Vector2_int p1, int radius,  Color color);

int main() {
    int width = 900, height = 600;
    Vector2_int p1 = {.x=300, .y=350};
    Vector2_int p2 = {.x=290, .y=400};
    int radius = 250;
    InitWindow(width, height, "Bresenham's line algorithm");
    SetTargetFPS(60);

    int i = 0;
    while (!WindowShouldClose()) {
        BeginDrawing();
        {
            my_drawcircle(p1,  radius, BLUE);
            // DrawCircleLines(p1.x, p1.y, (float) radius, RED);
            my_drawline(p1, p2, BLUE);
            DrawLine(p1.x, p1.y, p2.x, p2.y, RED);
        }
        EndDrawing();
    }
    return 0;
}

void my_drawline(Vector2_int p1, Vector2_int p2, Color color) {
    Vector2_int p0 = {.x=p1.x, .y=p1.y}; // curr point
    int dx = abs(p2.x-p1.x);
    int dy = abs(p2.y-p1.y);
    int e_xy = dx - dy;

    int sx = p1.x < p2.x ? 1: -1;
    int sy = p1.y < p2.y ? 1: -1;

    int e_x = -dy;
    int e_y = dx;

    while (p0.x != p2.x || p0.y != p2.y) {
        DrawPixel(p0.x, p0.y, color);
        int err2 = e_xy*2;

        if (err2 >= e_x) {
            p0.x += sx;
            e_xy += e_x;
        }

        if (err2 <= e_y) {
            p0.y += sy;
            e_xy += e_y;
        }
    }
    return;
}

void my_drawcircle(Vector2_int p1, int radius,  Color color) {
    Vector2_int p0 = {.x=0, .y=radius}; // curr point
    int err = 2 - 2 * radius;

    DrawPixel(p1.x, p1.y, color);
    while (p0.x <= p0.y) { // 1/8 of circle
        DrawPixel(p1.x+p0.x, p1.y+p0.y, color);
        DrawPixel(p1.x+p0.x, p1.y-p0.y, color);
        DrawPixel(p1.x-p0.x, p1.y+p0.y, color);
        DrawPixel(p1.x-p0.x, p1.y-p0.y, color);
        DrawPixel(p1.x+p0.y, p1.y+p0.x, color);
        DrawPixel(p1.x+p0.y, p1.y-p0.x, color);
        DrawPixel(p1.x-p0.y, p1.y+p0.x, color);
        DrawPixel(p1.x-p0.y, p1.y-p0.x, color);
        if (err > 0) {
            err += 4*(p0.x-p0.y) + 10;
            p0.y--;
        }
        else {
            err += 4*p0.x + 6;
        }
        p0.x++;
    }
    return;
}

