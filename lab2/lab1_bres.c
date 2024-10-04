#include <raylib.h>
#include <stdlib.h>
#include <raymath.h>
#include <math.h>

#if 1
Color setColor() {
    return BLUE;
}
#else
Color setColor() {
    int r = rand() % 255 + 1;
    int g = rand() % 255 + 1;
    int b = rand() % 255 + 1;

    return (Color) {
        r,g,b,255
    };
}
#endif

// Color randColor() {
//     int r = rand() % 255 + 1;
//     int g = rand() % 255 + 1;
//     int b = rand() % 255 + 1;
//
//     return (Color) {
//         r,g,b,255
//     };
// }
typedef struct Vector2_int {
    int x;
    int y;
} Vector2_int;

void my_drawline(Vector2_int p1, Vector2_int p2, Color color);

void pifagor(int N, Vector2 startPos, int baseLen, float baseAngle, float roofAngle);
Vector2 my_DrawSquare(Vector2 v1, int baseLen, float baseAngle, Color color);
Vector2 my_DrawTriangle(Vector2 v1, int baseLen, float baseAngle, float angle, Color color);

int main() {
    int width = 900, height = 600;
    InitWindow(width, height, "Pythagoras tree");
    SetTargetFPS(60);

    /* Initial values */
    int N = 15;
    Vector2 startPos = {400,500};
    int baseLen = 100; // base length
    float baseAngle = 0; // base angle
    float roofAngle = PI/6; // roof angle

    while (!WindowShouldClose()) {
        // baseAngle+=0.01;
        // roofAngle+=0.01;
        BeginDrawing();
        {
            // ClearBackground(BLACK);
            pifagor(N, startPos, baseLen, baseAngle, roofAngle);
        }
        EndDrawing();
    }
    CloseWindow();
    return 0;
}

void my_drawline(Vector2_int p1, Vector2_int p2, Color color) {
    Vector2_int p0 = {.x=p1.x, .y=p1.y}; // curr point
    int dx = abs(p2.x-p1.x);
    int dy = abs(p2.y-p1.y);
    // int D = 2*dy-dx;
    int e_xy = dx - dy;

    int sx = p1.x < p2.x ? 1: -1;
    int sy = p1.y < p2.y ? 1: -1;

    int e_x = -dy;
    int e_y = dx;

    // printf("Start:\n");
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
        // printf("x=%d, y=%d, err=%d\n", p0.x, p0.y, e_xy);
    }
    // printf("Stop:\n");
    return;
}
Vector2 my_DrawSquare(Vector2 v1, int baseLen, float baseAngle, Color color) {
    /* draws a rectangle starting from bottom-left vertex.
     * Vertices are calculated anti-clockwise.
     * returns last vertice v4 */
    Vector2 v2 = Vector2Add(v1, (Vector2) {
        baseLen*cosf(baseAngle), -baseLen*sinf(baseAngle)
    });

    Vector2 v3 = Vector2Add(v2, (Vector2) {
        baseLen*cosf(baseAngle+PI/2.), -baseLen*sinf(baseAngle+PI/2.)
    });

    Vector2 v4 = Vector2Add(v3, (Vector2) {
        baseLen*cosf(baseAngle+PI), -baseLen*sinf(baseAngle+PI)
    });

    // DrawLineV(v1, v2, color);
    // DrawLineV(v2, v3, color);
    // DrawLineV(v3, v4, color);
    // DrawLineV(v4, v1, color);
    my_drawline( (Vector2_int) {
        (int) v1.x,(int) v1.y
    }, (Vector2_int) {
        (int) v2.x,(int) v2.y
    }, color);
    my_drawline( (Vector2_int) {
        (int) v2.x,(int) v2.y
    }, (Vector2_int) {
        (int) v3.x,(int) v3.y
    }, color);
    my_drawline( (Vector2_int) {
        (int) v3.x,(int) v3.y
    }, (Vector2_int) {
        (int) v4.x,(int) v4.y
    }, color);
    my_drawline( (Vector2_int) {
        (int) v4.x,(int) v4.y
    }, (Vector2_int) {
        (int) v1.x,(int) v1.y
    }, color);
    return v4;
}

Vector2 my_DrawTriangle(Vector2 v1, int baseLen, float baseAngle, float angle, Color color) {
    /* draws a triangle starting from bottom-left vertex.
     * Vertices are calculated anti-clockwise.
     * returns last vertice v3*/
    Vector2 v2 = Vector2Add(v1, (Vector2) {
        baseLen*cosf(baseAngle), -baseLen*sinf(baseAngle)
    });
    Vector2 v3 = Vector2Add(v1, (Vector2) {
        // some fancy geometry shenanigans
        baseLen*cosf(angle)*cosf(angle+baseAngle), -baseLen*cosf(angle)*sinf(angle+baseAngle)
    });
    // DrawLineV(v2, v3, color);
    // DrawLineV(v3, v1, color);

    my_drawline( (Vector2_int) {
        (int) v2.x,(int) v2.y
    }, (Vector2_int) {
        (int) v3.x,(int) v3.y
    }, color);
    my_drawline( (Vector2_int) {
        (int) v3.x,(int) v3.y
    }, (Vector2_int) {
        (int) v1.x,(int) v1.y
    }, color);
    return v3;
}

void pifagor(int N, Vector2 startPos, int baseLen, float baseAngle, float roofAngle) {
    Vector2 triangleV1 = my_DrawSquare(startPos, baseLen, baseAngle, setColor());
    Vector2 rightSquareV1 = my_DrawTriangle(triangleV1, baseLen, baseAngle, roofAngle, setColor());

    if (--N > 0) {
        // left
        pifagor(N, triangleV1, baseLen*cosf(roofAngle), baseAngle+roofAngle, roofAngle);
        // right
        pifagor(N, rightSquareV1, baseLen*sinf(roofAngle), -PI/2+baseAngle+roofAngle, roofAngle);
    }

}
