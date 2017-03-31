#include <stdio.h>
#include <math.h>
#include "model.h"
#include "tga.h"

void swap(int *a, int *b);
int iabs(int a);
void line (tgaImage *image, int x0, int y0, int x1, int y1, tgaColor color);
int sign(int a);
void meshgrid(tgaImage *image, Model *model);

int main(int argc, char **argv){
    int rv = 0;
    if (argc < 3) {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        return -1;
    }
    Model *model = loadFromObj(argv[2]);

    int height = 800;
    int width = 800;
    tgaImage * image = tgaNewImage(height, width, RGB);
    
    meshgrid(image, model);

    if (-1 == tgaSaveToFile(image, argv[1])) {
        perror("tgaSateToFile");
        rv = -1;
    }

    tgaFreeImage(image);
    freeModel(model);    
    return rv;
}

void line (tgaImage *image, int x0, int y0, int x1, int y1, tgaColor color){
    int steep = 0;
    if (iabs(y1 - y0) > iabs(x1 - x0)) {
        steep = 1;
        swap(&x0, &y0);
        swap(&x1, &y1);
    }
    if (x0 > x1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
    }
    int dx = x1 - x0;
    int dy = y1 - y0;
    int de = 2 * iabs(dy);
    int e = 0;
    int y = y0;
    int x;
    for (x = x0; x <= x1; ++x){
         if (steep == 1) {
            tgaSetPixel(image, (unsigned int)y, (unsigned int)x, color);
         } else {
            tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);  
         }
         e = e + de;
         if (e > dx){
            y = y + sign(dy);
            e = e - 2 * dx;
         }
    } 
}

int sign(int a){
   int b;
   if (a > 0) {
      b = 1;
   }
   if (a == 0) {
      b = 0;
   }
   if (a < 0) {
      b = -1;
   }
   return b;
}

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

int iabs(int a) {
    return (a >= 0) ? a : -a;
}
void meshgrid(tgaImage *image, Model *model) {
	int i, j;
	tgaColor white = tgaRGB(255, 255, 255);
	for (i = 0; i < model->nface; ++i) {
		int screen_coords[3][2];
		for (j = 0; j < 3; ++j) {
			Vec3 *v = &(model->vertices[model->faces[i][3 * j]]);
			screen_coords[j][0] = ((*v)[0] + 1) * image->width / 2;
			screen_coords[j][1] = (1 - (*v)[1]) * image->height / 2;
		}
		for (j = 0; j < 3; ++j) {
			line(image, screen_coords[j][0], screen_coords[j][1], screen_coords[(j + 1) % 3][0], screen_coords[(j + 1) % 3][1], white);
		}
	}
}
