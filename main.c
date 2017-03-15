#include <stdio.h>
#include <math.h>
#include "model.h"
#include "tga.h"

void swap(int *a, int *b);
int iabs(int a);
void line (tgaImage *image, int x0, int y0, int x1, int y1, tgaColor color);
void meshgrid(tgaImage *image, Model *model);

int main(int argc, char **argv){
    int rv = 0;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        return -1;
    }

	Model *model = loadFromObj(argv[1]);

	int width = 800;
	int height = 800;
	tgaImage * image = tgaNewImage(height, width, RGB);

	meshgrid(image, model);

    /*int i;
    tgaColor white = tgaRGB(255, 255, 255);
    tgaColor red = tgaRGB(255, 0, 0);
    tgaColor blue = tgaRGB(0, 0, 255);
    for (i = 0; i < 1000000; ++i) {
        line(image, 13, 20, 90, 40, white);
        line(image, 20, 13, 40, 80, red);
        line(image, 80, 40, 13, 20, blue);
    }*/

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

    int x;
    double y;
    double k = ((double)(y1 - y0))/(x1 - x0);
    for (x = x0, y = y0; x <= x1; ++x, y += k) {
        if (steep) {
            tgaSetPixel(image, (unsigned int)y, (unsigned int)x, color);
        } else {
            tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);
        }
    }
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

