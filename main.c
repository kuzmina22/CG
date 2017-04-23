#include <stdio.h>
#include <math.h>
#include "model.h"
#include "tga.h"

void swap(int *a, int *b);
int iabs(int a);
int sign(int a);
void line (tgaImage *image, int x0, int y0, int x1, int y1, tgaColor color);
void meshgrid_1(tgaImage *image, Model *model);
void triangle(int coords[3][3], tgaColor color, tgaImage *image);
void normal(double coords[3][3], Vec3* n);
double intensity(Vec3 light, Vec3 n);
void meshgrid_2(tgaImage *image, Model *model);
void z_triangle(int coords[3][3], tgaColor color, tgaImage *image, int (*zbuffer)[700][700]);
void meshgrid_3(tgaImage *image, Model *model);

int main(int argc, char **argv){
    int rv = 0;
    if (argc < 3) {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        return -1;
    }
    Model *model = loadFromObj(argv[2]);
    int height = 700;
    int width = 700;
    tgaImage * image = tgaNewImage(height, width, RGB);
//  meshgrid_1(image, model);
//  meshgrid_2(image, model);
    meshgrid_3(image, model);
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
   if (a > 0) 
      b = 1;
   if (a == 0) 
      b = 0;
   if (a < 0) 
      b = -1;
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

void meshgrid_1(tgaImage *image, Model *model) {
	int i, j;
	tgaColor white = tgaRGB(255, 255, 255);
	for (i = 0; i < model->nface; ++i) {
		int screen_coords[3][3];
		for (j = 0; j < 3; ++j) {
			Vec3 *v = &(model->vertices[model->faces[i][3 * j]]);
                        screen_coords[j][0] = ((*v)[0] + 1) * image->width / 2;
			screen_coords[j][1] = (1 - (*v)[1]) * image->height / 2;
                        screen_coords[j][2] = (*v)[2];      
		}
	        for (j = 0; j < 3; ++j) {
			line(image, screen_coords[j][0], screen_coords[j][1], screen_coords[(j + 1) % 3][0], screen_coords[(j + 1) % 3][1], white);
		} 
        }
}

void triangle(int coords[3][3], tgaColor color, tgaImage *image){
    if (coords[0][1] == coords[1][0] && coords[0][1] == coords[2][0]) return;
    for (int i = 0; i < 3; i++){
        for (int j = 2; j > i; j--){
            if (coords[j - 1][1] > coords[j][1]){
                swap(&coords[j - 1][0], &coords[j][0]);
                swap(&coords[j - 1][1], &coords[j][1]);  
           }
        }
    } 
    int total_height = coords[2][1] - coords[0][1];
    for (int y = 0; y < total_height; y++) {
        _Bool second_half = y > coords[1][1] - coords[0][1] || coords[1][1] == coords[0][1];
        int segment_height = second_half ? coords[2][1] - coords[1][1] : coords[1][1] - coords[0][1];
        float alpha = (float)y/total_height;
        float beta  = (float)(y-(second_half ? coords[1][1] - coords[0][1] : 0))/segment_height;
        int xa = coords[0][0] + (coords[2][0] - coords[0][0])*alpha;
        int xb = second_half ? coords[1][0] + (coords[2][0] - coords[1][0])*beta : coords[0][0] + (coords[1][0] - coords[0][0])*beta;
        if (xa > xb){
            swap(&xa, &xb);
        }
        for (int x = xa; x <= xb; x++) {
            tgaSetPixel(image, (unsigned int)x, (unsigned int)coords[0][1] + y, color);
        }
    }
}

void normal(double coords[3][3], Vec3* n){
    double a[3], b[3];
    for (int j = 0; j < 3; ++j){
        a[j] = coords[1][j] - coords[0][j];
        b[j] = coords[2][j] - coords[0][j];
    }
    (*n)[0] = a[1] * b[2] - b[1] * a[2];
    (*n)[1] = - a[0] * b[2] + b[0] * a[2];
    (*n)[2] = a[0] * b[1] - b[0] * a[1];
    double absn = sqrt((*n)[0] * (*n)[0] + (*n)[1] * (*n)[1] + (*n)[2] * (*n)[2]);
    (*n)[0] = (*n)[0] / absn;
    (*n)[1] = (*n)[1] / absn;
    (*n)[2] = (*n)[2] / absn;
}

double intensity(Vec3 light, Vec3 n){
    double I;
    for (int i = 0; i < 3; ++i){
        I = light[i] * n[i];
    }
    return I;
}

void meshgrid_2(tgaImage *image, Model *model) {
	int i, j;
        Vec3 light = {0, 0, -1};
	for (i = 0; i < model->nface; ++i) {
                double coords[3][3];
		for (j = 0; j < 3; ++j) {
			Vec3 *v = &(model->vertices[model->faces[i][3 * j]]);            
                        coords[j][0] = (*v)[0];
	  		coords[j][1] = (*v)[1];
                        coords[j][2] = (*v)[2];       
		}        
                Vec3 n;
                normal(coords, &n);
                double I = intensity(light, n);
                if (I < 0){
                        I = (-1) * I;
                        int screen_coords[3][3];
                        for (j = 0; j < 3; ++j){
		        	screen_coords[j][0] = (coords[j][0] + 1) * image->width / 2;
	            		screen_coords[j][1] = (1 - coords[j][1]) * image->height / 2;
                    		screen_coords[j][2] = coords[j][2];    
            		}
            		tgaColor color = tgaRGB(I * 255, I * 255, I * 255);
            		triangle(screen_coords, color, image);         
        	}
    	}
}

void z_triangle(int coords[3][3], tgaColor color, tgaImage *image, int (*zbuffer)[700][700]){
    if (coords[0][1] == coords[1][0] && coords[0][1] == coords[2][0]) return;
    for (int i = 0; i < 3; i++){
        for (int j = 2; j > i; j--){
           if (coords[j - 1][1] > coords[j][1]){
              swap(&coords[j - 1][0], &coords[j][0]);
              swap(&coords[j - 1][1], &coords[j][1]);
              swap(&coords[j - 1][2], &coords[j][2]);  
           }
        }
    }
    int total_height = coords[2][1] - coords[0][1];
    for (int y = 0; y < total_height; y++) {
        _Bool second_half = y > coords[1][1] - coords[0][1] || coords[1][1] == coords[0][1];
        int segment_height = second_half ? coords[2][1] - coords[1][1] : coords[1][1] - coords[0][1];
        float alpha = (float)y/total_height;
        float beta  = (float)(y-(second_half ? coords[1][1] - coords[0][1] : 0))/segment_height;
        int xa = coords[0][0] + (float)(coords[2][0] - coords[0][0])*alpha;
        int za = coords[0][2] + (float)(coords[2][2] - coords[0][2])*alpha;
        int xb = second_half ? coords[1][0] + (float)(coords[2][0] - coords[1][0])*beta : coords[0][0] + (float)(coords[1][0] - coords[0][0])*beta;
        int zb = second_half ? coords[1][2] + (float)(coords[2][2] - coords[1][2])*beta : coords[0][2] + (float)(coords[1][2] - coords[0][2])*beta;
        if (xa > xb){
            swap(&xa, &xb);
            swap(&za, &zb);
        }
        for (int x = xa; x <= xb; x++) {
            float gamma = xa == xb ? 1.0 : (float)(x - xa)/(float)(xb - xa);
            int z = (float)za + (float)(zb - za)*gamma;
            if (z > (*zbuffer)[x][y]){
                (*zbuffer)[x][y] = z;
                tgaSetPixel(image, (unsigned int)x, (unsigned int)coords[0][1] + y, color);
            }
        }
    }    
}

void meshgrid_3(tgaImage *image, Model *model) {
	int i, j;
    	int zbuffer[700][700];
    	for (i = 0; i < 700; ++i){
        	for (j = 0; j < 700; ++j){
            		zbuffer[i][j] = -1;
        	}
    	}   
    	Vec3 light = {0, 0, -1};
	for (i = 0; i < model->nface; ++i) {
        	double coords[3][3];
		for (j = 0; j < 3; ++j) {
			Vec3 *v = &(model->vertices[model->faces[i][3 * j]]);            
            		coords[j][0] = (*v)[0];
			coords[j][1] = (*v)[1];
            		coords[j][2] = (*v)[2];       
		}        
       		Vec3 n;
        	normal(coords, &n);
        	double I = intensity(light, n);
        	if (I < 0){
            	I = (-1) * I;
            	int screen_coords[3][3];
 		for (j = 0; j < 3; ++j){
                	screen_coords[j][0] = (coords[j][0] + 1) * image->width / 2;
			screen_coords[j][1] = (1 - coords[j][1]) * image->height / 2;
                        screen_coords[j][2] = (coords[j][2] + 1) * 127.5;    
            	}
            	tgaColor color = tgaRGB(I * 255, I * 255, I * 255);
            	z_triangle(screen_coords, color, image, &zbuffer);         
        }
    }
}
