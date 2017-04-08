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

int main(int argc, char **argv){
    int rv = 0;
    if (argc < 3) {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        return -1;
    }
    Model *model = loadFromObj(argv[2]);
    int height = 1600;
    int width = 1600;
    tgaImage * image = tgaNewImage(height, width, RGB);
//  meshgrid_1(image, model);
    meshgrid_2(image, model);
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
   if (a > 0) {
      return 1;
   }
   if (a == 0) {
      return 0;
   }
   if (a < 0) {
      return -1;
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
  int i, j;
  for (i = 0; i < 3; i++){
      for (j = 2; j > i; j--){
         if (coords[j - 1][1] > coords[j][1]){
            swap(&coords[j - 1][0], &coords[j][0]);
            swap(&coords[j - 1][1], &coords[j][1]);  
         }
      }
  }
  int xa, xb, x, y;
  if ((coords[1][1] == coords[2][1]) && (coords[0][1] != coords[2][1])){     
      for (y = coords[0][1]; y <= coords[2][1]; ++y){
          xa = (y - coords[0][1]) * (coords[0][0] - coords[1][0]) / (coords[0][1] - coords[1][1]) + coords[0][0];
          xb = (y - coords[0][1]) * (coords[0][0] - coords[2][0]) / (coords[0][1] - coords[2][1]) + coords[0][0];
         if (xa > xb){
              swap(&xa, &xb);
         }    
         for (x = xa; x <= xb; ++x){
              tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);
         }           
      }
  } 
  if ((coords[0][1] == coords[1][1]) && (coords[0][1] != coords[2][1])){      
      for (y = coords[0][1]; y <= coords[2][1]; ++y){
          xa = (y - coords[1][1]) * (coords[1][0] - coords[2][0]) / (coords[1][1] - coords[2][1]) + coords[1][0];
          xb = (y - coords[0][1]) * (coords[0][0] - coords[2][0]) / (coords[0][1] - coords[2][1]) + coords[0][0];
         if (xa > xb){
              swap(&xa, &xb);
         }    
         for (x = xa; x <= xb; ++x){
              tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);
         }           
      }
  }
  if (coords[0][1] == coords[2][1]){
      for (i = 0; i < 3; i++){
          for (j = 2; j > i; j--){
              if (coords[j - 1][0] > coords[j][0]){
                  swap(&coords[j - 1][0], &coords[j][0]);
                  swap(&coords[j - 1][1], &coords[j][1]);  
              }
          }
      }
      y = coords[0][1];
      for (x = coords[0][0]; x <= coords[2][0]; ++x){
          tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);
      }    
  }
  if ((coords[0][1] != coords[1][1]) && (coords[1][1] != coords[2][1]) && (coords[0][1] != coords[2][1])){
      for (y = coords[0][1]; y <= coords[2][1]; ++y){
          if (y <= coords[1][1]){
              xa = (y - coords[0][1]) * (coords[0][0] - coords[1][0]) / (coords[0][1] - coords[1][1]) + coords[0][0];
              xb = (y - coords[0][1]) * (coords[0][0] - coords[2][0]) / (coords[0][1] - coords[2][1]) + coords[0][0];
          } else{
              xa = (y - coords[1][1]) * (coords[1][0] - coords[2][0]) / (coords[1][1] - coords[2][1]) + coords[1][0];
              xb = (y - coords[0][1]) * (coords[0][0] - coords[2][0]) / (coords[0][1] - coords[2][1]) + coords[0][0];
          }
          if (xa > xb){
              swap(&xa, &xb);
          }    
          for (x = xa; x <= xb; ++x){
              tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);
          } 
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
