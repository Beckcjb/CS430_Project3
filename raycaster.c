// Charles Beck
// CS 430
// 10/4/16
// PROJECT 2: RAYCASTER
// ==================================================================================
/* 	   
	   This program reads in five arguments from the comand line in the order of :
		
		program_name width height output.ppm input.json
   
	   The program will then read in the arguments and store the data where necessary.
   Then the program will parse the deisred JSON file while error checking for a 
   well made and valid JSON file. After the program gets through parsing, the data
   is stored into structures that will allow the program to test intersections of pixels 
   with a ray that is constructed and shot at every pixel in the scene. The program will
   store the intersection data and then write it to a P3 ppm file format with color. 
*/  
// ==================================================================================

									// Includes contain necessary libraries,
									//  header files, and other c files for the program

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vector_math.h"
#include "parser.c"
#include "ppmwrite.c"


Object objects[128];				// Maximum number of objects

static inline double sqr(double v){	// Inline Square method , returns squared value

  return v*v;
}


static inline void normalize(double* v){
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

// Plane intersection method returns
//   where the plane intersects the ray 
//   and returns a t value.
double plane_intersect(double* p, double* n, double* Rd, double* Ro){

    normalize(n);
    double deDot = v3_dot(n, Rd);
    double c[3];
    v3_subtract(p, Ro, c);


    if (fabs(deDot) < 0.0001f){
		return -1;
	}   						

    double t = v3_dot(c, n) / deDot;

									// no intersection
    if(t < 0.0){
		return -1;
	} 

    return t;
}


// Sphere intersection method returns
//   where the Sphere intersects the ray 
//   and returns a t value.
double sphere_intersect(double* p, double r, double* Rd, double* Ro){

    double a, b, c;
									// calculate quadratic formula
 
    a = sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]);
    b = 2 * (Rd[0]*(Ro[0]-p[0]) + Rd[1]*(Ro[1]-p[1]) + Rd[2]*(Ro[2]-p[2]));
    c = sqr(Ro[0]-p[0]) + sqr(Ro[1]-p[1]) + sqr(Ro[2]-p[2]) - sqr(r);


									// normalized
    if (a > 1.0001 || a < .9999){
        fprintf(stderr, "Ray direction was not normalized correctly");
        exit(-1);

    }

									// check discriminant
    double disc = sqr(b) - 4*a*c;

    double t0, t1;  				// t value solutions

									// no intersection
    if (disc < 0){
			return -1;		
	}


    else if (disc == 0){			
        t0 = -1*(b / (2*a));
        return t0;
    }

  
    else{
        t0 = (-1*b - sqrt(sqr(b) - 4*c))/2;
        t1 = (-1*b + sqrt(sqr(b) - 4*c))/2;
    
    }

									// no intersection
    if (t0 < 0 && t1 < 0){
		return -1;
	}

    else if (t0 < 0 && t1 > 0){
		return t1;
	}
    else if (t0 > 0 && t1 < 0){
		return t0;
	}
   
    else
    {

        if (t0 <= t1){
			return t0;
		}
        else{ 
			return t1;}
    }
 

}

// This method casts a ray that tests the intersections of josn objects
// with the ray and stores the information.
int ray_cast(Object objects[], Pixmap * buffer, double width, double height, int items){
	double cx, cy, h, w, pixelHeight, pixelWidth;
	int i, x, y;
    double Ro[3] = {0, 0, 0};
	double Rd[3] = {0, 0, 0};
	double point[3] = {0,0, 1};
	double view[2] = {0,0};

	cx = 0;
	cy = 0;
	// buffer
	buffer->width = width;
	buffer->height = height;
	buffer->color = 255;

												//get the size of the view plane
	for(i = 0; i < items; i++){
		if(strcmp(objects[i].type, "camera") == 0){
			h = objects[i].structures.camera.height;
			w = objects[i].structures.camera.width;
		}
	}

	pixelHeight = h / height;
    pixelWidth = w / width;

   // invert y and go through pixels to obtain intersections
	for (y = 0; y < width; y++){
        point[1] = -(view[1] - h/2.0 + pixelHeight*(y + 0.5));
		for (x = 0; x < height; x++){
            point[0] = view[0] - w/2.0 + pixelWidth*(x + 0.5);
            normalize(point);
			// normalize Rd
			Rd[0] = point[0];
            Rd[1] = point[1];
            Rd[2] = point[2];
            double best_t = INFINITY;
										// Go through each of the objects at each pixel 
										//   to carry out intersection calculaions.
			int best_i = 0;
			for (i = 0; i < items; i++){
				double t = 0;
				if(strcmp(objects[i].type, "sphere") == 0){
					// get intersections
					t = sphere_intersect(objects[i].structures.sphere.position, objects[i].structures.sphere.radius,Rd, Ro);
				} 
				else if(strcmp(objects[i].type, "plane") == 0){
					// get intersctions
					t = plane_intersect(objects[i].structures.plane.position, objects[i].structures.plane.normal, Rd, Ro);
				}

				if (t > 0 && t < best_t){
					// appearance on view plane
					best_t = t;
                    best_i = i;
				}

													// save color to buffer
				if(best_t > 0 && best_t != INFINITY){
					// store sphere data
                    if(strcmp(objects[best_i].type, "sphere") == 0){
                        buffer->image[y*3 * buffer->width + x*3].r = objects[best_i].structures.sphere.color[0] *255;
                        buffer->image[y*3 * buffer->width + x*3+1].g = objects[best_i].structures.sphere.color[1] *255;
                        buffer->image[y*3 * buffer->width + x*3+2].b = objects[best_i].structures.sphere.color[2] *255;
                    }
					// store plane data
                    else if(strcmp(objects[best_i].type, "plane") == 0){
                        buffer->image[y*3 * buffer->width + x*3].r = objects[best_i].structures.plane.color[0]*255 ;
                        buffer->image[y*3 * buffer->width + x*3+1].g = objects[best_i].structures.plane.color[1]*255;
                        buffer->image[y*3 * buffer->width + x*3+2].b = objects[best_i].structures.plane.color[2]*255;
                    }
				}
				else{ 
													// no intersection results in background color
                    buffer->image[y*3 * buffer->width + x*3].r = 0*255 ;
                    buffer->image[y*3 * buffer->width + x*3+1].g = 0*255;
                    buffer->image[y*3 * buffer->width + x*3+2].b = 0*255;
				}
			}

		}
	}

	return 0;
}

int main(int argc, char *argv[]){
	if (argc != 5) {										// Error is cmd arguments are not formated correctly
	fprintf(stderr, "Error: Usage width height output.ppm input.json");
	exit(1);
	}
	FILE *json;
	int items, i;
	int ppmFormat = 3;										// slected PPM format (forced)
	double width, height;
	int w, h;
	char* fileNameIn;
	fileNameIn = argv[5];
	
	width = atof(argv[1]);									// store desired width and height of image
	height = atof(argv[2]);

	
    Pixmap picbuffer;
    picbuffer.image = (PixelColor*)malloc(sizeof(PixelColor)*width* (height*3));

	json = fopen(argv[4], "r");
	if(json == NULL){											// if file DNE then end error
		fprintf(stderr, "Error: could not open file.\n");
		fclose(json);
		exit(1);
	}
	
	

	else{
		items = read_scene(json, objects);						// parse in objects
        ray_cast(objects, &picbuffer, width, height, items);	// call ray_cast and carry out intersections

        int size = height * width;

        ppmWriter(&picbuffer, argv[3], size , ppmFormat);		// send final fata to ppmWriter
	}

    fclose(json);						// close file
    free(picbuffer.image);				// free up malloced data

	return 0;
}