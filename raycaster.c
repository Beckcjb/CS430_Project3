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

// Light color clamp
double clamp(double colorVal){
    double max = 1;
    double min = 0;

	if(colorVal > max) return max;

    else if (colorVal < min) return min;

    else   return colorVal;
}


static inline double v3_len(double* x)
{
    return sqrt(sqr(x[0]) + sqr(x[1]) + sqr(x[2]));
}



void v3_reflect(double* x, double* y, double* z)
{
    double scalar = 2.0 * v3_dot(x, y);

    double tmp[3];

    //printf("not working");
    v3_scale(y, scalar, tmp);

    v3_subtract(x, tmp, z);

}
int shadows(Object objects[], double* newRd, double* newRo, int items, int closestObject, double maxDistance)
{
     // Do intersections with new Ron and Rdn of all other objects in the scene(planes or spheres)
    int k;
    int newBest_o = -1;
   // printf("%d", closestObject);
    normalize(newRd);
    double newBestT = INFINITY;

    for(k = 0; k < items; k++)
    {   // skip the closest object
        double newT = 0;
        if (k == closestObject)
        { //printf("hey");
            continue;
        }
        else if(strcmp(objects[k].type, "sphere") == 0){
            newT = sphereIntersection(objects[k].structures.sphere.position, objects[k].structures.sphere.radius, newRd, newRo);

        } else if(strcmp(objects[k].type, "plane") == 0){
            newT = planeIntersection(objects[k].structures.plane.position, objects[k].structures.plane.normal, newRd, newRo);
        }
        // Shade the pixel because there is another object in the way of the light source

        if (maxDistance != INFINITY && newT > maxDistance)

            continue;

        if (newT > 0 && newT < newBestT)
        {
            //best_t = t;

            newBest_o = k;

        }

    }
    // found no intersections
    return newBest_o;

}

// Find the diffuse for the closest object
// takes in the object normal, light vector, light color, and the objects diffuse color
// places the result into the outColor
void diffuseHandle(double *objNormal, double *light, double *illumColor, double *objDiffuse, double *outColor) {

    // K_a*I_a should be added to the beginning of this whole thing, which is a constant and ambient light
    double normDotLight = v3_dot(objNormal, light);

    if (normDotLight > 0)
    {
        double diffuseProduct[3];

        diffuseProduct[0] = objDiffuse[0] * illumColor[0];
        diffuseProduct[1] = objDiffuse[1] * illumColor[1];
        diffuseProduct[2] = objDiffuse[2] * illumColor[2];

        // multiply by n_dot_l and store in out_color
        v3_scale(diffuseProduct, normDotLight, outColor);
        /*printf("Color is: %lf", outColor[1]);
        printf("Color is: %lf", outColor[2]);*/
    }

    else
    {
        // would normally return K_a*I_a here...
        outColor[0] = 0;
        outColor[1] = 0;
        outColor[2] = 0;

    }
}

// Find the specular for the closest object
// Takes in the ns which as shown in class we can set to one for testing at least
// light vector, light reflection, object normal, ray direction(V), object specular, light color
// Returns the proper color in the out color
void specularHandle(double ns, double *light, double *lightRef, double *objNormal,
                        double *V, double *objSpecular, double *illumColor, double *outColor) {

    double rayDotLight = v3_dot(light, lightRef);
    double normDotLight = v3_dot(objNormal, light);

    if (rayDotLight > 0 && normDotLight > 0)
    {
        double rayPowerOfNS = pow(rayDotLight, ns);

        double specularProduct[3];
        specularProduct[0] = objSpecular[0] * illumColor[0];
        specularProduct[1] = objSpecular[1] * illumColor[1];
        specularProduct[2] = objSpecular[2] * illumColor[2];

        v3_scale(specularProduct, rayPowerOfNS, outColor);

    }

    else
    {
        outColor[0] = 0;
        outColor[1] = 0;
        outColor[2] = 0;
    }

}

double angular_attenuation(Object objects[], double objectray[3], int items, int light_on)
{
        //check if light is not a spotlight. if not, return 1.0
            double lightDirectionDP = v3_dot(objects[light_on].structures.light.direction, objectray);
            double fang = pow(lightDirectionDP, objects[light_on].structures.light.angular_a0);
            return 1.0;

}

double radial_attenuation(double a1, double a2, double a0, double distance)
{
    if(distance == INFINITY)
    {
        return 1;
    }

    return 1/(a2*pow(distance,2) + a1*distance + a0);


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
					t = sphere_intersect(objects[i].structures.sphere.position, objects[i].structures.sphere.radius,Rd, Ro);
				} 
				else if(strcmp(objects[i].type, "plane") == 0){
					t = plane_intersect(objects[i].structures.plane.position, objects[i].structures.plane.normal, Rd, Ro);
				}

				if (t > 0 && t < best_t){
								// location on view plane
					best_t = t;
                    best_i = i;
				}
				int l,k;
                double* color = malloc(sizeof(double)*3);

									// color of the current intersection into the image buffer
				if(best_t > 0 && best_t != INFINITY)
                {
                    if(strcmp(objects[best_i].type, "sphere") == 0)
                    {  
                        color[0] = 0;//ambientColor[0];
                        color[1] = 0;//ambientColor[1];
                        color[2] = 0;//ambientColor[2];
                        // In order to find a shadow
                        for (l = 0; l < items; l++)
                        {
                            // Look for a light to see if that object has a shadow casted on it by a light
                            if(strcmp(objects[l].type, "light") == 0)
                            {   
                                double temp[3];
                                double Ron[3];
                                double Rdn[3];
                                v3_scale(Rd, best_t, temp);
                                v3_add(temp, Ro, Ron);
                                v3_subtract(objects[l].structures.light.position, Ron, Rdn);

                                double light_dist = v3_len(Rdn);
                                normalize(Rdn);

                                double shadow_intersect = shadows(objects, Rdn, Ron, items, best_i, light_dist);
                                // there is an object in the way so ad shadow 
								if(shadow_intersect != -1)
                                {   
                                    continue;
                                }
                                // no object in between item and the light source
                                else
                                {
                                 
                                    double sphere_position[3] = {objects[best_i].structures.sphere.position[0],objects[best_i].structures.sphere.position[1],objects[best_i].structures.sphere.position[2]};

                                    double n[3] = {Ron[0] - sphere_position[0], Ron[1]-sphere_position[1], Ron[2]-sphere_position[2]}; 

                                    normalize(n);
                                    double vector_L[3] = {Rdn[0], Rdn[1], Rdn[2]}; 
                                    double distanceVector[3] = {Ron[0] - objects[best_i].structures.light.position[0], Ron[1] - objects[best_i].structures.light.position[1],Ron[2] - objects[best_i].structures.light.position[2]};

                                    double reflection_L[3];
                                    double V[3] = {Rd[0], Rd[1], Rd[2]}; 
									double object_light_range[3];
                                    double diffuseColor[3];
                                    double specularColor[3];
                                    double diffuseSpecular[3];
                                    v3_reflect(vector_L, n, reflection_L); 

                                    diffuseHandle(n, vector_L, objects[l].structures.light.color, objects[best_i].structures.sphere.diffuseColor, diffuseColor);
                                    specularHandle(1, vector_L, reflection_L, n, V, objects[best_i].structures.sphere.specularColor, objects[l].structures.light.color, specularColor);

                                    v3_add(diffuseColor, specularColor, diffuseSpecular);

                                    v3_scale(vector_L, -1, object_light_range);

                                    double fang = angular_attenuation(objects, object_light_range, items, l);
                                    double frad = radial_attenuation(objects[l].structures.light.radial_a1, objects[l].structures.light.radial_a2, objects[l].structures.light.radial_a0, light_dist);
                                    
                                    color[0] += frad * fang * diffuseSpecular[0];
                                    color[1] += frad * fang * diffuseSpecular[1];
                                    color[2] += frad * fang * diffuseSpecular[2];

                                    buffer->image[y*3 * buffer->width + x*3].r = clamp(color[0]) *255;
                                    buffer->image[y*3 * buffer->width + x*3+1].g = clamp(color[1]) *255;
                                    buffer->image[y*3 * buffer->width + x*3+2].b = clamp(color[2]) *255;
                                }

                            }

                        }
					else if(strcmp(objects[best_i].type, "plane") == 0)
                    {
                        color[0] = 0;						//ambientColor[0];
                        color[1] = 0;						//ambientColor[1];
                        color[2] = 0;						//ambientColor[2];
                        // In order to find a shadow
                        for (l = 0; l < items; l++)
                        {
                            // Look for a light to see if that object has a shadow casted on it by a light
                            if(strcmp(objects[l].type, "light") == 0)
                            {   // calc new ray origin and direction
                                double temp[3];
                                double Ron[3];
                                double Rdn[3];
                                v3_scale(Rd, best_t, temp);
                                v3_add(temp, Ro, Ron);
                                v3_subtract(objects[l].structures.light.position, Ron, Rdn);
                                normalize(Rdn);
                                double light_dist = v3_len(Rdn);

                                double shadow_intersect = shadows(objects, Rdn, Ron, items, best_i, light_dist);
                                // there is an object in the way so shade in
                                if(shadow_intersect != -1)
                                {   

                                    continue;
                                }
													//  no object in between item and the light source
                                else
                                {
                                    double plane_position[3] = {objects[best_i].structures.plane.position[0],objects[best_i].structures.plane.position[1],objects[best_i].structures.plane.position[2]};

                                    double vector_L[3] = {Rdn[0], Rdn[1], Rdn[2]}; 
                                    double reflection_L[3];
                                    double V[3] = {Rd[0], Rd[1], Rd[2]};
                                    double diffuseColor[3];
                                    double specularColor[3];
                                    double diffuseSpecular[3];
                                    double object_light_range[3];
                                    v3_reflect(vector_L, Ron, reflection_L); 

                                    diffuseHandle(Ron, vector_L, objects[l].structures.light.color, objects[best_i].structures.plane.diffuseColor, diffuseColor);
                                    specularHandle(1, vector_L, reflection_L, Ron, V, objects[best_i].structures.plane.specularColor, objects[l].structures.light.color, specularColor);

                                    v3_add(diffuseColor, specularColor, diffuseSpecular);

                                    v3_scale(vector_L, -1, object_light_range);

                                    double fang = angular_attenuation(objects, object_light_range, items, l);
                                    double frad = radial_attenuation(objects[l].structures.light.radial_a1, objects[l].structures.light.radial_a2, objects[l].structures.light.radial_a0, light_dist);
									color[0] += frad * fang * diffuseSpecular[0];
                                    color[1] += frad * fang * diffuseSpecular[1];
                                    color[2] += frad * fang * diffuseSpecular[2];

                                    buffer->image[y*3 * buffer->width + x*3].r = clamp(color[0]) *255;
                                    buffer->image[y*3 * buffer->width + x*3+1].g = clamp(color[1]) *255;
                                    buffer->image[y*3 * buffer->width + x*3+2].b = clamp(color[2]) *255;
                                }

                            }
                        }
                    }

				}
				else{ 
													// no intersection results in background color
                    buffer->image[y*3 * buffer->width + x*3].r = 0*255 ;
                    buffer->image[y*3 * buffer->width + x*3+1].g = 0*255;
                    buffer->image[y*3 * buffer->width + x*3+2].b = 0*255;
				}
					free(color);			// free color buffer
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