// Charles Beck
// Parser: Json files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "objectsStructs.h"

// Line number to report errors
int lineNumber = 1;


// Checks to see if the current character is new line
int next_c(FILE* json) {
    int c = fgetc(json);
    
    if (c == '\n')
    {
        lineNumber += 1;
    }

							// If end of file exit program and close the file
    if (c == EOF)
    {
        fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", lineNumber);
        fclose(json);
        exit(-1);
    }

    return c;
}


// Checks that the next character is what is expected
// If it is not then we will exit and return an error
void expect_c(FILE* json, int d) {
    int c = next_c(json);
    if (c == d) return;

    fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, lineNumber);
    fclose(json);
    exit(-1);
}


// Checks the current character in the file and if it is a space then we will skip it
void skip_ws(FILE* json) {
    int c = next_c(json);

    while (isspace(c))
    {
      c = next_c(json);
    }
    ungetc(c, json);
}


// Validation which checks the string for escape codes, if it is larger than ascii
// and if it is larger than 128
char* next_str(FILE* json)
{
    char buffer[129];
    int c = next_c(json), i = 0;

    if (c != '"')
    {
        fprintf(stderr, "Error: Expected string on line %d.\n", lineNumber);
        fclose(json);
        exit(-1);
    }
    c = next_c(json);

    while (c != '"')
    {
        if (i >= 128)
        {
          fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
          fclose(json);
          exit(-1);
        }
        if (c == '\\')
        {
          fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
          fclose(json);
          exit(-1);
        }
        if (c < 32 || c > 126)
        {
          fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
          fclose(json);
          exit(-1);
        }
        // Add the character to the buffer
        buffer[i] = c;
        i += 1;
        c = next_c(json);
    }
    // Null terminate the buffer and return the string for use
    buffer[i] = 0;
    return strdup(buffer);
}


double next_num(FILE* json)
{
    double value;
    if(fscanf(json, "%lf", &value) == 0)
    {
        fprintf(stderr, "Error, line number %d; expected numeric value.\n", lineNumber);
        fclose(json);
        exit(-1);
    }
    return value;
}

// Takes in an array containing a vector and returns a proper array of doubles
double* next_vec(FILE* json)
{
    double* vector = malloc(3*sizeof(double));

    expect_c(json, '[');
    skip_ws(json);
    vector[0] = next_num(json);

    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    vector[1] = next_num(json);

    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    vector[2] = next_num(json);

    skip_ws(json);
    expect_c(json, ']');

    return vector;
}

// Reads in a scene of objects taken from a JSON file
int read_scene(FILE *json, Object objects[]){
    int c, items;
	double *vector;
	char *key, *value;
	items = 0;
	skip_ws(json);
	c = next_c(json);
	
	if(c != '['){
		fprintf(stderr, "Error, line number %d; invalid scene definition '%c'\n", lineNumber, c);
		fclose(json);
		exit(-1);
	}

	skip_ws(json);
	c = next_c(json);

	// Check to see if the json file has no objects in it
	if(c != ']'){
		ungetc(c, json);
	}

	while(c != ']'){
		skip_ws(json);
		c = next_c(json);


		if(c != '{'){
			fprintf(stderr, "Error, line number %d; invalid object definition '%c'\n", lineNumber, c);
			fclose(json);
			exit(-1);
		}

		skip_ws(json);
		c = next_c(json);

		while(c != '}'){

			if(c == '"'){
				ungetc(c, json);
			}
			key = next_str(json);

			if(strcmp(key, "type") == 0){
				skip_ws(json);
				c = next_c(json);

				if(c != ':'){
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);

				}

				else{
					skip_ws(json);
					value = next_str(json);
					objects[items].type = value;
				}

			} 
			else if(strcmp(key, "width") == 0) {

				skip_ws(json);

				c = next_c(json);

				if(c != ':') {

					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(1);
				}
				
				else{
					skip_ws(json);
					objects[items].structures.camera.width = next_num(json);
				}

			} 
			else if(strcmp(key, "height") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':') {

					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(1);
				}

				else{
					skip_ws(json);
					objects[items].structures.camera.height = next_num(json);
				}

			} 
			else if(strcmp(key, "radius") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':') {

					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);

				}
				else{
					skip_ws(json);
					objects[items].structures.sphere.radius = next_num(json);
				}

			}
			else if(strcmp(key, "color") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':') {

					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(1);

				}
				else{
					skip_ws(json);
					vector = next_vec(json);

					if(strcmp(objects[items].type, "sphere") == 0){
						objects[items].structures.sphere.color[0] = vector[0];
						objects[items].structures.sphere.color[1] = vector[1];
						objects[items].structures.sphere.color[2] = vector[2];

					} 
					else if(strcmp(objects[items].type, "plane") == 0) {

						objects[items].structures.plane.color[0] = vector[0];
						objects[items].structures.plane.color[1] = vector[1];
						objects[items].structures.plane.color[2] = vector[2];
					}
				}

			} 
			else if(strcmp(key, "position") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':') {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(1);

				}
				else{
					skip_ws(json);
					vector = next_vec(json);

					if(strcmp(objects[items].type, "sphere") == 0){
						objects[items].structures.sphere.position[0] = vector[0];
						objects[items].structures.sphere.position[1] = vector[1];
						objects[items].structures.sphere.position[2] = vector[2];

					} 
					else if(strcmp(objects[items].type, "plane") == 0) {
						objects[items].structures.plane.position[0] = vector[0];
						objects[items].structures.plane.position[1] = vector[1];
						objects[items].structures.plane.position[2] = vector[2];
					}
				}

			} 
			else if(strcmp(key, "normal") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':') {
					fprintf(stderr, "Error, line number %d; unexpected character '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);

				}
				else{
					skip_ws(json);
					vector = next_vec(json);
					objects[items].structures.plane.normal[0] = vector[0];
					objects[items].structures.plane.normal[1] = vector[1];
					objects[items].structures.plane.normal[2] = vector[2];
				}

			}
			else{
				fprintf(stderr, "Error, line number %d; invalid type '%s'.\n", key);
				fclose(json);
				exit(-1);
			}

			skip_ws(json);
			c = next_c(json);

			if(c == ','){
				skip_ws(json);
                c = next_c(json);
			}
		}  // End of object parsing

		skip_ws(json);
		c = next_c(json);

		if(c == '{') {
			ungetc(c, json);
		}
		if(c == ','){
			skip_ws(json);
			c = next_c(json);

			if(c == '{') ungetc(c, json);
		}

													// Increment array index counter
		items += 1;

	} 												// End of JSON file

													// Return the items
	return items;
}