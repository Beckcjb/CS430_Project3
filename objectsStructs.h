#ifndef OBJECTSTRUCTS_H_INCLUDED
#define OBJECTSTRUCTS_H_INCLUDED
// object Camera struct 
typedef struct Camera
{
	double width;
	double height;
} Camera;

// Plane object
typedef struct Plane
{
	double diffuseColor[3];
	double specularColor[3];
    double position[3];
	double normal[3];

} Plane;

// Shpere object
typedef struct Sphere
{
	double diffuseColor[3];
	double specularColor[3];
    double position[3];
	double radius;

} Sphere;

// Light object
typedef struct Light
{
    double color[3];
    double position[3];
    double direction[3];
    double radial_a0, radial_a1, radial_a2;
    double angular_a0;
	double theta;

} Light;

typedef struct Object
{
	char *type;

	union structures
	{
		Camera camera;
        Light light;
		Plane plane;
		Sphere sphere;

	} structures;

} Object;

#endif // OBJECSTRUCTS_H_INCLUDED