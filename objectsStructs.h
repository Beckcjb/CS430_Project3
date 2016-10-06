
#ifndef OBJECTSTRUCTS_H_INCLUDED
#define OBJECTSTRUCTS_H_INCLUDED
// object Camer struct 
typedef struct Camera
{
	double width;
	double height;
} Camera;

typedef struct LightDiffuse{
	
	double color[3];
	double position[3];
	double direction[3];
	double radial_a2;
	double radial_a1;
	double radial_a0;
	double angular_a0;
} LightDiffuse;

typedef struct LightEmission{
	
	double color[3];
	double position[3];
	double direction[3];
	double radial_a2;
	double radial_a1;
	double radial_a0;
	double angular_a0;
} LightEmission;


typedef struct LightAmbient{
	
	double color[3];
	double position[3];
	double direction[3];
	double radial_a2;
	double radial_a1;
	double radial_a0;
	double angular_a0;
} LightAmbient;

// object plane struct
typedef struct Plane
{
	double color[3];
	double position[3];
	double normal[3];

} Plane;

// object sphere struct
typedef struct Sphere
{
	double diff_color[3];
	double specular_color[3];
	double position[3];
	double radius;

} Sphere;
// objects union
typedef struct Object
{
	char *type;

	union structures
	{
		Camera camera;
		LightAmbient lightAmbient;
		LightDiffuse lightDiffuse;
		LightEmission lightEmission;
		Plane plane;
		Sphere sphere;

	} structures;

} Object;

#endif // OBJECSTRUCTS_H_INCLUDED