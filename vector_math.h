#ifndef VECTORMATH_H_INCLUDED
#define VECTORMATH_H_INCLUDED
typedef double* v3;

// Vector Addition
static inline void v3_add(v3 x, v3 y, v3 c){

    c[0] = x[0] + y[0];
    c[1] = x[1] + y[1];
    c[2] = x[2] + y[2];


}
// vector subtraction
static inline void  v3_subtract(v3 x, v3 y, v3 c){

    c[0] = x[0] - y[0];
    c[1] = x[1] - y[1];
    c[2] = x[2] - y[2];

}

// vector scale
static inline void v3_scale(v3 x, double s, v3 z){    

	v3 c;

    c[0] = x[0] * z[0];
    c[1] = x[1] * z[1];
    c[2] = x[2] * z[2];


}

// vector dot plot
static inline double v3_dot(v3 a, v3 b){
    double c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return c;
}

// vector cross section
static inline void v3_cross(v3 a, v3 b, v3 c){

    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];

}


#endif // VECTORMATH_H_INCLUDED