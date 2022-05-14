//
// This code was written by Paul Bourke (1989) and is available from 
// http://paulbourke.net/papers/triangulate/
//
#ifndef Delaunay_H
#define Delaunay_H

#include <cmath>
#include <iostream>
#include <stdlib.h> // for C qsort
#include <time.h>   // for random

const double EPSILON = 0.000001;

struct ITRIANGLE {
  int p1, p2, p3;
};

struct IEDGE {
  int p1, p2;
};

struct XYZ {
  double x, y, z;
};

int XYZCompare(const void* v1, const void* v2);
int Triangulate(int nv, XYZ pxyz[], ITRIANGLE v[], int& ntri);
int CircumCircle(double, double, double, double, double, double, double, double,
                 double&, double&, double&);

#endif
