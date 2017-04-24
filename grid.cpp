#include "grid.h"

grid::~grid(){
  delete[] parameter;
  delete[] diff;
}

grid::grid(int x, int y, float sig):
  nx( x ),
  ny( y ),
  sigma( sig ),
  parameter ( new float[(nx+1)*(ny+1)] ),
  diff ( new float[(nx+1)*(ny+1)] )  {}

void grid::set(int index, float value) {
  parameter[index] = value;
}

float grid::get(int index) {
  return parameter[index];
}

void grid::d_set(int index, float value) {
  diff[index] = value;
}

float grid::d_get(int index) {
  return diff[index];
}

void grid::smooth() {
  //smoothing part 1
  for (int j = 1; j <= ny-1; j++) {
    for (int i = 1; i <= nx-1; i++) {
      d_set(j*(nx+1) + i,sigma*parameter[j*(nx+1) + i] +\
        ((1.0-sigma)/4.0)*(parameter[(j+1)*(nx+1) + i+1] +\
        parameter[(j-1)*(nx+1) + i+1] + parameter[(j+1)*(nx+1) + i-1] +\
        parameter[(j-1)*(nx+1) + i-1]));
    }
  }

  for (int j = 1; j <= ny-1; j++) {
    for (int i = 1; i <= nx-1; i++) {
      set(j*(nx+1) + i,d_get(j*(nx+1) + i));
    }
  }
}

//I wrote a copy constructor yay
/*
grid::grid(const grid &obj) :
  dx(obj.dx),
  dy(obj.dy),
  ddx(obj.ddx),
  ddy(obj.ddy),
  sigma(obj.sigma),
  parameter( new float[(int)(ddx*ddy/(dx*dy))] ),
  grid_150( new float[(int)(ddx*ddy/(dx*dy))] ),
  diff( new float[(int)(ddx*ddy/(dx*dy))] ) {
    for (int i = 0; i <= ddx; i++) {
      for (int j = 0; j <= ddy; j++) {
	parameter[(int)(j*(nx) + i)] = obj.parameter[(int)(j*(nx) + i)];
	grid_150[(int)(j*(nx) + i)] = obj.grid_150[(int)(j*(nx) + i)];
	diff[(int)(j*(nx) + i)] = obj.diff[(int)(j*(nx) + i)];
      }
    }
  }
*/

//Ill overload your operator=
grid& grid::operator= (grid &g) {
  for (int j = 0; j <= ny; j++) { 
    for (int i = 0; i <= nx; i++) {
      set(j*(nx+1) + i, g.get(j*(nx+1) + i));
    }
  }
  return *this;
}

