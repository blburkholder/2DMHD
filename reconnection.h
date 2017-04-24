#include "MHD_2D.h"

class reconnection: public MHD_2D {
public:
  reconnection(float deltax, float deltay, int width, int height,\
    float deltat, int number_of_steps, int output, float gam, float alph, float bg,\
    float bs, float sig, float resistivity, int c);
  ~reconnection() { }
  void initialize_grid();
  void bound();
  void j_bound();
  void mag_bound();
};
