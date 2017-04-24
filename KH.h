#include "MHD_2D.h"

class KH: public MHD_2D {
public:
  KH(float deltax, float deltay, int width, int height,\
    float deltat, int number_of_steps, int output, float gam, float alph, float bg,\
    float bs, float sig, float resistivity, int c);
  ~KH() { }
  void initialize_grid();
  void bound();
  void j_bound();
};
