#include "MHD_2D.h"

class bomb: public MHD_2D {
public:
  bomb(float deltax, float deltay, int width, int height,\
    float deltat, int number_of_steps, int output, float gam, float alph, float bg,\
    float bs, float sig, float resistivity, int c);
  ~bomb() { }
  void initialize_grid();
  void bound();
  void j_bound();
};
