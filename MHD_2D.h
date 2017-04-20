#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <algorithm>
#include "grid.h"

class MHD_2D {
public:

  MHD_2D (int init, float deltax, float deltay, int width, int height,\
    float deltat, int stps, int output, float gam, float alph,\
    float sig, float bs, float bg, float resistivity, int cl);
  ~MHD_2D ();
  void initialize_grid();
  void initialize_grid2();
  void initialize_grid_recon();
  void initialize_grid_KH();

  void integrate();

  float rho_f(int i) { return rho->get(i)*ux->get(i); }
  float rho_g(int i) { return rho->get(i)*uy->get(i); }

  float ux_f(int i) { return rho->get(i)*ux->get(i)*ux->get(i) +\
    0.5*(p->get(i) - bx->get(i)*bx->get(i) +\
      by->get(i)*by->get(i) + bz->get(i)*bz->get(i)); }
  float ux_g(int i) { return rho->get(i)*ux->get(i)*uy->get(i) -\
    bx->get(i)*by->get(i); }

  float uy_f(int i) { return rho->get(i)*ux->get(i)*uy->get(i) -\
    bx->get(i)*by->get(i); }
  float uy_g(int i) { return rho->get(i)*uy->get(i)*uy->get(i) +\
    0.5*(p->get(i) + bx->get(i)*bx->get(i) -\
      by->get(i)*by->get(i) + bz->get(i)*bz->get(i)); }

  float uz_f(int i) { return rho->get(i)*ux->get(i)*uz->get(i) -\
    bx->get(i)*bz->get(i); }
  float uz_g(int i) { return rho->get(i)*uy->get(i)*uz->get(i) -\
    by->get(i)*bz->get(i); }

  float bx_f(int i) { return 0; }
  float bx_g(int i) { return jz->get(i)*eta->get(i) +\
    uy->get(i)*bx->get(i) - ux->get(i)*by->get(i); }

  float by_f(int i) { return ux->get(i)*by->get(i) -\
    uy->get(i)*bx->get(i) - jz->get(i)*eta->get(i); }
  float by_g(int i) { return 0; }

  float bz_f(int i) { return eta->get(i)*jy->get(i) +\
    ux->get(i)*bz->get(i) - uz->get(i)*bx->get(i); }
  float bz_g(int i) { return uy->get(i)*bz->get(i) -\
    uz->get(i)*by->get(i) - eta->get(i)*jx->get(i); }

  float h_f(int i) { return h->get(i)*ux->get(i); }
  float h_g(int i) { return h->get(i)*uy->get(i); }

  float w_f(int i) { return ux->get(i)*(h->get(i) +\
    0.5*(p->get(i) + bx->get(i)*bx->get(i) +\
    by->get(i)*by->get(i) + bz->get(i)*bz->get(i))) -\
    bx->get(i)*(ux->get(i)*bx->get(i) +\
    uy->get(i)*by->get(i) + uz->get(i)*bz->get(i)) +\
    eta->get(i)*(jy->get(i)*bz->get(i) - jz->get(i)*by->get(i)); }
  float w_g(int i) { return uy->get(i)*(h->get(i) +\
    0.5*(p->get(i) + bx->get(i)*bx->get(i) +\
    by->get(i)*by->get(i) + bz->get(i)*bz->get(i))) -\
    by->get(i)*(ux->get(i)*bx->get(i) +\
    uy->get(i)*by->get(i) + uz->get(i)*bz->get(i)) +\
    eta->get(i)*(jz->get(i)*bx->get(i) - jx->get(i)*bz->get(i)); }

  //calculate resistivity and current density here
  //also do pressure
  void ampere();
  //smoothing is numerical diffusion
  void smooth();
  void bound();
  void fancy_bound();
  void mag_bound();
  void j_bound();
  void j_bound_fancy();
  void j_bound_KH();
  void bound_KH();
  void leap();
  void save_iteration();
  int index(int x, int y);
  void quit();

private:
  //dx,dy is grid spacing, ddx,ddy are grid lengths
  const int init_con;
  const float dx;
  const float dy;
  const int ddx;
  const int ddy;
  const int nx;
  const int ny;
  const float dt;
  int steps;
  const int out_steps;
  const float gamma;
  const float alpha;
  const float b_guide;
  float res;
  const int close;
  int here;
  int up;
  int down;
  int left;
  int right;

  grid * h;
  grid * rho;
  grid * ux;
  grid * uy;
  grid * uz;
  grid * bx;
  grid * by;
  grid * bz;

  grid *  p;
  grid * jx;
  grid * jy;
  grid * jz;
  grid * jj;
  grid * eta;

  grid * hv;
  grid * rhov;
  grid * uxv;
  grid * uyv;
  grid * uzv;
  grid * bxv;
  grid * byv;
  grid * bzv;

  grid * lp_h;
  grid * lp_rho;
  grid * lp_ux;
  grid * lp_uy;
  grid * lp_uz;
  grid * lp_bx;
  grid * lp_by;
  grid * lp_bz;

  //output files
  std::ofstream uxfile;
  std::ofstream uyfile;
  std::ofstream uzfile;
  std::ofstream bxfile;
  std::ofstream byfile;
  std::ofstream bzfile;
  std::ofstream rhofile;
  std::ofstream pfile;
};
