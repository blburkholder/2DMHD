#include "MHD_2D.h"
MHD_2D::~MHD_2D(){}

//energy2() sucks?
//more smoothing for pressure than velocity
//current dependent resistivity
MHD_2D::MHD_2D(int init, float deltax, float deltay, int width, int height,\
    float deltat, int number_of_steps, int output, float gam, float alph, float bg,\
    float sig, float resistivity, int cl) :
  init_con( init ),
  dx( deltax ),
  dy( deltay ),
  ddx( width ),
  ddy( height ),
  nx (ddx/dx ),
  ny ( ddy/dy ),
  dt( deltat ),
  steps( number_of_steps ),
  out_steps( output ),
  gamma( gam ),
  alpha( alph ),
  b_guide( bg ),
  res( resistivity ),
  close ( cl ),
  here ( 0 ),
  up ( 0 ),
  down ( 0 ),
  left ( 0 ),
  right ( 0 ),
  h( new grid(nx,ny,sig) ),
  rho ( new grid(nx,ny,sig) ),
  ux ( new grid(nx,ny,sig) ),
  uy ( new grid(nx,ny,sig) ),
  uz ( new grid(nx,ny,sig) ),
  bx ( new grid(nx,ny,sig) ),
  by ( new grid(nx,ny,sig) ),
  bz ( new grid(nx,ny,sig) ),

  p ( new grid(nx,ny,sig) ),
  jx ( new grid(nx,ny,sig) ),
  jy ( new grid(nx,ny,sig) ),
  jz ( new grid(nx,ny,sig) ),
  jj ( new grid(nx,ny,sig) ),
  eta ( new grid(nx,ny,sig) ),

  hv( new grid(nx,ny,sig) ),
  rhov ( new grid(nx,ny,sig) ),
  uxv ( new grid(nx,ny,sig) ),
  uyv ( new grid(nx,ny,sig) ),
  uzv ( new grid(nx,ny,sig) ),
  bxv ( new grid(nx,ny,sig) ),
  byv ( new grid(nx,ny,sig) ),
  bzv ( new grid(nx,ny,sig) ),

  lp_h( new grid(nx,ny,sig) ),
  lp_rho ( new grid(nx,ny,sig) ),
  lp_ux ( new grid(nx,ny,sig) ),
  lp_uy ( new grid(nx,ny,sig) ),
  lp_uz ( new grid(nx,ny,sig) ),
  lp_bx ( new grid(nx,ny,sig) ),
  lp_by ( new grid(nx,ny,sig) ),
  lp_bz ( new grid(nx,ny,sig) )
  { uxfile.open("ux.txt");
    uyfile.open("uy.txt");
    uzfile.open("uz.txt");
    bxfile.open("bx.txt");
    byfile.open("by.txt");
    bzfile.open("bz.txt");
    rhofile.open("rho.txt");
    pfile.open("p.txt"); 
    if ( init_con == 1 ) initialize_grid();
    else if ( init_con == 2 )initialize_grid2();
    else if ( init_con == 3 ) initialize_grid_recon();
    else initialize_grid_KH();
 }
//plasma bomb
//secondary closure relation not implemented
void MHD_2D::initialize_grid2() {
  for (int j = 0; j <= ny; j++) {
    for (int i = 0; i <= nx; i++) {
      here = index(i,j);
      rho->set(here,0.5);
      ux->set(here,0);
      uy->set(here,0);
      uz->set(here,0);

      if (pow(i - 0.5*nx,2) + pow(j - 0.5*ny,2) < 0.5*nx) {
        p->set(here,0.5 + exp(-(pow(i - 0.5*nx,2) + pow(j - 0.5*ny,2))/10));
      } 
      else p->set(here,0.5);
      h->set(here,pow(p->get(here)/2,1/gamma));
      bx->set(here,0);      
      by->set(here,1);
      bz->set(here,0);
      eta->set(here,res);

      jx->set(here,(0.5/dy)*(bz->get(index(i,j+1)) - bz->get(index(i,j-1))));
      jy->set(here,(-0.5/dx)*(bz->get(index(i+1,j)) - bz->get(index(i-1,j))));
      jz->set(here,(-0.5/dy)*(bx->get(index(i,j+1)) - bx->get(index(i,j-1))) +\
        (0.5/dx)*(by->get((index(i+1,j))) - by->get((index(i-1,j)))));
      jj->set(here,jx->get(here)*jx->get(here) + jy->get(here)*jy->get(here) + jz->get(here)*jz->get(here));
      j_bound();

      lp_h->set(here,h->get(here));
      lp_rho->set(here,rho->get(here));
      lp_ux->set(here,ux->get(here));
      lp_uy->set(here,uy->get(here));
      lp_uz->set(here,uz->get(here));
      lp_bx->set(here,bx->get(here));      
      lp_by->set(here,by->get(here));
      lp_bz->set(here,bz->get(here));

      hv->set(here,h->get(here));
      rhov->set(here,rho->get(here));
      uxv->set(here,ux->get(here));
      uyv->set(here,uy->get(here));
      uzv->set(here,uz->get(here));
      bxv->set(here,bx->get(here));
      byv->set(here,by->get(here));
      bzv->set(here,bz->get(here));
    }
  }
  bound(); 
  save_iteration();
}
//alfven wave
//secondary closure relation not implemented
void MHD_2D::initialize_grid() {
  for (int j = 0; j <= ny; j++) {
    for (int i = 0; i <= nx; i++) {
      here = index(i,j),
      rho->set(here,0.5);
      ux->set(here,0);
      uy->set(here,0);

      if (pow(i - 0.5*nx,2) + pow(j - 0.5*ny,2) < 0.5*nx) {
        uz->set(here,exp(-(pow(i - 0.5*nx,2) + pow(j - 0.5*ny,2))/10));
      } 
      else uz->set(here,0);
      p->set(here,0.5);
      h->set(here,pow(p->get(here)/2,1/gamma));
      bx->set(here,1);      
      by->set(here,0);
      bz->set(here,0);
      eta->set(here,res);

      jx->set(here,(0.5/dy)*(bz->get(index(i,j+1)) - bz->get(index(i,j-1))));
      jy->set(here,(-0.5/dx)*(bz->get(index(i+1,j)) - bz->get(index(i-1,j))));
      jz->set(here,(-0.5/dy)*(bx->get(index(i,j+1)) - bx->get(index(i,j-1))) +\
        (0.5/dx)*(by->get((index(i+1,j))) - by->get((index(i-1,j)))));
      jj->set(here,jx->get(here)*jx->get(here) + jy->get(here)*jy->get(here) + jz->get(here)*jz->get(here));
      j_bound();

      lp_h->set(here,h->get(here));
      lp_rho->set(here,rho->get(here));
      lp_ux->set(here,ux->get(here));
      lp_uy->set(here,uy->get(here));
      lp_uz->set(here,uz->get(here));
      lp_bx->set(here,bx->get(here));      
      lp_by->set(here,by->get(here));
      lp_bz->set(here,bz->get(here));

      hv->set(here,0);
      rhov->set(here,0);
      uxv->set(here,0);
      uyv->set(here,0);
      uzv->set(here,0);
      bxv->set(here,0);
      byv->set(here,0);
      bzv->set(here,0);
    }
  }
  bound(); 
  save_iteration();
}

void MHD_2D::initialize_grid_recon() {
  float xd = 70;
  float k1 = 2;
  float k2 = 1;
  float p_inf = 0.2;
  float px = 0;
  float c2 = (1 - p_inf)/(k2*pow(2,k1+1)/k1 - 1 - tanh(k2));
  float c1 = 1 - p_inf + c2*(1 + tanh(k2));
  float y = 0;

  for (int j = 0; j <= ny; j++) {
    for (int i = 0; i <= nx; i++) {
      here = index(i,j);
      //vertical dimension is symmetric about axis
      y = dy*j - ddy/2;
      //pressure in the middle
      px = c1*pow(1+i*dx/xd,-k1) + c2*tanh(k2*(i*dx/xd - 1)) - c2 + p_inf;

      p->set(here,0.25 + (1/(cosh(y*sqrt(px))*cosh(y*sqrt(px))))*px);
      rho->set(here,p->get(here));

      bx->set(here,-sqrt(px)*tanh(y*sqrt(px)));
      by->set(here,(-0.5/px)*(1 - y*sqrt(px)*tanh(y*sqrt(px)))*\
        (c2*k2*(1/(cosh(k2*(i*dx/xd-1))*cosh(k2*(i*dx/xd-1))))/xd -\
        k1*c1*pow(1+i*dx/xd,-k1-1)/xd));
      bz->set(here,b_guide);

      ux->set(here,0);
      uy->set(here,0);
      uz->set(here,0);
      //leap quantities are the same for first iteration
      lp_rho->set(here,rho->get(here));
      lp_ux->set(here,0);
      lp_uy->set(here,0);
      lp_uz->set(here,0);
      lp_bx->set(here,bx->get(here));
      lp_by->set(here,by->get(here));
      lp_bz->set(here,bz->get(here));

      if ( close == 1 ) lp_h->set(here,pow(p->get(here)/2.0,1.0/gamma));
      else lp_h->set(here,0.5*rho->get(here)*(ux->get(here)*ux->get(here) +\
        uy->get(here)*uy->get(here) + uz->get(here)*uz->get(here)) + \
        0.5*(bx->get(here)*bx->get(here) + by->get(here)*by->get(here) +\
        bz->get(here)*bz->get(here)) + 0.5*p->get(here)/(gamma-1));

      h->set(here,lp_h->get(here));
      eta->set(here,res);
      jx->set(here,(0.5/dy)*(bz->get(index(i,j+1)) - bz->get(index(i,j-1))));
      jy->set(here,(-0.5/dx)*(bz->get(index(i+1,j)) - bz->get(index(i-1,j))));
      jz->set(here,(-0.5/dy)*(bx->get(index(i,j+1)) - bx->get(index(i,j-1))) +\
        (0.5/dx)*(by->get((index(i+1,j))) - by->get((index(i-1,j)))));
      jj->set(here,jx->get(here)*jx->get(here) + jy->get(here)*jy->get(here) + jz->get(here)*jz->get(here));
      j_bound_fancy();

      hv->set(here,h->get(here));
      rhov->set(here,rho->get(here));
      uxv->set(here,ux->get(here));
      uyv->set(here,uy->get(here));
      uzv->set(here,uz->get(here));
      bxv->set(here,bx->get(here));
      byv->set(here,by->get(here));
      bzv->set(here,bz->get(here));
    }
  }
  fancy_bound();
  save_iteration();
}

void MHD_2D::initialize_grid_KH() {
  float y = 0;
  for (int j = 0; j <= ny; j++) {
    for (int i = 0; i <= nx; i++) {
      here = index(i,j);
      y = dy*j - ddy/2;     
      
      rho->set(here,1+0.1*tanh(y));
      ux->set(here,tanh(y));
      uy->set(here,-0.5*exp(-(float)(16*y*y))*sin((float)i*dx));
      uz->set(here,0);

      bx->set(here,0);      
      by->set(here,0);
      bz->set(here,tanh(i*dx));
      eta->set(here,res);
      p->set(here,1.5-bz->get(here)*bz->get(here));
      h->set(here,pow(p->get(here)/2,1/gamma));

      jx->set(here,(0.5/dy)*(bz->get(index(i,j+1)) - bz->get(index(i,j-1))));
      jy->set(here,(-0.5/dx)*(bz->get(index(i+1,j)) - bz->get(index(i-1,j))));
      jz->set(here,(-0.5/dy)*(bx->get(index(i,j+1)) - bx->get(index(i,j-1))) +\
        (0.5/dx)*(by->get((index(i+1,j))) - by->get((index(i-1,j)))));
      jj->set(here,jx->get(here)*jx->get(here) +\
        jy->get(here)*jy->get(here) + jz->get(here)*jz->get(here));
      bound_KH();

      lp_h->set(here,h->get(here));
      lp_rho->set(here,rho->get(here));
      lp_ux->set(here,ux->get(here));
      lp_uy->set(here,uy->get(here));
      lp_uz->set(here,uz->get(here));
      lp_bx->set(here,bx->get(here));      
      lp_by->set(here,by->get(here));
      lp_bz->set(here,bz->get(here));

      hv->set(here,h->get(here));
      rhov->set(here,rho->get(here));
      uxv->set(here,ux->get(here));
      uyv->set(here,uy->get(here));
      uzv->set(here,uz->get(here));
      bxv->set(here,bx->get(here));
      byv->set(here,by->get(here));
      bzv->set(here,bz->get(here));
    }
  }
}

int MHD_2D::index(int x, int y) {
  return y*(nx+1) + x;
}

void MHD_2D::leap() {
  for (int i = 1; i <= steps; i++) {
    integrate();

    lp_rho = rho;
    lp_ux = ux;
    lp_uy = uy;
    lp_uz = uz;
    lp_bx = bx;
    lp_by = by;
    lp_bz = bz;
    lp_h = h;

    if ( init_con < 3 ) bound();
    else if ( init_con == 3) {
      mag_bound();
      fancy_bound();
    }
    else bound_KH();
    //maybe no smoothing for KH?
    smooth(i);
    //do this to make sure del dot b = 0 after smoothing too
    if ( init_con == 3 ) mag_bound(); 
    ampere();
    if (steps == 0) std::cout << "quitting at iteration" << i << std::endl;

    rho = rhov;
    ux = uxv;
    uy = uyv;
    uz = uzv;
    bx = bxv;
    by = byv;
    bz = bzv;
    h = hv;

    if ((i % out_steps) == 0) {
      std::cout << "integration step:" << i << std::endl;
      save_iteration();
    }
  }
}

void MHD_2D::integrate() {
  for (int j = 1; j <= ny-1; j++) {
    for (int i = 1; i <= nx-1; i++) {
      here = index(i,j);
      up = index(i,j+1);
      down = index(i,j-1);
      left = index(i-1,j);
      right = index(i+1,j);
//continuity
      rhov->set(here,lp_rho->get(here) - (dt/(2*dx))*(rho_f(right) -\
        rho_f(left)) - (dt/(2*dy))*(rho_g(up) - rho_g(down)));
//momentum
      uxv->set(here,(1/rhov->get(here))*(lp_rho->get(here)*lp_ux->get(here) -\
        (dt/(2*dx))*(ux_f(right) - ux_f(left)) - (dt/(2*dy))*(ux_g(up) - ux_g(down))));
      uyv->set(here,(1/rhov->get(here))*(lp_rho->get(here)*lp_uy->get(here) -\
        (dt/(2*dx))*(uy_f(right) - uy_f(left)) - (dt/(2*dy))*(uy_g(up) - uy_g(down))));
      uzv->set(here,(1/rhov->get(here))*(lp_rho->get(here)*lp_uz->get(here) -\
        (dt/(2*dx))*(uz_f(right) - uz_f(left)) - (dt/(2*dy))*(uz_g(up) - uz_g(down))));
//faraday
      bxv->set(here,lp_bx->get(here) - (dt/(2*dx))*(bx_f(right) -\
        bx_f(left)) - (dt/(2*dy))*(bx_g(up) - bx_g(down)));  
      byv->set(here,lp_by->get(here) - (dt/(2*dx))*(by_f(right) -\
        by_f(left)) - (dt/(2*dy))*(by_g(up) - by_g(down)));
      bzv->set(here,lp_bz->get(here) - (dt/(2*dx))*(bz_f(right) -\
        bz_f(left)) - (dt/(2*dy))*(bz_g(up) - bz_g(down)));
//energy
      if ( close == 1 ) {
        hv->set(here,lp_h->get(here) - (dt/(2*dx))*(h_f(right) -\
          h_f(left)) - (dt/(2*dy))*(h_g(up) - h_g(down)) +\
          ((gamma-1)/gamma)*(pow(h->get(here),1-gamma)*eta->get(here)*\
          jj->get(here)));
      }
      else {
        hv->set(here,lp_h->get(here) - (dt/(2*dx))*(w_f(right) -\
          w_f(left)) - (dt/(2*dy))*(w_g(up) - w_g(down)));
      }
    }
  }
}

void MHD_2D::ampere() {
  if ( close == 1 ) {
    for (int j = 0; j <= ny; j++) {
      for (int i = 0; i <= nx; i++) {
        here = index(i,j);
        p->set(here, pow(hv->get(here),gamma)*2);
        if isnan(p->get(here)) quit();
      }
    }
  }
  else {
    for (int j = 0; j <= ny; j++) {
      for (int i = 0; i <= nx; i++) {
        here = index(i,j);
        p->set(here, 2*(gamma-1)*(hv->get(here) -\
          0.5*rho->get(here)*(ux->get(here)*ux->get(here) +\
          uy->get(here)*uy->get(here) + uz->get(here)*uz->get(here)) -\
          0.5*(bx->get(here)*bx->get(here) +\
          by->get(here)*by->get(here) + bz->get(here)*bz->get(here))));
        if isnan(p->get(here)) quit();
      }
    }
  }
  for (int j = 1; j <= ny-1; j++) {
    for (int i = 1; i <= nx-1; i++) {
      here = index(i,j);
      up = index(i,j+1);
      down = index(i,j-1);
      left = index(i-1,j);
      right = index(i+1,j);
      jx->set(here,(0.5/dy)*(bzv->get(up) - bzv->get(down)));
      jy->set(here,(-0.5/dx)*(bzv->get(right) - bzv->get(left)));
      jz->set(here,(-0.5/dy)*(bxv->get(up) - bxv->get(down)) +\
        (0.5/dx)*(byv->get(right) - byv->get(left)));
      jj->set(here,jx->get(here)*jx->get(here) +\
        jy->get(here)*jy->get(here) + jz->get(here)*jz->get(here)); 
    }
  }
  if ( init_con < 3 ) j_bound();
  else if (init_con == 3 ) j_bound_fancy();
}

void MHD_2D::j_bound() {
  for (int i = 1; i <= nx-1; i++) {
    jx->set(index(i,0),jx->get(index(i,ny-1)));
    jy->set(index(i,0),jy->get(index(i,ny-1)));
    jz->set(index(i,0),jz->get(index(i,ny-1)));
    jj->set(index(i,0),jj->get(index(i,ny-1)));

    jx->set(index(i,ny),jx->get(index(i,1)));
    jy->set(index(i,ny),jy->get(index(i,1)));
    jz->set(index(i,ny),jz->get(index(i,1)));
    jj->set(index(i,ny),jj->get(index(i,1)));
  }
  for (int j = 1; j <= ny-1; j++) {
    jx->set(index(0,j),jx->get(index(nx-1,j)));
    jy->set(index(0,j),jy->get(index(nx-1,j)));
    jz->set(index(0,j),jz->get(index(nx-1,j)));
    jj->set(index(0,j),jj->get(index(nx-1,j)));

    jx->set(index(nx,j),jx->get(index(1,j)));
    jy->set(index(nx,j),jy->get(index(1,j)));
    jz->set(index(nx,j),jz->get(index(1,j)));
    jj->set(index(nx,j),jj->get(index(1,j)));
  }
}

void MHD_2D::bound() {
  for (int i = 1; i <= nx-1; i++) {
    rhov->set(index(i,0),rhov->get(index(i,ny-1)));
    hv->set(index(i,0),hv->get(index(i,ny-1)));
    uxv->set(index(i,0),uxv->get(index(i,ny-1)));
    uyv->set(index(i,0),uyv->get(index(i,ny-1)));
    uzv->set(index(i,0),uzv->get(index(i,ny-1)));
    bxv->set(index(i,0),bxv->get(index(i,ny-1)));
    byv->set(index(i,0),byv->get(index(i,ny-1)));
    bzv->set(index(i,0),bzv->get(index(i,ny-1)));

    rhov->set(index(i,ny),rhov->get(index(i,1)));
    hv->set(index(i,ny),hv->get(index(i,1)));
    uxv->set(index(i,ny),uxv->get(index(i,1)));
    uyv->set(index(i,ny),uyv->get(index(i,1)));
    uzv->set(index(i,ny),uzv->get(index(i,1)));
    bxv->set(index(i,ny),bxv->get(index(i,1)));
    byv->set(index(i,ny),byv->get(index(i,1)));
    bzv->set(index(i,ny),bzv->get(index(i,1)));
  }
  for (int j = 1; j <= ny-1; j++) {
    rhov->set(index(0,j),rhov->get(index(nx-1,j)));
    hv->set(index(0,j),hv->get(index(nx-1,j)));
    uxv->set(index(0,j),uxv->get(index(nx-1,j)));
    uyv->set(index(0,j),uyv->get(index(nx-1,j)));
    uzv->set(index(0,j),uzv->get(index(nx-1,j)));
    bxv->set(index(0,j),bxv->get(index(nx-1,j)));
    byv->set(index(0,j),byv->get(index(nx-1,j)));
    bzv->set(index(0,j),bzv->get(index(nx-1,j)));

    rhov->set(index(nx,j),rhov->get(index(1,j)));
    hv->set(index(nx,j),hv->get(index(1,j)));
    uxv->set(index(nx,j),uxv->get(index(1,j)));
    uyv->set(index(nx,j),uyv->get(index(1,j)));
    uzv->set(index(nx,j),uzv->get(index(1,j)));
    bxv->set(index(nx,j),bxv->get(index(1,j)));
    byv->set(index(nx,j),byv->get(index(1,j)));
    bzv->set(index(nx,j),bzv->get(index(1,j)));
  }
}

void MHD_2D::j_bound_fancy() {
  for (int i = 1; i <= nx-1; i++) {
    jx->set(index(i,0),jx->get(index(i,2)) - alpha*(jx->get(index(i,3)) - jx->get(index(i,1))));
    jy->set(index(i,0),jy->get(index(i,2)) - alpha*(jy->get(index(i,3)) - jy->get(index(i,1))));
    jz->set(index(i,0),jz->get(index(i,2)) - alpha*(jz->get(index(i,3)) - jz->get(index(i,1))));
    jj->set(index(i,0),jx->get(index(i,0))*jx->get(index(i,0)) +\
      jy->get(index(i,0))*jy->get(index(i,0)) +\
      jz->get(index(i,0))*jz->get(index(i,0)));
    jx->set(index(i,ny),jx->get(index(i,ny-2)) + alpha*(jx->get(index(i,ny-1)) - jx->get(index(i,ny-3))));
    jy->set(index(i,ny),jy->get(index(i,ny-2)) + alpha*(jy->get(index(i,ny-1)) - jy->get(index(i,ny-3))));
    jz->set(index(i,ny),jz->get(index(i,ny-2)) + alpha*(jz->get(index(i,ny-1)) - jz->get(index(i,ny-3))));
    jj->set(index(i,ny),jx->get(index(i,ny))*jx->get(index(i,ny)) +\
      jy->get(index(i,ny))*jy->get(index(i,ny)) +\
      jz->get(index(i,ny))*jz->get(index(i,ny)));
  }
  for (int j = 1; j <= ny-1; j++) {
    jx->set(index(0,j),jx->get(index(2,j)) - alpha*(jx->get(index(3,j)) - jx->get(index(1,j))));
    jy->set(index(0,j),jy->get(index(2,j)) - alpha*(jy->get(index(3,j)) - jy->get(index(1,j))));
    jz->set(index(0,j),jz->get(index(2,j)) - alpha*(jz->get(index(3,j)) - jz->get(index(1,j))));
    jj->set(index(0,j),jx->get(index(0,j))*jx->get(index(0,j)) +\
      jy->get(index(0,j))*jy->get(index(0,j)) +\
      jz->get(index(0,j))*jz->get(index(0,j)));
    jx->set(index(nx,j),jx->get(index(nx-2,j)) + alpha*(jx->get(index(nx-1,j)) - jx->get(index(nx-3,j))));
    jy->set(index(nx,j),jy->get(index(nx-2,j)) + alpha*(jy->get(index(nx-1,j)) - jy->get(index(nx-3,j))));
    jz->set(index(nx,j),jz->get(index(nx-2,j)) + alpha*(jz->get(index(nx-1,j)) - jz->get(index(nx-3,j))));
    jj->set(index(nx,j),jx->get(index(nx,j))*jx->get(index(nx,j)) +\
      jy->get(index(nx,j))*jy->get(index(nx,j)) +\
      jz->get(index(nx,j))*jz->get(index(nx,j)));
  }
}

void MHD_2D::fancy_bound() {
  for (int j = 1; j <= ny-1; j++) {
  //left
    here = index(0,j);
    right = index(2,j);
    rhov->set(here,rhov->get(right) - alpha*(rhov->get(index(3,j)) - rhov->get(index(1,j))));
    hv->set(here,hv->get(right) - alpha*(hv->get(index(3,j)) - hv->get(index(1,j))));
    uxv->set(here,-uxv->get(right));
    uyv->set(here,-uyv->get(right));
    uzv->set(here,-uzv->get(right));
  //right
    here = index(nx,j);
    left = index(nx-2,j);
    rhov->set(here,rhov->get(left) + alpha*(rhov->get(index(nx-1,j)) - rhov->get(index(nx-3,j))));
    hv->set(here,hv->get(left) + alpha*(hv->get(index(nx-1,j)) - hv->get(index(nx-3,j))));
    uxv->set(here,uxv->get(left));
    uyv->set(here,uyv->get(left) + alpha*(uyv->get(index(nx-1,j)) - uyv->get(index(nx-3,j))));
    uzv->set(here,uzv->get(left) + alpha*(uzv->get(index(nx-1,j)) - uzv->get(index(nx-3,j))));
  }

  for (int i = 1; i <= nx-1; i++) {
  //bottom 
    here = index(i,0);
    up = index(i,2);
    rhov->set(here,rhov->get(up));
    hv->set(here,hv->get(up));
    uxv->set(here,uxv->get(up));
    uyv->set(here,-uyv->get(up));
    uzv->set(here,uzv->get(up) - alpha*(uzv->get(index(i,3)) - uzv->get(index(i,1))));
  //top
    here = index(i,ny);
    down = index(i,ny-2);
    rhov->set(here,rhov->get(down));
    hv->set(here,hv->get(down));
    uxv->set(here,uxv->get(down));
    uyv->set(here,-uyv->get(down));
    uzv->set(here,uzv->get(down) + alpha*(uzv->get(index(i,ny-1)) - uzv->get(index(i,ny-3))));
  }
    //dont forget the corners since they are included in smoothing!
    rhov->set(0,0.5*(rhov->get(index(0,1)) + rhov->get(index(1,0))));
    rhov->set((nx+1)*ny,0.5*(rhov->get(index(0,ny-1)) + rhov->get(index(1,ny))));
    rhov->set(nx,0.5*(rhov->get(index(nx-1,0)) + rhov->get(index(nx,1))));
    rhov->set((nx+1)*ny + nx,0.5*(rhov->get(index(nx,ny-1)) + rhov->get(index(nx-1,ny))));

    hv->set(0,0.5*(hv->get(index(0,1)) + hv->get(index(1,0))));
    hv->set((nx+1)*ny,0.5*(hv->get(index(0,ny-1)) + hv->get(index(1,ny))));
    hv->set(nx,0.5*(hv->get(index(nx-1,0)) + hv->get(index(nx,1))));
    hv->set((nx+1)*ny + nx,0.5*(hv->get(index(nx,ny-1)) + hv->get(index(nx-1,ny))));

    uxv->set(0,0.5*(uxv->get(index(0,1)) + uxv->get(index(1,0))));
    uxv->set((nx+1)*ny,0.5*(uxv->get(index(0,ny-1)) + uxv->get(index(1,ny))));
    uxv->set(nx,0.5*(uxv->get(index(nx-1,0)) + uxv->get(index(nx,1))));
    uxv->set((nx+1)*ny + nx,0.5*(uxv->get(index(nx,ny-1)) + uxv->get(index(nx-1,ny))));

    uyv->set(0,0.5*(uyv->get(index(0,1)) + uyv->get(index(1,0))));
    uyv->set((nx+1)*ny,0.5*(uyv->get(index(0,ny-1)) + uyv->get(index(1,ny))));
    uyv->set(nx,0.5*(uyv->get(index(nx-1,0)) + uyv->get(index(nx,1))));
    uyv->set((nx+1)*ny + nx,0.5*(uyv->get(index(nx,ny-1)) + uyv->get(index(nx-1,ny))));

    uzv->set(index(0,0),0.5*(uzv->get(index(0,1)) + uzv->get(index(1,0))));
    uzv->set((nx+1)*ny,0.5*(uzv->get(index(0,ny-1)) + uzv->get(index(1,ny))));
    uzv->set(nx,0.5*(uzv->get(index(nx-1,0)) + uzv->get(index(nx,1))));
    uzv->set((nx+1)*ny + nx,0.5*(uzv->get(index(nx,ny-1)) + uzv->get(index(nx-1,ny))));

}

void MHD_2D::mag_bound() {
  bx = bxv;
  by = byv;

  for (int j = 1; j <= ny-1; j++) {
    //left
    here = index(0,j);
    right = index(2,j);
    byv->set(here,byv->get(right));
    bxv->set(here,bx->get(right) + (dx/dy)*(by->get(index(1,j+1)) - by->get(index(1,j-1))));
    bzv->set(here,bzv->get(right) - alpha*(bzv->get(index(3,j)) - bzv->get(index(1,j))));
    //right
    here = index(nx,j);
    left = index(nx-2,j);
    byv->set(here,byv->get(left) + alpha*(byv->get(index(nx-1,j)) - byv->get(index(nx-3,j))));
    bxv->set(here,bx->get(left) - (dx/dy)*(by->get(index(nx-1,j+1)) - by->get(index(nx-1,j-1))));
    bzv->set(here,bzv->get(left) + alpha*(bzv->get(index(nx-1,j)) - bzv->get(index(nx-3,j))));
  }

  for (int i = 1; i <= nx-1; i++) {
    //bottom
    here = index(i,0);
    up = index(i,2);
    bxv->set(here,bxv->get(up));
    byv->set(here,by->get(up) + (dy/dx)*(bx->get(index(i+1,1)) - bx->get(index(i-1,1))));
    bzv->set(here,bzv->get(up) - alpha*(bzv->get(index(i,3)) - bzv->get(index(i,1)))); 
    //top
    here = index(i,ny);
    down = index(i,ny-2);
    bxv->set(here,bxv->get(down));
    byv->set(here,by->get(down) - (dy/dx)*(bx->get(index(i+1,ny-1)) - bx->get(index(i-1,ny-1))));
    bzv->set(here,bzv->get(down) + alpha*(bzv->get(index(i,ny-1)) - bzv->get(index(i,ny-3))));
  }
}

void MHD_2D::bound_KH() {
  for (int j = 1; j <= ny-1; j++) {
    //left
    rhov->set(index(0,j),rhov->get(index(nx-1,j)));
    hv->set(index(0,j),hv->get(index(nx-1,j)));
    uxv->set(index(0,j),uxv->get(index(nx-1,j)));
    uyv->set(index(0,j),uyv->get(index(nx-1,j)));
    uzv->set(index(0,j),uzv->get(index(nx-1,j)));
    bxv->set(index(0,j),bxv->get(index(nx-1,j)));
    byv->set(index(0,j),byv->get(index(nx-1,j)));
    bzv->set(index(0,j),bzv->get(index(nx-1,j)));
    //right
    rhov->set(index(nx,j),rhov->get(index(1,j)));
    hv->set(index(nx,j),hv->get(index(1,j)));
    uxv->set(index(nx,j),uxv->get(index(1,j)));
    uyv->set(index(nx,j),uyv->get(index(1,j)));
    uzv->set(index(nx,j),uzv->get(index(1,j)));
    bxv->set(index(nx,j),bxv->get(index(1,j)));
    byv->set(index(nx,j),byv->get(index(1,j)));
    bzv->set(index(nx,j),bzv->get(index(1,j)));
  }
  for (int i = 1; i <= nx-1; i++) {
    //bottom
    rhov->set(index(i,0),rhov->get(index(i,2)) - alpha*(rhov->get(index(i,3)) - rhov->get(index(i,1))));
    hv->set(index(i,0),hv->get(index(i,2)) - alpha*(hv->get(index(i,3)) - hv->get(index(i,1))));
    uxv->set(index(i,0),uxv->get(index(i,2)));
    uyv->set(index(i,0),-uyv->get(index(i,2)));
    uzv->set(index(i,0),-uzv->get(index(i,2)));
    bxv->set(index(i,0),bxv->get(index(i,2)));
    byv->set(index(i,0),byv->get(index(i,2)) + (dy/dx)*(bxv->get(index(i+1,1)) - bxv->get(index(i-1,1))));
    bzv->set(index(i,0),bzv->get(index(i,2)) - alpha*(bzv->get(index(i,3)) - bzv->get(index(i,1)))); 
    //top
    rhov->set(index(i,ny),rhov->get(index(i,ny-2)) + alpha*(rhov->get(index(i,ny-1)) - rhov->get(index(i,ny-3))));
    hv->set(index(i,ny),hv->get(index(i,ny-2)) + alpha*(hv->get(index(i,ny-1)) - hv->get(index(i,ny-3))));
    uxv->set(index(i,ny),uxv->get(index(i,ny-2)));
    uyv->set(index(i,ny),-uyv->get(index(i,ny-2)));
    uzv->set(index(i,ny),-uzv->get(index(i,ny-2)));
    bxv->set(index(i,ny),bxv->get(index(i,ny-2)));
    byv->set(index(i,ny),byv->get(index(i,ny-2)) - (dy/dx)*(bxv->get(index(i+1,ny-1)) - bxv->get(index(i-1,ny-1))));
    bzv->set(index(i,ny),bzv->get(index(i,ny-2)) + alpha*(bzv->get(index(i,ny-1)) - bzv->get(index(i,ny-3))));;
  }
  for (int i = 1; i <= nx-1; i++) {
    jx->set(index(i,0),jx->get(index(i,2)) - alpha*(jx->get(index(i,3)) - jx->get(index(i,1))));
    jy->set(index(i,0),jy->get(index(i,2)) - alpha*(jy->get(index(i,3)) - jy->get(index(i,1))));
    jz->set(index(i,0),jz->get(index(i,2)) - alpha*(jz->get(index(i,3)) - jz->get(index(i,1))));
    jj->set(index(i,0),jx->get(index(i,0))*jx->get(index(i,0)) +\
      jy->get(index(i,0))*jy->get(index(i,0)) +\
      jz->get(index(i,0))*jz->get(index(i,0)));

    jx->set(index(i,ny),jx->get(index(i,ny-2)) + alpha*(jx->get(index(i,ny-1)) - jx->get(index(i,ny-3))));
    jy->set(index(i,ny),jy->get(index(i,ny-2)) + alpha*(jy->get(index(i,ny-1)) - jy->get(index(i,ny-3))));
    jz->set(index(i,ny),jz->get(index(i,ny-2)) + alpha*(jz->get(index(i,ny-1)) - jz->get(index(i,ny-3))));
    jj->set(index(i,ny),jx->get(index(i,ny))*jx->get(index(i,ny)) +\
      jy->get(index(i,ny))*jy->get(index(i,ny)) +\
      jz->get(index(i,ny))*jz->get(index(i,ny)));
  }
  for (int j = 1; j <= ny-1; j++) {
    jx->set(index(0,j),jx->get(index(nx-1,j)));
    jy->set(index(0,j),jy->get(index(nx-1,j)));
    jz->set(index(0,j),jz->get(index(nx-1,j)));
    jj->set(index(0,j),jj->get(index(nx-1,j)));

    jx->set(index(nx,j),jx->get(index(1,j)));
    jy->set(index(nx,j),jy->get(index(1,j)));
    jz->set(index(nx,j),jz->get(index(1,j)));
    jj->set(index(nx,j),jj->get(index(1,j)));
  }
}

void MHD_2D::smooth(int k) {
  rhov->smooth(k);
  hv->smooth(k);
  uxv->smooth(k);
  uyv->smooth(k);
  uzv->smooth(k);
  //bxv->smooth(k);
  //byv->smooth(k);
  //bzv->smooth(k);
}

void MHD_2D::quit() {
  for(int j = 0; j <= ny; ++j) {
    for(int i = 0; i <= nx; ++i) {
      here = index(i,j);
      bxfile << lp_bx->get(here) << '\t';
      byfile << lp_by->get(here) << '\t';
      bzfile << lp_bz->get(here) << '\t';
      pfile << 2*pow(lp_h->get(here),gamma) << '\t';
      uzfile << lp_uz->get(here) << '\t';
      uxfile << lp_ux->get(here) << '\t';
      uyfile << lp_uy->get(here) << '\t';
      rhofile << lp_h->get(here) << '\t';
    }
  } 
  uzfile << std::endl;
  uxfile << std::endl;
  uyfile << std::endl;
  bxfile << std::endl;
  byfile << std::endl;
  bzfile << std::endl;
  pfile << std::endl;
  rhofile << std::endl;

  save_iteration();
  steps = 0;
}

void MHD_2D::save_iteration() {
  for(int j = 0; j <= ny; ++j) {
    for(int i = 0; i <= nx; ++i) {
      here = index(i,j);
      bxfile << bx->get(here) << '\t';
      byfile << by->get(here) << '\t';
      bzfile << bz->get(here) << '\t';
      pfile << p->get(here) << '\t';
      uzfile << uz->get(here) << '\t';
      uxfile << ux->get(here) << '\t';
      uyfile << uy->get(here) << '\t';
      rhofile << rho->get(here) << '\t';
    }
  } 
  uzfile << std::endl;
  uxfile << std::endl;
  uyfile << std::endl;
  bxfile << std::endl;
  byfile << std::endl;
  bzfile << std::endl;
  pfile << std::endl;
  rhofile << std::endl;
}


