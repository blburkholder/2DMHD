#include "reconnection.h"

reconnection::reconnection(float deltax, float deltay, int width, int height,\
    float deltat, int number_of_steps, int output, float gam, float alph, float bg,\
    float bs, float sig, float resistivity, int cl) :
  MHD_2D(deltax,deltay,width,height,deltat,number_of_steps,output,gam,alph,bg,bs,sig,resistivity,cl) 
{uxfile.open("ux.txt");
    uyfile.open("uy.txt");
    uzfile.open("uz.txt");
    bxfile.open("bx.txt");
    byfile.open("by.txt");
    bzfile.open("bz.txt");
    rhofile.open("rho.txt");
    pfile.open("p.txt");}

void reconnection::initialize_grid() {
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
      y = dy*(float)j - (float)ddy/2.0;
      //pressure in the middle
      px = c1*pow(1.0+(float)i*dx/xd,-k1) + c2*tanh(k2*((float)i*dx/xd - 1.0)) - c2 + p_inf;

      p->set(here,0.25 + (1.0/(cosh(y*sqrt(px))*cosh(y*sqrt(px))))*px);
      rho->set(here,p->get(here));

      bx->set(here,-sqrt(px)*tanh(y*sqrt(px)));
      by->set(here,(-0.5/px)*(1.0 - y*sqrt(px)*tanh(y*sqrt(px)))*\
        (c2*k2*(1.0/(cosh(k2*((float)i*dx/xd-1.0))*cosh(k2*((float)i*dx/xd-1.0))))/xd -\
        k1*c1*pow(1.0+(float)i*dx/xd,-k1-1.0)/xd));
      bz->set(here,b_guide);

      ux->set(here,0.0);
      uy->set(here,0.0);
      uz->set(here,0.0);
      //leap quantities are the same for first iteration

      lp_rho->set(here,rho->get(here));
      lp_ux->set(here,0.0);
      lp_uy->set(here,0.0);
      lp_uz->set(here,0.0);
      lp_bx->set(here,bx->get(here));
      lp_by->set(here,by->get(here));
      lp_bz->set(here,bz->get(here)); 

      if ( close == 1 ) lp_h->set(here,pow(p->get(here)/2.0,1.0/gamma));
      else lp_h->set(here,0.5*rho->get(here)*(ux->get(here)*ux->get(here) +\
        uy->get(here)*uy->get(here) + uz->get(here)*uz->get(here)) + \
        0.5*(bx->get(here)*bx->get(here) + by->get(here)*by->get(here) +\
        bz->get(here)*bz->get(here)) + 0.5*p->get(here)/(gamma-1.0));

      h->set(here,lp_h->get(here));
      eta->set(here,res);

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
  for (int j = 1; j <= ny-1; j++) {
    for (int i = 1; i <= nx-1; i++) {
      here = index(i,j);
      jx->set(here,(0.5/dy)*(bz->get(index(i,j+1)) - bz->get(index(i,j-1))));
      jy->set(here,(-0.5/dx)*(bz->get(index(i+1,j)) - bz->get(index(i-1,j))));
      jz->set(here,(-0.5/dy)*(bx->get(index(i,j+1)) - bx->get(index(i,j-1))) +\
        (0.5/dx)*(by->get((index(i+1,j))) - by->get((index(i-1,j)))));
      jj->set(here,jx->get(here)*jx->get(here) + jy->get(here)*jy->get(here) + jz->get(here)*jz->get(here));
    }
  }
  j_bound();
  save_iteration();
}
//combine mag bound into here
void reconnection::bound() {
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
    //add another one to this average? weight it differently? (more?)
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
    mag_bound();
}

void reconnection::mag_bound() {
  *bx = *bxv;
  *by = *byv;

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

    bxv->set(0,0.5*(bxv->get(index(0,1)) + bxv->get(index(1,0))));
    bxv->set((nx+1)*ny,0.5*(bxv->get(index(0,ny-1)) + bxv->get(index(1,ny))));
    bxv->set(nx,0.5*(bxv->get(index(nx-1,0)) + bxv->get(index(nx,1))));
    bxv->set((nx+1)*ny + nx,0.5*(bxv->get(index(nx,ny-1)) + bxv->get(index(nx-1,ny))));

    byv->set(0,0.5*(byv->get(index(0,1)) + byv->get(index(1,0))));
    byv->set((nx+1)*ny,0.5*(byv->get(index(0,ny-1)) + byv->get(index(1,ny))));
    byv->set(nx,0.5*(byv->get(index(nx-1,0)) + byv->get(index(nx,1))));
    byv->set((nx+1)*ny + nx,0.5*(byv->get(index(nx,ny-1)) + byv->get(index(nx-1,ny))));

    bzv->set(index(0,0),0.5*(bzv->get(index(0,1)) + bzv->get(index(1,0))));
    bzv->set((nx+1)*ny,0.5*(bzv->get(index(0,ny-1)) + bzv->get(index(1,ny))));
    bzv->set(nx,0.5*(bzv->get(index(nx-1,0)) + bzv->get(index(nx,1))));
    bzv->set((nx+1)*ny + nx,0.5*(bzv->get(index(nx,ny-1)) + bzv->get(index(nx-1,ny))));
}

void reconnection::j_bound() {
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
