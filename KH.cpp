#include "KH.h"

KH::KH(float deltax, float deltay, int width, int height,\
    float deltat, int number_of_steps, int output, float gam, float alph, float bg,\
    float bs, float sig, float resistivity, int cl) :
  MHD_2D(deltax,deltay,width,height,deltat,number_of_steps,output,gam,alph,bg,bs,sig,resistivity,cl) 
{uxfile.open("KH_ux.txt");
    uyfile.open("KH_uy.txt");
    uzfile.open("KH_uz.txt");
    bxfile.open("KH_bx.txt");
    byfile.open("KH_by.txt");
    bzfile.open("KH_bz.txt");
    rhofile.open("KH_rho.txt");
    pfile.open("KH_p.txt");}

void KH::initialize_grid() {
  float y = 0;
  for (int j = 0; j <= ny; j++) {
    for (int i = 0; i <= nx; i++) {
      here = index(i,j);
      y = dy*j - ddy/2;     
      
      rho->set(here,1.0);
      ux->set(here,0.5*tanh(y));
      uy->set(here,-0.1*exp(-(16.0*y*y))*sin((float)i*dx/2));
      uz->set(here,0);

      bx->set(here,1.0);      
      by->set(here,0);
      //bz->set(here,tanh(i*dx));
      bz->set(here,0);
      eta->set(here,res);
      p->set(here,1.0);

      if ( close == 1 ) h->set(here,pow(p->get(here)/2.0,1.0/gamma));
      else h->set(here,0.5*rho->get(here)*(ux->get(here)*ux->get(here) +\
        uy->get(here)*uy->get(here) + uz->get(here)*uz->get(here)) + \
        0.5*(bx->get(here)*bx->get(here) + by->get(here)*by->get(here) +\
        bz->get(here)*bz->get(here)) + 0.5*p->get(here)/(gamma-1.0));

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
  for (int j = 1; j <= ny-1; j++) {
    for (int i = 1; i <= nx-1; i++) {
      jx->set(here,(0.5/dy)*(bz->get(index(i,j+1)) - bz->get(index(i,j-1))));
      jy->set(here,(-0.5/dx)*(bz->get(index(i+1,j)) - bz->get(index(i-1,j))));
      jz->set(here,(-0.5/dy)*(bx->get(index(i,j+1)) - bx->get(index(i,j-1))) +\
        (0.5/dx)*(by->get((index(i+1,j))) - by->get((index(i-1,j)))));
      jj->set(here,jx->get(here)*jx->get(here) +\
        jy->get(here)*jy->get(here) + jz->get(here)*jz->get(here));
    }
  }
  bound();
}

void KH::bound() {
  *bx = *bxv;
  *by = *byv;
  for (int i = 1; i <= nx-1; i++) {
    //bottom
    here = index(i,0);
    up = index(i,2);
    rhov->set(left,rhov->get(up) - alpha*(rhov->get(index(i,3)) - rhov->get(index(i,1))));
    hv->set(left,hv->get(up) - alpha*(hv->get(index(i,3)) - hv->get(index(i,1))));
    uxv->set(left,uxv->get(up));
    uyv->set(left,-uyv->get(up));
    uzv->set(left,-uzv->get(up));
    bxv->set(left,bx->get(up));
    byv->set(left,by->get(up) + (dy/dx)*(bx->get(index(i+1,1)) - bx->get(index(i-1,1))));
    bzv->set(left,bzv->get(up) - alpha*(bzv->get(index(i,3)) - bzv->get(index(i,1)))); 
    //top
    here = index(i,ny);
    down = index(i,ny-2);
    rhov->set(index(i,ny),rhov->get(down) + alpha*(rhov->get(index(i,ny-1)) - rhov->get(index(i,ny-3))));
    hv->set(index(i,ny),hv->get(down) + alpha*(hv->get(index(i,ny-1)) - hv->get(index(i,ny-3))));
    uxv->set(index(i,ny),uxv->get(down));
    uyv->set(index(i,ny),-uyv->get(down));
    uzv->set(index(i,ny),-uzv->get(down));
    bxv->set(index(i,ny),bx->get(down));
    byv->set(index(i,ny),by->get(down) - (dy/dx)*(bx->get(index(i+1,ny-1)) - bx->get(index(i-1,ny-1))));
    bzv->set(index(i,ny),bzv->get(down) + alpha*(bzv->get(index(i,ny-1)) - bzv->get(index(i,ny-3))));;
  }

  for (int j = 1; j <= ny-1; j++) {
    //left
    here = index(0,j);
    right = index(nx-1,j);
    rhov->set(here,rhov->get(right));
    hv->set(here,hv->get(right));
    uxv->set(here,uxv->get(right));
    uyv->set(here,uyv->get(right));
    uzv->set(here,uzv->get(right));
    bxv->set(here,bxv->get(right));
    byv->set(here,byv->get(right));
    bzv->set(here,bzv->get(right));
    //right
    here = index(nx,j);
    left = index(1,j);
    rhov->set(here,rhov->get(left));
    hv->set(here,hv->get(left));
    uxv->set(here,uxv->get(left));
    uyv->set(here,uyv->get(left));
    uzv->set(here,uzv->get(left));
    bxv->set(here,bxv->get(left));
    byv->set(here,byv->get(left));
    bzv->set(here,bzv->get(left));
  }

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

void KH::j_bound() {
  for (int i = 1; i <= nx-1; i++) {
    //bottom
    jx->set(index(i,0),jx->get(index(i,2)) - alpha*(jx->get(index(i,3)) - jx->get(index(i,1))));
    jy->set(index(i,0),jy->get(index(i,2)) - alpha*(jy->get(index(i,3)) - jy->get(index(i,1))));
    jz->set(index(i,0),jz->get(index(i,2)) - alpha*(jz->get(index(i,3)) - jz->get(index(i,1))));
    jj->set(index(i,0),jx->get(index(i,0))*jx->get(index(i,0)) +\
      jy->get(index(i,0))*jy->get(index(i,0)) +\
      jz->get(index(i,0))*jz->get(index(i,0)));
    //top
    jx->set(index(i,ny),jx->get(index(i,ny-2)) + alpha*(jx->get(index(i,ny-1)) - jx->get(index(i,ny-3))));
    jy->set(index(i,ny),jy->get(index(i,ny-2)) + alpha*(jy->get(index(i,ny-1)) - jy->get(index(i,ny-3))));
    jz->set(index(i,ny),jz->get(index(i,ny-2)) + alpha*(jz->get(index(i,ny-1)) - jz->get(index(i,ny-3))));
    jj->set(index(i,ny),jx->get(index(i,ny))*jx->get(index(i,ny)) +\
      jy->get(index(i,ny))*jy->get(index(i,ny)) +\
      jz->get(index(i,ny))*jz->get(index(i,ny)));
  }
  for (int j = 1; j <= ny-1; j++) {
    //left
    jx->set(index(0,j),jx->get(index(nx-1,j)));
    jy->set(index(0,j),jy->get(index(nx-1,j)));
    jz->set(index(0,j),jz->get(index(nx-1,j)));
    jj->set(index(0,j),jj->get(index(nx-1,j)));
    //right
    jx->set(index(nx,j),jx->get(index(1,j)));
    jy->set(index(nx,j),jy->get(index(1,j)));
    jz->set(index(nx,j),jz->get(index(1,j)));
    jj->set(index(nx,j),jj->get(index(1,j)));
  }
}



