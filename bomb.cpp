#include "bomb.h"

bomb::bomb(float deltax, float deltay, int width, int height,\
    float deltat, int number_of_steps, int output, float gam, float alph, float bg,\
    float bs, float sig, float resistivity, int cl) :
  MHD_2D(deltax,deltay,width,height,deltat,number_of_steps,output,gam,alph,bg,bs,sig,resistivity,cl)
{uxfile.open("b_ux.txt");
    uyfile.open("b_uy.txt");
    uzfile.open("b_uz.txt");
    bxfile.open("b_bx.txt");
    byfile.open("b_by.txt");
    bzfile.open("b_bz.txt");
    rhofile.open("b_rho.txt");
    pfile.open("b_p.txt");}

//plasma bomb
//secondary closure relation not implemented
void bomb::initialize_grid() {
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
      bx->set(here,1);      
      by->set(here,0);
      bz->set(here,b_guide);
      eta->set(here,res);

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

void bomb::j_bound() {
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

void bomb::bound() {
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
