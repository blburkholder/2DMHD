/* Brandon Burkholder
   3-14-2017
   Phys 629 Numerical Simulation in Fluids and Plasmas 
   Plasmoid Project */

#include "MHD_2D.h"
#include <string>
int main(int argc, char *argv[]) {
    std::ifstream infile("input.txt");

    float deltax, deltay, deltat, gam, alph, sig, resistivity, b_guide;
    int width, height, number_of_steps, output, init, close;

    std::string b;

    infile >> init >> b;
    infile >> deltax >> b;
    infile >> deltay >> b;
    infile >> width >> b;
    infile >> height >> b;
    infile >> deltat >> b;
    infile >> number_of_steps >> b;
    infile >> output >> b;
    infile >> gam >> b;
    infile >> alph >> b;
    infile >> sig >> b;
    infile >> b_guide >> b;
    infile >> resistivity >> b;
    infile >> close >> b;


    MHD_2D plasmoid(init, deltax, deltay, width, height,\
      deltat, number_of_steps, output, gam, alph, b_guide, sig,\
      resistivity, close);

    plasmoid.leap();
}
