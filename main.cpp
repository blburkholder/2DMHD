/* Brandon Burkholder
   3-14-2017
   Phys 629 Numerical Simulation in Fluids and Plasmas 
   Plasmoid Project */

#include "MHD_2D.h"
#include "reconnection.h"
#include "KH.h"
#include "alfven.h"
#include "bomb.h"
#include <string>
int main(int argc, char *argv[]) {
    std::ifstream infile("input.txt");

    float deltax, deltay, deltat, gam, alph, sig, resistivity, b_guide, bs;
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
    infile >> bs >> b;
    infile >> b_guide >> b;
    infile >> resistivity >> b;
    infile >> close >> b;

    switch(init) {
	case 1 : {
		alfven ynwa(deltax, deltay, width, height,\
			deltat, number_of_steps, output, gam, alph, b_guide, sig,\
			bs, resistivity, close);
                 ynwa.leap();
		break; }

	case 2 : {
		bomb ynwa(deltax, deltay, width, height,\
			deltat, number_of_steps, output, gam, alph, b_guide, sig,\
			bs, resistivity, close);
                ynwa.leap();
		break; }

	case 3 : {
    		reconnection ynwa(deltax, deltay, width, height,\
      			deltat, number_of_steps, output, gam, alph, b_guide, sig,\
      			bs, resistivity, close);
    		ynwa.leap();
		break; }
	case 4 : {
    		KH ynwa(deltax, deltay, width, height,\
      			deltat, number_of_steps, output, gam, alph, b_guide, sig,\
      			bs, resistivity, close);
                ynwa.leap();
		break; }
        default:
	std::cout << "wrong" << std::endl;
    }
}
