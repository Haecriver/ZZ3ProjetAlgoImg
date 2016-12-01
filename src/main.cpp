#include "CImg.h"
#include "contAct.hpp"
#include "flotOpti.hpp"

int main(int argc, char* argv[]){

	int nbIter = cimg_option("-n", 2000, "Nombre d'itérations");
	float alpha = cimg_option("-a", 0.05, "Coefficient Alpha");
	float beta = cimg_option("-b", 0.5, "Coefficient Beta");
	float gamma = cimg_option("-g", 0.8, "Coefficient Gamma");
	
 	CImg<> seq = CImg<>("../rsc/taxi1.bmp").channel(0);
 	seq.append(CImg<>("../rsc/taxi2.bmp").channel(0),'z');
	executeContourActif(seq, nbIter, alpha, beta, gamma);
	return EXIT_SUCCESS;
}
