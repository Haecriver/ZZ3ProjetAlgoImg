#include <iostream>

#include "CImg.h"
#include "contAct.hpp"
#include "flotOpti.hpp"

int main(int argc, char* argv[]){

	int nbIter = cimg_option("-n", 10000, "Nombre d'itérations");
	float alpha = cimg_option("-a", 0.05, "Coefficient Alpha (propagation)");
	float beta = cimg_option("-b", 0.5, "Coefficient Beta (advection)");
	float gamma = cimg_option("-g", 0.5, "Coefficient Gamma (vitesse)");
	float delta_t = cimg_option("-dt", 2.0, "Pas temporel");
	float ballon = cimg_option("-ballon", -0.01, "Force ballon");
	float theta = cimg_option("-t", 1.0, "Force theta");

	std::cout << "teta ;" << theta << std::endl;

	const char *firstImage = cimg_option("-i1", "../rsc/taxi1.bmp", "Premiere image");
	const char *secondImage = cimg_option("-i2", "../rsc/taxi2.bmp", "Seconde image");

	CImg<> seq = CImg<>(firstImage).channel(0);
	seq.append(CImg<>(secondImage).channel(0),'z');

	executeContourActif(seq, nbIter, alpha, beta, gamma, delta_t, ballon, theta);

	return EXIT_SUCCESS;
}
