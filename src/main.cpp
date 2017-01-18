#include <iostream>

#include "CImg.h"
#include "contAct.hpp"
#include "flotOpti.hpp"

int main(int argc, char* argv[]){

/*
	VALEUR carre1.png carre2.png
	./mosuofe -a 0.3 -t 0.05
	
	VALEUR carreBois1.png carreBois2.png
	./mosuofe -a 0.5 -t 0.05 -g 0.3 -b 0.3

*/

	int nbIter = cimg_option("-n", 10000, "Nombre d'itérations");
	float alpha = cimg_option("-a", 0.1, "Coefficient Alpha (propagation)");
	float beta = cimg_option("-b", 0.5, "Coefficient Beta (advection)");
	float gamma = cimg_option("-g", 0.5, "Coefficient Gamma (vitesse)");
	float delta_t = cimg_option("-dt", 2.0, "Pas temporel");
	float ballon = cimg_option("-ballon", -0.01, "Force ballon");
	float theta = cimg_option("-t", 0.5, "Force theta");

	std::cout << "teta ;" << gamma << std::endl;

	const char *firstImage = cimg_option("-i1", "../rsc/carreBois1.bmp", "Premiere image");
	const char *secondImage = cimg_option("-i2", "../rsc/carreBois2.bmp", "Seconde image");

	CImg<> seq = CImg<>(firstImage).channel(0);
	seq.append(CImg<>(secondImage).channel(0),'z');

	executeContourActif(seq, nbIter, alpha, beta, gamma, delta_t, ballon, theta);

	return EXIT_SUCCESS;
}
