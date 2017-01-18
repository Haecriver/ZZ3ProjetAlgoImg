#include <iostream>

#include "CImg.h"
#include "contAct.hpp"
#include "flotOpti.hpp"

int main(int argc, char* argv[]){

/*
	VALEUR carre1.png carre2.png
	./mosuofe -img1 ../rsc/carre1.bmp -img2 ../rsc/carre2.bmp -a 0.3 -t 0.05
	
	VALEUR carreBois1.png carreBois2.png
	./mosuofe -img1 ../rsc/carreBois1.bmp -img2 ../rsc/carreBois2.bmp -a 0.5 -t 0.05 -g 0.3 -b 0.3

	VALEUR taxi1.png taxi2.png NE MARCHE PAS TROP TROP BIEN QUOI
	./mosuofe -a 0.1 -t 0.2 -g 1.0 -b 0.5

*/

	const char* image_path_1 = cimg_option("-img1", "../rsc/taxi1.bmp", "Image 1");
	const char* image_path_2 = cimg_option("-img2", "../rsc/taxi2.bmp", "Image 1");

	int nbIter = cimg_option("-n", 10000, "Nombre d'itérations");
	float alpha = cimg_option("-a", 0.1, "Coefficient Alpha (propagation)");
	float beta = cimg_option("-b", 0.5, "Coefficient Beta (advection)");
	float gamma = cimg_option("-g", 0.5, "Coefficient Gamma (vitesse)");
	float theta = cimg_option("-t", 0.5, "Coefficient theta (geodesique)");
	float delta_t = cimg_option("-dt", 2.0, "Pas temporel");
	float ballon = cimg_option("-ballon", -0.01, "Force ballon");


	const char *firstImage = cimg_option("-i1", image_path_1, "Premiere image");
	const char *secondImage = cimg_option("-i2", image_path_2, "Seconde image");

	CImg<> seq = CImg<>(firstImage).channel(0);
	seq.append(CImg<>(secondImage).channel(0),'z');

	executeContourActif(seq, nbIter, alpha, beta, gamma, delta_t, ballon, theta);

	return EXIT_SUCCESS;
}
