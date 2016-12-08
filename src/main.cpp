#include "CImg.h"
#include "contAct.hpp"
#include "flotOpti.hpp"

int main(int argc, char* argv[]){

	int nbIter = cimg_option("-n", 2000, "Nombre d'itérations");
	float alpha = cimg_option("-a", 0.05, "Coefficient Alpha");
	float beta = cimg_option("-b", 0.5, "Coefficient Beta");
	float gamma = cimg_option("-g", 0.8, "Coefficient Gamma");

	const char *firstImage = cimg_option("-i1", "../rsc/taxi1.bmp", "Premiere image");
	const char *secondImage = cimg_option("-i2", "../rsc/taxi2.bmp", "Seconde image");

	CImg<> seq = CImg<>(firstImage).channel(0);
	seq.append(CImg<>(secondImage).channel(0),'z');

	executeContourActif(seq, nbIter, alpha, beta, gamma);

	return EXIT_SUCCESS;
}
