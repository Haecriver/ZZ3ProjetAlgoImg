////////////////////////////////////////////////////////////////////////////////
//                          TP 6 : Flot optique                     //
////////////////////////////////////////////////////////////////////////////////

#ifndef FLOT_OPTI_HPP
#define FLOT_OPTI_HPP

#include "CImg.h"
#include <math.h>
#include <iostream>

using namespace cimg_library;

// Matrice de moyennage
static const CImg<double> moyenneur(3,3,1,1,1./9);
// Nb d'ite de HornSchunck
static const unsigned k = 50;
// Valeur du parametre lambda
static const double lambda = 10.0;

CImg<> HornSchunck(CImg<> seq)
{
	// On initialise le resultat
	// 2 channels et 1 de depth
	// initialise avec 0
	CImg<> field(seq.width(),seq.height(),1,2,0);

	// On calcul le gradient de l'image donnee
	CImgList<> grad = seq.get_gradient("xyz",0);

	// Pour un nombre k d'iteration
	for(unsigned i =0; i<k ; i++){
		// On calcul la moyenne des valeurs
		CImg<> moy = field.get_convolve(moyenneur);
		
		// On parcours le resultat pour le assigner les valeurs
		cimg_forXY(field,x,y){
			
			// On associe les valeurs de la formule
 			field(x,y,0) =
 				moy(x,y,0) - ( grad(0)(x,y) * (grad(0)(x,y) * moy(x,y,0)  + grad(1)(x,y) * moy(x,y,1) + grad(2)(x,y)) ) /
 				( (lambda * lambda) + (grad(0)(x,y) * grad(0)(x,y)) + (grad(1)(x,y) * grad(1)(x,y)) );

 			field(x,y,1) =
	 			moy(x,y,1) - (grad(1)(x,y) * (grad(0)(x,y) * moy(x,y,0) + grad(1)(x,y) * moy(x,y,1) + grad(2)(x,y))) /
 				( (lambda * lambda) + (grad(0)(x,y) * grad(0)(x,y)) + (grad(1)(x,y) * grad(1)(x,y)) );
		}

	}

	return field;
}

#endif
