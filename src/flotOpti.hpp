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
const CImg<double> moyenneur(3,3,1,1,1./9);
// Nb d'ite de HornSchunck
const unsigned k = 50;
// Valeur du parametre lambda
const double lambda = 10.0;

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

// Valeur seuille permettant le calcul de res en fonction des valeurs propres
const float tresholdVp = 10.0;

CImg<> LucasKanade(CImg<> seq, bool valeursPropres)
{
	// Taille des 
	int tailleV = 3;
	int tailleVCarre = tailleV * tailleV;
	int tailleVhalf = tailleV/2;

	// Res
	CImg<double> field(seq.width(),seq.height(),1,2);
	
	CImg<double> A(2,tailleVCarre,1,1);
	CImg<double> b(1,tailleVCarre,1,1);
	CImg<double> W(tailleVCarre,tailleVCarre,1,1);
	CImg<double> v;
	CImg<double> B;
	CImg<double> B_inv(2,2,1,1); // Inverse de B
	double detB;	// Determinant de B
	
	double proj;	// Coefficient pour trouver les valeurs H 

	CImgList<> vp;
	
  // Init de >
	W.identity_matrix();
	/*	
	// Ou on peut utiliser ce W (pour tailleV = 3 ):
	W(0,0)=1;
	W(1,1)=2;
	W(2,2)=3;
	W(3,3)=4;
	W(4,4)=5;
	W(5,5)=4;
	W(6,6)=3;
	W(7,7)=2;
	W(8,8)=1;
	*/

	// Calul du gradient de l'image
	CImgList<> grads =  seq.get_gradient("xyz",3);
	
	// On va parcourir toutes les valeurs d'image seq
	// Comme on remplit nos matrices avec ces valeurs, on a besoin de
	// faire attention au bordure
	// On parcourt donc de "1 à width -1"
	for(int i=tailleVhalf; i<seq.width()-tailleVhalf; i++){
		for(int j=tailleVhalf; j<seq.height()-tailleVhalf; j++){

			// Ici on calcul les matrices autour du point courant
			for(int x=-tailleVhalf; x<=tailleVhalf; x++){
				for(int y=-tailleVhalf; y<=tailleVhalf; y++){
					A(0,x+tailleVhalf + (y+tailleVhalf) * tailleV) = grads[0](x+i, y+j);
					A(1,x+tailleVhalf + (y+tailleVhalf) * tailleV) = grads[1](x+i, y+j);

					b(0,x+tailleVhalf + (y+tailleVhalf) * tailleV) = -grads[2](x+i, y+j);
				}
			}

			// Calcul de B
			B = A.get_transpose()*W*W*A;
			
			// Calcul de l'inverse de B
			//B_inv = B.get_invert();
			
			detB = 1 / (B(0,0)*B(1,1) - B(0,1)*B(1,0));
			B_inv(0,0) = detB * B(1,1);
			B_inv(0,1) = -detB * B(0,1);
			B_inv(1,0) = -detB * B(1,0);
			B_inv(1,1) = detB * B(0,0);
	
			// Calcul de v
			v = B_inv*A.get_transpose()*W*W*b;
			
			// Dans le cas ou on utilise valeurs propres
			if(valeursPropres){
				/*vp = B.get_eigen();
				double lambda1 = vp(0,0);
				double lambda2 = vp(0,1);*/
				
				double trB = B(0,0) + B(1,1);
				double lambda1 = trB/2 + sqrt((trB*trB) / (4-detB));
				double lambda2 = trB/2 - sqrt((trB*trB) / (4-detB));
						
				// Si lambda1, lambda2 grande
				if(lambda1>tresholdVp && lambda2>tresholdVp){
					field(i,j,0)=v[0];
					field(i,j,1)=v[1];
				// Si lambda1, lambda2 proche de 0
				}else if (lambda1<=tresholdVp && lambda2<=tresholdVp){
					field(i,j,0)=0.01;
					field(i,j,1)=0.01;
				// Si l'un des lambda est significativement non nulle
				}else{
					// On fait une projection de v sur le gradient en i,j
					proj = 
						((v[0])*grads[0](i,j) + (v[1])*grads[1](i,j)) 
						/ sqrt(grads[0](i,j)*grads[0](i,j)+grads[1](i,j)*grads[1](i,j));
					
					field(i,j,0)= ( (proj* grads[0](i,j)) / sqrt(grads[0](i,j)*grads[0](i,j)+grads[1](i,j)*grads[1](i,j)) );
					field(i,j,1)= ( (proj* grads[1](i,j)) / sqrt(grads[0](i,j)*grads[0](i,j)+grads[1](i,j)*grads[1](i,j)) );

				}
				
			// Sinon on utilise dans tous les cas les valeurs de v
			}else{
				field(i,j,0)=v[0];
				field(i,j,1)=v[1];
			}
		}
	}
	return field;
}

/*******************************************************************************

                                    Main

*******************************************************************************/
void executeFlotOptique()
{
 CImg<> seq = CImg<>("../rsc/taxi1.bmp").channel(0);
 seq.append(CImg<>("../rsc/taxi2.bmp").channel(0),'z');

 CImg<> displacementHS  = HornSchunck(seq);
 CImg<> displacementLK  = LucasKanade(seq,false);
 CImg<> displacementLK2 = LucasKanade(seq,true);

 // Affichage du champ r�sultat
 float color=500; unsigned int  sampling = 8; float  factor = 40; int  quiver_type = 0; float  opacity = 0.5;

 CImg<> imageHS = seq.get_slice(0).draw_quiver(displacementHS,&color,opacity,sampling,factor,quiver_type);
 CImgDisplay resHS_disp(imageHS,"Horn et Schunck");

 CImg<> imageLK = seq.get_slice(0).draw_quiver(displacementLK,&color,opacity,sampling,factor,quiver_type);
 CImgDisplay resLK_disp(imageLK,"Lucas et Kanade");

 CImg<> imageLK2 = seq.get_slice(0).draw_quiver(displacementLK2,&color,opacity,sampling,factor,quiver_type);
 CImgDisplay resLK2_disp(imageLK2,"Lucas et Kanade avec gestion des valeurs propres");

 while (!resHS_disp.is_closed() && !resLK_disp.is_closed() && !resLK2_disp.is_closed())
 {

 }
}

#endif
