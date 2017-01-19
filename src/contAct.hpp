////////////////////////////////////////////////////////////////////////////////
//                       TP 7 : Contours actifs implicites                    //
////////////////////////////////////////////////////////////////////////////////

#ifndef CONT_ACT_HPP
#define CONT_ACT_HPP

#include "CImg.h"
#include "flotOpti.hpp"
#include <iostream>

#define MAX(a,b) (((a)<(b)) ? (b) : (a) )
#define MIN(a,b) (((a)<(b)) ? (a) : (b) )


using namespace cimg_library;

/*******************************************************************************

     Calcul de la position approximative du contour à partir de la carte
   levelset (psi): isocontour de valeur 0. On recherche donc l'ensemble des
    pixels possédant un voisin de signe opposé (interface entre domaine <0
                                et domaine >0)

*******************************************************************************/
CImg<> ExtractContour(CImg<> LevelSet)
{
	CImg<> Contour(LevelSet.width(), LevelSet.height(), 1, 1);
	Contour.fill(0);

	CImg_3x3(I, double);
	cimg_for3x3(LevelSet, x, y, 0, 0, I, double)
	{
		// Comparaison avec les points voisin => si signe opposé alors on est sur un contour
		if(Icc*Icp <= 0 || Icc*Icn <= 0 || Icc*Ipc <= 0 || Icc*Inc <= 0)
			Contour(x, y) = 1;
	}
	
	return Contour;
}

/*******************************************************************************

                     Affichage d'un contour dans une image

*******************************************************************************/
void DrawContour(CImg<>* imgIn,CImg<> Contour)
{
	const double color = 2*imgIn->max();

	cimg_forXY(Contour,x,y)
	{
		if(Contour(x,y) == 1)
			imgIn->draw_point(x, y, 0, &color, 0.7);
	}
}

/*******************************************************************************

     Initialisation du LevelSet (psi) à l'aide de la distance euclidienne
    signée. Le contour initial est un cercle de centre (x0,y0) et de rayon r

*******************************************************************************/
void InitLevelSet(CImg<>* imgIn, int x0,int y0, int r)
{
	cimg_forXY(*imgIn,x,y)
	{
		(*imgIn)(x,y) = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))-r;
	}
}

/*******************************************************************************

      Algorithme de propagation d'un contour implicite (modéle géodésique)

*******************************************************************************/
void Propagate(
	CImg<> imgIn,
	CImg<>* LevelSet,
	int nbIter,						// Nombre d'itération
	float alpha, 					// Pondération du terme de propagation
	float beta, 						// Pondération du terme d'advection
	float gamma,
	float delta_t,					// Pas temporel
	float ballon,
	float theta) 				// Force ballon
{
	// Calcul des gradients spaciaux régularisés de l'image
	const CImgList<> gradImg = imgIn.get_channel(0).get_gradient("xy", 4 /* Filtre qui régularise */);
	const CImg<> squareGradImgNorme = (gradImg[0].get_sqr() + gradImg[1].get_sqr());
	
	// Calculer le module de chaque vecteur
	const CImg<> imgInHS = HornSchunck(imgIn);

/******************************************************************************/
// Affichage HornSchunck
{
	//float color = 500; unsigned int sampling = 8; float factor = 40; int quiver_type = 0; float opacity = 0.5;

	//CImg<> imageHS = imgIn.get_slice(0).draw_quiver(imgInHS, &color, opacity, sampling, factor, quiver_type);
    //CImgDisplay resHS_disp(imageHS, "Horn et Schunck");

	//while (!resHS_disp.is_closed()) { resHS_disp.wait(); }
}
/******************************************************************************/

	// Calcul du module des vecteurs vitesse
	const CImg<> imgInVitesse = ((imgInHS.get_channel(0).get_sqr() + imgInHS.get_channel(1).get_sqr()).get_sqrt());

	// Normalement pas util
	const CImgList<> gradimgInVitesse = imgInVitesse.get_gradient("xy",0);
	const CImg<> squareGradImgVitess = (gradimgInVitesse[0].get_sqr() + gradimgInVitesse[1].get_sqr());
	
	// Valeur quand ca commence a pete
	//gamma : 0.505139

	// Calcul de la fonction d'arrêt g du modèle géodésique
	CImg<> g(imgIn.width(), imgIn.height(), 1, 1);
	cimg_forXY(g, x, y)
	{
		g(x, y) = ( 1.0 / double(1.0 + squareGradImgNorme(x, y)*(float)theta ) ) + ballon - (float)imgInVitesse(x,y)*(float)gamma;
	}
	
	/////////////// AFFICHAGE DEBUG ///////////////
	
	// CImgDisplay resg_disp(g, "gedodesique");
	// CImgDisplay resVitesse_disp(imgInVitesse, "Vitesse");
	// CImgDisplay resGradVitesse_disp(squareGradImgVitess, "Gradient Vitesse");
	// CImgDisplay ressquareGradImgNorme_disp(squareGradImgNorme, "Norme grad");
	
	// while (!resVitesse_disp.is_closed() && 
	// !resGradVitesse_disp.is_closed() && 
	// !ressquareGradImgNorme_disp.is_closed() &&
	// !resg_disp.is_closed()) 
	// {
	// 	resVitesse_disp.wait(); 
	// }
	
	
	const CImgList<> gradG = g.get_gradient("xy");

	// Note : p = plus / m = moins
	for (int i = 0 ; i < nbIter ; ++i)
	{
		// Calcul des dérivées de LevelSet
		const CImgList<> Dm = LevelSet->get_gradient("xy", -1);
		//const CImgList<> Dc = LevelSet->get_gradient("xy", 0);
		const CImgList<> Dp = LevelSet->get_gradient("xy", 1);
				
		cimg_forXY(*LevelSet, x, y)
		{		
			// Valeurs des dérivées de LevelSet
			const float& Dxm = Dm[0](x, y);
			const float& Dxp = Dp[0](x, y);
			const float& Dym = Dm[1](x, y);
			const float& Dyp = Dp[1](x, y);
			
			// Calcul des nablas
			// Pour Dx
			const float maxDxm0 = MAX(Dxm, 0);
			const float minDxm0 = MIN(Dxm, 0);
			const float maxDxp0 = MAX(Dxp, 0);
			const float minDxp0 = MIN(Dxp, 0);
			
			// Pour Dy
			const float maxDym0 = MAX(Dym, 0);
			const float minDym0 = MIN(Dym, 0);
			const float maxDyp0 = MAX(Dyp, 0);
			const float minDyp0 = MIN(Dyp, 0);
			
			// Nablas
			const float nablap = std::sqrt(	maxDxm0 * maxDxm0
											+ minDxp0 * minDxp0
											+ maxDym0 * maxDym0
											+ minDyp0 * minDyp0);
			const float nablam = std::sqrt(	maxDxp0 * maxDxp0
											+ minDxm0 * minDxm0
											+ maxDyp0 * maxDyp0
											+ minDym0 * minDym0);
		
			// Calcul de la vitesse de propagation
			const float Fprop = - (MAX(g(x, y), 0) * nablap + MIN(g(x, y), 0) * nablam);

			// Calcul de la vitesse d'advection
			const float& u = gradG[0](x, y);
			const float& v = gradG[1](x, y);
			const float Fadv = - (MAX(u, 0) * Dxm + MIN(u, 0) * Dxp + MAX(v, 0) * Dym + MIN(v, 0) * Dyp);
			
			(*LevelSet)(x, y) += delta_t * (alpha * Fprop + beta * Fadv );
		}
	}
}

/*******************************************************************************

                                      main

*******************************************************************************/

void executeContourActif(	
	CImg<> img, 
	int nbIter,
	float alpha, 
	float beta, 
	float gamma,
	float delta_t, 
	float ballon,
	float theta)
{
	// Définition d'un contour initial circulaire
	int x0 = img.width()/2;
	int y0 = img.height()/2;
	int r  = (img.height() > img.width() ? img.width() / 2 - 10 : img.height() / 2 - 10);

	CImg<> levelset(img.width(), img.height(), 1, 1);
	InitLevelSet(&levelset, x0, y0, r);

	// Propagation du contour
	Propagate(img, &levelset, nbIter, alpha, beta, gamma, delta_t, ballon,theta);

	// Extraction du résultat
	CImg<> Contour = ExtractContour(levelset);

	// Tracé du résultat dans l'image
	DrawContour(&img, Contour);

	// Affichage finale
	CImgDisplay dispSpatial(img, "Image Segmentée");

	while (!dispSpatial.is_closed())
	{
		dispSpatial.wait();
	}
}

#endif
