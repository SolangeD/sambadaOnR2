/*************************************************************************
 * Copyright (©) 2011-2015 EPFL (Ecole Polytechnique fédérale de Lausanne)
 * Laboratory of Geographic information systems (LaSIG)
 * 
 * This file is part of Sambada.
 *  
 * Sambada is free software ; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation ; either version 3 of the License, or (at your option) any later version.
 * Sambada is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY ; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with Sambada ; if not, see <http://www.gnu.org/licenses/>.
 * 
 * Authors : Sylvie Stucki (sylvie.stucki@a3.epfl.ch), Stéphane Joost (stephane.joost@epfl.ch) 
 * Laboratory of Geographic information systems
 * EPFL ENAC IIE LASIG
 * Station 18
 * CH-1015 Lausanne
 * Web site : http://lasig.epfl.ch/sambada
 * 
 * Sambada includes two libraries: Scythe Statistical Library (under GPL 3) and Shapefile C Library (under LGPL 2.1, courtesy of Frank Warmerdam).
 * 
 * Scythe Statistical Library
 * Copyright (C) 2000-2002 Andrew D. Martin and Kevin M. Quinn;
 * 2002-2012 Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  All Rights Reserved.
 * 
 * Shapefile C Library
 * Copyright (c) 1999, Frank Warmerdam
 *************************************************************************/


#include "RegressionLogistique.h"
#include "Supervision2.h"
#include <ctime>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace scythe;

// [[Rcpp::export]]
int supervision(std::string nomFichier, int numBlock, int blockSize, int maxDimension, std::string selScore, double scoreThreshold, std::string sortScore, char wordDelim)
{	
	/*if (argc==1) 
	{
		RegressionLogistique::messageBienvenue(cout, true);
	}
	else
	{
		RegressionLogistique::messageBienvenue(cout);
	}
	*/
	
	Supervision regisseur;
	
	/*if(argc!=2 && argc!=5 && argc!=7 && argc !=8 && argc!=9)
	{
		cerr << "Nombre d'arguments incorrects." << endl;
		return 1;
	}*/
	
	time_t temps_start(time(NULL));
	try
	{
		//cerr << string(argv[1]) << endl;
		if (numBlock==0)
		{
			regisseur.preparationsCalculs(nomFichier);
		}
		else 
		{	
			typeScore scoreSel=Both;
			if(selScore=="G"){
				scoreSel=G;
			}
			else if(selScore=="Wald"){
				scoreSel=Wald;
			}


			typeScore scoreTri=Wald;
			if(sortScore=="G"){
				typeScore scoreTri=G;
			}
			else if(sortScore=="AIC"){
				scoreTri=AIC;
			}
			else if(sortScore=="BIC"){
				scoreTri=BIC;
			}
			
			reel seuilScore(scoreThreshold);
			if(scoreThreshold==0){
				reel seuilScore(0);
			}
			
			regisseur.fusionResultats(nomFichier, numBlock, blockSize, maxDimension, scoreSel, seuilScore, scoreTri, wordDelim);
		}

	}
	catch (const Erreur& err) 
	{
		cout << err.what() << endl;
		return(1); //Au lieu de exit(1) pour que R ne soit pas vexe
	}
	//time_t temps_interm(time(NULL));
	//cout << "Fin de la lecture : " << difftime(temps_interm, temps_start) << " s." << endl;
	
	
	time_t temps_fin_calculs(time(NULL));
	//cout << "Fin des calculs " << endl;
	cout << "Temps écoulé : " << difftime(temps_fin_calculs, temps_start) << " s." << endl;
	
	
	time_t temps_stop(time(NULL));
	
	cout << "Ecriture des résultats : " << difftime(temps_stop, temps_fin_calculs) << " s." << endl;
	
}


