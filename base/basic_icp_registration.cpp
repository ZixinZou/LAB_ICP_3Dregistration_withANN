#include "stdafx.h"
#include <stdio.h>
#include "Matrix.h"
// pour le tps d'exécution
#include <windows.h>

typedef float *TabFloat;

int main(int argc, char* argv[])
{
	// initialisation des fichiers
	FILE *in1 = fopen("foot.pts.neu","r"); // to rename files   bunny.pts.neu   PontchartrainVirageI.neu
	FILE *in2 = fopen("foot_perturbed.pts.neu","r"); //  bunny_perturbed.pts.neu  PontchartrainVirageII_prerecale.neu

	FILE *out1 = fopen("result_Foot.neu","w");
	FILE *out2 = fopen("error_Foot.txt","w");

	FILE *out3 = fopen("transformation_Foot.txt","w");
	FILE *out4= fopen("tempsdexecution_Foot.txt","w");
	// ******************************************************************
	/* Variables : compteurs*/
	int i,j, iteration = 0;
	float x,y,z;
	int n=3,NbrPt1=0,NbrPt2=0;//nombre de point

	/* Variables de temps d'exécution*/
	LARGE_INTEGER start, end, freq;
    double elapsed;

	/*tableaux de points*/
	TabFloat *Point1,*Point2, *Voisin;
	//-- Nos variables --
	//transformation T(R,t)
	Matrix R(n,n,false); // Matrice de rotation
	Vector t(3); //vecteur de translation

	//Calcul de la SVD {H.ComputeSVD(U,W,V);}
	Matrix H(n,n,false); //Matrice
	Matrix U(3,3,false);
	Matrix W(3,3,false);
	Matrix V(3,3,false);

	int NbrPtC, NbrPtV; // nombre des points des fichiers
	float seuil; // le seuil ?utilis?par ICP
	float erreur = 1000; // l'écart type moyen (~1000) valeur par défaut
	float max; //utilis?pour la distance maximale entre 2 points (distance cartésienne)
	float min; //utilisée pour la distance minimale entre 2 points (distance cartésienne)
	int k; // compteur des itération d'ICP
	int l; // indices utilisable pour la programmation
	float* PointG1 = new float[3]; //centre de gravit?nuage 1
	float* PointG2 = new float[3]; //centre de gravit?nuage 2

	// *********************** INITIALISATION *******************************************
	//--Init : centroide
	for (i=0;i<3;i++)
	{
		PointG1[i] = 0.;
		PointG2[i] = 0.;
	}
	//-- Init : nombre points : lecture pour comptage des points
	while (fscanf(in1,"%f %f %f",&x,&y,&z) != EOF)
		NbrPt1++;
	while (fscanf(in2,"%f %f %f",&x,&y,&z) != EOF)
		NbrPt2++;
	//-- Init allocation memoire :  allocations des tableaux des deux ensembles de points
	Point1 = new TabFloat[NbrPt1];
	for (i=0; i<NbrPt1; i++)
		Point1[i] = new float[3];

	Point2 = new TabFloat[NbrPt2];
	for (i=0;i<NbrPt2;i++)
		Point2[i] = new float[3];

	//Init : retour au début du fichier
	rewind(in1);
	rewind(in2);

	//Init :  lecture pour la mise mémoire
	for (i=0;i<NbrPt1;i++)
		fscanf(in1,"%f %f %f",&Point1[i][0],&Point1[i][1],&Point1[i][2]);
	for (i=0;i<NbrPt2;i++)
		fscanf(in2,"%f %f %f",&Point2[i][0],&Point2[i][1],&Point2[i][2]);

	//Init :  tableau du plus proche voisin
	Voisin = new TabFloat[NbrPt1];
	for (i=0; i<NbrPt1; i++)
		{
		Voisin[i] = new float[2];// distance, Index
		Voisin[i][1] = i ;
		}
	NbrPtC = NbrPt2;
	k=0,l=0;

	// ************************** ITERATION ICP ****************************************
	// Boucle d'iteration selon l'erreur ou/et le nombre d'itération
	while((erreur>0.1)&&(k<100))
	{
	  /* Début du temps d'exécution de l'itération*/
	  QueryPerformanceFrequency(&freq);
      /* Lancement de la mesure     */
	  QueryPerformanceCounter(&start);
	  /* Affichage du compteur K    */
	  printf("L'iteration %d \n",k); k++;
	//*************************************************************************************
	//*** STEP 0 : Calcul du seuil en fonction de la distance max ou bien fixation en dure
	//*************************************************************************************
		/*max = PT_PT_DIST(Point1[0],Point2[0]);
		for (i=0; i<NbrPt1; i++)
			{
			for (j=1; j<NbrPt2; j++)
				{
				if (max< PT_PT_DIST(Point1[i],Point2[j]))
						max= PT_PT_DIST(Point1[i],Point2[j]);
				}
			}
		seuil = max / 10;*/
		seuil =500;
	//*************************************************************************************
	//*** STEP 1 : Recherche des plus proches voisins : TO DO
	//*************************************************************************************
		for (i = 0; i < NbrPt1; i++)
		{
			min = PT_PT_DIST(Point1[i], Point2[0]);
			Voisin[i][1] = 0;
			for (j = 0; j < NbrPt2; j++)
			{
				if (PT_PT_DIST(Point1[i], Point2[j]) < min)
				{
					min= PT_PT_DIST(Point1[i], Point2[j]);
					Voisin[i][1] = (float)j;
				}	
			}
			if (min > seuil)
				Voisin[i][1] = -1;
			else
				Voisin[i][0] = min;
		}
			
	//*************************************************************************************
	//*** STEP 2: CALCUL DE LA TRANSFORMEE : TO DO
	//*************************************************************************************
	//-- STEP 2.1 : Recherche des centres de masse : TO DO 
		Vector pm1(3);
		Vector pm2(3);
		int N = 0;
		Vector somme1(3);
		Vector somme2(3);
		for (i = 0; i< NbrPt1; i++)
		{
			if (Voisin[i][1] != -1)
			{
				N++;
				somme1.operator +=(Point1[i]);
				l = (int)Voisin[i][1];//index de voisin
				somme2.operator +=(Point2[l]);
			}
		}

		pm1 = somme1 / N;
		pm2 = somme2 / N;

	//-- STEP 2.2 : Calcul de H
		Vector q1(3);
		Vector q2(3);

		H.SetToZero();
		for (i = 0; i< NbrPt1; i++)
		{
			if (Voisin[i][1] != -1)
			{
				q1[0] = Point1[i][0] - pm1[0];
				q1[1] = Point1[i][1] - pm1[1];
				q1[2] = Point1[i][2] - pm1[2];
				l = (int)Voisin[i][1];
				q2[0] = Point2[l][0] - pm2[0];
				q2[1] = Point2[l][1] - pm2[1];
				q2[2] = Point2[l][2] - pm2[2];
				H = H + abT(q2, q1);
			}
		}
	//-- STEP 2.3 : Estimation de la transformation rigide (R,t)
		H.ComputeSVD(U, W, V);
		R = V * Transpose(U);
		Vector t(3);
		t = pm1 - R * pm2;
	//*************************************************************************************
	//*** STEP 3: APPLICATION DE LA TRANSFORMEE : TO DO
	//*************************************************************************************
	//--STEP 3.1 Application de la transformation : TO DO
		Vector tampon(3);
		Vector tempPoint2(3);

		for (i = 0; i<NbrPt2; i++)
		{
			tempPoint2[0] = Point2[i][0];
			tempPoint2[1] = Point2[i][1];
			tempPoint2[2] = Point2[i][2];

			tampon = R * tempPoint2;

			for (j = 0; j<3; j++)
				Point2[i][j] = tampon[j] + t[j];
		}
	//--STEP 3.2 fin de temps d'exécution

	/* Arret de la mesure         */
    QueryPerformanceCounter(&end);
	/* Conversion en millsecondes */
    elapsed = (1000.0 * (end.QuadPart - start.QuadPart)) / freq.QuadPart;
    printf("Tps exécution itéraration : %d : %.0f millisecondes entre start et end.\n",k, elapsed);
	fprintf(out4,"%.0f\n", elapsed);

	//*************************************************************************************
	//*** STEP 4: CALCUL DE L'ERREUR ET AFFICHAGE DES RESULTATS : TO DO
	//*************************************************************************************

	erreur = 0;
	for (i = 0; i< NbrPt1; i++)
	{
		if (Voisin[i][1] != -1)
		{
			l = (int)Voisin[i][1];
			erreur = erreur + SQR(PT_PT_DIST(Point1[i], Point2[l]));
		}
	}
	erreur = erreur / N;
	fprintf(out2, "%f \n", erreur);
	fprintf(out3, "Itération %d\nT : %f %f %f\n", k, t[0], t[1], t[2]);
	fprintf(out3, "R : \n");

	for (i = 0; i<3; i++)
		fprintf(out3, "%f %f %f\n\n", R(i, 0), R(i, 1), R(i, 2));
	}//fin boucle while
	// ******************************************************************

	// sauvegarde du fichier
	printf("Sauvegarde du fichier\n");
	for (i=0;i<NbrPt2;i++)
		fprintf(out1,"%f %f %f\n",Point2[i][0],Point2[i][1],Point2[i][2]);

	printf("Sauvegarde du fichier\n");

	fclose(in1);
	fclose(in2);
	fclose(out1);
	fclose(out2);
   	fclose(out3);
	fclose(out4);
	return 0;
}
