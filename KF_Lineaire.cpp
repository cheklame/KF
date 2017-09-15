/*
  * KF_Lineaire.cpp
 *
 *  Created on: 9 mai 2017
 *      Author: cheklat
 */

#include "KF.h"



void  KF_Lineaire(TYPE_DATA Y[ROW_COL_M],
				  TYPE_DATA X_in[ROW_COL_M],
				  TYPE_DATA X_out[ROW_COL_M],
				  TYPE_DATA P_in[ROW_COL_M][ROW_COL_M],
				  TYPE_DATA P_out[ROW_COL_M][ROW_COL_M],
				  TYPE_DATA Err[ROW_COL_M],
				  TYPE_DATA Tx,
				  TYPE_DATA Ty,
				  TYPE_DATA sigma_capt,
				  TYPE_DATA sigma_sys){



	// X = AX + B*U + M*W

	TYPE_DATA A[ROW_COL_M][ROW_COL_M];
	TYPE_DATA C[ROW_COL_M][ROW_COL_M];
	TYPE_DATA K[ROW_COL_M][ROW_COL_M];
	TYPE_DATA R[ROW_COL_M][ROW_COL_M];
	TYPE_DATA Q[ROW_COL_M][ROW_COL_M];
	TYPE_DATA M[ROW_COL_M][ROW_COL_M];
	TYPE_DATA w[ROW_COL_M];
	TYPE_DATA v[ROW_COL_M];



	a1:for(int i = 0; i<ROW_COL_M; i++){
		a2:for(int j= 0; j<ROW_COL_M; j++){
			if(i==j){
				R[i][j] = sigma_sys*sigma_sys;   //bruit de mesure
				Q[i][j] = sigma_capt;               // bruit de capteurs
				K[i][j] = 0;
				A[i][j] = 1;
				C[i][j] = 1;
				M[i][j] = 1;

			}
			else{
			    R[i][j] = 0;
				Q[i][j] = 0;
				K[i][j] = 0;
				A[i][j] = 0;
				C[i][j] = 0;
				M[i][j] = 0;}
		}
	}


	A[0][2] = Tx;
	A[1][3] = Ty;
	A[4][4] = 0;
	A[5][5] = 0;

	//*****begining of algorithme of KF

// phase de prédiction

TYPE_DATA X_inter[ROW_COL_M];

 mvmult(A,X_in,X_inter);


// -----  P = A*P*A' + M*Q*M^T   -----------

//#ifndef __SYNTHESIS__
	//cette partie ne sera pas syntethisé mettre les affichages
//#endif

 TYPE_DATA ATran[ROW_COL_M][ROW_COL_M];
 TYPE_DATA CTran[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Res_AP[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Res_APAtran[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Inv[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Res_CX[ROW_COL_M];
 TYPE_DATA Res_Y_CX[ROW_COL_M];
 TYPE_DATA Res_K_Y_CX[ROW_COL_M];
 TYPE_DATA Res_PC[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Res_CtranPC[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Res_CtranPC_R[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Res_KPC[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Res_P_KCP[ROW_COL_M][ROW_COL_M];
 TYPE_DATA P_inter[ROW_COL_M][ROW_COL_M];
 TYPE_DATA MTran[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Res_MQ[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Res_MQMtran[ROW_COL_M][ROW_COL_M];
 TYPE_DATA Res_CtranInv[ROW_COL_M][ROW_COL_M];
 int success;


 mtran(A,ATran);
 mmult(A,P_in,Res_AP);
 mmult(ATran,Res_AP,Res_APAtran);

 mtran(M,MTran);
 mmult(M,Q,Res_MQ);
 mmult(MTran,Res_MQ,Res_MQMtran);

 madd(Res_APAtran,Res_MQMtran,P_inter);


//phase de mise à jour

	//K = P*C^t[C*P*C^t +  R]^-1;
	 mtran(C,CTran);                   // trasposée de c

	 mmult(C,P_inter,Res_PC);                // Ctran * p
	 mmult(CTran,Res_PC,Res_CtranPC);  //// Ctran * p * C
	 madd(Res_CtranPC,R,Res_CtranPC_R); // Ctran * p * C + R


	 cout<<"matrice a inversé  " << endl;
	 success = KF_cholesky_inverse(Res_CtranPC_R,Inv);  //(Ctran * p * C + R)^-1
	 cout<<"inv ca marche" << endl;

	 mmult(CTran,Inv,Res_CtranInv);   // Ctran*(Ctran * p * C + R)^-1
	 mmult(Res_CtranInv,P_inter,K);  // P*Ctran*(Ctran * p * C)^-1


    // X = X + K(Y - C*X)
    mvmult(C,X_inter,Res_CX);        //c*X
    sousvec(Y,Res_CX,Res_Y_CX); // Y - C*X
    for(int r=0;r<ROW_COL_M;r++)
    Err[r] = Res_Y_CX[r];
    mvmult(K,Res_Y_CX,Res_K_Y_CX);   // K(Y - C*X)
    addvec(X_inter,Res_K_Y_CX,X_out);          // X + K(Y - C*X)


    //P =  P - K*C*P
    mmult(K,Res_PC,Res_KPC); //K*C*P
    msous(P_inter,Res_KPC,P_out);     //P - K*C*P

}





//*****************KF_Lineaire avec des entrées en streaming**************//

void KF_Lineaire_HLS(  hls::stream<AXI_VALUE> &in_stream,
		               hls::stream<AXI_VALUE> &out_stream){

#pragma HLS INTERFACE axis port=in_stream
#pragma HLS INTERFACE axis port=out_stream
#pragma HLS INTERFACE s_axilite port=return bundle=CRTL_BUS



	TYPE_DATA P[ROW_COL_M][ROW_COL_M];
	TYPE_DATA X[ROW_COL_M];
	TYPE_DATA P_out[ROW_COL_M][ROW_COL_M];
	TYPE_DATA X_out[ROW_COL_M];
	TYPE_DATA Y[ROW_COL_M];
	TYPE_DATA Err[ROW_COL_M];
	TYPE_DATA Ts[2];
	TYPE_DATA Sigmas[2];


	AXI_VALUE aValue1;
	AXI_VALUE aValue2;
	AXI_VALUE aValue3;
	AXI_VALUE aValue4;
	AXI_VALUE aValue5;
	AXI_VALUE aValue6;



	//récupération du vecteur de mesures Y
	a1: for(int i=0; i<ROW_COL_M; i++){
#pragma HLS PIPELINE
		in_stream.read(aValue1);
		Y[i] = (TYPE_DATA) aValue1.data >> BIT_SHIFT;
		}

	//récupération du vecteur d'état précédent X
	b1: for(int i=0; i<ROW_COL_M; i++){
#pragma HLS PIPELINE
		in_stream.read(aValue2);
		X[i] = (TYPE_DATA) aValue2.data >> BIT_SHIFT;
		}

	//récuperation de la matrice P précédente

	c1: for(int i=0; i<ROW_COL_M; i++){
		c2: for(int j=0; j<ROW_COL_M; j++){
#pragma HLS PIPELINE
			in_stream.read(aValue3);
			P[i][j] = (TYPE_DATA) aValue3.data >> BIT_SHIFT;

		}
	}

	//récupération du vecteur qui contient les valeurs des transitions Tx et Ty
	e1: for(int i=0; i<2; i++){
#pragma HLS PIPELINE
		in_stream.read(aValue5);
		Ts[i] = (TYPE_DATA) aValue5.data >> BIT_SHIFT;
	    }

	//récupération du vecteur qui contient les valeurs des sigmas Sigma_capt et sigma_sys
	f1: for(int i=0; i<2; i++){
#pragma HLS PIPELINE
		in_stream.read(aValue6);
		Sigmas[i] =  (TYPE_DATA) aValue6.data >> BIT_SHIFT;
		}


	//appel de la fonction KF_Liniaire

	KF_Lineaire(Y,X,X_out,P,P_out,Err,Ts[0],Ts[1],Sigmas[0],Sigmas[1]);



	//mettre le vecteur X en sortie

	 R1: for( int i=0; i<ROW_COL_M; i++){
#pragma HLS PIPELINE
		       int  temp = (int)(X_out[i] << BIT_SHIFT);
			   aValue2.data = temp;
			   aValue2.user = (i==0) ? 1 : 0;
			   aValue2.strb = -1;
			   aValue2.keep = 15;
			   aValue2.id = 0;
			   aValue2.dest = 0;
			   aValue2.last = (i==ROW_COL_M -1) ? 1 : 0;
			   out_stream.write(aValue2);
		   }

	//mettre P en sortie

	Res1: for(int l=0; l<ROW_COL_M; l++){
	   Res2: for( int h=0; h<ROW_COL_M; h++){
#pragma HLS PIPELINE
		   int  temp = (int) (P_out[l][h] << BIT_SHIFT);
		   aValue3.data = temp;
		   aValue3.user = ((l==0) && (h==0))  ? 1 : 0;
		   aValue3.strb = -1;
		   aValue3.keep = 15;
		   aValue3.id = 0;
		   aValue3.dest = 0;
		   aValue3.last = ((l==ROW_COL_M -1) && (h==ROW_COL_M-1)) ? 1 : 0;
		   out_stream.write(aValue3);

	   }
	}

	//mettre le vecteur Err en sortie

	L1: for( int i=0; i<ROW_COL_M; i++){
#pragma HLS PIPELINE
		     int  temp = (int)(Err[i] << BIT_SHIFT);
	           aValue4.data = temp;
			   aValue4.user = (i==0) ? 1 : 0;
			   aValue4.strb = -1;
			   aValue4.keep = 15;
			   aValue4.id = 0;
			   aValue4.dest = 0;
			   aValue4.last = (i==ROW_COL_M -1) ? 1 : 0;
			   out_stream.write(aValue4);

}

}
