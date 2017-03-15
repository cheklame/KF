#include <stdio.h>
#include <stdlib.h>


#define Nb_element 3





void Inv_Mat(float Mat[Nb_element][Nb_element], float Id_Mat[Nb_element][Nb_element]){

int i,j,k;
float p;



for(k = 0; k<Nb_element; k++){
	if(Mat[k][k] == 0){
		printf("error \n");
	}
    else{
	  p = Mat[k][k];

	//normalisation

    for(j = 0; j<Nb_element; j++){
    	Mat[k][j] = Mat[k][j] / p;
        Id_Mat[k][j] = Id_Mat[k][j] / p;
    }


    //reduction
     for(i = 0; i<Nb_element; i++){
     	if(i!=k){
     		p = Mat[i][k];
     		for(j=0; j<Nb_element; j++){
     			Mat[i][j] -= p*Mat[k][j];
     		    Id_Mat[i][j] -= p*Id_Mat[k][j];}
     	}
     }

}
}

}




void Mult_Mat(float M1[Nb_element][Nb_element], float M2[Nb_element][Nb_element], float Res[Nb_element][Nb_element]){

int i,j,k;

for(i = 0; i<Nb_element;i++){

	for(j = 0; j<Nb_element; j++){

      Res[i][j] = 0;

	}
}

for(i = 0; i<Nb_element;i++){

	for(j = 0; j<Nb_element; j++){

		for(k = 0; k<Nb_element; k++){

		Res[i][j] += M1[i][k] * M2[k][j];

        }

    }
}
}



void Tran_Mat(float Mat[Nb_element][Nb_element]){

  int i,j;
  float temp[Nb_element][Nb_element];

  for(i=0; i<Nb_element; i++){
    for(j=0; j<Nb_element; j++){

     temp[j][i] = Mat[i][j];

    }
}
   for(i = 0; i<Nb_element;i++){
	 for(j = 0; j<Nb_element; j++){

      Mat[i][j] = temp[i][j];
	}
}

}


void Add_Mat(float M1[Nb_element][Nb_element],float M2[Nb_element][Nb_element],float Res[Nb_element][Nb_element]){

int i,j;

     for(i = 0; i<Nb_element;i++){
     	for(j = 0; j<Nb_element; j++){

     		Res[i][j] = M1[i][j] + M2[i][j];
     	}
     }
}






void Display(float M[Nb_element][Nb_element]){

for(i = 0; i<Nb_element; i++){
	for(j = 0; j<Nb_element; j++){

      printf("  %f  ",M[i][j]);

	}

	printf("\n");
}
printf("\n\n");


}


void Sous_Mat(float M1[Nb_element][Nb_element],float M2[Nb_element][Nb_element],float Res[Nb_element][Nb_element]){

int i,j;

     for(i = 0; i<Nb_element;i++){
     	for(j = 0; j<Nb_element; j++){

     		Res[i][j] = M1[i][j] - M2[i][j];
     	}
     }
}







int main(){
int i,j;


float Mat1[Nb_element][Nb_element] = {{2,1,-4} , {3, 3, -5} , {4, 5, -2}};
float Mat2[Nb_element][Nb_element] = {{2,1,-4} , {3, 3, -5} , {4, 5, -2}};
float Res[Nb_element][Nb_element];


float Id_Mat[Nb_element][Nb_element] = {{1,0,0},{0,1,0},{0,0,1}};


Inv_Mat(Mat2, Id_Mat);
Display(Inv_Mat);

//----------------------------------------------------------------------

Mult_Mat(Mat1,Mat2,Res);
Display(Res);


//-----------------------------------------------

Tran_Mat(Mat1);
Display(Mat1);

//------------------------------------------

Add_Mat(Mat1,Mat2,Res);
Display(Res);

//-------------------------------------------

Sous_Mat(Mat1,Mat2,Res);
Display(Res);


	return 0;
}
