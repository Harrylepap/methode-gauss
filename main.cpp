#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <fstream>
using namespace std;

void ** getData(float**& A, float *& b, float *&x, int &dim);
void displayResult(float**A, float* b, float* x,int dim);
void erreur (string message);
void displayMat(float** A, int dimL, int dimC);
void displayVec(float *b, int dim);
void gauss(float**&A, float*& b, float*& x, int dim);
void triangMat(float **&A, float *&b, int dim);
void solveTriangsup(float** A, float* b, float*& x, int dim);
void permuter(float**& A, float*& b, int dim);
void getFileData(float**& A, float *& b, float *&x, int &dim);

int main(){
    cout << "Programme de test resolution equation" << endl;
    /// Donnees
    float **A(0), *b(0), *x(0);
    int dim(4);
    //Allocation
    getFileData(A,b,x,dim);

    /// Calcul
    gauss(A,b,x,dim);

    /// Resultat
    displayResult(A,b,x,dim);
    return 0;
}
void gauss(float**&A, float*& b, float*& x, int dim){
    /// Triangularisation du systeme
    triangMat(A,b,dim);
    /// Resolution du systeme triangulaire superieur
    solveTriangsup(A,b,x,dim);
}
void solveTriangsup(float** A, float* b, float*& x, int dim){
    int i(0),j(0);
    float s(0);
    for(i=dim-1; i>=0; i--){
        for(s=0, j=i+1; j<dim ;j++){
            s+=A[i][j]*x[j];
        }
        x[i]=(b[i]-s)/A[i][i];
    }
}
void triangMat(float **&A, float *&b, int dim){
    int i(0),j(0),k(0),lp(0);
    float p(0),*tp(0),t(0);
    for(k=0;k<dim;k++){
    //Recherche d'un plus grand pivot
    lp=k;
    p=fabs(A[k][k]);
    for(i=k+1;i<dim;i++){
        if(fabs(A[i][k])>p){
            p=fabs(A[i][k]); lp=i;
        }
    }

    // Permutation des lignes k et lp
    tp=A[k]; A[k]=A[lp]; A[lp]=tp;
    // Permutation du second membre
    t=b[k]; b[k]=b[lp]; b[lp]=t;
    /// Elimination
        for(i=k+1;i<dim;i++){
            for(j=k+1;j<dim;j++){
                A[i][j]-=(A[i][k]/A[k][k]*A[k][j]);
            }
            b[i]-=(A[i][k]/A[k][k]*b[k]);
            A[i][k]=0;
        }
    }
}
void erreur (string message){
    cout << "Erreur: " << message<< endl;
    exit(1);
}
void displayResult(float**A, float* b, float* x,int dim){
    cout << "\nApres triangularisation" << endl;
    cout << "\nLa matrice du probleme : " << endl;
    displayMat(A,dim,dim);
    cout << "\nLe second membre : " << endl;
    displayVec(b,dim);
    cout << "\nLa solution trouvee : " << endl;
    displayVec(x,dim);
}
void displayMat(float** A, int dimL, int dimC){
    for(int i(0);i<dimL;i++){
        for(int j(0);j<dimC;j++){
            cout << "\t" << A[i][j];
        }
        cout << endl;
    }
}
void displayVec(float *b, int dim){
    for(int i(0);i<dim;i++){
        cout << "\t" << b[i] << endl;
    }
}
void ** getData(float**& A, float *& b, float *&x, int &dim){
    /// Creation des matrices et vecteurs
    A=new (nothrow) float *[dim];
    if(A){
        for(int i(0);i<dim;i++){
            A[i]=new (nothrow) float [dim];
            if(!A[i]) erreur("Probleme d'allocation de A[i]...");
        }
    }
    else erreur("Probleme d'allocation de A...");

    b=new (nothrow) float[dim];
    if(!b) erreur("Probleme d'allocation de b...");
    x=new (nothrow) float[dim];
    if(!x) erreur("Probleme d'allocation de x...");
    /// Remplissage des donnees
    /*A[0][0]=1;A[0][1]=2;A[0][2]=3;
    A[1][0]=7;A[1][1]=8;A[1][2]=9;
    A[2][0]=2;A[2][1]=4;A[2][2]=5;
    b[0]=4;b[1]=29;b[2]=9;
    x[0]=0;x[1]=1;x[2]=2;*/

    A[0][0]=2;A[0][1]=1;A[0][2]=0;A[0][3]=4;
    A[1][0]=-4;A[1][1]=-2;A[1][2]=3;A[1][3]=-7;
    A[2][0]=4;A[2][1]=1;A[2][2]=-2;A[2][3]=8;
    A[3][0]=0;A[3][1]=-3;A[3][2]=-12;A[3][3]=-1;
    b[0]=2;b[1]=-9;b[2]=2;b[3]=2;
    x[0]=0;x[1]=0;x[2]=0;x[3]=0;

}
void getFileData(float**& A, float *& b, float *&x, int &dim){
    ifstream pfile(0);
    int i(0),j(0);
    // Creation des matrices et vecteurs
    pfile.open("data.txt");
    if(!pfile) erreur("Probleme d'ouverture du fichier...");
    pfile >> dim;
    A=new (nothrow) float *[dim];
    if(A){
        for(int i(0);i<dim;i++){
            A[i]=new (nothrow) float [dim];
            if(!A[i]) erreur("Probleme d'allocation de A[i]...");
        }
    }
    else erreur("Probleme d'allocation de A...");

    b=new (nothrow) float[dim];
    if(!b) erreur("Probleme d'allocation de b...");
    x=new (nothrow) float[dim];
    if(!x) erreur("Probleme d'allocation de x...");

    for(i=0;i<dim;i++) x[i]=0;
    /// Remplissage des donnees
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            pfile >> A[i][j];
        }
    }
    for(i=0;i<dim;i++) pfile >> b[i];
    pfile.close();
}
