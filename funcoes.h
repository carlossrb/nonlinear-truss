#ifndef FUNCOES_H
#define FUNCOES_H

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include"nrutil.h"

#define PI 4.0 * atan(1.0)
#define MINORADOR 1.15 //Aço............................................................................../
//#define smt 300 //Suavizador da função de terceiro grau que define as tensões.........................../
#define DELTA1 0.5 //Parâmetro delta do método de newmark................................................./
#define ALFA1 0.25 //Parâmetro alfa do método de newmark................................................../
#define BEGIN 0 //Iteração inicial para os métodos de newton............................................../
/*************************Parâmetros constantes na integração numérica de newmark*************************/
#define C0 (1.0/(ALFA1*pow(deltaT,2)))
#define C1 (DELTA1/(ALFA1*deltaT))
#define C2 (1.0/(ALFA1*deltaT))
#define C3 ((1.0/(2.0*ALFA1))-1.0)
#define C4 ((DELTA1/ALFA1)-1.0)
#define C5 deltaT*((DELTA1/(2.0*ALFA1))-1.0)

/***********************************PROTÓTIPOS DAS FUNÇÕES CRIADAS****************************************/

////Calcula a matriz de rigidez geométrica................................................................/
void matriz_geom(double, double); 

////Calcula matriz de rigidez constitutiva................................................................/
void matriz_const(double);

////calculos da matriz de rigidez........................................................................./
void matrizrigidez(void);

////Matriz de rigidez global [K]........................................................................../
void matririglob(void);

////Matriz massa global e local [M]......................................................................./
void matrizmassa(void);

////Matriz de amortecimento de Rayleigh [C]=alfa[M] + beta[K]............................................./
void matrizamortecimento(void);

////Calcula as freqências naturais e seus respectivos modos de vibração (para matriz em repouso apenas).../
void frequencias(void);

////Integração no tempo por Newmark......................................................................./
void newmark(void);

////Chama as funções responsáveis pela análise dinâmica.................................................../
void dinamic(void);

////Aloca matriz de massa em coordenadas globais de cada elemento........................................./
void mmass(void);

////Chama funções, procede as iterações de newton para análise estática.................................../
void newtonestatico(void);

////desalocar memória de variáveis......................................................................../
void free_var(void);

////vetor de forças externas............................................................................../
void vectorfext(void); 

////Vetor de forças internas............................................................................../
void vectorfint(double); 

////Método da penalidade em matrizes de rigidez.........................................................../
void penal(double**); 

////Define a força residual para a análise dinâmica....................................................../
void forcaresid(void);

////Define a matriz efetiva de rigidez do sistema para análise dinâmica.................................../
void rigefetiva(void);

////Atualiza coordenadas dos vetores de coordenadas......................................................./
void refreshcord(void); 

////Zera vetores double.................................................................................../
void empvectord(double *, long , long); 

////Limpa matriz double.................................................................................../
void empmatrixd(double **, long, long, long, long);

////Atribuir valores de matriz 1 em matriz 2 (matrizes double)............................................/
void equalmatrixd(double**, double **, long, long, long, long); 

////Atribuir valores de vetor 1 wm vetor 2 (vetores double).............................................../
void equalvectord(double*, double *, long , long ); 
				
////Leitura de arquivo..................................................................................../
void leitura(void);

//// Função para chamar a condição de parada dos métodos................................................../
void condparada(void);

//// Função de saída de dados para pós processamento no GiD.............................................../
void saidadados(void);
/*********************************************************************************************************/

/*****************************************PARA O MÉTODO PBCG**********************************************/

////Chama o método PBCG em função para cálculo............................................................/

void callpbcg(void);

//Para a montagem dos vetores ija e sa da matriz de rigidez ou pseudorigidez (esparsa)..................../

extern double *sa;
extern unsigned long *ija;

//parâmetro de convergencia e números de iterações do PBCG................................................/
extern int Conv, n_iter;
/*********************************************************************************************************/

/***********************************PROTÓTIPOS DAS FUNÇÕES CHOLESKY***************************************/
//chama o método Cholesky para resolver o sistema........................................................./
void callcholesky(void); 

//Protótipos para cálculo pelo cholesky.................................................................../
void choldc(double **a, int n, double p[]);
void cholsl(double **a, int n, double p[], double b[], double x[]);
/*********************************************************************************************************/

/*************************PROTÓTIPOS DAS FUNÇÕES PBCG/VETOR DE MATRIZ ESPARSA*****************************/
////Separa a matriz esparsa em dois vetores para otimizar memória........................................./
void sprsin(double **a, unsigned long n, double thresh, unsigned long nmax, double sa[],
	unsigned long ija[]);

////Cálcula usando o PBCG a partir dos vetores de matriz esparsa........................................../
void linbcg(unsigned long n, double b[], double x[], int itol, double tol,
	int itmax, int *iter, double *err);
/*********************************************************************************************************/

/***********************************PROTÓTIPOS DAS FUNÇÕES JACOBI*****************************************/
//chama a função de rotação para encontrar os autovalores (freq. natural) e autovetores (modos de vibração)
void jacobi(double **a, int n, double d[], double **v, int *nrot);

//organiza as frequencias naturais em ordem decrescente.................................................../
void eigsrt(double d[], double **v, int n);
/*********************************************************************************************************/

/**************************DECLARAÇÃO DE VARIÁVEIS GLOBAIS PRESENTES EM var.h*****************************/
extern double **K, **kaux,**M, **mass, **C, **Kf, *Ff;
extern double *Cx, *Cy, *Cz, *auxCx, *auxCy, *auxCz, *F_x, *F_y, *F_z;
extern int n_barras, n_nos, n_mat, n_dim, *pos_i, *pos_j, *material, *R_x, *R_y, *R_z, *ft;
extern double *area, *Eac, *rho, *fy, *ko, penalty, Erro, Ti, Tf, deltaT, Csi;
extern double *Fi, *Fe, *Fr;
extern int je, countw2, psi[4][7], psi1[3][5], amort;
extern double *Fd, *sigma, *lt, *lambda, *epsilon;
extern char solve[9], deformation[11], analise[10], increm[5];
extern double Gc, GI, *l0, *bv;
extern double *V, *V1, *A, *A1, *D, *D1, alfa, beta, *w2, **autovector; //
extern unsigned long int ne, r;
extern unsigned int inc, it;
/*********************************************************************************************************/
#endif // !FUNCOES_H

