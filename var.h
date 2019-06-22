#ifndef VAR_H
#define VAR_H

/****************************************************VARIÁVEIS GLOBAIS*********************************************************/
/*********************************************Variáveis da entrada de dados****************************************************/

int n_barras; //Número de elementos............................................................................................/
int n_nos; //Número de nós...................................................................................................../
int n_mat; //Número de materiais.............................................................................................../
int n_dim; //Número de dimensões do problema.................................................................................../
int *material; //Condição de material em cada barra............................................................................/
double *Cx; //Coordenadas x no nó............................................................................................../
double *Cy; //Coordenadas y no nó............................................................................................../
double *Cz; //Coordenadas z no nó............................................................................................../
int *pos_i; //Posição do nó (inicial)........................................................................................../
int *pos_j; //Posição do nó (final)............................................................................................/
double *F_x; //Força na direção x............................................................................................../
double *F_y; //Força na direção y............................................................................................../
double *F_z; //Força na direção z............................................................................................../
int *ft; //tipo de comportamento da carga aplicada (1,2,3,4,5)................................................................./
double *ko; //Constante multiplicativa........................................................................................./
int *R_x; //Restrição na direção x............................................................................................./
int *R_y; //Restrição na direção y............................................................................................./
int *R_z; //Restrição na direção z............................................................................................./
double *area; //Área da seção transversal....................................................................................../
double *Eac; //Módulo de elasticidade........................................................................................../
double *fy; //Tensão de escoamento (sem minorador)............................................................................./
double penalty; //Função penalidade............................................................................................/
double Erro; //Limites de Erros................................................................................................/
int Conv, n_iter, amort; /*Parâmetro de convergência para o PBCG, número de iterações mínimas para o método,
definição se terá ou não amortecimento*///...................................................................................../
char solve[9], deformation[11], analise[10], increm[5]; //solve(PBCG ou Cholesky), análise (din ou est), incremental-iterativo?/
unsigned int inc; //Quantidade de incrementos para análise incremental-iterativa.............................................../
/********************************************Variáveis gerais usadas em funções************************************************/

double *auxCx, *auxCy, *auxCz;//Vetores auxiliares para coordenadas............................................................/
double **K; //Matriz de rigidez global........................................................................................./
double **kaux; //Matrizes auxiliares para cálculo............................................................................../
int je; //Importante variável para interação (p/ fazer matriz de rig de cada barra)............................................/
double *Fe, *Fi, *Fd, *Fr; //Vetores forças internas e externas e diferença de forças........................................../
double *D, *sigma, *lt,*l0, *lambda; //Incremento de deslocamento, tensão, comprimento corrente, comprimento inicial, estiramento
unsigned long int r;//Contador aux.........................................................................................................../
unsigned long int ne=0; //contador de newmark;................................................................................./
unsigned int it; //contador do newton incremental iterativo..................................................................../
double *epsilon; //deformações nas barras....................................................................................../
double *D1; //Vetor deslocamento.............................................................................................../
double Gc, GI; //Condição de parada (corrente e inicial)......................................................................./
double *bv; //vetor das coordenadas (xi-xf,yi-yf,zi-zf)......................................................................../

//Matriz estática psi (posição 0 desconsiderada para cálculos)................................................................./
//para estruturas 3D
int psi[4][7] = { 
{ 0,0,0,0,0,0,0 },
{ 0,-1,0,0,1,0,0 }, 
{ 0,0,-1,0,0,1,0 },
{ 0,0,0,-1,0,0,1 } 
};
//para estruturas 2D
int psi1[3][5] = {
	{ 0,0,0,0,0 },
	{ 0,-1,0,1,0 },
	{ 0,0,-1,0,1}
	};
/*************************************Variáveis unicamente usadas na análise dinâmica******************************************/


double **M; //matriz de massa global [M]......................................................................................./
double **mass; //Matriz de massa elementar em coordenadas globais............................................................../
double **C; //matriz de amortecimento global [C].............................................................................../
double **Kf; //pseudomatriz de rigidez [K*]..................................................................................../
double *Ff;//pseudoforçainterna................................................................................................/
double *rho; //densidade do material.........................................................................................../
double Ti, Tf, deltaT, Csi; //tempo inicial, tempo final, intervalo de integração, taxa de amortecimento crítico.............../
double *V,*V1, *A, *A1; //incrementos de velocidade e aceleração e vetores respectivos........................................./
double alfa, beta, *w2,**autovector ; /*parâmetros para encontrar matriz de amortecimento, frequencias naturais e 
                                       modos de vibração*/ //................................................................../
int countw2; //contador da quantidade de frequências naturais................................................................../

/**********************************************PARA RESOLUÇÃO DE SISTEMAS LINEARES*********************************************/

//vetores ija e sa - para resolução de sistema (PBCG)........................................................................../
double *sa; //aloca valores desta matriz nestas posições referentes à *ija...................................................../
unsigned long *ija; //aloca posições da matriz a ser resolvida................................................................./
#endif // !VAR_H


