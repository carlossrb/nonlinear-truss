#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_
/*DIFERENTE DO NUMERICAL RECIPES. FUNÇŐES USADAS PARA ALOCAÇŐES DE MEMÓRIA DINÂMICA
ÚTIL NO PBCG, ENTRE OUTROS*/

double *dvector(long nl, long nh);
double *dvector_zeros(long nl, long nh);
unsigned int *uivector(unsigned nl, unsigned nh);
unsigned int *uivector_zeros(unsigned nl, unsigned nh);
unsigned long *lvector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
int *ivector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix_zeros(long nrl, long nrh, long ncl, long nch);//matrix of zeros
int **imatrix(long nrl, long nrh, long ncl, long nch);
void nrerror(char[]);
void free_dvector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_uivector(unsigned  *v, unsigned nl, unsigned nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);

#endif /* _NR_UTILS_H_ */
