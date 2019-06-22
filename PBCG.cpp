#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#define EPS 1.0e-14
//floats mudados para doubles
/*-----------------------------FUNÇÃO PRINCIPAL DO PBCG--------------------------------*/

/*x[1..n], given b[1..n], by the iterative biconjugate gradient method.
On input x[1..n] should be set to an initial guess of the solution (or all zeros); itol is 1,2,3,
or 4, specifying which convergence test is applied; itmax is the maximum number
of allowed iterations; and tol is the desired convergence tolerance. On output, x[1..n] is
reset to the improved solution, iter is the number of iterations actually taken, and err is the
estimated error. The matrix A is referenced only through the user-supplied routines atimes,
which computes the product of either A or its transpose on a vector; and asolve*/

void linbcg(unsigned long n, double b[], double x[], int itol, double tol,
	int itmax, int *iter, double *err)
{
	/*MUDEI: entra-se com sa e ija também na função para não usar Extern*/
	/*em tol independente da tolernancia dada, o que ditará a convergência é o 
	numero máx (nmax) de iterações*/
	
	void asolve(unsigned long n, double b[], double x[], int itrnsp);
	void atimes(unsigned long n, double x[], double r[], int itrnsp);
	double snrm(unsigned long n, double sx[], int itol);
	unsigned long j;
	double ak, akden, bk, bkden, bknum, bnrm, dxnrm, xnrm, zm1nrm, znrm;
	double *p, *pp, *r, *rr, *z, *zz;
	p = dvector(1, n);
	pp = dvector(1, n);
	r = dvector(1, n);
	rr = dvector(1, n);
	z = dvector(1, n);
	zz = dvector(1, n);
	*iter = 0;
	atimes(n, x, r, 0);
	for (j = 1; j <= n; j++) {
		r[j] = b[j] - r[j];
		rr[j] = r[j];
	}
	/* atimes(n,r,rr,0); */
	if (itol == 1) {
		bnrm = snrm(n, b, itol);
		asolve(n, r, z, 0);
	}
	else if (itol == 2) {
		asolve(n, b, z, 0);
		bnrm = snrm(n, z, itol);
		asolve(n, r, z, 0);
	}
	else if (itol == 3 || itol == 4) {
		asolve(n, b, z, 0);
		bnrm = snrm(n, z, itol);
		asolve(n, r, z, 0);
		znrm = snrm(n, z, itol);
	}
	else nrerror("illegal itol in linbcg");
	while (*iter <= itmax) {
		++(*iter);
		asolve(n, rr, zz, 1);
		for (bknum = 0.0, j = 1; j <= n; j++) bknum += z[j] * rr[j];
			if (*iter == 1) {
				for (j = 1; j <= n; j++) {
					p[j] = z[j];
					pp[j] = zz[j];
				}
			}
			else {
				bk = bknum / bkden;
				for (j = 1; j <= n; j++) {
					p[j] = bk*p[j] + z[j];
					pp[j] = bk*pp[j] + zz[j];
				}
			}
			bkden = bknum;
			atimes(n, p, z, 0);
			for (akden = 0.0, j = 1; j <= n; j++) akden += z[j] * pp[j];
			ak = bknum / akden;
			atimes(n, pp, zz, 1);
			for (j = 1; j <= n; j++) {
				x[j] += ak*p[j];
				r[j] -= ak*z[j];
				rr[j] -= ak*zz[j];
			}
			asolve(n, r, z, 0);
			if (itol == 1)
				*err = snrm(n, r, itol) / bnrm;
			else if (itol == 2)
				*err = snrm(n, z, itol) / bnrm;
			else if (itol == 3 || itol == 4) {
				zm1nrm = znrm;
				znrm = snrm(n, z, itol);
				if (fabs(zm1nrm - znrm) > EPS*znrm) {
					dxnrm = fabs(ak)*snrm(n, p, itol);
					*err = znrm / fabs(zm1nrm - znrm)*dxnrm;
				}
				else {
					*err = znrm / bnrm;
					continue;
				}
				xnrm = snrm(n, x, itol);
				if (*err <= 0.5*xnrm) *err /= xnrm;
				else {
					*err = znrm / bnrm;
					continue;
				}
			}
			printf("iter=%4d err=%12.6f\n", *iter, *err);
			if (*err <= tol) break;
	}
	free_dvector(p, 1, n);
	free_dvector(pp, 1, n);
	free_dvector(r, 1, n);
	free_dvector(rr, 1, n);
	free_dvector(z, 1, n);
	free_dvector(zz, 1, n);
}


/*--------------------------PARA TRANSFORM MATRIZ EM VETOR-----------------------------*/

/*Converts a square matrix a[1..n][1..n] into row-indexed sparse storage mode.
Only elements of a with magnitude ≥thresh are retained. Output is in two linear arrays
with dimension nmax (an input parameter): sa[1..] contains array values, indexed by ija[1..]. The
number of elements filled of sa and ija on output are both ija[ija[1]-1]-1.
Location 1 of ija is always equal to N + 2. (It can be read to determine N.). ija[1]-1=ija[1]+n*/


/*IMPORTANTE P/ nmax: Para tamanho da matriz definido, tem-se tamanho de ija e sa dado como 2*A(tamanho)-1.
Para numero de elementos definidos como tamanho da matriz, tem-se 2*N-1*/

//mudei os floats para doubles em **a, thresh e sa[]. 
//Mudei int n para unsigned long

void sprsin(double **a, unsigned long n, double thresh, unsigned long nmax, double sa[],
	unsigned long ija[])
{
	void nrerror(char error_text[]); //protótipo - função em nrutil.cpp
	int i, j;
	unsigned long k;
	for (j = 1; j <= n; j++) sa[j] = a[j][j];
	ija[1] = n + 2;
	k = n + 1;
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) {
			if (fabs(a[i][j]) >= thresh && i != j) {
				if (++k > nmax) nrerror("sprsin: nmax too small");
				sa[k] = a[i][j];
				ija[k] = j;
			}
		}
		ija[i + 1] = k + 1;
	}

}


/*-----------------------------FUNÇÕES INTERNAS DO PBCG--------------------------------*/

/*So that the specifications for the routines atimes and asolve are clear, we list here
simple versions that assume a matrix A stored somewhere in row-index sparse format.*/

extern double *sa;
extern unsigned long *ija;

double snrm(unsigned long n, double sx[], int itol)
{
	/*Compute one of two norms for a vector sx[1..n], as signaled by itol. Used by linbcg*/
	unsigned long i, isamax;
	double ans;
	if (itol <= 3) {
		ans = 0.0;
		for (i = 1; i <= n; i++) ans += sx[i] * sx[i];
		return sqrt(ans);
	}
	else {
		isamax = 1;
		for (i = 1; i <= n; i++) {
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax = i;
		}
		return fabs(sx[isamax]);
	}

}

void atimes(unsigned long n, double x[], double r[], int itrnsp)
{

	void dsprsax(double sa[], unsigned long ija[], double x[], double b[],
		unsigned long n);
	void dsprstx(double sa[], unsigned long ija[], double x[], double b[],
		unsigned long n);

	if (itrnsp) dsprstx(sa, ija, x, r, n);
	else dsprsax(sa, ija, x, r, n);
}

void asolve(unsigned long n, double b[], double x[], int itrnsp)
{


	/*The matrix Ae is the diagonal part of A, stored in the first n elements of sa.Since the
	transpose matrix has the same diagonal, the flag itrnsp is not used.*/

	unsigned long i;
	for (i = 1; i <= n; i++) x[i] = (sa[i] != 0.0 ? b[i] / sa[i] : b[i]);
}


/*-----------------------------RESOLVER EM FUNÇÃO ATIMES--------------------------------*/

void dsprsax(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n)
{
	/*Para utilização em atimes. Convertido manualmente de float para double*/
	/*Multiply a matrix in row-index sparse storage arrays sa and ija by a vector x[1..n], giving
	a vector b[1..n].*/

	void nrerror(char error_text[]); //protótipo - função em nrutil.cpp
	unsigned long i, k;
	if (ija[1] != n + 2) nrerror("sprsax: mismatched vector and matrix");
	for (i = 1; i <= n; i++) {
		b[i] = sa[i] * x[i];
		for (k = ija[i]; k <= ija[i + 1] - 1; k++)
			b[i] += sa[k] * x[ija[k]];
	}
}

void dsprstx(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n)
{
	/*Para utilização em atimes. Convertido manualmente de float para double*/
	/*Multiply the transpose of a matrix in row - index sparse storage arrays sa and ija by a vector
	x[1..n], giving a vector b[1..n].*/

	void nrerror(char error_text[]); //protótipo - função em nrutil.cpp
	unsigned long i, j, k;
	if (ija[1] != n + 2) nrerror("mismatched vector and matrix in sprstx");
	for (i = 1; i <= n; i++) b[i] = sa[i] * x[i];
	for (i = 1; i <= n; i++) {
		for (k = ija[i]; k <= ija[i + 1] - 1; k++) {
			j = ija[k];
			b[j] += sa[k] * x[i];
		}
	}
}