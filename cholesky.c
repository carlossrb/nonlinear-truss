#include <math.h>

//OBS: todos os 'doubles' são dados com 'floats' no numerical recipes. Mudei!

void choldc(double **a, int n, double p[]) {
	/* recebe matriz a, dimensão n, vetor vazio p.
	calcula a matriz triangular inferior (a) e os
	elementos da diagonal principal são retornados
	*/
	int i, j, k;
	double sum;
	void nrerror(char[]);

	for (i = 1; i <= n; i++) {
		for (j = i; j <= n; j++) {
			for (sum = a[i][j], k = i - 1; k >= 1; k--) sum -= a[i][k] * a[j][k];
			if (i == j) {
				if (sum <= 0.0)
					nrerror("choldc failed");
				p[i] = sqrt(sum);
			}
			else a[j][i] = sum / p[i];
		}
	}

	/* A typical use of choldc and cholslis in the inversion of covariancematrices describing the fit
	of data to a model, see, e.g., x15.6. In this, and many other applications, oe often needs L^-1. The
	Lower triangle of this matriz can be efficiently found from the output of choldc*/

	/*for (i = 1; i <= n; i++) {
	a[i][i] = 1.0 / p[i];
	for (j = i + 1; j <= n; j++) {
	sum = 0.0;
	for (k = i; k<j; k++) sum -= a[j][k] * a[k][i];
	a[j][i] = sum / p[j];
	}
	}*/
}

void cholsl(double **a, int n, double p[], double b[], double x[]) {
	/* recebe matriz a, dimensão n, vetor vazio p,
	vetor força b e vetor x que calcula o deslocamento e retorna
	*/
	int i, k;
	double sum;
	for (i = 1; i <= n; i++) {
		for (sum = b[i], k = i - 1; k >= 1; k--) sum -= a[i][k] * x[k];
		x[i] = sum / p[i];
	}
	for (i = n; i >= 1; i--) {
		for (sum = x[i], k = i + 1; k <= n; k++) sum -= a[k][i] * x[k];
		x[i] = sum / p[i];
	}

}
