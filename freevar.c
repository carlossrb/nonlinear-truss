#include"var.h"
#include"nrutil.h"
#include <string.h>
/******************LIBERA ESPAÇO DE MEMÓRIA DOS VETORES/MATRIZES ALOCADOS DINAMICAMENTE******************/

void free_var(void) {

	//Caso seja análise dinâmica desaloca essas variáveis
	if (strcmp(analise, "Dinâmica") == 0) {
		free_dvector(bv, 1, 3);

		free_ivector(R_x, 1, n_nos);
		free_ivector(R_y, 1, n_nos);
		free_ivector(R_z, 1, n_nos);
		free_dvector(F_x, 1, n_nos);
		free_dvector(F_y, 1, n_nos);
		free_dvector(F_z, 1, n_nos);

		free_ivector(pos_i, 1, n_barras);
		free_ivector(pos_j, 1, n_barras);
		free_ivector(material, 1, n_barras);

		free_dvector(Cx, 1, n_nos);
		free_dvector(Cy, 1, n_nos);
		free_dvector(Cz, 1, n_nos);
		free_dvector(auxCx, 1, n_nos);
		free_dvector(auxCy, 1, n_nos);
		free_dvector(auxCz, 1, n_nos);
		free_dvector(Fd, 1, n_dim * n_nos);
		free_dvector(Fi, 1, n_dim * n_nos);
		free_dvector(D, 1, n_dim * n_nos);
		free_dvector(V, 1, n_dim * n_nos);
		free_dvector(V1, 1, n_dim * n_nos);
		free_dvector(A, 1, n_dim * n_nos);
		free_dvector(A1, 1, n_dim * n_nos);
		free_dvector(D1, 1, n_dim * n_nos); // ficar atento a isso aqui
		free_dvector(w2, 1, n_dim * n_nos);
		free_dvector(Fe, 1, n_dim * n_nos);

		free_dvector(epsilon, 1, n_barras);
		free_dvector(sigma, 1, n_barras);
		free_dvector(area, 1, n_barras);
		free_dvector(rho, 1, n_barras);
		free_dvector(fy, 1, n_barras);
		free_dvector(Eac, 1, n_barras);
		free_dvector(lt, 1, n_barras);
		free_dvector(l0, 1, n_barras);
		free_dvector(lambda, 1, n_barras);


		free_dmatrix(kaux, 1, 2 * n_dim, 1, 2 * n_dim);
		free_dmatrix(K, 1, n_dim * n_nos, 1, n_dim*n_nos);
		free_dmatrix(M, 1, n_dim * n_nos, 1, n_dim* n_nos);
		free_dmatrix(C, 1, n_dim * n_nos, 1, n_dim * n_nos);
		free_dmatrix(autovector, 1, n_dim * n_nos, 1, n_dim * n_nos);

	}
	//Caso seja análise Estática desaloca essas variáveis
	else {
		free_dmatrix(kaux, 1, 2 * n_dim, 1, 2 * n_dim);
		free_dmatrix(K, 1, n_dim * n_nos, 1, n_dim*n_nos);

		free_dvector(epsilon, 1, n_barras);
		free_dvector(sigma, 1, n_barras);
		free_dvector(area, 1, n_barras);
		free_dvector(rho, 1, n_barras);
		free_dvector(fy, 1, n_barras);
		free_dvector(Eac, 1, n_barras);
		free_dvector(lt, 1, n_barras);
		free_dvector(l0, 1, n_barras);
		free_dvector(lambda, 1, n_barras);

		free_dvector(Cx, 1, n_nos);
		free_dvector(Cy, 1, n_nos);
		free_dvector(Cz, 1, n_nos);
		free_dvector(auxCx, 1, n_nos);
		free_dvector(auxCy, 1, n_nos);
		free_dvector(auxCz, 1, n_nos);
		free_dvector(Fd, 1, n_dim * n_nos);
		free_dvector(Fi, 1, n_dim * n_nos);
		free_dvector(D, 1, n_dim * n_nos);

		free_dvector(bv, 1, 3);

		free_ivector(R_x, 1, n_nos);
		free_ivector(R_y, 1, n_nos);
		free_ivector(R_z, 1, n_nos);
		free_dvector(F_x, 1, n_nos);
		free_dvector(F_y, 1, n_nos);
		free_dvector(F_z, 1, n_nos);

		free_ivector(pos_i, 1, n_barras);
		free_ivector(pos_j, 1, n_barras);
		free_ivector(material, 1, n_barras);
	}
}

