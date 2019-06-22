#include"funcoes.h"


void saidadados(void) {

	FILE *output1;
	char filedat[1024], filelog[1024];
	extern char projname[1024];

	strcpy(filedat, projname);
	strcat(filedat, ".post.res");
	
	output1 = fopen(filedat, "a");
	if (strcmp(analise, "Estática") == 0) {
		if (it == 1) {
			fprintf(output1, "GiD Post Results File 1.0\n\n");
			// GaussPoints

			fprintf(output1, "GaussPoints ");
			fprintf(output1, "\"Gauss Points\" ");
			fprintf(output1, "ElemType Linear\n");
			fprintf(output1, "  Number Of Gauss Points: 4\n");
			fprintf(output1, "  Nodes included\n");
			fprintf(output1, "  Natural Coordinates: Internal\n");
			fprintf(output1, "End gausspoints\n\n");
		}

		// RESULTS GROUP deslocamentos e forças
		fprintf(output1, "ResultGroup ");
		fprintf(output1, " \"Analise trelica - nós\"");
		fprintf(output1, " %i ", it);
		fprintf(output1, "OnNodes \n");

		// Impressao dos deslocamentos
		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Deslocamento (m)\"");
		fprintf(output1, " Vector\n");
		fprintf(output1, "                     ComponentNames ");
		fprintf(output1, "\"Desloc X (m)\", \"Desloc Y (m)\", \"Desloc Z (m)\"\n");

		// Impressao das forças Externas
		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Forças Externas (N)\"");
		fprintf(output1, " Vector\n");
		fprintf(output1, "                     ComponentNames ");
		fprintf(output1, "\"Força X (N)\", \"Força Y (N)\", \"Força Z (N)\"\n");

		// Impressao das forças de reação
		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Reações (N)\"");
		fprintf(output1, " Vector\n");
		fprintf(output1, "                     ComponentNames ");
		fprintf(output1, "\"Força X (N)\", \"Força Y (N)\", \"Força Z (N)\"\n");

		//Impressao das forças internas
		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Forças Internas (N)\"");
		fprintf(output1, " Vector\n");
		fprintf(output1, "                     ComponentNames ");
		fprintf(output1, "\"Força X (N)\", \"Força Y (N)\", \"Força Z (N)\"\n");

		fprintf(output1, "Values\n");

		//Caso seja uma estrutura 3D ou 2D
		if (n_dim == 3) {
			for (int i = 1; i <= n_nos; i++) {
				fprintf(output1, "%i\t  %30.15f\t  %30.15f\t  %30.15f\t", i, D[n_dim * i - (n_dim - 1)], D[n_dim * i - (n_dim - 2)], D[n_dim * i]);
				fprintf(output1, "%30.15f\t", Fe[n_dim * i - (n_dim - 1)]);
				fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\t", Fe[n_dim * i - (n_dim - 2)], Fe[n_dim * i], Fd[n_dim * i - (n_dim - 1)], Fd[n_dim * i - (n_dim - 2)], Fd[n_dim * i]);
				fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\n", Fi[n_dim * i - (n_dim - 1)], Fi[n_dim * i - (n_dim - 2)], Fi[n_dim * i]);
			}
			fprintf(output1, "End Values\n");
		}
		else {
			for (int i = 1; i <= n_nos; i++) {
				fprintf(output1, "%i\t  %30.15f\t  %30.15f\t  %30.15f\t", i, D[n_dim * i - (n_dim - 1)], D[n_dim * i - (n_dim - 2)], 0.0);
				fprintf(output1, "%30.15f\t", Fe[n_dim * i - (n_dim - 1)]);
				fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\t", Fe[n_dim * i - (n_dim - 2)], 0.0, Fd[n_dim * i - (n_dim - 1)], Fd[n_dim * i - (n_dim - 2)], 0.0);
				fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\n", Fi[n_dim * i - (n_dim - 1)], Fi[n_dim * i - (n_dim - 2)], 0.0);
			}
			fprintf(output1, "End Values\n");
		}

		fprintf(output1, "\nResultGroup ");
		fprintf(output1, " \"Analise trelica - elementos\"");
		fprintf(output1, " %lu ", it);
		fprintf(output1, "OnGaussPoints \"Gauss Points\" \n");

		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Tamanho das barras (m)\"");
		fprintf(output1, " Scalar\n");

		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Tensão normal (Pa)\"");
		fprintf(output1, " Scalar\n");

		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Deformações (m)\"");
		fprintf(output1, " Scalar\n");

		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Forças nas barras (N)\"");
		fprintf(output1, " Scalar\n");

		fprintf(output1, "Values\n");
		for (int i = 1; i <= n_barras; i++) {

			fprintf(output1, "%i\t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\n", i, lt[i], sigma[i], epsilon[i], sigma[i] * area[i]);
			fprintf(output1, "  \t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\n", lt[i], sigma[i], epsilon[i], sigma[i] * area[i]);
			fprintf(output1, "  \t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\n", lt[i], sigma[i], epsilon[i], sigma[i] * area[i]);
			fprintf(output1, "  \t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\n", lt[i], sigma[i], epsilon[i], sigma[i] * area[i]);
		}

		fprintf(output1, "End Values\n");
	}

	//PARA ANÁLISE DINÂMICA IMPRIME AS TRÊS MENORES FREQUENCIAS VIBRACIONAIS E SEUS RESPECTIVOS MODOS DE VIBRAR

	if (strcmp(analise, "Dinâmica") == 0) {
		
		if(ne == 0) {
			fprintf(output1, "GiD Post Results File 1.0\n\n");
			// GaussPoints

			fprintf(output1, "GaussPoints ");
			fprintf(output1, "\"Gauss Points\" ");
			fprintf(output1, "ElemType Linear\n");
			fprintf(output1, "  Number Of Gauss Points: 4\n");
			fprintf(output1, "  Nodes included\n");
			fprintf(output1, "  Natural Coordinates: Internal\n");
			fprintf(output1, "End gausspoints\n\n");
		}
		// RESULTS GROUP deslocamentos e forças
		fprintf(output1, "ResultGroup ");
		fprintf(output1, " \"Analise trelica - nós\"");
		fprintf(output1, " %i ", ne);
		fprintf(output1, "OnNodes \n");

		// Impressao dos deslocamentos
		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Deslocamento (m)\"");
		fprintf(output1, " Vector\n");
		fprintf(output1, "                     ComponentNames ");
		fprintf(output1, "\"Desloc X (m)\", \"Desloc Y (m)\", \"Desloc Z (m)\"\n");

		// Impressao das forças Externas
		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Forças Externas (N)\"");
		fprintf(output1, " Vector\n");
		fprintf(output1, "                     ComponentNames ");
		fprintf(output1, "\"Força X (N)\", \"Força Y (N)\", \"Força Z (N)\"\n");

		// Impressao das forças de reação
		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Reações (N)\"");
		fprintf(output1, " Vector\n");
		fprintf(output1, "                     ComponentNames ");
		fprintf(output1, "\"Força X (N)\", \"Força Y (N)\", \"Força Z (N)\"\n");

		//Impressao das forças internas
		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Forças Internas (N)\"");
		fprintf(output1, " Vector\n");
		fprintf(output1, "                     ComponentNames ");
		fprintf(output1, "\"Força X (N)\", \"Força Y (N)\", \"Força Z (N)\"\n");

		fprintf(output1, "Values\n");

		//Caso seja uma estrutura 3D ou 2D
		if (n_dim == 3) {
			for (int i = 1; i <= n_nos; i++) {
				fprintf(output1, "%i\t  %30.15f\t  %30.15f\t  %30.15f\t", i, D[n_dim * i - (n_dim - 1)], D[n_dim * i - (n_dim - 2)], D[n_dim * i]);
				fprintf(output1, "%30.15f\t", Fe[n_dim * i - (n_dim - 1)]);
				fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\t", Fe[n_dim * i - (n_dim - 2)], Fe[n_dim * i],-Fd[n_dim * i - (n_dim - 1)], -Fd[n_dim * i - (n_dim - 2)], -Fd[n_dim * i]);
				fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\n", Fi[n_dim * i - (n_dim - 1)], Fi[n_dim * i - (n_dim - 2)], Fi[n_dim * i]);
			}
			fprintf(output1, "End Values\n");
		}
		else {
			for (int i = 1; i <= n_nos; i++) {
				fprintf(output1, "%i\t  %30.15f\t  %30.15f\t  %30.15f\t", i, D[n_dim * i - (n_dim - 1)], D[n_dim * i - (n_dim - 2)], 0.0);
				fprintf(output1, "%30.15f\t", Fe[n_dim * i - (n_dim - 1)]);
				fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\t", Fe[n_dim * i - (n_dim - 2)], 0.0, -Fd[n_dim * i - (n_dim - 1)], -Fd[n_dim * i - (n_dim - 2)], 0.0);
				fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\n", Fi[n_dim * i - (n_dim - 1)], Fi[n_dim * i - (n_dim - 2)], 0.0);
			}
			fprintf(output1, "End Values\n");
		}

		fprintf(output1, "\nResultGroup ");
		fprintf(output1, " \"Analise trelica - elementos\"");
		fprintf(output1, " %lu ", ne);
		fprintf(output1, "OnGaussPoints \"Gauss Points\" \n");

		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Tamanho das barras (m)\"");
		fprintf(output1, " Scalar\n");

		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Tensão normal (Pa)\"");
		fprintf(output1, " Scalar\n");

		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Deformações (m)\"");
		fprintf(output1, " Scalar\n");

		fprintf(output1, "  ResultDescription ");
		fprintf(output1, "\"Forças nas barras (N)\"");
		fprintf(output1, " Scalar\n");

		fprintf(output1, "Values\n");
		for (int i = 1; i <= n_barras; i++) {

			fprintf(output1, "%i\t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\n", i, lt[i], sigma[i], epsilon[i], sigma[i] * area[i]);
			fprintf(output1, "  \t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\n", lt[i], sigma[i], epsilon[i], sigma[i] * area[i]);
			fprintf(output1, "  \t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\n", lt[i], sigma[i], epsilon[i], sigma[i] * area[i]);
			fprintf(output1, "  \t  %30.15f\t  %30.15f\t  %30.15f\t  %30.15f\n", lt[i], sigma[i], epsilon[i], sigma[i] * area[i]);
		}

		fprintf(output1, "End Values\n");
		//A única finalidade em cálcular os modos de vibração é caso haja amortecimento no sistema
		if (amort == 1) {
			// RESULTS GROUP PARA DINÂMICO AUTOVETORES PARA AS 3 MENORES FREQUÊNCIAS DE VIBRAÇÃO
			fprintf(output1, "\nResultGroup ");
			fprintf(output1, " \"Analise Modal - nós\"");
			fprintf(output1, " %i ", ne);
			fprintf(output1, "OnNodes \n");

			// Impressao do modo de deformação associado à primeira menor frequência natural
			fprintf(output1, "  ResultDescription ");
			fprintf(output1, "\"Frequência %6.2f", sqrt(w2[n_dim*n_nos - countw2]) / (2 * PI));
			fprintf(output1, " Hz\"");
			fprintf(output1, " Vector\n");
			fprintf(output1, "                     ComponentNames ");
			fprintf(output1, "\"Modo X\", \"Modo Y\", \"Modo Z\"\n");

			// Impressao do modo de deformação associado à segunda menor frequência natural
			fprintf(output1, "  ResultDescription ");
			fprintf(output1, "\"Frequência %6.2f", sqrt(w2[n_dim*n_nos - countw2 - 1]) / (2 * PI));
			fprintf(output1, " Hz\"");
			fprintf(output1, " Vector\n");
			fprintf(output1, "                     ComponentNames ");
			fprintf(output1, "\"Modo X\", \"Modo Y\", \"Modo Z\"\n");

			// Impressao do modo de deformação associado à terceira menor frequência natural
			fprintf(output1, "  ResultDescription ");
			fprintf(output1, "\"Frequência %6.2f", sqrt(w2[n_dim*n_nos - countw2 - 2]) / (2 * PI));
			fprintf(output1, " Hz\"");
			fprintf(output1, " Vector\n");
			fprintf(output1, "                     ComponentNames ");
			fprintf(output1, "\"Modo X\", \"Modo Y\", \"Modo Z\"\n");

			fprintf(output1, "Values\n");
			//Caso seja uma estrutura 3D ou 2D
			if (n_dim == 3) {
				for (int i = 1; i <= n_nos; i++) {
					fprintf(output1, "%i\t  %30.15f\t  %30.15f\t  %30.15f\t", i, autovector[n_dim * i - (n_dim - 1)][n_dim*n_nos - countw2], autovector[n_dim * i - (n_dim - 2)][n_dim*n_nos - countw2], autovector[n_dim * i][n_dim*n_nos - countw2]);
					fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\t", autovector[n_dim * i - (n_dim - 1)][n_dim*n_nos - countw2 - 1], autovector[n_dim * i - (n_dim - 2)][n_dim*n_nos - countw2 - 1], autovector[n_dim * i][n_dim*n_nos - countw2 - 1]);
					fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\n", autovector[n_dim * i - (n_dim - 1)][n_dim*n_nos - countw2 - 2], autovector[n_dim * i - (n_dim - 2)][n_dim*n_nos - countw2 - 2], autovector[n_dim * i][n_dim*n_nos - countw2 - 2]);
				}
				fprintf(output1, "End Values\n");
			}
			else {
				for (int i = 1; i <= n_nos; i++) {
					fprintf(output1, "%i\t  %30.15f\t  %30.15f\t  %30.15f\t", i, autovector[n_dim * i - (n_dim - 1)][n_dim*n_nos - countw2], autovector[n_dim * i - (n_dim - 2)][n_dim*n_nos - countw2], 0.0);
					fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\t", autovector[n_dim * i - (n_dim - 1)][n_dim*n_nos - countw2 - 1], autovector[n_dim * i - (n_dim - 2)][n_dim*n_nos - countw2 - 1], 0.0);
					fprintf(output1, "%30.15f\t  %30.15f\t  %30.15f\n", autovector[n_dim * i - (n_dim - 1)][n_dim*n_nos - countw2 - 2], autovector[n_dim * i - (n_dim - 2)][n_dim*n_nos - countw2 - 2], 0.0);
				}
				fprintf(output1, "End Values\n");
			}
		}
	}
	fclose(output1);

}
