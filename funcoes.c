#include"funcoes.h"

//////////////////////////////////FUNÇÕES BÁSICAS DO PROGRAMA DE TRELIÇAS///////////////////////////////////
// Protótipos presentes em funcoes.h																	  //
// Funções Básicas de resolução do sistema linear para cada iteração de Newton							  //																	  
// Variáveis presentes em var.h																			  //
// Alocação de memória dinâmica e desalocação presentes em nrutil.h										  //
// Chamada da função newmark() ou newtonestatico() em origem.cpp. Isto varia de acordo com solicitação GiD// 
// A partir delas, chama-se as demais funções referentes à cada análise									  //				 
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*****************************************************************************************************/
/***************************************INTEGRAÇÃO NO TEMPO*******************************************/
void newmark() {
	/*****************************************************************************************************/
	/******************************ALOCAÇÃO DINÂMICA DE VARIÁVEIS UTILIZADAS******************************/
	K = dmatrix_zeros(1, n_dim * n_nos, 1, n_dim * n_nos);
	M = dmatrix_zeros(1, n_dim * n_nos, 1, n_dim * n_nos);
	Kf = dmatrix_zeros(1, n_dim * n_nos, 1, n_dim * n_nos);
	C = dmatrix_zeros(1, n_dim * n_nos, 1, n_dim * n_nos);
	D = dvector_zeros(1, n_dim * n_nos); 
	D1 = dvector_zeros(1, n_dim * n_nos);
	A = dvector_zeros(1, n_dim * n_nos);
	A1 = dvector_zeros(1, n_dim * n_nos);
	V = dvector_zeros(1, n_dim * n_nos);
	V1 = dvector_zeros(1, n_dim * n_nos);
	l0 = dvector(1, n_barras);
	kaux = dmatrix_zeros(1, 2 * n_dim, 1, 2 * n_dim);
	bv = dvector(1, 3);
	Fi = dvector_zeros(1, n_dim * n_nos);
	Fe = dvector_zeros(1, n_dim * n_nos);
	Ff = dvector_zeros(1, n_dim * n_nos);
	Fr = dvector_zeros(1, n_dim * n_nos);
	Fd = dvector_zeros(1, n_dim * n_nos);
	sigma = dvector_zeros(1, n_barras);
	lambda = dvector_zeros(1, n_barras);
	lt = dvector_zeros(1, n_barras);
	epsilon = dvector_zeros(1, n_barras);
	mass = dmatrix_zeros(1, n_dim * 2, 1, n_dim * 2);
	w2 = dvector_zeros(1, n_dim * n_nos);
	autovector = dmatrix_zeros(1, n_dim * n_nos, 1, n_dim * n_nos);
	
	//Vetores auxiliares de coordenadas
	auxCx = dvector_zeros(1, n_nos);
	auxCy = dvector_zeros(1, n_nos);
	auxCz = dvector_zeros(1, n_nos);

	//Copia coordenadas  para referência do lagrangeano total
	equalvectord(Cx, auxCx, 1, n_nos);
	equalvectord(Cy, auxCy, 1, n_nos);
	equalvectord(Cz, auxCz, 1, n_nos);

	/*****************************************************************************************************/
	//matriz de rigidez passo inicial
	r == BEGIN;
	matrizrigidez();
	//A matriz massa [M]
	matrizmassa();
	//As frequências naturais e modos de vibração respectivos (apenas se houver amortecimento)
	if (amort == 1) frequencias();
	//define a matriz amortecimento[C]  (apenas se houver amortecimento)
	if (amort == 1) matrizamortecimento();
	//matriz efetiva
	rigefetiva();

	do
	{	
		//Imprime arquivos resultantes da análise dinâmica
		saidadados();

		//Zera vetor
		empvectord(Fd, 1, n_dim*n_nos);

		//Aloca vetor força externa
		vectorfext();

		//Cálculo das forças efetivas a partir da força interna [K] e de [M] e [C]
		//Caso não haja nenhuma força aplicada na estrutura no passo, não é executado o método de newmark. Pula o passo
		unsigned long int count = 0;
		for (int i = 1; i <= n_dim * n_nos; i++) {
			for (int j = 1; j <= n_dim * n_nos; j++) { Fd[i] += M[i][j] * (C2*V[j] + C3 * A[j]) + C[i][j] * (C4* V[j] + C5 * A[j]); }
			Fd[i] += +Fe[i] - Fi[i];
			if (fabs(Fd[i]) >= Erro*deltaT) count++;
		}
		if (count == 0) { 
			//atualização do intervalo de tempo
			Ti += deltaT; 
			// Atualiza o contador
			ne++;
			continue; 
		}

		//Aplicação da penalidade
		penal(Kf);

		//Cálculo do Dinicial
		strcmp(solve, "PBCG") == 0 ? callpbcg() : callcholesky();

		//Inicie a energia - Condição de parada
		condparada();
		
		// Atualiza o contador
		r++;

		//Newton aplicado ao processo dinâmico
		dinamic();

		//atualização do intervalo de tempo
		Ti += deltaT;
		
		// Atualiza o contador
		ne++;
	} while (Ti <= Tf);

}
/*****************************************************************************************************/
/********************************MÉTODO DE NEWTON (ANÁLISE ESTÁTICA)**********************************/
void newtonestatico() {
	/******************************ALOCAÇÃO DINÂMICA DE VARIÁVEIS UTILIZADAS******************************/
	K = dmatrix_zeros(1, n_dim * n_nos, 1, n_dim * n_nos);
	Fd = dvector_zeros(1, n_dim * n_nos);
	D = dvector_zeros(1, n_dim * n_nos); // Atribui deslocamentos iniciais aos nós 
	D1 = dvector_zeros(1, n_dim * n_nos);
	l0 = dvector(1, n_barras);
	kaux = dmatrix_zeros(1, 2 * n_dim, 1, 2 * n_dim);
	bv = dvector(1, 3);
	Fi = dvector_zeros(1, n_dim * n_nos);
	Fe = dvector_zeros(1, n_dim * n_nos);
	Fr = dvector_zeros(1, n_dim * n_nos);
	sigma = dvector_zeros(1, n_barras);
	lambda = dvector_zeros(1, n_barras);
	lt = dvector_zeros(1, n_barras);
	epsilon = dvector_zeros(1, n_barras);
	//Vetores auxiliares de coordenadas
	auxCx = dvector_zeros(1, n_nos);
	auxCy = dvector_zeros(1, n_nos);
	auxCz = dvector_zeros(1, n_nos);

	//Aloca vetor força externa
	vectorfext();
	/*****************************************************************************************************/
	/**********FAZ AS ITERAÇÕES DE NEWTON PARA ACHAR O ACRESCIMO DE DESLOCAMENTO A CADA UMA DELAS*********/

	//iguala coordenadas para ter uma referência lagrangeana inicial
	equalvectord(Cx, auxCx, 1, n_nos);
	equalvectord(Cy, auxCy, 1, n_nos);
	equalvectord(Cz, auxCz, 1, n_nos);

	//Vetor auxiliar para cálculo incremental
	double *Fei;
	Fei = dvector_zeros(1, n_dim * n_nos);
	equalvectord(Fe, Fei, 1, n_dim*n_nos);

	//Inicia o procedimento d Newton incremental-iterativo. Se NÃO for incremental, este loop só executa 1x
	for (it = 1; it<=inc ; it++) {

		r = BEGIN;
		double aux1;
		do {  //dispara newton. Em r=0 faz o cálculo para a matriz linear

			  //Escolha do tipo de análise (Estática ou dinâmica)
			matrizrigidez();

			//aplica a penalidade referente aos nós restringidos para que possa resolver o sistema

			penal(K);

			//Diferenças de forças (internas-externas) grad1
			//Caso não haja nenhuma força aplicada na estrutura, não é executado o método de newton
			unsigned long int count = 0;
			for (int i = 1; i <= n_dim * n_nos; i++) {
				Fd[i] = Fi[i] - Fe[i];
				if (fabs(Fd[i]) > Erro) count++;
			}
			if (count == 0) break;

			/*Cálculo do deslocamento*///Cálculo da direção inicial do sistema
			strcmp(solve, "PBCG") == 0 ? callpbcg() : callcholesky();

			//aplica decremento (deslocamento inicial menos o deslocamento corrente)
			for (int i = 1; i <= n_dim * n_nos; i++) D[i] -= D1[i];

			//Inicie a energia - Condição de parada
			condparada();
			
			//Caso não haja convergência, escreve arquivo de saída

			if (r > BEGIN && (fabs(Gc) / aux1) > 1.0) {
				//Arquivo de saida de dados para erros
				char fileerr[1024];
				extern char projname[1024];
				FILE *ferr;
				strcpy(fileerr, projname);
				strcat(fileerr, ".err");
				ferr = fopen(fileerr, "w");
				fprintf(ferr, "\n \n \n *****  ERRO: Convergência não atingida ***** \n");
				fprintf(ferr, "\n \n \n Iteração de Newton passo: %i", r);
				fclose(ferr);
				exit(1);
			}
			aux1 = fabs(Gc);

			//atualiza valores de coordenadas
			refreshcord();

			//Verifica condição de equilíbrio: Se Gc<=TOL*GI, para a iteração
			r++;
		} while (fabs(Gc) >= Erro * fabs(GI));

		//escreve arquivo para análise estática
		saidadados();

		//Quebra o laço caso não seja incremental
		if (strcmp(increm, "Não") == 0) break;

		// Iguala para começar novo procedimento incremental
		equalvectord(auxCx, Cx, 1, n_nos);
		equalvectord(auxCy, Cy, 1, n_nos);
		equalvectord(auxCz, Cz, 1, n_nos);

		// Incrementa a carga externa

		for (int k = 1; k <= n_dim * n_nos; k++) Fe[k] += Fei[k];
		// Zera vetor dos deslocamentos

		empvectord(D, 1, n_dim*n_nos);
		
	}

	//desaloca variável local
	free_dvector(Fei, 1, n_dim*n_nos);
}
/*****************************************************************************************************/
/*****************************************FORÇA RESIDUAL**********************************************/
void forcaresid() {
	//Zera vetor
	empvectord(Fd, 1, n_dim*n_nos);

	//Cálculo de Fd (força residual) a partir da força interna [K] e de [M] e [C]
	for (int i = 1; i <= n_dim*n_nos; i++) {
		for (int j = 1; j <= n_dim*n_nos; j++) {
			Fd[i] += -M[i][j] * A[j] - C[i][j] * V[j];
		}
		Fd[i] += -Fi[i] + Fe[i];
	}

}
/*****************************************************************************************************/
/*************************************RIGIDEZ EFETIVA DINÂMICA****************************************/
void rigefetiva() {
	//Cálculo de Mf (matriz efetiva de rigidez) a partir da matriz de rigidez([K]) e de [M] e [C]
	for (int i = 1; i <= n_dim*n_nos; i++) {
		for (int j = 1; j <= n_dim*n_nos; j++) {
			Kf[i][j] = K[i][j] + C0 * M[i][j] + C1 * C[i][j];
		}
	}

}
/*****************************************************************************************************/
/************************************MATRIZ DE MASSA DO SISTEMA***************************************/
void matrizmassa() {
	empmatrixd(M, 1, n_dim*n_nos, 1, n_dim*n_nos);
	/*Cálculo da matriz massa diferenciando para casos em 2D e 3D, já que o processo de triangulização feito por cholesky
	não funciona com matrizes 3D completas, já que o determinante é igual a zero nos casos em que a estrutura é 2d*/
	
	//calcula para cada elemento a matriz massa em coordenadas globais
	mmass();

	//aloca matriz massa global 
	for (int p = 1; p <= n_barras; p++) {

		//Aloca na matriz massa global do sistema
		for (int i = n_dim*pos_i[p] - (n_dim - 1), i0 = n_dim*pos_j[p] - (n_dim - 1), k = 1;
			i <= n_dim*pos_i[p], i0 <= n_dim*pos_j[p], k <= n_dim; i++, i0++, k++) {
			for (int j = n_dim*pos_i[p] - (n_dim - 1), j0 = n_dim*pos_j[p] - (n_dim - 1), k0 = 1;
				j <= n_dim*pos_i[p], j0 <= n_dim*pos_j[p], k0 <= n_dim; j++, j0++, k0++) {
				M[i][j] += ((rho[p] * area[p] * lt[p]) / 6)*mass[k][k0];
				M[i][j0] += ((rho[p] * area[p] * lt[p]) / 6)*mass[k][k0 + n_dim];
				M[i0][j] += ((rho[p] * area[p] * lt[p]) / 6)*mass[k + n_dim][k0];
				M[i0][j0] += ((rho[p] * area[p] * lt[p]) / 6)*mass[k + n_dim][k0 + n_dim];

			}
		}
	}

}
/*****************************************************************************************************/
/********************************MATRIZ DE AMORTECIMENTO DO SISTEMA***********************************/
void matrizamortecimento() {

	/*Monta-se a matriz de amortecimento para cada nova matriz de rigidez [K] do processo de newton.
	Conserva-se a matriz massa [M] e os parãmetros como alfa e beta, além das frequências e modos de
	vibração (são calculados para sistema em repouso, apenas com matriz massa e rigidez linear).
	*/
	for (int i = 1; i <= n_dim * n_nos; i++) {
		for (int j = 1; j <= n_dim * n_nos; j++) {
			C[i][j] = alfa*M[i][j] + beta*K[i][j];
		}
	}

}
/*****************************************************************************************************/
/****************************CALCULO DAS FREQUENCIAS E MODOS DE VIBRAR********************************/
void frequencias() {

	//Variáveis
	double **mdinamic, *diag, sum, **sum1, **mmassglob2;
	int *nrot;

	//Conta a quantidade de frequências naturais da estrutura considerando suas restrições
	countw2 = 0;
	for (int i = 1; i <= n_nos; i++) {
		if (R_x[i] == 1) countw2++;
		if (R_y[i] == 1) countw2++;
		if (R_z[i] == 1) countw2++;
	}

	//Alocação dinâmica

	diag = dvector(1, n_dim * n_nos);
	mdinamic = dmatrix_zeros(1, n_dim * n_nos, 1, n_dim * n_nos);
	sum1 = dmatrix_zeros(1, n_dim* n_nos, 1, n_dim * n_nos);
	nrot = ivector(1, 1);
	mmassglob2 = dmatrix(1, n_dim * n_nos, 1, n_dim * n_nos);

	//iguala matrizes para prosseguir cálculo
	equalmatrixd(M, mmassglob2, 1, n_dim * n_nos, 1, n_dim * n_nos);

	//Calcula a matriz triangular inferior da matriz massa [M]
	choldc(mmassglob2, n_dim * n_nos, diag);

	/*A typical use of choldc and cholslis in the inversion of covariancematrices describing the fit
	of data to a model, see, e.g., x15.6. In this, and many other applications, oe often needs L^-1. The
	Lower triangle of this matriz can be efficiently found from the output of choldc (PRESS,199) fsafasd*/ 

	for (int i = 1; i <= n_dim * n_nos; i++) {
		mmassglob2[i][i] = 1.0 / diag[i];
		for (int j = i + 1; j <= n_dim * n_nos; j++) {
			sum = 0.0;
			for (int k = i; k<j; k++) sum -= mmassglob2[j][k] * mmassglob2[k][i];
			mmassglob2[j][i] = sum / diag[j];
		}
	}

	//Calcula a matriz dinâmica [D] para posterior cálculo dos autovalores e autovetores [D]=[L]^-1 * [K] * [L]^-1^t
	for (int i = 1; i <= n_dim * n_nos; i++) {
		for (int j = 1; j <= n_dim * n_nos; j++) {
			for (int k = 1; k <= i; k++) mdinamic[i][j] += mmassglob2[i][k] * K[k][j];
		}
		for (int j = 1; j <= n_dim * n_nos; j++) {
			for (int k = 1; k <= j; k++) sum1[i][j] += mdinamic[i][k] * mmassglob2[j][k];
		}
	}

	//iguala matrizes para achar a matriz dinâmica [D]
	equalmatrixd(sum1, mdinamic, 1, n_dim * n_nos, 1, n_dim * n_nos);

	/*zera as linhas e colunas onde há restrições na matriz dinâmica para que os
	autovalores calculados sejam considerando estes amarramentos*/
	for (int i = 1; i <= n_nos; i++) {
		if (R_x[i] == 1) for (int j = 1; j <= n_dim*n_nos; j++) {
			mdinamic[n_dim * i - (n_dim - 1)][j] = 0.0;
			mdinamic[j][n_dim * i - (n_dim - 1)] = 0.0;
		}
		if (R_y[i] == 1) for (int j = 1; j <= n_dim * n_nos; j++) {
			mdinamic[n_dim * i - (n_dim - 2)][j] = 0.0;
			mdinamic[j][n_dim * i - (n_dim - 2)] = 0.0;
		}
		if (R_z[i] == 1) for (int j = 1; j <= n_dim * n_nos; j++) {
			mdinamic[n_dim * i][j] = 0.0;
			mdinamic[j][n_dim * i] = 0.0;
		}
	}

	//------------------------------------------------------------------------------------------//
	//calcula e ordena decrescentemente os autovalores e autovetores por Jacobi
	jacobi(mdinamic, n_dim * n_nos, w2, autovector, nrot);
	//ordena de forma decrescente os autovalores e autovetores
	eigsrt(w2, autovector, n_dim * n_nos);
	//------------------------------------------------------------------------------------------//

	//Seleciona as duas menores frequências para cálculo de alfa e beta (utilizáveis para matriz de amortecimento de Rayleigh)
	//1) beta:
	beta = 2 * (Csi / 100) / (sqrt(w2[n_dim * n_nos - countw2]) + sqrt(w2[n_dim * n_nos - countw2 - 1]));
		
	//2) alfa:
	alfa = 2 * (Csi / 100)*sqrt(w2[n_dim * n_nos - countw2])*sqrt(w2[n_dim * n_nos - countw2-1])/
		(sqrt(w2[n_dim * n_nos - countw2]) + sqrt(w2[n_dim * n_nos - countw2 - 1]));

	//Libera memória
	free_dmatrix(mmassglob2, 1, n_dim * n_nos, 1, n_dim * n_nos);
	free_dmatrix(mdinamic, 1, n_dim * n_nos, 1, n_dim * n_nos);
	free_dmatrix(sum1, 1, n_dim * n_nos, 1, n_dim * n_nos);
	free_dvector(diag, 1, n_dim * n_nos);
	free_ivector(nrot, 1, 1);
}
/*****************************************************************************************************/
/**************************************NEWTON ANÁLISE DINÂMICA****************************************/
void dinamic() {
	//vetores auxiliares
	double *D2, *D3;
	D2 = dvector_zeros(1, n_dim*n_nos);
	D3 = dvector_zeros(1, n_dim*n_nos);
	//Iguala velocidades, acelerações e deslocamentos com passoa anterior de tempo
	equalvectord(V, V1, 1, n_dim*n_nos); equalvectord(A, A1, 1, n_dim*n_nos); equalvectord(D, D3, 1, n_dim*n_nos);
	
	do
	{
		//Avalia as novas velocidades, acelerações e deslocamentos
		for (int i = 1; i <= n_dim * n_nos; i++) {
			D[i] = D1[i]+D3[i];
			V[i] = C1 * D1[i] - C4 * V1[i] - C5 * A1[i];
			A[i] = C0 * D1[i] - C2 * V1[i] - C3 * A1[i];
		}
		//Atualiza as coordenadas para cada newton (atualiza matriz de rigidez e vetor de forças internas)
		refreshcord();

		//matriz de rigidez
		matrizrigidez();

		//A matriz massa [M]
		matrizmassa();

		//As frequências naturais e modos de vibração respectivos para incremento (apenas se houver amortecimento)
		if (amort == 1) frequencias();

		//define a matriz amortecimento[C] para incremento
		if (amort == 1) matrizamortecimento();

		//Rigidez efetiva 
		rigefetiva();

		//Força residual
		forcaresid();

		//Aplicação da penalidade
		penal(Kf);

		//Cálculo do Dinicial+1
		equalvectord(D1, D2, 1, n_dim*n_nos);
		strcmp(solve, "PBCG") == 0 ? callpbcg() : callcholesky();

		//atualiza deslocamento
		for (int i = 1; i <= n_dim * n_nos; i++) D1[i] += D2[i];
		//Condição de parada
		condparada();
	} while (fabs(Gc) >= Erro * fabs(GI));

	//Avalia as novas velocidades, acelerações e deslocamentos
	for (int i = 1; i <= n_dim * n_nos; i++) {
		D[i] = D1[i] + D3[i];
		V[i] = C1 * D1[i] - C4 * V1[i] - C5 * A1[i];
		A[i] = C0 * D1[i] - C2 * V1[i] - C3 * A1[i];
	}
	
	//desaloca memoria dinamica
	free_dvector(D2, 1, n_dim*n_nos);
	free_dvector(D3, 1, n_dim*n_nos);
	
}
/*****************************************************************************************************/
/**************************************MATRIZ DE MASSA GLOBAL*****************************************/
void mmass() {

	//Matriz massa triangular para uma barra em 3D
	if (n_dim == 3) {
		mass[1][1] = 2;
		mass[2][1] = 0; mass[2][2] = 2;
		mass[3][1] = 0; mass[3][2] = 0; mass[3][3] = 2;
		mass[4][1] = 1; mass[4][2] = 0; mass[4][3] = 0; mass[4][4] = 2;
		mass[5][1] = 0; mass[5][2] = 1; mass[5][3] = 0; mass[5][4] = 0; mass[5][5] = 2;
		mass[6][1] = 0; mass[6][2] = 0; mass[6][3] = 1; mass[6][4] = 0; mass[6][5] = 0; mass[6][6] = 2;

	}
	//Matriz massa triangular para uma barra em 2D
	else {
		mass[1][1] = 2;
		mass[2][1] = 0; mass[2][2] = 2;
		mass[3][1] = 1; mass[3][2] = 0;  mass[3][3] = 2;
		mass[4][1] = 0; mass[4][2] = 1; mass[4][3] = 0; mass[4][4] = 2;

	}
	//Completa os outros elementos restantes
	for (int i = 2 * n_dim; i >= 1; i--) {
		for (int j = 1; j <= i; j++) {
			mass[j][i] = mass[i][j];
		}
	}

}
/*****************************************************************************************************/
/*************************MATRIZ DE RIGIDEZ DO SISTEMA(PROCESSOS INICIAIS)****************************/
void matrizrigidez() {

	//declaração de variáveis locais
	double Eac1;
	double const_geom1, const_geom2, const_constit;

	//Zera valores para cada newton
	empvectord(Fi, 1, n_dim * n_nos);
	empmatrixd(K, 1, n_dim * n_nos, 1, n_dim * n_nos);


	/********************FAZ A MATRIZ DE RIGIDEZ LOCAL DE CADA BARRA E ALOCA NA GLOBAL********************/

	for (je = 1; je <= n_barras; je++) {

		//calc módulo do vetor b
		bv[1] = Cx[pos_j[je]] - Cx[pos_i[je]];
		bv[2] = Cy[pos_j[je]] - Cy[pos_i[je]];
		bv[3] = Cz[pos_j[je]] - Cz[pos_i[je]];
		//Cálculo do módulo
		lt[je] = sqrt(pow(bv[1], 2) + pow(bv[2], 2) + pow(bv[3], 2));
		//se a iteração de newton for a primeira, l0=lt;(lambda =1)
		if (r == BEGIN) {
			l0[je] = lt[je];
		}
		// Cálculo do estiramento (alongamento linear)
		lambda[je] = lt[je] / l0[je];

		//para calculo da deformação de engenharia
		epsilon[je] = lambda[je] - 1;
		
		//condições para calculo da tensão no aço (ou material a definir)
		//MINORADOR é o coeficiente de minoração do aço, normalmente 1.1
		if (epsilon[je] >= -(fy[je] / MINORADOR) / Eac[je] && epsilon[je] <= 0) {
			Eac1 = Eac[je];
			sigma[je] = epsilon[je] * Eac1;
		}
		else if (epsilon[je] <= (fy[je] / MINORADOR) / Eac[je] && epsilon[je] >= 0) {
			Eac1 = Eac[je];
			sigma[je] = epsilon[je] * Eac1;
		}
		else if (epsilon[je] > (fy[je] / MINORADOR) / Eac[je]) {
			Eac1 = 0;
			sigma[je] = fy[je] / MINORADOR;
		}
		else {
			Eac1 = 0;
			sigma[je] = -fy[je] / MINORADOR;
		}

		//calculo da constante da matriz geom e da constante da matriz constitutiva
		const_geom1 = area[je] * l0[je] * sigma[je] * lambda[je] * (-pow(lt[je], -4));
		const_geom2 = area[je] * l0[je] * sigma[je] * lambda[je] * pow(lt[je], -2);
		const_constit = area[je] * l0[je] * Eac1*pow(lambda[je], 2)*(pow(lt[je], -4));

		//Cálculo da matriz de rigidez geométrica
		matriz_geom(const_geom1, const_geom2);
		//acrescenta matriz de rigidez geométrica na matriz de rigidez global
		matririglob();
		//Cálculo da matriz de rigidez constitutiva
		matriz_const(const_constit);
		//acrescenta matriz de rigidez constitutiva na matriz de rigidez global
		matririglob();
		//cálculo da força interna
		vectorfint(const_geom2);

	}

}
/*****************************************************************************************************/
/***********************************CÁLCULO DA ENERGIA DO SISTEMA*************************************/
void condparada() {
	int i;
	double aux = 0;

	//Trabalho das forças internas e externas p/ cada iteração de newton
	for (i = 1; i <= n_dim * n_nos; i++) {
		aux += D1[i] * (Fd[i]);
	}
	Gc = aux;
	if (r == BEGIN) GI = Gc;


}
/*****************************************************************************************************/
/*************************************VETOR DAS FORÇAS EXTERNAS***************************************/
void vectorfext() {

	//Zera para cada iteração de newmark
	empvectord(Fe, 1, n_dim*n_nos);

	for (int i = 1; i <= n_nos; i++) {

		//Para toda análise estática, as cargas hão de ser constantes ao longo do tempo. 
		
		//Para cargas constantes ao longo do tempo
		if (ft[i] == 1) {
			if (fabs(F_x[i]) > 0) Fe[n_dim * i - (n_dim - 1)] = ko[i] * F_x[i];
			if (fabs(F_y[i]) > 0) Fe[n_dim * i - (n_dim - 2)] = ko[i] * F_y[i];
			if (n_dim == 3) {
				if (fabs(F_z[i]) > 0) Fe[n_dim * i] = ko[i] * F_z[i];
			}
		}
		//para cargas senoidais ao longo do tempo
		else if (ft[i] == 2) {
			if (fabs(F_x[i]) > 0) Fe[n_dim * i - (n_dim - 1)] = F_x[i] * sin(ko[i] * Ti);
			if (fabs(F_y[i]) > 0) Fe[n_dim * i - (n_dim - 2)] = F_y[i] * sin(ko[i] * Ti);
			if (n_dim == 3) { 
				if (fabs(F_z[i]) > 0)  Fe[n_dim * i] = F_z[i] * sin(ko[i] * Ti);
			}

		}
		//para cargas lineares (função 1 grau) no tempo
		else if (ft[i] == 3) {
			if (fabs(F_x[i]) > 0) Fe[n_dim * i - (n_dim - 1)] = F_x[i] * (ko[i] * Ti);
			if (fabs(F_y[i]) > 0) Fe[n_dim * i - (n_dim - 2)] = F_y[i] * (ko[i] * Ti);
			if (n_dim == 3) {
				if (fabs(F_z[i]) > 0) Fe[n_dim * i] = F_z[i] * (ko[i] * Ti);
			}

		}
		//para cargas cossenoidais no tempo
		else if (ft[i] == 4) {
			if (fabs(F_x[i]) > 0) Fe[n_dim * i - (n_dim - 1)] = F_x[i] * cos(ko[i] * Ti);
			if (fabs(F_y[i]) > 0) Fe[n_dim * i - (n_dim - 2)] = F_y[i] * cos(ko[i] * Ti);
			if (n_dim == 3) {
				if (fabs(F_z[i]) > 0) Fe[n_dim * i] = F_z[i] * cos(ko[i] * Ti);
			}

		}
		//para cargas quadráticas (função 2 grau) no tempo
		else if (ft[i] == 5) {
			if (fabs(F_x[i]) > 0)Fe[n_dim * i - (n_dim - 1)] = F_x[i] * ko[i] * pow(Ti, 2);
			if (fabs(F_y[i]) > 0)Fe[n_dim * i - (n_dim - 2)] = F_y[i] * ko[i] * pow(Ti, 2);
			if (n_dim == 3) {
				if (fabs(F_z[i]) > 0) Fe[n_dim * i] = F_z[i] * ko[i] * pow(Ti, 2);
			}
		}
	}

}
/*****************************************************************************************************/
/*************************************MATRIZ DE RIGIDEZ GLOBAL****************************************/
void matririglob() {
	//aloca matriz auxiliar (geométrica e constitutiva) na matriz de rigidez global
	for (int i = (n_dim * pos_i[je] - (n_dim - 1)), j = (n_dim * pos_j[je] - (n_dim - 1)), k = 1;
		i <= n_dim * pos_i[je], j <= n_dim * pos_j[je], k <= n_dim; i++, j++, k++) {
		for (int i0 = (n_dim * pos_i[je] - (n_dim - 1)), j0 = (n_dim * pos_j[je] - (n_dim - 1)), k0 = 1;
			i0 <= n_dim * pos_i[je], j0 <= n_dim * pos_j[je], k0 <= n_dim; i0++, j0++, k0++) {
			K[i][i0] += kaux[k][k0];
			K[i][j0] += kaux[k][k0 + n_dim];
			K[j][j0] += kaux[k + n_dim][k0 + n_dim];
			K[j][i0] += kaux[k + n_dim][k0];

		}
	}


}
/*****************************************************************************************************/
/********************************ATUALIZAÇÃO DAS COORDENADAS NODAIS***********************************/
void refreshcord() {
	//atualiza as coordenadas dos nós (sistema vetorial)
	for (int i = 1; i <= n_nos; i++) {//ok

		Cx[i] = auxCx[i] + D[n_dim * i - (n_dim - 1)];
		Cy[i] = auxCy[i] + D[n_dim * i - (n_dim - 2)];
		if (n_dim == 3) Cz[i] = auxCz[i] + D[n_dim * i];

	}
}
/*****************************************************************************************************/
/*****************************PENALIDADE APLICADA À MATRIZ [K] ou [Kf]********************************/
void penal(double**matrix) {

	//alocar penalidades apenas na matriz de rigidez pra nós restringidos
	for (int m0 = 1; m0 <= n_nos; m0++) {

		if (R_x[m0] == 1) matrix[n_dim * m0 - (n_dim - 1)][n_dim * m0 - (n_dim - 1)] += penalty;
		if (R_y[m0] == 1) matrix[n_dim * m0 - (n_dim - 2)][n_dim * m0 - (n_dim - 2)] += penalty;
		if (R_z[m0] == 1) {
			matrix[n_dim * m0][n_dim * m0] += penalty;
		}
		else continue;
	}

}
/*****************************************************************************************************/
/*********************************FUNÇÃO DE CHAMADA DO MÉTODO PBCG************************************/
void callpbcg() {
	//variáveis e vetores utilizados
	unsigned int aux = 0;
	int *iter;
	double *err;

	// Contagem de elementos não nulos fora da matriz principal para armazenagem em ija e sa
	for (int m = 1; m <= n_dim * n_nos; m++) {
		for (int n = 1; n <= n_dim * n_nos; n++) {
			if (m != n && fabs(strcmp(analise, "Estática") == 0 ? K[m][n] : Kf[m][n]) >= Erro) {
				aux++;
			}
		}
	}
	// Alocação de memória para cálculo
	sa = dvector(1, aux + n_dim * n_nos + 1);
	ija = lvector(1, aux + n_dim * n_nos + 1);
	iter = ivector(1, 1);
	err = dvector(1, 1);

	// Partição da matriz esparsa em sa e ija
	sprsin(strcmp(analise, "Estática") == 0 ? K : Kf, n_dim * n_nos, Erro, aux + n_dim * n_nos + 1, sa,
		ija);

	//Aplicação do PBCG
	linbcg(n_dim * n_nos, Fd, D1, Conv, Erro, n_iter, iter, err);

	//Liberação de memória para nova contagem
	free_dvector(sa, 1, aux + n_dim * n_nos + 1);
	free_lvector(ija, 1, aux + n_dim * n_nos + 1);
	free_dvector(err, 1, 1);
	free_ivector(iter, 1, 1);
}
/*****************************************************************************************************/
/*******************************FUNÇÃO DE CHAMADA DO MÉTODO CHOLESKY**********************************/
void callcholesky() {
	//variável local
	double *vectoraux1;
	//alocação dinâmica de vetor auxiliar
	vectoraux1 = dvector_zeros(1, n_dim * n_nos);

	//Método de Cholesky
	choldc(strcmp(analise, "Estática") == 0 ? K : Kf, n_dim * n_nos, vectoraux1);

	cholsl(strcmp(analise, "Estática") == 0 ? K : Kf, n_dim * n_nos, vectoraux1, Fd, D1);


	//liberação de memória
	free_dvector(vectoraux1, 1, n_dim * n_nos);

}
/*****************************************************************************************************/
/*************************************VETOR DAS FORÇAS INTERNAS***************************************/
void vectorfint(double const_geom2) {

	//auxiliares
	double aux = 0, *vectoraux, *vectoraux2;

	vectoraux = dvector_zeros(1, 2 * n_dim);
	vectoraux2 = dvector_zeros(1, 2 * n_dim);

	//calcula psi transposto * psi e acrescenta na matriz auxiliar
	for (int m0 = 1; m0 <= 2 * n_dim; m0++) {//ok
		for (int m1 = 1; m1 <= 2 * n_dim; m1++) {
			for (int m2 = 1; m2 <= n_dim; m2++) {
				n_dim == 3 ? aux += psi[m2][m0] * psi[m2][m1] : aux += psi1[m2][m0] * psi1[m2][m1];

			}
			kaux[m0][m1] = 0;//
			kaux[m0][m1] += aux;
			aux = 0;
		}
	}
	//Aloca as coordenadas nodais iniciais e finais em um vetor auxiliar
	if (n_dim == 3) {
		vectoraux2[1] = Cx[pos_i[je]];
		vectoraux2[2] = Cy[pos_i[je]];
		vectoraux2[3] = Cz[pos_i[je]];
		vectoraux2[4] = Cx[pos_j[je]];
		vectoraux2[5] = Cy[pos_j[je]];
		vectoraux2[6] = Cz[pos_j[je]];
	}
	else {
		vectoraux2[1] = Cx[pos_i[je]];
		vectoraux2[2] = Cy[pos_i[je]];
		vectoraux2[3] = Cx[pos_j[je]];
		vectoraux2[4] = Cy[pos_j[je]];
	}

	//Calcula força interna para cada barra (localmente)
	for (int m0 = 1; m0 <= 2 * n_dim; m0++) {//ok
		for (int m1 = 1; m1 <= 2 * n_dim; m1++) {

			vectoraux[m0] += kaux[m0][m1] * vectoraux2[m1] * const_geom2;
		}
	}

	//Calculo da força interna (globalmente) //ok
	for (int m0 = n_dim * pos_i[je] - (n_dim - 1), m1 = n_dim * pos_j[je] - (n_dim - 1), m2 = 1;
		m0 <= n_dim * pos_i[je], m1 <= n_dim * pos_i[je], m2 <= n_dim; m0++, m1++, m2++) {

		Fi[m0] += vectoraux[m2];
		Fi[m1] += vectoraux[m2 + n_dim];

	}

	free_dvector(vectoraux, 1, 2 * n_dim);
	free_dvector(vectoraux2, 1, 2 * n_dim);
}
/*****************************************************************************************************/
/***********************************MATRIZ DE RIGIDEZ GEOMÉTRICA**************************************/
void matriz_geom(double const_geom1, double const_geom2) {
	double  *vector1, aux = 0;

	//alocação de memória
	vector1 = dvector_zeros(1, 2 * n_dim);

	//vector1 é psi transposto * vectorb
	for (int m0 = 1; m0 <= 2 * n_dim; m0++) { //ok
		for (int m1 = 1; m1 <= n_dim; m1++) {
			n_dim == 3 ? vector1[m0] += psi[m1][m0] * bv[m1] : vector1[m0] += psi1[m1][m0] * bv[m1];
		}
	}
	//matril é vector1*vector1 transposto
	for (int m0 = 1; m0 <= 2 * n_dim; m0++) { //ok
		for (int m1 = 1; m1 <= 2 * n_dim; m1++) {
			kaux[m0][m1] = vector1[m0] * vector1[m1] * const_geom1;
		}
	}

	//calcula psi transposto * psi, multiplica pela const_geom2 e acrescenta na matriz de rigidez geometrica
	for (int m0 = 1; m0 <= 2 * n_dim; m0++) { //ok
		for (int m1 = 1; m1 <= 2 * n_dim; m1++) {
			for (int m2 = 1; m2 <= n_dim; m2++) {
				n_dim == 3 ? aux += psi[m2][m0] * psi[m2][m1] * const_geom2 : aux += psi1[m2][m0] * psi[m2][m1] * const_geom2;

			}

			kaux[m0][m1] += aux;
			aux = 0;
		}

	}

	free_dvector(vector1, 1, 2 * n_dim);//limpa vetor auxiliar

}
/*****************************************************************************************************/
/**********************************MATRIZ DE RIGIDEZ CONSTITUTIVA*************************************/
void matriz_const(double const_constit) {
	double  *vector1;
	vector1 = dvector_zeros(1, 2 * n_dim);


	//vector1 é psi transposto * vectorb
	for (int m0 = 1; m0 <= 2 * n_dim; m0++) {//ok
		for (int m1 = 1; m1 <= n_dim; m1++) {
			n_dim == 3 ? vector1[m0] += psi[m1][m0] * bv[m1] : vector1[m0] += psi1[m1][m0] * bv[m1];
		}
	}


	//matril é vector1*vector1 transposto
	for (int m0 = 1; m0 <= 2 * n_dim; m0++) {//ok
		for (int m1 = 1; m1 <= 2 * n_dim; m1++) {
			kaux[m0][m1] = vector1[m0] * vector1[m1] * const_constit;

		}
	}

	free_dvector(vector1, 1, 2 * n_dim);//limpa vetor auxiliar

}

/*****************************************************************************************************/

/*************************************OUTRAS FUNÇÕES(SECUNDÁRIAS)*************************************/
// Funções responsáveis por igualar matrizes e vetores, além de zerá-los em condição necessária
// Funções úteis para economizar espaço no código principal acim
//
/*****************************************************************************************************/
/****************************************ZERA VETOREs DOUBLE******************************************/
void empvectord(double *vector, long n1, long nh) {
	/*n1 posição incial do vetor, nh posição final*/

	for (int i = n1; i <= nh; i++) vector[i] = 0;

}
/*****************************************************************************************************/
/***************************************ZERA MATRIZ DOUBLE********************************************/
void empmatrixd(double **matrix, long nrl, long nrh, long ncl, long nch) {
	/*zera matriz double de tamanho indeterminado*/

	for (int i = nrl; i <= nrh; i++) {
		for (int j = ncl; j <= nch; j++) {
			matrix[i][j] = 0;
		}
	}
}
/*****************************************************************************************************/
/***********************************IGUALA DUAS MATRIZES DOUBLE***************************************/
void equalmatrixd(double**matrix1, double ** matrix2, long nrl, long nrh, long ncl, long nch) {
	/*Atribuir valores double de uma matriz 1 para uma matriz 2*/

	for (int i = nrl; i <= nrh; i++) {
		for (int j = ncl; j <= nch; j++) {

			matrix2[i][j] = matrix1[i][j];
		}
	}

}
/*****************************************************************************************************/
/************************************IGUALA DOIS VETORES DOUBLE***************************************/
void equalvectord(double*vectorx1, double *vectorx2, long nrl, long nrh) {
	/*Atribuir valores double de um vetor 1 para uma vetor 2*/

	for (int i = nrl; i <= nrh; i++) vectorx2[i] = vectorx1[i];

}
