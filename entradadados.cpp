#include"funcoes.h"


void jumpline(FILE*); //pula linhas nos arquivos de leitura

void jumpline(FILE*jump) {

	char buffer[1024];

	fgets(buffer, 1023, jump);

}

void leitura(void) {

	int aux;
	char aux2[12];

	//criação de arquivos para leitura do .dat e criação de arquivo .err
	char filedat[1024], fileerr[1024];
	extern char projname[1024];
	FILE *fp, *ferr;

	strcpy(filedat, projname); //copia nome inicial para arquivo de leitura (filedat)
	strcat(filedat, ".dat"); //atribui uma extensão para o arquivo anteriormente copiado

	//abertura do arquivo para leitura
	fp = fopen(filedat, "r");

	//Leitura da quantidade de barras, nós e materiais
	for (int i = 0; i < 4; i++) jumpline(fp);

	fscanf(fp, "%i %i %i %i", &n_barras, &n_nos, &n_mat, &n_dim); 
	if (n_mat < 1) {
		strcpy(fileerr, projname);
		strcat(fileerr, ".err");
		ferr = fopen(fileerr, "w");
		fprintf(ferr, "\n \n \n *****  ERRO: Nenhum material alocado ***** \n");
		fprintf(ferr, "(Tente fazê-lo)\n");
		exit(1);
	}

	//Atribuição de variáveis (alocação dinâmica)
	Cx = dvector(1, n_nos);
	Cy = dvector(1, n_nos);
	Cz = dvector(1, n_nos);
	F_x = dvector(1, n_nos);
	F_y = dvector(1, n_nos);
	F_z = dvector(1, n_nos);
	R_x = ivector(1, n_nos);
	R_y = ivector(1, n_nos);
	R_z = ivector(1, n_nos);
	ft = ivector(1, n_nos);
	ko = dvector(1, n_nos);
	area = dvector(1, n_barras);
	Eac = dvector(1, n_barras);
	fy = dvector(1, n_barras);
	material = ivector(1, n_barras);
	pos_i = ivector(1, n_barras);
	pos_j = ivector(1, n_barras);
	rho = dvector(1, n_barras);
	
	//Penalidade inicial, erros, parâmetros de convergencias do PBGC definidos
	for (int i = 0; i < 4; i++) jumpline(fp);
	fscanf(fp, "%s %le %le %s %s %u %s %i %i %lf %lf %lf %lf %i", &analise, &penalty,&Erro, &deformation,&increm,&inc, &solve, &Conv,
		&n_iter, &Ti, &Tf, &deltaT, &Csi, &amort);

	//Coordenadas e forças nos nós
	for (int i = 0; i < 4; i++) jumpline(fp);
	for (int i = 1; i <= n_nos; i++) {
		fscanf(fp, "%i %lf %lf %lf %lf %lf %lf %s %lf", &aux, &Cx[i], &Cy[i], &Cz[i], &F_x[i], &F_y[i], &F_z[i], &aux2, &ko[i]);
		//Seleção de cada tipo de carregamento utilizado na análise (Constante, Senoidal, Linear,cossenoidal,decrescente e segunda ordem)
		if (strcmp(aux2, "k.A") == 0) ft[i] = 1;
		else if (strcmp(aux2, "A.sin(k.T)") == 0) ft[i] = 2;
		else if (strcmp(aux2, "A.(k.T)") == 0) ft[i] = 3;
		else if (strcmp(aux2, "A.cos(k.T)") == 0) ft[i] = 4;
		else ft[i] = 5; //A.(k.T²)

	}
	
	//Restrições nodais
	for (int i = 0; i < 4; i++) jumpline(fp);
	for (int i = 1; i <= n_nos; i++) fscanf(fp, "%i %i %i %i", &aux, &R_x[i], &R_y[i], &R_z[i]);

	//Elementos e suas propriedades
	for (int i = 0; i < 4; i++) jumpline(fp);
	for (int i = 1; i <= n_barras; i++) fscanf(fp, "%i %i %i %i %lf %le %lf %lf", &aux, &pos_i[i], &pos_j[i], &material[i], 
		&area[i], &Eac[i], &fy[i], &rho[i]);

	fclose(fp);

	}