#include"funcoes.h"


//arquivo recebido inicialmente com nome qualquer que seja dado pelo save do GiD
char projname[1024]; 

int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("Erro: Nome do projeto não fornecido corretamente.\n");
		return 1;
	}
	strcpy(projname, argv[1]);

	//FAZ A LEITURA DE DADOS DO ARQUIVO GERADO PELO GiD
	leitura();
	//CHAMA OS CÁLCULOS NA ESTRUTURA
	strcmp(analise, "Dinâmica") == 0 ? newmark() : newtonestatico();
	//LIMPA AS VARIÁVEIS DINÂMICAS USADAS NO FINAL DO PROCESSO
	free_var();
	return 0;
}

