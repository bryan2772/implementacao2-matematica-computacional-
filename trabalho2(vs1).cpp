/*Implemente o Método Algébrico de Gauss-Jordan e o Método Iterativo de Gauss-Seidel para a resolução de
Sistemas de Equações Lineares com ordem n ≤ 15 | n = m. Ambos os métodos deverão tratar a ocorrência de
coeficientes-pivôs nulos ou pequenos por meio do Método de Pivotamento Parcial.

Método de Gauss-Jordan:
O programa deve ler arquivos de texto com os dados de entrada digitados da seguinte maneira:
𝑛
𝑎11 ␣ 𝑎12 ␣ 𝑎13 ␣ … ␣ 𝑎1𝑛
𝑎21 ␣ 𝑎22 ␣ 𝑎23 ␣ … ␣ 𝑎2𝑛
𝑎31 ␣ 𝑎32 ␣ 𝑎33 ␣ … ␣ 𝑎3𝑛
⋮
𝑎𝑚1 ␣ 𝑎𝑚2 ␣ 𝑎𝑚3 ␣ … ␣ 𝑎𝑚𝑛
𝑏1 ␣ 𝑏2 ␣ 𝑏3 ␣ … ␣ 𝑏𝑚
Obs.: o caractere “␣” representa um espaço em branco.
A partir destes dados, o programa deverá imprimir a matriz aumentada original [𝐴 ⋮ 𝑏] e a matriz equivalente
[𝐴′ ⋮ 𝑏′] junto à solução do sistema com os valores obtidos para cada 𝑥𝑖 após a diagonalização. Caso não seja possível
determinar a solução do sistema, o programa deverá exibir uma mensagem informativa.

Método de Gauss-Seidel:
O programa deve ler arquivos de texto com os dados de entrada digitados da seguinte maneira:
𝑛
𝑎11 ␣ 𝑎12 ␣ 𝑎13 ␣ … ␣ 𝑎1𝑛
𝑎21 ␣ 𝑎22 ␣ 𝑎23 ␣ … ␣ 𝑎2𝑛
𝑎31 ␣ 𝑎32 ␣ 𝑎33 ␣ … ␣ 𝑎3𝑛
⋮
𝑎𝑚1 ␣ 𝑎𝑚2 ␣ 𝑎𝑚3 ␣ … ␣ 𝑎𝑚𝑛
𝑏1 ␣ 𝑏2 ␣ 𝑏3 ␣ … ␣ 𝑏𝑚
𝑘
𝜀
Obs.: o caractere “␣” representa um espaço em branco.
A partir destes dados, o programa deverá calcular o Critério de Convergência de Sassenfeld e imprimir se há
ou não a certeza de que o Método de Gauss-Seidel convergirá para a solução do sistema. Em seguida, o programa
deverá imprimir o sistema 𝑥 = 𝐹𝑥 + 𝑑 gerado e, para toda equação 𝑖 resolvida durante cada iteração 𝑘, deverá
imprimir o 𝑥𝑖 obtido. Ao final de cada iteração 𝑘, o programa deverá analisar se a condição do critério de parada 𝜀 foi
satisfeita. Caso afirmativo, o programa deverá parar e apresentar a solução obtida. Caso negativo, o programa deverá
parar apenas quando chegar à iteração 𝑘 e apresentar a solução aproximada obtida.*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

int maxtam= 0;

void pause (){//funçao de pausar o sistema linux
    int ch;
    while((ch = fgetc(stdin)) != EOF && ch != '\n');//ja limpa o buffer antes
    printf ("\nPressione qualquer tecla para continuar...");
    scanf("%*c");//não PRECISO LIMPAR O BUFFER porque O USUARIO não VAI DIGITAR NADA
}

void fflush_stdin(){//funçao que limpa o buff
    int ch;
    while ((ch = getchar()) != '\n' && ch != EOF);
}

//1-Método Algébrico de Gauss-Jordan
//2-Método Iterativo de Gauss-Seidel
//3-Ambos os métodos deverão tratar a ocorrência de coeficientes-pivôs nulos 
//ou pequenos por meio do Método de Pivotamento Parcial.

int main(){//funcao principal do programa
    int  i= 0,j= 0,cont= 0,op= 1;
    char string[200];
    
	FILE *file; //declaracao do ponteiro arquivo para o arquivo 1
    file= fopen("Inputs3.txt","r");//abre o arquivo 
    
    if(file==NULL){//verifica se o file esta abrindo corretamente
        printf("nao foi possivel abrir o arquivo.\n");
        getchar();//pausa
        exit(0);//sai do programa
    }
   
    fscanf(file,"%d",&maxtam);//le um valor de uma variavel do arquivo como se o usuario estivesse digitado
    double matriz[maxtam][maxtam],vetor[maxtam];//declara a matriz e o vetor b no tamanho lido
    
    for(i=0; i<maxtam; i++){//preenche a matriz com os valores do arquivo
        for(j=0; j<maxtam; j++){
            fscanf(file,"%lf",&matriz[i][j]);
        }     
    }   

    for(i=0; i<maxtam; i++){//preenche o vetor b com os valores do arquivo
        fscanf(file,"%lf",&vetor[i]);
    }

    fclose(file);//fecha o arquivo para evitar erros

    while(op!=0){//laço para gerir a hora de sair do menu 
        
        printf("\nMENU\n0- Sair\n1-opçao para mostrar a matriz \n2-opçao para executar o metodo de Gauss-Jordan\n3-opçao para executar o metodo de Gauss-Seidel\n");
        scanf("%d",&op);
        system("clear||cls");//limpar a tela

        switch (op){//switch para seleçao de qual opçao seguir
            case 1:        
                for(i=0; i<maxtam; i++){//mostra a matriz na tela
                    for(j=0; j<maxtam; j++){
                        printf(" %.2lf",matriz[ i ][ j ]);
                    }
                    printf ("= %.2lf \n",vetor[i]);
                }
                pause();
            break;

            case 2: 
                printf("\n2-\n");
                
                pause();
            break;

            case 3: 
                printf("\n3-\n");
                pause();
            break;

            default:
                if(op!=0){
                    printf("\nopçao invalida...\n");
                    pause();
                }
            break;
        }
        system("clear||cls");//limpa a tela tanto no window quanto no linux
    }

    return 0;//encerra o programa
}