//meta exploid
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
#include <math.h>

int maxtam= 0;

void pause (){//funçao de pausar o sistema linux
    int ch;
   // while((ch = fgetc(stdin)) != EOF && ch != '\n');//ja limpa o buffer antes
    printf ("\nPressione qualquer tecla para continuar...");
    scanf("%*c");//não PRECISO LIMPAR O BUFFER porque O USUARIO não VAI DIGITAR NADA
}

void fflush_stdin(){//funçao que limpa o buff
    int ch;
    while ((ch = getchar()) != '\n' && ch != EOF);
}

void pivotamento(double matrizaux[maxtam][maxtam+1],int linha){
    int j=0;

    double maior=fabs(matrizaux[linha][linha]), aux[maxtam+1];
    int posicao=linha;
    for(int i=linha;i<maxtam;i++){
        j=linha;
        if(maior < fabs(matrizaux[i][j])){
            maior=fabs(matrizaux[i][j]);
            posicao=i;
        }
    }
    if(posicao!=linha){
        for (int j =0; j <=maxtam; j++){
            aux[j]=matrizaux[posicao][j];
        }
        for (int j = 0; j <=maxtam; j++){
            matrizaux[posicao][j]=matrizaux[linha][j];
            matrizaux[linha][j]=aux[j];
        }   
    }
}

void imprime(double matriz[maxtam][maxtam+1]){
    for(int i=0; i<maxtam; i++){//mostra a matriz na tela
        for(int j=0; j<=maxtam; j++){
            if(j==maxtam){
                printf(" =");
            }
            printf(" %.2lf",matriz[ i ][ j ]);
        }
        printf ("\n");
    }
}
/*
matriz[3][4]=
     j   j   j =  j
i    2   3  -1 = -7
i    1   1   1 = -1
i   -1  -2   3 = 15

*/
void gauss_jordan(double matrizaux[maxtam][maxtam+1]){
    int i=0, j=0,k=0;
    double v[maxtam+1],ajj,aij;
    for (j=0;j<maxtam;j++){//linhas
        pivotamento(matrizaux,j);
        ajj=matrizaux[j][j];
        for (k= 0; k <=maxtam; k++){
            if(matrizaux[j][k]==0||ajj==0){
                matrizaux[j][k]=0;
            }else{
                matrizaux[j][k]=matrizaux[j][k]/ajj;;//𝐿𝑗 ← 𝐿𝑗/𝑎𝑗𝑗 ;
            }
            v[k]=matrizaux[j][k];//𝑉 ← 𝐿𝑗 ;
        }

        for (i=0;i<maxtam;i++) {
            if (i!=j){
                aij=matrizaux[i][j];
                for(k=0;k<=maxtam;k++){
                    matrizaux[j][k]=matrizaux[j][k]*aij;// 𝐿𝑗 ← 𝐿𝑗 ∗ 𝑎𝑖𝑗 ;
                    matrizaux[i][k]=matrizaux[i][k]-matrizaux[j][k];//𝐿𝑖 ← 𝐿𝑖 − 𝐿𝑗 ;
                    matrizaux[j][k]=v[k];//𝐿𝑗 ← 𝑉;
                }
            }
        }
    }
}

//1-Método Algébrico de Gauss-Jordan
//2-Método Iterativo de Gauss-Seidel
//3-Ambos os métodos deverão tratar a ocorrência de coeficientes-pivôs nulos 
//ou pequenos por meio do Método de Pivotamento Parcial.

int main(){//funcao principal do programa
    int  i= 0,j= 0;
    
	FILE *file; //declaracao do ponteiro arquivo para o arquivo 1
    file= fopen("Inputs3-2.txt","r");//abre o arquivo 
    
    if(file==NULL){//verifica se o file esta abrindo corretamente
        printf("nao foi possivel abrir o arquivo.\n");
        getchar();//pausa
        exit(0);//sai do programa
    }
   
    fscanf(file,"%d",&maxtam);//le um valor de uma variavel do arquivo como se o usuario estivesse digitado
    double matriz[maxtam][maxtam+1],vetor[maxtam],matrizaux[maxtam][maxtam+1];//declara a matriz e o vetor b no tamanho lido
    
    for(i=0; i<maxtam; i++){//preenche a matriz com os valores do arquivo
        for(j=0; j<maxtam; j++){
            fscanf(file,"%lf",&matriz[i][j]);
        }     
    }   

    for(i=0; i<maxtam; i++){//preenche o vetor b com os valores do arquivo
        fscanf(file,"%lf",&vetor[i]);
    }

    fclose(file);//fecha o arquivo para evitar erros

    j=maxtam;
    for(i=0; i<maxtam; i++){//adiciona o vetor b a matriz na ultima coluna
       matriz[ i ][ j ]=vetor[i];
    }

    for(i=0; i<maxtam; i++){//mostra a matriz na tela
        for(j=0; j<=maxtam; j++){
            matrizaux[i][j]=matriz[ i ][ j ];
        }
    }
       
    imprime(matriz);
    printf ("\n");
    gauss_jordan(matrizaux);
    imprime(matrizaux);
    printf ("\n");

    pause();

    return 0;//encerra o programa
}