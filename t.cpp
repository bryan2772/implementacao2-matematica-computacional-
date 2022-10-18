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

//1-Método Algébrico de Gauss-Jordan
//2-Método Iterativo de Gauss-Seidel
//3-Ambos os métodos deverão tratar a ocorrência de coeficientes-pivôs nulos 
//ou pequenos por meio do Método de Pivotamento Parcial.

int main(){//funcao principal do programa
    int  i=0,j=0,cont=0,maxtam=0;
    char  string[200];
    char *test=NULL;
    
    

	FILE *file; //declaracao do ponteiro arquivo para o arquivo 1
    file= fopen("Inputs1.txt","r");
    
    if(file==NULL){//verifica se o file esta abrindo corretamente
        printf("nao foi possivel abrir o arquivo.\n");
        getchar();//pausa
        exit(0);//sai do programa
    }

    fgets(string,200,file) != NULL;
    test = strtok(string," "); 
    maxtam=atoi(test); 
    double matriz[maxtam][maxtam],vetor[maxtam];

    while((fgets(string,200,file) != NULL) && (cont<=maxtam)){//fgets copia a string do arquivo para o programa

            if(cont<maxtam){
            j=0;
            test = strtok(string," "); 
            matriz[cont][j]= atoi(test); 
            for(j=1;j<maxtam;j++){
                test = strtok(NULL," "); 
                matriz[cont][j]= atoi(test); 
            }
            cont++;
            }else{
            j=0;
            test = strtok(string," "); 
            vetor[j]= atoi(test); 
            for(j=1;j<maxtam;j++){
                test = strtok(NULL," "); 
                vetor[j]= atoi(test); 
            }
            }
           
    }

    fclose(file); //fecha o arquivo para evitar erros
    
    for ( i=0; i<maxtam; i++ ){
        for ( j=0; j<maxtam; j++ ){
            printf (" %lf ", matriz[ i ][ j ]);
        }
        printf ("\n");
    }

    for ( j=0; j<maxtam; j++ ){
            printf (" %lf ", vetor[j]);
        }
   
    return 0;
}