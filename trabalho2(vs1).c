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

void pivotamento(int maxtam, long double matrizaux[maxtam][maxtam+1],int linha){
    int j=0;//funçao que faz o pivotamento ela identifica a maior linha e faz a troca 

    long double maior=fabs(matrizaux[linha][linha]), aux[maxtam+1];//maior recebe o primeiro elemento em modulo
    int posicao=linha;//salva a posiçao da linha inicial
    for(int i=linha;i<maxtam;i++){//percorre as linhas em busca da posiçao maior
        j=linha;
        if(maior < fabs(matrizaux[i][j])){//se a posiçao for maior ele salva ela e o valor
            maior=fabs(matrizaux[i][j]);
            posicao=i;
        }
    }
    if(posicao!=linha){//faz a troca das linhas
        for (int j =0; j <=maxtam; j++){//salva no auxiliar para fazer a substituiçao
            aux[j]=matrizaux[posicao][j];
        }
        for (int j = 0; j <=maxtam; j++){
            matrizaux[posicao][j]=matrizaux[linha][j];
            matrizaux[linha][j]=aux[j];
        }   
    }
}

void pivotamentocompleto(int maxtam,long double matrizaux[maxtam][maxtam+1]){
    long double maior=0,aux=0,multiplicador=0;
    int k,i,j,l,c;

    for(k=0;k<maxtam-1;k++){
        l = c = k;
        maior = fabs(matrizaux[k][k]);
        for(i=k;i<maxtam;i++){// Encontra o maior elemento em modulo  
            for(j=k;j<maxtam;j++){
                if(fabs(matrizaux[i][j]) > maior){
                    maior = fabs(matrizaux[i][j]);
                    l = i;
                    c = j;
                }
            }
        }
        if(l != k){// Se a linha do maior elemento encontrado for maior que o pivo original troca as linhas
            for(i=0;i<=maxtam;i++){
                aux = matrizaux[k][i];
                matrizaux[k][i] = matrizaux[l][i];
                matrizaux[l][i] = aux;
            }
        }
        if(c != k){ //troca a coluna 
            for(i=0;i<maxtam;i++){
                aux = matrizaux[i][k];
                matrizaux[i][k] = matrizaux[i][c];
                matrizaux[i][c] = aux;
            }
        }
        // Zera os elementos abaixo da diagonal principal (escalonamento)
        for(i=k+1;i<maxtam;i++){
            multiplicador = matrizaux[i][k]/matrizaux[k][k];
			matrizaux[i][k] = 0;
			for(j=k+1;j<=maxtam;j++){
				matrizaux[i][j] -= multiplicador * matrizaux[k][j];
			}
		}	
    }
}

void imprime(int maxtam, long double matriz[maxtam][maxtam+1]){//funçao para imprimir a matriz
    for(int i=0; i<maxtam; i++){//mostra a matriz na tela
        for(int j=0; j<=maxtam; j++){
            if(j==maxtam){
                printf(" = %.4Lf",matriz[ i ][ j ]);
               
            }else{
                printf(" %.2Lf",matriz[ i ][ j ]);
            }
            
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
void gauss_jordan(int maxtam,long double matrizaux[maxtam][maxtam+1]){
    int i=0, j=0,k=0;
    long double v[maxtam+1],ajj,aij;
    for (j=0;j<maxtam;j++){//linhas
        //imprime(matrizaux);
        //printf ("\n");
        pivotamento(maxtam,matrizaux,j);//troca a linha maior
        if(matrizaux[j][j]==0){//verifica se a diagonal é 0 mesmo depois da troca
            //if(matrizaux[j][maxtam]==0)
                //printf("\né um SPI\n");
           // else
                //printf("\né um SI\n");
            return;//caso seja 0 ele encerra as execuçoes
        }
        ajj=matrizaux[j][j];//salva a posiçao ajj da matriz pois sera alterada
        for (k= 0; k <=maxtam; k++){//for para percorrer cada elemento da linha j
            if(matrizaux[j][k]==0||ajj==0){//verificaçao se esta fazendo  divisao por 0;
                matrizaux[j][k]=0;
            }else{
                matrizaux[j][k]=matrizaux[j][k]/ajj;;//𝐿𝑗 ← 𝐿𝑗/𝑎𝑗𝑗 ;
            }
            v[k]=matrizaux[j][k];//𝑉 ← 𝐿𝑗 ;]
        }

        for (i=0;i<maxtam;i++) {
            if (i!=j){
                aij=matrizaux[i][j];//salva a posiçao i j da matriz
                for(k=0;k<=maxtam;k++){//percorre cada elemento das linhas 
                    matrizaux[j][k]=matrizaux[j][k]*aij;// 𝐿𝑗 ← 𝐿𝑗 ∗ 𝑎𝑖𝑗 ;
                    matrizaux[i][k]=matrizaux[i][k]-matrizaux[j][k];//𝐿𝑖 ← 𝐿𝑖 − 𝐿𝑗 ;
                    matrizaux[j][k]=v[k];//𝐿𝑗 ← 𝑉;
                }
            }
        }
    }
}

void solucaogaussjordan(int maxtam,long double matrizaux[maxtam][maxtam+1]){

    if(matrizaux[maxtam-1][maxtam-1]==0){//verifica se a ultima diagonal principal e 0
        if(matrizaux[maxtam-1][maxtam]==0){//verifica o ultimo vetor b
            printf("o sistema e um SPI\n");
        }else{
            printf("o sistema e um SI\n");
        }
    }else{
        printf("o sistema e um SPD:\n");
        for (int i=0; i < maxtam; i++){//imprime o resultado do vetor b (Xi)
           printf("x%d: %.10Lf ",i,matrizaux[i][maxtam]);
        }
        printf("\n");
    }

}

void gauss_seidel(int maxtam,long double matrizaux[maxtam][maxtam+1],int K,long double epsilon){
    /*A partir destes dados, o programa deverá calcular o Critério de Convergência de Sassenfeld 
    e imprimir se há ou não a certeza de que o Método de Gauss-Seidel convergirá para a solução 
    do sistema.*/
    int i=0,j=0,k=0;
    long double beta[maxtam],Sassenfeld=0,soma=0,X[maxtam],Xanterior[maxtam],elementomaior=epsilon;
    printf("\nmatriz apos o pivotamento completo:\n");
    pivotamentocompleto(maxtam,matrizaux);
    imprime(maxtam,matrizaux);

    beta[0]=0;
    for(i=1;i<maxtam;i++){//salva o beta como 1 para na atrapalhar na divisao
        beta[i]=1;
    }

    for(j=0;j<maxtam;j++){//operacionaliza linha 1
        beta[0]=beta[j]+matrizaux[0][j];
        if(matrizaux[0][0]!=0){//verifica divisao por zero
            beta[0]=beta[0]/matrizaux[0][0];
        }else{
            beta[0]=0;
        }
    }

    for ( i = 1; i < maxtam; i++){//percorre as linhas
        soma=0;
        if(matrizaux[0][0]!=0){//verifica divisao por zero
            for ( j = 0; j < maxtam; j++){//percorre as colunas
                if(i!=j){//verifica se nao esta na diagonal principal
                    soma+=(fabs(matrizaux[i][j])*beta[j]);//soma cada eleento da linha
                }
                if(matrizaux[i][i]!=0){//verifica divisao por zero
                    beta[i]=soma/fabs(matrizaux[i][i]);
                }else{
                    beta[i]=0;
                }
            }
        }else{
            beta[i]=0;
        }
    }
        
    Sassenfeld=beta[0];
	for(i=1;i<maxtam;i++){
		if(beta[i]>Sassenfeld){
			Sassenfeld=beta[i];
		}
	}
    if(Sassenfeld<1){
        printf("\nhá a certeza de que o Método de Gauss-Seidel convergirá para a solução do sistema.\n beta = %Lf < 1\n\n",Sassenfeld);
    }else{
        printf("\nnão a certeza de que o Método de Gauss-Seidel convergirá para a solução do sistema.\nbeta = %Lf > 1\n\n",Sassenfeld);
    }
    /* Em seguida, o programa deverá imprimir o sistema 𝑥 = 𝐹𝑥 + 𝑑 gerado*/
    
    //imprime(maxtam,matrizaux);
    printf("sistema 𝑥 = 𝐹𝑥 + 𝑑: \n");
    for(i=0;i<maxtam;i++){
        printf("x%d= ",i+1);
        printf(" (%Lf",matrizaux[i][maxtam]);
        
        for(j=0;j<maxtam;j++){
            if(i!=j){
                if(matrizaux[i][j]!=0){
                    matrizaux[i][j]=matrizaux[i][j]*-1;
                    if(matrizaux[i][j]>0){
                        printf(" +%Lfx%d",matrizaux[i][j],j+1);
                    }else{
                        printf(" %Lfx%d",matrizaux[i][j],j+1);
                    }
                }
            }
        }
        printf(") / %Lf \n",matrizaux[i][i]);
    }
     //imprime(maxtam,matrizaux);




/*e, para 
    toda equação 𝑖 resolvida durante cada iteração 𝑘, deverá imprimir o 𝑥𝑖 obtido.
    Ao final de cada iteração 𝑘, o programa deverá analisar se a condição do critério
    de parada 𝜀 foi satisfeita. Caso afirmativo, o programa deverá parar e 
    apresentar a solução obtida. Caso negativo, o programa deverá parar apenas
    quando chegar à iteração 𝑘 e apresentar a solução aproximada obtida.*/




    for(i=0;i<maxtam;i++){
        X[i]=0;
        Xanterior[i]=0;
    }
    soma=0;
    for(k=1;k<=K && elementomaior >= epsilon;k++){
        printf("\nk=%d ",k);
        for(i=0;i<maxtam;i++){ 
            soma=matrizaux[i][maxtam];
            for(j=0;j<maxtam;j++){
                if(i!=j){
                   soma=soma+matrizaux[i][j]*X[j];
                }
            }

            X[i]=soma/matrizaux[i][i];
            printf(" X%d = %Lf  ",X[i],i+1);
            soma=0;
        }

        for(i=0;i<maxtam;i++){
            Xanterior[i]=X[i];
        }

    }
}

//11 15
//1-Método Algébrico de Gauss-Jordan =9,10,13,14--5 7 12
//2-Método Iterativo de Gauss-Seidel
//3-Ambos os métodos deverão tratar a ocorrência de coeficientes-pivôs nulos 
//ou pequenos por meio do Método de Pivotamento Parcial.

int main(){//funcao principal do programa
    int  i= 0,j= 0,K=0;
    int maxtam= 0;
    long double epsilon=0.005;
    //float epsilon=0.005;
	FILE *file; //declaracao do ponteiro arquivo para o arquivo 1
    file= fopen("Inputs3.txt","r");//abre o arquivo 
    
    if(file==NULL){//verifica se o file esta abrindo corretamente
        printf("nao foi possivel abrir o arquivo.\n");
        getchar();//pausa
        exit(0);//sai do programa
    }
   
    fscanf(file,"%d",&maxtam);//le um valor de uma variavel do arquivo como se o usuario estivesse digitado
    long double matriz[maxtam][maxtam+1],vetor[maxtam],matrizaux[maxtam][maxtam+1];//declara a matriz e o vetor b no tamanho lido
    
    for(i=0; i<maxtam; i++){//preenche a matriz com os valores do arquivo
        for(j=0; j<maxtam; j++){
            fscanf(file,"%Lf",&matriz[i][j]);
        }     
    }   

    for(i=0; i<maxtam; i++){//preenche o vetor b com os valores do arquivo
        fscanf(file,"%Lf",&vetor[i]);
    }
    fscanf(file,"%d",&K);
    fscanf(file,"%Lf",&epsilon);

    fclose(file);//fecha o arquivo para evitar erros
    
    j=maxtam;
    for(i=0; i<maxtam; i++){//adiciona o vetor b a matriz na ultima coluna
       matriz[ i ][ j ]=vetor[i];
    }

    for(i=0; i<maxtam; i++){//cria uma matriz auxiliar
        for(j=0; j<=maxtam; j++){
            matrizaux[i][j]=matriz[ i ][ j ];
        }
    }

    printf ("Matriz aumentada original [𝐴 ⋮ 𝑏] : \n");
        imprime(maxtam,matriz);
        printf ("\n");

    gauss_jordan(maxtam, matrizaux);

    printf ("Matriz equivalente [𝐴′ ⋮ 𝑏′] : \n");
        imprime(maxtam,matrizaux);
        printf ("\n");

    solucaogaussjordan(maxtam,matrizaux);

    for(i=0; i<maxtam; i++){//cria uma matriz auxiliar
        for(j=0; j<=maxtam; j++){
            matrizaux[i][j]=matriz[ i ][ j ];
        }
    }
    pause();
    gauss_seidel(maxtam,matrizaux,K,epsilon);
    //printf(" k=%d ,𝜀=%Lf ",K,epsilon);
    pause();

    return 0;//encerra o programa
}