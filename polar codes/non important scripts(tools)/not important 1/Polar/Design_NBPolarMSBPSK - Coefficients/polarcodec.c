#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "./include/polarcodec.h"

void Encoder (code_t *code, table_t *table, int *KSYMB, int *NSYMB, int **NBIN, decoder_t *decoder){
    int k,n,d,m,i;
    int K,N, logGF;
    int depth;

    K=code->K;
    N=code->N;
    logGF=code->logGF;

    depth=code->n;

    int u_symb[N];
    int temp_symb[N];

    for(n=0;n<code->N;n++){
        u_symb[n]=0;
        decoder->encoded_symbols[depth][n]=0;
    }
    for(k=0;k<code->K;k++){
//        printf("%d ", KSYMB[k]);
        u_symb[code->reliability_sequence[code->N-K+k]]=KSYMB[k];
        decoder->encoded_symbols[depth][code->reliability_sequence[code->N-K+k]]=KSYMB[k];
    }

    /*for(int i=0;i<depth;i++){
        for(int j=0;j<N/2;j++){
            code->polar_coeff[i][j]=table->DECGF[code->polar_coeff[i][j]];
        }
    }*/

    m=1; // Size of symbols to combine
    int nb_oper, group_n, group_index,ind1, ind2;
    for(d = 1; d <=depth; d++){
         m = 1 << d;
         nb_oper = m >> 1;
        for(i=0; i<N/2; i++){
            group_n = i/nb_oper;
            group_index = i - nb_oper*group_n;
            ind1 = group_n*m + group_index;
            ind2 = group_n*m + nb_oper + group_index;
            temp_symb[ind1]=table->ADDGF[u_symb[ind1]][u_symb[ind2]];
            temp_symb[ind2]=table->MULGF[code->polar_coeff[d-1][i]][u_symb[ind2]];
        }

//        printf("\n");
        for(n=0;n<N;n++){
            u_symb[n]=temp_symb[n];
            decoder->encoded_symbols[depth-d][n]=temp_symb[n];
        }
//                    printf("%d ", depth-d);
    }
//            getchar();

    for(int n=0;n<N;n++){
        NSYMB[n]=u_symb[n];
//        printf("%d ", NSYMB[n]);
        for(int q=0;q<logGF;q++){
//            printf("%d ", table->BINGF[NSYMB[n]][q]);
            NBIN[n][q]= table->BINGF[NSYMB[n]][q];
        }
    }

//    printf("\n");
//    for(n=0;n<code->N;n++)
//    printf("%d ", NSYMB[n]);
//    getchar();

    return;
}


void DecoderSC (decoder_t *decoder, code_t *code, table_t *table, float err[code->N/2][code->GF], int dep_ind, int node_ind){
    int N,K,GF,n;

    N=code->N;
    K=code->K;
    GF=code->GF;
    n=code->n;

    int node, depth;
    int decoding_status=0; // indicates whether the decoder has finished.

    int check_vector[N];
    for(int i=0;i<N;i++){
        if(i<N-K){
            check_vector[code->reliability_sequence[i]]=1; //frozen position
        }
        else{
            check_vector[code->reliability_sequence[i]]=0; // data position
        }
    }

    float LLRs[n+1][N][GF]; // Used for computing LLR beliefs.
    int states[2*N-1]; // Used for Node states.
    int ucap[n+1][N]; // Used for computing the decisions

    for(int i=0;i<2*N-1;i++){
        states[i]=0;
    }

    node=0; depth=0; // start at root
    int bit_len=0;
    int  node_position=0;

    float **L_node;
    L_node=calloc((size_t)N,sizeof(float *));
    L_node[0]= calloc((size_t)N*GF,sizeof(float));
    for (int k=1; k<N; k++)
        L_node[k] = L_node[0] + k*GF;

    float **temp_func;
    temp_func=calloc((size_t)(N),sizeof(float *));
    temp_func[0]= calloc((size_t)(N)*GF,sizeof(float));
    for (int k=1; k<N; k++)
        temp_func[k] = temp_func[0] + k*GF;

    int left_node;
    int right_node;

    int left_depth;
    int left_len;
    int child_depth;
    int child_len;
    int *temp_coeff=calloc(N, sizeof(int));

    int *temp_ucap=calloc(N, sizeof(int));

    int *real_ucap=calloc(N, sizeof(int));

    for(int i=0; i<=n;i++){
        for(int j=0;j<N;j++){
            ucap[i][j]=-1;
            for(int q=0;q<GF;q++)
                LLRs[i][j][q]=10e5;
        }
    }

    for(int i=0;i<N;i++){
        for(int q=0;q<GF;q++)
            LLRs[0][i][q]=decoder->intrinsic_LLR[i][q]; // Beliefs of the root;
    }

    float min;

    while(decoding_status==0){ //Traversal Loop. Keeps Traversing until all bits are decoded.
//            printf("%d %d", depth, node);
//    getchar();
        if(depth==n){ //leaf procedure
//             if frozen procedure ucap is always 0
            if(check_vector[node]==1){
                ucap[depth][node]=0;
            }
            else{ //message position so check LLR
                min=10e5;
                for(int q=0;q<GF;q++){
                    if(LLRs[depth][node][q]<min){
                        min=LLRs[depth][node][q];
                        ucap[depth][node]=q;
                    }
                }
            }
//            printf("u[%d]= %d", node, ucap[depth][node]);
//            getchar();
            if(node==N-1){
                decoding_status =1;
            }
            else{
                node=floor(node/2);
                depth=depth-1;
            }
        }
        else{// Non-leaf procedure
            bit_len=pow(2, n-depth); // Bit length is halved over layers
            node_position=pow(2,depth)-1+node; // position of node in node state array.
            if(states[node_position]==0){ //Left-node processing
//                                printf("\nU' U\n ");
                for(int i=0;i<bit_len;i++){
                    for(int q=0;q<GF;q++){
                        L_node[i][q]=LLRs[depth][node*bit_len+i][q];
//                        if(L_node[i][q]==0){
//                            printf("%d ", q);
//                        }
                    }
                    if(i<bit_len/2){
                        temp_coeff[i]=code->polar_coeff[code->n-depth-1][node*(bit_len/2)+i];
                    }
                    real_ucap[i]=decoder->encoded_symbols[depth][node*bit_len+i];
                }
//                getchar();
                SC_CN(L_node, bit_len, temp_coeff, table, GF, temp_func);
                node=2*node;
                depth=depth+1;
                bit_len=bit_len/2;
                for(int i=0;i<bit_len;i++){
                    for(int q=0;q<GF;q++){
                        LLRs[depth][node*bit_len+i][q]=temp_func[i][q];
                    }
//                    if(depth==dep_ind && node == node_ind){
//                        err[node*bit_len+i][temp_coeff[i]]=fabs(computeErrProb(L_node[i], real_ucap[i], GF)-computeErrProb(L_node[bit_len+i], real_ucap[bit_len+i], GF));
////                        printf("Result saved at ERR[%d][%d][%d] = %g\n", depth-1, node*bit_len+i, temp_coeff[i], err[node*bit_len+i][temp_coeff[i]]);
//                    }

                }
                states[node_position]=1;
            }
            else if (states[node_position]==1){//Right leaf processing
                for(int i=0;i<bit_len;i++){
                    for(int q=0;q<GF;q++){
                        L_node[i][q]=LLRs[depth][node*bit_len+i][q];
                    }
                }
                left_node=2*node;
                left_depth=depth+1;
                left_len=bit_len/2;

                for(int i=0;i<left_len;i++){
                    temp_ucap[i]= ucap[left_depth][left_node*left_len+i];
                    temp_coeff[i]=code->polar_coeff[code->n-depth-1][node*left_len+i];
                    real_ucap[i]=decoder->encoded_symbols[left_depth][left_node*left_len+i];
                }
                node=node*2+1;
                depth=depth+1;
                SC_VN(L_node, temp_ucap, bit_len, temp_coeff, table, GF, temp_func);
                for(int i=0;i<left_len;i++){
                    for(int q=0;q<GF;q++){
                        LLRs[depth][(node*left_len)+i][q]=temp_func[i][q];
                    }
                }
                if(depth-1==dep_ind && (node/2)==node_ind){
//                    printf("D: %d len:%d combined to %d\n", depth-1, bit_len, bit_len/2);
                    for(int i=0;i<bit_len;i++){
//                        printf("LLR[%d][%d] with U[%d][%d] processed with LLR[%d][%d] with U[%d][%d]!\t", depth, (node/2)*bit_len+i, depth, (node/2)*bit_len+i, depth, node*(bit_len/2)+i, depth, node*(bit_len/2)+i);
                        real_ucap[i]=decoder->encoded_symbols[depth][(node/2*bit_len)+i];
                        for(int q=0;q<GF;q++)
                            temp_func[i][q]=LLRs[depth][(node/2*bit_len)+i][q];
                    }
                    for(int i=0;i<bit_len/2;i++){
                        //err[i][temp_coeff[i]]=fabs(computeErrProb(temp_func[i], real_ucap[i], GF)-computeErrProb(temp_func[bit_len/2+i], real_ucap[bit_len/2+i], GF));
                        err[i][temp_coeff[i]]=fabs(computeEntropy(temp_func[i], GF)-computeEntropy(temp_func[bit_len/2+i],GF));
//                        if(temp_coeff[i]==1)
//                            printf("ERR[%d][%d] is obtained by ERR[%d][%d] with ERR[%d][%d]\n", depth-1, i, depth-1, (node/2*bit_len)+i, depth-1, (node/2*bit_len)+bit_len/2+i);
                    }
                    return;
                }

//                if(depth==l)
//                getchar();
                bit_len=bit_len/2;
                states[node_position]=2;
            }

            else{
                left_node=2*node;
                right_node=2*node+1;
                child_depth=depth+1;
                child_len=bit_len/2;
//                printf("Depth: %d, node:%d\n", depth, node);
                for(int i=0;i<child_len;i++){
                    temp_ucap[i]=table->ADDGF[ucap[child_depth][left_node*child_len+i]][ucap[child_depth][right_node*child_len+i]];
                    temp_ucap[child_len+i]=table->MULGF[code->polar_coeff[code->n-1-depth][node*child_len+i]][ucap[child_depth][right_node*child_len+i]];
//                    temp_ucap[child_len+i]=ucap[child_depth][right_node*child_len+i];
//                    printf("Coef[%d][%d] divided by u[%d][%d]\n", code->n-1-depth, node*child_len+i, child_depth, right_node*child_len+i);
                }
//                getchar();
                for(int i=0;i<bit_len;i++){
                    ucap[depth][(node*bit_len)+i]=temp_ucap[i];
                }
                node=floor(node/2);
                depth=depth-1;
            }
        }
    }
    free(L_node[0]);
    free(L_node);
    free(temp_func[0]);
    free(temp_func);
    free(temp_ucap);
    free(real_ucap);
    free(temp_coeff);
    return;
}

void SC_VN(float **L, int *c, int arr_size, int *coeff, table_t *table, int GF, float **res){
    float a[GF];
    float b[GF];
    for(int i=0;i<arr_size/2;i++){
        for(int q=0;q<GF;q++){
            a[q]=L[i][q];
            b[q]=L[arr_size/2+i][q];
        }
        G_MS(a, b, c[i], coeff[i],table,GF, res[i]);
    }
    return;
}

void SC_CN(float **L, int arr_size, int* coeff, table_t *table, int GF, float **result){
    float a[GF];
    float b[GF];
    for(int i=0;i<arr_size/2;i++){
        for(int q=0;q<GF;q++){
            a[q]=L[i][q];
            b[q]=L[arr_size/2+i][q];
        }
        F_MS(a, b, coeff[i], table, GF, result[i]);
    }
    return;
}

