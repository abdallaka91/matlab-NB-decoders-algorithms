#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "./include/polarcodec.h"

void Encoder (code_t *code, table_t *table, int *KSYMB, int *NSYMB, int **NBIN){
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
    }
    for(k=0;k<code->K;k++){
//        printf("%d ", KSYMB[k]);
        u_symb[code->reliability_sequence[code->N-K+k]]=KSYMB[k];
        //printf("%d ", KSYMB[k]);
    }
//    for(n=0;n<N;n++)
//        printf("%d ", u_symb[n]);
    //getchar();

    m=1; // Size of symbols to combine
    int nb_oper, group_n, group_index,ind1, ind2;
    for(d = 1; d <=depth; d++){
         m = 1 << d;
//         printf("m=%d\n", m);
         nb_oper = m >> 1;
        for(i=0; i<N/2; i++){
            group_n = i/nb_oper;
            group_index = i - nb_oper*group_n;
            ind1 = group_n*m + group_index;
            ind2 = group_n*m + nb_oper + group_index;
            temp_symb[ind1]=table->ADDGF[u_symb[ind1]][u_symb[ind2]];
            temp_symb[ind2]=table->MULGF[code->polar_coeff[d-1][i]][u_symb[ind2]];
        }

        for(n=0;n<N;n++){
            u_symb[n]=temp_symb[n];
//            printf("%d ", u_symb[n]);
        }
//        getchar();

    }

    for(n=0;n<N;n++){
        NSYMB[n]=u_symb[n];
        for(int q=0;q<logGF;q++)
            NBIN[n][q]= table->BINGF[NSYMB[n]][q];
    }

//    printf("\n");
//    for(n=0;n<code->N;n++)
//    printf("%d ", NSYMB[n]);
//    getchar();

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

void GenieAidedDecoderSC (decoder_t *decoder, code_t *code, table_t *table, int* KSYMB){
    int N,K,GF,n;

    N=code->N;
    K=code->K;
    GF=code->GF;
    n=code->n;

    int node, depth;
    int decoding_status=0; // indicates whether the decoder has finished.

    int data_positions[K];
    for(int i=0;i<K;i++){
        data_positions[i]=code->reliability_sequence[N-K+i];
    }

    float LLRs[n+1][N][GF]; // Used for computing LLR beliefs.
    int states[2*N-1]; // Used for Node states.
    int real_u[n+1][N];

    for(int i=0;i<2*N-1;i++){
        states[i]=0;
    }

    node=0; depth=0; // start at root
    int bit_len=0;
    int  node_position=0;

    float **L_node;
    L_node=calloc((size_t)N,sizeof(float *));
    L_node[0]= calloc((size_t)N*GF,sizeof(float));
    for (int k=1; k<N; k++){
        L_node[k] = L_node[0] + k*GF;
    }
    float **temp_func;
    temp_func=calloc((size_t)(N/2),sizeof(float *));
    temp_func[0]= calloc((size_t)(N/2)*GF,sizeof(float));
    for (int k=1; k<N/2; k++){
        temp_func[k] = temp_func[0] + k*GF;
    }

    int left_node;
    int right_node;

    int left_depth;
    int left_len;
    int child_depth;
    int child_len;

    int *temp_ucap=calloc(N, sizeof(int));
    int temp_coeff[code->N];


    for(int i=0; i<=n;i++){
        for(int j=0;j<N;j++){
            real_u[i][j]=-1;
            for(int q=0;q<GF;q++)
                LLRs[i][j][q]=10e5;
        }
    }

    for(int i=0;i<N;i++){
        real_u[n][i]=0;
    }

    for(int i=0;i<K;i++){
        real_u[n][data_positions[i]]=KSYMB[i];
    }

    for(int i=0;i<N;i++){
        for(int q=0;q<GF;q++)
            LLRs[0][i][q]=decoder->intrinsic_LLR[i][q]; // Beliefs of the root;
    }

    while(decoding_status==0){ //Traversal Loop. Keeps Traversing until all symbols are decoded.
        if(depth==n){ //leaf procedure
            for(int q=0;q<GF;q++){
                decoder->extrinsic_LLR[node][q]=LLRs[depth][node][q];
            }
            if(node==N-1){
                decoding_status =1;
                break;
            }
            else{
                node=floor(node/2);
                depth=depth-1;
            }
        }
        else{// Non-leaf procedure
            bit_len=pow(2, n-depth); // Bit length is halved over layers
            node_position=pow(2,depth)-1+node; // position of node in node state array.
            if(states[node_position]==0){ // At state 0, perform min-sum approx for beliefs.
                for(int i=0;i<bit_len;i++){
                    for(int q=0;q<GF;q++){
                        L_node[i][q]=LLRs[depth][node*bit_len+i][q];
                    }
                    if(i<bit_len/2){
                        temp_coeff[i]=code->polar_coeff[code->n-depth-1][node*(bit_len/2)+i];
//                    printf("Coef[%d][%d] permutes L[%d][%d]\n", code->n-1-depth, node*bit_len/2+i, depth,node*bit_len+bit_len/2+i);
                    }
                }
                SC_CN(L_node, bit_len, temp_coeff, table, GF, temp_func);
                node=2*node;
                depth=depth+1;
                bit_len=bit_len/2;
                for(int i=0;i<bit_len;i++){
                    for(int q=0;q<GF;q++){
                        LLRs[depth][node*bit_len+i][q]=temp_func[i][q];
                    }
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
                    temp_ucap[i]= real_u[left_depth][left_node*left_len+i];
                    temp_coeff[i]=code->polar_coeff[code->n-depth-1][node*left_len+i];
//                    printf("%d ", temp_ucap[i]);
                }
//                getchar();
                node=node*2+1;
                depth=depth+1;
                SC_VN(L_node, temp_ucap, bit_len, temp_coeff, table, GF, temp_func);
                for(int i=0;i<left_len;i++){
                for(int q=0;q<GF;q++){
                        LLRs[depth][(node*left_len)+i][q]=temp_func[i][q];
                    }
                }
                bit_len=bit_len/2;
                states[node_position]=2;
            }

            else{
                left_node=2*node;
                right_node=2*node+1;
                child_depth=depth+1;
                child_len=bit_len/2;
                for(int i=0;i<child_len;i++){
                    temp_ucap[i]=table->ADDGF[real_u[child_depth][left_node*child_len+i]][real_u[child_depth][right_node*child_len+i]];
                    temp_ucap[child_len+i]=table->MULGF[code->polar_coeff[code->n-1-depth][node*child_len+i]][real_u[child_depth][right_node*child_len+i]];
                }

                for(int i=0;i<bit_len;i++){
                    real_u[depth][(node*bit_len)+i]=temp_ucap[i];
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
    return;
}





