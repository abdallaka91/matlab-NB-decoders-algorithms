// initialization
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "./include/struct.h"
#include "./include/init.h"
#include "./include/tools.h"

#define STR_MAXSIZE 350

void AllocateDecoder (code_t *code, decoder_t *decoder)
{
    const int N = code->N;
    int nbRow, nbCol;
    int k;

    nbRow = N;
    nbCol = code->GF;
//    /* APP [N][GF] */
//    decoder->APP =calloc((size_t)nbRow,sizeof(softdata_t *));
//    //if (decoder->APP  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
//    decoder->APP [0] = calloc((size_t)nbRow*nbCol,sizeof(softdata_t));
//    //if (decoder->APP [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
//    for (k=1; k<nbRow; k++) decoder->APP[k] = decoder->APP[0] + k*nbCol;

    /* intrinsic_LLR [N][GF] */
    decoder->intrinsic_LLR =calloc((size_t)nbRow,sizeof(softdata_t *));
    decoder->intrinsic_LLR [0] = calloc((size_t)nbRow*nbCol,sizeof(softdata_t));
    for (k=1; k<nbRow; k++) decoder->intrinsic_LLR[k] = decoder->intrinsic_LLR[0] + k*nbCol;

    /* intrinsic_GF [N][GF] */
    decoder->intrinsic_GF 	= calloc((size_t)nbRow,sizeof(int *));
    decoder->intrinsic_GF [0] 	= calloc((size_t)nbRow*nbCol,sizeof(int));
    for (k=1; k<nbRow; k++) decoder->intrinsic_GF[k] = decoder->intrinsic_GF[0] + k*nbCol;

    decoder->extrinsic_LLR =calloc((size_t)nbRow,sizeof(softdata_t *));
    decoder->extrinsic_LLR [0] = calloc((size_t)nbRow*nbCol,sizeof(softdata_t));
    for (k=1; k<nbRow; k++) decoder->extrinsic_LLR[k] = decoder->extrinsic_LLR[0] + k*nbCol;

}

void FreeDecoder (decoder_t *decoder)
{
    free(decoder->intrinsic_LLR[0] );
    free(decoder->intrinsic_LLR);
    free(decoder->intrinsic_GF[0] );
    free(decoder->intrinsic_GF);
    free(decoder->extrinsic_LLR[0]);
    free(decoder->extrinsic_LLR);

}


void Table_Add_GF(table_t *table, int GF, int logGF)
{
    int i,j,k;
    int temp[12];

    for(j=0; j<GF; j++)
    {
        for(k=0; k<GF; k++)
        {
            for(i=0; i<logGF; i++)
            {
                temp[i] = (table->BINGF[j][i])^(table->BINGF[k][i]);
            }
            table->ADDGF[j][k] = Bin2GF(temp,GF,logGF,table->BINGF);
        }
//        printf("%d",j);
    }
}

// multiply GF values and output decimal
void Table_Mul_DEC(table_t *table, int GF)
{

    int i,j;
    for(i=0; i<GF; i++)
    {
        for(j=0; j<GF; j++)
        {
        //table->MULDEC[table->DECGF[i]][table->DECGF[j]]=table->DECGF[table->MULGF[i][j]];
        table->MULDEC[i][j]=table->DECGF[table->MULGF[i][j]];
        }
    }

//    for(i=0; i<GF; i++)
//    {
//        for(j=0; j<GF; j++)
//        {
//          printf("%d ",table->MULDEC[i][j]);
//        }
//        printf(" \n ");
//    }
//    getchar();

}

//divide decimal by GF and output GF
void Table_Div_DEC(table_t *table, int GF)
{

    int i,j;
    for(i=0; i<GF; i++)
    {
        for(j=0; j<GF; j++)
        {
        //table->DIVDEC[table->DECGF[i]][table->DECGF[j]]=table->DECGF[table->DIVGF[i][j]];
        table->DIVDEC[table->DECGF[i]][j]=table->DIVGF[i][j];
        }
    }

//    for(i=0; i<GF; i++)
//    {
//        for(j=0; j<GF; j++)
//        {
//          printf("%d ",table->DIVDEC[i][j]);
//        }
//        printf(" \n ");
//    }
//    getchar();

}




void Table_dec_GF(table_t *table, int GF, int logGF)
{
    int i,j;

    //bin2dec
    int sum;
    int tmp;
    for (j=0; j<GF; j++)
    {
        sum = 0;
        for (i=0; i<logGF; i++)
        {
            tmp = table->BINGF[j][i];
            //printf("%d",tmp);
            sum =sum + (tmp<<i);
        }
        table->DECGF[j]=sum;
        table->GFDEC[sum]=j;
        //printf(" \n bin2dec of GF %d is %d \n",j,sum);
    }
    //getchar();
}



/*!
 * \fn Table_Mul_GF
 * \brief compute the multiplication table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
void Table_Mul_GF(int **MULGF, int GF)
{
    int i,j,temp;
    for(i=0; i<GF; i++)
    {
        for(j=0; j<GF; j++)
        {
            if (i==0 || j==0)
                MULGF[i][j] = 0;
            else if (i==1)
                MULGF[i][j] = j;
            else if (j==1)
                MULGF[i][j] = i;
            else
            {
                temp=i+j-2;
                if(temp<GF-1)
                    MULGF[i][j] = temp+1;
                else
                    MULGF[i][j]=(temp%(GF-1))+1;
            }
        }
    }
//    for(i=0; i<GF; i++)
//    {
//        for(j=0; j<GF; j++)
//        {
//          printf("%d ",MULGF[i][j]);
//        }
//        printf(" \n ");
//    }
//    getchar();

}




/*!
 * \fn Table_Div_GF
 * \brief compute the division table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
void Table_Div_GF(int **DIVGF, int GF)
{
    int i,j,nb;
    nb=GF-1;
    for(i=0; i<GF; i++)
    {
        for(j=0; j<GF; j++)
        {
            if(j==0)
            {
                DIVGF[i][j]=0;
            }
            else if(i==0)
            {
                DIVGF[i][j]=0;
            }
            else if(j==1)
            {
                DIVGF[i][j]=i;
            }
            else
            {
                DIVGF[i][j]=nb--;
            };
            if(nb<1)
            {
                nb=GF-1;
            }
        }
    }
//    for(i=0; i<GF; i++)
//    {
//        for(j=0; j<GF; j++)
//        {
//          printf("%d ",DIVGF[i][j]);
//        }
//        printf(" \n ");
//    }
//    getchar();
}

void LoadCode (code_t *code, float SNR)
{
    int N=code->N;
    code->K=N;
    if(ceil(log2(N)) != floor(log2(N)) || ceil(log2(code->GF)) != floor(log2(code->GF))){
        printf("\n\t\t------N should be power of 2---------\n");
        exit(-1001);
    }
    code->logGF=log2(code->GF);
    code->n=log2(N);
    code->reliability_sequence=calloc(N, sizeof(int));
    code->ErrorProb=calloc(N, sizeof(float));
    code->Entropy=calloc(N, sizeof(float *));
    for(int i=0;i<N;i++){
        code->reliability_sequence[i]=i;
        code->ErrorProb[i]=0;
        code->Entropy[i]=0;
    }

    FILE *opfile;
    char fname[70];
    sprintf (fname,"./matrices/N%d/mat_N%d_GF%d_SNR%.2f.txt",code->N, code->N,code->GF, SNR);
    printf("%s\n", fname);
    opfile=fopen(fname,"r");
    if(opfile==NULL){
        printf("Sequence Unavailable!!\n");
        exit(-1010);
    }

    code->polar_coeff=calloc((size_t)code->n,sizeof(int *));
    code->polar_coeff[0] = calloc((size_t)code->n*code->N/2,sizeof(int));
    for (int k=1; k<code->n; k++){
        code->polar_coeff[k] = code->polar_coeff[0] + k*code->N/2;
    }

    int tmp;
    for(int i=0;i<code->n;i++){
        for(int j=0;j<code->N/2;j++){
            fscanf(opfile,"%d",&tmp);
            code->polar_coeff[i][j]=tmp;
        }
    }
	fclose(opfile);




}

void FreeCode(code_t *code)
{
    free(code->reliability_sequence);
    free(code->ErrorProb);
    free(code->Entropy);
    free(code->polar_coeff[0]);
    free(code->polar_coeff);
}

/**
 * \fn void LoadTables (table_t *table, int GF, int logGF)
 * \brief Memory allocation for the tables and Initialization of the tables.
 * Inputs        :
 * 	- table_t table : Structure that contains pointers to the tables
 * 	- int GF    : order of the field
 * 	- int logGF : logGF = log2(GF)
 * Outputs       :
 */
void LoadTables (table_t *table, int GF, int logGF)
{
    int nbRow, nbCol, g,k,l;

    if(GF!=16 && GF!=32 && GF!=64 && GF!=256 && GF!=4096)
    {
        printf("The binary image of GF(%d) is not available in this version of the program. Please try GF(64) or GF(256)\n",GF);
        exit(EXIT_FAILURE);
    }

    nbRow = GF;
    nbCol = logGF;
    /* BINGF [GF][logGF] */
    table->BINGF =calloc((size_t)nbRow,sizeof(int *));
    //if (table->BINGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    table->BINGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
    //if (table->BINGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) table->BINGF[k] = table->BINGF[0] + k*nbCol;

    nbRow = GF;
    nbCol = GF;
    /* ADDGF [GF][GF] */
    table->ADDGF =calloc((size_t)nbRow,sizeof(int *));
    //if (table->ADDGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    table->ADDGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
    //if (table->ADDGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) table->ADDGF[k] = table->ADDGF[0] + k*nbCol;

    /*DECGF [GF] */
    table->DECGF = calloc(GF,sizeof(int));
    /*GFDEC [GF] */
    table->GFDEC = calloc(GF,sizeof(int));

    /* MULGF [GF][GF] */
    table->MULGF =calloc((size_t)nbRow,sizeof(int *));
    //if (table->MULGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    table->MULGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
    //if (table->MULGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) table->MULGF[k] = table->MULGF[0] + k*nbCol;

    /* DIVGF [GF][GF] */
    table->DIVGF =calloc((size_t)nbRow,sizeof(int *));
    //if (table->DIVGF  == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    table->DIVGF [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
    //if (table->DIVGF [0] == NULL) err(EXIT_FAILURE,"%s:%d > malloc failed !",__FILE__,__LINE__);
    for (k=1; k<nbRow; k++) table->DIVGF[k] = table->DIVGF[0] + k*nbCol;

    /* MULDEC [GF][GF] */
    table->MULDEC =calloc((size_t)nbRow,sizeof(int *));
    table->MULDEC [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
    for (k=1; k<nbRow; k++) table->MULDEC[k] = table->MULDEC[0] + k*nbCol;

    /* DIVDEC [GF][GF] */
    table->DIVDEC =calloc((size_t)nbRow,sizeof(int *));
    table->DIVDEC [0] = calloc((size_t)nbRow*nbCol,sizeof(int));
    for (k=1; k<nbRow; k++) table->DIVDEC[k] = table->DIVDEC[0] + k*nbCol;

    if(GF==16)
    {
        for(g=0; g<GF; g++)
            for(l=0; l<logGF; l++)
                table->BINGF[g][l] = BinGF_16[g][l];
        //printf("Loading of the binary image of GF(64): Success\n");
        //fflush(stdout);
    }


    if(GF==32)
    {
        for(g=0; g<GF; g++)
            for(l=0; l<logGF; l++)
                table->BINGF[g][l] = BinGF_32[g][l];
        //printf("Loading of the binary image of GF(64): Success\n");
        //fflush(stdout);
    }

    if(GF==64)
    {
        for(g=0; g<GF; g++)
            for(l=0; l<logGF; l++)
                table->BINGF[g][l] = BinGF_64[g][l];
        //printf("Loading of the binary image of GF(64): Success\n");
        //fflush(stdout);
    }

    if(GF==256)
    {
        for(g=0; g<GF; g++)
            for(l=0; l<logGF; l++)
                table->BINGF[g][l] = BinGF_256[g][l];
        //printf("Loading of the binary image of GF(256): Success\n");
        //fflush(stdout);
    }


    if(GF==4096)
    {
        for(g=0; g<GF; g++)
            for(l=0; l<logGF; l++)
                table->BINGF[g][l] = BinGF_4096[g][l];
        //printf("Loading of the binary image of GF(256): Success\n");
        //fflush(stdout);
    }


    /*
     * Build the addition, multiplication and division tables (corresponding to GF[q])
     */
   // Table_Add_GF(table,GF,logGF);
    Table_dec_GF(table, GF,logGF);
    Table_Mul_GF (table->MULGF,GF);
    Table_Div_GF (table->DIVGF,GF);
    Table_Mul_DEC(table,GF);
    Table_Div_DEC(table,GF);
    Table_Add_GF(table, GF, logGF);
}


/*!
 * \fn FreeTable(table_t *table)
 * \brief Free the memory in a table_t structure
 * Inputs        :
 * 	- table_t table : Structure that contains pointers to the tables
 * Outputs       :
 */
void FreeTable(table_t *table)
{
    free(table->ADDGF[0]);
    free(table->MULGF[0]);
    free(table->DIVGF[0]);
    free(table->BINGF[0]);
    free(table->ADDGF);
    free(table->MULGF);
    free(table->DIVGF);
    free(table->BINGF);
}




