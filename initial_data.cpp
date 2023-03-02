#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <vector>
#include <time.h>

#include "../parameters.h"


//----------------------------------
// FUNCTION initialdata
//----------------------------------
void initialdata(double *protein, int protein_col, double *pre_weight, int pre_weight_num)
{
    /* protein */
    for (int j = 0; j < protein_col; j++)
    {
        if (j < 2)
        {
            protein[j] = 1.0;
        }
        else
        {
            protein[j] = 0.0;
        }
    }

    /* weight */
    for (int i = 0; i < pre_weight_num; i++)
    {  
        pre_weight[i] = ((double)rand() / RAND_MAX);
    }
}

//----------------------------------
// Allocation
//----------------------------------
double** mallocDouble(int row, int col)
{
    double **d;

    d = (double**)malloc(sizeof(double)*row);
    for(int i = 0; i < row; i++)
    {
        d[i] = (double*)malloc(sizeof(double*)*col);
        for(int j = 0; j < col; j++)
        {
            d[i][j] = 0;
        }
    }

    return d;
}

double* mallocSingle(int col)
{
    double *dd;

    dd = (double*)malloc(sizeof(double)*col);
    for(int j = 0; j < col; j++)
    {
        dd[j] = 0;
    }

    return dd;
}


void freeMallocDouble(double **data, int row)
{
    int i = 0;

    for(i = 0; i < row; i++)
    free(data[i]);
    free(data);
}

void freeMallocSingle(double *data)
{
    free(data);
}
