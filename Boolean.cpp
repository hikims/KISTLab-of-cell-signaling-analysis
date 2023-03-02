#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <vector>
#include <time.h>
#include <unistd.h>

#include "Boolean.hpp"
#include "../parameters.h"

//------------------------
// Hyper Tangent function
//------------------------
double T(double x)
{
    return tanh( x / taneps);
}

//---------------------------------
// Euler-Maruyama numerical method
//---------------------------------
double protein_calc(double a, double b, double c, double w, double para1, double para2, double para3)
{
    return a + T(para1 * b - para2 * c * a) * dt + a * para3 * w;
}

//----------------------------------------------
// Algorithm for maximum and minimum of weights
//----------------------------------------------
void get_min_max(int _size, double *_min, double *_max, double _num1, double _num2, double _num3, double _num4, double _num5)
{
    double w[5] = {_num1, _num2, _num3, _num4, _num5};

    *_min = *_max = _num1;

    for (int i = 0; i < _size; i++)
    {
        if (w[i] < *_min)
        {
            *_min = w[i];
        }
        if (*_max < w[i])
        {
            *_max = w[i];
        }
    }
}

//--------------------
// Upstream algorithm
//--------------------
double get_input(int _size, double p_num1, double b_num1, double p_num2, double b_num2, double p_num3, double b_num3, double p_num4, double b_num4)
{
    double pro[4] = {p_num1, p_num2, p_num3, p_num4};
    double wei[4] = {b_num1, b_num2, b_num3, b_num4};
    double temp = 0.0;

    for (int i = 0; i < _size; i++)
    {
        temp += pro[i] * wei[i];
    }
    return temp;
}

//-------------------
// Boolean algorithm
//-------------------
void first_Boolean_gate(double *boolean, double *output, double *new_protein, double protein, double para1, double para2, double para3, double z, 
double input, int size, double out_w1, double out_w2, double out_w3, double out_w4, double out_w5)
{
    double min = 0.0, max = 0.0;

    get_min_max(size, &min, &max, out_w1, out_w2, out_w3, out_w4, out_w5);

    double boo = rb * min + (1.0 - rb) * max;
    *boolean = boo;
    *output = boo * protein;
    *new_protein = protein_calc(protein, input, *output, z, para1, para2, para3);
}

