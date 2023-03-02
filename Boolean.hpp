#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <vector>
#include <time.h>
#include <unistd.h>

#include "../parameters.h"

//--------------------------------------------
// FUNCTION PROTOTYPES
//--------------------------------------------

//double protein_calc(double a, double b, double c, double w);

//double get_boolean(int size, double out_w1 = 0.0, double out_w2 = 0.0, double out_w3 = 0.0, double out_w4 = 0.0, double out_w5 = 0.0);

//void first_Boolean_gate(double boo, double *output, double *new_protein, double protein, double z, double input);

void get_min_max(int _size, double *_min, double *_max, double _num1 = 0.0, double _num2 = 0.0, double _num3 = 0.0, double _num4 = 0.0, double _num5 = 0.0);

void first_Boolean_gate(double *boolean, double *output, double *new_protein, double protein, double para1, double para2, double para3, double z, double input, int size,
						double out_w1 = 0.0, double out_w2 = 0.0, double out_w3 = 0.0, double out_w4 = 0.0, double out_w5 = 0.0);

double get_input(int _size, double p_num1 = 0.0, double b_num1 = 0.0, double p_num2 = 0.0, double b_num2 = 0.0,
				 double p_num3 = 0.0, double b_num3 = 0.0, double p_num4 = 0.0, double b_num4 = 0.0);