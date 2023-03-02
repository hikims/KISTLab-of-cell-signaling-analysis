#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <vector>
#include <time.h>
#include "../parameters.h"

//----------------------------------
// Amplitude
//----------------------------------
double amplitude(double *data, double max, double min, int startsize, int endsize)
{		
    double temp = 0.0;
    for(int i = startsize; i < endsize; i++)
    {
        if(data[i] > max)
        {
            max = data[i];
        }
        if(data[i] < min)
        {
            min = data[i];
        }
    } 
    
    temp = fabs(max - min);
    return temp;
}

//----------------------------------
// Root Mean Squar Error
//----------------------------------
double RMSE(double *data1, double mean, int size)
{
    double sum = 0.0;

    for(int i =0; i < size; i++)
    {
        sum += pow((data1[i] - mean), 2);
    }
    double temp = sum / simul_time;
    double rmse = sqrt(temp);
    return rmse;
}

//----------------------------------
// Inhibition Average
//----------------------------------
double avg_function(double *data, int inhi_size)
{
    double avg_temp = 0.0;
    
    for(int j = 0; j < inhi_size; j++)
    {
        avg_temp += data[j];     
    }
    double avg = avg_temp / inhi_size;
    return avg;
}

double inhi_avg_function(double *data, int inhi_size, int size)
{
    double inhi_avg_temp = 0.0;

    for(int j = inhi_size; j < size; j++)
    {
        inhi_avg_temp += data[j];
    } 
    double inhi_avg = inhi_avg_temp / inhi_size;
    return inhi_avg;
}

//----------------------------------
// Mean
//----------------------------------
double MEAN(double *data, int size)
{
    double avg;
    double sum = 0.0;
    for(int i = 0; i < size; i++)
    {
        sum += data[i];
    }
    avg = sum / size;
    return avg;
}

//----------------------------------
// Standard Variation
//----------------------------------
double VAR(double *data, double mean)
{
    double var, std;
    double v_sum = 0.0;
    for(int j = 0; j < simul_time; j++)
    {
        v_sum += pow( (data[j] - mean), 2);
    }
    var = v_sum / simul_time;
    std = sqrt(var);
    return std;
}

//----------------------------------
// Reduction rate
//----------------------------------
double reduction(double *data1, double *data2, int size)
{    
    double diff_avg[size];
    for(int i = 0; i < size; i++)
    {
        diff_avg[i] = 0.0;
        diff_avg[i] = data1[i] - data2[i]; 
    }
    double avg1 = MEAN(data1, size);
    double avg2 = MEAN(data2, size);
    double temp = 1.0 - (avg2 / avg1);
    
    return temp;
}    
    
//----------------------------------
// Bliss index
//----------------------------------    
double bliss(double **data2, double *data1, int pro_size)
{       
    for(int i = 0; i < pro_size; i++)
    {
        for(int j = 0; j < pro_size; j++)
        {
            data2[i][j] = 0.0;
           if(j > i)
           {
              int z = j + 6 * i - i * (i - 1) * 0.5;
              double predicted_value = data1[i] + data1[j] - (data1[i] * data1[j]);
              data2[i][j] = predicted_value / data1[z];
           }         
        }
    }
    return 0;
}

//----------------------------------
// Smooth spline
//----------------------------------       
double spline(double *data1, double *avg, double *std, double *tm, int size, int interval)
{    
    for(int z = 0; z < size; z++)
    {
        int st_point = z * interval;
        int end_point = (z + 1) * interval;

        for(int i = st_point; i < end_point; i++)
        {
            avg[z] += data1[i] / interval;
        } 

        double temp_std1 = 0.0;
        double temp_std2 = 0.0;
        for(int i = st_point; i < end_point; i++)
        {   
            temp_std1 = data1[i] - avg[z];
            temp_std2 += pow(temp_std1, 2) / interval;
        }
        std[z] = sqrt(temp_std2);
        tm[z] = st_point;
    }
    return 0;
}   

//----------------------------------
// Coefficient variation
//----------------------------------    
double CV(double **data, double* final_avg, double* std, int size, int col)
{        
    double avg[size];
    
    for(int i = 0; i < size; i++)
    {
        avg[i] = 0.0;
        for (int j = 0; j < col; j++)
        {
            avg[i] += data[i][j] / col;
        }
    }

    *final_avg = 0.0;
    for(int i = 0; i < size; i++)
    {
        *final_avg += avg[i] / size;
    }

    *std = 0.0;
    double temp_std1 = 0.0;
    double temp_std2 = 0.0;
    for (int i = 0; i < size; i++)
    {
        temp_std1 = avg[i] - *final_avg;
        temp_std2 += pow(temp_std1, 2) / size;
    }
    *std = sqrt(temp_std2) / *final_avg;

    return 0;
}
    
//----------------------------------
// Interval coefficient variation
//----------------------------------    
double data_sort(double **data4, double* data1, double* data2, int size, int row_size)
{
    for(int j = 0; j < row_size; j++)
    {
        data1[j] = 0.0;
        data2[j] = 0.0;
        double row_temp_std1 = 0.0;
        double row_temp_std2 = 0.0;
        for (int i = 0; i < size; i++)
        {
            data1[j] += data4[i][j] / size;
        }
        for (int i = 0; i < size; i++)
        {
            row_temp_std1 = data4[i][j] - data1[j];
            row_temp_std2 += pow(row_temp_std1, 2) / size;
        }
        data2[j] = sqrt(row_temp_std2);
    }
    return 0;
}    

double interCV(double *data1, double *avg, double *std, double *temp_std1, double *temp_std2, double *tm, int size, int interval)
{    
    for(int z = 0; z < size; z++)
    {
        int st_point = z * interval;
        int end_point = (z + 1) * interval;

        for(int i = st_point; i < end_point; i++)
        {
            avg[z] += data1[i] / interval;
        } 

        for(int i = st_point; i < end_point; i++)
        {   
            temp_std1[z] = data1[i] - avg[z];
            temp_std2[z] += pow(temp_std1[z], 2) / interval;
        }
        std[z] = sqrt(temp_std2[z]) / avg[z];
        tm[z] = st_point;
    }
    return 0;
}   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
