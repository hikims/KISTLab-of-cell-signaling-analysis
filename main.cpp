#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <vector>
#include <time.h>
#include "../parameters.h"
#include "Boolean.hpp"
#include "initial_data.hpp"
#include "functions.hpp"

#define PATH "USER PATH"

//----------------------------------
// Normal distribution
//----------------------------------
double gaussianRandom(double average, double stdev)
{
    double v1, v2, s, temp;
    do
    {
        v1 = 2 * ((double)rand() / RAND_MAX) - 1;
        v2 = 2 * ((double)rand() / RAND_MAX) - 1;
        s = v1 * v1 + v2 * v2;
    } while (s >= 1 || s == 0);

    s = sqrt((-2 * log(s)) / s);

    temp = v1 * s;
    temp = (stdev * temp) + average;
    return temp;
}

int main(int argc, char **argv)
{
    double *protein = mallocSingle(Np);
    double *new_protein = mallocSingle(Np);
    double *input = mallocSingle(Np);
    double *output = mallocSingle(Np);

    double *row_protein0 = mallocSingle(iter);
    double *row_protein1 = mallocSingle(iter);
    double *row_protein2 = mallocSingle(iter);
    double *row_protein3 = mallocSingle(iter);
    double *row_protein4 = mallocSingle(iter);
    double *row_protein5 = mallocSingle(iter);
    
    double *row_protein6 = mallocSingle(iter);
    double *row_protein7 = mallocSingle(iter);
    double *row_protein8 = mallocSingle(iter);
    double *row_protein9 = mallocSingle(iter);
    double *row_protein10 = mallocSingle(iter);
    
    double *row_protein11 = mallocSingle(iter);
    double *row_protein12 = mallocSingle(iter);
    double *row_protein13 = mallocSingle(iter);
    double *row_protein14 = mallocSingle(iter);
    double *row_protein15 = mallocSingle(iter);
    
    double **erk12 = mallocDouble(simul_time, iter);
    
    double *avg = mallocSingle(simul_time);
    double *inhibi_avg = mallocSingle(simul_time);
    double *erk_amp = mallocSingle(simul_time);
    
    double **bliss_value = mallocDouble(p,p);

    int spinter_time = iter / spinterval; 

    double *sppro_avg = mallocSingle(spinter_time);
    double *sppro_std = mallocSingle(spinter_time);
    double *sppro_tm = mallocSingle(spinter_time);  

    int cvinter_time = iter / cvinterval;

    double *row_avg = mallocSingle(iter);
    double *row_std = mallocSingle(iter);
    
    double *std_avg = mallocSingle(cvinter_time);
    double *std_std = mallocSingle(cvinter_time);
    double *std_temp_std1 = mallocSingle(cvinter_time);
    double *std_temp_std2 = mallocSingle(cvinter_time);
    double *std_tm = mallocSingle(cvinter_time);  

    double *protein_avg = mallocSingle(cvinter_time);
    double *protein_std = mallocSingle(cvinter_time);
    double *protein_temp_std1 = mallocSingle(cvinter_time);
    double *protein_temp_std2 = mallocSingle(cvinter_time);
    double *protein_tm = mallocSingle(cvinter_time);  

    
    int num_pro = p + (p * (p - 1)) / 2; // combination order by chosen inhibition proteins
    double *avg_per = mallocSingle(num_pro);
    double *inhi_index = mallocSingle(num_pro);

    int inhibi_point = iter / 2; // set inhibition start point
    double pre_weight[Nw] = {
        0,
    };
    double boolean_boo[Np];

    int i, j, k, num;

    srand((unsigned)time(NULL) + (unsigned)getpid()); // random wiener value for each simulation

    char filename[1024];

    //double ku = atof(argv[1]);
    //printf("ku : %.2lf\n", ku);
    double ku = 4.0;
    //double kd = atof(argv[2]);
    double kd = 3.5;
    //double sigma = atof(argv[3]);
    double sigma = 0.01;

    double phi = -1.0; // Negative feedback index

    double erk_mean = 0.978204; // ERK1/2 amplitude of non-stochastic case.    

for(int z = 0; z < simul_time; z++)
{
    //printf("%dth simulation\n", z);

    // Initial data
    initialdata(protein, Np, pre_weight, Nw);
    
    // Fixed weights
    FILE *fwei = fopen(PATH "./code/weight.txt", "r");
    if (!fwei)
    {
        printf("stochastic: Failed to open!\n");
    }
    for(int k = 0; k < Nw; k++)
    {
        fscanf(fwei, "%lf\n", &pre_weight[k]);
    }
    fclose(fwei);

    // Numerical Calculation
    for (k = 0; k < iter; k++)
    {
            for (j = 0; j < Np; j++)
            {
                boolean_boo[j] = 0.0;
            }

        if (k % cut == 0)
        {
            printf("calc time: %.3lf\n", k * dt);
        }

                for (i = 0; i < Np; i++)
                {
                    input[i] = 0.0;
                    output[i] = 0.0;

                    if (i == 0)
                    {
                        row_protein0[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        first_Boolean_gate(&(boolean_boo[0]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, 0.0, 1, pre_weight[0]);
                    }
                    if (i == 1)
                    {
                        row_protein1[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        first_Boolean_gate(&(boolean_boo[1]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, 0.0, 1, pre_weight[1]);
                    }
                    if (i == 2)
                    {   
                        row_protein2[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        first_Boolean_gate(&(boolean_boo[0]), &(output[0]), &(new_protein[0]), protein[0], ku, kd, sigma, wiener, 0.0, 1, pre_weight[0]);
                        input[i] = get_input(1, protein[0], boolean_boo[0]);
                        first_Boolean_gate(&(boolean_boo[2]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[2]);
                    }
                    if (i == 3)
                    {                  
                        row_protein3[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        first_Boolean_gate(&(boolean_boo[1]), &(output[1]), &(new_protein[1]), protein[1], ku, kd, sigma, wiener, 0.0, 1, pre_weight[1]);
                        input[i] = get_input(1, protein[1], boolean_boo[1]);
                        first_Boolean_gate(&(boolean_boo[3]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 3, pre_weight[3], pre_weight[4], pre_weight[5]);
                    }
                    if (i == 4)
                    {
                        row_protein4[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[2] = get_input(1, protein[0], boolean_boo[0]);
                        first_Boolean_gate(&(boolean_boo[2]), &(output[2]), &(new_protein[2]), protein[2], ku, kd, sigma, wiener, input[2], 1, pre_weight[2]);
                        input[15] = get_input(1, protein[12], boolean_boo[12]);
                        first_Boolean_gate(&(boolean_boo[15]), &(output[15]), &(new_protein[15]), protein[15], ku, kd, sigma, wiener, input[15], 2, pre_weight[20], pre_weight[21]);
                        input[i] = get_input(2, protein[2], boolean_boo[2], protein[15], phi * boolean_boo[15]);
                        first_Boolean_gate(&(boolean_boo[4]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[6]);

                        // Non-negative feedback
                        /*row_protein4[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[2] = get_input(1, protein[0], boolean_boo[0]);
                        first_Boolean_gate(&(boolean_boo[2]), &(output[2]), &(new_protein[2]), protein[2], ku, sigma, wiener, input[2], 1, pre_weight[2]);
                        input[i] = get_input(1, protein[2], boolean_boo[2]);
                        first_Boolean_gate(&(boolean_boo[4]), &(output[i]), &(new_protein[i]), protein[i], ku, sigma, wiener, input[i], 1, pre_weight[6]);*/

                    }
                    if (i == 5)
                    {
                        row_protein5[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[3] = get_input(1, protein[1], boolean_boo[1]);
                        first_Boolean_gate(&(boolean_boo[3]), &(output[3]), &(new_protein[3]), protein[3], ku, kd, sigma, wiener, input[3], 3, pre_weight[3], pre_weight[4], pre_weight[5]);
                        input[i] = get_input(1, protein[3], boolean_boo[3]);
                        first_Boolean_gate(&(boolean_boo[5]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 2, pre_weight[9], pre_weight[10]);
                    }
                    if (i == 6)
                    {    
                        row_protein6[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[3] = get_input(1, protein[1], boolean_boo[1]);
                        first_Boolean_gate(&(boolean_boo[3]), &(output[3]), &(new_protein[3]), protein[3], ku, kd, sigma, wiener, input[3], 3, pre_weight[3], pre_weight[4], pre_weight[5]);
                        input[4] = get_input(2, protein[2], boolean_boo[2], protein[21], boolean_boo[21]);
                        first_Boolean_gate(&(boolean_boo[4]), &(output[4]), &(new_protein[4]), protein[4], ku, kd, sigma, wiener, input[4], 1, pre_weight[6]);
                        input[i] = get_input(2, protein[3], boolean_boo[3], protein[4], boolean_boo[4]);
                        first_Boolean_gate(&(boolean_boo[6]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[7]);
                    }
                    if (i == 7)
                    {    
                        row_protein7[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[3] = get_input(1, protein[1], boolean_boo[1]);
                        first_Boolean_gate(&(boolean_boo[3]), &(output[3]), &(new_protein[3]), protein[3], ku, kd, sigma, wiener, input[3], 3, pre_weight[3], pre_weight[4], pre_weight[5]);
                        input[i] = get_input(1, protein[3], boolean_boo[3]);
                        first_Boolean_gate(&(boolean_boo[7]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[8]);
                    }
                    if (i == 8)
                    {
                        row_protein8[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[5] = get_input(1, protein[3], boolean_boo[3]);
                        first_Boolean_gate(&(boolean_boo[5]), &(output[5]), &(new_protein[5]), protein[5], ku, kd, sigma, wiener, input[5], 2, pre_weight[9], pre_weight[10]);
                        input[i] = get_input(1, protein[5], boolean_boo[5]);
                        first_Boolean_gate(&(boolean_boo[8]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[11]);
                    }
                    if (i == 9)
                    {
                        row_protein9[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[5] = get_input(1, protein[3], boolean_boo[3]);
                        first_Boolean_gate(&(boolean_boo[5]), &(output[5]), &(new_protein[5]), protein[5], ku, kd, sigma, wiener, input[5], 2, pre_weight[9], pre_weight[10]);
                        input[8] = get_input(1, protein[5], boolean_boo[5]);
                        first_Boolean_gate(&(boolean_boo[8]), &(output[8]), &(new_protein[8]), protein[8], ku, kd, sigma, wiener, input[8], 1, pre_weight[11]);
                        input[i] = get_input(2, protein[5], boolean_boo[5], protein[8], boolean_boo[8]);
                        first_Boolean_gate(&(boolean_boo[9]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 3, pre_weight[13], pre_weight[14], pre_weight[15]);
                    }
                    if (i == 10)
                    {               
                        row_protein10[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[6] = get_input(2, protein[3], boolean_boo[3], protein[4], boolean_boo[4]);
                        first_Boolean_gate(&(boolean_boo[6]), &(output[6]), &(new_protein[6]), protein[6], ku, kd, sigma, wiener, input[6], 1, pre_weight[7]);
                        input[9] = get_input(2, protein[5], boolean_boo[5], protein[8], boolean_boo[8]);
                        first_Boolean_gate(&(boolean_boo[9]), &(output[9]), &(new_protein[9]), protein[9], ku, kd, sigma, wiener, input[9], 3, pre_weight[13], pre_weight[14], pre_weight[15]);
                        input[i] = get_input(2, protein[6], boolean_boo[6], protein[9], boolean_boo[9]);
                        first_Boolean_gate(&(boolean_boo[10]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[12]);
                    }
                    if (i == 11)
                    {
                        row_protein11[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[9] = get_input(2, protein[5], boolean_boo[5], protein[8], boolean_boo[8]);
                        first_Boolean_gate(&(boolean_boo[9]), &(output[9]), &(new_protein[9]), protein[9], ku, kd, sigma, wiener, input[9], 3, pre_weight[13], pre_weight[14], pre_weight[15]);
                        input[i] = get_input(1, protein[9], boolean_boo[9]);
                        first_Boolean_gate(&(boolean_boo[11]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[16]);
                    }
                    if (i == 12)
                    {                     
                        row_protein12[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[10] = get_input(2, protein[6], boolean_boo[6], protein[9], boolean_boo[9]);
                        first_Boolean_gate(&(boolean_boo[10]), &(output[10]), &(new_protein[10]), protein[10], ku, kd, sigma, wiener, input[10], 1, pre_weight[12]);
                        input[i] = get_input(1, protein[10], boolean_boo[10]);
                        first_Boolean_gate(&(boolean_boo[12]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[17]);
                    }
                    if (i == 13)
                    {                     
                        row_protein13[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[7] = get_input(1, protein[3], boolean_boo[3]);
                        first_Boolean_gate(&(boolean_boo[7]), &(output[7]), &(new_protein[7]), protein[7], ku, kd, sigma, wiener, input[7], 1, pre_weight[8]);
                        input[i] = get_input(1, protein[7], boolean_boo[7]);
                        first_Boolean_gate(&(boolean_boo[13]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[18]);
                    }
                    if (i == 14)
                    {                   
                        row_protein14[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[9] = get_input(2, protein[5], boolean_boo[5], protein[8], boolean_boo[8]);
                        first_Boolean_gate(&(boolean_boo[9]), &(output[9]), &(new_protein[9]), protein[9], ku, kd, sigma, wiener, input[9], 3, pre_weight[13], pre_weight[14], pre_weight[15]);
                        input[11] = get_input(1, protein[9], boolean_boo[9]);
                        first_Boolean_gate(&(boolean_boo[11]), &(output[11]), &(new_protein[11]), protein[11], ku, kd, sigma, wiener, input[11], 1, pre_weight[16]);
                        input[i] = get_input(2, protein[9], boolean_boo[9], protein[11], phi * boolean_boo[11]);
                        first_Boolean_gate(&(boolean_boo[14]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[19]);
                    }
                    if (i == 15)
                    {                     
                        row_protein15[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[12] = get_input(1, protein[10], boolean_boo[10]);
                        first_Boolean_gate(&(boolean_boo[12]), &(output[12]), &(new_protein[12]), protein[12], ku, kd, sigma, wiener, input[12], 1, pre_weight[17]);
                        input[i] = get_input(1, protein[12], boolean_boo[12]);
                        first_Boolean_gate(&(boolean_boo[15]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 2, pre_weight[20], pre_weight[21]);

                        // Non-negative feedback
                        /*row_protein15[k] = protein[i];
                        double wiener = gaussianRandom(0, dt);
                        input[12] = get_input(1, protein[10], boolean_boo[10]);
                        first_Boolean_gate(&(boolean_boo[12]), &(output[12]), &(new_protein[12]), protein[12], ku, kd, sigma, wiener, input[12], 1, pre_weight[17]);
                        input[i] = get_input(1, protein[12], boolean_boo[12]);
                        first_Boolean_gate(&(boolean_boo[15]), &(output[i]), &(new_protein[i]), protein[i], ku, kd, sigma, wiener, input[i], 1, pre_weight[21]);*/
                    }                   
                }
            for(i = 0; i < Np; i++)
            {
                protein[i] = new_protein[i];
                
                // control stimuli  
                protein[0] = 1.0;
                protein[1] = 1.0;
                
                // choose inhibitor 
                if(k >= inhibi_point)
                {
                    //num = 9; // order of inhibition cases	
                    //protein[4] = 0.0;
                    //protein[6] = 0.0;
                    //protein[10] = 0.0;
                    //protein[12] = 0.0;
                    //protein[3] = 0.0;
                    //protein[5] = 0.0;
                    //protein[9] = 0.0;
                }
            }
            
       // multiple simulation data 
       erk12[z][k] = row_protein15[k];
    }
        
    //  amplitude   
    /*int st_point = iter / 2;
    double erk_max = row_protein15[st_point];
    double erk_min = row_protein15[st_point];
    double erk_amplitude = amplitude(row_protein15, erk_max, erk_min, st_point, iter);     
    erk_amp[z] = erk_amplitude;

    double sos_max = row_protein4[st_point];
    double sos_min = row_protein4[st_point];
    double sos_amplitude = amplitude(row_protein4, sos_max, sos_min, st_point, iter);     

    FILE *famp = fopen(PATH "./data/amplitude.txt", "a");
    fprintf(famp, "%.2lf %lf\n", ku, sos_amplitude);
    fclose(famp);*/

    /*if(z == 0) // amplitude for chosen simulation
    {    
        FILE *fprotein = fopen(PATH "./data/trajectory.txt", "w");
        for(k = 0; k < iter; k++)
        {
        fprintf(fprotein, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", k * dt, row_protein0[k], row_protein1[k], row_protein2[k], row_protein3[k], row_protein4[k], 
        row_protein5[k], row_protein6[k], row_protein7[k], row_protein8[k], row_protein9[k], row_protein10[k], row_protein11[k], row_protein12[k], row_protein13[k], row_protein14[k], row_protein15[k]);
        }
        fclose(fprotein);
    }*/
    
    // Smooth spline    
    /*if(z == 0) // smooth spline for chosen simulation
    {
        spline(row_protein15, sppro_avg, sppro_std, sppro_tm, spinter_time, spinterval);
        FILE *fspline;
        fspline = fopen(PATH "./data/spline.txt", "w");
        for(int k = 0; k < spinter_time; k++)
        {
            fprintf(fspline, "%.2lf %lf %lf\n", sppro_tm[k], sppro_avg[k], sppro_std[k]);
        }
        fclose(fspline);
    }*/

    // inhibition average 	
    /*avg[z] = avg_function(row_protein15, inhibi_point);
    inhibi_avg[z] = inhi_avg_function(row_protein15, inhibi_point, iter);*/
}    

    // Interval coefficient variation 
    /*data_sort(erk12, row_avg, row_std, simul_time, iter);
    interCV(row_std, std_avg, std_std, std_temp_std1, std_temp_std2, std_tm, cvinter_time, cvinterval);
    interCV(row_avg, protein_avg, protein_std, protein_temp_std1, protein_temp_std2, protein_tm, cvinter_time, cvinterval);
    FILE *fintercv;
    fintercv = fopen(PATH "./data/interCV.txt", "w");
    for(int k = 0; k < cvinter_time; k++)
    {
        fprintf(fintercv, "%.2lf %lf %lf\n", protein_tm[k], protein_avg[k], std_avg[k]);
    }
    fclose(fintercv);*/

	
    // Coefficient variation 
    /*double erk_avg;
    double erk_std;
    const char* pro_name = "ERK1/2";
    CV(erk12, &erk_avg, &erk_std, simul_time, iter);     
    FILE *ferk;
    ferk = fopen(PATH "./data/erk12.txt", "w");
    fprintf(ferk, "%s %lf %lf\n", pro_name, erk_avg, erk_std);
    fclose(ferk);*/
        
    
    // RMSE 
    /*double erk_rmse;
    erk_rmse = RMSE(erk_amp, erk_mean, simul_time);
    FILE *frmse = fopen(PATH "./data/rmse.txt", "a");
    fprintf(frmse, "%.2lf %.2lf %lf\n", kd, sigma, erk_rmse);
    fclose(frmse);*/
    
    
    // bliss index
    /*double avg_reduc = 0.0;
    avg_reduc = reduction(avg, inhibi_avg, simul_time);
    FILE *freduc = fopen(PATH "./data/reduction_rate.txt", "a");
    fprintf(freduc, "%d %lf\n", num, avg_reduc);
    fclose(freduc);

    if (num == (num_pro - 1))
    {
        //Reduction rate data
        FILE *freducdata = fopen(PATH "./data/reduction_rate.txt", "r");
        if(!freducdata)
        {
            printf("Reduction rate : Failed to open!\n");
        }
        for(int d = 0; d < num_pro; d++)
        {
            fscanf(freducdata, "%lf %lf\n", &inhi_index[d], &avg_per[d]);
        }
        fclose(freducdata);

    	bliss(bliss_value, avg_per, p);
        FILE *fbliss = fopen(PATH "./data/bliss.txt", "w");
        for(int i = 0; i < p; i++)
        {
            fprintf(fbliss, "\n");
            for(int j = 0; j < p; j++)
            {
                fprintf(fbliss, "%lf\t", bliss_value[i][j]);
            }
            fprintf(fbliss, "\n");
        }
        fclose(fbliss);
    }*/    
    
        
    freeMallocSingle(protein);
    freeMallocSingle(new_protein);
    freeMallocSingle(input);
    freeMallocSingle(output);

    freeMallocSingle(row_protein0);
    freeMallocSingle(row_protein1);
    freeMallocSingle(row_protein2);
    freeMallocSingle(row_protein3);
    freeMallocSingle(row_protein4);
    freeMallocSingle(row_protein5);
    freeMallocSingle(row_protein6);

    freeMallocSingle(row_protein7);
    freeMallocSingle(row_protein8);
    freeMallocSingle(row_protein9);
    freeMallocSingle(row_protein10);
    freeMallocSingle(row_protein11);
    freeMallocSingle(row_protein12);
    freeMallocSingle(row_protein13);
    freeMallocSingle(row_protein14);
    freeMallocSingle(row_protein15);
    
    freeMallocDouble(bliss_value, p);
    freeMallocDouble(erk12, simul_time);
    
    freeMallocSingle(sppro_avg);
    freeMallocSingle(sppro_std);
    freeMallocSingle(sppro_tm);  

    freeMallocSingle(row_avg);
    freeMallocSingle(row_std);  
  
    freeMallocSingle(std_avg);
    freeMallocSingle(std_std);
    freeMallocSingle(std_temp_std1);
    freeMallocSingle(std_temp_std2);
    freeMallocSingle(std_tm);  
    
    freeMallocSingle(protein_avg);
    freeMallocSingle(protein_std);
    freeMallocSingle(protein_temp_std1);
    freeMallocSingle(protein_temp_std2);
    freeMallocSingle(protein_tm);  
    
    freeMallocSingle(avg);
    freeMallocSingle(inhibi_avg);
    freeMallocSingle(erk_amp);
    freeMallocSingle(avg_per);
    freeMallocSingle(inhi_index);


    printf("\n");
    printf("Complete Simulations!!\n");
    printf("\n");

    return 0;
}
