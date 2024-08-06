# KISTLab-of-cell-signaling-analysis
README

These files are open-source for computing cell signaling pathways and are made up of C/C++.

Compiling

	In Linux and macOS, the exec.sh can compile the main.cpp executables, 
 	as well as Boolean.cpp, functions.cpp, and initial_data.cpp with

	chmod +x exec.sh

	./exec.sh

According to the user computer system, -std=c++11 in exec.sh can be changeable.
If the user wants to calculate variation by a specific parameter, 
use cat filename.txt | xargs -L 1 ./a.out in exec.sh. (Some examples are written in exec.sh.)


Parameters

	Np: the number of proteins in the signaling pathway
	Nw: the number of weights in the signaling pathway
	rb: the parameter of the Waller-Kraft operator
	p: the number of inhibitors
	taneps: the parameter of the hyper tangent function
	dt: time shift
	iter: calculation time
	cut: a specific interval of simulation time
	simul_time: simulation time
	cvinterval: coefficient variation period on simulation time
	spinterval: smooth spline period on calculation time 

 	In this study, we set arbitrary parameters (parameters.h) as follow:
    		
      	#define Np (16)
		#define Nw (22)
		#define rb (0.75)
		#define p (4)
		#define taneps  (1.0)
		#define dt (1.0)
		#define iter (100000)
		#define cut (1)
		#define simul_time  (1)
		#define cvinterval (1000)
		#define spinterval (1000)

The upstream (k_u), downstream (k_d), and stochastic rate (σ) imply the control parameters of the protein’s activation, and they are indicated in the main script (main.cpp).


Boolean operator

	When the user calculates the Boolean value of each protein, use the get_min_max function in Boolean.cpp. 

	In this study, our open source considers a maximum of four weights for each protein, 
 	and they have to be inserted sequentially. 

 		double get_input(int _size, double p_num1, double b_num1, double p_num2, double b_num2, double p_num3, 
  		double b_num3, double p_num4, double b_num4){
    
    		double pro[4] = {p_num1, p_num2, p_num3, p_num4};
   		 	double wei[4] = {b_num1, b_num2, b_num3, b_num4};
   	 		double temp = 0.0;
    		for (int i = 0; i < _size; i++){temp += pro[i] * wei[i];}
    	return temp;}

    For the upstream value, the user can use the get_input function and the protein and 
    its boolean value has to be inserted crosswise. 
    For example, if one protein is activated by two proteins, the upstream value is calculated:

		input = get_input(2, protein1, boolean_value1, protein2, boolean_value2)

	Finally, using the Boolean value and input, the activation of protein can be computed by first_Boolean_gate, T, 
 	and protein_calc functions.

		void first_Boolean_gate(double *boolean, double *output, double *new_protein, double protein, double para1, 
  		double para2, double para3, double z, double input, int size, double out_w1, double out_w2, double out_w3, 
  		double out_w4, double out_w5){
   
   			double min = 0.0, max = 0.0;
			get_min_max(size, &min, &max, out_w1, out_w2, out_w3, out_w4, out_w5);
 			# get_min_max function is to simply find the maximum and minimim of weights (in Boolean.cpp).

    		double boo = rb * min + (1.0 - rb) * max;
    		*boolean = boo;
    		*output = boo * protein;
    		*new_protein = protein_calc(protein, input, *output, z, para1, para2, para3);}

		double T(double x){return tanh( x / taneps);}

		double protein_calc(double a, double b, double c, double w, double para1, double para2, double para3){
    		return a + T(para1 * b - para2 * c * a) * dt + a * para3 * w;
      		# This function is based on Euler-Maruyama numerical method.}



Calculation Functions

Single simulation: Amplitude and smooth spline of each protein are calculated via a single simulation (simul_time=1). 
	
	Amplitude is estimated as the difference between the maximum and minimum of protein’s activation 
 	from a specific time point (startsize, endsize) without stochastic effects (amplitude in functions.cpp).
  
		double amplitude(double *data, double max, double min, int startsize, int endsize){		
    		
      		double temp = 0.0;
    		for(int i = startsize; i < endsize; i++){
        		if(data[i] > max){max = data[i];}
        		if(data[i] < min){min = data[i];}} 
    		temp = fabs(max - min);
   		return temp;}

	Smooth spline is estimated at a constant period (splinterval) without stochastic effects (spline in functions.cpp). 

 		double spline(double *data1, double *avg, double *std, double *tm, int size, int interval){    
    		for(int z = 0; z < size; z++){
        		int st_point = z * interval;
        		int end_point = (z + 1) * interval;
        			for(int i = st_point; i < end_point; i++){avg[z] += data1[i] / interval;} 

        		double temp_std1 = 0.0;
        		double temp_std2 = 0.0;
        			for(int i = st_point; i < end_point; i++){   
            		temp_std1 = data1[i] - avg[z];
            		temp_std2 += pow(temp_std1, 2) / interval;}
       		std[z] = sqrt(temp_std2);
        		tm[z] = st_point;}
    	return 0;}   


Multiple simulations: Root mean square error, coefficient variation, and bliss index of each protein are calculated via multiple simulations (simul_time=1000). 

	The root mean square error (RMSE in functions.cpp) is to find suitable parameters set based on a specific parameter. 

 		double RMSE(double *data1, double mean, int size){
    
    		double sum = 0.0;
			for(int i =0; i < size; i++){sum += pow((data1[i] - mean), 2);}
    		double temp = sum / simul_time;
    		double rmse = sqrt(temp);
    	return rmse;}

  
	Coefficient variation (CV in functions.cpp) is estimated as the mean (MEAN in functions.cpp) and 
 	standard deviation (VAR in functions.cpp) of multiple simulations with stochastic effects. 

		double CV(double **data, double* final_avg, double* std, int size, int col){        
   			double avg[size];
    		for(int i = 0; i < size; i++){
        		avg[i] = 0.0;
        			for (int j = 0; j < col; j++){avg[i] += data[i][j] / col;}}
   			*final_avg = 0.0;
    			for(int i = 0; i < size; i++){*final_avg += avg[i] / size;}
    		*std = 0.0;
    		double temp_std1 = 0.0;
    		double temp_std2 = 0.0;
    		for (int i = 0; i < size; i++){
        		temp_std1 = avg[i] - *final_avg;
        		temp_std2 += pow(temp_std1, 2) / size;}
    		*std = sqrt(temp_std2) / *final_avg;
    	return 0;}
  
	For interval coefficient variation (interCV in functions.cpp) is calculated through the data array (data_sort in functions.cpp)  
 	with a constant period (cvlinterval). We set period as 1000 in this study.

  		double interCV(double *data1, double *avg, double *std, double *temp_std1, double *temp_std2, 
   		double *tm, int size, int interval){    
    		for(int z = 0; z < size; z++){
        		int st_point = z * interval;
        		int end_point = (z + 1) * interval;
        			for(int i = st_point; i < end_point; i++){avg[z] += data1[i] / interval;} 
       				for(int i = st_point; i < end_point; i++){   
            			temp_std1[z] = data1[i] - avg[z];
            			temp_std2[z] += pow(temp_std1[z], 2) / interval;}
        		std[z] = sqrt(temp_std2[z]) / avg[z];
        		tm[z] = st_point;}
    	return 0;}   

   	The reduction rate (reduction in functions.cpp) is calculated using the estimated mean (MEAN in functions.cpp) and 
    standard deviation (VAR in functions.cpp) and makes a file for each inhibitor. 
    The user has to be careful about the order of inhibitors. 
    This open-source considers the order from single to combination inhibitor. 
	After making a file of inhibition, the bliss index can be calculated using the file (bliss in functions.cpp).

       	double reduction(double *data1, double *data2, int size){    
    		double diff_avg[size];
    		for(int i = 0; i < size; i++){
        		diff_avg[i] = 0.0;
        		diff_avg[i] = data1[i] - data2[i];}
    		double avg1 = MEAN(data1, size);
   			double avg2 = MEAN(data2, size);
    		double temp = 1.0 - (avg2 / avg1);
   	 	return temp;} 

	Since the bliss index (bliss in functions.cpp) is calculated via inhibitors, the user has to choose the number of inhibitors (p). 
 	A protein is set as zero at a specific simulation time and estimates the mean (MEAN in functions.cpp) and 
  	standard deviation (VAR in functions.cpp) for before and after inhibition while multiple simulations. 

   		double bliss(double **data2, double *data1, int pro_size){       
    		for(int i = 0; i < pro_size; i++){
        			for(int j = 0; j < pro_size; j++){
            		data2[i][j] = 0.0;
          			if(j > i){
             	 	int z = j + 6 * i - i * (i - 1) * 0.5;
              		double predicted_value = data1[i] + data1[j] - (data1[i] * data1[j]);
             		data2[i][j] = predicted_value / data1[z];}}}
    	return 0;}
   


Main code

	Protein has to be set individually because of different weights, upstream, and downstream. 

	The proteins including a negative feedback loop indicate the phi (=-1) index in main.cpp. 

	The inhibi_point in main.cpp can be changeable by the user (It is set as half of the calculation time in open-source). 

	The num implies the order of inhibition for each inhibitor. 

	Each inhibition has to be compiled individually to make a file. 

	When the user wants to do multiple stimulations, set 1 to the receptor protein, in contrast, set zero to the receptor protein for single stimulation.

	Non-negative feedback in main.cpp represents the case without a negative feedback loop. 
