//----------------------------------
// FUNCTION PROTPTYPES
//----------------------------------
double amplitude(double *data, double max, double min, int startsize, int endsize);
double RMSE(double *data1, double mean, int size);
double avg_function(double *data, int inhi_size);
double inhi_avg_function(double *data, int inhi_size, int size);
double MEAN(double *data, int size);
double VAR(double *data, double mean);
double reduction(double *data1, double *data2, int size);
double bliss(double **data2, double *data1, int pro_size);
double spline(double *data1, double *avg, double *std, double *tm, int size, int interval);
double CV(double **data, double* final_avg, double* std, int size, int col);
double data_sort(double **data4, double *data1, double *data2, int size, int row_size);
double interCV(double *data1, double *avg, double *std, double *temp_std1, double *temp_std2, double *tm, int size, int interval);
