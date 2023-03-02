//----------------------------------
// FUNCTION PROTPTYPES
//----------------------------------
void initialdata(double *protein, int protein_col, double *pre_weight, int pre_weight_num);
double** mallocDouble(int row, int col);
double* mallocSingle(int col);
void freeMallocDouble(double **data, int row);
void freeMallocSingle(double *data);
