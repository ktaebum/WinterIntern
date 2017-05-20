#ifndef _MATRIX_H
#define _MATRIX_H

double **make_matrix(int row, int col);
double **multi(double **matrix1, int row1, int col1, double **matrix2, int row2, int col2);
double **transp(double **matrix, int row, int col);
double determinant(double **matrix, int size);
double **inverse(double **matrix, int size);
void print_matrix(double **matrix, int row, int col);
void free_matrix(double **matrix, int row);





#endif
