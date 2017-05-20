#include "pyarray.h"

double **pymatrix_to_Cmatrixptr(PyArrayObject *array)
{
    /* Creates C array that point to the NUMPY data.
     * Memory allocate(must free) */
    double **matrix, *numpy;
    int row, col;

    row = array->dimensions[0]; 
    col = array->dimensions[1]; 
	matrix = (double **)malloc(row * sizeof(double *));
    numpy = (double *) array->data;

    for(int i = 0; i < row; i++)
        matrix[i] = numpy + i * col;

    return matrix;
}

void free_Cmatrixptr(double **v)
{
    free(v);
}
double ***pytensor_to_Ctensorptr(PyArrayObject *array)
{
    double ***tensor, *numpy;
    int row, col1, col2;
    row = array->dimensions[0];
    col1 = array->dimensions[1];
    col2 = array->dimensions[2];
	tensor = (double ***)malloc(row * sizeof(double **));
	for(int i = 0; i < row; i++)
		tensor[i] = (double **)malloc(col1 * sizeof(double *));
    numpy = (double *)array->data;


    for(int i = 0; i < row; i++)
    {
        for(int j = 0; j < col1; j++)
            tensor[i][j] = numpy + i * col1 * col2  + j * col2;
    }

    return tensor;
}


void free_Ctensorptr(double ***t, int row)
{
	for(int i = 0; i < row; i++)
	{
		free(t[i]);
	}
	free(t);

}

double *pyvector_to_Cvectorptr(PyArrayObject *array)
{
    return (double *) array->data;  /* pointer to arrayin data as double */
}

long **pyintmatrix_to_Clongmatrixptr(PyArrayObject *array)
{
    long **matrix, *numpy;
    int row, col;

    row = array->dimensions[0];
    col = array->dimensions[1];
	matrix = (long **)malloc(row * sizeof(long *));
    numpy = (long *) array->data;

    for (int i = 0; i < row; i++)
        matrix[i] = numpy + i * col;
    return matrix;
}
void free_Clongmatrixptr(long **v)
{
    free(v);
}
