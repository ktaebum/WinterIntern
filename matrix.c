#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "matrix.h"

double **make_matrix(int row, int col)
{
	double **matrix = (double**)malloc(row * sizeof(double*));
	for(int i = 0; i < row; i++)
		matrix[i] = (double*)malloc(col * sizeof(double));
	return matrix;
}

double **multi(double **matrix1, int row1, int col1, double **matrix2, int row2, int col2)
{
    double **result = NULL;

    /* Error Handling; Cannot Define Matrix Multiplication */
    if(col1 != row2)
        return NULL;

    /* result matrix's size is (row1 by col2) */
    result = (double **) malloc(row1 * sizeof(double *));

    for(int i = 0; i < row1; i++)
    {
        result[i] = (double *) malloc(col2 * sizeof(double));

        for(int j = 0; j < col2; j++)
            result[i][j] = 0;
    }


    for(int k = 0; k < col1; k++)
    {
        for(int i = 0; i < row1; i++)
        {
            for(int j = 0; j < col2; j++)
                result[i][j] += matrix1[i][k] * matrix2[k][j];
        }
    }

    return result;
}


double **transp(double **matrix, int row, int col)
{
    double **result = (double **) malloc(col * sizeof(double *));

    for(int i = 0; i < col; i++)
    {
        result[i] = (double *) malloc(row * sizeof(double));

        for(int j = 0; j < row; j++)
            result[i][j] = matrix[j][i];
    }

    return result;
}

double determinant(double **matrix, int size)
{
    /* Get Determinant
     * Base Case 1: size == 2
     * det(A) = ad-bc
     * Base Case 2: size == 1
     * The Only element is deteminant
     * Using recursion */
    if(size == 2)
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    else if(size == 1)
        return matrix[0][0];
    else
    {
        double sum = 0;
        double **cofactor = NULL;

        for(int i = 0; i < size; i++)
        {
            cofactor = (double **) malloc((size - 1) * sizeof(double *));

            for(int j = 0; j < size - 1; j++)
            {
                cofactor[j] = (double *) malloc((size - 1) * sizeof(double));

                for(int k = 0; k < size; k++)
                {
                    if(k == i)
                        continue;
                    else if(k > i)
                        cofactor[j][k - 1] = matrix[j + 1][k];
                    else
                        cofactor[j][k] = matrix[j + 1][k];
                }
            }

            sum += pow(-1, i) * matrix[0][i] * determinant(cofactor, size - 1);

            for(int k = 0; k < size - 1; k++)
                free(cofactor[k]);

            free(cofactor);
        }

        return sum;
    }
}

double **inverse(double **matrix, int size)
{
    /* Get Inverse Matrix
     * Mathmatical Idea: for matrix A, A_inverse = (C^T)/det(A)
     * C is cofactor of A */
    double **result = (double **) malloc(size * sizeof(double *));

    for(int i = 0; i < size; i++)
        result[i] = (double *) malloc(size * sizeof(double));

    /* This is a real return value */
    double **inverse = NULL;
    double detA = determinant(matrix, size);

    /* If Determinant is 0, there is no inverse matrix */
    if(detA == 0)
        return NULL;

    double **cofactor = NULL;
    int row, col;

    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            row = 0;
            cofactor = (double **) malloc((size - 1) * sizeof(double *));

            for(int l = 0; l < size; l++)
            {
                if(l == i)
                    continue;

                col = 0;
                cofactor[row] = (double *) malloc((size - 1) * sizeof(double));

                for(int m = 0; m < size; m++)
                {
                    if(m == j)
                        continue;

                    cofactor[row][col] = matrix[l][m];
                    col++;
                }

                row++;
            }

            result[i][j] = pow(-1, i + j) * determinant(cofactor, size - 1) / detA;

            for(int l = 0; l < size - 1; l++)
                free(cofactor[l]);

            free(cofactor);

        }
    }

    inverse = transp(result, size, size);

    for(int i = 0; i < size; i++)
        free(result[i]);

    free(result);


    return inverse;
}


void print_matrix(double **matrix, int row, int col)
{
	for(int i = 0; i < row; i++)
	{
		printf("[");
		for(int j = 0; j < col; j++)
		{
			printf(" %lf ", matrix[i][j]);
		}
		printf("]\n");
	}
}


void free_matrix(double **matrix, int row)
{
	for(int i = 0; i < row; i++)
		free(matrix[i]);
	free(matrix);
}



