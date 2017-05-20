#ifndef _PYARRAY_H
#define _PYARRAY_H

#include <python2.7/Python.h>
#include <numpy/arrayobject.h>

double **pymatrix_to_Cmatrixptr(PyArrayObject *array);
void free_Cmatrixptr(double **v);

double ***pytensor_to_Ctensorptr(PyArrayObject *array);
void free_Ctensorptr(double ***t, int n);

double *pyvector_to_Cvectorptr(PyArrayObject *array);

long **pyintmatrix_to_Clongmatrixptr(PyArrayObject *array);
void free_Clongmatrixptr(long **v);


#endif
