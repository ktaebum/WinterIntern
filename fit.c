#include "fit.h"
#include "matrix.h"
#include "pyarray.h"


double getdist(double x, double y, double z)
{
	return pow(x, 2) + pow(y, 2) + pow(z, 2);
}

static PyObject *setmat(PyObject *self, PyObject *args)
{
	printf("\n\n========= Using C-API =========\n");
	printf("Set Matrix Equation: E * lambda = I * lambda\n");
	PyArrayObject *candidate_in;
	int len;
    if(!PyArg_ParseTuple(args, "O!", &PyArray_Type, &candidate_in))
        return NULL;

	len = candidate_in->dimensions[1];
	double **candidate = pymatrix_to_Cmatrixptr(candidate_in);
	double **D = make_matrix(len, 4);
	for(int i = 0; i < len; i++)
	{
		D[i][0] = pow(candidate[0][i], 2);
		D[i][1] = candidate[0][i] * candidate[1][i];
		D[i][2] = pow(candidate[1][i], 2);
		D[i][3] = 1;
	}
	double **D_T = transp(D, len, 4);
	double **S = multi(D_T, 4, len, D, len, 4);
	
	double **S12 = make_matrix(3, 1);
	for(int i = 0; i < 3; i++)
		S12[i][0] = S[i][3];
	
	double **S21 = make_matrix(1, 3);
	memcpy(S21[0], S[3], 3 * sizeof(double));
	
	double S22;
	S22 = 1./S[3][3];
	
	double **C_inv = make_matrix(4, 4);
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			if(i + j == 2)
			{
				if(i == 2 || j == 2)
					C_inv[i][j] = 0.5;
				else
					C_inv[i][j] = -1;
			}
			else
				C_inv[i][j] = 0;
		}
	}

	double **S_temp;
	double **E_copy;

	S_temp = multi(S12, 3, 1, S21, 1, 3);
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			S[i][j] -= S_temp[i][j] * S22;
		}
	}
	E_copy = multi(C_inv, 4, 4, S, 4, 4);
	printf("Finished Get \'E\'\n");
	printf("Now Save in to PyArray and Return\n");
	int dim2[2] = {3, 3};
	int dim1[1] = {3};
	PyArrayObject *E = (PyArrayObject*) PyArray_FromDims(2, dim2, NPY_DOUBLE);
	double **E_C = pymatrix_to_Cmatrixptr(E);
	PyArrayObject *s21 = (PyArrayObject*) PyArray_FromDims(1, dim1, NPY_DOUBLE);
	double *s21_C = pyvector_to_Cvectorptr(s21);
	for(int i = 0; i < 3; i++)
	{
		s21_C[i] = S21[0][i];
		for(int j = 0; j < 3; j++)
			E_C[i][j] = E_copy[i][j];
	}

	free_matrix(S12, 3);
	free_matrix(S21, 1);
	free_matrix(S_temp, 3);
	free_matrix(D_T, 4);
	free_matrix(D, len);
	free_matrix(S, 4);
	free_matrix(E_copy, 4);
	free_matrix(C_inv, 4);
	free_Cmatrixptr(candidate);
	free_Cmatrixptr(E_C);
	printf("========= Exit C-API ==========\n\n\n");
	return Py_BuildValue("NNd", E, s21, S22);
}
double gsquare(double x, double y, double z, double major, double minor)
{
    return (pow(x / major, 2) + (pow(y , 2) + pow(z, 2)) / pow(minor, 2));
}
static PyObject *compRho(PyObject *self, PyObject *args)
{
    PyArrayObject *rho_in;
    PyArrayObject *pos_in;
	PyArrayObject *index_in;
    PyArrayObject *vector_in;
    double rho_bar;

    if(!PyArg_ParseTuple(args, "O!O!O!O!d", &PyArray_Type, &vector_in, &PyArray_Type, &pos_in, 
				&PyArray_Type, &index_in, &PyArray_Type, &rho_in, &rho_bar))
        return NULL;
	double ***rho = pytensor_to_Ctensorptr(rho_in);
	double **pos = pymatrix_to_Cmatrixptr(pos_in);
	long **index = pyintmatrix_to_Clongmatrixptr(index_in);
	double *vector = pyvector_to_Cvectorptr(vector_in);
	int npart = pos_in->dimensions[0];
    double least = 0;
	double expect;
	double real;
	double gsq;
	int count = 0;
    for(int i = 0; i < npart; i++)
    {
		gsq = gsquare(pos[i][0], pos[i][1], pos[i][2], vector[1], vector[2]);
		if(gsq < 1)
        {
			count++;
			real = rho[index[i][0]][index[i][1]][index[i][2]];
			expect = rho_bar * pow(1 - gsq, vector[0]);
			least += pow(real - expect, 2);
        }
    }
	free_Cmatrixptr(pos);
	free_Clongmatrixptr(index);
	free_Ctensorptr(rho, rho_in->dimensions[0]);
    return Py_BuildValue("d", least / (double) count);
}

static PyObject *bulgefit(PyObject *self, PyObject *args)
{
	PyArrayObject *vector_in;
	PyArrayObject *pos_in;
	PyArrayObject *index_in;
	PyArrayObject *rho_in;
	double rho_bul;

	if(!PyArg_ParseTuple(args, "O!O!O!O!d", &PyArray_Type, &vector_in, &PyArray_Type, 
				&pos_in, &PyArray_Type, &index_in, &PyArray_Type, &rho_in, &rho_bul))
		return NULL;

	double **pos = pymatrix_to_Cmatrixptr(pos_in);
	double *vector = pyvector_to_Cvectorptr(vector_in);
	double ***rho = pytensor_to_Ctensorptr(rho_in);
	long **index = pyintmatrix_to_Clongmatrixptr(index_in);
	
	int count = 0;
	double rhocomp = 0.;
	double dist;
	for(int i = 0; i < pos_in->dimensions[0]; i++)
	{
		dist = getdist(pos[i][0], pos[i][1], pos[i][2]);
		/* If inside the bulge */
		if(dist < vector[0])
		{
			count++;
			rhocomp += pow(rho[index[i][0]][index[i][1]][index[i][2]] - 
					rho_bul * pow(1 + dist / vector[0] / vector[0], -1.5), 2);
		}
	}
	free_Cmatrixptr(pos);
	free_Clongmatrixptr(index);
	free_Ctensorptr(rho, rho_in->dimensions[0]);
	return Py_BuildValue("d", rhocomp / (double)count);
}

static PyObject *withbulge(PyObject *self, PyObject *args)
{

	PyArrayObject *vector_in;
	PyArrayObject *pos_in;
	PyArrayObject *index_in;
	PyArrayObject *rho_in;
	double rho_c;

	if(!PyArg_ParseTuple(args, "O!O!O!O!d", &PyArray_Type, &vector_in, &PyArray_Type, &pos_in,
				&PyArray_Type, &index_in, &PyArray_Type, &rho_in, &rho_c))
		return NULL;
	
	double ***rho = pytensor_to_Ctensorptr(rho_in);
	long **index = pyintmatrix_to_Clongmatrixptr(index_in);
	double **pos = pymatrix_to_Cmatrixptr(pos_in);
	/* vector[0] = rho_bulge
	 * rho_bar = rho_c - vector[0] 
	 * vector[1] = radius_bulge
	 * vector[2] = n
	 * vector[3] = major
	 * vector[4] = minor */
	double *vector = pyvector_to_Cvectorptr(vector_in);
	double dist, gsq;
	double rhocomp;
	double rho_bar = rho_c - vector[0];
	double rho_chi = 0;
	int count = 0;
	for(int i = 0; i < pos_in->dimensions[0]; i++)
	{
		rhocomp = 0;
		dist = getdist(pos[i][0], pos[i][1], pos[i][2]);
		gsq = gsquare(pos[i][0], pos[i][1], pos[i][2], vector[3], vector[4]);

		if(gsq < 1)
		{
			count++;
			rhocomp += rho_bar * pow(1 - gsq, vector[2]);
		}
		rhocomp += vector[0] * pow(1 + dist / vector[1] / vector[1], -1.5);

		rho_chi += pow(rho[index[i][0]][index[i][1]][index[i][2]] - rhocomp, 2);

	}

	free_Cmatrixptr(pos);
	free_Clongmatrixptr(index);
	free_Ctensorptr(rho, rho_in->dimensions[0]);
	return Py_BuildValue("d", rho_chi / (double)pos_in->dimensions[0]);
}



