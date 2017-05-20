#ifndef _CAL_H
#define _CAL_H

#include <python2.7/Python.h>
#include <numpy/arrayobject.h>

static PyObject *calCOM(PyObject *self, PyObject *args);
static char calCOM_docstring[] =
	"Calculate Center of Mass in the System.";
static PyObject* calRho1d(PyObject* self, PyObject* args);
static char calRho1d_docstring[] =
    "Caculate Sliced Line Density of Galaxy.";
static PyObject* calRho2d(PyObject* self, PyObject* args);
static char calRho2d_docstring[] =
    "Caculate Sliced Surface Density of Galaxy.";
static PyObject* calRho3d(PyObject* self, PyObject* args);
static char calRho3d_docstring[] =
    "Caculate Total Density of Galaxy.";

static PyObject *calRho3d_removebulge(PyObject *self, PyObject *args);
static char calRho3d_removebulge_docstring[] =
	"Calculate Total Density smoothing bulge region.";
static PyObject *calInertia(PyObject *self, PyObject *args);
static char calInertia_docstring[] =
	"Get Inertia Tensor.";
static PyObject *findMax(PyObject *self, PyObject *args);
static char findMax_docstring[] =
	"find Maximum index of matrix.";
static PyObject *caltheta(PyObject *self, PyObject *args);
static char caltheta_docstring[] =
	"find theta between vel and pos to get radial velocity";
double getdist(double x, double y, double z);
int findWhere(double *grid, int n, double pos, double ds);


static struct PyMethodDef methods[] =
{
    {"calrho1d", calRho1d, METH_VARARGS, calRho1d_docstring},
    {"calrho2d", calRho2d, METH_VARARGS, calRho2d_docstring},
    {"calrho3d", calRho3d, METH_VARARGS, calRho3d_docstring},
    {"calrho3d_nobulge", calRho3d_removebulge, METH_VARARGS, calRho3d_removebulge_docstring},
    {"calcom", calCOM, METH_VARARGS, calCOM_docstring},
    {"calinertia", calInertia, METH_VARARGS, calInertia_docstring},
    {"findmax", findMax, METH_VARARGS, findMax_docstring},
    {"caltheta", caltheta, METH_VARARGS, caltheta_docstring},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC inittbcal(void)
{
    PyObject* m = Py_InitModule("tbcal", methods);
    if (m == NULL)
    {
        printf("Failed To Load Module\n");
        return;
    }
	import_array();
}
#endif
