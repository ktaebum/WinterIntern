#ifndef _FIT_H
#define _FIT_H

#include <python2.7/Python.h>
#include <numpy/arrayobject.h>


double getdist(double x, double y, double z);
static PyObject *setmat(PyObject *self, PyObject *args);
static char setmat_docstring[] = 
	"Set Matrix for solve Eigenvalue Problem";
static PyObject *compRho(PyObject *self, PyObject *args);
static char compRho_docstring[] =
	"Chi-square function of ferrers bar";

static PyObject *bulgefit(PyObject *self, PyObject *args);
static char bulgefit_docstring[] = 
	"fitting bulge of the galaxy";

static PyObject *withbulge(PyObject *self, PyObject *args);
static char withbulge_docstring[] =
	"fitting summation of bulge and bar";

static struct PyMethodDef methods[] =
{
    {"setmat", setmat, METH_VARARGS, setmat_docstring},
    {"compbar", compRho, METH_VARARGS, compRho_docstring},
    {"compbulge", bulgefit, METH_VARARGS, bulgefit_docstring},
	{"fitbulge", withbulge, METH_VARARGS, withbulge_docstring},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC inittbfit(void)
{
    PyObject *m = Py_InitModule("tbfit", methods);

    if(m == NULL)
 	{   
        printf("Failed to Load Module\n");
        return;
    }

    import_array();

}

#endif
