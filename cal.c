#include "cal.h"
#include "pyarray.h"
double getdist(double x, double y, double z)
{
    return pow(x, 2) + pow(y, 2) + pow(z, 2);
}
int findWhere(double *grid, int n, double pos, double ds)
{
    int index;
    int first = 0, last = n - 1, mid;
    double grid_val;

    while(1)
    {
        mid = (first + last) / 2;
		grid_val = grid[mid];
        if(pos < grid_val + ds / 2. && pos > grid_val - ds / 2.)
        {
            index = mid;
            break;
        }
        else
        {
            if(pos >= grid_val + ds / 2.)
            {
                if(mid == (n - 2))
                {
                    index = n - 1;
                    break;
                }
                else
                    first = mid;
            }
            else if(pos <= grid_val - ds / 2.)
            {
                if(mid == 0)
                {
                    index = 0;
                    break;
                }
                else
                    last = mid;
            }
        }
    }

    return index;

}

static PyObject *calCOM(PyObject *self, PyObject *args)
{
    /* Get Center of Mass of System
     * _p is position vector
     * _m is mass array
     * tot_mass is sum of all mass*/
    PyArrayObject *position_in;
    PyArrayObject *mass_in;
    double tot_mass;

    if(!PyArg_ParseTuple(args, "O!O!d", &PyArray_Type, &position_in, &PyArray_Type, &mass_in, &tot_mass))
        return NULL;
	printf("\n\n========= Using C-API =========\n");
	printf("Calculating Center of Mass of System...\n");
    double x_c = 0, y_c = 0, z_c = 0;
	double **position = pymatrix_to_Cmatrixptr(position_in);
	double *mass = pyvector_to_Cvectorptr(mass_in);
    for(int i = 0; i < mass_in->dimensions[0]; i++)
    {
		x_c += position[i][0] * mass[i];
		y_c += position[i][1] * mass[i];
		z_c += position[i][2] * mass[i];
    }
	printf("x_center = %.10lf\n", x_c / tot_mass);
	printf("y_center = %.10lf\n", y_c / tot_mass);
	printf("z_center = %.10lf\n", z_c / tot_mass);
	printf("Calculating Finished\n");
	printf("========= Exit C-API ==========\n\n\n");

	free_Cmatrixptr(position);
    return Py_BuildValue("ddd", x_c / tot_mass, y_c / tot_mass, z_c / tot_mass);
}

static PyObject *calRho1d(PyObject *self, PyObject *args)
{
    PyArrayObject *pos_in;
    PyArrayObject *mass_in;
    PyArrayObject *grid_in;
    double ds;

    if(!PyArg_ParseTuple(args, "O!O!O!d", &PyArray_Type, &pos_in, &PyArray_Type, &mass_in, &PyArray_Type, &grid_in, &ds))
        return NULL;

	double *pos = pyvector_to_Cvectorptr(pos_in);
	double *mass = pyvector_to_Cvectorptr(mass_in);
	double *grid = pyvector_to_Cvectorptr(grid_in);
	printf("\n\n========= Using C-API =========\n");
	printf("Calculating 1-Dimensional Density\n");
	printf("Smoothing: Divide by 3\n");
    int n = grid_in->dimensions[0];
    npy_intp dim[1] = {n};
    PyArrayObject *rho = (PyArrayObject *)PyArray_ZEROS(1, dim, NPY_DOUBLE, 0);
	double *c_rho = pyvector_to_Cvectorptr(rho);
    int index = 0;

    for(int i = 0; i < mass_in->dimensions[0]; i++)
    {
        index = findWhere(grid, n, pos[i], ds);
    	for(int j = index - 1; j <= index + 1; j++)
		{
			if(j == -1 || j == n)
				continue;
			c_rho[j] += mass[i] / 3.;
		}
	}
	printf("Calculating Finished\n");
	printf("========= Exit C-API ==========\n\n\n");
    return Py_BuildValue("N", rho);
}

static PyObject *calRho2d(PyObject *self, PyObject *args)
{
    PyArrayObject *x_pos_in;
    PyArrayObject *y_pos_in;
    PyArrayObject *mass_in;
    PyArrayObject *grid_in;
	double *x_pos, *y_pos, *mass, *grid;
    double ds;

    if(!PyArg_ParseTuple(args, "O!O!O!O!d", &PyArray_Type, &x_pos_in, &PyArray_Type, &y_pos_in, 
				&PyArray_Type, &mass_in, &PyArray_Type, &grid_in, &ds))
        return NULL;

	double **rho_out;
	int n;
	printf("\n\n========= Using C-API =========\n");
	printf("Calculating 2-Dimensional Density\n");
	printf("Smoothing: Divide by 9\n");
    npy_intp dim[2];
	n = dim[0] = dim[1] = grid_in->dimensions[0];
    PyArrayObject *xy_rho = (PyArrayObject *)PyArray_ZEROS(2, dim, NPY_DOUBLE, 0);
	x_pos = pyvector_to_Cvectorptr(x_pos_in);
	y_pos = pyvector_to_Cvectorptr(y_pos_in);
	mass = pyvector_to_Cvectorptr(mass_in);
	grid = pyvector_to_Cvectorptr(grid_in);
	rho_out = pymatrix_to_Cmatrixptr(xy_rho);
	

    int x = 0, y = 0;

    for(int i = 0; i < mass_in->dimensions[0]; i++)
    {
        x = findWhere(grid, n, x_pos[i], ds);
        y = findWhere(grid, n, y_pos[i], ds);
      
		for(int j = y - 1; j <= y + 1; j++)
		{
			if(j == - 1 || j == n)
				continue;
			for(int k = x - 1; k <= x + 1; k++)
			{
				if(k == -1 || k == n)
					continue;
				rho_out[j][k] += mass[i] / 9.;
			}
		}
	}
	printf("Calculating Finished\n");
	printf("========= Exit C-API ==========\n\n\n");
	free_Cmatrixptr(rho_out);
    return Py_BuildValue("N", xy_rho);
}

static PyObject *calRho3d(PyObject *self, PyObject *args)
{
    PyArrayObject *pos_in;
    PyArrayObject *mass_in;
    PyArrayObject *grid_in;
    double ds;
	int type = 0;
    if(!PyArg_ParseTuple(args, "O!O!O!d|i", &PyArray_Type, &pos_in, &PyArray_Type, 
				&mass_in, &PyArray_Type, &grid_in, &ds, &type))
        return NULL;
	
	printf("\n\n========= Using C-API =========\n");
	printf("Calculating 3-Dimensional Density\n");
	printf("Smoothing: Divide by 27\n");
	if(type != 0)
		printf("Return Index Matrix: True\n");
	else
		printf("Return Index Matrix: False\n");
    int n = grid_in->dimensions[0];
	int npart = mass_in->dimensions[0];
    npy_intp dim[3] = {n, n, n};
    PyArrayObject *rho = (PyArrayObject *)PyArray_ZEROS(3, dim, NPY_DOUBLE, 0);
	PyArrayObject *id = NULL;
	double ***rho_out = pytensor_to_Ctensorptr(rho);
	long **id_out = NULL;
	double **pos = pymatrix_to_Cmatrixptr(pos_in);
	double *mass = pyvector_to_Cvectorptr(mass_in);
	double *grid = pyvector_to_Cvectorptr(grid_in);
	if(type != 0)
	{
		int dim2[2] = {npart, 3};
		id = (PyArrayObject *)PyArray_FromDims(2, dim2, NPY_LONG);
		id_out = pyintmatrix_to_Clongmatrixptr(id);
	}
    int x = 0, y = 0, z = 0;

    for(int i = 0; i < npart; i++)
    {
        x = findWhere(grid, n, pos[i][0], ds);
        y = findWhere(grid, n, pos[i][1], ds);
        z = findWhere(grid, n, pos[i][2], ds);
		for(int j = z - 1; j <= z + 1; j++)
		{
			if(j == -1 || j == n)
				continue;
			for(int k = y - 1; k <= y + 1; k++)
			{
				if(k == -1 || k == n)
					continue;
				for(int l = x - 1; l <= x + 1; l++)
				{
					if(l == -1 || l == n)
						continue;
					rho_out[j][k][l] += mass[i] / 27.;
				}
			}
		}

		if(type != 0)
		{
			int copy = 0;
			for(int j = 0; j < 3; j++)
			{
				switch(j){
					case 0:
						copy = z;
						break;
					case 1:
						copy = y;
						break;
					case 2:
						copy = x;
						break;
					default:
						break;
				}
				id_out[i][j] = copy;
			}
		}

    }
	printf("Calculating Finished\n");
	printf("========= Exit C-API ==========\n\n\n");
	free_Cmatrixptr(pos);
	free_Ctensorptr(rho_out, n);
	if(type == 0)
    	return Py_BuildValue("N", rho);
	else
	{
		free_Clongmatrixptr(id_out);
		return Py_BuildValue("NN", rho, id);
	}
}


static PyObject *calRho3d_removebulge(PyObject *self, PyObject *args)
{
	PyArrayObject *grid_in;
	PyArrayObject *pos_in;
	PyArrayObject *mass_in;
	double bulge;
	double percent;
	if(!PyArg_ParseTuple(args, "O!O!O!dd", &PyArray_Type, &grid_in, &PyArray_Type, &pos_in, 
				&PyArray_Type, &mass_in, &bulge, &percent))
		return NULL;

	double *grid = pyvector_to_Cvectorptr(grid_in);
	double **pos = pymatrix_to_Cmatrixptr(pos_in);
	double *mass = pyvector_to_Cvectorptr(mass_in);

	double ds = grid[1] - grid[0];
	int n = grid_in->dimensions[0];
	int npart = mass_in->dimensions[0];
    npy_intp dim[3] = {n, n, n};
    PyArrayObject *rho = (PyArrayObject *)PyArray_ZEROS(3, dim, NPY_DOUBLE, 0);
	double ***rho_out = pytensor_to_Ctensorptr(rho);
	int x = 0, y = 0, z = 0;

    for(int i = 0; i < npart; i++)
    {
        x = findWhere(grid, n, pos[i][0], ds);
        y = findWhere(grid, n, pos[i][1], ds);
        z = findWhere(grid, n, pos[i][2], ds);
		for(int j = z - 1; j <= z + 1; j++)
		{
			if(j == -1 || j == n)
				continue;
			for(int k = y - 1; k <= y + 1; k++)
			{
				if(k == -1 || k == n)
					continue;
				for(int l = x - 1; l <= x + 1; l++)
				{
					if(l == -1 || l == n)
						continue;
					if(sqrt(getdist(grid[x], grid[y], grid[z])) < bulge)
					{
						/* If Inside The Bulge */
						rho_out[j][k][l] += (mass[i] * percent) / 27.;
					}
					else
						rho_out[j][k][l] += mass[i] / 27.;
				}
			}
		}
	}
	free_Cmatrixptr(pos);
	free_Ctensorptr(rho_out, n);
	return Py_BuildValue("N", rho);
}


static PyObject *calInertia(PyObject *self, PyObject *args)
{
    PyArrayObject *disk_in;
    PyArrayObject *mass_in;

    if(!PyArg_ParseTuple(args, "OO", &disk_in, &mass_in))
        return NULL;
	

    int dim[2] = {3, 3};
    PyArrayObject *I = (PyArrayObject *)PyArray_FromDims(2, dim, NPY_DOUBLE);
	double **disk = pymatrix_to_Cmatrixptr(disk_in);
	double *mass = pyvector_to_Cvectorptr(mass_in);
	double **Iptr = pymatrix_to_Cmatrixptr(I);
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
			Iptr[i][j] = 0;
            for(int k = 0; k < disk_in->dimensions[1]; k++)
            {
				Iptr[i][j] += mass[i] * (getdist(disk[0][k], disk[1][k], disk[2][k]) * (i == j) 
						- disk[j][k] * disk[i][k]);
            }
        }
    }
	free_Cmatrixptr(disk);
	free_Cmatrixptr(Iptr);
    return Py_BuildValue("N", I);
}

static PyObject *findMax(PyObject *self, PyObject *args)
{
	PyArrayObject *rho_in;
	
	if(!PyArg_ParseTuple(args, "O!", &PyArray_Type, &rho_in))
		return NULL;

	double **rho = pymatrix_to_Cmatrixptr(rho_in);
	double max = 0;
	int y = 0, x = 0;

	for(int i = 0; i < rho_in->dimensions[0]; i++)
	{
		for(int j = 0; j < rho_in->dimensions[1]; j++)
		{
			if(rho[i][j] > max)
			{
				max = rho[i][j];
				y = i;
				x = j;
			}
		}
	}

	free_Cmatrixptr(rho);
	return Py_BuildValue("ii", y, x);
}

double vector_mag(double x, double y)
{
	return sqrt(pow(x, 2) + pow(y, 2));
}

static PyObject *caltheta(PyObject *self, PyObject *args)
{
	PyArrayObject *pos_in;
	PyArrayObject *vel_in;
	
	if(!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &pos_in, &PyArray_Type, &vel_in))
		return NULL;

	double **pos = pymatrix_to_Cmatrixptr(pos_in);
	double **vel = pymatrix_to_Cmatrixptr(vel_in);
	
	int npart[] = {pos_in->dimensions[0]};
	PyArrayObject *theta_out = (PyArrayObject*)PyArray_FromDims(1, npart, NPY_DOUBLE);
	PyArrayObject *vel_out = (PyArrayObject*)PyArray_FromDims(1, npart, NPY_DOUBLE);
	double *theta = pyvector_to_Cvectorptr(theta_out);
	double *scalar_vel = pyvector_to_Cvectorptr(vel_out);
	for(int i = 0; i < npart[0]; i++)
	{
		scalar_vel[i] = vector_mag(vel[i][0], vel[i][1]);
		theta[i] = acos((pos[i][0] * vel[i][0] + pos[i][1] * vel[i][1])
			/ (vector_mag(pos[i][0], pos[i][1]) * scalar_vel[i]));
		
	}

	free_Cmatrixptr(vel);
	free_Cmatrixptr(pos);
	return Py_BuildValue("OO", vel_out, theta_out);
}

