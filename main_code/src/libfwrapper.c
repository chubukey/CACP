/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* libwrapper = fortran wrappers for LAMMPS library functions.
   See README for compilation instructions */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "stdint.h"
#include "library.h"        /* this is a LAMMPS include file */

/* wrapper for creating a lammps instance from fortran.
   since fortran has no simple way to emit a C-compatible
   argument array, we don't support it. for simplicity,
   the address of the pointer to the lammps object is
   stored in a 64-bit integer on all platforms. */

void lammps_open_(MPI_Fint *comm, int64_t *ptr)
{
    void *obj;
    MPI_Comm ccomm;

    /* convert MPI communicator from fortran to c */
    ccomm = MPI_Comm_f2c(*comm);

    lammps_open(0,NULL,ccomm,&obj);
    *ptr = (int64_t) obj;
}

/* no-MPI version of the wrapper from above. */

void lammps_open_no_mpi_(int64_t *ptr)
{
    void *obj;

    lammps_open_no_mpi(0,NULL,&obj);
    *ptr = (int64_t) obj;
}

/* wrapper for shutting down a lammps instance from fortran. */

void lammps_close_(int64_t *ptr)
{
    void *obj;
    obj = (void *) *ptr;

    lammps_close(obj);
}

/* wrapper for passing an input file to lammps from fortran.
   since fortran strings are not zero terminated, we have
   to pass the length explicitly and make a copy that is. */

void lammps_file_(int64_t *ptr, char *fname, MPI_Fint *len)
{
    void *obj;
    char *cpy;

    obj = (void *) *ptr;

    cpy = (char *)calloc(*len + 1,sizeof(char));
    memcpy(cpy,fname,*len);

    lammps_file(obj,cpy);
    free(cpy);
}

/* wrapper for passing a line input to lammps from fortran.
   since fortran strings are not zero terminated, we have
   to pass the length explicitly and make a copy that is. */

void lammps_command_(int64_t *ptr, char *line, MPI_Fint *len)
{
    void *obj;
    char *cpy;

    obj = (void *) *ptr;
    cpy = (char *)calloc(*len + 1,sizeof(char));
    memcpy(cpy,line,*len);
//    printf("%s\n",cpy);
    lammps_command(obj,cpy);
    free(cpy);
}

/* fortran wrapper to get the number of atoms from lammps.
   return values require an interface in fortran, so we 
   make the wrapper into a procedure. */

void lammps_get_natoms_(int64_t *ptr, MPI_Fint *natoms)
{
    void *obj;
    obj = (void *) *ptr;
    
    *natoms = lammps_get_natoms(obj);
}

/* wrapper to copy coordinates from lammps to fortran */

void lammps_extract_global_(int64_t *ptr, char *name, double *data, MPI_Fint *len, MPI_Fint *size_count)
{
   void *obj;
   obj = (void *) *ptr;
   char *cpy;
   double *vector;
   cpy = (char *)calloc(*len+1,sizeof(char));
   memcpy(cpy,name,*len);
   vector = (double*)lammps_extract_global(obj,cpy);
   int i;
   int count = *size_count;
   for(i=0;i<count;i++)   data[i]  = vector[i];
   free(cpy);
}
//     call lammps_extract_global(ptr_lammps,'boxxlo',MD_boxdim(1),6)


// NOTE: this is now out-of-date, needs to be updated to lammps_gather_atoms()

//void lammps_gather_atoms_(int64_t *ptr, MPI_Fint *natoms, MPI_Fint *len,  char *name, double *data)

// CALL lammps_gather_atoms(ptr,'fext',1,3,fext,5);

void lammps_gather_atoms_(int64_t *ptr, char *name,MPI_Fint *type, MPI_Fint *count, void *data,MPI_Fint *len)
{
    void *obj;
    obj = (void *) *ptr;
    char *cpy;
    cpy = (char *)calloc(*len+1,sizeof(char));
    memcpy(cpy,name,*len);
    lammps_gather_atoms(obj,cpy,*type,*count,data);
    free(cpy);
}

//void lammps_scatter_atoms(void *,char *, int, int, void *);

//void lammps_gather_atoms_(int64_t *ptr, MPI_Fint *natoms, MPI_Fint *len,  char *name, double *data)
//{
//    void *obj;
//    obj = (void *) *ptr;
//    char *cpy;
//    cpy = (char *)calloc(*len+1,sizeof(char));
//    memcpy(cpy,name,*len);
//    *natoms = lammps_get_natoms(obj);
//    lammps_gather_atoms(obj,cpy,1,3,data);
//    free(cpy);
//}

void lammps_scatter_atoms_(int64_t *ptr,char *name, MPI_Fint *type, MPI_Fint *count,void *data, MPI_Fint *len)
{
    void *obj;
    obj = (void *) *ptr;
    char *cpy;
    cpy = (char *)calloc(*len+1,sizeof(char));
    memcpy(cpy,name,*len);
    
    lammps_scatter_atoms(obj,cpy,*type,*count,data);
    free(cpy);
}

void lammps_scatter_global_(int64_t *ptr,char *name, MPI_Fint *type, MPI_Fint *count,void *data, MPI_Fint *len)
{
    // because in marcc 64 bit is used, and here assuming all integer used 32 bit is enough, thus for 64 bit data, first 32 bit = original value, next 32 bit is always zero
   
    void *obj;
    obj = (void *) *ptr;
    char *cpy; 
    cpy = (char *)calloc(*len+1,sizeof(char));
    memcpy(cpy,name,*len);
    lammps_scatter_global(obj,cpy,*type,*count,data);
    free(cpy);
}

// only for fortran without C_binding
// default data structure is double
void lammps_extract_variable_(int64_t *ptr, char *name, char *group, MPI_Fint *len_name, MPI_Fint *len_group,MPI_Fint *dimension, double *data)
{
    void *obj;
    obj = (void *) *ptr;
    char *cpy_name;
    char *cpy_group;
    double *vector;
    cpy_name = (char *)calloc(*len_name+1,sizeof(char));
    cpy_group = (char *)calloc(*len_group+1,sizeof(char));
    memcpy(cpy_name,name,*len_name);
    memcpy(cpy_group,group,*len_group);
    vector = (double *) lammps_extract_variable(obj,cpy_name,cpy_group);
//    printf("vector value is%f %f\n",vector[0],vector[1]);
    int i;
    for(i=0;i<*dimension;i++)
    {
       data[i] = vector[i];
    }
//    printf("data value is%f %f\n",data[0],data[1]);
    free(vector);
    free(cpy_name);
    free(cpy_group);
}

void lammps_gather_variable_(int64_t *ptr, char *name, char *group, MPI_Fint *len_name, MPI_Fint *len_group, MPI_Fint *type, MPI_Fint *count, double *data)
{
    void *obj;
    obj = (void *) *ptr;
    char *cpy_name;
    char *cpy_group;
    cpy_name = (char *)calloc(*len_name+1,sizeof(char));
    cpy_group = (char *)calloc(*len_group+1,sizeof(char));
    memcpy(cpy_name,name,*len_name);
    memcpy(cpy_group,group,*len_group);
    lammps_gather_variable(obj,cpy_name,cpy_group,*type,*count,data);
 //   printf("data value is%f %f\n",data[0],data[1]);
    free(cpy_name);
    free(cpy_group);
}

void lammps_extract_compute_(int64_t *ptr, char *name, MPI_Fint *len_name, MPI_Fint *type, MPI_Fint *count, double *data)
{
    void *obj;
    obj = (void *) *ptr;
    char *cpy_name;
    int style = 0;         // currently this only works for global type of compute, i.e. not per-atom or local
    cpy_name = (char *)calloc(*len_name+1,sizeof(char));
    memcpy(cpy_name,name,*len_name);
//    lammps_extract_compute(obj,cpy_name,cpy_group,*type);
//    printf("before use lammps_extract, %d %d\n",*style,*type);
//    void *happy = lammps_extract_compute(obj,cpy_name,*style,*type);
//    double *vector = (double *)lammps_extract_compute(obj,cpy_name,*style,*type);
//    void *vector = (double *)lammps_extract_compute(obj,cpy_name,0,1);
    int i;
    void *happy = lammps_extract_compute(obj,cpy_name,style,*type);
    double *vector = (double *)happy;
    int n = (int) *count;
    for(i=0;i<n;i++)
    {
       data[i] = vector[i];
    }
//   printf("n is %d\n",n);
//    printf("data value is:%f %f %f %f\n",data[0],data[1],data[n-2],data[n-1]);
    free(cpy_name);
}



/* wrapper to copy coordinates from fortran to lammps */

/* NOTE: this is now out-of-date, needs to be updated to lammps_scatter_atoms()

void lammps_put_coords_(int64_t *ptr, double *coords)
{
    void *obj;
    obj = (void *) *ptr;

    lammps_put_coords(obj,coords);
}

*/
