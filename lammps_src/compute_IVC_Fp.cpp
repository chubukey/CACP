/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
     K-space terms added by Stan Moore (BYU)
------------------------------------------------------------------------- */
#include "mpi.h"
#include "string.h"
#include "compute_IVC_Fp.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "math.h"
#include "comm.h"
#include "math_const.h"
#include "stdio.h"
#include "domain.h"
#include "modify.h"
#include "math.h"
#include "memory.h"
#include "force.h"
#include "stdlib.h"

using namespace LAMMPS_NS;
using namespace MathConst;
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PERATOM 8
#define INVOKED_LOCAL 16

#define SMALL 0.00001

/* ---------------------------------------------------------------------- */


ComputeIVCFp::ComputeIVCFp(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{

// e.g. compute 1 lower group/group upper pair yes
  scalar_flag = vector_flag = 1;
  extscalar = 1;
  extvector = 1;

  if (narg < 5){
     if(me==0) printf("narg is less than 4\n");
     error->all(FLERR,"Illegal compute IVC_Fp command");
  }

// e.g. compute 1 lower group/group upper pair yes
  nnode_inter = force->inumeric(FLERR,arg[3]);
  compute_id = arg[4];  
  int nplane = 4;
  int nFp = 9;

  size_vector = nnode_inter*nplane*nFp;     // need to modify later

//  char* compute_id ="FeFp" ;  
  int icompute = modify->find_compute(compute_id);
  if(icompute==-1){
     error->all(FLERR,"cannot find FeFp compute");
  }


  me = comm->me;
//  printf("compute start icompute %d at proc %d\n",icompute,me);


  vector = new double[size_vector];
//  printf("vector allocated with size %d at proc %d\n",size_vector,me);
}

/* ---------------------------------------------------------------------- */

ComputeIVCFp::~ComputeIVCFp()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

double *ComputeIVCFp::compute_one()
{
  int i;
  me = comm->me;
  
  compute_id = "Ft";
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *IVC_id = atom->IVC_id;
  onecount = 0.0;
  allcount = 0.0;
//  printf("in proc %d, nlocal: %d; IVC_id: %d %d %d %d %d\n",me,nlocal,IVC_id[0],IVC_id[1],IVC_id[2],IVC_id[3],IVC_id[4]);

  nvalues = size_vector;
  double *onevec = new double[size_vector];
  for(i=0;i<nvalues;i++){
     onevec[i] = 0.0;
  }
  int icompute = modify->find_compute(compute_id);
  if(icompute==-1){
     error->all(FLERR,"cannot find FeFp compute");
  }
  int length = modify->compute[icompute]->size_vector;
  Compute *computeFeFp = modify->compute[icompute];
  
  double *Fp_local;
  double **Fp_local_plane;
  memory->create(Fp_local,nFp,"IVCFp:Fp_local");
  memory->create(Fp_local_plane,nplane,nFp,"IVCFp:Fp_local_plane");
  int iFp;
  double distFp;

  if(!(computeFeFp->invoked_flag & INVOKED_PERATOM)) {
     computeFeFp->compute_peratom();
     computeFeFp->invoked_flag |= INVOKED_PERATOM;
  }
  double **carray_atom = computeFeFp->array_atom;
  int n = nlocal;
//  if(n>100) n = 100;
  double I3[] = {1,0,0,0,1,0,0,0,1};
//  printf("Check 2, in proc %d icompute is %d, size_vector: %d, carray: %f\n",me, icompute,nvalues,carray_atom[0][1]);
  double min_distFp = 0.15*3;
//  if(me<4) printf("in proc %d, in compute_one before loop, check IVC_id %d %d %d %d\n",me,IVC_id[0],IVC_id[1],IVC_id[n-5],IVC_id[n-1]);
  for(i = 0; i < n; i++){
//     if(mask[i] & groupbit){
     //do things here
        distFp = 0.0;
        flag_Fp_yes = 0;
//        IVC_id[i] = -1;
        if(IVC_id[i]>0) continue;
        else{
           // extract the Fp part of FeFp calculation, index 10-18: Fp, 1-9: Fe, 0: i_elastic_bound
           for(iFp=0;iFp<9;iFp++){
             Fp_local[iFp] = carray_atom[i][10+iFp];   
             distFp+= fabs(Fp_local[iFp]-I3[iFp]);
           }
           if(distFp>min_distFp){
              /*
              printf("atom %d in proc %d value: %f\n",i,me,distFp);
              printf("the Fp value is %f %f %f %f %f %f %f %f %f\n",
                      Fp_local[0],Fp_local[1],Fp_local[2],
                      Fp_local[3],Fp_local[4],Fp_local[5],
                      Fp_local[6],Fp_local[7],Fp_local[8]);
              printf("atom %d in proc %d decompose done: %f\n",i,me,Fp_local_plane[0][0]);

              */
              decompose_into_plane(Fp_local,Fp_local_plane,natoms_plane_one,natoms_allslip_one,flag_Fp_yes);
              
              plane_summation(Fp_local_plane,IVC_id[i],onevec);
           }
        }
//     }
  }
  
//  if(me==0)  printf("Check 4.5, before mpi_reduce proc %d val: %f \n",me,onevec[0]);
  for(int m = 0; m < nvalues; m++){
    MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);
  }
  for(int iplane = 0; iplane < nplane ;iplane++){
    MPI_Allreduce(&natoms_plane_one[iplane],&natoms_plane_all[iplane],1,MPI_DOUBLE,MPI_SUM,world);
  }
  MPI_Allreduce(&natoms_allslip_one,&natoms_allslip_all,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&onecount,&allcount,1,MPI_DOUBLE,MPI_SUM,world);
//  printf("Check 5 in proc %d, calculation done %f \n",me,vector[0]);
  if(me==0) printf("total number atoms in IVC is %f\n",allcount);
  memory->destroy(Fp_local);
  memory->destroy(Fp_local_plane);
  delete [] onevec;
  return vector;
}


double ComputeIVCFp::compute_scalar()
{
//  printf("compute IVC_Fp scalar called for proc %d \n",me);
  invoked_scalar = update->ntimestep;
  vector = compute_one();
  double scalar = vector[0];
  return scalar;
}


void ComputeIVCFp::compute_vector()
{
//  printf("compute IVC_Fp vector called for proc %d \n",me);
  invoked_vector = update->ntimestep;
  nplane = 4;
  nFp = 9;
  for(int iplane=0;iplane<nplane;iplane++){
      natoms_plane_one[iplane] = 0.0;
      natoms_plane_all[iplane] = 0.0;
   }
  vector = compute_one();
  int i,j;
  if(me==0){
//    printf("IVC value for 4 planes is called at timstep: %d\n",invoked_vector);
    for(i=0;i<4;i++){ 
      for(j=0;j<9;j++){
         if(i*9+j<nvalues) printf("%f ",vector[i*9+j]);
      }
      printf("\n");
    }
    printf("natoms_plane is %f %f %f %f, tot_slip %f\n",
    natoms_plane_all[0],natoms_plane_all[1],natoms_plane_all[2],natoms_plane_all[3],
    natoms_allslip_all);
  }
//  for (int m = 0; m < nvalues; m++)
//      MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);
//  }
}

void ComputeIVCFp::init()
{
//   printf("compute IVC_Fp init() called for proc %d \n",me);
   int i = 0;
   int iFp,iplane;
   nFp = 9;
   nplane = 4;
   memory->create(Fp_tot_plane,nplane,nFp,"IVCFp:Fp_tot_plane");
   memory->create(Fp_tot,nFp,"IVCFp:Fp_tot");
   memory->create(natoms_plane_one,nplane,"IVCFp:natoms_plane_one");
   memory->create(natoms_plane_all,nplane,"IVCFp:natoms_plane_all");

   for(iFp=0;iFp<9;iFp++){
      for(iplane=0;iplane<4;iplane++){
         Fp_tot_plane[iplane][iFp] = 0.0;
      }
      Fp_tot[iFp] = 0.0;
   }
   natoms_allslip_one = 0.0;
   natoms_allslip_all = 0.0;
   for(iplane=0;iplane<nplane;iplane++){
      natoms_plane_one[iplane] = 0.0;
      natoms_plane_all[iplane] = 0.0;
   }
//   printf("compute IVC_Fp init() over for proc %d \n",me);
}


void ComputeIVCFp::decompose_into_plane(double *Fp_local,double **Fp_local_plane,double *natoms_plane,double natoms_allslip_one, int flag_Fp_yes)
{
  
   double I3[] = {1,0,0,0,1,0,0,0,1};
   int nsystem = 12;
   int i,j,iplane,isystem,iFp;
   int randint = rand() % 100;
//   printf(" in decompose now nplane: %d nFp:%d seed: %d\n",nplane,nFp,randint);
   double maxdF = 1.0;
   double rmin = 0.1;
   double min_angledev = 0.08;
   double min_Fpdot = 0.15;
   double min_distFp = 0.15*3;
   double crit_c = 0.03;
   double check_plane_crit = 0.15;
   double Fp[] = {1,0,0,0,1,0,0,0,1};
   double Fpdot[3];

//   printf(" in decompose check 1\n");
   double s_plane[4][3] = {1,1,1, -1,1,1, 1,-1,1, 1,1,-1};

   double dotplane[4];
   for(iplane=0;iplane<nplane;iplane++){
      dotplane[iplane] = 0.0;
      for(j=0;j<3;j++){
         dotplane[iplane] += s_plane[iplane][j]*s_plane[iplane][j];
      }
      dotplane[iplane] = sqrt(dotplane[iplane]);
   }
   double check_planex[4];
   double check_planey[4];
   double check_planez[4];

//   printf(" in decompose check 2\n");

   double distFp = 0.0;
   int flag_bad_atom = 0;
   int atom_planecount = 0;

   int combination[6];
   int combination_count = 0;
   int icomb = -1;
   int whichplane = -1;
   double Fsep1[] = {1,0,0,0,1,0,0,0,1};
   double Fsep2[] = {1,0,0,0,1,0,0,0,1};
   double Fsep1_ori[9],Fsep2_ori[9];
   double Fpa,Fpb;
   int pindex1 = 0; 
   int pindex2 = 0;
   
   for(iFp=0;iFp<9;iFp++){
      for(iplane=0;iplane<4;iplane++){
         Fp_local_plane[iplane][iFp] = 0.0;
      }
   }

   

   for(iFp=0;iFp<9;iFp++){
      distFp+= fabs(Fp_local[iFp]-I3[iFp]);
   }

   if(distFp>3*min_Fpdot){
      for(iFp=0;iFp<9;iFp++){
         Fp[iFp] = Fp_local[iFp] - I3[iFp];
         if(fabs(Fp[iFp])>maxdF){
            flag_bad_atom = 1;
         }
      } 
      
 
      // only proceed if this is not bad atom
      if(flag_bad_atom == 0){
         natoms_allslip_one =  natoms_allslip_one + 1;
         /*
         printf("in decompose in plane with F: %f %f %f %f %f %f %f %f %f \n",
         Fp[0],Fp[1],Fp[2],
         Fp[3],Fp[4],Fp[5],
         Fp[6],Fp[7],Fp[8]);
         */
         
         for(i=0;i<3;i++){
            Fpdot[i] = sqrt(Fp[i*3]*Fp[i*3]+Fp[i*3+1]*Fp[i*3+1]+Fp[i*3+2]*Fp[i*3+2]);
         }
         /*
         printf("check Fpdot: %f %f %f, Fp: %f %f %f %f %f %f %f %f %f \n",
         Fpdot[0],Fpdot[1],Fpdot[2],
         Fp[0],Fp[1],Fp[2],
         Fp[3],Fp[4],Fp[5],
         Fp[6],Fp[7],Fp[8]);
         */
      // first determine which slip plane this atom belongs
         for(iplane=0;iplane<nplane;iplane++){
            if(fabs(Fpdot[0]) < 0.001){
               check_planex[iplane] = 1.0;
            }else{
               check_planex[iplane] = fabs(s_plane[iplane][0]*Fp[0]+s_plane[iplane][1]*Fp[1]+s_plane[iplane][2]*Fp[2])/(dotplane[iplane]*Fpdot[0]); 
            }
            if(fabs(Fpdot[1]) < 0.001){
               check_planey[iplane] = 1.0;
            }else{
               check_planey[iplane] = fabs(s_plane[iplane][0]*Fp[3]+s_plane[iplane][1]*Fp[4]+s_plane[iplane][2]*Fp[5])/(dotplane[iplane]*Fpdot[1]); 
            }
            if(fabs(Fpdot[2]) < 0.001){
               check_planez[iplane] = 1.0;
            }else{
               check_planez[iplane] = fabs(s_plane[iplane][0]*Fp[6]+s_plane[iplane][1]*Fp[7]+s_plane[iplane][2]*Fp[8])/(dotplane[iplane]*Fpdot[2]); 
            }
            if(((check_planex[iplane]>1-min_angledev)||(Fpdot[0]<min_Fpdot))&&\
               ((check_planey[iplane]>1-min_angledev)||(Fpdot[1]<min_Fpdot))&&\
               ((check_planez[iplane]>1-min_angledev)||(Fpdot[2]<min_Fpdot)))  
            {
               whichplane = iplane;
               atom_planecount = atom_planecount+1;
            }
        
         }
         if(atom_planecount == 1){
            // atom belong to only one slip plane
            natoms_plane[whichplane] = natoms_plane[whichplane]+1;
//            printf("plane %d add one atom to %f\n",whichplane,natoms_plane[whichplane]);
            flag_Fp_yes = 1;
            for(iFp=0;iFp<9;iFp++){
               Fp_local_plane[whichplane][iFp] = Fp_local_plane[whichplane][iFp]+Fp_local[iFp]-I3[iFp];
//             Fp_known[iFp] = Fp_known[iFp]+Fp_local[iFp]-I3[iFp]
            }  
         }else if(atom_planecount == 0){
            // atom belong to non of known slip plane
            whichplane = -1;
         }else{
            // atom belong to more than one slip plane
            whichplane = 4;  
         }
//         printf("which plane is %d\n",whichplane);
         /*
         printf("check planex: %f %f %f %f\n",check_planex[0],check_planex[1],check_planex[2],check_planex[3]);
         printf("check planey: %f %f %f %f\n",check_planey[0],check_planey[1],check_planey[2],check_planey[3]);
         printf("check planez: %f %f %f %f\n",check_planez[0],check_planez[1],check_planez[2],check_planez[3]);
         */
         // try to apply plane decompose on -1 type
         if(whichplane == -1){
            for(i=0;i<6;i++){
               combination[i] = 0;      // 0+1,0+2,0+3,1+2,1+3,2+3
            }
            //0+1 => Fp: 1=2,4=5,7=8 
            if((fabs(Fp[1]-Fp[2])<crit_c)&&(fabs(Fp[4]-Fp[5])<crit_c)&&(fabs(Fp[7]-Fp[8])<crit_c)){
               combination[0] = 1;
               for(iFp=0;iFp<3;iFp++){
                  Fpa = +0.5*Fp[iFp*3+0]+0.25*(Fp[iFp*3+1]+Fp[iFp*3+2]);
                  Fpb = -0.5*Fp[iFp*3+0]+0.25*(Fp[iFp*3+1]+Fp[iFp*3+2]);
                  Fsep1[iFp*3+0] =  Fsep1[iFp*3+1] = Fsep1[iFp*3+2] = Fpa;
                  Fsep2[iFp*3+0] =  Fsep2[iFp*3+1] = Fsep2[iFp*3+2] = Fpa;
                  Fsep2[iFp*3+0] = -Fsep2[iFp*3+0];
               }
            }
            // 0+2 => Fp: 0=2,3=5,6=8
            if((fabs(Fp[0]-Fp[2])<crit_c)&&(fabs(Fp[3]-Fp[5])<crit_c)&&(fabs(Fp[6]-Fp[8])<crit_c)){
               combination[1] = 1;
               for(iFp=0;iFp<3;iFp++){
                   Fpa = +0.5*Fp[iFp*3+1]+0.25*(Fp[iFp*3+0]+Fp[iFp*3+2]);
                   Fpb = -0.5*Fp[iFp*3+1]+0.25*(Fp[iFp*3+0]+Fp[iFp*3+2]);
                   Fsep1[iFp*3+0] = Fsep1[iFp*3+1] = Fsep1[iFp*3+2] = Fpa;
                   Fsep2[iFp*3+0] = Fsep2[iFp*3+1] = Fsep2[iFp*3+2] = Fpb;
                   Fsep2[iFp*3+1] = -Fsep2[iFp*3+1];
               }
            }
            //0+3 => Fp: 0=1,3=4,6=7
            if((fabs(Fp[0]-Fp[1])<crit_c)&&(fabs(Fp[3]-Fp[4])<crit_c)&&(fabs(Fp[6]-Fp[7])<crit_c)){
               combination[2] = 1;
               for(iFp=0;iFp<3;iFp++){
                   Fpa = +0.5*Fp[iFp*3+2]+0.25*(Fp[iFp*3+0]+Fp[iFp*3+1]);
                   Fpb = -0.5*Fp[iFp*3+2]+0.25*(Fp[iFp*3+0]+Fp[iFp*3+1]);
                   Fsep1[iFp*3+0] = Fsep1[iFp*3+1] = Fsep1[iFp*3+2] = Fpa;
                   Fsep2[iFp*3+0] = Fsep2[iFp*3+1] = Fsep2[iFp*3+2] = Fpb;
                   Fsep2[iFp*3+2] = -Fsep2[iFp*3+2];
               }
            }
            //1+2 => Fp: 0=-1,3=-4,6=-7
            if((fabs(Fp[0]+Fp[1])<crit_c)&&(fabs(Fp[3]+Fp[4])<crit_c)&&(fabs(Fp[6]+Fp[7])<crit_c)){
               combination[3] = 1;
               for(iFp=0;iFp<3;iFp++){
                   Fpa = +0.5*Fp[iFp*3+2]+0.25*(Fp[iFp*3+1]-Fp[iFp*3+0]);
                   Fpb = +0.5*Fp[iFp*3+2]-0.25*(Fp[iFp*3+1]+Fp[iFp*3+0]);
                   Fsep1[iFp*3+0] = Fsep1[iFp*3+1] = Fsep1[iFp*3+2] = Fpa;
                   Fsep2[iFp*3+0] = Fsep2[iFp*3+1] = Fsep2[iFp*3+2] = Fpb;
                   Fsep1[iFp*3+0] = -Fsep1[iFp*3+0];
                   Fsep2[iFp*3+1] = -Fsep2[iFp*3+1];
               }
            }
            //1+3 => Fp: 0=-2,3=-5,6=-8
            if((fabs(Fp[0]+Fp[2])<crit_c)&&(fabs(Fp[3]+Fp[5])<crit_c)&&(fabs(Fp[6]+Fp[8])<crit_c)){
               combination[4] = 1;
               for(iFp=0;iFp<3;iFp++){
                   Fpa = +0.5*Fp[iFp*3+1]+0.25*(Fp[iFp*3+2]-Fp[iFp*3+0]);
                   Fpb = +0.5*Fp[iFp*3+1]-0.25*(Fp[iFp*3+2]+Fp[iFp*3+0]);
                   Fsep1[iFp*3+0] = Fsep1[iFp*3+1] = Fsep1[iFp*3+2] = Fpa;
                   Fsep2[iFp*3+0] = Fsep2[iFp*3+1] = Fsep2[iFp*3+2] = Fpb;
                   Fsep1[iFp*3+0] = -Fsep1[iFp*3+0];
                   Fsep2[iFp*3+2] = -Fsep2[iFp*3+2];
               }
            }
            //2+3 => Fp: 1=-2,4=-5,7=-8
            if((fabs(Fp[1]+Fp[2])<crit_c)&&(fabs(Fp[4]+Fp[5])<crit_c)&&(fabs(Fp[7]+Fp[8])<crit_c)){
               combination[5] = 1;
               for(iFp=0;iFp<3;iFp++){
                   Fpa = +0.5*Fp[iFp*3+0]+0.25*(Fp[iFp*3+2]-Fp[iFp*3+1]);
                   Fpb = +0.5*Fp[iFp*3+0]-0.25*(Fp[iFp*3+2]+Fp[iFp*3+1]);
                   Fsep1[iFp*3+0] = Fsep1[iFp*3+1] = Fsep1[iFp*3+2] = Fpa;
                   Fsep2[iFp*3+0] = Fsep2[iFp*3+1] = Fsep2[iFp*3+2] = Fpb;
                   Fsep1[iFp*3+1] = -Fsep1[iFp*3+1];
                   Fsep2[iFp*3+2] = -Fsep2[iFp*3+2];
               }
            }

               
            for(i=0;i<6;i++){
                combination_count += combination[i];
                if(combination[i] == 1){
                   icomb = i;
                }
            }
            if(combination_count == 1){
               flag_Fp_yes = 1;
               whichplane = 10+icomb;
               if(icomb == 0){
                  pindex1 = 0;
                  pindex2 = 1;
               }
               else if(icomb == 1){
                  pindex1 = 0;
                  pindex2 = 2;
               }
               else if(icomb == 2){
                  pindex1 = 0;
                  pindex2 = 3;
               }
               else if(icomb == 3){
                  pindex1 = 1;
                  pindex2 = 2;
               }
               else if(icomb == 4){
                  pindex1 = 1;
                  pindex2 = 3;
               }else if(icomb == 5){
                  pindex1 = 2;
                  pindex2 = 3;
               }

               natoms_plane[pindex1] += 0.5;
               natoms_plane[pindex2] += 0.5;
//               printf("plane %d and %d with icomb %d add one atom to %f and %f\n",pindex1,pindex2,icomb,natoms_plane[pindex1],natoms_plane[pindex2]);
               
               /*** rotate
               for(iFp=0;iFp<0;iFp++){
                   Fsep1[iFp] += I3[iFp];
                   Fsep2[iFp] += I3[iFp];
               }
               Frotate(Fsep1_ori,Fsep1,Rt)
               Frotate(Fsep2_ori,Fsep2,Rt)

               ***/               
               

               for(iFp=0;iFp<9;iFp++){
                   Fsep1_ori[iFp] = Fsep1[iFp];
                   Fsep2_ori[iFp] = Fsep2[iFp];
                   Fp_local_plane[pindex1][iFp] += Fsep1_ori[iFp]-I3[iFp];
                   Fp_local_plane[pindex2][iFp] += Fsep2_ori[iFp]-I3[iFp];
               }
           }

           
           /*
           if(combination_count == 1){
              printf("icomb is %d\n",icomb);
           }
           else{
              printf("combination is %d %d %d %d %d %d\n",
              combination[0],combination[1],
              combination[2],combination[3],
              combination[4],combination[5]);
           }*/
           
         }
      /*
      printf("in docompose got distFp %f  F[1:9] %f %f %f %f %f %f %f %f %f, iplane: %d, icomb:%d\n",distFp,
      Fp[0],Fp[1],Fp[2],
      Fp[3],Fp[4],Fp[5],
      Fp[6],Fp[7],Fp[8],
      whichplane,icomb,flag_Fp_yes);
      */
      }      
      
   }
//   printf("decompose over for seed %d\n",randint); 
}

void ComputeIVCFp::plane_summation(double **Fp_local_plane, int IVC_id_me, double *onevec)
{
//  printf("in plane summation\n");
  int id = -IVC_id_me;
  onecount = 0;
  id = id - 1;     // IVC_id_me start from -1 to -n
  int iFp,iplane;
  if((id>=0)&&(id<nnode_inter))
  {
    for(iplane=0;iplane<4;iplane++)
    {
      for(iFp=0;iFp<9;iFp++)
      {
         onevec[id*9*4+iplane*9+iFp] += Fp_local_plane[iplane][iFp];
      }
    }
    onecount ++;
  }
}

void ComputeIVCFp::Frotate(double *F0, double *Frot, double *R)
{

}


double ComputeIVCFp::memory_usage()
{
  double bytes = maxatom * sizeof(double);
  return bytes;
}
