      INTEGER:: flag_lammps              ! 1 if coupling is enabled, 0 otherwise
      INTEGER,PARAMETER :: fp_lammps = 20
      INTEGER::ptr_lammps
      integer:: mpi_status(MPI_STATUS_SIZE)
      INTEGER:: narg, nprocs_cpfem, nprocs_lammps, comm_split 
      INTEGER:: iproc_split,nprocs_split,flag_job                
      INTEGER:: natoms, i_lammps,j_lammps,n_lammps
      INTEGER:: couple_iter
      integer:: n_couple_iter,n_couple_relax,n_MD_run
      INTEGER:: flag_inside_iter,flag_mdskip,flag_load_status
      INTEGER,PARAMETER :: maxatom_of_one_node = 2500
      REAL(kind=8),ALLOCATABLE:: x_atoms(:),v_atoms(:),f_atoms(:)    ! per_atom array that may be used
      REAL(kind=8),ALLOCATABLE:: zu(:)                               ! per_atom array to correct dispz
      REAL(kind=8),ALLOCATABLE:: fa_ext(:),fa_ext_check(:)
      REAL(kind=8),ALLOCATABLE:: avedisp_x(:),avedisp_y(:),avedisp_z(:)
      real(8):: maxfa_ext,maxfa_ext_check
      REAL(kind=8),PARAMETER  :: epsilon = 0.1
      real(kind=8):: dxlo,dxhi,dylo,dyhi,dzlo,dzhi,ratio_deform
      CHARACTER(len=64)  :: lammps_arg
      CHARACTER(len=1024):: line_lammps
      CHARACTER(len=32)  :: tmpline1,tmpline2,tmpline3,tmpline4,tmpline5
      CHARACTER(len=12)  :: fmt
      real(kind=8)::FEM_zlo,FEM_zhi,MD_zlo,MD_zhi
      real(kind=8)::pin_x,pin_y,pin_z,pall_x,pall_y,pall_z,pin_z0
      real(kind=8)::npt_px,npt_py,npt_pz
      INTEGER:: nnode_inter,id_atom,atom_id     ! # of node at atom interface
      INTEGER:: id_node_interface(maxnode)          ! 1~nnode_interface: exact node id
      INTEGER:: internode_pbc(maxnode)          ! the pbc pair of a given interface node
      INTEGER:: atom_at_node(maxnode)

      !---------------------------------------------------------------------
      ! these variables involved with lammps_scatter into a 32-bit integer
      INTEGER(kind=4):: natoms_IVC(maxnode)            ! # of atoms in cutoff of the interface node
      INTEGER(kind=4):: natoms_SVC(maxnode)            ! # of atoms in cutoff of the interface node
      INTEGER(kind=4):: natoms_FVC(maxnode)	           ! FVC : where deformation gradient is calculated
      integer(kind=4),ALLOCATABLE:: group_tag(:)
      !---------------------------------------------------------------------
 
      INTEGER,ALLOCATABLE::atom_IVC(:,:)                   ! 2d array,  (i,j) ith atom id of node j 
      INTEGER,ALLOCATABLE::atom_SVC(:,:)                   ! 2d array,  (i,j) ith atom id of node j 
      INTEGER,ALLOCATABLE::atom_FVC(:,:)                   ! 2d array,  (i,j) ith atom id of node j 

      REAL(kind=8),allocatable:: weight_atom_IVC(:)		! this is for ATC to get average disp in fem procs
      REAL(kind=8),allocatable:: weight_atom_SVC(:)             ! this is for spread of surface compensation in fem procs
      REAL(kind=8),allocatable:: weight_atom_FVC(:)             ! this is for spread of surface compensation in fem procs
      real(kind=8),allocatable:: dist0(:,:,:),distT(:,:,:)
      real(kind=8),allocatable:: dist0_all(:,:),distT_all(:,:)
      real(kind=8),allocatable::DG_FA(:,:,:),DG_FC(:,:,:),refvector(:,:)
      real(kind=8):: xc_IVC(maxnode*3)
      real(kind=8):: xc_SVC(maxnode*3)
      real(kind=8):: xc_FVC(maxnode*3)
      real(kind=8):: u_CTA(3,maxnode)
      real(kind=8)::MD_boxdim(6),MD_dim0(3)
      real(kind=8)::boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi,xa_tmp(3)
      INTEGER:: flag_inter_node(maxnode),flag_havingLU
      INTEGER:: flag_surf_compen_check,flag_neighlist_read
      real(kind=8)::FEMtime,timepoint(10),CTAtime,ATCtime,Fcalctime
      real(kind=8)::tot_FEMtime,tot_CTAtime,tot_ATCtime,tot_Fcalctime
      real(kind=8)::MDtime,tot_MDtime
      INTEGER:: npadding_atoms,file_status,flag_slu                  ! # of padding atoms < # of tot atoms
      INTEGER,ALLOCATABLE:: padding_atoms(:),elem_of_atom(:)
      REAL(kind=8),ALLOCATABLE::shape_coeff_atom(:)
      REAL(kind=8),PARAMETER:: cutoff = 5.0
      real(kind=8):: thick1,MD_r1,thick2,MD_r2
      integer::LU_flag
      integer::n_checkpoint
      integer::tmp_neigh(12),loopi,loopj
      integer,allocatable:: neighlist(:,:),n_neigh(:)
      REAL(kind=8),PARAMETER:: scale_coeff = 10
      integer:: elem_of_node(maxnodelem,maxnode),n_elem_node(maxnode)
      integer:: nodeid2interid(maxnode)
      real(kind=8):: f_react(mdofx),maxf_react
      character(len=20)::runtime_md

!//-------------following is parameters for self-consistent simulation ------
      integer::flag_SC=0		! 1 if SC process is on
      integer::SC_startstep             ! from which step the self-consistent start to calculate and change the stiffness in media
      real(kind=8)::SC_BOX(6)			! xlo,xhi,ylo,yhi,zlo,zhi  outside of which will be defined by media
      integer:: flag_in_media
      real(kind=8)::sum_stress_FE(6),sum_strain_FE(6)   !sum of stress/strain in each procs in FE, Integral(stress dV)
      real(kind=8)::sum_stress_MD(6),sum_strain_MD(6)   !sum of stress/strain in MD
      real(kind=8)::sum_strain_tot(6),sum_stress_tot(6) !sum of stress/strain by all procs in FE
      real(kind=8)::sum_energy_FE(6),sum_energy_tot(6) !sum of energy
      real(kind=8)::V_FE,V_FE_tot,V_MD,V0_MD                 !sum of volume by FE, MD
      real(kind=8)::ave_stress_FE(6),ave_strain_FE(6)   !sum of FE divided by FE volume
      real(kind=8)::ave_stress_MD(6),ave_strain_MD(6)   !sum of MD divided by MD volume

      real(kind=8)::V_all
      real(kind=8)::sum_strain_all(6),sum_stress_all(6)  !sum of stress/strain by FE+MD
      real(kind=8)::ave_strain_all(6),ave_stress_all(6)  !sum of stress/strain by FE+MD devided by V_all

      real(kind=8)::ave_strain_prev(6),ave_stress_prev(6) !record the value of equilibrium state at last loadstep
      real(kind=8)::ave_strain_incre(6),ave_stress_incre(6) !record the value of equilibrium state at last loadstep

      real(kind=8)::stiff_media(6,6),stiff_media_tan(6,6)                   
      
      real(kind=8)::strain_all,strain_MD,strain_FE,strain_incre
      real(kind=8)::stress_all,stress_MD,stress_FE,stress_incre


! 3 blocks determines interface node, MD region,  Padding Elem bound one by one
! dimension: Padding_Elem > MD_Box > Interface_node
!      real(kind=8)::Inter_node_lo(3),Inter_node_hi(3)            ! to determine the interface node position
!      real(kind=8)::MD_box_lo(3),MD_box_hi(3)                    ! to determine the mobile atom region, outside which will be padding atoms
!      real(kind=8)::Padding_region_lo(3),Padding_region_hi(3)    ! to determine the region having elements having padding elements



      ! this part is for plastic deformation
      integer:: flag_IVCFp,iflag_statusoutput
      REAL(kind=8),allocatable:: data_IVCFp(:)          ! size = nnode_inter*9*nplane, to receive compute_IVC_Fp from lammps
      real(kind=4)::data_print
      
