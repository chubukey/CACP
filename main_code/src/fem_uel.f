      program Couple_Atomistic_Continnum_Program
      
      use superlu_mod   
     
      implicit none

      include 'pardis.h'
      include 'mpif.h'
      include 'ESCM_couple.h'
      parameter(maxnrhs = 1)


!*********Without Constraint************************************
      real(8):: vtau(mdofx),g0xyz(maxcrd),vbc1(mbound1),vbc3(mbound3)
      real(8):: vtau0(mdofx),f(mdofx),vt(mdofx)
      real(8):: vbcs1(mbound1),vbcs2(mbound2),xund(maxcrd) 
      real(8):: svars1(maxsv),svars2(maxsv),trac(mdofx)
      real(8):: promat2d(maxprops,maxgrp),time1(2)
      real(8):: vbc2(mbound2),gpxyz(maxcrd)
      
      integer:: imid(3*maxel),nbc1(mbound1)
      integer:: nbc2(mbound2),idf(mbound1)
      integer:: ijk(mnelx),nbc3(mbound3)
      integer:: ibct1(mdofx),npropsg(maxgrp)
      integer:: ibc1(mdofx),ibct(mdofx)
      integer:: ibc(mdofx)

      integer:: ielemno(maxnodelem*maxnode),nindx(maxnode)
      integer:: ibelem(maxel),ilf(maxel),iface_tet(12)
      integer:: iface_brick(24)
      integer:: nboud1,nboud2,nboud3,nupnch
      integer:: nx,nelx,neq,nplane,node
      integer:: nsvars1,ndofel1,mdload1,npredf1,mcrd1 
      integer:: nelst,nelend,nelnox
      integer:: ierror,ierror2
      integer:: iproc,nproc,nprow,npcol,MPIINFSIZE,id_split
      integer:: mpiinf(mpiinfsizemax),iprocinfo(maxproc)
      integer:: isendorrecv(maxproc)
      integer:: nprocs,icyc,file_length,ngauss,nnode
      integer:: ndf,mdim, nostp,ngroup
      integer:: nstate,i1,i2,ia,ib,itod
      integer:: nseg,i,j,icut,ist_dt,ipoly_write
      integer:: n,neq1, iq

      integer:: ii,n1,n_update,in_st,isize,nstep
      integer:: iflag_over,iwrite_flag
      integer:: idisp_write,maxref,niter_ref
      integer:: jacflag,ibroyflag,ncyc
      integer:: idt_stop,ite,itr,iteref,numupd
      integer:: numrhs,kinc,iflag
      integer(4)::iam
      integer:: iga_id,m
      integer:: nnz_loc,is,ie
      integer:: nrhs,nrow_loc
      integer:: iii,nmpd1,init,nrow,ncol
      integer:: niter,npost1,jj,istv_start3,kk
      integer:: nshift,itag,istv_start
      integer:: istv_start1,istv_start2,islip,ist1
      integer:: ikin_state,ishift1,ishift,ishift2
      integer:: indx2,indx,info,kstep
      
      real(8):: pressb(maxel) !,gtmp(maxel/minproc,maxnp,nslip)
!      real(8)::fptmp(maxel/minproc,maxnp,mphase,3,3)
!      real(8):: chitmp(maxel/minproc,maxnp,nslip)
!      real(8):: gammatmp(maxel/minproc,maxnp,nslip)
      real(8):: solveinitialized
      real(8):: t1,totme,xidt,xmxdt,xmindt
      real(8):: t_nodwell,t_dwell
      real(8):: samp,smin,uamp,umin,tseg,tramp
      real(8):: tperiod,sload,st_rate

      integer:: flag_multirate
      real(8):: list_st_rate(100),rate_last
      integer:: list_rate_steps(100),nrate,stage_rate,nextstep

      real(8):: deltme,tmp,fintme,deltme0,deltme_this
      real(8):: twrite,pnewdt
      real(8):: ectl,xeul(3,maxel)
      real(8):: step,t2,telapsed
      real(8):: xmaxv,w2,w1,berr
      real(8):: colperm,rowperm,grainsize(maxel)
      real(8):: g0val,pnewdt1, b1, b2
      real(8):: b3, c1, c2, c3, dgammadt
      real(8):: dchidt,dfp2,dfp1,deltme_prev
      real(8):: dfp3,dgdt,ect1,ftol,gcval
      
      
      
      
      common/bounda/nboud1,nboud2,nboud3,nupnch
      common/elemnt/nx,nelx,neq,nplane,node
      common/stateupdate/svars2
      common/abq/time1,nsvars1,ndofel1,mdload1,npredf1,svars1,mcrd1
      common/press/pressb,ibelem,ilf,iface_tet,iface_brick
      common/el_par/nelst,nelend,nelnox
******************************************************************
c     Added by Pritam
      integer,allocatable:: ibindx(:)
      integer,allocatable:: irowptr(:)
      integer,allocatable:: ibindx_loc(:)
      integer,allocatable:: irowptr_loc(:)
      real(8),allocatable:: resid(:)
      real(8),allocatable:: val(:)
      real(8),allocatable:: resid_loc(:)
      real(8),allocatable:: var(:), war(:)
      real(8),allocatable:: dvtau(:)
      real(8),allocatable:: dvtau_recv(:), dvtau0(:)
!      real(8),allocatable:: svars_out(:),svars_out_recv(:)
      real(8)::vbc_elim(mdofx)
      common/elimination/vbc_elim
      
      
      
      
!****************************************************************************      
!     Added by Jiahao
!     output part
      integer:: ielcur      

! Added by Jiahao :  GND part 
      integer:: nsize_GND, choose   
      integer:: ishift3, k, grainID
      integer:: elem_grain(maxel)
      parameter (nsize_GND=70)
      real(8):: Fp_dgamma_elem(maxel*(3+9))
      real(8):: Fp_dgamma_node(maxnode*(3+9),ngrains)
      real(8):: Fp_dgam_send(maxel*(3+9))       
      integer:: nodelem(nsize_GND,maxnode),nodeset(4)
      real(8):: A_GND(3,3), A_inv(3,3), det_A, volume(maxel)
      real(8):: xset(4), yset(4),zset(4),totvolume(ngrains)
      real(8):: nodalvalue(ngrains), GND_weight 
      real(8):: xpos,ypos,zpos,distance,center(maxel,3)
      real(4):: r4array(10)
      
!****************************************************************************    


      common/processorinfo/iproc,nprocs,mpiinfsize,mpiinf,
     &                     iprocinfo,isendorrecv
      common/nodeconnec/ielemno,nindx 

      common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin,uamp,
     &umin
      common/load_mono/sload,st_rate
      common/icyc_flag/icyc

      dimension ipoly_write(maxel)
      character*10 fname
      character*50 filname
      allocatable::tseg(:)
      dimension berr(maxnrhs)

!------------Super LU------------------
      integer(superlu_ptr) :: grid
      integer(superlu_ptr) :: options
      integer(superlu_ptr) :: ScalePermstruct
      integer(superlu_ptr) :: LUstruct
      integer(superlu_ptr) :: SOLVEstruct
      integer(superlu_ptr) :: A
      integer(superlu_ptr) :: stat
!--------------------------------------

      
!-----------------------------------------------------------------------
      call MPI_INIT(ierror)
!      write(*,*) 'after init ierror:',ierror
      narg = command_argument_count()
      if(narg/= 2) then
          print *, 'Wrong input argument for mpirun'
          call MPI_ABORT(mpi_comm_world,1,ierror)
      endif 
      
      call mpi_comm_rank(mpi_comm_world,iproc,ierror)
      call mpi_comm_size(mpi_comm_world,nprocs,ierror)
      !-------For Slu Dist--------------------------------

      open(102,file='slu_proc_distrib.inp')
      read(102,*)nprow,npcol
      close(102)

      nprocs_cpfem = nprow*npcol
      !--------------------------------------------------

      call GET_COMMAND_ARGUMENT(1,lammps_arg)
      if(nprocs_cpfem > nprocs) then
          if(iproc == 0) then
              print *,'ERROR:CPFEM cant use more procs than available'
              call mpi_abort(mpi_comm_world,2,ierror)
          endif
      elseif(nprocs_cpfem == nprocs) then
          flag_lammps = 0
      else
          flag_lammps = 1
      endif
      
      nprocs_lammps = nprocs - nprocs_cpfem
      if(iproc < nprocs_cpfem) then
          flag_job = 1    ! flag_job == 1 => running cpfem
          id_split = iproc
      else
          flag_job = 0    ! flag_job == 0 => running lammps
          id_split = iproc-nprocs_cpfem
      endif
!      write(*,*) 'before split iproc',iproc,'flag',flag_job
      comm_split = 0
      call mpi_comm_split(mpi_comm_world,flag_job,0,
     & comm_split,ierror)
      call GET_COMMAND_ARGUMENT(2,lammps_arg)
      OPEN(UNIT=fp_lammps,FILE=lammps_arg,ACTION='READ',STATUS='OLD',
     &     IOSTAT=ierror)
      if(ierror /= 0) then
          print *,'ERROR:Could not open Lammps input'
          call mpi_abort(comm_split,4,ierror)
          call mpi_abort(mpi_comm_world,4,ierror)
      endif
      
  
      iproc_split = 0
      nprocs_split = 0
      call mpi_comm_rank(comm_split,iproc_split,ierror)
      call mpi_comm_size(comm_split,nprocs_split,ierror)
      

      if(flag_job==0) then
         call lammps_open(comm_split,ptr_lammps)
      endif
      
      ! added by Jiaxi, initiate flag_inter_node,id_node_interface,natoms_node
      n_checkpoint = 0
      flag_neighlist_read = 1
      tot_FEMtime = 0.d0
      tot_CTAtime = 0.d0
      tot_ATCtime = 0.d0
      tot_Fcalctime = 0.d0
      tot_MDtime = 0.d0
      ! Inter_node_lo/hi in unit of continuum (nm), so there is 10 factor to AA
 
      flag_slu = 0
      do i_lammps = 1,maxnode
         flag_inter_node(i_lammps) = 0
         id_node_interface(i_lammps) = 0
         xc_IVC(i_lammps*3-2) = 0.0
         xc_IVC(i_lammps*3-1) = 0.0
         xc_IVC(i_lammps*3-0) = 0.0
         xc_SVC(i_lammps*3-2) = 0.0
         xc_SVC(i_lammps*3-1) = 0.0
         xc_SVC(i_lammps*3-0) = 0.0
         xc_FVC(i_lammps*3-2) = 0.0
         xc_FVC(i_lammps*3-1) = 0.0
         xc_FVC(i_lammps*3-0) = 0.0
         natoms_IVC(i_lammps) = 0.0
         natoms_SVC(i_lammps) = 0.0
         natoms_FVC(i_lammps) = 0.0
      enddo
      !       from now on, only procs belong to CPFEM running group
      !--------(flag_CPFEM == 1) will excute the CPFEM part
      if((iproc_split==0).and.(flag_job==0))then
         open(134,position='append',file='fext_data')
      endif

      open(322,file='couple_run.inp')
      read(322,*) tmpline1
      read(322,*) iflag_statusoutput
      read(322,*) tmpline1
      read(322,*) n_couple_iter,n_couple_relax
      read(322,*) tmpline1
      read(322,*) n_MD_run
      read(322,*) tmpline1
      read(322,*) MD_r1,thick1,thick2
      MD_r2 = MD_r1-thick1-thick2
      read(322,*) tmpline1
      read(322,*) SC_startstep
      read(322,*) SC_BOX(1:6)
      close(322)

      nnode_inter = 0
      nrate = 0
      st_rate = 0.0
      list_st_rate(:) = 0.0
      list_rate_steps(:) = 0
      stage_rate = 0
      nextstep = 0

      if(flag_job==1) then         !if #A1

      t1=MPI_WTIME()

      if(nprocs_cpfem.lt.minproc)then
         write(*,*)'nprocs_cpfem less than minproc',nprocs_cpfem,minproc
         stop
      endif

      mpiinfsize=0
      call proc_comb(nprocs_cpfem,mpiinf,mpiinfsize)

      call arrange_proc(iprocinfo,isendorrecv,mpiinf,
     &   mpiinfsize,iproc_split)
     
     
      
!---------------------------------------------------------
      open(201,file='loadtime.inp')
      
      call readstr(201)
      read(201,*)icyc

      call readstr(201)
      read(201,*)fname

      close(201)
            
      file_length=len_trim(fname)
      filname=fname(1:file_length)//'.inp'

      open(unit=lr,file=filname)

      if(iproc_split.eq.0)then


         filname=fname(1:file_length)//'.out'

         open(unit=131,file=filname)
         open(unit=100,file='disp.out')

      endif

      call get_subscript(iproc_split,fname)

      if(icyc.eq.1)then

         filname='gfs'//fname
         open(500,file=filname)
         filname='fpfs'//fname
         open(501,file=filname)
         filname='chifs'//fname
         open(502,file=filname)
         filname='dgdtfs'//fname
         open(503,file=filname)
         filname='dfpdtfs'//fname
         open(504,file=filname)
         filname='dchidtfs'//fname
         open(505,file=filname)
         filname='gammafs'//fname
         open(506,file=filname)
         filname='dgammadtfs'//fname
         open(507,file=filname)
         filname='stressfs'//fname
         open(508,file=filname)
         filname='time.out'
         open(132,file=filname)
         
      else

            filname='stressfs'//fname
            open(500,file=filname)
            
!            filname='fpfs'//fname
!            open(501,file=filname)

!            filname='wpfs'//fname
!            open(502,file=filname)
            
            filname='time.out'
            open(132,file=filname)
            
!            filname='GND'//fname
!            open(510,file=filname)
            
!            filname='totgamma'//fname
!            open(511,file=filname)
            
!            filname='Nu_twin'//fname
!            open(520,file=filname)
            
!            filname='twin_f'//fname
!            open(521,file=filname)
            

      endif


      !-----------if 1 then brick, if 2 then tet-------------

      if(ibrickortet.eq.1)then

         ngauss=8
         nnode=8
         node=8

      elseif(ibrickortet.eq.2)then

         ngauss=1
         nnode=4
         node=4

      endif
!------------------------------------------------------



      pressb(1:nelx)=0.d0
      ibelem(1:nelx)=0
      ilf(1:nelx)=0

      if((iproc_split.eq.0).and.(iflag_statusoutput==1))then
      write(*,*)'start reading'
      endif
      
      call readgm(ndf,g0xyz,imid,ijk,mdim,iproc_split,nprocs_cpfem,
     &            pbc_box)
      if((iproc.eq.0).and.(iflag_statusoutput==1))then
      write(*,*)'gm reading done'
      endif
      if(iproc_split.eq.0)then
      write(*,*)'ndf  nx, nelx', ndf ,nx, nelx
      endif
      
      call readnset(iproc_split)
      call readelset()
      call tmstep(totme,xidt,xmxdt,xmindt,nostp)
      
      if(iproc_split.eq.0)then
      write(*,*)'totme xidt,xmxdt,xmindt,nostp', 
     & totme, xidt,xmxdt,xmindt,nostp
      endif
      
      call readmt(ngroup,promat2d,npropsg,nstate,grainsize,nelx,
     & iproc_split)
      
      ! added by Jiaxi for SC model to initialize Stiff_media
      stiff_media(1,1) = promat2d(1,1)    !c11
      stiff_media(1,2) = promat2d(2,1)    !c12
      stiff_media(1,3) = promat2d(3,1)    !c13
      stiff_media(3,3) = promat2d(4,1)    !c33
      stiff_media(4,4) = promat2d(5,1)    !c55

!------------Jiahao initiation of GND related variavles-------------------
      if(GND_switch == 1) then 
         do n=1,nx*12
           do i=1,ngrains
              Fp_dgamma_node(n,i)=0.0d0
           enddo
         enddo
      endif
!***************************************************************************
     
      
      i1=nelx/nprocs_cpfem+1
      i2=nelx/nprocs_cpfem
      
      ia=nelx-i2*nprocs_cpfem
      ib=nprocs_cpfem-ia
      
      if((iproc_split+1).le.ia)then
         nelst=iproc_split*i1+1
         nelnox=i1
      else
         nelst=ia*i1+(iproc_split-ia)*i2+1
         nelnox=i2
      endif
      
      nelend=nelst+nelnox-1
      
!-------------------------------------------------------
      nnode_inter = 0
      npairs_pbc(1:3) = 0
      st_rate = 0.0


      call initabq(nstate,mdim)
      call readbc(nbc1,vbc1,idf,nbc2,vbc2,nbc3,vbc3,ndf,ibc1,ibct1,
     & nprocs_cpfem,
     & flag_inter_node,nnode_inter,id_node_interface,internode_pbc,
     & npairs_pbc,id_node_pbc)             !modified by Jiaxi to add pbc
      if(npairs_pbc(1)+npairs_pbc(2)+npairs_pbc(3) > 0)then
        iflag_pbc = 1
        iflag_gotw = 0
        if(iproc_split == 0)then 
          write(*,*) 'PBC on, BOX:',pbc_box(1:3)
        endif
      else
        iflag_pbc = 0
      endif

!      write(*,*) 'proc after ', iproc
********************************************************************
c     Added by Pritam
      call node_elem_connec(ijk,ndf)

c     Added by Jiaxi
      call  node_elem_connec_new(ijk,n_elem_node,elem_of_node,
     & maxnode,maxnodelem,iproc,nx,nelx,node,mnelx)
!      write(*,*) 'check node_elem point', iproc
********************************************************************

!--------------Reads traction or displacement flag----------------------

      if(icyc.eq.1)then

         open(201,file='loadtime_cyc.inp')
         
         call readstr(201)
         read(201,*)t_nodwell,t_dwell

         call readstr(201)
         read(201,*)itod 

         if(itod.eq.1)then

            call readstr(201)
            
            
            read(201,*)samp,smin

            uamp=0.d0
            umin=0.d0

         elseif(itod.eq.2)then

            call readstr(201)
            read(201,*)uamp,umin

            samp=0.d0
            smin=0.d0

         endif

         call readstr(201)
         read(201,*)nseg

         allocate(tseg(nseg+1))
         tseg=0.d0
         
         read(201,*)(tseg(i),i=1,nseg+1)
         
         close(201)
            
            
         ! write(787,*)'tseg is',tseg   
            
         if(itod.eq.1)then
            tramp=t_nodwell*smin/samp
         else
            tramp=t_nodwell*umin/uamp
         endif
         
      
         tperiod=2.d0*t_nodwell+t_dwell

         sload=0.d0
         st_rate=0.d0


      else

         open(201,file='loadtime_creep.inp')
      
         call readstr(201)
         read(201,*)tramp

         call readstr(201)
         read(201,*)itod 

         call readstr(201)
         
         if(itod.eq.1)then
            st_rate=0.d0
            read(201,*)sload
         elseif(itod.eq.2)then
            flag_multirate = 0
            sload=0.d0
            read(201,*)st_rate
         elseif(itod.eq.3)then
            flag_multirate = 1
            sload=0.0
            read(201,*) nrate
            read(201,*) list_rate_steps(1:nrate)
            read(201,*) list_st_rate(1:nrate)
         endif
            
         close(201)

         tperiod=0.d0
         t_nodwell=0.d0
         t_dwell=0.d0
         samp=0.d0
         smin=0.d0
         uamp=0.d0
         umin=0.d0

         nseg=0.d0

      endif

      if(itod.eq.2)then
        if(iproc == 0) write(*,*) 'loadcreep done,1 rate=',st_rate
      elseif(itod.eq.3)then
        if(iproc == 0) then
           write(*,*) 'loadcreep done,nrate=',nrate,
     & 'rate =',list_st_rate(1:nrate),
     & 'steps=',list_rate_steps(1:nrate-1)
        endif
      endif
!----------------------------------------------------------------------------------

      deltme=xidt
      icut=0

      ist_dt=0
     
      open(unit=201,file='texture.dat')                   ! read texture

      do i=1,nelx
        read(201,*) xeul(1,i),xeul(2,i),xeul(3,i),tmp
        ipoly_write(i)=tmp
      enddo
        
      close(201)
      
!-------------------------------------------------------------------
      open(unit=202, file='elem_grain.dat')
           
      do i=1,nelx
        read(202,*) elem_grain(i)
!     	elem_grain(i)=1
      enddo
      
      close(202)  
!-------------------------------------------------------------------      
         
      do i=1,mcrd1*nx
         gpxyz(i)=g0xyz(i)
         xund(i)=g0xyz(i)
      enddo

      do i=1,neq
         vtau(i)=0.0d0
         vtau0(i)=0.d0
         vt(i)=0.d0
         trac(i)=0.d0
      enddo
                        


      neq1=neq


      do n=1,nboud1
         if(idf(n).eq.0)then
            do ii=1,ndf
               n1=ndf*(n-1)+ii      
               vbcs1(n1)=vbc1(n1)/totme
               vbc1(n1)=0.d0
            enddo
         endif
      enddo

      do n=1,nboud2
         do ii=1,ndf
            n1=ndf*(n-1)+1
            vbcs2(n1)=vbc2(n1)/totme
         enddo
      enddo
      
      ! write(787,*)'totme is',totme

!*************************************************************************
c     AZTEC routine performing a linear partitioning of the nodes

      call node_partition(N_update,in_st,iproc_split, nprocs_cpfem,
     & nx,ndf)
!      write(*,*) 'N_update is',N_update,'in_st',in_st,'iproc',iproc

      isize=0
      if(iflag_pbc == 0)then
        call init_size_slu(ijk,ndf,iproc_split,N_update,in_st,isize)
      elseif(iflag_pbc == 1)then
        call init_size_slu_pbc(ijk,ndf,iproc_split,N_update,in_st,isize,
     &           npairs_pbc,id_node_pbc)
      endif
      allocate(resid(N_update))
      allocate(irowptr(N_update+1))
      allocate (val(isize))
      allocate (ibindx(isize))

      val=0.d0
      resid=0.d0
      if(iflag_pbc == 0)then
         call init_msr_slu(ijk,ndf,iproc_split,N_update,in_st,
     &           irowptr,ibindx,isize)
      elseif(iflag_pbc == 1)then
         call init_msr_slu_pbc(ijk,ndf,iproc_split,N_update,in_st,
     &           irowptr,ibindx,isize,
     &           npairs_pbc,id_node_pbc)
      endif

      write(*,*) 'after init_msr_slu, flagpbc:',iflag_pbc
!      write(*,*)'N',N_update,'isize',isize,irowptr(:)
!      write(*,*)'ibindx',ibindx(:)
      endif          ! endif # A1
!*************************************************************************
    
! ----- system init complete before running ------------------------------
! ----- now build node-atom connectivity and element-atom connectivity---

      call mpi_barrier(mpi_comm_world,ierror)
      if(flag_lammps == 1)then
        call mpi_bcast(nnode_inter,1,MPI_INTEGER,0,
     & mpi_comm_world, ierror)
        call mpi_bcast(id_node_interface,maxnode,MPI_INTEGER,0,
     & mpi_comm_world, ierror)
      endif
!      call mpi_bcast(st_rate1,1,MPI_DOUBLE_PRECISION,0,
!     & mpi_comm_world, ierror)
!      call mpi_bcast(st_rate2,1,MPI_DOUBLE_PRECISION,0,
!     & mpi_comm_world, ierror)
      if(flag_job==0) then
         n_lammps = 0
         do
           if(iproc_split==0) then
             read(UNIT=fp_lammps,FMT='(A)',IOSTAT=ierror) line_lammps
!           write(*,*) 'iproc is',iproc,iproc_split
             n_lammps = 0
             if(ierror == 0) then
                n_lammps = LEN(TRIM(line_lammps))
                if(n_lammps==0) then
                    line_lammps = ' '
                    n_lammps = 1
                endif
             endif
           endif

           
           call mpi_bcast(n_lammps,1,MPI_INTEGER,0,COMM_SPLIT,ierror)
           if(n_lammps==0)   exit
           call mpi_bcast(line_lammps,n_lammps,MPI_CHARACTER,0,
     &    comm_split,ierror)
!           write(*,*) 'proc',iproc_split,'on line',line_lammps(1:20)
           call lammps_command(ptr_lammps,line_lammps,n_lammps)
         enddo

         close(unit=fp_lammps)
          
!        line_lammps = 'minimize 1e-11 1e-11 100 100'
!        n_lammps = LEN(TRIM(line_lammps))
!        call lammps_command(ptr_lammps,line_lammps,n_lammps)

      endif
!      write(*,*) 'MD relax done',iproc,flag_job
      call mpi_barrier(mpi_comm_world,ierror)
      natoms = 0
      if (flag_job == 0) then
          CALL lammps_get_natoms(ptr_lammps,natoms)
      endif
      call MPI_BARRIER(mpi_comm_world,ierror)

      if(flag_lammps == 1)then
         call mpi_bcast(natoms,1,MPI_INTEGER,nprocs_cpfem,
     &      mpi_comm_world,ierror)
      endif
!      if(iproc == nprocs_cpfem) then
!         call MPI_SEND(natoms,1,MPI_INTEGER,
!     &                   0,100,MPI_COMM_WORLD,ierror)
!      elseif(iproc == 0) then
!         call MPI_RECV(natoms,1,MPI_INTEGER,
!     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
!      endif
!      if(flag_job == 1) then
!      call mpi_bcast(natoms,1,MPI_INTEGER,0,
!     &     comm_split,ierror)
!      endif
!      if(iproc == 0) then
!     write(*,*) 'natoms and iproc',natoms,iproc
!      endif     
 
!-------- allocate necessary arrays for all procs---------------!


      call mpi_barrier(mpi_comm_world,ierror)

! ********* data structure interact with LAMMPS *******!
      if(flag_lammps == 1)then
        if(natoms>20000000)then
          write(*,*) 'error, not enough memory for allocation'
          call mpi_abort(mpi_comm_world,3,ierror)
        else
          ALLOCATE(x_atoms(3*natoms))
          ALLOCATE(zu(natoms))
          ALLOCATE(fa_ext(3*natoms))
          ALLOCATE(avedisp_x(natoms))
          ALLOCATE(avedisp_y(natoms))
          ALLOCATE(avedisp_z(natoms))
          ALLOCATE(group_tag(natoms))
        endif
        if(iproc_split == 0)then
          write(*,*) 'before all IVCFp',flag_job,iproc,nnode_inter
        endif
     
        flag_IVCFp = 0
        if(flag_IVCFp == 1)then
           ALLOCATE(data_IVCFp(nnode_inter*36))
        endif
!********************************************************

        fa_ext(1:3*natoms) = 0.d0
      
        if(flag_job == 1) then
           ALLOCATE(atom_IVC(maxatom_of_one_node,nnode_inter))
           ALLOCATE(atom_SVC(maxatom_of_one_node,nnode_inter))
           ALLOCATE(atom_FVC(maxatom_of_one_node,nnode_inter))
           ALLOCATE(weight_atom_IVC(nnode_inter))
           ALLOCATE(weight_atom_SVC(nnode_inter))
           ALLOCATE(weight_atom_FVC(nnode_inter))


           do i_lammps = 1,nnode_inter
              natoms_IVC(i_lammps) = 0
              natoms_SVC(i_lammps) = 0
              natoms_FVC(i_lammps) = 0
           enddo
           call find_inter_node_index(nnode_inter,
     & id_node_interface,nodeid2interid,nx)
        endif
      endif
!-------- allocate over ---------------------------------------!
       call mpi_barrier(mpi_comm_world,ierror)



!-------- initiate on the output files -----------------------!
      if(flag_lammps == 1)then  ! if B1
      if((flag_job == 1).and.(iproc_split==0))then
         open(283,file='ATC_u.out')
         write(283,*) 'This file record the ATC disp for 
     & each couple iterations'
         close(283)
         open(284,file='ATC_du.out')
         write(284,*) 'This file record the ATC disp change
     & each couple iterations'
         close(284)
!         open(291,file='Stiff_media.out')
!         write(291,*) 'This file record Stiff_media based on S-C
!     & iterations'
!         close(291)

!         open(264,file='Fc_react_column.dat')
!         write(264,*) 'This file record the CTA f_react for
!     & each couple iterations in column'
!         close(264)

         open(265,file='Fc_react_row.out')
         write(265,*) 'This file record the CTA f_react for
     & each couple iterations in row'
         close(265)

         open(266,file='Deformation_Grad_A.out')
         write(266,*) 'This file record the DG_FA data'
         close(266)

         open(267,file='Deformation_Grad_C.out')
         write(267,*) 'This file record the DG_FC data'
         close(267)

         open(319,file ='Ave_FA.out',action='write')
         write(319,*) 'Average DG_FA per step'
         open(320,file ='Ave_FC.out',action='write')
         write(320,*) 'Average DG_FC per step'
         
         close(319)
         close(320)
      endif
      if((flag_job==0).and.(iproc_split==0))then
         open(325,file='IVCFp.out',action='write')
         write(325,*) 'This file records IVCFp value for each IVC cell'
         close(325)
      endif
      endif
!-------------------------------------------------------------!

      call mpi_barrier(mpi_comm_world,ierror)
      if(iproc_split == 0)then
        write(*,*) 'before gather xatoms'
      endif
      if (flag_job == 0) then
           call lammps_gather_atoms(ptr_lammps,'x',2,3,x_atoms,1)
           call lammps_gather_variable
     & (ptr_lammps,'zu','all',2,3,1,1,zu)
           call lammps_extract_variable
     & (ptr_lammps,'pallz','NULL',5,4,1,pall_z)
           call lammps_extract_variable
     & (ptr_lammps,'pin3','NULL',4,4,1,pin_z0)
      endif

      
!-----------  added by Jiaxi, build the atom-node connectivity
!      write(*,*) 'after gather'
      call mpi_barrier(mpi_comm_world,ierror)
      if(flag_lammps == 1)then
      if(iproc == nprocs_cpfem) then
         call MPI_SEND(x_atoms,3*natoms,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
         call MPI_SEND(zu,natoms,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
      elseif(iproc == 0) then
         call MPI_RECV(x_atoms,3*natoms,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
         call MPI_RECV(zu,natoms,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
      endif
      if(flag_job == 1) then             
         call mpi_bcast(x_atoms,3*natoms,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
         call mpi_bcast(zu,natoms,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
      endif
      if(iproc_split == 0) then
         write(*,*) 'before Voronoi' 
      endif
      call mpi_barrier(mpi_comm_world,ierror)
      if(flag_job==1) then ! calculate the connectivity on 1st proc in cpfem_group
!         write(*,*) 'id_n_inter,natoms',id_node_interface,natoms
!         call build_atom_IVC(g0xyz,x_atoms,maxatom_of_one_node,
!     &   nnode_inter,id_node_interface,atom_IVC,natoms,cutoff,
!     &   maxnode,maxcrd,natoms_node,weight_atom_IVC,scale_coeff,
!     &   xc_IVC,atom_at_node)
!          write(*,*) 'start voronoi partition',MD_r,thickness
         call Atom_Voronoi_Partion(g0xyz,x_atoms,zu,maxatom_of_one_node,
     &  nnode_inter,id_node_interface,thick1,thick2,natoms,iproc_split,
     &  natoms_IVC,atom_IVC,xc_IVC,weight_atom_IVC,
     &  natoms_SVC,atom_SVC,xc_SVC,weight_atom_SVC,
     &  natoms_FVC,atom_FVC,xc_FVC,weight_atom_FVC,
     &  maxnode,maxcrd,scale_coeff,MD_r1,MD_r2,group_tag)
!          write(*,*) 'weight check',xc_IVC(1:nnode_inter*3)
      endif
      call mpi_barrier(mpi_comm_world,ierror)
      if(iproc == 0)then
         write(*,*) 'voronoi finished'
!     & ,group_tag(1),group_tag(13232),group_tag(3214),group_tag(422:450)
      endif

       call mpi_barrier(mpi_comm_world,ierror)
       call mpi_bcast(group_tag,natoms,MPI_INTEGER4,0,
     & mpi_comm_world, ierror)
       call mpi_bcast(natoms_IVC,nnode_inter,MPI_INTEGER4,0,
     & mpi_comm_world, ierror)
       call mpi_bcast(natoms_SVC,nnode_inter,MPI_INTEGER4,0,
     & mpi_comm_world, ierror)
       call mpi_bcast(natoms_FVC,nnode_inter,MPI_INTEGER4,0,
     & mpi_comm_world, ierror)

      if (flag_job == 0) then
         call lammps_extract_global(ptr_lammps,'boxxlo',boxxlo,6,1)
         call lammps_extract_global(ptr_lammps,'boxxhi',boxxhi,6,1)
         call lammps_extract_global(ptr_lammps,'boxylo',boxylo,6,1)
         call lammps_extract_global(ptr_lammps,'boxyhi',boxyhi,6,1)
         call lammps_extract_global(ptr_lammps,'boxzlo',boxzlo,6,1)
         call lammps_extract_global(ptr_lammps,'boxzhi',boxzhi,6,1)
         MD_boxdim(1) = boxxlo  
         MD_boxdim(2) = boxxhi
         MD_boxdim(3) = boxylo
         MD_boxdim(4) = boxyhi
         MD_boxdim(5) = boxzlo
         MD_boxdim(6) = boxzhi
         MD_dim0(1) = boxxhi - boxxlo
         MD_dim0(2) = boxyhi - boxylo
         MD_dim0(3) = boxzhi - boxzlo
      endif



      call mpi_barrier(mpi_comm_world,ierror)
      if(iproc == nprocs_cpfem) then
         call MPI_SEND(MD_boxdim,6,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
      elseif(iproc == 0) then
         call MPI_RECV(MD_boxdim,6,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
      endif

      if(flag_job == 1) then             
         call mpi_bcast(MD_boxdim,6,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
      endif


      V0_MD = MD_r1*MD_r1*3.1416*(MD_boxdim(6)-MD_boxdim(5))
     & /(scale_coeff*scale_coeff*scale_coeff)

      call mpi_barrier(mpi_comm_world,ierror)
      if(flag_job==0) then
         call lammps_scatter_atoms(ptr_lammps,'IVC_id',0,1,
     &      group_tag,6,natoms)
         call lammps_scatter_global(ptr_lammps,'natoms_SVC',0,
     & nnode_inter,natoms_SVC,10)
         call lammps_scatter_global(ptr_lammps,'natoms_IVC',0,
     & nnode_inter,natoms_IVC,10)
      endif
      call mpi_barrier(mpi_comm_world,ierror)
            
!      if(iproc_split == 0)then
!         write(*,*) 'before allocate DG_calc variable'
!      endif
      ! notice know that the deforamtion gradient DG_FC:continuum / DG_FA:atomic 
      ! are allocated in different procs, later mpi_send/cast is needed
      allocate(DG_FA(3,3,nnode_inter))
      allocate(DG_FC(3,3,nnode_inter))

      allocate(dist0(3,maxatom_of_one_node,nnode_inter))
      allocate(distT(3,maxatom_of_one_node,nnode_inter))

      allocate(dist0_all(3,natoms))
      allocate(distT_all(3,natoms))
      
      
      if(flag_SC == 1)then
        ave_stress_prev(1:6) = 0.d0      
        ave_strain_prev(1:6) = 0.d0      
        ave_stress_all(1:6) = 0.d0      
        ave_strain_all(1:6) = 0.d0      
        ave_stress_incre(1:6) = 0.d0      
        ave_strain_incre(1:6) = 0.d0      
        sum_strain_all(1:6) = 0.d0
        sum_stress_all(1:6) = 0.d0
        sum_strain_tot(1:6) = 0.0d0
        sum_stress_tot(1:6) = 0.0d0
      endif
      if (flag_job == 0) then
           call lammps_gather_atoms(ptr_lammps,'xu',2,3,x_atoms,1)
      endif
      
      call mpi_barrier(mpi_comm_world,ierror)
      if(iproc == nprocs_cpfem) then
         call MPI_SEND(x_atoms,3*natoms,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
      elseif(iproc == 0) then
         call MPI_RECV(x_atoms,3*natoms,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
      endif

      if(flag_job == 1) then             
         call mpi_bcast(x_atoms,3*natoms,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
      endif

      if(flag_job == 1)then
           if(iproc_split == 0)then
              write(*,*) 'calculate dist0'
           endif
           call f_calc_dist0(nnode_inter,id_node_interface,natoms,
     & dist0,natoms_FVC,atom_FVC,maxnode,g0xyz,x_atoms,
     & maxatom_of_one_node,maxcrd,scale_coeff,xc_FVC,
     & group_tag,dist0_all)

!            write(*,*) 'dist0 check1',iproc,dist0_all(1:3,1:10)
      endif
       
      call mpi_barrier(mpi_comm_world,ierror)

      if(flag_job == 0) then

          flag_neighlist_read = 0
          if((flag_neighlist_read == 1)) then      ! if 1 neighlist.dat will be read for calculating deformation gradient
             if(iproc_split == 0)then
                write(*,*) 'use neigh data'
             endif
          allocate(neighlist(12,natoms))
          allocate(n_neigh(natoms))
      ! fist read the neighborlist from file 'neighlist.dat'
 
          open(unit = 361,file = 'neighlist.dat')
          do loopi = 1,9
             read(361,*)
          enddo
          do loopi = 1,natoms
             read(361,*) atom_id,tmp_neigh(1),tmp_neigh(2),tmp_neigh(3),
     & tmp_neigh(4),tmp_neigh(5),tmp_neigh(6),tmp_neigh(7),tmp_neigh(8),
     & tmp_neigh(9),tmp_neigh(10),tmp_neigh(11),tmp_neigh(12)
             n_neigh(atom_id) = 12
             do loopj = 1,12
                if(tmp_neigh(loopj)==-1) then
                   n_neigh(atom_id) = loopj-1
                   exit
                else
                   neighlist(loopj,atom_id) = tmp_neigh(loopj)
                endif
             enddo
          enddo
          write(*,*) 'neighlist read finish',n_neigh(1),n_neigh(100)
          close(361)
      ! second calculate the initial eta on all node_interface position
!        call f_dist_calc(x_atoms,nnode_inter,id_node_interface,
!     &  neighlist,n_neigh,natoms,maxnode,atom_at_node,dist0)
  
        endif
      endif

      !check if padding atoms is already built, if yes, read files, if not call build_padding_atoms
!     file_status = 0

      
      call mpi_barrier(mpi_comm_world,ierror)
!------------ before loading, check init surf force -------------!
       flag_surf_compen_check = 1
       if(flag_surf_compen_check == 1)then
         if(flag_job == 0) then             
             call lammps_gather_variable
     & (ptr_lammps,'avedispx','all',8,3,1,1,avedisp_x)
             call lammps_gather_variable
     & (ptr_lammps,'avedispy','all',8,3,1,1,avedisp_y)
             call lammps_gather_variable
     & (ptr_lammps,'avedispz','all',8,3,1,1,avedisp_z)
         endif
         call mpi_barrier(mpi_comm_world,ierror)
         if(iproc == nprocs_cpfem) then
           call MPI_SEND(avedisp_x,natoms,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
           call MPI_SEND(avedisp_y,natoms,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
           call MPI_SEND(avedisp_z,natoms,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
         elseif(iproc == 0) then
           call MPI_RECV(avedisp_x,natoms,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
           call MPI_RECV(avedisp_y,natoms,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
           call MPI_RECV(avedisp_z,natoms,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
         endif
         if(flag_job == 1) then             
          call mpi_bcast(avedisp_x,natoms,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
          call mpi_bcast(avedisp_y,natoms,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
          call mpi_bcast(avedisp_z,natoms,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
         endif

         if((flag_job == 1).and.(iproc_split == 0))then
           call Surf_f_check(natoms,avedisp_x,avedisp_y,   !node1 is the current node id
     &  avedisp_z,atom_IVC,maxatom_of_one_node,maxnode,nnode_inter,
     &  weight_atom_IVC,natoms_IVC,scale_coeff,xc_IVC,iproc) 
         endif
         if(iproc_split == 0)then
            write(*,*) 'checking init surf compen force'
         endif
      endif

      endif      ! if B1
!------------ before loading, check init surf force -------------!

      ! here check whether elements contain padding atoms include the interface node

      call mpi_barrier(mpi_comm_world,ierror)
      nstep= 1
      niter= 0
      fintme=0.d0
      iflag_over=0
            
      iwrite_flag=0
      idisp_write=0

      LU_flag = 0 
      if(icyc.eq.2)then
         iwrite_flag=1
         idisp_write=1
      endif


      maxref=24
      niter_ref=14

      jacflag=0
      ibroyflag=0



      if(dabs(tramp-0.d0).lt.1d-10)iwrite_flag=1
      ncyc=0

      if(flag_job==1)then
        allocate(dvtau0(neq),dvtau(neq),dvtau_recv(neq))
        allocate(var(N_update*maxref),war(N_update*maxref))
      endif

      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       deltme0 = deltme
!       dvtau0(1:neq) = 0.0d0

      if(iproc_split == 0) then
          write(*,*) 'start loading'
      endif
      
      do while(iflag_over.eq.0)
        flag_mdskip = 0
        do couple_iter = 1,n_couple_iter+n_couple_relax
           if(couple_iter>n_couple_iter)then
              flag_load_status = 0
           else
              flag_load_status = 1
           endif

           rate_last = st_rate

           if(flag_multirate == 1) then
             if(abs(rate_last) < 1.0e-10)then
               rate_last = list_st_rate(1)
             endif

             if(iproc == 0)then
               if((nstep > nextstep).and.(stage_rate<nrate))then
                 stage_rate = stage_rate+1
                 st_rate = list_st_rate(stage_rate)
                 nextstep = nextstep + list_rate_steps(stage_rate)
               else
                 st_rate = list_st_rate(stage_rate)
               endif
             endif
           endif 
           call mpi_barrier(mpi_comm_world,ierror)
           call mpi_bcast(st_rate,1,MPI_DOUBLE_PRECISION,0,
     &    mpi_comm_world, ierror)
           
           if(iproc == 0)then
        write(*,*) 'step=',nstep,'stage=',stage_rate,'rate=',st_rate,
     & list_st_rate(stage_rate),'last=',rate_last,
     & 'next_step=',nextstep
           endif
           call mpi_barrier(mpi_comm_world,ierror)


           if(iproc_split == 0)then
              write(*,*) 'this is in couple_iter,flagover and iproc',
     &        couple_iter,iflag_over,iproc
           endif
           if(couple_iter == 1) then
              flag_inside_iter = 0
           else
              flag_inside_iter = 1
           endif
        ! if(iproc==0) then
        !   write(*,*) 'this is iteration # ',i_lammps
        !    i_lammps = i_lammps+1
        ! endif
         
           if(flag_job == 1)then 
              if(iproc_split==0)then       
                 write(*,*)'strat increment', nstep  !, iproc
              endif
            
              if(icut.eq.0.and.ibroyflag.eq.0)then
                idt_stop=0
              endif
           endif
         
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

! we start from evaluating the difference in the handshake region between atomic evalutation and continuum evalutation
         if((iproc_split.eq.0).and.(iflag_statusoutput==1))then
            write(*,*) 'before CTA,flag_skip,iproc',flag_mdskip,iproc
         endif
         call cpu_time(timepoint(3))
! here we do the CtA then run MD then AtC two way weak coupling
       
! ----------------CTA part -----------------------
!        if(flag_job == 1)then        
!           call update_MD_boundary(FEM_zlo,FEM_zhi,g0xyz,gpxyz,nx,
!     &        iproc_split)
!           if(iproc_split == 0)then
!              write(*,*) 'update_MD_zb at pos 1 is called' 
!           endif
!         endif
      
      ! 2. do CTA_update on cpfem_proc
      if(flag_mdskip == 0)then
         if((iproc_split.eq.0).and.(iflag_statusoutput==1))then
              write(*,*) 'before CTA_f update',iproc
         endif
!         if((flag_job == 1).AND.(iproc_split ==0))then
         if((flag_job == 1) )then
            call CTA_f_update(f_react,fa_ext,natoms,mdofx,
     &    maxatom_of_one_node,maxnode,weight_atom_IVC,natoms_IVC,
     &    maxcrd,nnode_inter,id_node_interface,atom_IVC,
     &    x_atoms,MD_r1,thick1,natoms_SVC,atom_SVC)
           if(iproc_split==0) then
              open(265,position = 'Append',file='Fc_react_row.out')
!              write(265,*) time(1)
              do i = 1,nnode_inter
                 write(265,*) i,id_node_interface(i),
     &   f_react(id_node_interface(i)*3-2),
     &    f_react(id_node_interface(i)*3-1),
     &   f_react(id_node_interface(i)*3),
     &   g0xyz(id_node_interface(i)*3-2),g0xyz(id_node_interface(i)*3-1)

              enddo
              close(265)
           endif
         endif
         call mpi_barrier(mpi_comm_world,ierror)
      ! old one send the new x_atoms back to lammps_proc
      ! 3. send the the external force accordint to FEM solution back to lammps
      ! need to use scatter_atoms(ptr_lammps,natoms,1,'f',f_atoms)
         if((iproc_split.eq.0).and.(iflag_statusoutput==1))then
            write(*,*) 'before send fa_ext',iproc
!            write(*,*) 'before FE stress/strain',iproc,
!     & sum_stress_tot(1),sum_strain_tot(1),flag_job
!            write(*,*) 'before tot stress/strain',iproc,
!     & sum_stress_all(1),sum_strain_all(1),flag_job
         endif
         if(flag_lammps == 1)then
         if(iproc == 0) then
            call MPI_SEND(fa_ext,3*natoms,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
            call MPI_SEND(FEM_zlo,1,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
            call MPI_SEND(FEM_zhi,1,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
            if(flag_SC == 1)then
                 call MPI_SEND(sum_strain_all,6,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
                 call MPI_SEND(sum_stress_all,6,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
                 call MPI_SEND(V_FE_tot,1,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
            endif
         elseif(iproc == nprocs_cpfem) then
            call MPI_RECV(fa_ext,3*natoms,MPI_DOUBLE_PRECISION,
     &           0,100,MPI_COMM_WORLD,mpi_status,ierror)
            call MPI_RECV(FEM_zlo,1,MPI_DOUBLE_PRECISION,
     &           0,100,MPI_COMM_WORLD,mpi_status,ierror)
            call MPI_RECV(FEM_zhi,1,MPI_DOUBLE_PRECISION,
     &           0,100,MPI_COMM_WORLD,mpi_status,ierror)
            if(flag_SC == 1)then
               call MPI_RECV(sum_strain_all,6,MPI_DOUBLE_PRECISION,
     &              0,100,MPI_COMM_WORLD,mpi_status,ierror)
               call MPI_RECV(sum_stress_all,6,MPI_DOUBLE_PRECISION,
     &              0,100,MPI_COMM_WORLD,mpi_status,ierror)
               call MPI_RECV(V_FE_tot,1,MPI_DOUBLE_PRECISION,
     &              0,100,MPI_COMM_WORLD,mpi_status,ierror)
            endif
         endif
         endif
         call MPI_BARRIER(mpi_comm_world,ierror)
!         if(iproc == 0) then
!            write(*,*) 'before bcast fa_ext'
!         endif

         if(flag_job == 0) then             
           call mpi_bcast(fa_ext,3*natoms,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
           call mpi_bcast(FEM_zlo,1,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
           call mpi_bcast(FEM_zhi,1,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
           if(flag_SC == 1)then
             call mpi_bcast(sum_strain_all,6,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
             call mpi_bcast(sum_stress_all,6,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
             call mpi_bcast(V_FE_tot,1,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
           endif
         endif
         
         call MPI_BARRIER(mpi_comm_world,ierror)
         if(flag_job==0) then
            CALL lammps_scatter_atoms(ptr_lammps,'fext',1,3,fa_ext,4)
         endif
         call MPI_BARRIER(COMM_SPLIT,ierror)


         if((iproc.eq.0).and.(iflag_statusoutput==1))then
            write(*,*) 'CTA done'
         endif
        
!        if(flag_job==0) then
!         CALL lammps_gather_atoms(ptr_lammps,'v',1,3,fa_ext_check,1)
!           maxfa_ext_check = 0.d0
!           do i = 1,natoms
!             maxfa_ext_check = max(maxfa_ext_check,abs(fa_ext_check(i)))
!           enddo
!           if(iproc_split == 0)then
!           write(*,*) 'regather v is',maxfa_ext_check
!           endif
!         endif
! --------------- CTA over -------------------

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call CPU_TIME(timepoint(4))
         if(flag_mdskip == 1)then
            timepoint(4) = timepoint(3)
         endif
!         write(*,*) 'report couuple_iter value',couple_iter,iproc
! --------------- run MD --------------------
         if((iproc.eq.0).and.(iflag_statusoutput==1))then
            write(*,*) 'after CTA before run MD'
         endif
         if(flag_job==0) then

            line_lammps =  'unfix affine'
            n_lammps = LEN(TRIM(line_lammps))
            call lammps_command(ptr_lammps,line_lammps,n_lammps)
            if(flag_load_status == 1)then
               ratio_deform = 1.0
            elseif(flag_load_status == 0)then
               ratio_deform = 0.0
            endif
            dxlo = +st_rate*1.0*MD_boxdim(1)*ratio_deform
            dxhi = +st_rate*1.0*MD_boxdim(2)*ratio_deform
            dylo = +st_rate*1.0*MD_boxdim(3)*ratio_deform
            dyhi = +st_rate*1.0*MD_boxdim(4)*ratio_deform
            fmt = '(F10.7)'
            write(tmpline2,fmt) dxlo
            write(tmpline3,fmt) dxhi
            write(tmpline4,fmt) dylo
            write(tmpline5,fmt) dyhi
            line_lammps = 'fix affine all deform 1 '//
     &  'x delta '//trim(tmpline2)//' '//trim(tmpline3)//
!     &  ' y delta '//trim(tmpline4)//' '//trim(tmpline5)//
     &  ' units box remap x'
            n_lammps = LEN(TRIM(line_lammps))
            call lammps_command(ptr_lammps,line_lammps,n_lammps)

            line_lammps = 'unfix pz'
            n_lammps = LEN(TRIM(line_lammps))
            call lammps_command(ptr_lammps,line_lammps,n_lammps)
!            line_lammps = 'unfix tcontrol'
!            n_lammps = LEN(TRIM(line_lammps))
!            call lammps_command(ptr_lammps,line_lammps,n_lammps)

!            line_lammps = 'fix tcontrol CTA langevin 1 1 0.5 3361'
!            n_lammps = LEN(TRIM(line_lammps))
!            call lammps_command(ptr_lammps,line_lammps,n_lammps)
           
            npt_pz = pall_z + (pin_z-pin_z0)/4
            fmt = '(F10.2)'
            write(tmpline2,fmt) npt_pz
            line_lammps ='fix pz all press/berendsen
     &         z'// tmpline2//' '//tmpline2//'1 modulus 2000'
            n_lammps = LEN(TRIM(line_lammps))
            call lammps_command(ptr_lammps,line_lammps,n_lammps)

            fmt = '(I5)'
            write(tmpline1,fmt) n_MD_run
            line_lammps = 'run '//tmpline1(1:10)
            if(iproc_split == 0) then
               write(*,*) 'run lammps using',line_lammps(1:100)
            endif
           
            n_lammps = LEN(TRIM(line_lammps))
            call lammps_command(ptr_lammps,line_lammps,n_lammps)
         endif
 
 
! -------------  MD finish ------------------
!         write(*,*) 'before gather',iproc
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         if((iproc_split.eq.0).and.(iflag_statusoutput==1))then
            write(*,*) 'MD finished'
         endif

         call CPU_TIME(timepoint(5))
!----------------ATC part ---------------------
         if((iproc_split.eq.0).and.(iflag_statusoutput==1))then
           write(*,*) 'ATC start'
         endif
         if(flag_job == 0) then             
             call lammps_gather_variable
     & (ptr_lammps,'avedispx','all',8,3,1,1,avedisp_x)
             call lammps_gather_variable
     & (ptr_lammps,'avedispy','all',8,3,1,1,avedisp_y)
             call lammps_gather_variable
     & (ptr_lammps,'avedispz','all',8,3,1,1,avedisp_z)
             call lammps_extract_variable
     & (ptr_lammps,'pallz','NULL',5,4,1,pall_z)
             call lammps_extract_variable
     & (ptr_lammps,'pin3','NULL',4,4,1,pin_z)
             call lammps_extract_global(ptr_lammps,'boxxlo',boxxlo,6,1)
             call lammps_extract_global(ptr_lammps,'boxxhi',boxxhi,6,1)
             call lammps_extract_global(ptr_lammps,'boxylo',boxylo,6,1)
             call lammps_extract_global(ptr_lammps,'boxyhi',boxyhi,6,1)
             call lammps_extract_global(ptr_lammps,'boxzlo',boxzlo,6,1)
             call lammps_extract_global(ptr_lammps,'boxzhi',boxzhi,6,1)
            
             pbc_strain(1:3) = 0.0
             pbc_strain(3) = (boxzhi - boxzlo)/MD_dim0(3) - 1
             if(iproc_split == 0)then
               write(*,*) 'got pbc strain',pbc_strain(3)
             endif
             MD_boxdim(1) = boxxlo
             MD_boxdim(2) = boxxhi
             MD_boxdim(3) = boxylo
             MD_boxdim(4) = boxyhi
             MD_boxdim(5) = boxzlo
             MD_boxdim(6) = boxzhi

        ! this part is for Fp transfer, extract compute_IVCFp
      if(flag_IVCFp == 1)then
             call lammps_extract_compute(ptr_lammps,'IVCFp',5,1,
     & nnode_inter*4*9,data_IVCFp)
         
         if(iproc_split == 0)then
            open(325,file='IVCFp.out',action='write',position='append')
            write(325,*) 'ndata: ',nnode_inter*4*9,'nstep: ',nstep
            do i=1,nnode_inter*4*9
               data_print = data_IVCFp(i)
               write(325,*) data_print
            enddo
            close(325)
         endif
      endif

!             call lammps_extract_variable
!     & (ptr_lammps,'avestrain','all',9,3,6,sum_strain_MD)
!             call lammps_extract_variable
!     & (ptr_lammps,'pin1','all',10,3,1,ave_stress_MD(1))
!             call lammps_extract_variable
!     & (ptr_lammps,'pin2','all',10,3,1,ave_stress_MD(2))
!             call lammps_extract_variable
!     & (ptr_lammps,'pin3','all',10,3,1,ave_stress_MD(3))
!             call lammps_extract_variable
!     & (ptr_lammps,'pin4','all',10,3,1,ave_stress_MD(4))
!             call lammps_extract_variable
!     & (ptr_lammps,'pin5','all',10,3,1,ave_stress_MD(5))
!             call lammps_extract_variable
!     & (ptr_lammps,'pin6','all',10,3,1,ave_stress_MD(6))
         endif
         call mpi_barrier(mpi_comm_world,ierror)
         if((iproc_split.eq.0).and.(iflag_statusoutput==1))then
            write(*,*) 'before send',iproc
         endif
         if(flag_lammps == 1)then
           if(iproc == nprocs_cpfem) then
           call MPI_SEND(avedisp_x,natoms,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
           call MPI_SEND(avedisp_y,natoms,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
           call MPI_SEND(avedisp_z,natoms,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
           call MPI_SEND(pbc_strain,3,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
!           call MPI_SEND(sum_strain_MD,6,MPI_DOUBLE_PRECISION,
!     &                   0,100,MPI_COMM_WORLD,ierror)
!           call MPI_SEND(DG_FA,3*3*nnode_inter,MPI_DOUBLE_PRECISION,
!     &                   0,100,MPI_COMM_WORLD,ierror)
           elseif(iproc == 0) then
           call MPI_RECV(avedisp_x,natoms,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
           call MPI_RECV(avedisp_y,natoms,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
           call MPI_RECV(avedisp_z,natoms,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
           call MPI_RECV(pbc_strain,3,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
!           call MPI_RECV(sum_strain_MD,6,MPI_DOUBLE_PRECISION,
!     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
!           call MPI_RECV(DG_FA,3*3*nnode_inter,MPI_DOUBLE_PRECISION,
!     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
           endif
         else
           pbc_strain(1:3) = 0.0     ! in case MD not provide strain
         endif
         if((iproc_split.eq.0).and.(iflag_statusoutput==1))then
            write(*,*) 'before bcast',iproc
         endif

         if(flag_job == 1) then             
          call mpi_bcast(avedisp_x,natoms,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
          call mpi_bcast(avedisp_y,natoms,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
          call mpi_bcast(avedisp_z,natoms,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
          call mpi_bcast(pbc_strain,3,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
!          call mpi_bcast(sum_strain_MD,6,MPI_DOUBLE_PRECISION,
!     &           0,COMM_SPLIT,ierror)          
!          call mpi_bcast(DG_FA,3*3*nnode_inter,MPI_DOUBLE_PRECISION,
!     &           0,COMM_SPLIT,ierror)          
         endif
         if(flag_SC == 1) then
         if(iproc == nprocs_cpfem) then
           call MPI_SEND(ave_stress_MD,6,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
         elseif(iproc == 0) then
           call MPI_RECV(ave_stress_MD,6,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
         endif
         if(flag_job == 1) then             
           call mpi_bcast(ave_stress_MD,6,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
         endif
         endif

         if((iproc_split.eq.0).and.(iflag_statusoutput==1))then
            write(*,*) 'ATC finished'
         endif
      endif         ! flag_mdskip if
          
! -------------- ,ATC over --------------------
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call CPU_TIME(timepoint(6))
      if(flag_mdskip == 1)then
         timepoint(6) = timepoint(5)
      endif

! -------------- compare F ------------------- !
      call mpi_barrier(mpi_comm_world,ierror)
! if using ave_disp, (since this is how xc_IVC is calculated), no need to find x_atoms

      if(flag_lammps == 1)then ! if B2

      if(iproc == nprocs_cpfem) then
!         call MPI_SEND(x_atoms,3*natoms,MPI_DOUBLE_PRECISION,
!     &                   0,100,MPI_COMM_WORLD,ierror)
         call MPI_SEND(MD_boxdim,6,MPI_DOUBLE_PRECISION,
     &                   0,100,MPI_COMM_WORLD,ierror)
      elseif(iproc == 0) then
!         call MPI_RECV(x_atoms,3*natoms,MPI_DOUBLE_PRECISION,
!     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
         call MPI_RECV(MD_boxdim,6,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,mpi_status,ierror)
      endif

      if(flag_job == 1) then             
!         call mpi_bcast(x_atoms,3*natoms,MPI_DOUBLE_PRECISION,
!     &           0,COMM_SPLIT,ierror)          
         call mpi_bcast(MD_boxdim,6,MPI_DOUBLE_PRECISION,
     &           0,COMM_SPLIT,ierror)          
      endif
      if((flag_job == 1).AND.(ibrickortet==2)) then
            call f_calc_distT(nnode_inter,id_node_interface,natoms,
     & dist0,distT,natoms_FVC,atom_FVC,maxnode,gpxyz,
     & avedisp_x,avedisp_y,avedisp_z,
     & maxatom_of_one_node,maxcrd,scale_coeff,xc_FVC,MD_boxdim, 
     & group_tag,dist0_all,distT_all)
!            write(*,*) 'dist0check',iproc,dist0_all(1:3,1:10)
!            write(*,*) 'distTcheck',iproc,distT_all(1:3,1:10)
           call compareF(gpxyz,g0xyz,maxcrd,ijk,nx,nelx,mnelx,
     & nnode_inter,id_node_interface,maxatom_of_one_node,maxnodelem,
     & elem_of_node,n_elem_node,node,DG_FC,iproc_split,fintme,
     & dist0,distT,natoms,DG_FA,natoms_FVC,atom_FVC)
!           if(couple_iter > 5)then
!              call newstiff_LeastSquare(DG_FA,DG_FC,iproc_split,maxcrd,
!     &  g0xyz,nnode_inter,id_node_interface,natoms_IVC,promat2d(1:5,1))
!           endif
           
           if(couple_iter >=2) then
              call f_calc_Atomic_all(dist0_all,distT_all,
     &             ave_strain_MD,natoms,group_tag)
!              if(iproc_split == 0)then
!                write(*,*) 'ave_MD_strain is ',ave_strain_MD(1:2),iproc
!                write(*,*) 'ave_MD_stress is ',ave_stress_MD(1:2),iproc
!              endif
!              write(*,*) 'ave_MD_strain is calculated',ave_strain_MD(:)
           endif
       endif
       endif  ! if B2

!---------------------- SC_calc ------------------------!
       if((flag_SC == 1).and.(couple_iter>1)) then
! -------------- ave_value in both MD and FEM for diagnosis ---------!
         if(iproc == 0)      write(*,*) 'SC start'
        !  MD_proc has all correct value, FEM does not have correct MD value
          call mpi_barrier(mpi_comm_world,ierror)
          if(iproc == 0) then
             call MPI_SEND(ave_strain_MD,6,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
             call MPI_SEND(ave_stress_MD,6,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
             call MPI_SEND(sum_strain_tot,6,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
             call MPI_SEND(sum_stress_tot,6,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
             call MPI_SEND(V_MD,1,MPI_DOUBLE_PRECISION,
     &           nprocs_cpfem,100,MPI_COMM_WORLD,ierror)
          elseif(iproc == nprocs_cpfem) then
             call MPI_RECV(ave_strain_MD,6,MPI_DOUBLE_PRECISION,
     &           0,100,MPI_COMM_WORLD,mpi_status,ierror)
             call MPI_RECV(ave_stress_MD,6,MPI_DOUBLE_PRECISION,
     &           0,100,MPI_COMM_WORLD,mpi_status,ierror)
             call MPI_RECV(sum_strain_tot,6,MPI_DOUBLE_PRECISION,
     &           0,100,MPI_COMM_WORLD,mpi_status,ierror)
             call MPI_RECV(sum_stress_tot,6,MPI_DOUBLE_PRECISION,
     &           0,100,MPI_COMM_WORLD,mpi_status,ierror)
             call MPI_RECV(V_MD,1,MPI_DOUBLE_PRECISION,
     &           0,100,MPI_COMM_WORLD,mpi_status,ierror)
          endif
          if(flag_job == 0) then
             call mpi_bcast(ave_strain_MD,6,MPI_DOUBLE_PRECISION,
     &       0,COMM_SPLIT,ierror)
             call mpi_bcast(ave_stress_MD,6,MPI_DOUBLE_PRECISION,
     &       0,COMM_SPLIT,ierror)
             call mpi_bcast(sum_strain_tot,6,MPI_DOUBLE_PRECISION,
     &       0,COMM_SPLIT,ierror)
             call mpi_bcast(sum_stress_tot,6,MPI_DOUBLE_PRECISION,
     &       0,COMM_SPLIT,ierror)
             call mpi_bcast(V_MD,1,MPI_DOUBLE_PRECISION,
     &       0,COMM_SPLIT,ierror)
          endif
          call mpi_barrier(mpi_comm_world,ierror)
!          if(iproc_split == 0)then
!            write(*,*) 'ave_MD_strain bcasted',iproc,ave_strain_MD(1:3),
!     &      ave_stress_MD(1:2),V0_MD
!            write(*,*) 'sum_FE_strain/stress',iproc,sum_strain_tot(1:2),
!     &      sum_stress_tot(1:2),V_FE_tot
!          endif
!         write(*,*) 'ave_strain_MD is',ave_strain_MD(1:6),iproc,V_MD
!         write(*,*) 'ave_strain_FE is',ave_strain_FE(1:6),iproc,V_FE_tot
!         write(*,*) 'ave_stress_MD is',ave_stress_MD(1:6)/10,iproc,V_MD
!         write(*,*) 'ave_stress_FE is',ave_stress_FE(1:6),iproc,V_FE_tot
!--------care the calculation of V_MD, stress_MD needs to be consistent, need revisit------
!------- also the unit system is different, FEM unit: 1M Pa  MD unit: bar = 0.1 MPa -------
          
!           V_MD = V0_MD*(1+ave_strain_MD(1))
!     &        *(1+ave_strain_MD(2))*(1+ave_strain_MD(3))
          V_MD = V0_MD
!          if(iproc_split == 0)then
!             write(*,*) 'check V first',iproc,V_MD,ave_strain_MD(1:3)
!          endif
          do i_lammps = 1,6
           sum_strain_MD(i_lammps) = ave_strain_MD(i_lammps)*V_MD
           sum_stress_MD(i_lammps) = ave_stress_MD(i_lammps)*V_MD / 10
           ave_strain_FE(i_lammps) = sum_strain_tot(i_lammps)/V_FE_tot
           ave_stress_FE(i_lammps) = sum_stress_tot(i_lammps)/V_FE_tot
          enddo
!         if(iproc_split ==0)then
!         write(*,*) 'ave_strain_MD is',ave_strain_MD(1:6),iproc,V_MD
!         write(*,*) 'ave_strain_FE is',ave_strain_FE(1:6),iproc,V_FE_tot
!         write(*,*) 'ave_stress_MD is',ave_stress_MD(1:6),iproc,V_MD
!         write(*,*) 'ave_stress_FE is',ave_stress_FE(1:6),iproc,V_FE_tot
!          endif

! ------------------------------------------------------------------- !

          V_all = V_MD + V_FE_tot
          if(iproc_split == 0) then
           write(*,*)  'check all value',iproc,V_all,
     & sum_stress_tot(1),sum_stress_MD(1),
     & sum_strain_tot(1),sum_strain_MD(1)
          endif
          call mpi_barrier(mpi_comm_world,ierror)
          do i_lammps = 1,6
             sum_stress_all(i_lammps) = sum_stress_tot(i_lammps) + 
     &                                  sum_stress_MD(i_lammps)
             sum_strain_all(i_lammps) = sum_strain_tot(i_lammps) +
     &                                  sum_strain_MD(i_lammps)

             ave_stress_all(i_lammps) = sum_stress_all(i_lammps)/V_all
             ave_strain_all(i_lammps) = sum_strain_all(i_lammps)/V_all
          
             ave_stress_incre(i_lammps) = ave_stress_all(i_lammps) -
     &                                    ave_stress_prev(i_lammps)
             ave_strain_incre(i_lammps) = ave_strain_all(i_lammps) - 
     &                                    ave_strain_prev(i_lammps)
          enddo

!          flag_SC = 0
!          if(flag_SC == 1) then
!          if(iproc_split == 0) then
!             write(*,*) 'strain part sum_all,V,ave_all,sum_FE,sum_MD',
!     & sum_strain_all(1),V_all,ave_strain_all(1),
!     & sum_strain_tot(1),sum_strain_MD(1)
             
!             write(*,*) 'stress part',
!     & sum_stress_all(1),V_all,ave_stress_all(1)
!             write(*,*) 'V_tot,stress_tot,strain_tot',
!     & V_all,ave_stress_all(1:2),ave_strain_all(1:2)
!             write(*,*) 'ave_stress_all',ave_stress_all(1:2)
!             write(*,*) 'ave_strain_all',ave_strain_all(1:2)
!          endif
          
          call norm_vector(ave_stress_incre,stress_incre,6)
          call norm_vector(ave_strain_incre,strain_incre,6)
          call norm_vector(ave_stress_all,stress_all,6)
          call norm_vector(ave_strain_all,strain_all,6)
 
          if(couple_iter == n_couple_iter) then
            ave_stress_prev(1:6) = ave_stress_all(1:6)
            ave_strain_prev(1:6) = ave_strain_all(1:6)
            if(iproc == 0) then
                write(*,*) 'stress/strain recorded at step:',couple_iter
            endif
          endif
          if(iproc_split == 0) then
             write(*,*) 'check norm',iproc,strain_incre,strain_all
          endif
          if(couple_iter > SC_startstep) then
             open(291,file='Stiff_media.out')
             if(strain_all > 1e-7) then
                 call stiff_media_calc
     &            (ave_stress_all,ave_strain_all,stiff_media)
                if(iproc ==0 ) then
                   write(*,*) 'stiff_media_all',
     &             stiff_media(1,1),stiff_media(1,2)
                   write(291,*) 'stiff_all', stiff_media(1,1),
     &             stiff_media(1,2),stiff_media(4,4)
                endif
             endif
             if(strain_incre > 1e-7) then
                call stiff_media_calc
     &             (ave_stress_incre,ave_strain_incre,stiff_media_tan)
                if(iproc ==0 ) then
                   write(*,*) 'stiff_media_tan',
     &             stiff_media_tan(1,1),stiff_media_tan(1,2)
                   write(291,*) 'stiff_tan',stiff_media_tan(1,1),
     &             stiff_media_tan(1,2),stiff_media_tan(4,4)
                endif
                stiff_media = stiff_media_tan
             endif
             close(291)
           endif
       endif
       
!-------------------------------------------------------!
 
! -------------- compare F over ------------------
       if(iproc == 0)then
          write(*,*) 'Compare F finished',iproc
       endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierror)
 
! -------------- run CPFEM---------------------

       if(flag_job == 1)then  !if # A2
         if(iproc_split == 0) then
           write(*,*) 'CPFEM start' 
         endif
         call cpu_time(timepoint(1))
         if(icyc.eq.1)then
               
           call tdwell(fintme,deltme,xidt,xmindt,iwrite_flag,nseg,tseg) 
         else
               
           call t_creep_cstrainrate(fintme,deltme,xidt,xmindt,xmxdt)

         endif


         if(idt_stop.eq.1)then
            if(iproc.eq.0)then
          write(131,*)'Require smaller increment than specified minimum'
            endif
            stop
         endif

         
!         if(flag_inside_iter == 0)then                  !old way, pulse loading
!            fintme=fintme+deltme
!         endif
!          fintme = fintme+deltme                         !new way gradual loading


         if(flag_load_status == 1)then                   ! 1      < n < n_iter
            deltme_this = deltme
         elseif(flag_load_status == 0)then               ! n_iter <  < n_iter+n_relax
            deltme_this = 0.0
         endif
            

         fintme = fintme +deltme_this
         time1(1)=fintme
         time1(2)=fintme

 


         if(iproc_split.eq.0)then
            write(*,*)'fintme,tramp,deltme',fintme,tramp,deltme_this
         endif
            
            
         twrite=time1(1)-tramp
         

         if(iproc_split.eq.0)then
      !      write(*,*)'Attempting nstep=',nstep
            write(131,*)'Attempting nstep=',nstep
         endif
         
         kstep=1
         pnewdt=1.d0

        ! vtau0(1:neq) = vt(1:neq)+dvtau0(1:neq)*deltme/deltme0
        ! deltme0 = deltme
*****************************************************************
c         kstep=nstep
         
*******************************************************************
         
      ! here we need to cast the arrays needed from iproc_split = 0
      ! to other proc having flag_job = 1


 
!         vtau0(1:neq) = vt(1:neq)+0.49*dvtau0(1:neq)*deltme/deltme0
!         deltme0 = deltme

!         if(iproc_split == 0)then
!             write(*,*) 'check gpxyz before disinc_C',gpxyz(1:6)
!             write(*,*) 'check g0xyz before disinc_C',g0xyz(1:6)
!         endif
          call disinc_C(vtau0,vbcs1,nbc1,deltme_this,idf,nstep,kstep,
     &  g0xyz,gpxyz,time1,ndf,mcrd1,nx,nodeid2interid,
     &  nnode_inter,internode_pbc,maxatom_of_one_node,scale_coeff,
     &  natoms,avedisp_x,avedisp_y,avedisp_z,
     &  pbc_box,pbc_strain,
     &  flag_mdskip,flag_inside_iter,couple_iter,
     &  natoms_IVC,atom_IVC,weight_atom_IVC,xc_IVC,
     &  natoms_FVC,atom_FVC,weight_atom_FVC,xc_FVC)
         call trainc(nbc2,vbc2,vbcs2,deltme,ndf,trac) 
         if((iproc_split == 0).and.(iflag_statusoutput==1))then
         write(*,*)  'after disinc_C'
         write(*,*)  'deltme',deltme,'nstep',nstep,'kstep',kstep
         write(*,*)  'time1',time1(1:2),'couple_iter',couple_iter
         write(*,*) 'check dvtau0'
         write(*,*) '4 xyz',vtau0(10:12)
         write(*,*) '3582 xyz',vtau0(10744:10746)
         endif
!-----------------------------------------qsint---------------------------------------------
       !  ftol=1.d-1
         ftol=1.0d-3
         ite=0
         itr=0
         iteref=0
         numupd=0
         numrhs=0 
         step=1.d0
         ectl=1d-3
         kinc=nstep


         dvtau0(1:neq)=0.d0

         var(1:N_update*maxref)=0.d0
         war(1:N_update*maxref)=0.d0

         vtau(1:neq)=vtau0(1:neq)
         dvtau0(1:neq)=vtau(1:neq)-vt(1:neq)

         iflag=1

         if(jacflag.eq.0)then
            iflag=0
            val=0.d0
            ibroyflag=0

         endif


         f(1:neq)=0.d0
         resid(1:N_update)=0.d0

         do i =1,neq              
           vbc_elim(i) = vtau(i) - vt(i)
         enddo
         
        
         if(iproc_split.eq.0) then        
            write(*,*)'check main program before first loadvc',iflag
         endif
       
!         call loadvc(f,vtau,dvtau0,xund,ijk,imid,npropsg,
         call loadvc_C(f,f_react,vtau,dvtau0,xund,ijk,imid,npropsg,
     & promat2d,deltme,nbc1,ibc1,kinc,ndf,N_update,in_st,irowptr,
     & ibindx,val,resid,isize,iflag,ngauss,pnewdt,elem_grain,xeul,
     & grainsize,Fp_dgamma_node,iproc_split,nprocs_cpfem,comm_split,
     & sum_stress_tot,sum_strain_tot,sum_energy_tot,
     & V_FE,V_FE_tot,stiff_media,SC_BOX,flag_SC,g0xyz,
     & npairs_pbc,id_node_pbc,pbc_strain,pbc_box,iflag_PBC,w_penalty)

!      if (iproc_split.eq.0) then        
!       write(*,*)'check main program -after loadvc'
!      endif


         if(pnewdt.ge.1.d0)then

            xmaxv=dabs(f(1))
            do i=2,neq
               if(dabs(f(i)).gt.xmaxv)then
                  xmaxv=dabs(f(i))
               endif
            enddo
            
            if(iproc_split.eq.0)then
               write(*,*)'kinc=',kinc,'Resid:',xmaxv
            endif

         endif

         !--------------------Super  LU------------------------
         if(iproc_split ==0)then
           write(*,*) 'before start superlu couple_iter is',couple_iter
!     & ,jacflag,(xmaxv.gt.ftol),(couple_iter.lt.n_couple_iter),iproc
         endif

         if((jacflag.eq.0.and.xmaxv.gt.ftol
!     & .or.couple_iter<=2
     & .or.couple_iter.eq.1).AND.LU_flag.eq.0)  then
            if(iproc_split == 0) then
               write(*,*) 'LU_flag is',LU_flag,'create LU'
            endif
            call f_create_gridinfo_handle(grid) 
            call f_create_options_handle(options)  
            call f_create_ScalePerm_handle(ScalePermstruct)
            call f_create_LUstruct_handle(LUstruct)
            call f_create_SOLVEstruct_handle(SOLVEstruct)
            call f_create_SuperMatrix_handle(A)
            call f_create_SuperLUStat_handle(stat)
            
            call f_superlu_gridinit(COMM_SPLIT, nprow, npcol, grid)
            call get_GridInfo(grid, iam=iam)
            LU_flag = 1
            if ( iam >= nprow * npcol ) then
               
               if(iga_id.eq.0)then
                  write(*,*)'iam>=nprow*npcol'
                  
               endif
               
               call mpi_finalize(ierror)
               stop
               
            endif
            

            allocate(ibindx_loc(irowptr(N_update+1)))
            allocate(irowptr_loc(N_update+1))
            flag_slu = 1
            ibindx_loc=0
            irowptr_loc=0
            
            n=neq
            m=neq
            
            nnz_loc=0
            do i=1,N_update
               is=irowptr(i)
               ie=irowptr(i+1)-1
               if(i.eq.N_update)ie=irowptr(i+1)
               do j=is,ie
                  ibindx_loc(j)=ibindx(j)-1
                  nnz_loc=nnz_loc+1
               enddo
            enddo
            
            
            
            do i=1,N_update
               irowptr_loc(i)=irowptr(i)-1
            enddo
            irowptr_loc(N_update+1)=irowptr(N_update+1)
            
            call f_dCreate_CompRowLoc_Mat_dist(A, m, n, nnz_loc,
     &      N_update,in_st,val, ibindx_loc, irowptr_loc, SLU_NR_loc,
     &      SLU_D, SLU_GE)
            
            
            nrhs = 1
                        
            call get_CompRowLoc_Matrix(A, nrow_loc=N_update)

            call f_set_default_options(options)
            call set_superlu_options(options,ColPerm=mmd_at_plus_a)
            call set_superlu_options(options,RowPerm=largediag)
            call get_SuperMatrix(A,nrow=m,ncol=n)
            call f_ScalePermstructInit(m, n, ScalePermstruct)
            call f_LUstructInit(m, n, LUstruct)
            call f_PStatInit(stat)

            
         endif

!------------------------------------------------------
         
         if(xmaxv.gt.ftol)then
                  
            dvtau(1:neq)=0.d0
            
            allocate(resid_loc(N_update))
            
            resid_loc(1:N_update)=resid(1:N_update)
            
            call mpi_barrier(COMM_SPLIT,ierror)
            
            if(jacflag.eq.0)then
            
!            if (iproc_split.eq.0)then   
!            write(*,*)'check bf slu factorization',iproc_split
!            endif
               
            call f_pdgssvx(options, A, ScalePermstruct, resid_loc,
     &    N_update, nrhs, grid, LUstruct, SOLVEstruct, berr, stat, info)
               
               
               
            jacflag=1
         
               
            else
               
               call set_superlu_options(options,Fact=FACTORED)
               
               call f_pdgssvx(options, A, ScalePermstruct, resid_loc,
     &   N_update, nrhs, grid, LUstruct, SOLVEstruct, berr, stat, info)
               
               
            endif
            dvtau(1+in_st:N_update+in_st)=resid_loc(1:N_update)
            
            deallocate(resid_loc)
         
         
         
       call mpi_barrier(COMM_SPLIT,ierror)
       call mpi_allreduce(dvtau,dvtau_recv,neq,
     & mpi_double_precision,
     & mpi_sum,COMM_SPLIT,ierror)
         dvtau(1:neq)=dvtau_recv(1:neq)
            
            
            g0val=0.d0
            gcval=dabs(g0val)

         endif

        

         do while(xmaxv.gt.ftol.and.pnewdt.ge.1.d0)       !<--------------------


          if (iproc.eq.0) then
             write(*,*)'enter qsint, niter',time1(1),iteref
             write(*,*)'residue',xmaxv
          endif
        
        
            iflag=2

              
            ite=ite+1
            itr=itr+1
            step=1.d0
            vtau(1:neq)=vtau(1:neq)+step*dvtau(1:neq)
            dvtau0(1:neq)=vtau(1:neq)-vt(1:neq)
             
!            write(*,*)'before loadvc2 dvtau is',iflag
            f(1:neq)=0.d0
            resid(1:N_update)=0.d0
  
!         call loadvc(f,vtau,dvtau0,xund,ijk,imid,npropsg,
         call loadvc_C(f,f_react,vtau,dvtau0,xund,ijk,imid,npropsg,
     & promat2d,deltme,nbc1,ibc1,kinc,ndf,N_update,in_st,irowptr,
     & ibindx,val,resid,isize,iflag,ngauss,pnewdt,elem_grain,xeul,
     & grainsize,Fp_dgamma_node,iproc_split,nprocs_cpfem,comm_split,
     & sum_stress_tot,sum_strain_tot,sum_energy_tot,
     & V_FE,V_FE_tot,stiff_media,SC_BOX,flag_SC,g0xyz,
     & npairs_pbc,id_node_pbc,pbc_strain,pbc_box,iflag_PBC,w_penalty)

!         if(iproc_split == 0)then
!            maxf_react = 0.d0
!            do i = 1,neq
!               if(abs(f_react(i)) > abs(maxf_react))then
!                  maxf_react = f_react(i)
!               endif
!            enddo
!            write(*,*) 'f_react is',maxf_react,iflag
!            write(*,*) 'f_react is',f_react(1:neq)
!         endif 
     


            if(pnewdt.ge.1.d0)then

               xmaxv=dabs(f(1))
               do i=2,neq
                  if(dabs(f(i)).gt.xmaxv)then
                     xmaxv=dabs(f(i))
                  endif
               enddo

               
               iteref=iteref+1
               
               if(iteref.gt.maxref)then 
                  
                  
                  if(iproc_split.eq.0)then
                     write(*,*)'The broy. iteration did not converge 
     &                even with reformations'
                  endif
                  
                  pnewdt=0.5d0
                  ibroyflag=1
                  
                  
               endif
            

               if(ibroyflag.eq.0)then
            
            
                  dvtau(1:neq)=step*dvtau(1:neq)
                  nmpd1=numupd*N_update
                     
              war(nmpd1+1:nmpd1+N_update)=dvtau(in_st+1:in_st+N_update)
                  
                  allocate(resid_loc(N_update))
                  
                  resid_loc(1:N_update)=resid(1:N_update)
                  
                  call mpi_barrier(COMM_SPLIT,ierror)
                  
                  call set_superlu_options(options,Fact=FACTORED)
     
                  
            call f_pdgssvx(options, A, ScalePermstruct, resid_loc,
     &    N_update, nrhs, grid, LUstruct, SOLVEstruct, berr, stat, info)
     
                  
                  dvtau(1:neq)=0.d0
                  dvtau(1+in_st:N_update+in_st)=resid_loc(1:N_update)
                  
                  deallocate(resid_loc)
                  
         call mpi_barrier(COMM_SPLIT,ierror)
         call mpi_allreduce(dvtau,dvtau_recv,neq,
     &         mpi_double_precision,
     &         mpi_sum,COMM_SPLIT,ierror)
                  dvtau(1:neq)=dvtau_recv(1:neq)
                  
          call broy_parallel(dvtau,var,war,numupd,nmpd1,maxref,in_st,
     &         N_update,neq,step,COMM_SPLIT)
                  
                  
                  g0val=0.d0

                
               endif

            endif
            

         enddo                                                       !<--------------------
         
         
  
            if(iproc_split ==0)then
               write(*,*) 'before end superlu',ibroyflag,jacflag,
     & flag_havingLU,iproc
!    & couple_iter
            endif
         if(ibroyflag.eq.1.and.jacflag.eq.1
!     & .or.couple_iter.eq.1       
     & .or.couple_iter.eq.n_couple_iter)then
            if(iproc_split ==0)then
               write(*,*) 'destroy superLU'
            endif
            deallocate(ibindx_loc,irowptr_loc)
                                ! Deallocate the storage allocated by SuperLU_DIST
            call f_PStatFree(stat)
            call f_Destroy_SuperMat_Store_dist(A)
            call f_ScalePermstructFree(ScalePermstruct)
            call f_Destroy_LU(n, grid, LUstruct)
            call f_LUstructFree(LUstruct)
            call get_superlu_options(options, SolveInitialized=init)
            if (init == YES) then
               call f_dSolveFinalize(options, SOLVEstruct)
            endif
         
         
         
            call f_superlu_gridexit(grid)
            
!     Deallocate the C structures pointed to by the Fortran handles
            call f_destroy_gridinfo_handle(grid)
            call f_destroy_options_handle(options)
            call f_destroy_ScalePerm_handle(ScalePermstruct)
            call f_destroy_LUstruct_handle(LUstruct)
            call f_destroy_SOLVEstruct_handle(SOLVEstruct)
            call f_destroy_SuperMatrix_handle(A)
            call f_destroy_SuperLUStat_handle(stat)
            LU_flag = 0
         else
         
            niter=iteref+1
      
         endif

         if(couple_iter.eq.n_couple_iter)then
            jacflag = 0
         endif    

         


!--------------------------------qsint end-----------------------------------------------------------

         if(pnewdt.lt.1.d0)then

            icut=icut+1

            if(iproc_split.eq.0)then
               write(*,*)'Cutback',icut
            endif

         else
            icut=0
         endif
         
         if(ibroyflag.eq.1)then

            jacflag=0

         endif
         
         
!         write(*,*) 'check point jiaxi 1',iproc
         if(icut.ne.0)then


            if(icut.gt.maxcut)then
               if(iproc_split.eq.0)then
                  write(131,*)'Too many cutbacks. Terminating...'
               endif
               stop
            endif
            
            fintme=fintme-deltme            
            deltme=pnewdt*deltme

            if(deltme.lt.xmindt)then
               if(iproc_split.eq.0)then
                  write(131,*)'Too small increment required',deltme
               endif
               stop
            endif

            if(iproc_split.eq.0)then
               write(131,*)'User Routine Requests smaller increment'
            endif
            
         endif
         
!         write(*,*) 'check point jiaxi 2',iproc
         
         

         if(icut.eq.0.and.ibroyflag.eq.0)then
            

            deltme_prev=deltme
            
            
            
!        write(*,*) 'check point jiaxi 3',iproc
            
            if(iproc_split.eq.0)then
               write(*,*)'step=',nstep,'time=',fintme
               write(*,*)'niter=',niter
            endif

            
            
            if(iproc_split.eq.0)then
          write(131,*)'step=',nstep,'time=',fintme
          write(131,*)'frac completed=',fintme/totme,'time incr=',deltme
          write(131,*)'niter=',niter
            endif



            do iii=1,neq
               gpxyz(iii)=g0xyz(iii)+vtau(iii)
               vt(iii)=vtau(iii)
               vtau0(iii)=vtau(iii)
            enddo
!            if((flag_job == 1).and.(couple_iter>5))then        
!               call update_MD_boundary(FEM_zlo,FEM_zhi,g0xyz,gpxyz,nx,
!     &           iproc_split)
!           write(*,*) 'update_MD_zb at pos 2 is called' 
!            endif

            do ii=1,nsvars1*nelnox*ngauss
               svars1(ii)=svars2(ii)
            enddo



            pnewdt1=(1.d0*niter_ref/(1.d0*niter))**0.5d0


            if(icyc.eq.1)then
               if(pnewdt1.ge.2.d0)pnewdt=2.d0
            else
               if(pnewdt1.gt.pnewdt)pnewdt=pnewdt1
            endif

            
            

            deltme=pnewdt*deltme
            
            

            if(iwrite_flag.eq.1.or.idisp_write.eq.1)then

               idisp_write=1
               if(iproc_split.eq.0)then
                  npost1=nstep-nstep/npost*npost
                  if (npost1.eq.0) then
!                     write(100,*) nstep,fintme,couple_iter
                     do ii=1,nx
                        is=(ii-1)*3
                        write(100,*)vtau(is+1),vtau(is+2),vtau(is+3)
                     enddo
                  end if
               endif


               if(icyc.eq.1)then

                  
                  do ii=1,nelnox
                     do jj=1,ngauss
                    istv_start3=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+210
                        
                        
                        do kk=1,6
                           write(508,*)svars1(istv_start3+kk)
                        enddo
                        
                     enddo
                  enddo

               endif


               if(iproc_split.eq.0)then
                  npost1=nstep-nstep/npost*npost
                  if (npost1.eq.0) then
                  
                     if(icyc.eq.1)then
                        
                        write(132,*)fintme-tramp
                     
                     endif

               endif               

            endif
            
         endif
            
            

            call mpi_barrier(COMM_SPLIT,ierror)



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
           
            
            if((iwrite_flag.eq.1.or.dabs(ncyc*tperiod-twrite).lt.1d-10)
     &.and.icyc.eq.1)then

!               b1=4.d0/3.d0
!               b2=1.d0/3.d0
!               b3=2.d0/3.d0

!               c1=b1/(1.d0+b3)
!               c2=b2/(1.d0+b3)
!               c3=b3/(1.d0+b3)

               

!               iwrite_flag=0

!               ncyc=ncyc+1

!               do ii=1,nelnox
!                  do jj=1,ngauss
                     
!                     nshift=0
                     
!                     do itag=1,ipoly_write(ii-1+nelst)
                        
!                 istv_start=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+nshift
!                        istv_start1=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1
!                        istv_start2=(ii-1)*ngauss+jj
                        
!                        islip=30
!                        ist1=0
!                        if(itag.eq.2)islip=48
!                        if(itag.eq.2)ist1=30
                        
!                        do kk=1,islip
!                           if(ist_dt.eq.1)then
!                        dgdt=(svars1(istv_start+kk)-gtmp(ii,jj,ist1+kk))
!                              write(503,*)dgdt
!                           endif
!                           write(500,*)svars1(istv_start+kk)
!                           gtmp(ii,jj,ist1+kk)=svars1(istv_start+kk)
!                        enddo


                        
!                        do kk=1,islip
!                           if(ist_dt.eq.1)then
!          dgammadt=(svars1(istv_start+islip+kk)-gammatmp(ii,jj,ist1+kk))
!                              write(507,*)dgammadt
!                           endif
                           
!                           write(506,*)svars1(istv_start+islip+kk)
!                     gammatmp(ii,jj,ist1+kk)=svars1(istv_start+islip+kk)
!                        enddo
                        

!                        nshift=2*islip+6

                        
!                        do kk=1,3

!                           indx2=istv_start+nshift+(kk-1)*3

!                           if(ist_dt.eq.1)then

!                           dfp1=(svars1(indx2+1)-fptmp(ii,jj,itag,kk,1))
!                           dfp2=(svars1(indx2+2)-fptmp(ii,jj,itag,kk,2))
!                           dfp3=(svars1(indx2+3)-fptmp(ii,jj,itag,kk,3))
!                              write(504,*)dfp1,dfp2,dfp3


!                           endif
!             write(501,*)svars1(indx2+1),svars1(indx2+2),svars1(indx2+3)
!                           fptmp(ii,jj,itag,kk,1)=svars1(indx2+1)
!                           fptmp(ii,jj,itag,kk,2)=svars1(indx2+2)
!                           fptmp(ii,jj,itag,kk,3)=svars1(indx2+3)
!                        enddo
                        
!                        nshift=nshift+9
                        
!                        ikin_state=220
!                        if(itag.eq.2)ikin_state=250
                        
!                        do kk=1,islip
!                           if(ist_dt.eq.1)then
!                              dchidt=(svars1(istv_start1+ikin_state+kk)
!     &                        -chitmp(ii,jj,ist1+kk))
!                              write(505,*)dchidt
!                           endif
!                           write(502,*)svars1(istv_start1+ikin_state+kk)
!                 chitmp(ii,jj,ist1+kk)=svars1(istv_start1+ikin_state+kk)
!                        enddo                     
!                     enddo                  
!                  enddo
!              enddo

!               ist_dt=1

               


               write(*,*)'Write Output'

            else
            
        
         
!------------------------------------------------------------------
!-----------------------cst or creep ------------------------------

      if (iproc_split.eq.0) then
            write(*,*)'start to Write Output' , iproc_split
      endif
              
               
!----------------- write time.out -------------------------------                
                  npost1=nstep-nstep/npost*npost
                  if (npost1.eq.0) then
                    
                        write(132,*)fintme
              
                  endif
!----------------------------------------------------------------

                  w1=0.93
                  w2=0.07
                  
                  npost1=nstep-nstep/npost*npost
                  if (npost1.eq.0) then
                  write(500,*),nelst,nelend
                  endif
                  
                  do ii=1,nelnox
                     do jj=1,ngauss
                     
                        ishift=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+210

                        npost1=nstep-nstep/npost*npost
                        if (npost1.eq.0) then

                           do kk=1,6
                              
                              indx=indx+1
      !                        write(500,*)svars_out_recv(ishift+kk)
      
                               ielcur=nelst+ii-1
                               write(500,*)svars1(ishift+kk), ielcur
                               
                           enddo
                        endif
                     enddo
                  enddo

        
!-------------------------------------------------------------------                 
                  
                  npost1=nstep-nstep/npost*npost
!                  if (npost1.eq.0) then
!                  write(501,*),nelst,nelend 
!                  endif
                  
                   
                  do ii=1,nelnox
                     do jj=1,ngauss
                        
                        ishift1=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+66
                        ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+177
                        ishift3=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+150             ! Jiahao delta_gamma
                        
                        npost1=nstep-nstep/npost*npost
                        
                        
!                        if (npost1.eq.0) then

!                           do kk=1,9
                              
!                              ielcur=nelst+ii-1
                              
!                              if(ipoly_write(ielcur).eq.1)then
                                 
!                              write(501,*)svars1(ishift1+kk), ielcur
!                              else
                                 
!                             write(501,*)svars1(ishift1+kk), ielcur
!                              endif

!                           enddo

!                        endif
                        
!------------------------------------------------------------------------------                         
!                        do kk=1,9
!                            ielcur=nelst+ii-1                                     ! here we assume only 1 gauss point per elem
!                Fp_dgamma_elem((ielcur-1)*12+kk+3)=svars1(ishift1+kk)  
!                Fp_dgam_send((ielcur-1)*12+kk+3)=svars1(ishift1+kk)
!                        enddo
                        
                            
!                          do kk=1,3    
!                             ielcur=nelst+ii-1   
!                Fp_dgamma_elem((ielcur-1)*12+kk)=svars1(ishift3+kk)
!                 Fp_dgam_send((ielcur-1)*12+kk)=0.0d0        !svars1(ishift3+kk)
!                          enddo                         
!------------------------------------------------------------------------------  
                    
                     enddo
                  enddo   

                    
                  npost1=nstep-nstep/npost*npost
!                  if (npost1.eq.0) then
!                  write(502,*),nelst,nelend
!                  write(510,*),nelst,nelend
!                  write(511,*),nelst,nelend
!                  write(520,*),nelst,nelend
!                  write(521,*),nelst,nelend
!                  endif
                  
!                  do ii=1,nelnox
!                     do jj=1,ngauss
                        
!                        ishift1=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+299
                        
!                        npost1=nstep-nstep/npost*npost
!                        ielcur=nelst+ii-1
!                        if (npost1.eq.0) then
                           
!                  ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+124               ! S_critical
!                            write(502,*)svars1(ishift2), ielcur
!                  ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+125               ! g_size
!                            write(502,*)svars1(ishift2), ielcur
!                  ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+289               ! twin_nu   
!                            write(502,*)svars1(ishift2), ielcur

!                  ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+288                ! tot_GND
!                       write(510,*)svars1(ishift2), ielcur
                       
!                        do kk=1,12
!                  ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+250+kk             ! rho_GND
!                       write(510,*)svars1(ishift2), ielcur
!                        enddo   


!                        do kk=1,12
!                  ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+30+kk               ! tot_gamma
!                       write(511,*)svars1(ishift2), ielcur
!                        enddo  
                      
!                      do kk=1,12
!                 ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+0+kk                ! g_tau     
!                 write(520,*)svars1(ishift2), ielcur
!                 ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+75+kk               ! chi_tau      
!                 write(520,*)svars1(ishift2), ielcur                  
!                       enddo 
                        
!                       do kk=1,1
!                ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+289               ! curl_fp      ! twin   
!                 write(520,*)svars1(ishift2), ielcur
!                       enddo
                        
                 
!                ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+123                   !  mobile_dislocation_tt
!                 write(521,*)svars1(ishift2), ielcur
                        
!                       do kk=1,6
!                 ishift2=(ii-1)*ngauss*nsvars1+(jj-1)*nsvars1+280+kk                ! twin_fraction   
!                  write(521,*)svars1(ishift2), ielcur
!                       enddo
                  
                              
!                        endif

!                     enddo
!                  enddo
                  
            endif                ! if write_flag =1

!*****************************   GND calculation      ********************************  

      if (GND_cal.eq.1) then         

!     First we already find out the element value of delta_gamma and Fp and stored in Fp_dgamma_elem
     

         call mpi_barrier(COMM_SPLIT,ierror)
         call mpi_allreduce(Fp_dgam_send,Fp_dgamma_elem,nelx*ngauss*12,
     &     mpi_double_precision,
     &     mpi_sum,COMM_SPLIT,ierror)


!     next we find the nodal value of delta_gamma, Fp by taking volume average.
!     now calculate the volume of each element  
      ! gpxyz is the updated nodal coordinates 
      ! ijk is connectivity

!  !  find for each node, which element it belongs to 
         do n=1,nx
            do i=1,nsize_GND
               nodelem(i,n)=0
            enddo
         enddo
      
         do n=1,nx
            k=0
            kk=0
            do i=1,nelx
               do j=1,4
                  kk=kk+1
                  if(ijk(kk).eq.n)then
                     k=k+1
                     nodelem(k,n)=i
                  endif
               enddo
            enddo
        
        
           if (k.eq.0)then
              write(*,*)'this node is not connected to any element', n
           endif
        
           if (k.gt.nsize_GND) then
              write(*,*)'maxnodelem is not large enough', k
           endif
        
        enddo
        
                 
                 
! calculte the position of each elem's gauss point
!      if(GND_switch == 1) then        ! if # GND 2


      do k=1,nelx
           nodeset(1)=ijk(k*4-3)
           nodeset(2)=ijk(k*4-2)
           nodeset(3)=ijk(k*4-1)
           nodeset(4)=ijk(k*4)
            do j=1,4
            xset(j)=gpxyz(nodeset(j)*3-2)
            yset(j)=gpxyz(nodeset(j)*3-1)
            zset(j)=gpxyz(nodeset(j)*3)
            enddo

       center(k,1)=1.0/4.0*(xset(1)+xset(2)+xset(3)+xset(4))
       center(k,2)=1.0/4.0*(yset(1)+yset(2)+yset(3)+yset(4))
       center(k,3)=1.0/4.0*(zset(1)+zset(2)+zset(3)+zset(4)) 
          
       enddo    
      
      
      if(iproc_split.eq.0)then
       write(*,*)'finish calculating the volume of each elem'
      endif
      

!  take the first 3 value to be position of gauss point (to calculate x_star)
      
      do k=1,nelx
      Fp_dgamma_elem((k-1)*12+1)=center(k,1)
      Fp_dgamma_elem((k-1)*12+2)=center(k,2)
      Fp_dgamma_elem((k-1)*12+3)=center(k,3)
      enddo
      
      
       
!   choose which variable (from delta_gamma, Fp) you want to calculate for nodal value      
     

      do choose =1,12                                               !  1
      
! calculate nodal value of Fp      

      do n=1,nx
        
        do ii=1,ngrains
        nodalvalue(ii)=0.d0
        totvolume(ii)=0.d0
        enddo
        
        xpos=gpxyz(n*3-2)
        ypos=gpxyz(n*3-1)
        zpos=gpxyz(n*3)
        
        do iq=1,nsize_GND                            ! nsize_GND=50
        
          k=nodelem(iq,n)                  !  k is one of the elements that node i belongs to 

          grainID=elem_grain(k)
          
          if (k.ne.0) then
      
            distance=((xpos-center(k,1))**2.0 
     &  + (ypos-center(k,2))**2.0  +(zpos-center(k,3))**2.0)**0.5  
       
           GND_weight= exp(-2.0*distance)
           
               nodalvalue(grainID)=nodalvalue(grainID) +
     &            GND_weight*Fp_dgamma_elem((k-1)*12+choose)
            
           totvolume(grainID)=totvolume(grainID)+GND_weight
          
          endif
           
        enddo
         
         
         do ii=1,ngrains    
         if (totvolume(ii).gt.0.0) then
         nodalvalue(ii)=nodalvalue(ii)/totvolume(ii)
         Fp_dgamma_node((n-1)*12+choose,ii)=nodalvalue(ii)
         else
         Fp_dgamma_node((n-1)*12+choose,ii)=0.0d0
         endif
         enddo
                      
      enddo   
           
      enddo                                                      !  1
      
      
      if(iproc_split.eq.0)then
       write(*,*)'finish calculating nodal value of Fp'
      endif

      endif              ! endif #GND 2
    
    
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
!-------------------------------------------------------------------------------------------            

            call mpi_barrier(COMM_SPLIT,ierror)
            
            nstep=nstep+1
            if(iproc_split==0)then
               write(*,*) 'current step is',nstep,couple_iter
            endif
            if(nstep.gt.nostp)then
               if(iproc_split==0)then
                  write(*,*) 'step# larger than cutoff',nstep,nostp
               endif
               iflag_over=1
            endif
            if(dabs(fintme-totme).lt.1d-19.or.fintme.gt.totme)then
               if(iproc_split==0)then
                  write(*,*) 'runtime larger than cutoff'
               endif
               iflag_over=2
            endif
          endif
      elseif (flag_job == 0)then
          nstep = nstep + 1
      endif      !endif # A2
     
      call cpu_time(timepoint(2))
      FEMtime = timepoint(2)-timepoint(1)
      ATCtime = timepoint(6)-timepoint(5)
      CTAtime =  timepoint(4)-timepoint(3)
      Fcalctime = timepoint(1) - timepoint(6)
      MDtime = timepoint(5) - timepoint(4)

      if(FEMtime>1e3.or.FEMtime<0) FEMtime = 0
      if(ATCtime>1e3.or.ATCtime<0) ATCtime = 0
      if(CTAtime>1e3.or.CTAtime<0) CTAtime = 0
      if(Fcalctime>1e3.or.Fcalctime<0) Fcalctime = 0
      if(MDtime>1e3.or.MDtime<0) MDtime = 0

      tot_FEMtime = tot_FEMtime + FEMtime
      tot_ATCtime = tot_ATCtime + ATCtime
      tot_CTAtime = tot_CTAtime + CTAtime
      tot_Fcalctime = tot_Fcalctime + Fcalctime
      tot_MDtime = tot_MDtime + MDtime

      
      if(iproc_split < 0)then
         write(*,*) 'FEM time is',FEMtime,tot_FEMtime,iproc
         write(*,*) 'CTA time is',CTAtime,tot_CTAtime,iproc
         write(*,*) 'MD  time is',Fcalctime,tot_Fcalctime,iproc
         write(*,*) 'ATC time is',ATCtime,tot_ATCtime,iproc
         write(*,*) 'Fclac time is',Fcalctime,tot_Fcalctime,iproc
      endif

      if(flag_mdskip == 1) then
         flag_mdskip = 0
      endif
      call mpi_bcast(iflag_over,1,MPI_INTEGER,0,
     &        mpi_comm_world,ierror)
      if(iflag_over>0)then
         stop
      endif         

      enddo
!      write(*,*) 'before enddo iflag_over is',iflag_over
      enddo                                 ! step n, do while(iflag)
      
! %%%%%%%%%%%%%%%%%%%%%%%%%%%  Finish all steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      

      
      if(flag_job==1) then        !if # A3
      t2=MPI_WTIME()
      
      telapsed=(t2-t1)

      deallocate(val,resid,ibindx,irowptr)
!      deallocate(svars_out,svars_out_recv)

      if(icyc.eq.1)then
         deallocate(tseg)
      endif


      deallocate(dvtau0,dvtau,dvtau_recv)
      deallocate(var,war)
      
      if(jacflag.eq.1)then
         deallocate(ibindx_loc,irowptr_loc)

                                ! Deallocate the storage allocated by SuperLU_DIST
         call f_PStatFree(stat)
         call f_Destroy_SuperMat_Store_dist(A)
         call f_ScalePermstructFree(ScalePermstruct)
         call f_Destroy_LU(n, grid, LUstruct)
         call f_LUstructFree(LUstruct)
         call get_superlu_options(options, SolveInitialized=init)
         if (init == YES) then
            call f_dSolveFinalize(options, SOLVEstruct)
         endif
         
         
         
         call f_superlu_gridexit(grid)
         
!     Deallocate the C structures pointed to by the Fortran handles
         call f_destroy_gridinfo_handle(grid)
         call f_destroy_options_handle(options)
         call f_destroy_ScalePerm_handle(ScalePermstruct)
         call f_destroy_LUstruct_handle(LUstruct)
         call f_destroy_SOLVEstruct_handle(SOLVEstruct)
         call f_destroy_SuperMatrix_handle(A)
         call f_destroy_SuperLUStat_handle(stat)
         
         
      endif

      if(iproc_split.eq.0)then
         write(131,*)'Time of operation=',telapsed,'secs'
         write(131,*)'Program Ends'
         close(131)
         close(100)         
      endif

!      write(*,*) 'check point 2 fem proc',iproc,icyc
      if(icyc.eq.1)then

         close(500)
         close(501)
         close(502)
         close(503)
         close(504)
         close(505)
         close(506)
         close(507)
         close(508)

      else

         close(500)
         close(501)
         close(502)

      endif

!     write(*,*) 'check point 3 fem proc',iproc

      if (errorcalornot.eq.1) then
      close(787)
      close(788)
      endif
      

      if(iproc_split.eq.0)then
         close(132)
      endif
      
      endif    ! endif # A3
      
!      write(*,*) 'check point 5 all proc',iproc
      if(flag_job == 1) then
         DEALLOCATE(atom_IVC)
         DEALLOCATE(atom_SVC)
         DEALLOCATE(atom_FVC)
         DEALLOCATE(weight_atom_IVC)
         DEALLOCATE(weight_atom_SVC)
         DEALLOCATE(weight_atom_FVC)
      endif
      deallocate(dist0)
      deallocate(distT)
      deallocate(dist0_all)
      deallocate(distT_all)
      if(flag_job == 0) then
!         DEALLOCATE(dist0)
!         DEALLOCATE(distT)


         if(flag_neighlist_read == 1) then
            DEALLOCATE(n_neigh)
            DEALLOCATE(neighlist) 
         endif
         if(iproc_split == 0) then
           close(361)
         endif
      endif

      DEALLOCATE(fa_ext)
      DEALLOCATE(x_atoms)
      DEALLOCATE(zu)
      DEALLOCATE(avedisp_x)
      DEALLOCATE(avedisp_y)
      DEALLOCATE(avedisp_z)
      if(flag_IVCFp == 1)then
         DEALLOCATE(data_IVCFp)
      endif
!      DEALLOCATE(fa_ext_check)
      if(flag_job == 0) then
         call lammps_close(ptr_lammps)
      endif
      if(iproc.eq.0)then 
        call  PostProcessing(nx,nelx,nprocs_cpfem,nstep,gpxyz,g0xyz,ijk)
      endif
 
      call MPI_FINALIZE(ierror)
      
      
      
      !
      
     

      end  











*********************************************************

      subroutine readgm(ndf,gxyz,imid,ijk,mdim,iproc_split,nprocs_cpfem,
     & pbc_box)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension gxyz(maxcrd),imid(3*maxel)
      dimension ijk(mnelx)
      dimension ijknel(maxeln),xc(3),xmin(3),xmax(3)
      integer::mpiinf(mpiinfsizemax),iprocinfo(maxproc)
      integer::isendorrecv(maxproc)
      integer::iproc_split,nprocs_cpfem
      common/elemnt/nx,nelx,neq,nplane,node
      common/processorinfo/iproc,nprocs,mpiinfsize,mpiinf,
     &                     iprocinfo,isendorrecv

      call readstr(lr)
      read(lr,*)mdim,ndfadd
      ndf=mdim+ndfadd     
      if(ndf.gt.ndofelx)then
         if(iproc_split.eq.0)then
            write(131,*)'Insufficient Mem-ndofelx'
         endif
         stop
      endif
      call readstr(lr)
      read(lr,*) nx
      neq=ndf*nx
      if(ndf.gt.mdofx)then
         if(iproc_split.eq.0)then
            write(131,*)'Insufficient Memory-mdofx'
         endif
         stop
      endif

      if(nx*mdim.gt.maxcrd)then
         write(131,*)"Insufficient mem -maxcrd",nx,ndim,maxcrd
         stop
      endif

      if(mdim.gt.maxdim)then
         write(131,*)'Insufficient mem-maxdim'
         stop
      endif
      
      xlarge = 1.0e6
      xmin(1) = +xlarge 
      xmax(1) = -xlarge
      xmin(2) = +xlarge 
      xmax(2) = -xlarge
      xmin(3) = +xlarge 
      xmax(3) = -xlarge
      zmd = 12*sqrt(1.0)*3.52          ! orientation 4
      zshift = -0.0001/zmd/10*0.0
      z0 = -5.0
      zmesh = 100.0
      do i=1,nx
         read(lr,*) nnum,(xc(ii),ii=1,mdim)
         xc(3) = (xc(3)-z0)*zmd/zmesh+zshift
         do j = 1,3
            if(xc(j)<xmin(j)) then
               xmin(j) = xc(j)
            elseif(xc(j)>xmax(j))then
               xmax(j) = xc(j)
            endif
         enddo
               
         do ii=1,mdim
            gxyz(mdim*(nnum-1)+ii)=xc(ii)
         enddo
      enddo
      pbc_box(1) = xmax(1) - xmin(1)
      pbc_box(2) = xmax(2) - xmin(2)
      pbc_box(3) = xmax(3) - xmin(3)
      if(iproc==0)then
        write(*,*) 'FEM box info'
        write(*,*) xmin(1:3)
        write(*,*) xmax(1:3)
      endif
c     ( read element connectivity )

      call readstr(lr)

      if(ibrickortet.eq.1)then

         nnode=8
         node=8
         ngauss=8

      else
         
         ngauss=1
         nnode=4
         node=4

      endif


      mpe=1
      iep=1
      read(lr,*) nelx
      if(nelx.gt.maxel)then
         if(iproc_split.eq.0)then
            write(131,*)'Insufficient Memory -mnel'
         endif
         stop
      endif
      jj=0
      do n=1,nelx
         read(lr,*) nel,(ijknel(j),j=1,nnode)
         imid(nel)=nnode
         imid(nelx+nel)=mpe
         imid(2*nelx+nel)=iep
         do j=1,nnode
            jj=jj+1
            ijk(jj)=ijknel(j)
         enddo
      enddo
      return
      end

************************************************************      

      subroutine readmt(ngroup,promat2d,npropsg,nstate,grainsize,nelx,
     & iproc_split)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension promat2d(maxprops,maxgrp)
      dimension npropsg(maxgrp),grainsize(maxel)
      integer::mpiinf(mpiinfsizemax),iprocinfo(maxproc)
      integer::isendorrecv(maxproc)
      common/processorinfo/iproc,nprocs,mpiinfsize,mpiinf,
     &                     iprocinfo,isendorrecv

      call readstr(lr)
      read(lr,*)ngroup
      if(ngroup.gt.maxgrp)then
         if(iproc_split.eq.0)then
            write(131,*)'Insufficient Mem -maxgrp'
         endif
         stop
      endif
      
      read(lr,*)(npropsg(n),n=1,ngroup)
      do i=1,ngroup
         if(npropsg(i).gt.maxprops)then
            if(iproc_split.eq.0)then
               write(131,*)'Insufficient Mem - maxprops'
            endif
            stop
         endif
      enddo
      
      do n=1,ngroup
         matgrp=npropsg(n)
         ncount=matgrp/8
         nremain=mod(matgrp,8)
         nshift=0
         if(iproc_split.eq.0)then
            write(131,*)'num of props: ',matgrp
         endif
         do j=1,ncount
            read(lr,*)(promat2d(nshift+i,n),i=1,8)
            nshift=nshift+8
         enddo
         read(lr,*)(promat2d(nshift+i,n),i=1,nremain)
      enddo
      read(lr,*)nstate
      if(nstate.gt.msv)then
         if(iproc_split.eq.0)then
            write(131,*)'Insufficient Memory -msv'
         endif
         stop
      endif
      
      
      
!------------Modified by Jiahao 05/02/2013-----------      
      
      
      npropsg(1)=300
      
      do i=1,npropsg(1)
      promat2d(i,1)=0.0
      enddo
      
      open(201,file='mat1.dat')
      nprops1=209
      ic=0
      nlines=nprops1/8
      do i=1,nlines
         read(201,*)(promat2d(j,1),j=ic+1,ic+8)
         ic=ic+8
      enddo
      
      if(8*nlines.lt.nprops1) then
         read(201,*)(promat2d(j,1), j=8*nlines+1,nprops1)
      endif
      close(201) 
  
!--------------------------------------- 
      open(unit=201,file='mat_singlet.dat')


c     G_hcp------shear modulus 
c     g0_ab-----hardness of basal
c     g0_ap-----hardness of prism
c     g0_bcc---hardness of bcc
c     gi_ab----internal g of basal
c     gi_ap----internal g of prism
c     gi_b1,gi_b2,gi_b3-----internal g of three families of bcc

c     Material constants and hardness
c     primary alpha
      call readstr(201)
      read(201,*)(promat2d(nprops1+i,1), i=1,8)
c     shear modulus and poissons ratio
      call readstr(201)
      read(201,*)G_hcp_a,xnu_hcp_a
c     burger vector
      call readstr(201)
      read(201,*)bur_hcp
      call readstr(201)
      read(201,*)son   ! single or not
      call readstr(201)
      read(201,*)nGNDswitch, sample_size
      close(201)
         
         
      promat2d(nprops1+9,1)=G_hcp_a
      promat2d(nprops1+10,1)=xnu_hcp_a
      promat2d(nprops1+11,1)=bur_hcp
      promat2d(nprops1+12,1)=son
      promat2d(nprops1+13,1)=nGNDswitch
      promat2d(nprops1+14,1)=sample_size                                       !223
         
         
      open(unit=201,file='size.dat')
 
      do ii=1,nelx
        read(201,*)iel,grainsize(ii)
      enddo
      close(201)
         
!      if(iproc == 0)then
!          write(*,*) 'read stiff:',promat2d(1:5,1) 
!      endif
!---------------------------------------------------
      
      return
      end

*******************************************************************

      subroutine disinc(v,vbcs1,nbc1,delt,idf,kinc,kstep,
     &             g0xyz,time,ndf,mcrd,nodeid2interid,
     &  natoms_IVC,i_atom_IVC,maxatom_of_one_node,
     &  avedisp_x,avedisp_y,avedisp_z,natoms,nnode_inter,
     &  weight_atom_IVC,scale_coeff,xc_IVC)

      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension v(mdofx),vbcs1(mbound1),nbc1(mbound1),u(3),time(2)
      dimension coords(3),idf(mbound1),g0xyz(maxcrd)
      dimension xc_IVC(maxnode*3),xn_tmp(3)
      integer::mpiinf(mpiinfsizemax),iprocinfo(maxproc)
      integer::isendorrecv(maxproc),natom,nnode_inter
      real(kind=8)::avedisp_x(natoms)
      real(kind=8)::avedisp_y(natoms)
      real(kind=8)::avedisp_z(natoms)
      INTEGER::i_atom_IVC(maxatom_of_one_node,maxnode)  ! 2d array,  (i,j) ith atom id of node j 
      integer::nodeid2interid(maxnode),index_inter
      integer::natoms_IVC(nnode_inter),maxatom_of_one_node
      common/bounda/nboud1,nboud2,nboud3,nupnch
      common/processorinfo/iproc,nprocs,mpiinfsize,mpiinf,
     &                     iprocinfo,isendorrecv

      ndf1=ndf+1
      do n=1,nboud1
         node1=nbc1(ndf1*(n-1)+1) 
         do  i=1,ndf
            index=nbc1(ndf1*(n-1)+1+i)
            
            if(index.eq.1) then
               nn=ndf*(node1-1)+i
!               nn=ieqno(nn)
               
               if(idf(n).eq.0)then
                  n1=ndf*(n-1)+i
                  v(nn)=v(nn)+vbcs1(n1)*delt
               endif
               
               if(idf(n).eq.1)then
                  icn=mcrd*(node1-1)
                  do ii=1,mcrd
                     coords(ii)=g0xyz(icn+ii)
                  enddo
                  call disp(u,time)
                  v(nn)=u(1)
!                  if(nn<=10) then
!                  write(*,*)'u value is ',u(1),'node', nn,'dof is',i
!                  endif
               endif
               if(idf(n).eq.6)then
                  icn = mcrd*(node1-1)
                  do ii = 1,mcrd
                     coords(ii) = g0xyz(icn+ii)
                  enddo
                  call dispRadical(u,time,coords)
                  if(i<4) then
                    v(nn) = u(i)
                  endif
               endif
               if(idf(n).eq.5)then
                  icn=mcrd*(node1-1)
                  do ii=1,mcrd
                     coords(ii)=g0xyz(icn+ii)
                  enddo
                  call disp2(u,time)
                  v(nn)=u(1)
!                  if(nn<=10) then
!                  write(*,*)'u value is ',u(1),'node', nn,'dof is',i
!                  endif
               endif
!--------------------------------------------------------------------------------          
               if(idf(n).eq.2)then                    !  Jiahao, added user-defined displacement field
                  icn=mcrd*(node1-1)
                  do ii=1,mcrd
                     coords(ii)=g0xyz(icn+ii)
                  enddo
                  call dispX(u,time,coords)                 ! call dispX to define displacement in x 
                  v(nn)=u(1)
               endif          
               
               if(idf(n).eq.3)then                    !  Jiahao, added user-defined displacement field
                  icn=mcrd*(node1-1)
                  do ii=1,mcrd
                     coords(ii)=g0xyz(icn+ii)
                  enddo
                  call dispY(u,time,coords)              ! call dispY to define displacement in x 
                  v(nn)=u(1)
               endif     
               
              if(idf(n)==4) then                    ! Jiaxi, nodal disp b.c. by atoms average
                 icn=mcrd*(node1-1)
                 do ii = 1,mcrd
                    coords(ii) = g0xyz(icn+ii)
                 enddo
! here bcast all arrays needed to all procs with flag_job = 1
                if(i == 1) then 
                 index_inter = nodeid2interid(node1)  !node1 is the current node id
             call ATC_update(u,index_inter,natoms,scale_coeff,
     &  avedisp_x,avedisp_y,avedisp_z,
     &  maxatom_of_one_node,maxnode,iproc_split,
     &  weight_atom_IVC,natoms_IVC,xc_IVC,i_atom_IVC,
     &  weight_atom_FVC,natoms_FVC,xc_FVC)
!                 write(*,*) 'disp by ATC',u(1),u(2),u(3)
!                 write(*,*) 'g0xyz',g0xyz(node1*3-2),g0xyz(node1*3-1),
!     &   g0xyz(node1*3)
                 endif
                 if(i<4) then
                    v(nn) = u(i)
                 endif
              endif  


!----------------------------------------------------------------------
            endif
         enddo
      enddo
      return
      end

*************************************************

      subroutine readstr(lr)
      implicit real*8(a-h,o-z)
      dimension flag(20)
      read(lr,505) (flag(i),i=1,20)
505   format(20a4)
      return
      end
*************************************************

      subroutine initabq(nstate,mdim)
      implicit real*8 (a-h,o-z)
      include 'pardis.h'
      dimension time1(2),svars1(maxsv),svars2(maxsv)
      common/abq/time1,nsvars1,ndofel1,mdload1,npredf1,svars1,mcrd1    
      common/elemnt/nx,nelx,neq,nplane,node
      common/dof/ndof
      common/el_par/nelst,nelend,nelnox
      common/stateupdate/svars2
      
      mcrd1=mdim
      ndof=neq/nx
      numn=4
      if(mcrd1.eq.2)numn=4
      ndofel1=ndof*numn
      mdload1=1
      npredf1=1
      time1(1)=0.d0
      time1(2)=0.d0
      nsvars1=nstate
      do i=1,maxsv
         svars1(i)=0.d0
         svars2(i)=0.d0
      enddo
      return
      end

******************************************************

      subroutine tmstep(totme,xidt,xmxdt,xmindt,nostp)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension flag(20)
      read(lr,505) (flag(i),i=1,20)
      read(lr,*) totme,xidt,xmxdt,xmindt,nostp
505   format(20a4)     
      return
      end

*****************************************************************


      subroutine readbc(nbc1,vbc1,idf,nbc2,vbc2,nbc3,vbc3,ndf,ibc,ibct,
     & nprocs_cpfem,
     & flag_inter_node,nnode_inter,id_node_interface,internode_pbc,           !modified by Jiaxi to calculate nnode_inter
     & npairs_pbc,id_node_pbc)             !modified by Jiaxi to add pbc
      implicit none
      
      include 'pardis.h'
      
      real(8):: vbc1(mbound1),vbc2(mbound2)
      real(8):: vbc3(mbound3)
      real(8):: dk(maxdim),pressb(maxel)
      integer:: nbc3(mbound3),ik(maxdim)
      integer:: nbc1(mbound1),nbc2(mbound2)
      integer:: nnset(maxnset),nlist(maxns,maxnset)
      integer:: idf(mbound1), ibc(mdofx),ibct(mdofx)
      integer:: ipress_set(maxpress),ipress_sf(maxpress)
      integer:: ilf(maxel),iface_tet(12),iface_brick(24)
      integer:: ibelem(maxel)
      integer:: nelset(maxelset),nellist(maxels,maxelset)
      integer:: mpiinf(mpiinfsizemax),iprocinfo(maxproc)
      integer:: isendorrecv(maxproc)
      integer:: ndf,i,ifunc,iend,ii
      integer:: iproc,index,ipf, j,iproc_split
      integer:: istart,iset,mpiinfsize,nconstr
      integer:: n,nbou2,nbn,nbou1
      integer:: nboud1,nboud2,nboud3
      integer:: nconst, nplane, node1,nn,neset
      integer:: ndf1,nel, neq,nelx,node
      integer:: nnod, nodes,nx,nupnch,nset
      integer:: nprocs,npress,nprocs_cpfem
      real(8):: press_val

      integer:: flag_inter_node(maxnode)
      integer:: inter_node_count,nnode_inter
      integer:: id_node_interface(max_internode)
      integer:: internode_pbc(max_internode)

      integer:: flag_pbc_surfx
      integer:: flag_pbc_surfy
      integer:: flag_pbc_surfz
      integer:: id_pbc1,id_pbc2
      integer:: xpbc_count,ypbc_count,zpbc_count,ispbc
      
      common/cbelset/neset,nelset,nellist
      common/elemnt/nx,nelx,neq,nplane,node
      common/bounda/nboud1,nboud2,nboud3,nupnch
      common/cbnset/nset,nnset,nlist
      common/press/pressb,ibelem,ilf,iface_tet,iface_brick
      common/processorinfo/iproc,nprocs,mpiinfsize,mpiinf,
     &                     iprocinfo,isendorrecv
c
c     ( read the 1-st type boundary condition
      if(iproc_split == 0)then
         write(*,*) 'reading bc'
      endif
      call readstr(lr)

      read(lr,*) nbou1
      
      inter_node_count = 0         ! init the inter_node_count

      xpbc_count = 0               ! number of pairs of xpbc
      ypbc_count = 0
      zpbc_count = 0
      flag_pbc_surfx = 0           ! if this is face 1 or face 2
      flag_pbc_surfy = 0
      flag_pbc_surfz = 0
      ndf1=ndf+1                                    ! ndf: # of dof per node
      nboud1=0                                      
      do n=1,nbou1                                  ! nbou1: 4

         ispbc = 0

         read(lr,*) nbn,iset,ifunc,istart,iend
         do i=1,ndf
            ik(i)=0
            dk(i)=0.d0
         enddo
         
         if(ifunc.eq.0)then
            read(lr,*)(dk(i),i=istart,iend)
         endif
         do i=istart,iend
            ik(i)=1
         enddo
         nnod=1
         if(iset.eq.1)nnod=nnset(nbn)
         do i=1,nnod
            nboud1=nboud1+1
            if(nboud1.gt.mbound1)then
               if(iproc_split.eq.0)then
                  write(131,*)'Insufficient Memory-mbound1'
               endif
               stop
            endif
            node1=nbn
            if(iset.eq.1)node1=nlist(i,nbn)
            idf(nboud1)=ifunc
            if(ifunc==4) then
               flag_inter_node(node1) = 1
               inter_node_count = inter_node_count+1
               id_node_interface(inter_node_count) = node1
            elseif(ifunc==11)then
               ispbc = 1
               xpbc_count = xpbc_count+1
               if(flag_pbc_surfx==0)then
                  id_node_pbc(xpbc_count,1,1) = node1
               elseif(flag_pbc_surfx==1)then
                  id_node_pbc(xpbc_count,2,1) = node1
               endif
               if(nnod == xpbc_count) then
                  flag_pbc_surfx = flag_pbc_surfx + 1
                  npairs_pbc(1) = xpbc_count
                  xpbc_count = 0 
               endif
            elseif(ifunc==12)then
               ispbc = 1
               ypbc_count = ypbc_count+1
               if(flag_pbc_surfy==0)then
                  id_node_pbc(ypbc_count,1,2) = node1
               elseif(flag_pbc_surfy==1)then
                  id_node_pbc(ypbc_count,2,2) = node1
               endif
               if(nnod == ypbc_count) then
                  flag_pbc_surfy = flag_pbc_surfy + 1
                  npairs_pbc(2) = ypbc_count
                  ypbc_count = 0
               endif
            elseif(ifunc==13)then
               ispbc = 1
               zpbc_count = zpbc_count+1
               if(flag_pbc_surfz==0)then
                  id_node_pbc(zpbc_count,1,3) = node1
               elseif(flag_pbc_surfz==1)then
                  id_node_pbc(zpbc_count,2,3) = node1
               endif
               if(nnod == zpbc_count) then
                  flag_pbc_surfz = flag_pbc_surfz + 1
                  npairs_pbc(3) = zpbc_count
                  zpbc_count = 0
               endif
            endif
            do ii=1,ndf1 
              if(ispbc == 0)then
               if(ii.eq.1)then
                  nbc1(ndf1*(nboud1-1)+ii)=node1
               else
                  if(ndf1*(nboud1-1)+ii.gt.mbound1)then
                     if(iproc_split.eq.0)then
                        write(131,*)'Insufficient Memeory-mbound1'
                     endif
                     stop
                  endif
                  nbc1(ndf1*(nboud1-1)+ii)=ik(ii-1)
                  if(ifunc.eq.0)vbc1(ndf*(nboud1-1)+ii-1)=dk(ii-1)
               endif
              endif
            enddo
         enddo
      enddo
      nnode_inter = inter_node_count

      internode_pbc(:) = 0
      do n = 1,npairs_pbc(3)
         id_pbc1 = id_node_pbc(n,1,3) 
         id_pbc2 = id_node_pbc(n,2,3) 
         do nn = 1,nnode_inter
            if(id_node_interface(nn) == id_pbc1)then
                internode_pbc(nn) = id_pbc2
            elseif(id_node_interface(nn) == id_pbc2)then
                internode_pbc(nn) = id_pbc1
            endif
         enddo
      enddo
!      do n = 1,npairs_pbc(3)
!        write(*,*) 'n',n,id_node_pbc(n,1:2,3)
!      enddo
      if(iproc == 0)then
        write(*,*) 'n_zpbc:',npairs_pbc(3),'nnode_inter',nnode_inter
!        do n = 1,nnode_inter 
!           write(*,*) id_node_interface(n),'my pair is',internode_pbc(n)
!        enddo
      endif

      do n=1,neq
         ibc(n)=0
      enddo
      do n=1,nboud1
         nodes=nbc1(ndf1*(n-1)+1)
         do i=1,ndf
            index=nbc1(ndf1*(n-1)+1+i)
            if(index.eq.1)then
               nn=ndf*(nodes-1)+i
               ibc(nn)=1
            endif
         enddo         
      enddo
c
c     ( read the 2-nd type boundary condition : point forces )
c
      call readstr(lr)      
      read(lr,*) nbou2
      nboud2=0
      do n=1,nbou2
         read(lr,*) nbn,iset,ifunc,istart,iend       ! nbn:node/nodeset #,
                                                     ! iset:nbn is node or nodeset
                                                     ! ifunc: 0 if value is given, 1234 if through func

         do i=1,ndf
            ik(i)=0
            dk(i)=0.d0
         enddo

         if(ifunc.eq.0)then
            read(lr,*)(dk(i),i=istart,iend)
         endif
         do i=istart,iend
            ik(i)=1
         enddo
         nnod=1
         if(iset.eq.1)nnod=nnset(nbn)

         do i=1,nnod
            nboud2=nboud2+1
            if(nboud2.gt.mbound2)then
               if(iproc_split.eq.0)then
                  write(131,*)'Insufficient Memory-mbound2'
               endif
               stop
            endif
            node=nbn
            if(iset.eq.1) then
            node1=nlist(i,nbn)
            endif
            idf(nboud1)=ifunc
            do ii=1,ndf1        
               if(ii.eq.1)then
!                  nbc2(ndf1*(nboud1-1)+ii)=node1                        <-----------------------------------
                   nbc2(ndf1*(nboud1-1)+ii)=0 
               else
                  if(ndf1*(nboud2-1)+ii.gt.mbound2)then
                     if(iproc_split.eq.0)then
                        write(131,*)'Insufficient Memeory-mbound1'
                     endif
                     stop
                  endif
!                  nbc2(ndf1*(nboud2-1)+ii)=ik(ii-1)
                  nbc2(ndf1*(nboud2-1)+ii)=0
!                  if(ifunc.eq.0)vbc2(ndf*(nboud2-1)+ii-1)=dk(ii-1)
                  if(ifunc.eq.0)vbc2(ndf*(nboud2-1)+ii-1)=0.0d0
               endif
            enddo
         enddo
      enddo
      do n=1,neq
         ibct(n)=0
      enddo
      do n=1,nboud2
         nodes=nbc2(ndf1*(n-1)+1) 
         do i=1,ndf
            index=nbc2(ndf1*(n-1)+1+i)
            if(index.eq.1)then
               nn=ndf*(nodes-1)+i
               ibct(nn)=1
            endif
         enddo         
      enddo
c     
c      read the third type boundary condition (pressure)
c
      call readstr(lr)
      read(lr,*) nboud3

!-----------------------------------
      iface_tet(1)=1
      iface_tet(2)=2
      iface_tet(3)=3

      iface_tet(4)=1
      iface_tet(5)=2
      iface_tet(6)=4

      iface_tet(7)=2
      iface_tet(8)=3
      iface_tet(9)=4
     
      iface_tet(10)=1
      iface_tet(11)=3
      iface_tet(12)=4

!------------------------------------
      
      iface_brick(1)=1
      iface_brick(2)=2
      iface_brick(3)=3
      iface_brick(4)=4

      iface_brick(5)=5
      iface_brick(6)=6
      iface_brick(7)=7
      iface_brick(8)=8

      iface_brick(9)=2
      iface_brick(10)=3
      iface_brick(11)=7
      iface_brick(12)=6

      iface_brick(13)=1
      iface_brick(14)=4
      iface_brick(15)=8
      iface_brick(16)=5

      iface_brick(17)=4
      iface_brick(18)=3
      iface_brick(19)=7
      iface_brick(20)=8

      iface_brick(21)=1
      iface_brick(22)=2
      iface_brick(23)=6
      iface_brick(24)=5


!------------------------------------
      do n=1,nboud3
         read(lr,*)npress,ipf,(ipress_set(i),i=1,npress),
     &        (ipress_sf(i),i=1,npress)
         
         if(ipf.eq.0)then
            read(lr,*)press_val
         endif
         
         do i=1,npress
            iset=ipress_set(i)
            do j=1,nelset(iset)
               nel=nellist(j,iset)
               ibelem(nel)=ipress_sf(i)
               if(ipf.eq.0)then
                  ilf(nel)=0
                  pressb(nel)=press_val
               else
                  ilf(nel)=1
               endif

            enddo
         enddo
      enddo

      call readstr(lr)
C     READ CONSTRAINTS
      read(lr,*)nconst
      nconstr=0

      if(nconst.gt.0)then
         write(*,*)'Constraints not in code'
         stop
      endif
      
      
      do i=1,mbound2
      nbc2(i)=0
      vbc2(i)=0.0d0
      enddo
      
      do i=1,mbound3
      nbc3(i)=0
      vbc3(i)=0.0d0
      enddo
      
!-------------------------------------------------------------      

      return
      end 

************************************************************

      subroutine trainc(nbc2,vbc2,vbcs2,delt,ndf,trac)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      integer::mpiinf(mpiinfsizemax),iprocinfo(maxproc)
      integer::isendorrecv(maxproc)
      common/bounda/nboud1,nboud2,nboud3,nupnch
      common/processorinfo/iproc,nprocs,mpiinfsize,mpiinf,
     &                     iprocinfo,isendorrecv
      dimension vbc2(mbound2),vbcs2(mbound2),nbc2(mbound2),
     &     trac(mdofx)

      ndf1=ndf+1
      do n=1,nboud2
         node=nbc2(ndf1*(n-1)+1)
         do i=1,ndf
            index=nbc2(ndf1*(n-1)+i)
            if(index.eq.1)then
               nn=ndf*(node-1)+i
!               nn=ieqno(nn)
               n1=ndf*(n-1)+i
               trac(nn)=trac(nn)+vbcs2(n1)*delt
            endif
         enddo
      enddo
      return
      end

*******************************************************************

      subroutine loadvc(f,vtau,dvtau0,xund,ijk,imid,npropsg,
     & promat2d,deltme,nbc1,ibc,kinc,ndf,N_update,in_st,irowptr,
     & ibindx,val,resid,isize,iflag,ngauss,pnewdt,elem_grain,xeul,
     & grainsize,Fp_dgamma_node,iproc_split,nprocs_cpfem,comm_split)
      
      implicit none
      include 'pardis.h'
      include 'mpif.h'
      
      integer:: N_update
      real(8)::f(mdofx),ske(mdofe,mdofe),fe(mdofe)
      real(8)::xeun(maxdim,maxeln),xund(maxcrd),ske1(mdofe,mdofe)
      real(8)::f_recv(mdofx),vtau(mdofx),dvtau0(mdofx)
      real(8)::promat2d(maxprops,maxgrp),vtau_e(maxeln*ndofelx)
      real(8)::props(maxprops),time1(2),svars1(maxsv),time(2)
      real(8)::dvtau0_e(maxeln*ndofelx),grainsize(maxel)
      real(8)::extrnval(maxdim*8,mdofproc),extrnrhs(mdofproc,maxproc)
      real(8)::rhsarr(mdofproc),extrnmpi(mdofproc*maxdim*maxeln)
      real(8)::recvval(mdofproc*8*maxdim),recvrhs(mdofproc)
      real(8)::resid(N_update),val(isize),xeul(3,maxel)
      real(8)::deltme,pnewdt,pnewdt_loc,pnewdt_recv
      
      INTEGER::iproc_split,nprocs_cpfem,comm_split         !defined by jiaxi, for coupling
      
      integer::imid(maxdim*maxel),ijk(mnelx),nbc1(mbound1)
      integer::npropsg(maxgrp),ibc(mdofx)
      integer::numrowsproc(maxproc),irowno(mdofproc,maxproc)
      integer::ielno(mdofproc,maxproc),indxcolm(mdofproc,maxproc)
      integer::irowindx(mdofproc,maxproc),ieldof(mdofx*maxeln)
      integer::irowarr(mdofproc),ielarr(mdofproc),irindx(mdofproc)
      integer::kinc,ndf,in_st,isize
      integer::ngauss,iflag,irow,iib,iacolm
      integer::ie,icolm,iel,ierror,iend
      integer::iextrn,iia,ii,iii
      
      integer:: irecvrow(mdofproc),irecvelno(mdofproc)
      integer:: irecvrindx(mdofproc),istatus(MPI_STATUS_SIZE)
      integer:: nelx,neq,nplane,node
      integer:: nx,nsvars1,ndofel1,mdload1
      integer:: npredf1,mcrd1,mpiinf(mpiinfsizemax)
      integer:: iprocinfo(maxproc), isendorrecv(maxproc)
      integer:: nelst,nelend,nelnox
      integer:: ibindx(isize)
      integer:: irowptr(N_update+1)
      integer:: indx1,i,j,irno,isend
      integer:: is,jelem,itemp,jkia,jj
      integer:: jindx,jkia1,k,ll,kk
      integer:: mcrd,mdload,mlvarx,mpiprocno
      integer:: ndofel, ngroup,nno,nn
      integer:: nodeinproc,npredf,iproc
      integer:: mpiinfsize, nel, nrows
      integer:: nprocs,nrecvrows,nprops,nsvars
      integer:: itmp_indx
      integer:: nodeset(4)
      
      real(8):: Fp_dgamma_node(maxnode*12,ngrains)
      real(8):: Fp_dgamma_local(4,12)
      integer:: elem_grain(maxel), grainID

      
       
      common/elemnt/nx,nelx,neq,nplane,node
      common/abq/time1,nsvars1,ndofel1,mdload1,npredf1,svars1,mcrd1
      common/processorinfo/iproc,nprocs,mpiinfsize,mpiinf,
     &                     iprocinfo,isendorrecv
      common/el_par/nelst,nelend,nelnox

      
      ! write(*,*)'irowptr ---1 is', irowptr   
!      write(*,*)'N_update is', N_update 
!      write(*,*)'loadvc check 1 ', iproc


      
      
      do i=1,N_update
         resid(i)=0.0d0
      enddo
      
!      write(*,*)'irowptr ---1-0-0 is', irowptr
      
      itmp_indx=N_update + 1
      
      iend=irowptr(itmp_indx)
      
!       write(*,*)'loadvc check 2 ', iproc   


      time(1)=time1(1)
      time(2)=time1(2)
      nsvars=nsvars1
      ndofel=ndofel1
      mdload=mdload1
      npredf=npredf1
      mcrd=mcrd1
      mlvarx=ndofel

!      write(*,*)'irowptr ---1-0-1 is', irowptr   
****************************linear partitioning of elements***************************

      do i=1,nprocs_cpfem-1
         numrowsproc(i)=0
      enddo

      iextrn=1
***************************************************************************************
      
!      write(*,*)'loadvc check  3 ', iproc   
      
      pnewdt_loc=10.d0
     
      do nel=nelst,nelend


         write(*,*) 'this is elem #',nel
         node=imid(nel)
         ngroup=imid(nelx+nel)
         nprops=npropsg(ngroup)

         
         do ii=1,nprops
            props(ii)=promat2d(ii,ngroup)
         enddo
         
         iib=0
         do iia=1,node
            jkia=ijk(node*(nel-1)+iia) 
            jkia1=mcrd*(jkia-1) 
            do iii=1,mcrd
               xeun(iii,iia)=xund(jkia1+iii)
            enddo
            do ii=1,ndf
               nn=ndf*(jkia-1)+ii
               iib=iib+1
               vtau_e(iib)=vtau(nn)
               dvtau0_e(iib)=dvtau0(nn)
            enddo     
         enddo

         jelem=nel
         ndofel=node*mcrd
         
         
         props(209+15)=grainsize(jelem)                                           !224
         
!-----------------------------------------------------         
         do i=1,4         
         nodeset(i)=ijk((nel-1)*4+i)
         enddo
!----------------------------------------------------         
         
         
       grainID=elem_grain(jelem)
       
          props(209+16)=grainID                                                   ! 225
       
!        if(nel>=1255) then
!        write(*,*)'nel at point 4 is', nel
!        endif      
       do i=1,4
       do j=1,12
      
!       Fp_dgamma_local(i,j)=Fp_dgamma_node((nodeset(i)-1)*12+j,grainID)
       Fp_dgamma_local(i,j)= 0.0d0
       
       enddo
       enddo
!        if(nel>=1255) then
!        write(*,*)'nel at point 5 is', nel
!        endif      
         
!        if(jelem>=1255) then
!        write(*,*)'nel is', nel
!        write(*,*)'loadvc check  5 ', jelem,nodeset(1),nodeset(4) 
!        endif      
         call elst1(jelem,node,xeun,ske1,fe,time,deltme,props,nprops,
     &    nsvars,ndofel,mdload,mlvarx,npredf,mcrd,vtau_e,dvtau0_e,kinc,
     &    iflag,ngauss,pnewdt,nodeset,xeul,Fp_dgamma_local)
     
           
         if(pnewdt.lt.pnewdt_loc)pnewdt_loc=pnewdt


           write(*,*)'check after elst1', iproc , jelem 

c Creating transpose matrix

         if(iflag.eq.0)then
            do ii=1,ndf*node
               do jj=1,ndf*node
                  ske(jj,ii)=ske1(ii,jj)
               enddo
            enddo
         endif
         
         
***************************************************************************************
         do ii=1,node
            nno=ijk((nel-1)*node+ii)
            do jj=1,ndf
               ieldof((ii-1)*ndf+jj)=(nno-1)*ndf+jj
            enddo
         enddo
!          write(*,*)'irowptr ---1-4 is', irowptr   
!     write(*,*)'check in loadvc --1'

         do j=1,node
            nno=ijk((nel-1)*node + j)
            call findproc(nno,nprocs_cpfem,nodeinproc,ndf)

            do k=1,ndf

               irow=(nno-1)*ndf+k
               
!      write(*,*)'irow is', irow
!      write(*,*)'in_st is', in_st   
!      write(*,*)'irowptr ---2  is', irowptr   
               if(nodeinproc.eq.iproc_split)then
                  
                  indx1=irow-in_st
                  is=irowptr(indx1)
                  ie=irowptr(indx1+1)-1

                  if(indx1.eq.N_update)ie=irowptr(N_update+1)



                  if(iflag.eq.0)then

!        write(*,*)'check in loadvc --1-2'
!        write(*,*)'val is',val
!        write(*,*)'is is', is
!        write(*,*)'ie is', ie
!        write(*,*)'ibindx is', ibindx
                     do kk=1,node*ndf
                        iacolm=ieldof(kk)          
                        do ll=is,ie
                           if(ibindx(ll).eq.iacolm)then
                              val(ll)=val(ll)+ske((j-1)*ndf+k,kk)
                              goto 110
                           endif
                        enddo
 110                 enddo
 
                  endif

            ! write(*,*)'check in loadvc --2'
!           write(*,*) 'originial resid is',resid(indx1),indx1
           resid(indx1)=resid(indx1) + fe((j-1)*ndf + k)

!           write(*,*) 'new resid is',resid(indx1),(j-1)*ndf+k,
!     & fe((j-1)*ndf+k)
               else
               
            ! write(*,*)'check in loadvc --2-else'

                  if(nodeinproc.lt.iproc_split)then
                     indx1 = nodeinproc+1
                  else
                     indx1 = nodeinproc
                  endif

                  jindx=numrowsproc(indx1)+1
                  numrowsproc(indx1)=jindx
                  irowno(jindx,indx1)=irow
                  ielno(jindx,indx1)=nel
                  indxcolm(jindx,indx1)=iextrn
                  extrnrhs(jindx,indx1)=fe((j-1)*ndf+k)
                  irowindx(jindx,indx1)=(j-1)*ndf+k

                  if(iflag.eq.0)then

                     do kk=1,node*ndf
                        extrnval(kk,iextrn)=ske((j-1)*ndf+k,kk)
                     enddo
                     
                     iextrn=iextrn+1
                  endif


                  if(jindx.gt.mdofproc.or.iextrn.gt.mdofproc)then
                     write(*,*)'External communication variable 
     &                          mdofproc should be increased'
                     write(*,*)jindx,iextrn,mdofproc
                     stop
                  endif

               endif
        
            enddo
         enddo
         
         
      enddo
!********************************************************************

!        write(*,*)'loadvc check  5 ', iproc  

      call mpi_barrier(comm_split,ierror)


!        write(*,*)'check after barrier ', iproc  
      
      call mpi_allreduce(pnewdt_loc,pnewdt_recv,1,
     & mpi_double_precision,
     & mpi_min,comm_split,ierror)
      
!        write(*,*)'check after allreduce ', iproc 
      
      pnewdt=pnewdt_recv

      if(pnewdt.lt.1.d0)return
      

      do i=1,nprocs_cpfem-1
      
         
   
         mpiprocno=iprocinfo(i)


         isend=isendorrecv(i)

         if(iproc_split.gt.mpiprocno) then 
            itemp=mpiprocno+1
         else 
            itemp=mpiprocno
         endif

         if(isend.eq.1)then

            
         
            nrows=numrowsproc(itemp)
            
!            write(*,*)'loadvc check  6 ', iproc  

            call MPI_SSEND(nrows,1,MPI_INTEGER,mpiprocno,100,
     &                   COMM_SPLIT,ierror)

            do j=1,nrows
               irno=irowno(j,itemp)
               irowarr(j)=irno
               icolm=indxcolm(j,itemp)
               iel=ielno(j,itemp)
               ielarr(j)=iel
               irindx(j)=irowindx(j,itemp)
               rhsarr(j)=extrnrhs(j,itemp)

               if(iflag.eq.0)then
                  do k=1,node*ndf
                     extrnmpi((j-1)*node*ndf+k)=extrnval(k,icolm)
                  enddo
               endif
            enddo

            if(nrows.gt.0)then

!               write(*,*)'send1',iproc,mpiprocno


               call MPI_SSEND(irowarr,nrows,MPI_INTEGER,mpiprocno,100,
     &                   COMM_SPLIT,ierror)

               call MPI_SSEND(ielarr,nrows,MPI_INTEGER,mpiprocno,100,
     &                   COMM_SPLIT,ierror)

               call MPI_SSEND(irindx,nrows,MPI_INTEGER,mpiprocno,100,
     &                   COMM_SPLIT,ierror)


               call MPI_SSEND(rhsarr,nrows,MPI_DOUBLE_PRECISION,
     &                   mpiprocno,100,COMM_SPLIT,ierror)

               if(iflag.eq.0)then
            call MPI_SSEND(extrnmpi,nrows*node*ndf,MPI_DOUBLE_PRECISION,
     &                   mpiprocno,100,COMM_SPLIT,ierror)
               endif
            endif

            call MPI_RECV(nrecvrows,1,MPI_INTEGER,mpiprocno,100,
     &                   COMM_SPLIT,istatus,ierror)

            if(nrecvrows.gt.0)then

!               write(*,*)'recv1',iproc,mpiprocno

               call MPI_RECV(irecvrow,nrecvrows,MPI_INTEGER,mpiprocno,
     &                   100,COMM_SPLIT,istatus,ierror)

               call MPI_RECV(irecvelno,nrecvrows,MPI_INTEGER,mpiprocno,
     &                  100,COMM_SPLIT,istatus,ierror)

              call MPI_RECV(irecvrindx,nrecvrows,MPI_INTEGER,mpiprocno,
     &                  100,COMM_SPLIT,istatus,ierror)

         call MPI_RECV(recvrhs,nrecvrows,MPI_DOUBLE_PRECISION,mpiprocno,
     &                  100,COMM_SPLIT,istatus,ierror)
               
               if(iflag.eq.0)then
          call MPI_RECV(recvval,nrecvrows*ndf*node,MPI_DOUBLE_PRECISION,
     &                   mpiprocno,100,COMM_SPLIT,istatus,ierror)
               endif

         call update_local_array(ijk,nrecvrows,recvval,recvrhs,irecvrow,
     & irecvelno,irecvrindx,ndf,iproc,N_update,in_st,irowptr,ibindx,val,
     &  resid,isize,iflag)
     
     

     

            endif
            
!            write(*,*)'loadvc check  7 ', iproc  

         else

            call MPI_RECV(nrecvrows,1,MPI_INTEGER,mpiprocno,100,
     &                   COMM_SPLIT,istatus,ierror)


            if(nrecvrows.gt.0)then

!               write(*,*)'recv2',iproc,mpiprocno

               call MPI_RECV(irecvrow,nrecvrows,MPI_INTEGER,mpiprocno,
     &                   100,COMM_SPLIT,istatus,ierror)

               call MPI_RECV(irecvelno,nrecvrows,MPI_INTEGER,mpiprocno,
     &                  100,COMM_SPLIT,istatus,ierror)

              call MPI_RECV(irecvrindx,nrecvrows,MPI_INTEGER,mpiprocno,
     &                  100,COMM_SPLIT,istatus,ierror)


         call MPI_RECV(recvrhs,nrecvrows,MPI_DOUBLE_PRECISION,mpiprocno,
     &                  100,COMM_SPLIT,istatus,ierror)

               if(iflag.eq.0)then
         call MPI_RECV(recvval,nrecvrows*node*ndf,MPI_DOUBLE_PRECISION,
     &                   mpiprocno,100,COMM_SPLIT,istatus,ierror)
               endif

         call update_local_array(ijk,nrecvrows,recvval,recvrhs,irecvrow,
     & irecvelno,irecvrindx,ndf,iproc,N_update,in_st,irowptr,ibindx,val,
     & resid,isize,iflag)
     
     

            endif
            
            nrows=numrowsproc(itemp)

            call MPI_SSEND(nrows,1,MPI_INTEGER,mpiprocno,100,
     &                   COMM_SPLIT,ierror)
         
         
                 ! write(*,*)'check in loadvc -2'
         
            do j=1,nrows
               irno=irowno(j,itemp)
               irowarr(j)=irno
               icolm=indxcolm(j,itemp)
               iel=ielno(j,itemp)
               ielarr(j)=iel
               rhsarr(j)=extrnrhs(j,itemp)
               irindx(j)=irowindx(j,itemp)
               if(iflag.eq.0)then
                  do k=1,node*ndf
                     extrnmpi((j-1)*node*ndf+k)=extrnval(k,icolm)
                  enddo
               endif
            enddo


            if(nrows.gt.0)then

!               write(*,*)'send2',iproc,mpiprocno

               call MPI_SSEND(irowarr,nrows,MPI_INTEGER,mpiprocno,100,
     &                   COMM_SPLIT,ierror)

               call MPI_SSEND(ielarr,nrows,MPI_INTEGER,mpiprocno,100,
     &                   COMM_SPLIT,ierror)

               call MPI_SSEND(irindx,nrows,MPI_INTEGER,mpiprocno,100,
     &                   COMM_SPLIT,ierror)

         call MPI_SSEND(rhsarr,nrows,MPI_DOUBLE_PRECISION,mpiprocno,100,
     &                   COMM_SPLIT,ierror)

               if(iflag.eq.0)then
            call MPI_SSEND(extrnmpi,nrows*node*ndf,MPI_DOUBLE_PRECISION,
     &                   mpiprocno,100,COMM_SPLIT,ierror)
               endif

            endif

         endif


      enddo

!      write(*,*)'loadvc check  8 ', iproc  

      call remove_dispdof(ibc,N_update,in_st,irowptr,ibindx,val,
     &                     resid,isize,iflag)


      do i=1,N_update
         f(i+in_st)=resid(i)
      enddo

               ! write(*,*)'check in loadvc -3'

      call mpi_barrier(COMM_SPLIT,ierror)
      call mpi_allreduce(f,f_recv,neq,
     & mpi_double_precision,
     & mpi_sum,COMM_SPLIT,ierror)


      do i=1,neq
         f(i)=f_recv(i)
!         if(f(i)>1e4) then
!            write(*,*) 'f(i) ',i, f(i)
!          endif
      enddo
      
!             write(*,*)'loadvc check  9 ', iproc  
             
      return
      end 

************************************************
      subroutine exchn(a,b,n)   
      implicit real*8(a-h,o-z)
      dimension a(1),b(1)
      do 10 i=1,n
10    b(i)=a(i)
      return
      end

**************************************************
      
      subroutine readnset(iproc_split)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      integer:: nnset(maxnset),nlist(maxns,maxnset)
      integer::mpiinf(mpiinfsizemax),iprocinfo(maxproc)
      integer::isendorrecv(maxproc)
      integer:: n
      common/cbnset/nset,nnset,nlist
      common/processorinfo/iproc,nprocs,mpiinfsize,mpiinf,
     &                     iprocinfo,isendorrecv

      call readstr(lr)
      read(lr,*)nset
      
      if(nset.gt.maxnset)then
         if(iproc_split.eq.0)then
            write(131,*)'Insufficient Memory-maxnset'
         endif
         stop
      endif
     
      do i=1,nset
         read(lr,*)nsno,ngen
         if(ngen.eq.0)then
            read(lr,*)numn,(nlist(j,nsno),j=1,numn)
         endif
         if(ngen.eq.1)then
            read(lr,*)istart,iend,inc
            numn=0
            do n=istart,iend,inc
               numn=numn+1
               nlist(numn,nsno)=n
            enddo
         endif
         nnset(nsno)=numn
      enddo
      return
      end

*******************************************************
      subroutine readelset()
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension nelset(maxelset),nellist(maxels,maxelset)
      common/cbelset/neset,nelset,nellist
      call readstr(lr)
      read(lr,*)neset
      
      if(neset.gt.maxelset)then
         write(131,*)'Insufficient Memory-maxelset'
         stop
      endif
     
      do i=1,neset
         read(lr,*)nsno,ngen
         if(ngen.eq.0)then
            read(lr,*)numn,(nellist(j,nsno),j=1,numn)
         endif
         if(ngen.eq.1)then
            read(lr,*)istart,iend,inc
            numn=0
            do n=istart,iend,inc
               numn=numn+1
               nellist(numn,nsno)=n
            enddo
         endif
         nelset(nsno)=numn
      enddo
      return
      end

!**********************************************************
      function dotprd(v1,v2,nbc1,ieqno)
      implicit real*8(a-h,o-z)
      include 'pardis.h'     
      dimension v1(1),v2(1),nbc1(1),nfac(1),ieqno(1)
      common/bounda/nboud1,nboud2,nboud3,nupnch
       common/dof/ndof

      ndf=ndof
c
      dotprd=fdot(v1,v2)
      ndf1=ndf+1
      if(nboud1.eq.0)go to 15
      do 10 n=1,nboud1
      node=nbc1(ndf1*(n-1)+1)
      do 12 i=1,ndf
      index=nbc1(ndf1*(n-1)+1+i)
      if(index.ne.1)go to 12
      nn=ndf*(node-1)+i
      nn=ieqno(nn)
      dotprd=dotprd-v1(nn)*v2(nn)
12    continue
10    continue   
15    continue
      return
      end

***************************************************************
 
      subroutine elst1(jelem,nnode,xeun,ske,fe,time,deltme,props1,nprops
     &    ,nsvars,ndofel,mdload,mlvarx,npredf,mcrd,vtau_e,dvtau0_e,kinc,
     &    iflag,ngauss,pnewdt,nodeset,xeul,Fp_dgamma_local,flag_SC,
     &    elem_strain,elem_stress,elem_V,elem_energy)

!      implicit real *8(a-h,o-z)
      implicit none
      include 'pardis.h'

      
     
      
      integer::  ntens, nrhs, jelem, nnode
      integer::  nsvars, nprops
      integer::  ndofel, mdload, ngauss
      integer::  mlvarx,npredf,mcrd,kinc,iflag
      integer::  nodeset(4) 
      parameter (ntens=6)
      real(8):: deltme,pnewdt
      real(8)::  fe(mdofe),ske(mdofe,mdofe)
      real(8)::  xeun(maxdim,maxeln)
      real(8)::  rhs(mlvarx,1),svars(msv),energy(maxnp)
      real(8)::  u(ndofel),du(mlvarx,1),v(ndofel),a(ndofel)
      real(8)::  time(2),params(1),adlmag(mdload,1)
      real(8)::  ddlmag(mdload,1),predef(2,npredf,nnode)
      real(8)::  props(maxprops),props1(maxprops)
      real(8)::  svars1(maxsv),coords(mcrd,nnode),time1(2)
      real(8)::  svars2(maxsv),amatrx(ndofel,ndofel)
      real(8)::  vtau_e(maxeln*ndofelx),dvtau0_e(maxeln*ndofelx)
      real(8):: xeul(3,maxel), euler(3)      
      integer:: jprops(1),jdltyp(mdload,1),lflags(10)
      integer:: jtype,ndload
      
!-----------------------------------------------------------------------      
      real(8):: period
      integer:: njprop
!-----------------------------------------------------------------------
      
      integer::  i,j,ndofel1,nx
      integer:: igauss_init,iconv1,ii
      integer::  iinp,inumj,iwr
      integer:: nelend,nelnox,neq
      integer::  nelx,node, nplane,npredf1  
      integer::  nelst,k,kstep,mcrd1,mdload1
      integer::  nstatv,nstr,nsv,nsvars1 
      integer::  ibelem(maxel),ilf(maxel),iface_tet(12)     
      integer::  iface_brick(24),ibnode(4)
      real(8)::  statev(maxstatev),pressb(maxel), dtime


! ------------ SC elem stress strain energy volume ------
      integer:: flag_SC
      real(8):: elem_stress(6),elem_strain(6),elem_energy,elem_V
! --------------------------------------------------------    


!------------GND-----------------------------------------      
      integer:: choose
      real(8):: Fp_dgamma_local(4,12),B_GND(3,4)
      real(8):: G_GND(3,4),tmp_GND(4),Cmpen_factor1
      real(8):: x_star(3,4), Cmpen_factor2,Cmpen_factor3
!------------ Umat ------------------------------------------
        
      real(8):: dtemp, celent, drpldt,sse
      real(8):: layer,  scd, rpl, spd       
      integer:: ndi, nshr, kspt    


      character*80 cmname 
!----------------------------------------------------------------------
     

      common/press/pressb,ibelem,ilf,iface_tet,iface_brick      
      common/abq/time1,nsvars1,ndofel1,mdload1,npredf1,svars1,mcrd1     
      common/stateupdate/svars2
      common/elemnt/nx,nelx,neq,nplane,node
      common/el_par/nelst,nelend,nelnox
      common/inits/igauss_init
      common/conv/iconv1
      common/numj/inumj
      common/wr/iwr

      data igauss_init/0/
      data iinp/0/

      

      
      
      euler(1)=xeul(1,jelem)
      euler(2)=xeul(2,jelem)
      euler(3)=xeul(3,jelem)
       

       
       
       


      do i=1,10
         lflags(i)=0
      enddo

      do i=1,nprops
         props(i)=props1(i)
      enddo

      nsv=nsvars*ngauss*(jelem-1)-nsvars*ngauss*(nelst-1)
      do i=1,nsvars*ngauss
         svars(i)=svars1(nsv+i)
      enddo

      do i=1,mcrd
         do j=1,nnode
            coords(i,j)=xeun(i,j)
         enddo
      enddo
      
      dtime=deltme
      lflags(1)=1

********************************
      do i=1,ndofel
         du(i,1)=0.0
      enddo
*******************************
      
      do i=1,ndofel
         u(i)=vtau_e(i)
         du(i,1)=dvtau0_e(i)
      enddo 

      kstep=1
      pnewdt=1.0 
     
!****************************************************************************      
      call uel_disp_t(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     & props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     & kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     & lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period,iflag,
     & ngauss,Fp_dgamma_local,euler,flag_SC,
     & elem_strain,elem_stress,elem_V,elem_energy)
!****************************************************************************

      if(pnewdt.lt.1.d0)then
         write(131,*)'User Routine Requests smaller increment'
      
         return
      endif

      do i=1,ndofel
         fe(i)=rhs(i,1)
         do j=1,ndofel
            ske(i,j)=amatrx(i,j)
         enddo
      enddo
      
      do i=1,nsvars*ngauss
         svars2(nsv+i)=svars(i)
      enddo

!      fpdot_el(jelem-nelst+1,2,1:3,1:3,1:2)=fpdot(1:3,1:3,1:2)

      return
      end
           
****************************************************
      function fdot(a,b)
      implicit real*8(a-h,o-z)
      dimension a(1),b(1)
      common/elemnt/nx,nelx,neq,nplane,node
      fdot=0.d0
      do 100 i=1,neq
  100 fdot=fdot+a(i)*b(i)
      return
      end
****************************************************
      subroutine get_subscript(iproc,fname)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      character*10 fname
      dimension in(10),intmp(10)                                  !  <---------------------
      
      iproc1=iproc
      i=0

      do while(iproc1.gt.9)
         i1=iproc1/10
         i2=iproc1-i1*10
         i=i+1
         iproc1=i1
         in(i)=i2
      enddo

      i=i+1
      in(i)=iproc1

      do j=1,i
         intmp(i-j+1)=in(j)
      enddo

      do j=1,i
         in(j)=intmp(j)
      enddo

      if(i.eq.1)then
         fname='_'//achar(48+in(1))//'.out'
      elseif(i.eq.2)then
         fname='_'//achar(48+in(1))//achar(48+in(2))//'.out'
      else
         fname='_'//achar(48+in(1))//achar(48+in(2))//achar(48+in(3))
     &   //'.out'
      endif
      
      
      return
      end
!---------------------------------------------------------------
! Returns time based on dwell load such that a time point is present 
! at the change points
      subroutine tdwell(fintme,deltme,xidt,xmindt,iwrite_flag,nseg,tseg)
      implicit real*8(a-h,o-z)

      dimension tseg(nseg+1)

      common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin,uamp,
     &       umin
      common/load_mono/sload,st_rate
      common/icyc_flag/icyc

       
      xtol=1d-6

      
      if(tramp-fintme.gt.xtol)then

         fintme1=fintme+deltme+xmindt
         if(fintme1.ge.tramp)then
            deltme=tramp-fintme
            iwrite_flag=1
         endif



      else

         fintme2=fintme-tramp

         rem=dmod(fintme2,tperiod)

         if(dabs(rem).lt.xtol.or.dabs(rem-tperiod).lt.xtol)then
            xncyc=fintme/tperiod
         else
            xncyc=(fintme2-rem)/tperiod
         endif
 
         t=fintme2-xncyc*tperiod


         do i=1,nseg

            tst=tseg(i)
            tend=tseg(i+1)

            teq=t-tst


         if((teq.gt.xtol.or.dabs(teq).lt.xtol).and.(tend-t).gt.xtol)then


               t1=t+deltme+xmindt
               t1eq=t1-tend

               if(t1eq.gt.xtol.or.dabs(t1eq).lt.xtol)deltme=tend-t
               if(dabs(teq).lt.xtol)deltme=xidt


            endif


         enddo

      endif
      
      return
      end

!----------------------------------------------------------------------------------
      subroutine broy_parallel(dvtau,var,war,numupd,nmpd1,maxref,in_st,
     &             N_update,neq,step,COMM_SPLIT)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      include 'mpif.h'
      dimension dvtau(mdofx),var(*),war(*)
      allocatable::v(:),w(:),v_recv(:),w_recv(:)
 
      allocate(v(neq),w(neq),v_recv(neq),w_recv(neq))
      v(1:neq)=0.d0
      w(1:neq)=0.d0
      v_recv(1:neq)=0.d0
      w_recv(1:neq)=0.d0

      step1=1.d0/step


      if(numupd.eq.0)goto 30

      do i=1,numupd
         nn=(i-1)*N_update

         v(1:neq)=0.d0
         w(1:neq)=0.d0

         v(in_st+1:in_st+N_update)=var(nn+1:nn+N_update)
         w(in_st+1:in_st+N_update)=war(nn+1:nn+N_update)
      
         call mpi_barrier(COMM_SPLIT,ierror)
         call mpi_allreduce(v,v_recv,neq,
     &    mpi_double_precision,
     &    mpi_sum,COMM_SPLIT,ierror)
         call mpi_allreduce(w,w_recv,neq,
     &    mpi_double_precision,
     &    mpi_sum,COMM_SPLIT,ierror)

         v(1:neq)=v_recv(1:neq)
         w(1:neq)=w_recv(1:neq)

         dw=dot_product(dvtau(1:neq),w)

         dvtau(1:neq)=dvtau(1:neq)+dw*v(1:neq)
      enddo
 30   continue

      numupd=numupd+1
      
      w(1:neq)=0.d0
      w(in_st+1:in_st+N_update)=war(nmpd1+1:nmpd1+N_update)

      call mpi_barrier(COMM_SPLIT,ierror)
      call mpi_allreduce(w,w_recv,neq,
     & mpi_double_precision,
     & mpi_sum,COMM_SPLIT,ierror)

      w(1:neq)=w_recv(1:neq)

      v(1:neq)=step1*w(1:neq)-dvtau(1:neq)

      fac1=dot_product(v,w)

      fac=1.d0/fac1


      v(1:neq)=(w(1:neq)-v(1:neq))*fac

      var(nmpd1+1:nmpd1+N_update)=v(in_st+1:in_st+N_update)

      dw=dot_product(dvtau(1:neq),w)

      dvtau(1:neq)=dvtau(1:neq)+dw*v(1:neq)

      deallocate(v,w,v_recv,w_recv)

      return
      end
!--------------------------------------------------------------------

c$$$      subroutine tpredict(deltme_prev,deltme_pred,nstep,ipoly_write,xmaxdt)
c$$$      implicit real*8(a-h,o-z)
c$$$      include 'pardis.h'
c$$$      common/tpred_crit/fpdot_el(maxsv,2,3,3,2)
c$$$      common/el_par/nelst,nelend,nelnox
c$$$      dimension time1(2),svars1(maxsv)
c$$$      common/abq/time1,nsvars1,ndofel1,mdload1,npredf1,svars1,mcrd1
c$$$      dimension ipoly_write(maxel)
c$$$      dimension fp_n(3,3),d2fpdt2_n(3,3)
c$$$      
c$$$
c$$$      deltme_min=xmaxdt
c$$$      eta=1d-3
c$$$
c$$$      if(nstep.eq.1)then
c$$$         
c$$$         do nel=nelst,nelend
c$$$            iph=ipoly_write(nel)
c$$$
c$$$            do itag=1,ph
c$$$               fpdot_el(nel-nelst+1,1,1:3,1:3,itag)=fpdot_el(nel-nelst+1,2,1:3,1:3,itag)
c$$$            enddo
c$$$            
c$$$         enddo
c$$$
c$$$         deltme_pred=1000.d0
c$$$
c$$$      else
c$$$
c$$$         do nel=nelst,nelend
c$$$
c$$$            iph=ipoly_write(nel)
c$$$
c$$$            do itag=1,iph
c$$$
c$$$               if(itag.eq.1)then
c$$$
c$$$                  ist=(nel-nelst)*maxstatev+66
c$$$                  
c$$$                  do i=1,3
c$$$                     do j=1,3
c$$$                        fp_n(i,j)=svars1(ist+(i-1)*3+j)
c$$$                     enddo
c$$$                  enddo
c$$$
c$$$               else
c$$$
c$$$                  ist=(nel-nelst)*maxstatev+177
c$$$                                    
c$$$                  do i=1,3
c$$$                     do j=1,3
c$$$                        fp_n(i,j)=svars1(ist+(i-1)*3+j)
c$$$                     enddo
c$$$                  enddo
c$$$
c$$$               endif
c$$$
c$$$
c$$$               d2fpdt2_n(1:3,1:3)=(fpdot_el(nel-nelst+1,2,1:3,1:3,itag)-fpdot_el(nel-nelst+1,1,1:3,1:3,itag))/deltme_prev
c$$$
c$$$
c$$$               fp_n(1,1)=fp_n(1,1)-1.d0
c$$$               fp_n(2,2)=fp_n(2,2)-1.d0
c$$$               fp_n(3,3)=fp_n(3,3)-1.d0
c$$$
c$$$               call frob_norm(fp_n,3,fpnorm)
c$$$               call frob_norm(d2fpdt2_n,3,d2fpnorm)
c$$$
c$$$               fpdot_el(nel-nelst+1,1,1:3,1:3,itag)=fpdot_el(nel-nelst+1,2,1:3,1:3,itag)
c$$$
c$$$               if(fpnorm.lt.1d-8.or.d2fpnorm.lt.1d-10)then
c$$$                  deltme=xmaxdt
c$$$               else
c$$$                  deltme=(2.d0*eta*fpnorm/d2fpnorm)**0.5d0
c$$$
c$$$                  if(deltme.gt.xmaxdt)deltme=xmaxdt
c$$$               endif
c$$$
c$$$               if(deltme.lt.deltme_min)deltme_min=deltme
c$$$
c$$$            enddo
c$$$            
c$$$         enddo
c$$$         
c$$$         
c$$$         
c$$$
c$$$      endif
c$$$
c$$$      deltme_pred=deltme_min
c$$$
c$$$      end


!-----------------------------

      subroutine frob_norm(x,n,xnorm)
      implicit real*8(a-h,o-z)
      dimension x(n,n)
      
      xnorm=0.d0
      do i=1,n
         do j=1,n
            xnorm=xnorm+x(i,j)**2.d0
         enddo
      enddo

      xnorm=dsqrt(xnorm)
      return
      end 


! Returns time based for creep load or constant strain rate problem
 
      subroutine t_creep_cstrainrate(fintme,deltme,xidt,xmindt,xmaxdt)
      implicit real*8(a-h,o-z)

      common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin,uamp,
     &       umin
      common/load_mono/sload,st_rate
      common/icyc_flag/icyc

      if(deltme.gt.xmaxdt)then
         deltme=xmaxdt
      endif

       
      xtol=1d-6

      
      if(tramp-fintme.gt.xtol)then

         fintme1=fintme+deltme+xmindt
         if(fintme1.ge.tramp)then
            deltme=tramp-fintme
            iwrite_flag=1
         endif

      endif

      return      
      end


!----------------------------------------------------------------------------------

!*******************************************************************
!     Creates a sequence of processor pairs communicating(send/recv)
!     Objective: Reduces waiting time for communicating processor
!*******************************************************************
      subroutine proc_comb(nproc,mpiinf,mpiinfsize)
      implicit real*8 (a-h,o-z)
      include 'pardis.h'

      dimension mpiinf(mpiinfsizemax)
      dimension iparr1(maxcomm,2),iparr(2)

      allocatable::igrparr(:,:),igrpndx(:)
      
!*******************************************************************
!     nproc(I)-number of processors
!     mpiinf(O)-processor numbers in pairs. Ex - (0123..), means 
!     01,23,.. are different communicating pairs
!     mpiinfsize(O)-length of array mpiinf
!*******************************************************************

      nsize=nproc*(nproc-1)

      if(nsize.gt.maxcomm)then
        write(*,*)'Increse maxcomm to',nsize
        stop
      endif

      
      ilevel=0
      i2=nproc
      do while(i2.gt.1)
         ilevel=ilevel+1
         i1=i2/2
         i2=i2-i1
      enddo


      allocate(igrparr(maxcomm,ilevel),igrpndx(ilevel))
         
      do i=1,ilevel
            
         if(i.eq.1)then
           i1=nproc/2
           i2=nproc-i1
           igrpndx(i)=2
           igrparr(1,i)=i2
           igrparr(2,i)=i1
               
         else
           iloop=igrpndx(i-1)
           icnt=0
               
           do j=1,iloop
              ip1=igrparr(j,i-1)
              if(ip1.gt.1)then
                 i1=ip1/2
                 i2=ip1-i1
                 icnt=icnt+1
                 igrparr(icnt,i)=i2
                 icnt=icnt+1
                 igrparr(icnt,i)=i1
              else
                 icnt=icnt+1
                 igrparr(icnt,i)=0
              endif
           enddo

           igrpndx(i)=icnt

         endif

      enddo

      igndx=0
      do i=1,ilevel
         itemp=1
         iprocbegin=0

         do j=1,igrpndx(i)

            if(itemp.le.2)then
               ip1=igrparr(j,i)
               if(ip1.eq.0) then
                  iprocbegin=iprocbegin+1
                  goto 130
               else
                  iparr(itemp)=igrparr(j,i)
                  do k=0,iparr(itemp)-1
                     iparr1(k+1,itemp)=iprocbegin
                     iprocbegin=iprocbegin+1
                  enddo
                  itemp=itemp+1
               endif
 130        endif
   
            if(itemp.gt.2)then
               itemp=1
               indx1=1
               indx2=2
               if(iparr(1).gt.iparr(2)) then
                  indx1=2
                  indx2=1
               endif

               do k=0,iparr(indx2)-1
                  do l=1,iparr(indx1)
                     igndx=igndx+1
                     mpiinf(igndx)=iparr1(l,indx1)
                     igndx=igndx+1
                     i1=iparr(indx2)
                     i2=l+k
                     i3=(i2-1)/i1
                     iloc=i2-i3*i1
                     mpiinf(igndx)=iparr1(iloc,indx2)
                  enddo
               enddo
                  
            endif

         enddo
      enddo

      mpiinfsize=igndx

      if(mpiinfsize.gt.mpiinfsizemax)then
         write(*,*)'mpiinfsize',mpiinfsize,mpiinfsizemax
         stop
      endif

      deallocate(igrparr,igrpndx)

      end
!***************************END*****************************************
!***********************************************************************
!     Each processor identifies its sequence of communication with
!     other processors form mpiinf array
!***********************************************************************
      subroutine arrange_proc(iprocinfo,isendorrecv,mpiinf,           
     &    mpiinfsize,iproc_split)

      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension  mpiinf(mpiinfsizemax),iprocinfo(maxproc)
      dimension isendorrecv(maxproc)
!***********************************************************************
!     iprocinfo(O)-array containing processor numbers arranged in best
!     communication sequence
!     isendorrecv(O)-array containg tag which decides whether proc a would
!     send and then receive from b, or vice versa
!     Tag 1 - Send First
!     Tag 2 - Receive First
!     mpiinf(I)
!      mpiinfsize(I)
!     iproc(I)
!***********************************************************************
      indx=0

      do i=1,mpiinfsize,2
         ip1=mpiinf(i) 
         ip2=mpiinf(i+1)
         if(ip1.eq.iproc_split)then 
            indx=indx+1
            iprocinfo(indx)=ip2 
            isendorrecv(indx)=1
         elseif(ip2.eq.iproc_split)then
            indx=indx+1
            iprocinfo(indx)=ip1
            isendorrecv(indx)=2
         endif
      enddo

      end
!***********************END*********************************************

!***********************************************************************
!c     Creates linear array of element numbers for each node
!***********************************************************************      
      subroutine node_elem_connec(ijk,ndf)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      common/elemnt/nx,nelx,neq,nplane,node
      common/nodeconnec/ielemno,nindx
      dimension ijk(1),ielemno(maxnodelem*maxnode),nindx(maxnode)
!***********************************************************************
!c     ielemno: stores element numbers
!c     nindx: stores start and end index in ielemno array for
!c            each node
!***********************************************************************

      
      itemp=1
      do i=1,nx
         do j=1,nelx
            do k=1,node
               if(ijk((j-1)*node+k).eq.i)then
                  ielemno(itemp)=j
                  itemp=itemp+1
                  goto 20
               endif
            enddo
   20    enddo
         nindx(i)=itemp
      enddo 
      end

!*************************END*******************************************

!*************************END*******************************************
!***********************************************************************
!c     Finds proc no associated with a node considering linear
!c     partitioning of grid points
!***********************************************************************
      subroutine findproc(nno,nprocs_cpfem,nodeinproc,ndf)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      common/elemnt/nx,nelx,neq,nplane,node
!***********************************************************************
!c     nno(I)-node number
!c     nprocs-number of processors
!c     nodeinproc(O)-processor number to which node is associated
!***********************************************************************

      i1=nx/nprocs_cpfem+1
      i2=nx/nprocs_cpfem

      ia=nx-i2*nprocs_cpfem
      ib=nprocs_cpfem-ia


      if(nno.le.i1*ia)then
         nodeinproc=(nno-1)/i1
      else
         nno1=nno-i1*ia
         nodeinproc=ia+(nno1-1)/i2
      endif


      end
!*************************END*******************************************
      subroutine update_local_array(ijk, nrows,recvval,recvrhs,irowarr, 
     &   ielarr,irecvrindx,ndf,iproc,N_update,in_st,irowptr,ibindx,val, 
     &   resid,isize,iflag)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension ijk(mnelx)
      dimension irowarr(mdofproc),ielarr(mdofproc),                
     &          recvval(mdofproc*maxdim*8),recvrhs(mdofproc),        
     &          irecvrindx(mdofproc)
      dimension icolmarr(mdofe)
      common/elemnt/nx,nelx,neq,nplane,node
      dimension ibindx(isize),val(isize),resid(N_update),          
     & irowptr(N_update+1)



      do i=1,nrows
         irowno=irowarr(i)
         ielno=ielarr(i)


         do j=1,node
            nno=ijk((ielno-1)*node+j)
            do k=1,ndf
               itemp=(nno-1)*ndf+k
               icolmarr((j-1)*ndf+k)=itemp
            enddo
         enddo

         irow1=irowno-in_st
         is=irowptr(irow1)
         ie=irowptr(irow1+1)-1
               
         if(irow1.eq.N_update)ie=irowptr(N_update+1)

         if(iflag.eq.0)then
               
            do k=is,ie
               icolm=ibindx(k)
               do kk=1,node*ndf
                  if(icolm.eq.icolmarr(kk))then
                     val(k)=val(k)+recvval((i-1)*ndf*node+kk)
                     goto 120
                  endif
               enddo
 120        enddo
         endif


         resid(irow1)=resid(irow1)+recvrhs(i)

      enddo
         


      end

!*************************END*******************************************
!***********************************************************************


!*************************END*******************************************
!******Performs division of dof to each processor
      subroutine node_partition(N_update,in_st,iproc_split,nprocs_cpfem,
     & nx,ndf)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      



      i1=nx/nprocs_cpfem+1
      i2=nx/nprocs_cpfem
      
      ia=nx-i2*nprocs_cpfem
      ib=nprocs_cpfem-ia

      i1=i1*ndf
      i2=i2*ndf

      if((iproc_split+1).le.ia)then
         in_st=iproc_split*i1
         N_update=i1

      else
         in_st=ia*i1+(iproc_split-ia)*i2
         N_update=i2
      endif

      return
      end
!***********************end*********************************************
















      subroutine makegrad_tet4(xb_st,xbtr_st,xg_st,shptau_st,dvoltau_st,
     &             nnode,mcrd,ndofel,ngauss)
      implicit real*8(a-h,o-z)
      dimension xb_st(2*3,ndofel,ngauss),dvoltau_st(ngauss),               
     &      shptau_st(mcrd,nnode,ngauss),                               
     &     nshr(3),xb_av(ndofel),xbtr_st(ndofel,ngauss),               
     &     xg_st(9,ndofel,ngauss),nshr0(3)


      nshr(1)=2
      nshr(2)=3
      nshr(3)=3

      nshr0(1)=1
      nshr0(2)=1
      nshr0(3)=2
      
      if (mcrd.eq.3)then
         nstr=6
      endif

      if (mcrd.eq.2)then
         nstr=3
      endif
      xb_av(1:ndofel)=0.d0

      do ipt=1,ngauss
         do i=1,nstr
            do j=1,ndofel
               xb_st(i,j,ipt)=0.d0
               jtemp1=(j-1)/mcrd
               jtemp2=j-jtemp1*mcrd
               if(i.eq.jtemp2.and.i.le.mcrd) then
                  xb_st(i,j,ipt)=shptau_st(i,jtemp1+1,ipt)
               endif
               if(i.gt.mcrd)then
                  is=i-mcrd
                  ishr=nshr(is) 
                  is=nshr0(is)
                  if(jtemp2.eq.is) then
                     xb_st(i,j,ipt)=shptau_st(ishr,jtemp1+1,ipt)
                  endif
                  if(jtemp2.eq.ishr) then
                     xb_st(i,j,ipt)=shptau_st(is,jtemp1+1,ipt)
                  endif
               endif
            enddo
         enddo
         
  
         do j=1,ndofel
            xbtr_st(j,ipt)=0.d0
            do i=1,mcrd
              xbtr_st(j,ipt)=xbtr_st(j,ipt)+xb_st(i,j,ipt)
            enddo
            do i=1,mcrd
               xb_st(i,j,ipt)=xb_st(i,j,ipt)-1.d0/3.d0*xbtr_st(j,ipt)
            enddo
            xb_av(j)=xb_av(j)+xbtr_st(j,ipt)*dvoltau_st(ipt)
         enddo
         do j=1,ndofel
            jtemp=(j-1)/mcrd
            jtemp1=jtemp+1
            jtemp2=j-jtemp*mcrd
            do i=1,mcrd**2
               itemp=(i-1)/mcrd
               itemp1=itemp+1
               itemp2=i-itemp*mcrd
               xg_st(i,j,ipt)=0.d0
               if(jtemp2.eq.itemp1) then
                  xg_st(i,j,ipt)=shptau_st(itemp2,jtemp1,ipt)
               endif
               if(itemp2.eq.1) then
                  xg_st(i,j,ipt)=xg_st(i,j,ipt)-1.d0/3.d0*xbtr_st(j,ipt)
               endif
            enddo
         enddo    
      enddo


      tv=sum(dvoltau_st(1:ngauss))
      xb_av(1:ndofel)=xb_av(1:ndofel)/tv
      do ipt=1,ngauss
         do j=1,ndofel
            do i=1,mcrd
               xb_st(i,j,ipt)=xb_st(i,j,ipt)+1.d0/3.d0*xb_av(j)
               itemp=(i-1)*mcrd+1
               xg_st(itemp,j,ipt)=xg_st(itemp,j,ipt)+1.d0/3.d0*xb_av(j)
            enddo
         enddo
      enddo
      return
      end


      subroutine calc_dfgrd_tet4(ipt,u,mcrd,nnode,shp_st,dfgrd)
      implicit real*8 (a-h,o-z)
      dimension u(mcrd,nnode),dfgrd(3,3),shp_st(3,nnode,nnode)

      do i=1,mcrd
         do j=1,mcrd
            dfgrd(i,j)=0.d0
            if(i.eq.j)dfgrd(i,j)=1.d0
            do in=1,nnode
               dfgrd(i,j)=dfgrd(i,j)+shp_st(j,in,ipt)*u(i,in)
            enddo
         enddo
      enddo
      if(mcrd.eq.2)then
         dfgrd(1:3,3)=0.d0
         dfgrd(3,1:3)=0.d0
         dfgrd(3,3)=1.d0
      endif
      return
      end




      subroutine shape_tet4(xel,ipt,shp,dvol,xc)
      implicit none
      
      real(8):: shp1(3,4),shp(3,4),xjac(3,3),xjacinv(3,3)
      real(8):: xel(3,4),xc(3)
      real(8):: x(4),y(4),z(4),xx(3),yy(3),zz(3)
      real(8):: dvol
      integer:: i,j, inode, ipt

      do i=1,3
      do j=1,4
      shp1(i,j)=0.d0
      shp(i,j)=0.d0
      enddo
      enddo

      shp1(1,1)=1.d0
      shp1(2,2)=1.d0
      
      shp1(3,3)=1.d0
      shp1(1,4)=-1.d0
      shp1(2,4)=-1.d0

      shp1(3,4)=-1.d0

      
!		write(*,*) 'xel   is:  ', xel

      do i=1,4
         x(i)=xel(1,i)
         y(i)=xel(2,i)
         z(i)=xel(3,i)
      enddo

      do i=1,3
         xx(i)=xel(1,i)-xel(1,4)
         yy(i)=xel(2,i)-xel(2,4)
         zz(i)=xel(3,i)-xel(3,4)
      enddo

      xjac(1:3,1)=xx(1:3)
      xjac(1:3,2)=yy(1:3)
      xjac(1:3,3)=zz(1:3)


      call matinv3_uel(xjac,xjacinv,dvol)
      if(dvol.lt.1d-16)then
         write(*,*)'Negative Volume in shapet4'
         stop
      endif

      dvol=1.d0/6.d0*dabs(dvol)

      do inode=1,4
         shp(:,inode)=matmul(xjacinv,shp1(:,inode))
      enddo
      
      return
      end







      subroutine calc_dfg_tet4(ue,nnode,ndofel,mcrd,dfgrd_st,dvol_st,   
     & xjac_st,shp_st,xjbar,ibar,ngauss)
     
      implicit none 

      real(8)::  xjbar,tv,dvol,fac,xjac
      integer::  mcrd, nnode,ngauss, ibar, ndofel,ipt
      real(8)::  dfgrd0(3,3),dfgrd1(3,3),ue(mcrd,nnode),      
     &   xjac_st(ngauss),dfgrd_st(3,3,ngauss),dfgrd_inv(3,3),         
     &   dvol_st(ngauss),shp_st(mcrd,nnode,ngauss),dfgrd(3,3)
     
!      write(*,*)'check inside calc_dfg_tet4 -0'
!	write(*,*)'xjbar is ',xjbar
      xjbar= 0.0
!	write(*,*)'xjbar is ',xjbar
      tv=0.0
!	write(*,*)'tv is ',tv
!	write(*,*)'ngauss is ',ngauss
      
      do ipt=1,ngauss
!	   write(*,*)'check inside calc_dfg_tet4'
         call calc_dfgrd_tet4(ipt,ue,mcrd,nnode,shp_st,dfgrd)
!    write(*,*)'check inside calc_dfg_tet4 -1'
         call matinv3_uel(dfgrd,dfgrd_inv,xjac)
         xjac_st(ipt)=xjac
         dvol=dvol_st(ipt)
         xjbar=xjbar+xjac*dvol
         tv=tv+dvol
         dfgrd_st(1:3,1:3,ipt)=dfgrd(1:3,1:3)
      enddo
      xjbar=xjbar/tv
      if(ibar.eq.1)then
         do ipt=1,ngauss
            fac=xjbar**(1.d0/3.d0)*xjac_st(ipt)**(-1.d0/3.d0)
            dfgrd_st(1:3,1:3,ipt)=dfgrd_st(1:3,1:3,ipt)*fac
         enddo
      endif
      return
      end




      subroutine makespin(w,v)
      implicit real*8(a-h,o-z)
      dimension w(3,3),v(3)
      w=0.d0
      w(1,2)=-v(3)
      w(1,3)=v(2)
      w(2,1)=v(3)
      w(2,3)=-v(1)
      w(3,1)=-v(2)
      w(3,2)=v(1)
      return
      end
!------------------------------------------------
      subroutine matinv3_uel(a,ai,det)
      implicit double precision (a-h,o-z)
!      implicit real*4 (a-h,o-z)
      dimension a(3,3), ai(3,3)
      
       det=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(2,1)*a(1,2)   
     &      *a(3,3)+a(2,1)*a(1,3)*a(3,2)+a(3,1)*a(1,2)*a(2,3)-a(3,1)    
     &     *a(1,3)*a(2,2))
      ai(1,1) =  ( a(2,2)*a(3,3)-a(2,3)*a(3,2))/det
      ai(1,2) = -( a(1,2)*a(3,3)-a(1,3)*a(3,2))/det
      ai(1,3) = -(-a(1,2)*a(2,3)+a(1,3)*a(2,2))/det
      ai(2,1) = -( a(2,1)*a(3,3)-a(2,3)*a(3,1))/det
      ai(2,2) =  ( a(1,1)*a(3,3)-a(1,3)*a(3,1))/det
      ai(2,3) = -( a(1,1)*a(2,3)-a(1,3)*a(2,1))/det
      ai(3,1) =  ( a(2,1)*a(3,2)-a(2,2)*a(3,1))/det
      ai(3,2) = -( a(1,1)*a(3,2)-a(1,2)*a(3,1))/det
      ai(3,3) =  ( a(1,1)*a(2,2)-a(1,2)*a(2,1))/det
      return
      end
!---------------------------------------------------------------
      subroutine dot3(dot,u,v)
      implicit real*8 (a-h,o-z)
      dimension u(3),v(3),w(3)
      w=u*v
      dot=sum(w)
      return
      end
!----------------------------------------------
      subroutine cross(w,u,v)
      implicit real*8 (a-h,o-z)
      dimension w(3),u(3),v(3)
      w(1)=u(2)*v(3)-u(3)*v(2)
      w(2)=u(3)*v(1)-u(1)*v(3)
      w(3)=u(1)*v(2)-u(2)*v(1)
      return
      end
!-----------------------------------------------
      subroutine rudcmp_uel(f,rd,ud) 
      implicit real*8 (a-h,o-z)                                        
!      implicit integer*8 (i-n)                                          
!      real *8 li1,li2,li3,lamda1,lamda2,lamda3
!                 
      dimension rd(3,3),f(3,3),ud(3,3) 

!    
!	write(6,*)'entering ru'
      o3=1.0d0/3.0
      root3=dsqrt(3.d0)
      c11=f(1,1)*f(1,1)+f(2,1)*f(2,1)+f(3,1)*f(3,1)
      c12=f(1,1)*f(1,2)+f(2,1)*f(2,2)+f(3,1)*f(3,2)
      c13=f(1,1)*f(1,3)+f(2,1)*f(2,3)+f(3,1)*f(3,3)
      c23=f(1,2)*f(1,3)+f(2,2)*f(2,3)+f(3,2)*f(3,3)
      c22=f(1,2)*f(1,2)+f(2,2)*f(2,2)+f(3,2)*f(3,2)
      c33=f(1,3)*f(1,3)+f(2,3)*f(2,3)+f(3,3)*f(3,3)
      c1212=c12*c12
      c1313=c13*c13
      c2323=c23*c23
      c2313=c23*c13
      c1223=c12*c23
      c1213=c12*c13
      s11=c22*c33-c2323
      ui1=o3*(c11+c22+c33)
      ui2=s11+c11*c22+c33*c11-c1212-c1313
      ui3=c11*s11+c12*(c2313-c12*c33)+c13*(c1223-c22*c13)
      ui1s=ui1*ui1
      q    =dsqrt(-dmin1(o3*ui2-ui1s,0.d0))                     
      r    =0.5*(ui3-ui1*ui2)+ui1*ui1s
      xmod =q*q*q
      scl1 =.5d0+dsign(.5d0,xmod-1.d-30)                          
      scl2 =.5d0+dsign(.5d0,xmod-dabs(r))                   
      scl0 =dmin1(scl1,scl2)                               
      scl1 =1.-scl0
      sdetm=dacos(r/(xmod+scl1))*o3
      q  =scl0*q
      ct3=q*dcos(sdetm)
      st3=q*root3*dsin(sdetm)
      sdetm=scl1*dsqrt(dmax1(0.0d0,r))
      aa=2.000*(ct3+sdetm)+ui1 
      bb=-ct3+st3-sdetm+ui1
      cc=-ct3-st3-sdetm+ui1                         
      xlamda1=dsqrt(dmax1(aa,0.d0))
      xlamda2=dsqrt(dmax1(bb,0.d0))
      xlamda3=dsqrt(dmax1(cc,0.d0))
      sdetm=xlamda1*xlamda2
      xli1=xlamda1+xlamda2+xlamda3
      xli2= sdetm+xlamda2*xlamda3+xlamda3*xlamda1
      xli3= sdetm*xlamda3/xli1
      s11= c11+xli3
      s22= c22+xli3
      s33= c33+xli3
      s12= c2313-c12*s33
      s13= c1223-s22*c13
      s23=-c2323+s22*s33
      sdetm=1./(xli1*(s11*s23+c12*s12+c13*s13))
      c11=c11+xli2
      c22=c22+xli2
      c33=c33+xli2
      si11=sdetm*s23
      si12=sdetm*s12
      si13=sdetm*s13
      si22=sdetm*( s11*s33-c1313)
      si23=sdetm*(-s11*c23+c1213)
      si33=sdetm*( s11*s22-c1212)
      s12=c12*si12
      s13=c13*si13
      s23=c23*si23
      ui11=c11*si11+s12+s13
      ui22=s12+c22*si22+s23
      ui33=s13+s23+c33*si33
      ui12=c11*si12+c12*si22+c13*si23
      ui13=c11*si13+c12*si23+c13*si33
      ui23=c12*si13+c22*si23+c23*si33
      rd(1,1)=f(1,1)*ui11+f(1,2)*ui12+f(1,3)*ui13
      rd(1,2)=f(1,1)*ui12+f(1,2)*ui22+f(1,3)*ui23
      rd(1,3)=f(1,1)*ui13+f(1,2)*ui23+f(1,3)*ui33
      rd(2,1)=f(2,1)*ui11+f(2,2)*ui12+f(2,3)*ui13
      rd(2,2)=f(2,1)*ui12+f(2,2)*ui22+f(2,3)*ui23
      rd(2,3)=f(2,1)*ui13+f(2,2)*ui23+f(2,3)*ui33
      rd(3,1)=f(3,1)*ui11+f(3,2)*ui12+f(3,3)*ui13
      rd(3,2)=f(3,1)*ui12+f(3,2)*ui22+f(3,3)*ui23
      rd(3,3)=f(3,1)*ui13+f(3,2)*ui23+f(3,3)*ui33 



      do i=1,3
         do j=1,3
            ud(i,j)=0.d0
            do k=1,3
               ud(i,j)=ud(i,j)+rd(k,i)*f(k,j)
            enddo
         enddo
      enddo
!	write(6,*)'exiting ru'
      return
      end
!--------------------------------------------------
      subroutine drotmat(r,a_2d)
      implicit real*8(a-h,o-z)
      dimension r(3,3),a_2d(6,6),a_4d(3,3,3,3)
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a_4d(i,j,k,l)=r(k,i)*r(l,j)
               enddo
            enddo
         enddo
      enddo
      do i=1,3
         a_2d(i,1)=a_4d(i,i,1,1)
         a_2d(i,2)=a_4d(i,i,2,2)
         a_2d(i,3)=a_4d(i,i,3,3)
         a_2d(i,4)=a_4d(i,i,1,2)+a_4d(i,i,2,1)
         a_2d(i,5)=a_4d(i,i,1,3)+a_4d(i,i,3,1)
         a_2d(i,6)=a_4d(i,i,2,3)+a_4d(i,i,3,2)
      enddo
      
      do i=1,3
         a_2d(4,i)=a_4d(1,2,i,i)
         a_2d(5,i)=a_4d(1,3,i,i)
         a_2d(6,i)=a_4d(2,3,i,i)
      enddo
      
      a_2d(4,4)=a_4d(1,2,1,2)+a_4d(1,2,2,1)
      a_2d(4,5)=a_4d(1,2,1,3)+a_4d(1,2,3,1)
      a_2d(4,6)=a_4d(1,2,2,3)+a_4d(1,2,3,2)
      
      a_2d(5,4)=a_4d(1,3,1,2)+a_4d(1,3,2,1)
      a_2d(5,5)=a_4d(1,3,1,3)+a_4d(1,3,3,1)
      a_2d(5,6)=a_4d(1,3,2,3)+a_4d(1,3,3,2)

      a_2d(6,4)=a_4d(2,3,1,2)+a_4d(2,3,2,1)
      a_2d(6,5)=a_4d(2,3,1,3)+a_4d(2,3,3,1)
      a_2d(6,6)=a_4d(2,3,2,3)+a_4d(2,3,3,2)

      a_2d(1:6,4:6)=a_2d(1:6,4:6)/2

      a_2d(4:6,1:6)=2.d0*a_2d(4:6,1:6)
         
      return
      end


      
      



      subroutine uel_disp0(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     & props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     & kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     & lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period,iflag,
     & ngauss,Fp_dgamma_local,euler)
      
      implicit none
      include 'pardis.h'
      
!----------------- variables passed in and out-------------------      
      integer::  ntens, ndofel, nrhs, nsavrs, nprops
      integer::  nnode,jtype,kstep, kinc, jelem
      integer::  ndload, mdload, ngauss, njprop
      integer::  mlvarx,npredf,mcrd,iflag
      integer:: jprops(1),jdltyp(mdload,1),lflags(10)
      parameter (ntens=6)
      real(8)::  rhs(mlvarx,1),svars(msv),energy(maxnp)
      real(8)::  amatrx(ndofel,ndofel),props(maxprops)
      real(8):: coords(mcrd,nnode),u(ndofel)
      real(8):: du(mlvarx,1),v(ndofel),a(ndofel)
      real(8):: time(2), dtime, params(1),adlmag(mdload,1)
      real(8):: ddlmag(mdload,1),predef(2,npredf,nnode)
      real(8):: pnewdt,Fp_dgamma_local(4,12)
      real(8):: period, euler(3)
      
!-----------------other variables used in this code------------------
      real(8):: ske(mdofe,mdofe)
      real(8):: stress(6),ddsdde(ntens,ntens),ddsddt(ntens)
      real(8):: drot(3,3),stran(ntens),dstran(ntens)
      real(8):: predef1(1),dpred(1)
      real(8):: coords_ipt(3),drplde(ntens)
      
      real(8):: ue(mcrd,nnode),dfgrd1(3,3),dfgrd0(3,3)                  
      real(8):: ue_t(maxdim,maxeln),dfgrd1_inv(3,3),pk1_2d(3,3)
      real(8):: pk1_1d(9)  
      real(8):: xb(9,ndofel),f1(ndofel),sigma(3,3)

      real(8):: xe_t(mcrd,nnode),xe_tau(mcrd,nnode)                 
      real(8):: shp0(3,4),shptau(3,4)                              
      real(8):: shp0_st(3,4,ngauss),shptau_st(3,4,ngauss)       
      real(8):: dvol0_st(ngauss),dvoltau_st(ngauss),xb_av(12)

      real(8)::  dfgrd1_st(3,3,ngauss) 
      real(8)::  dfgrd0_st(3,3,ngauss)
      real(8)::  xcmid(maxdim)
      real(8)::  xbar_st(2*maxdim,ndofel,ngauss)
      real(8)::  xc0(maxdim)
      real(8)::  xctau(maxdim)
      real(8)::  xjac_st(ngauss)                  
      real(8)::  xbtr_st(ndofel,ngauss)
      real(8)::  ske1(ndofel,ndofel)
                      
      real(8)::  ske3(ndofel,ndofel)  
      real(8)::    ske2(ndofel,ndofel), stmat(maxdim**2,maxdim**2)
      real(8)::  stmat1(2*maxdim,2*maxdim),xmtemp1(maxdim**2,ndofel)   
      real(8)::    rt_tau(3,3),ut_tau(3,3),dfgrdt(3,3)
      real(8):: crot(2*maxdim,2*maxdim), ftemp(ndofel) 
      real(8)::    xtm(2*maxdim,ndofel),drotm(2*maxdim,2*maxdim)  
      real(8)::    xg_st(maxdim**2,ndofel,ngauss),dfgrd0_inv(3,3)
      real(8)::    xbar(ntens,ndofel),ct(ndofel,ndofel)
      real(8)::  statev(maxstatev),pressb(maxel)

      integer::  ibelem(maxel),ilf(maxel),iface_tet(12)     
      integer::     iface_brick(24),ibnode(4)
      real(8):: xcg(3),uv1(3),uv2(3),uv3(1:3),uv4(1:3)
      real(8):: shp_surf(3,12),fload(12),shp_surf_xi(3,12)               
      real(8)::     shp_surf_eta(3,12),dx_xi(3),dx_eta(3)      
      real(8)::     xel_surf(12),da(3),tempv(3),tempv1(3)
      real(8):: ske_load(12,12),tempv2(3)


      real(8):: xcyc(4,4,3),xjac_tet(3,3),xjac_tet_inv(3,3)              
      real(8):: grad_tet_b(3,4),grad_tet(9,12),dugrad(9)        
      real(8):: dfgrd1_tet(3,3),dfgrd0_tet(3,3),xjac_tet0(3,3)         
      real(8)::    xjac_tet0_inv(3,3),grad_tet0_b(3,4)         
      real(8)::  xcyc0(4,4,3),xbar_tet(6,12),temp2(3,3)
      real(8):: temp2d(3,4),grad_tet0(9,12) 
      

      real(8)::  xjti(3,3),xbtet(6,12),w3(3,3)

      real(8)::  x1(3),x2(3),x3(3),w1(3,3),w2(3,3)
      real(8):: dfgrd_store(3,3)
      real(8)::  FTFjiahao(3,3)
      real(8)::  dfgrd1_t(3,3)
      real(8)::  ckrone_2d(3,3)
      real(8)::  strain_2d(3,3), strain_1d(6)
      integer::  i,j        
 
      
      
      real(8):: deltme,det,a3
      integer:: nsvars

      integer:: jn, ista,iinp
      integer:: ifs,iconv1,ii
      real(8):: a1,a2,dvoltau, dvol_tet
      real(8):: dvol0,dvol0_tet
      integer:: ist1,ipt,inumj,in
      integer:: ist,is,ist2,istv_start
      integer:: nelst,k
      integer::nstatv,nstr,nsv
      real(8):: pr,temp,tvol0,tvoltau
      real(8):: xjbar      
      
      character*80 cmname 
     
!------------GND-----------------------------------------      
      integer:: nodeset(4),choose
      real(8):: B_GND(3,4)
      real(8):: G_GND(3,4),tmp_GND(4),Cmpen_factor1
      real(8):: x_star(3,4), Cmpen_factor2,Cmpen_factor3
!------------ Umat ------------------------------------------
        
      real(8):: dtemp, celent, drpldt,sse
      real(8):: layer,  scd, rpl, spd       
      integer:: ndi, nshr, kspt   
      real(8):: gradient_GND(3,12)
      
      
!----------------------------------------------------------     
      common/press/pressb,ibelem,ilf,iface_tet,iface_brick
      common/conv/iconv1
      common/numj/inumj
         
      
      if(lflags(3).eq.4)then
         amatrx(1:ndofel,1:ndofel)=0.d0
         return
      endif

      inumj=0



   



      ue(1:mcrd,1:nnode)=reshape(u,(/mcrd,nnode/))
      ue_t(1:mcrd,1:nnode)=ue(1:mcrd,1:nnode)-                          
     &    reshape(du(1:mcrd*nnode,1),(/mcrd,nnode/))
        
      xe_tau(1:mcrd,1:nnode)=coords(1:mcrd,1:nnode)+ue(1:mcrd,1:nnode)
      xe_t(1:mcrd,1:nnode)=coords(1:mcrd,1:nnode)+ue_t(1:mcrd,1:nnode)



      tvol0=0.d0
      tvoltau=0.d0

      nstr=6


      ske(1:ndofel,1:ndofel)=0.d0
      f1(1:ndofel)=0.d0
      ngauss=1



      do ipt=1,ngauss
         call shape_tet4(coords,ipt,shp0,dvol0,xc0)
         call shape_tet4(xe_tau,ipt,shptau,dvoltau,xctau)
   
! write(*,*)'coords is ', coords
! write(*,*)'xe_tau is ', xe_tau
   
         dvol0_st(ipt)=dvol0
         dvoltau_st(ipt)=dvoltau
         tvol0=tvol0+dvol0
         tvoltau=tvoltau+dvoltau
      
         shp0_st(1:mcrd,1:nnode,ipt)=shp0(1:mcrd,1:nnode)
         shptau_st(1:mcrd,1:nnode,ipt)=shptau(1:mcrd,1:nnode)

      enddo

        ! write(*,*) 'uel_disp check 3'
        
      do i=1,3
         do j=1,4
            do k=1,4
               xcyc(j,k,i)=xe_tau(i,j)-xe_tau(i,k)
               xcyc0(j,k,i)=coords(i,j)-coords(i,k)
            enddo
         enddo
      enddo


      xjac_tet=xcyc(1:3,4,1:3)

      xjac_tet0=xcyc0(1:3,4,1:3)

      call matinv3_uel(xjac_tet,xjac_tet_inv,dvol_tet)

      call matinv3_uel(xjac_tet0,xjac_tet0_inv,dvol0_tet)


!-------------------------- GND_part--------------------------------------
      if(GND_switch == 1) then              ! if # GND 3
!  compensation factor
      do i=1,4
      do j=1,3
      x_star(j,i)=Fp_dgamma_local(i,j)
      enddo
      enddo
      
      
!      write(4003,*)'x_star is'
!      do i=1,4
!      write(4003,*) (x_star(j,i),j=1,3)
!      enddo
      
      
      if(time(1).le.1.0d0) then
      Cmpen_factor1=1.0
      Cmpen_factor2=1.0
      Cmpen_factor3=1.0
      
      else
      Cmpen_factor1=dsqrt((xe_tau(1,1)-xe_tau(1,4))**2.0 +
     &(xe_tau(2,1)-xe_tau(2,4))**2.0+(xe_tau(3,1)-xe_tau(3,4))**2.0)
     &/dsqrt((x_star(1,1)-x_star(1,4))**2.0 +
     &(x_star(2,1)-x_star(2,4))**2.0+(x_star(3,1)-x_star(3,4))**2.0)
      
      Cmpen_factor2=dsqrt((xe_tau(1,2)-xe_tau(1,4))**2.0 +
     &(xe_tau(2,2)-xe_tau(2,4))**2.0+(xe_tau(3,2)-xe_tau(3,4))**2.0)
     &/dsqrt((x_star(1,2)-x_star(1,4))**2.0 +
     &(x_star(2,2)-x_star(2,4))**2.0+(x_star(3,2)-x_star(3,4))**2.0)
     
      Cmpen_factor3=dsqrt((xe_tau(1,3)-xe_tau(1,4))**2.0 +
     &(xe_tau(2,3)-xe_tau(2,4))**2.0+(xe_tau(3,3)-xe_tau(3,4))**2.0)
     &/dsqrt((x_star(1,3)-x_star(1,4))**2.0 +
     &(x_star(2,3)-x_star(2,4))**2.0+(x_star(3,3)-x_star(3,4))**2.0)
     
      endif 

      
!      write(4003,*) 'Cmpen_factor is', Cmpen_factor1,Cmpen_factor2,
!     &Cmpen_factor3, jelem 
      
      
      
      if (Cmpen_factor1.le.0.0) then
      write(*,*) 'Cmpen_factor1 is not right', Cmpen_factor1, jelem
      stop
      endif
      if (Cmpen_factor2.le.0.0) then
      write(*,*) 'Cmpen_factor2 is not right', Cmpen_factor1, jelem
      stop
      endif
      if (Cmpen_factor3.le.0.0) then
      write(*,*) 'Cmpen_factor3 is not right', Cmpen_factor1, jelem
      stop
      endif
      
      
      
      do i=1,3
      do j=1,4
      G_GND(i,j)=0.0d0
      enddo
      enddo
      
      G_GND(1,1)=1.0d0*Cmpen_factor1
      G_GND(2,2)=1.0d0*Cmpen_factor2
      G_GND(3,3)=1.0d0*Cmpen_factor3
      
      G_GND(1,4)=-1.0d0*Cmpen_factor1
      G_GND(2,4)=-1.0d0*Cmpen_factor2
      G_GND(3,4)=-1.0d0*Cmpen_factor3
      
      
      B_GND=matmul(xjac_tet_inv,G_GND)
      
      
      do choose=1,12
      do i=1,4
      tmp_GND(i)=Fp_dgamma_local(i,choose)
      enddo
      
      do i=1,3
      gradient_GND(i,choose)=0.0d0
      do j=1,4
      gradient_GND(i,choose)=gradient_GND(i,choose)+
     &                                             B_GND(i,j)*tmp_GND(j)
      enddo
      enddo
      enddo
      
      
      endif       ! endif # GND 3
!-------------------------------------------------------------------------	


        ! write(*,*) 'uel_disp check 4'

      temp2=matmul(xjac_tet,xjac_tet_inv)

      dvol_tet=1.d0/6.d0*dvol_tet

      temp2d=0

      temp2d(1,1)=1.d0
      temp2d(2,2)=1.d0
      temp2d(3,3)=1.d0
      temp2d(1:3,4)=-1.d0

      grad_tet_b=matmul(xjac_tet_inv,temp2d)
      grad_tet0_b=matmul(xjac_tet0_inv,temp2d)


      grad_tet=0.0d0
      grad_tet0=0.0d0
      
      
      do i=1,3
         ist1=(i-1)*3
         do j=1,4
            ist2=(j-1)*3+i
            grad_tet(ist1+1:ist1+3,ist2)=grad_tet_b(1:3,j)
            grad_tet0(ist1+1:ist1+3,ist2)=grad_tet0_b(1:3,j)
         enddo
      enddo
        
        ! write(*,*) 'uel_disp check 5'

      xjti=xjac_tet_inv
      xbtet=0.d0

      a1=xjti(1,1)+xjti(1,2)+xjti(1,3)
      a2=xjti(2,1)+xjti(2,2)+xjti(2,3)
      a3=xjti(3,1)+xjti(3,2)+xjti(3,3)

      xbtet(1,1)=xjti(1,1)
      xbtet(1,4)=xjti(1,2)
      xbtet(1,7)=xjti(1,3)
      xbtet(1,10)=-a1


      xbtet(2,2)=xjti(2,1)
      xbtet(2,5)=xjti(2,2)
      xbtet(2,8)=xjti(2,3)
      xbtet(2,11)=-a2

      xbtet(3,3)=xjti(3,1)
      xbtet(3,6)=xjti(3,2)
      xbtet(3,9)=xjti(3,3)
      xbtet(3,12)=-a3

      xbtet(4,1)=xjti(2,1)
      xbtet(4,2)=xjti(1,1)
      xbtet(4,4)=xjti(2,2)
      xbtet(4,5)=xjti(1,2)
      xbtet(4,7)=xjti(2,3)
      xbtet(4,8)=xjti(1,3)
      xbtet(4,10)=-a2
      xbtet(4,11)=-a1

      xbtet(6,2)=xjti(3,1)
      xbtet(6,3)=xjti(2,1)
      xbtet(6,5)=xjti(3,2)
      xbtet(6,6)=xjti(2,2)
      xbtet(6,8)=xjti(3,3)
      xbtet(6,9)=xjti(2,3)
      xbtet(6,11)=-a3
      xbtet(6,12)=-a2

      xbtet(5,1)=xjti(3,1)
      xbtet(5,3)=xjti(1,1)
      xbtet(5,4)=xjti(3,2)
      xbtet(5,6)=xjti(1,2)
      xbtet(5,7)=xjti(3,3)
      xbtet(5,9)=xjti(1,3)
      xbtet(5,10)=-a3
      xbtet(5,12)=-a1

      dugrad=matmul(grad_tet0,u(1:12))

      dfgrd1_tet=reshape(dugrad,(/3,3/))
      
      dugrad=matmul(grad_tet0,u(1:12)-du(1:12,1))
      dfgrd0_tet=reshape(dugrad,(/3,3/))

      do i=1,3
         dfgrd0_tet(i,i)=1.d0+dfgrd0_tet(i,i)
         dfgrd1_tet(i,i)=1.d0+dfgrd1_tet(i,i)
      enddo

      xbar_tet(1,1:12)=grad_tet(1,1:12)
      xbar_tet(2,1:12)=grad_tet(5,1:12)
      xbar_tet(3,1:12)=grad_tet(9,1:12)
      xbar_tet(4,1:12)=grad_tet(2,1:12)+grad_tet(4,1:12)
      xbar_tet(5,1:12)=grad_tet(3,1:12)+grad_tet(7,1:12)
      xbar_tet(6,1:12)=grad_tet(6,1:12)+grad_tet(8,1:12)

      ! write(*,*) 'uel_disp check 6'
      
!	xjbar=0.0
!	! write(*,*)'xjbar is ', xjbar
      call calc_dfg_tet4(ue_t,nnode,ndofel,mcrd,dfgrd0_st,dvol0_st,
     & xjac_st, shp0_st,xjbar,1,ngauss)

      ! write(*,*) 'uel_disp check 6-1'
      
!      xjbar=0.0
!      write(*,*)'xjbar is ', xjbar
      call calc_dfg_tet4(ue,nnode,ndofel,mcrd,dfgrd1_st,dvol0_st,
     & xjac_st,shp0_st,xjbar,1,ngauss)
     
!      write(*,*) 'uel_disp check 6-2'
      call makegrad_tet4(xbar_st,xbtr_st,xg_st,shptau_st,dvoltau_st,
     & nnode,mcrd,ndofel,ngauss)


      ! write(*,*) 'uel_disp check 7'

      xbar_st(1:6,1:12,1)=xbar_tet(1:6,1:12)

      f1(1:ndofel)=0.d0

      nstatv=nsvars
      do ipt=1,ngauss

         dfgrd1(1:3,1:3)=dfgrd1_st(1:3,1:3,ipt)
         dfgrd0(1:3,1:3)=dfgrd0_st(1:3,1:3,ipt)

         if(ipt.eq.1.and.jelem.eq.1)then
            dfgrd_store(1:3,1:3)=dfgrd1(1:3,1:3)
         endif
   
         kstep=1
   
         istv_start=nsvars*(ipt-1)
         do ista=1,nstatv
            statev(ista)=svars(istv_start+ista)
         enddo
   
         ! write(*,*) 'uel_disp check 8'

         call umat(stress,statev,ddsdde,sse,spd,scd,                 
     &        rpl,ddsddt,drplde,drpldt,                                
     &        stran,dstran,time,dtime,temp,dtemp,predef1,dpred,      
     &        cmname,ndi,nshr,ntens,nstatv,props,nprops,coords_ipt,   
     &        drot,pnewdt,celent,dfgrd0,dfgrd1,jelem,ipt,layer,      
     &        kspt,kstep,kinc,euler,gradient_GND)

       
       
!       open(unit=401,file='outtest.txt')
!       write (401,*),('check 1')        
!       write(401,*) 'stress=',stress
!       write(401,*) 'stress is fucking =',stress
       
!       write(401,*) 'stran=', stran

!	    EE=9.6d2
!	    mumu=0.30
!	    zo=0.0d0

!      DDSDDE=reshape((/1.-mumu,mumu,mumu,zo,zo,zo,mumu,1.-mumu,mumu,zo, 
!     &  zo,zo,mumu,mumu,1.-mumu,zo,zo,zo,zo,zo,zo,(1.-2.*mumu)/2.,zo,zo,
!     &  zo,zo,zo,zo,(1.-2.*mumu)/2.,zo,zo,zo,zo,zo,zo,(1.-2.*mumu)/2./),
!     &   (/6,6/))
  

!	 DDSDDE=EE/((1.+mumu)*(1.-2.*mumu))*DDSDDE  

!	 call trans(dfgrd1,dfgrd1_t)
!  	 call mat33(FTFjiahao, dfgrd1_t, dfgrd1,3)
!	 call cal_ckrone_2d(ckrone_2d)

!      strain_2d=0.5*(FTFjiahao-ckrone_2d)	
!      call tr2to1(strain_2d,strain_1d)

!	  do i=1,6
!	  stress(i)=0.0d0
!	  do j=1,6 
!        stress(i)=stress(i)+DDSDDE(i,j)*strain_1d(j)
!        enddo
!       enddo
             
             
             
! -------------------------------------------------------------------------             
       
!       write(*,*) 'time=',  time
!       write(*,*) 'stress=',stress
!       write(*,*) 'strain=', strain_1d
!       write(401,*) 'dfgrd0=', dfgrd0
!       write(401,*) 'dfgrd1=', dfgrd1
      

            ! write(*,*) 'uel_disp check 9'

         if(pnewdt.lt.1.d0)return

         do ista=1,nstatv
            svars(istv_start+ista)=statev(ista)
         enddo

!	   if(iconv1.eq.1)return

         
         xbar(1:nstr,1:ndofel)=xbar_st(1:nstr,1:ndofel,ipt)           
     &     *dvoltau_st(ipt)
         ftemp=matmul(transpose(xbar(1:nstr,1:ndofel)),stress(1:nstr))
         f1(1:ndofel)=f1(1:ndofel)+ftemp(1:ndofel)


         if(inumj.eq.0)then
            do i=1,mcrd
               sigma(i,i)=stress(i)
            enddo
            if(mcrd.eq.2)then
               sigma(1,2)=stress(3)
               sigma(2,1)=stress(3)
            endif
            if(mcrd.eq.3)then
               sigma(1,2)=stress(4)
               sigma(1,3)=stress(5)
               sigma(2,3)=stress(6)
               do i=1,3
                  do j=i,3
                     sigma(j,i)=sigma(i,j)
                  enddo
               enddo
            endif

            stmat(1:mcrd**2,1:mcrd**2)=0.d0

            do ii=1,mcrd
               ist=(ii-1)*mcrd
               do i=1,mcrd
                  do j=1,mcrd
                     stmat(ist+i,ist+j)=sigma(i,j)
                  enddo
               enddo
            enddo
 
         ! write(*,*) 'uel_disp check 10'
   
            stmat1=0.0d0
     
            do i=1,mcrd
               stmat1(i,i)=2*sigma(i,i)
            enddo

            if(mcrd.eq.3)then
               stmat1(1,4)=sigma(1,2)
               stmat1(4,1)=sigma(1,2)
               stmat1(1,5)=sigma(3,1)
               stmat1(5,1)=sigma(3,1)
               stmat1(2,4)=sigma(1,2)
               stmat1(4,2)=sigma(1,2)
               stmat1(2,6)=sigma(2,3)
               stmat1(6,2)=sigma(2,3)
               stmat1(3,6)=sigma(2,3)
               stmat1(6,3)=sigma(2,3)
               stmat1(3,5)=sigma(1,3)
               stmat1(5,3)=sigma(1,3)
            else
               stmat1(1,3)=sigma(1,2)
               stmat1(2,3)=sigma(2,1)
               stmat1(3,3)=0.5d0*(sigma(1,1)+sigma(2,2))
            endif
            if(mcrd.eq.3)then
               stmat1(4,4)=0.5d0*(sigma(1,1)+sigma(2,2))
               stmat1(6,6)=0.5d0*(sigma(2,2)+sigma(3,3))
               stmat1(5,5)=0.5d0*(sigma(1,1)+sigma(3,3))
               stmat1(4,6)=0.5d0*sigma(3,1)
               stmat1(6,4)=0.5d0*sigma(3,1)
               stmat1(4,5)=0.5d0*sigma(3,2)
               stmat1(5,4)=0.5d0*sigma(3,2)
               stmat1(6,5)=0.5d0*sigma(1,2)
               stmat1(5,6)=0.5d0*sigma(1,2)
            endif
            do i=1,ndofel
               do j=1,ndofel
                  ske1(i,j)=ftemp(i)*xbtr_st(j,ipt)
               enddo
            enddo

            ! write(*,*) 'uel_disp check 11'

            xmtemp1(1:mcrd**2,1:ndofel)=
     & matmul(stmat(1:mcrd**2,1:mcrd**2),xg_st(1:mcrd**2,1:ndofel,ipt))
            xmtemp1(1:mcrd**2,1:ndofel)=
     & xmtemp1(1:mcrd**2,1:ndofel)*dvoltau_st(ipt)
     
            ske2(1:ndofel,1:ndofel)=matmul(transpose(xg_st(1:mcrd**2,1:
     &     ndofel,ipt)),xmtemp1(1:mcrd**2,1:ndofel))

            call matinv3_uel(dfgrd0,dfgrd0_inv,det)      

            call mat33(dfgrdt,dfgrd1,dfgrd0_inv,3)

            call rudcmp_uel(dfgrdt,rt_tau,ut_tau)

            call drotmat(rt_tau,drotm)

            if(mcrd.eq.2)then
               do i=1,nstr-1
                  drotm(3,i)=drotm(4,i)
                  drotm(i,3)=drotm(i,4)
               enddo
               drotm(3,3)=drotm(4,4)
            endif


            do i=1,nstr
               do j=1,ndofel
                  ct(i,j)=0.d0
                  do k=1,nstr
                     ct(i,j)=ct(i,j)+ddsdde(i,k)*xbar_st(k,j,ipt)
                  enddo
               enddo
            enddo
 
        ! write(*,*) 'uel_disp check 12'

            crot(1:nstr,1:nstr)=matmul(ddsdde(1:nstr,1:nstr),          
     &     drotm(1:nstr,1:nstr))
            crot(1:nstr,1:nstr)=
     &     (crot(1:nstr,1:nstr)-stmat1(1:nstr,1:nstr))*dvoltau_st(ipt)  

            xtm(1:nstr,1:ndofel)=matmul(crot(1:nstr,1:nstr),  
     &      xbar_st(1:nstr,1:ndofel,ipt))
     
            
            
            ske3(1:ndofel,1:ndofel)=
     &      matmul(transpose(xbar_st(1:nstr,1:ndofel,ipt)),
     &      xtm(1:nstr,1:ndofel))
          
            ske(1:ndofel,1:ndofel)=ske(1:ndofel,1:ndofel)+
     &      ske3(1:ndofel,1:ndofel)+ske2(1:ndofel,1:ndofel)
            
         endif
      enddo

        ! write(*,*) 'uel_disp check 13'


      if(ibelem(jelem).ne.0)then

          ! write(*,*) 'uel_disp check 13-1'
    
         ifs=ibelem(jelem)
         ist=3*(ifs-1)
         ibnode(1:3)=iface_tet(ist+1:ist+3)

         do i=1,3
            xel_surf(3*i-2:3*i)=xe_tau(1:3,ibnode(i))
         enddo
         
         xcg(1:3)=1.d0/4.d0*sum(xe_tau(1:3,1:4),2)
         uv1(1:3)=xcg-xe_tau(1:3,ibnode(1))
         uv2(1:3)=xe_tau(1:3,ibnode(2))-xe_tau(1:3,ibnode(1))
         uv3(1:3)=xe_tau(1:3,ibnode(3))-xe_tau(1:3,ibnode(1))
   
         ! write(*,*) 'uel_disp check 13-2'
   
         call cross(uv4,uv3,uv2)
         call dot3(temp,uv1,uv4)
         is=1
         if(temp.le.0)is=-1
         if(ilf(jelem).eq.1)then
            call dload(jelem,pr,time)
         else
!       Implement Ramping Here
         endif
         pr=is*pr
         fload(1:3)=uv4*pr/3.d0
         fload(4:6)=uv4*pr/3.d0
         fload(7:9)=uv4*pr/3.d0
         fload(1:9)=0.5d0*fload(1:9)
         ske_load=0.d0

         x1=xe_tau(:,ibnode(1))
         x2=xe_tau(:,ibnode(2))
         x3=xe_tau(:,ibnode(3))


        ! write(*,*) 'uel_disp check 14'

         call makespin(w1,x3-x2)

         call makespin(w2,x1-x3)
         call makespin(w3,x2-x1)


         ske_load(1:9,1:9)=0.d0
   
         do i=0,2
            ske_load(3*i+1:3*i+3,1:3)=w1
            ske_load(3*i+1:3*i+3,4:6)=w2
            ske_load(3*i+1:3*i+3,7:9)=w3
         enddo
   
         ske_load(1:9,1:9)=0.5d0*ske_load(1:9,1:9)*pr/3.d0


         do i=1,3
            in=ibnode(i)
            f1(3*in-2:3*in)=f1(3*in-2:3*in)-fload(3*i-2:3*i)
            do j=1,3
               jn=ibnode(j)
               ske(3*in-2:3*in,3*jn-2:3*jn)=ske(3*in-2:3*in,3*jn-2:3*jn)
     &                       -ske_load(3*i-2:3*i,3*j-2:3*j)
            enddo
         enddo

      endif

            ! write(*,*) 'uel_disp check 15'
                
         ske(1:ndofel,1:ndofel)=0.5d0*(ske(1:ndofel,1:ndofel)+       
     &   transpose(ske(1:ndofel,1:ndofel)))

            ! write(*,*) 'uel_disp check 15-1'
            
            
            
         if(pnewdt.lt.1.d0)return


      rhs(1:ndofel,1)=-f1(1:ndofel)
      amatrx(1:ndofel,1:ndofel)=ske(1:ndofel,1:ndofel)

       ! write(*,*) 'uel_disp check main 2'
    
   
      return
      end
 
!----------------------------------------------
      subroutine   PostProcessing(nx,nelx,nprocs_cpfem,nstep,
     & gpxyz,g0xyz, ijk)
      implicit none
      include 'pardis.h'
      integer:: nx,nelx,nprocs_cpfem,nstep, nstp, dim2
      integer:: elemstr, elemend, dummy, ijk(mnelx)
      integer:: i,j,k,ii,jj,kk,n,t
      real(8):: gpxyz(maxcrd),g0xyz(maxcrd)
      real(8):: displacement(nx,3),fpout(nelx,6)
      integer:: nodeset(4), nodelem(maxnodelem,nx)
      real(8):: nodalvalue(nx,6), totvolume, weight,distance
      real(8):: center(nelx,3), xset(4),yset(4),zset(4)
      real(8):: xpos,ypos,zpos
      character*10 fname
      character (LEN=20) ::  filname
      
      dim2=6
      
      
      
      write(*,*)('****************************************')
      write(*,*)('****************************************')
      write(*,*)('Now start to write tecplot output file')
      
      
      open(601,file='disp.out')
      
      
      nstp=(nstep-1)/npost
      write(*,*) 'in post processing',nstp,nstep,npost
      do t=1,nstp
      do i=1,nx
         read(601,*)(displacement(i,j), j=1,3)           !only last step
      enddo
      enddo
      close(601)
      
      write(*,*)('finish reading displacement')
      
      
      do ii=0,nprocs_cpfem-1
      call get_subscript(ii,fname)      
      filname='stressfs'//fname
      open(602,file=filname)
         do t=1,nstp
            read(602,*)elemstr, elemend
            do i=elemstr,elemend
               do j=1,dim2 
                    read(602,*)fpout(i,j), dummy
                    if (dummy /=i)write(*,*)('error code 1')     !still only record last step
               enddo
            enddo
          enddo  
         close(602)
       enddo
       
      write(*,*)('finish reading stress')
      
      
      
! find for each node, which element it belongs to 
      do n=1,nx
      do i=1,maxnodelem
      nodelem(i,n)=0
      enddo
      enddo
      
      do n=1,nx
       k=0
       kk=0
       do i=1,nelx
            do j=1,4
            kk=kk+1
              if(ijk(kk).eq.n)then
              k=k+1
              nodelem(k,n)=i
              endif
            enddo
        enddo
        
        
        if (k.eq.0)then
        write(*,*)'this node is not connected to any element', n
        endif
        
        if (k.gt.maxnodelem) then
        write(*,*)'maxnodelem is not large enough', k
        endif
        
        enddo



       
! calculte the position of each elem's gauss point
      do k=1,nelx
           nodeset(1)=ijk(k*4-3)
           nodeset(2)=ijk(k*4-2)
           nodeset(3)=ijk(k*4-1)
           nodeset(4)=ijk(k*4)
            do j=1,4
            xset(j)=gpxyz(nodeset(j)*3-2)
            yset(j)=gpxyz(nodeset(j)*3-1)
            zset(j)=gpxyz(nodeset(j)*3)
            enddo

       center(k,1)=1.0/4.0*(xset(1)+xset(2)+xset(3)+xset(4))
       center(k,2)=1.0/4.0*(yset(1)+yset(2)+yset(3)+yset(4))
       center(k,3)=1.0/4.0*(zset(1)+zset(2)+zset(3)+zset(4)) 
          
       enddo    
      
      
       write(*,*)'finish calculating the center of each elem'
      
      
! now calculate nodal value of stress  
        do n=1,nx
        xpos=gpxyz(n*3-2)
        ypos=gpxyz(n*3-1)
        zpos=gpxyz(n*3)
        
        
        do j=1,dim2
        nodalvalue(n,j)=0.0d0
        totvolume=0.0
        do i=1,maxnodelem
          k=nodelem(i,n)
          if (k /= 0) then
         distance=sqrt((xpos-center(k,1))**2+(ypos-center(k,2))**2
     &   +(zpos-center(k,3))**2)
         weight=exp(-3*distance)
         nodalvalue(n,j)=nodalvalue(n,j) + weight*fpout(k,j)
         totvolume=totvolume+weight
          endif
        enddo
         nodalvalue(n,j)=nodalvalue(n,j)/totvolume
         enddo 
      enddo  

        write(*,*) ('finish calculating nodal value') 
      
      
       ! at last, output to tecplot format
      
      open(603,position='Append', file='TecplotInput_stress.dat')
        
      write(603,*)('VARIABLES=   "X","Y","Z","dispX",
     & "dispY","dispZ","stress11","stress22",
     & "stress33","stress12","stress13","stress23" ')

      write(603,*) 'ZONE N=  ', nx ,'  E=  ',nelx, 
     &  'F=FEPOINT, ET=TETRAHEDRON'

      do i=1,nx
      write(603,*) g0xyz(i*3-2),g0xyz(i*3-1),g0xyz(i*3),
     & gpxyz(i*3-2)-g0xyz(i*3-2),
     & gpxyz(i*3-1)-g0xyz(i*3-1),
     & gpxyz(i*3-0)-g0xyz(i*3-0),
     & (nodalvalue(i,j),j=1,dim2)
      enddo
      
      do i=1,nelx
      write(603,*) ijk(i*4-3),ijk(i*4-2),ijk(i*4-1),ijk(i*4)
      enddo

      close(603)
      
      write(*,*) ('finish writing output file') 
      
      
      return
      end




