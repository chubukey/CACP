      subroutine loadvc_C(f,f_react,vtau,dvtau0,xund,ijk,imid,npropsg,
     & promat2d,deltme,nbc1,ibc,kinc,ndf,N_update,in_st,irowptr,
     & ibindx,val,resid,isize,iflag,ngauss,pnewdt,elem_grain,xeul,
     & grainsize,Fp_dgamma_node,iproc_split,nprocs_cpfem,comm_split,
     & sum_stress_tot,sum_strain_tot,sum_energy_tot,
     & V_FE,V_FE_tot,stiff_media,SC_BOX,flag_SC,g0xyz,
     & npairs_pbc,id_node_pbc,pbc_strain,pbc_box,iflag_PBC,w_penalty)
      
      implicit none
      include 'pardis.h'
      include 'mpif.h'
      
      integer:: N_update
      real(8)::f_react(mdofx),g0xyz(mdofx)
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

      ! variables related to PBC      
      integer:: iflag_exam
      real(8):: w_penalty_local
      ! variables related to SC_calc
      integer:: flag_SC
      real(8):: sum_stress_FE(6),sum_strain_FE(6)
      real(8):: sum_stress_tot(6),sum_strain_tot(6)
      real(8):: V_FE,V_FE_tot
      real(8):: energy_FE,sum_energy_FE,sum_energy_tot
      real(8):: stiff_media(6,6)
      real(8):: elem_stress(6),elem_strain(6),elem_energy,elem_V
      real(8):: xec(3)
      real(8):: SC_BOX(6)
      integer:: flag_in_media

      common/elemnt/nx,nelx,neq,nplane,node
      common/abq/time1,nsvars1,ndofel1,mdload1,npredf1,svars1,mcrd1
      common/processorinfo/iproc,nprocs,mpiinfsize,mpiinf,
     &                     iprocinfo,isendorrecv
      common/el_par/nelst,nelend,nelnox

      real(8)::vbc_elim(mdofx)
      common/elimination/vbc_elim
      
      ! write(*,*)'irowptr ---1 is', irowptr   

! ---- initiation of SC_variables --- !
      do i=1,6
         sum_stress_FE(i) = 0.0
         sum_strain_FE(i) = 0.0
         sum_stress_tot(i) = 0.0
         sum_strain_tot(i) = 0.0
      enddo     
      V_FE = 0.0
      V_FE_tot = 0.0
      ! at end of the loop sum the sum_stress_FE, sum_strain_FE, V_FE
! ----------------------------------- !
!      write(*,*) 'val is',val(1),iflag
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
!      if(iproc_split == 0)then
!         write(*,*) 'loadvc check before loop'
!      endif
      do nel=nelst,nelend

         ngroup=imid(nelx+nel)
         nprops=npropsg(ngroup)

         
         do ii=1,nprops
            props(ii)=promat2d(ii,ngroup)
         enddo
         
         xec(1:3) = 0.d0
         iib=0
         do iia=1,node
            jkia=ijk(node*(nel-1)+iia) 
            jkia1=mcrd*(jkia-1) 
            do iii=1,mcrd
               xeun(iii,iia)=xund(jkia1+iii)
      ! add one line by jiaxi to calc center of elem
               xec(iii) = xec(iii)+g0xyz(jkia1+iii)/node
            enddo
            do ii=1,ndf
               nn=ndf*(jkia-1)+ii
               iib=iib+1
               vtau_e(iib)=vtau(nn)
               dvtau0_e(iib)=dvtau0(nn)
            enddo     
         enddo
           

          if(flag_SC == 1) then
      ! calc flag_in_media based on xcenter(1:3) and SC_BOX(6)
             if((xec(1)>SC_BOX(1)).and.(xec(1)<SC_BOX(2))
     &    .and.(xec(2)>SC_BOX(3)).and.(xec(2)<SC_BOX(4))
     &    .and.(xec(3)>SC_BOX(5)).and.(xec(3)<SC_BOX(6))) then
               flag_in_media = 0
            else
               flag_in_media = 1
            endif

            if(flag_in_media==1)then
              props(1) = stiff_media(1,1)     !c11
              props(2) = stiff_media(1,2)     !c12
              props(3) = stiff_media(1,3)     !c13  in fcc = c12
              props(4) = stiff_media(3,3)     !c33  in fcc = c11
              props(5) = stiff_media(4,4)     !c44  if iso=1/2*(c11-c12)
            endif
         endif


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
       
       do i=1,4
         do j=1,12
      
!       Fp_dgamma_local(i,j)=Fp_dgamma_node((nodeset(i)-1)*12+j,grainID)
           Fp_dgamma_local(i,j)= 0.0d0
       
          enddo
       enddo

!          write(*,*)'check before elst1', iproc , jelem 
          call elst1(jelem,node,xeun,ske1,fe,time,deltme,props,nprops,
     &    nsvars,ndofel,mdload,mlvarx,npredf,mcrd,vtau_e,dvtau0_e,kinc,
     &    iflag,ngauss,pnewdt,nodeset,xeul,Fp_dgamma_local,flag_SC,
     &    elem_strain,elem_stress,elem_V,elem_energy)
!          write(*,*)'check after elst1', iproc , jelem 
           if(pnewdt.lt.pnewdt_loc)pnewdt_loc=pnewdt



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
         
      ! at end of the loop sum the sum_stress_FE, sum_strain_FE, V_FE
         if((flag_in_media == 0).and.(flag_SC ==1)) then
            sum_stress_FE = sum_stress_FE + elem_stress
            sum_strain_FE = sum_strain_FE + elem_strain
            V_FE          = V_FE          + elem_V
         endif
          sum_energy_FE = sum_energy_FE + elem_energy
!         write(*,*) 'elem_E is',elem_energy,sum_energy_FE
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
     
       

!            write(*,*) 'check res NAN error',resid(1:3)
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
!      write(*,*) 'sum_V is',V_FE
!      write(*,*) 'sum_strain_FE is',sum_strain_FE
!      write(*,*) 'sum_stress_FE is',sum_stress_FE
!      write(*,*)'loadvc check  8 ', iproc  
     
      do i=1,N_update
         f(i+in_st)=resid(i)
      enddo
      call mpi_barrier(COMM_SPLIT,ierror)
      call mpi_allreduce(f,f_recv,neq,
     & mpi_double_precision,mpi_sum,COMM_SPLIT,ierror)
      do i=1,neq
         f_react(i)=f_recv(i)
      enddo
!      write(*,*) 'check f_react NAN error',f_react(1:3)

!      if(iproc_split == 0) then
!         do i = 1,20
!            write(*,*) 'ibc(n) is',i,
!     &      ibc(i),vtau(i)
!         enddo
!      endif
  
      iflag_exam = 0
      if(iflag_exam == 1)then
      write(*,*) 'exam vtau and dvtau0'
      write(*,*) 'node 1 z',vtau(3),dvtau0(3),g0xyz(3)
      write(*,*) 'node 7 x',vtau(19),dvtau0(19),g0xyz(19)
      write(*,*) 'node 8 x',vtau(22),dvtau0(22),g0xyz(22)
      write(*,*) 'node 8 y',vtau(23),dvtau0(23),g0xyz(23)
      write(*,*) 'node 8 z',vtau(24),dvtau0(24),g0xyz(24)
      write(*,*) 'node 3588 x',vtau(10762),dvtau0(10762),g0xyz(10762)
      write(*,*) 'node 3588 y',vtau(10763),dvtau0(10763),g0xyz(10763)
      write(*,*) 'node 3588 z',vtau(10764),dvtau0(10764),g0xyz(10764)
      endif


      if(iflag_PBC == 0)then
        call remove_dispdof(ibc,N_update,in_st,irowptr,ibindx,
     &  val,resid,isize,iflag)
       
      elseif(iflag_PBC == 1)then
        if(iflag==0) then
           w_penalty_local = 1.0d4*maxval(val(1:isize))
           call mpi_allreduce(w_penalty_local,w_penalty,1,
     & mpi_double_precision,mpi_max,comm_split,ierror)
!           write(*,*) 'me:',iproc_split,'w:',w_penalty
        endif
           
!        write(*,*) 'before apply PBC iflag',iflag,'w',w_penalty
!        do i=in_st+1,in_st+9
!        do i=in_st+1,in_st+N_update
!          write(*,*) 'i:',i,'irowptr',irowptr(i),irowptr(i+1)-1
!          write(*,*) 'ibindx',ibindx(irowptr(i):irowptr(i+1)-1)
!          write(*,*) 'val',val(irowptr(i):irowptr(i+1)-1)
!        enddo
!      write(*,*) ' res is'
!        write(*,*) 'bf me',iproc,'maxres:',maxval(resid(1:N_update))
        call apply_PBC(N_update,in_st,
     &  irowptr,ibindx,val,resid,isize,iflag,
     &  npairs_pbc,id_node_pbc,PBC_strain,pbc_box,w_penalty,vtau)
        if(iflag_exam == 1)then
          do i=in_st+22,in_st+30
!        do i=in_st+1,in_st+N_update
            write(*,*) 'mid i:',i,'irowptr',irowptr(i),irowptr(i+1)-1
            write(*,*) 'mid ibindx',ibindx(irowptr(i):irowptr(i+1)-1)
            write(*,*) 'mid val',val(irowptr(i):irowptr(i+1)-1)
          enddo
          write(*,*) 'mid res is'
          write(*,*) resid(22:30)
        endif
!        write(*,*) 'mid me',iproc,'maxres:',maxval(resid(1:N_update))
        call remove_dispdof_penalty(ibc,N_update,in_st,
     &  irowptr,ibindx,val,resid,isize,iflag,w_penalty)

!       write(*,*) 'af me',iproc,'maxres:',maxval(resid(1:N_update))
        if(iflag_exam == 1)then
          do i=in_st+22,in_st+30
            write(*,*) 'af i:',i,'irowptr',irowptr(i),irowptr(i+1)-1
            write(*,*) 'af ibindx',ibindx(irowptr(i):irowptr(i+1)-1)
            write(*,*) 'af val',val(irowptr(i):irowptr(i+1)-1)
          enddo
          write(*,*) 'af res is'
          write(*,*) resid(22:30)
        endif
      endif


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


! ------ sum up the sum_stress, sum_strain, V_FE values through all procs ------ !
      call mpi_barrier(COMM_SPLIT,ierror)
      if(flag_SC == 1) then
        call mpi_allreduce(sum_stress_FE,sum_stress_tot,6,
     &   mpi_double_precision,mpi_sum,COMM_SPLIT,ierror)
        call mpi_allreduce(sum_strain_FE,sum_strain_tot,6,
     &   mpi_double_precision,mpi_sum,COMM_SPLIT,ierror)
        call mpi_allreduce(V_FE,V_FE_tot,1,
     &   mpi_double_precision,mpi_sum,COMM_SPLIT,ierror)
      endif
       call mpi_allreduce(sum_energy_FE,sum_energy_tot,1,
     &   mpi_double_precision,mpi_sum,COMM_SPLIT,ierror)
!      if(iproc_split==0)then
!       write(*,*) 'sum_energy_tot is',sum_energy_tot
!      endif
      call mpi_barrier(COMM_SPLIT,ierror)
!      write(*,*) 'sum_V is',V_FE_tot
!      write(*,*) 'sum_strain_FE is',sum_strain_tot,iproc  ! strain value is wrong!!!
!      write(*,*) 'sum_stress_FE is',sum_stress_tot,iproc  

! ------------------------------------------------------------------------------ !

             
      return
      end 


      subroutine apply_PBC(N_update,in_st,
     & irowptr,ibindx,val,resid,isize,iflag,
     & npairs_pbc,id_node_pbc,pbc_strain,pbc_box,weight_pbc,vtau)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension ibc(mdofx)
      common/elemnt/nx,nelx,neq,nplane,node

      dimension ibindx(isize),val(isize),resid(N_update),               
     &            irowptr(N_update+1)
      dimension vtau(neq)
      dimension pbc_du(3)     ! pbc method: x,y,z * dux,duy,duz
      dimension valcheck(2,2),rescheck(2)
      itest = 0                 ! if test mode is on
      kk = 3
      pbc_du(1:3) = 0.0
      pbc_du(3) = pbc_strain(3)*pbc_box(3)
!      write(*,*) 'pbc_du',pbc_du(1:3),pbc_strain(:),pbc_box(:)
!      write(*,*) 'in_st',in_st,'N_update',N_update
      do i_pbc = 1,npairs_pbc(kk)
        do i_case = 1,2
          if(i_case == 1)then
            id_node1 = id_node_pbc(i_pbc,1,kk)
            id_node2 = id_node_pbc(i_pbc,2,kk) 
          elseif(i_case == 2)then
            id_node1 = id_node_pbc(i_pbc,2,kk)
            id_node2 = id_node_pbc(i_pbc,1,kk)
          endif
!          write(*,*) 'id_node1,2',id_node1,id_node2
          do i1 = 1,3
            igot = 0
            idgot1 = 0
            idgot2 = 0
            id_dof1 = id_node1*3-3+i1
            id_dof2 = id_node2*3-3+i1
            du_star = (-vtau(id_dof1)+vtau(id_dof2))
!            if(i1 == 2) then
!               write(*,*) 'dof1,2',id_dof1,id_dof2,du_star,i_case
!            endif
! if node1 = 1, node2 = 2, du_star = -u(1) + u(2), target: u2 - u1 = pbc_du
! resid(1) = resid(1) - pbc_du * w -> resid(1)-(pbc_du-u(2)+u(1))*w 
! resid(1) = resid(1)-(pbc_du-du_star)*w 
! if node1 = 2, node2 = 1, du_star =  u(1) - u(2), target: u2 - u1 = pbc_du
! resid(2) = resid(2) + pbc_du * w -> resid(+)+(pbc_du-u(2)+u(1))*w 
! resid(2) = resid(2)-(pbc_du+du_star)*w 
            if(iflag < 2)then
              if((id_dof1>in_st).and.(id_dof1<=in_st+N_update))then
                if(i_case == 1)then
!                  if((id_dof1-in_st <25).and.(id_dof1-in_st>21))then
!                      write(*,*) 'case1 id',id_dof1,id_dof2
!                      write(*,*) 'vtau',vtau(id_dof1),vtau(id_dof2)
!                      write(*,*) 'du_pbc and star',pbc_du(i1),du_star
!                      write(*,*) 'val',-(pbc_du(i1)-du_star)*weight_pbc
!                  endif

                  resid(id_dof1-in_st) = resid(id_dof1-in_st)
     &                    -(pbc_du(i1)-du_star)*weight_pbc
!                   if(id_dof1-in_st == 11)then
!                      write(*,*) 'resid(11) now is',
!     &                     resid(11),pbc_du(i1),weight_pbc
!                   endif
                elseif(i_case == 2)then
!                  if((id_dof1-in_st <25).and.(id_dof1-in_st>21))then
!                      write(*,*) 'case2 id',id_dof1,id_dof2
!                      write(*,*) 'vtau',vtau(id_dof1),vtau(id_dof2)
!                      write(*,*) 'du_pbc and star',pbc_du(i1),du_star
!                      write(*,*) 'val',-(pbc_du(i1)-du_star)*weight_pbc
!                  endif
                  resid(id_dof1-in_st) = resid(id_dof1-in_st)
     &                    +(pbc_du(i1)+du_star)*weight_pbc
                endif

                if(iflag == 0)then
                 do j = 1,N_update
                  irowno = in_st + j
                  idval = irowno
                  if(irowno == id_dof1)then
                    is = irowptr(j)
                    ie = irowptr(j+1)-1


                    if(j==N_update)then
                      ie = irowptr(N_update+1)
                    endif
                    do jj = is,ie
                      if(ibindx(jj)==id_dof1)then
                        val(jj) = val(jj)+weight_pbc
!                        if(irowno == 1)then      ! for test 
!                          write(*,*) 'id_dof1',id_dof1,jj,val(jj)
!                          write(*,*) 'pair d2:',i_pbc,id_node1,id_node2
!                        endif
                      elseif(id_dof2 == ibindx(jj)) then
                        val(jj) = val(jj)-weight_pbc
!                        if(irowno == 1)then      ! for test 
!                          write(*,*) 'id_dof2',id_dof2,jj,val(jj)
!                          write(*,*) 'pair d1:',i_pbc,id_node1,id_node2
!                        endif
                      endif
                    enddo
                  endif
                 enddo
                endif

              endif    

            elseif(iflag == 2)then
              if((id_dof1>in_st).and.(id_dof1<=in_st+N_update))then
                resid(id_dof1-in_st) = 0.0
              endif
            endif
          enddo
        enddo
      enddo

      end

      subroutine remove_dispdof(ibc,N_update,in_st,irowptr,ibindx,val,  
     &             resid,isize,iflag)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension ibc(mdofx)
      common/elemnt/nx,nelx,neq,nplane,node

      dimension ibindx(isize),val(isize),resid(N_update),               
     &            irowptr(N_update+1)
      real(8)::vbc_elim(mdofx)
      common/elimination/vbc_elim

      xmax_reac=0.d0
      stop
      do i=1,N_update
         irowno=in_st+i

         is=irowptr(i)
         ie=irowptr(i+1)-1
               
         if(i.eq.N_update)ie=irowptr(N_update+1)
         
         if(ibc(irowno).eq.1)then
            if(iflag.eq.0)then
               do j=is,ie
                  val(j)=0.d0
                  if(ibindx(j).eq.irowno)then
                     val(j)=1.d0
                  endif
               enddo
            endif
            resid(i)=0.d0
            goto 140
         endif

         if(iflag.eq.0)then
            do j=is,ie
               icolm=ibindx(j)
               if(ibc(icolm).eq.1)then
                  val(j)=0.d0
               endif
            enddo
         endif

 140  continue
      enddo
            
      end

      subroutine remove_dispdof_penalty(ibc,N_update,in_st,
     &  irowptr,ibindx,val,resid,isize,iflag,weight_penalty)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension ibc(mdofx)
      common/elemnt/nx,nelx,neq,nplane,node

      dimension ibindx(isize),val(isize),resid(N_update),               
     &            irowptr(N_update+1)
      real(8)::vbc_elim(mdofx)
      common/elimination/vbc_elim

      xmax_reac=0.d0
      do i=1,N_update
         irowno=in_st+i

         is=irowptr(i)
         ie=irowptr(i+1)-1
               
         if(i.eq.N_update)ie=irowptr(N_update+1)
         
         if(ibc(irowno).eq.1)then
            if(iflag.eq.0)then
               do j=is,ie
                  if(ibindx(j).eq.irowno)then
                     val(j) = val(j)+weight_penalty
                  endif
               enddo
!               resid(i)=resid(i)+weight_penalty*vbc_elim(irowno)
            elseif(iflag == 2)then      ! if(iflag.eq.2) ?
              resid(i) = 0.0
!              resid(i)=resid(i)+weight_penalty*vbc_elim(irowno)
            endif
            goto 140
         endif

 140  continue
      enddo
      
      ! when not constructing K, need to make residue of
      ! dof related to pbc or dispbc zero

      end

      subroutine remove_dispdof_old(ibc,N_update,in_st,irowptr,ibindx,
     &  val,resid,isize,iflag)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension ibc(mdofx)
      common/elemnt/nx,nelx,neq,nplane,node

      dimension ibindx(isize),val(isize),resid(N_update),               
     &            irowptr(N_update+1)

      xmax_reac=0.d0

      do i=1,N_update
         irowno=in_st+i

         is=irowptr(i)
         ie=irowptr(i+1)-1
               
         if(i.eq.N_update)ie=irowptr(N_update+1)
         
         if(ibc(irowno).eq.1)then
            if(iflag.eq.0)then
               do j=is,ie
                  val(j)=0.d0
                  if(ibindx(j).eq.irowno)then
                     val(j)=1.d0
                  endif
               enddo
            endif
            resid(i)=0.d0
            goto 140
         endif

         if(iflag.eq.0)then
            do j=is,ie
               icolm=ibindx(j)
               if(ibc(icolm).eq.1)then
                  val(j)=0.d0
               endif
            enddo
         endif

 140  continue
      enddo
            
      end
