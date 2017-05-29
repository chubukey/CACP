      subroutine disinc_C(v,vbcs1,nbc1,delt,idf,kinc,kstep,
     &  g0xyz,gpxyz,time,ndf,mcrd,nx,nodeid2interid,
     &  nnode_inter,internode_pbc,maxatom_of_one_node,scale_coeff,
     &  natoms,avedisp_x,avedisp_y,avedisp_z,
     &  pbc_box,pbc_strain,
     &  flag_mdskip,flag_inside_iter,couple_iter,
     &  natoms_IVC,atom_IVC,weight_atom_IVC,xc_IVC,
     &  natoms_FVC,atom_FVC,weight_atom_FVC,xc_FVC)

      implicit real*8(a-h,o-z)
      include 'pardis.h'
      dimension v(mdofx),vbcs1(mbound1),nbc1(mbound1),u(3),time(2)
      dimension uold(3),gpxyz(maxcrd)
      dimension coords(3),idf(mbound1),g0xyz(maxcrd)
      dimension xc_IVC(maxnode*3),natoms_IVC(nnode_inter)
      dimension xc_FVC(maxnode*3),natoms_FVC(nnode_inter)
      dimension atom_IVC(maxatom_of_one_node,nnode_inter)  ! 2d array,  (i,j) ith atom id of node j 
      dimension atom_FVC(maxatom_of_one_node,nnode_inter)  ! 2d array,  (i,j) ith atom id of node j 
      dimension weight_atom_IVC(nnode_inter)  ! 2d array,  (i,j) ith atom id of node j 
      dimension weight_atom_FVC(nnode_inter)  ! 2d array,  (i,j) ith atom id of node j 
      integer::mpiinf(mpiinfsizemax),iprocinfo(maxproc)
      integer::isendorrecv(maxproc),natoms,couple_iter
      real(kind=8)::avedisp_x(natoms)
      real(kind=8)::avedisp_y(natoms)
      real(kind=8)::avedisp_z(natoms)
      real(kind=8)::du_ave(3)
      integer::nx,nodeid2interid(nx),internode_pbc(nnode_inter)
      integer::index_inter
      integer::nnode_inter,maxatom_of_one_node
      integer::flag_inside_iter,flag,flag_mdskip                       ! flag = 1 -> no far field disp, only ATC disp, used inside equilibirum iteration
      integer::iflag_topbot  
      real(kind=8) :: ztol
      integer:: iflag_zflat
      integer:: nzlevel,iflag_newzlevel,id_zlevel
      dimension:: nnode_zlevel(10)
      dimension:: dispz_zlevel(10),avedispz_zlevel(10)
      dimension:: zlevel(10),couple_box(3,2)
!      integer:: node_zlevel(100,10)
!      integer:: nnode_zlevel(maxzlevel)
!      integer:: node_zlevel(maxnode_perlevel,maxzlevel)
!      real(kind=8):: dispz_zlevel(maxzlevel),avedispz_zlevel(maxzlevel)
!      real(kind=8)::zlevel(maxzlevel),ztol
      
      common/bounda/nboud1,nboud2,nboud3,nupnch
      common/processorinfo/iproc,nprocs,mpiinfsize,mpiinf,
     &                     iprocinfo,isendorrecv
     

      ztol = 0.01
      iflag_zflat = 0
      nzlevel = 0
      nnode_zlevel(:) = 0
      dispz_zlevel(:) = 0.d0      
      avedispz_zlevel(:) = 0.d0
      zlevel(:) = 0.d0
      couple_box(1,1) = -100.0
      couple_box(1,2) = +100.0
      couple_box(2,1) = -100.0
      couple_box(2,2) = -100.0
      couple_box(3,1) =  0.0
      couple_box(3,2) = +pbc_box(3)

      if(iproc == 0)then
         write(*,*) 'in disinc check flag_mdskip',flag_mdskip
         open(283,position = 'Append',file='ATC_u.out')
         write(283,*) 'at step:',time(1)
 
         open(284,position = 'Append',file='ATC_du.out')
         write(284,*) 'at step:',time(1)

!         write(*,*) 'in disinc_C check gpxyz', gpxyz(1:6)
!         write(*,*) 'in disinc_C check g0xyz', g0xyz(1:6)
      endif
      du_ave(1:3) = 0.0

      flag= flag_inside_iter
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
                    
                 if(abs(coords(3)- couple_box(3,1))<1e-4)then
                    iflag_topbot = -1
                 elseif(abs(coords(3)-couple_box(3,2))<1e-4)then
                    iflag_topbot = 1
                 else
                    iflag_topbot = 0
                 endif
                 uold(1) = gpxyz(node1*3-2) - g0xyz(node1*3-2)
                 uold(2) = gpxyz(node1*3-1) - g0xyz(node1*3-1)
                 uold(3) = gpxyz(node1*3-0) - g0xyz(node1*3-0)
                 if(i == 1) then
                   if(flag_mdskip >= 0)then
                      index_inter = nodeid2interid(node1)
!                 if((node1 == 4).or.(node1 == 3582))then
!                  write(*,*) 'disinc,node1:',node1,'index',index_inter
!                  write(*,*) 'coord:',coords(3),couple_box(3,1:2)
!                  write(*,*) 'flag:',iflag_topbot
!                 endif
        call ATC_update(u,index_inter,natoms,
     &  avedisp_x,avedisp_y,avedisp_z,nx,nodeid2interid,
     &  pbc_box,pbc_strain,iflag_topbot,
     &  maxatom_of_one_node,nnode_inter,internode_pbc,scale_coeff,
     &  natoms_IVC,atom_IVC,weight_atom_IVC,xc_IVC,
     &  natoms_FVC,atom_FVC,weight_atom_FVC,xc_FVC)
                      du_ave(1) = du_ave(1)+abs(u(1)-uold(1))
                      du_ave(2) = du_ave(2)+abs(u(2)-uold(2))
                      du_ave(3) = du_ave(3)+abs(u(3)-uold(3))
                      if(iproc == 0)then
                        write(283,*) index_inter,node1,u(1),u(2),u(3),
     &                  g0xyz(node1*3-2),g0xyz(node1*3-1),g0xyz(node1*3)
                        write(284,*) index_inter,node1,
     &                  u(1)-uold(1),u(2)-uold(2),u(3)-uold(3),
     &                  uold(1),uold(2),uold(3),
     &                  g0xyz(node1*3-2),g0xyz(node1*3-1),g0xyz(node1*3)
                      endif
                    endif
                 endif
                 if(i<4) then
                    v(nn) = u(i) 
                 endif
                 
                 if((i == 3).and.(iflag_zflat==1))then
                    if(nzlevel > 0)then
                       iflag_newzlevel = 1
                       do j = 1,nzlevel
                          if(abs(coords(3)-zlevel(j))<ztol)then
                             nnode_zlevel(j) = nnode_zlevel(j)+1
                             id_zlevel = j
                             iflag_newzlevel = 0
                          endif
                       enddo
                       if(iflag_newzlevel == 1)then
                          nzlevel = nzlevel+1
                          zlevel(nzlevel) = coords(3)
                          nnode_zlevel(nzlevel) = 1
                          id_zlevel = nzlevel
                       endif
                    else
                       nzlevel = 1
                       zlevel(1) = coords(3)
                       nnode_zlevel(1) = 1
                       id_zlevel = 1
                    endif

                    dispz_zlevel(id_zlevel)=dispz_zlevel(id_zlevel)+u(3)
!                    write(*,*) 'dispz updated','level',id_zlevel,
!     &  'n',n,'node',node1,'u',u(3),'newdispz',dispz_zlevel(id_zlevel)
                 endif
              endif              
  
              if(idf(n).eq.5)then        ! load according to position
                 icn=mcrd*(node1-1)
                 do ii=1,mcrd
                    coords(ii)=g0xyz(icn+ii)
                 enddo
!                  call disp_affine(u,time,coords)
                  call disp_affine_dt(u,delt,coords)
                 v(nn)=v(nn)+u(1)
!                write(*,*) 'bc 5 nn',nn,'x,u',coords(1),v(nn),u(1)
              endif
              if(idf(n).eq.6)then        ! load according to position
                 icn=mcrd*(node1-1)
                 do ii=1,mcrd
                    coords(ii)=g0xyz(icn+ii)
                 enddo
!                  call disp_affine(u,time,coords)
                  call disp_affine_dt(u,delt,coords)
                 v(nn)=v(nn)+u(2)
!                write(*,*) 'bc type 6 nn',nn,'coord v',coords(1:3),v(nn)
              endif
            endif
         enddo 
      enddo
      
      if(iproc == 0)then
       close(283)

       write(284,*) 'summary using du_tot',time(1),couple_iter,
     & du_ave(1),du_ave(2),du_ave(3)
       write(284,*) 'summary using du_ave',time(1),couple_iter,
     & du_ave(1)/nnode_inter,du_ave(2)/nnode_inter,du_ave(3)/nnode_inter
       close(284)
      endif

      if(iflag_zflat == 1)then
!      write(*,*) 'nnode_zlevel',nnode_zlevel(1:nzlevel)
!      write(*,*) 'dispz_zlevel',dispz_zlevel(1:nzlevel)
      do i = 1,nzlevel
         avedispz_zlevel(i) = dispz_zlevel(i)/nnode_zlevel(i)
      enddo

      do n=1,nboud1
         node1=nbc1(ndf1*(n-1)+1) 

         do  i=3,3
            index=nbc1(ndf1*(n-1)+1+i)
            
            if(index.eq.1) then
               nn=ndf*(node1-1)+i
               
              if((idf(n).eq.6).or.(idf(n).eq.4))then       ! enforce flag surface
                 icn=mcrd*(node1-1)
                 do ii=1,mcrd
                    coords(ii)=g0xyz(icn+ii)
                 enddo 
                 do j = 1,nzlevel
                    if(abs(coords(3)-zlevel(j))<ztol)then
                       if(ndf.ne.3) then
                          write(*,*) 'avez not applied on z dof'
                       endif
                       v(nn) = avedispz_zlevel(j)
                    endif
                 enddo
              endif

            endif
         enddo
      enddo
      endif
 
      return
      end



      subroutine disp_affine(u,time,coords)
      implicit real*8(a-h,o-z)
      dimension::u(3),time(2),coords(3)
      common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin,uamp,
     &umin
      common/load_mono/sload,st_rate
      common/icyc_flag/icyc
      u(1) = (dexp(st_rate*time(1))-1.d0)*coords(1)
      u(2) = (dexp(st_rate*time(1))-1.d0)*coords(2)
!      u(3) = (dexp(st_rate*time(1))-1.d0)*coords(3)
      u(3) = 0.0
      return
      end


      subroutine disp_affine_dt(u,dt,coords)
      implicit real*8(a-h,o-z)
      dimension::u(3),time(2),coords(3)
      common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin,uamp,
     &umin
      common/load_mono/sload,st_rate
      common/icyc_flag/icyc
      u(1) = (st_rate*dt)*coords(1)
      u(2) = (st_rate*dt)*coords(2)
!      u(3) = (dexp(st_rate*time(1))-1.d0)*coords(3)
      u(3) = 0.0
      return
      end
