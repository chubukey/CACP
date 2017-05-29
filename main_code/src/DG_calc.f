      subroutine f_calc_dist0(nnode_interface,id_node_interface,natoms,
     & dist0,natoms_IVC,atom_IVC,maxnode,g0xyz,x_atoms,
     & maxatom_of_node,maxcrd,scale_coeff,xc_IVC,
     & group_tag,dist0_all)
      implicit none
      real(kind=8)::scale_coeff,xc_IVC(nnode_interface*3)
      real(kind=8)::dist0(3,maxatom_of_node,nnode_interface)
      integer:: maxnode,i,j,k,l,natoms,maxatom_of_node,maxcrd
      integer:: nnode_interface, atom_id, node_id
      integer:: atom_IVC(maxatom_of_node,nnode_interface)
      integer:: id_node_interface(nnode_interface)
      real(kind=8)::x_atoms(natoms*3)
      real(kind=8)::node_coord(3),atom_coord(3)
      real(kind=8)::g0xyz(maxcrd)
      real(kind=8)::center(3),xc_check(3)

      real(kind=8)::xc_all(3)
      real(kind=8)::dist0_all(3,natoms)
      integer::natoms_mobile

      integer(kind=4)::group_tag(natoms)
      integer(kind=4):: natoms_IVC(nnode_interface)
      xc_check = 0.0
!      write(*,*) 'maxatom_of_node is',maxatom_of_node
      do i = 1,nnode_interface
!         write(*,*) 'node_id and natoms_IVC',node_id,natoms_IVC(i)
         node_coord(1:3) = xc_IVC((i*3-2):(i*3))
         do k = 1,3
            center(k) = node_coord(k)*scale_coeff
         enddo
!         write(*,*) 'center pos is',center(:)
         do j = 1,natoms_IVC(i)
            atom_id = atom_IVC(j,i)
            atom_coord(1:3) = x_atoms((atom_id*3-2):(atom_id*3))
!            write(*,*) '2atom pos is',atom_coord(:),atom_id
            do k = 1,3
               dist0(k,j,i) = atom_coord(k)-center(k)
               xc_check(k) = xc_check(k) + dist0(k,j,i)
!           write(*,*) 'dist0(k,j,i) is',dist0(k,j,i),k,j,i,atom_coord(:)
            enddo
         enddo
         if(abs(xc_check(k))>0.1)then
           write(*,*) 'xc_check not pass'
         endif
      enddo
      
      xc_all(1:3) = 0.0
      natoms_mobile = 0
      do j = 1,natoms
         if(group_tag(j) == 0) then
            natoms_mobile = natoms_mobile + 1
            atom_coord(1:3) = x_atoms((j*3-2):(j*3))
            do k = 1,3
               xc_all(k) = xc_all(k)+atom_coord(k)
            enddo
         endif
      enddo      
      
      xc_all(1) = xc_all(1)/natoms_mobile
      xc_all(2) = xc_all(2)/natoms_mobile
      xc_all(3) = xc_all(3)/natoms_mobile

      do j = 1,natoms
         if(group_tag(j) == 0) then
            atom_coord(1:3) = x_atoms((j*3-2):(j*3))
            do k = 1,3
               dist0_all(k,j) = atom_coord(k)-xc_all(k)
            enddo
         endif
      enddo

!       write(*,*) 'nmobile by tag = 0',natoms_mobile,xc_all(1:3)
!      write(*,*) 'check dist0(:,1,1)',dist0(:,1,1)
      return 
      end

      subroutine f_calc_distT(nnode_interface,id_node_interface,natoms,
     & dist0,distT,natoms_IVC,atom_IVC,maxnode,gpxyz,
     & avedisp_x,avedisp_y,avedisp_z,
     & maxatom_of_node,maxcrd,scale_coeff,xc_IVC,MD_boxdim,
     & group_tag,dist0_all,distT_all)
      implicit none
      real(kind=8)::scale_coeff,xc_IVC(nnode_interface*3)
      real(kind=8)::dist0(3,maxatom_of_node,nnode_interface)
      real(kind=8)::distT(3,maxatom_of_node,nnode_interface)
      integer:: maxnode,i,j,k,l,natoms,maxatom_of_node,maxcrd
      integer:: nnode_interface, atom_id, node_id
      integer:: atom_IVC(maxatom_of_node,nnode_interface)
      integer:: id_node_interface(nnode_interface)
      real(kind=8)::avedisp_x(natoms)
      real(kind=8)::avedisp_y(natoms)
      real(kind=8)::avedisp_z(natoms)
      real(kind=8)::node_coord(3),atom_coord(3)
      real(kind=8)::gpxyz(maxcrd)
      real(kind=8)::center(3),MD_boxdim(6),box(3)
      real(kind=8)::boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi
      real(kind=8)::Ddist(3),dist_sign(3),xc_check(3,nnode_interface)

      real(kind=8)::xc_all(3)
      real(kind=8)::dist0_all(3,natoms),distT_all(3,natoms)
      integer::natoms_mobile
      integer(kind=4)::group_tag(natoms)
      integer(kind=4)::natoms_IVC(nnode_interface)
!      write(*,*) 'nnode_interface is',nnode_interface
!      write(*,*) 'check gpxyz',gpxyz(1:10)

! warning -> now not using current center of mass yet!

      boxxlo = MD_boxdim(1) 
      boxxhi = MD_boxdim(2) 
      boxylo = MD_boxdim(3) 
      boxyhi = MD_boxdim(4) 
      boxzlo = MD_boxdim(5) 
      boxzhi = MD_boxdim(6) 
      box(1) = boxxhi-boxxlo
      box(2) = boxyhi-boxylo
      box(3) = boxzhi-boxzlo
!      write(*,*) 'box(1:3) is',box(1:3),MD_boxdim(1:6)
!      write(*,*) 'check distT(:,1,1)',distT(:,1,1)
      do i = 1,nnode_interface
         node_coord(1:3) = xc_IVC((i*3-2):(i*3))
         do k = 1,3
            center(k) = node_coord(k)*scale_coeff
         enddo
         do j = 1,natoms_IVC(i)
            atom_id = atom_IVC(j,i)
            Ddist(1) = avedisp_x(atom_id)
            Ddist(2) = avedisp_y(atom_id)
            Ddist(3) = avedisp_z(atom_id)
!            if(abs(Ddist(3))>0.5*box(3))then
!                write(*,*) 'in dist calc for atom ',atom_id,'in node',i,
!     & distT(3,j,i),dist0(3,j,i),distT(3,j,i)-box(3)*dist_sign(3),box(3)
!            endif
!            Ddist(3) = distT(3,j,i) - dist0(3,j,i)
!            if(Ddist(3)==0)then
!               dist_sign(3) = 0
!            else 
!               dist_sign(3) = Ddist(3)/abs(Ddist(3))
!            endif
!                write(*,*) 'in dist calc for atom ',atom_id,'in node',i,
!     & distT(3,j,i),dist0(3,j,i),distT(3,j,i)-box(3)*dist_sign(3),box(3)
! in dist calc for atom      1536760 in node           8  -100.888969852623       -2.32903668411082        4.71103014737702        105.600000000000     
!                distT(3,j,i) = distT(3,j,i)-box(3)*dist_sign(3)
!            endif
            distT(1,j,i) = dist0(1,j,i)+Ddist(1)
            distT(2,j,i) = dist0(2,j,i)+Ddist(2)
            distT(3,j,i) = dist0(3,j,i)+Ddist(3)
         enddo
      enddo
!       write(*,*) 'check distT(:,1,1)',distT(1,1,1),dist0(1,1,1)
      xc_all(1:3) = 0.0
      natoms_mobile = 0
      do j = 1,natoms
         if(group_tag(j) == -1) then
            natoms_mobile = natoms_mobile + 1
            atom_coord(1) = avedisp_x(j)
            atom_coord(2) = avedisp_y(j)
            atom_coord(3) = avedisp_z(j)
            do k = 1,3
               xc_all(k) = xc_all(k)+atom_coord(k)
            enddo
         endif
      enddo      
      
      xc_all(1) = xc_all(1)/natoms_mobile
      xc_all(2) = xc_all(2)/natoms_mobile
      xc_all(3) = xc_all(3)/natoms_mobile

      do j = 1,natoms
         if(group_tag(j) == 0) then
            atom_coord(1) = avedisp_x(j)
            atom_coord(2) = avedisp_y(j)
            atom_coord(3) = avedisp_z(j)
            do k = 1,3
               distT_all(k,j) = atom_coord(k)-xc_all(k)+dist0_all(k,j)
            enddo
         endif
      enddo
      

      return 
      end 

      subroutine f_calc_Atomic_all(dist0_all,distT_all,ave_strain_MD,
     & natoms,group_tag)
      implicit none
      integer:: i,j,k,m
      integer:: natoms
      real(8):: dist0_all(3,natoms),distT_all(3,natoms)
      real(8):: DG_FA_all(3,3),eta(3,3),omega(3,3),eta_inv(3,3)
      real(8):: dg_tmp(3,3),DG_FA_T(3,3),det_eta
      real(8):: ave_strain_MD(6),C(3,3),E(3,3)
      integer(kind=4):: group_tag(natoms)
      do j = 1,3
         do k = 1,3
            eta(j,k) = 0.0
            omega(j,k) = 0.0
         enddo
      enddo
      do m = 1,natoms
         if(group_tag(m) == 0)then
         do j = 1,3
            do k = 1,3
               eta(j,k) = eta(j,k) + dist0_all(j,m)*dist0_all(k,m)
               omega(j,k) = omega(j,k) + dist0_all(j,m)*distT_all(k,m)
            enddo
         enddo
         endif
      enddo
      call matinv3_uel(eta,eta_inv,det_eta)
!      write(*,*) 'eta is',eta(:,:)
!      write(*,*) 'det_eta is',det_eta
      call matmult(omega,eta_inv,dg_tmp,3,3,3)
      DG_FA_all(:,:) = dg_tmp(:,:)
      do j = 1,3
         do k = 1,3
            DG_FA_T(j,k) = DG_FA_all(k,j)
         enddo
      enddo
      
      call matmult(DG_FA_T,DG_FA_all,C,3,3,3)
      do j = 1,3
         do k = 1,3
            if(j==k) then
               E(j,k) = C(j,k) - 1
            else
               E(j,k) = C(j,k)
            endif
         enddo
      enddo
      ave_strain_MD(1) = E(1,1)
      ave_strain_MD(2) = E(2,2)
      ave_strain_MD(3) = E(3,3)
      ave_strain_MD(4) = E(2,3)
      ave_strain_MD(5) = E(3,1)
      ave_strain_MD(6) = E(1,2)
 
      return
      end

      subroutine f_calc_Atomic(nnode_interface,natoms_IVC,
     & dist0,distT,DG_FA,maxatom_of_node)
      implicit none
      integer:: i,j,k,nnode_interface
      integer(kind=4)::natoms_IVC(nnode_interface)
      integer:: m,n,atom_id,maxatom_of_node
      real(kind=8)::dist0(3,maxatom_of_node,nnode_interface)
      real(kind=8)::distT(3,maxatom_of_node,nnode_interface)
      real(kind=8):: eta(3,3),omega(3,3),eta_inv(3,3)
      real(kind=8):: DG_FA(3,3,nnode_interface)
      real(kind=8):: dg_tmp(3,3)
      real(kind=8):: det_eta
      
      do i = 1,nnode_interface
         do j = 1,3
            do k = 1,3
               eta(j,k) = 0.0
               omega(j,k) = 0.0
            enddo
         enddo
         do m = 1,natoms_IVC(i)
            do j = 1,3
               do k = 1,3
                  eta(j,k) = eta(j,k) + dist0(j,m,i)*dist0(k,m,i)
                  omega(j,k) = omega(j,k) + dist0(j,m,i)*distT(k,m,i)
               enddo
            enddo
         enddo
         call matinv3_uel(eta,eta_inv,det_eta)
!         if(i == 1)then
!           write(*,*) 'eta is',eta(:,:)
!            write(*,*) 'omega is',omega(:,:)
!            write(*,*) 'eta_inv is',det_eta,eta_inv(:,:)
!         endif
         call matmult(omega,eta_inv,dg_tmp,3,3,3)
         DG_FA(:,:,i) = dg_tmp(:,:)
      enddo
      return
      end



 

      subroutine f_calc_Continuum(gxyz,g0xyz,maxcrd,ijk,nx,nelx,mnelx,
     & nnode_interface,id_node_interface,maxnodelem,
     & elem_of_node,n_elem_node,node,DG_FC)
      implicit none
      integer:: i,j,k,node_id,elem_id
      integer:: nx,nelx,mnelx,maxcrd,nnode_interface
      integer:: maxnodelem,node
      integer:: id_node_interface(nnode_interface)
      integer:: ijk(mnelx),nodeset(node)
      real(kind=8)::gxyz(maxcrd),g0xyz(maxcrd),coords_tmp(3,4)
      real(kind=8)::xc(3),xnodal(3)
      integer:: elem_of_node(maxnodelem,nnode_interface)
      integer::n_elem_node(nnode_interface)
      real(kind=8)::dist_v(3),dist1,dist2
      real(kind=8)::weight_gauss,totweight
      real(kind=8)::shapef(3,4),dvol,ue(3,4),dfgrd(3,3),shp_st(3,4,4)
      real(kind=8)::DG_FC(3,3,nnode_interface)
      real(kind=8)::refvector(3,nnode_interface),refv_cutoff
      integer:: nelem_used
      refv_cutoff = 20.0

      if(node.NE.4) then
         write(*,*)'error, not tetra elem'
      endif

      do i = 1,nnode_interface
         do j = 1,3
            do k = 1,3
               DG_FC(j,k,i) = 0.0d0
            enddo
         enddo
         
!         refvector(1,i) = 0.0
!         refvector(2,i) = 0.0
!         refvector(3,i) = 0.0
         node_id = id_node_interface(i)
         totweight = 0
         nelem_used = 0
         xnodal(1) = gxyz(node_id*3-2)
         xnodal(2) = gxyz(node_id*3-1)
         xnodal(3) = gxyz(node_id*3)
         
         do elem_id = 1,nelx
            do k = 1,4
               nodeset(k) = ijk(4*(elem_id-1)+k)
               coords_tmp(1,k) = gxyz(nodeset(k)*3-2)
               coords_tmp(2,k) = gxyz(nodeset(k)*3-1)
               coords_tmp(3,k) = gxyz(nodeset(k)*3)
               ue(1,k) = gxyz(nodeset(k)*3-2)-g0xyz(nodeset(k)*3-2)
               ue(2,k) = gxyz(nodeset(k)*3-1)-g0xyz(nodeset(k)*3-1)
               ue(3,k) = gxyz(nodeset(k)*3)-g0xyz(nodeset(k)*3)
            enddo
            call shape_tet4(coords_tmp,1,shapef,dvol,xc)
            xc(1)=(coords_tmp(1,1)+coords_tmp(1,2)+coords_tmp(1,3)+
     & coords_tmp(1,4))/4
            xc(2)=(coords_tmp(2,1)+coords_tmp(2,2)+coords_tmp(2,3)+
     & coords_tmp(2,4))/4
            xc(3)=(coords_tmp(3,1)+coords_tmp(3,2)+coords_tmp(3,3)+
     & coords_tmp(3,4))/4
            dist_v(1) = -xnodal(1)+xc(1)
            dist_v(2) = -xnodal(2)+xc(2)
            dist_v(3) = -xnodal(3)+xc(3)
            dist2 = dist_v(1)*dist_v(1)+dist_v(2)*dist_v(2)+
     & dist_v(3)*dist_v(3)
            dist1 = sqrt(dist2)
            if(dist1<refv_cutoff)then
               weight_gauss = 1/exp(dist1/10.0)
               totweight = weight_gauss+totweight
               nelem_used = nelem_used + 1
!               refvector(1,i) = refvector(1,i) + dist_v(1)*weight_gauss
!               refvector(2,i) = refvector(2,i) + dist_v(2)*weight_gauss
!               refvector(3,i) = refvector(3,i) + dist_v(3)*weight_gauss
               shp_st(:,:,1) = shapef(:,:)
               call calc_dfgrd_tet4(1,ue,3,4,shp_st,dfgrd)
               DG_FC(:,:,i) = DG_FC(:,:,i)+dfgrd*weight_gauss
            endif
         enddo
         DG_FC(:,:,i) = DG_FC(:,:,i)/totweight
!         refvector(:,i) = refvector(:,i)/totweight

      enddo
      
      return
      end



      subroutine compareF(gpxyz,g0xyz,maxcrd,ijk,nx,nelx,mnelx,
     & nnode_interface,id_node_interface,maxatom_of_one_node,maxnodelem,
     & elem_of_node,n_elem_node,node,DG_FC,iproc,time,
     & dist0,distT,natoms,DG_FA,natoms_IVC,atom_IVC)
      implicit none
      integer:: iproc
      real(kind=8):: time
      integer:: nx,nelx,mnelx,maxcrd,nnode_interface,maxatom_of_one_node
      integer:: maxnodelem,node
      integer:: id_node_interface(nnode_interface),ijk(mnelx)
      real(kind=8)::gpxyz(maxcrd),g0xyz(maxcrd)
      integer:: elem_of_node(maxnodelem,nnode_interface)
      integer::n_elem_node(nnode_interface)
      real(kind=8)::DG_FC(3,3,nnode_interface)
      integer:: natoms
      INTEGER(kind=4):: natoms_IVC(nnode_interface)
      INTEGER:: atom_IVC(maxatom_of_one_node,nnode_interface)  ! 2d array,  (i,j) ith atom id of node j   
      real(kind=8)::x_atoms(natoms*3)
      integer:: dist0(3,maxatom_of_one_node,nnode_interface)
      integer:: distT(3,maxatom_of_one_node,nnode_interface)
       
      integer::i,j,k,node_id,atom_id
      real(kind=8)::DG_FA(3,3,nnode_interface)
      real(kind=8)::ave_FA(3,3),ave_FC(3,3)
      integer:: w_IVC,w2_IVC,zcount
      
      ! scheme: for each interface node, allocate value Fc(n_inter_node),Fa(n_inter_node)
      !         1. find continuum Fc(i) =  sum_n{w(n)*Fc(n)} n is all elems having node i
      !         2. find atomic Fa(i) = sum_n{w(n)*Fa(n)} n is all atoms belong to Voronoi partition of node i
      
      ! notes: try to not store the Fc nor Fa
      if(iproc == 0)then
         write(*,*) 'compareF is called'
         open(264,position = 'Append',file='Deformation_Grad_A.out')
         write(264,*) 'DG_F at current step is ',time
         write(264,*) 'nnode_internode,  DG_F(1,:), DG_F(2,:), DG(3,:)'
         open(266,position = 'Append',file='Deformation_Grad_C.out')
         write(266,*) 'DG_F at current step is ',time
         write(266,*) 'nnode_internode,  DG_F(1,:), DG_F(2,:), DG(3,:)'
      endif
      call f_calc_Continuum(gpxyz,g0xyz,maxcrd,ijk,nx,nelx,mnelx,
     & nnode_interface,id_node_interface,maxnodelem,
     & elem_of_node,n_elem_node,node,DG_FC)
      if(iproc == 0)then
         write(*,*) 'FC is calculated'
      endif
        
      ave_FA(:,:) = 0.0
      ave_FC(:,:) = 0.0
      w_IVC = 0
      w2_IVC = 0
      zcount = 0
      call f_calc_Atomic(nnode_interface,natoms_IVC,
     & dist0,distT,DG_FA,maxatom_of_one_node)
      if(iproc == 0)then
         write(*,*) 'FA is calculated'
      endif
      if(iproc == 0) then
         write(*,*) 'DG_FA(1) is',DG_FA(1,1,1),DG_FA(1,1,4),DG_FA(1,1,2)
         write(*,*) 'DG_FC(1) is',DG_FC(1,1,1),DG_FC(1,1,4),DG_FC(1,1,2)
         do i = 1,nnode_interface
            write(264,*) 'A: ',i,DG_FA(1,:,i),DG_FA(2,:,i),DG_FA(3,:,i)
!            write(*,*) 'A: ',i,DG_FA(1,1,i),DG_FA(2,2,i),DG_FA(3,3,i)
            write(266,*) 'C: ',i,DG_FC(1,:,i),DG_FC(2,:,i),DG_FC(3,:,i)
         enddo
         close(264)
         close(266)
      endif
      do i = 1,nnode_interface
         ave_FA(:,1:2) = ave_FA(:,1:2) + natoms_IVC(i)*DG_FA(:,1:2,i)
         ave_FC(:,:) = ave_FC(:,:) + natoms_IVC(i)*DG_FC(:,:,i)
         w_IVC = w_IVC + natoms_IVC(i)
         if(abs(g0xyz(id_node_interface(i)*3)-0)>1e-2)then
           ave_FA(:,3) = ave_FA(:,3) + natoms_IVC(i)*DG_FA(:,3,i)
           w2_IVC = w2_IVC + natoms_IVC(i)
         else
           zcount = zcount + 1
         endif
      enddo
      ave_FA(:,1:2) = ave_FA(:,1:2)/w_IVC
      ave_FA(:,3) = ave_FA(:,3)/w2_IVC
      ave_FC(:,:) = ave_FC(:,:)/w_IVC
      if(iproc == 0) then
         open(319,position = 'Append',file ='Ave_FA.out')
         open(320,position = 'Append',file ='Ave_FC.out')
         write(*,*) 'ave_Fa and ave_Fc is', ave_FA(1,1),ave_FC(1,1),
     &  ave_FA(2,2),ave_FC(2,2),ave_FA(3,3),ave_FC(3,3),zcount
         write(319,'(100F9.5)') ave_FA(:,:)
         write(320,'(100F9.5)') ave_FC(:,:)
         close(319)
         close(320)
      endif

      return
      end



      subroutine newstiff_LeastSquare(DG_FA,DG_FC,iproc,maxcrd,
     & g0xyz,ninter,id_inter,natoms_IVC,promat)
      implicit none
      integer::ninter,i,j,maxcrd
      integer(kind=4)::natoms_IVC(ninter)
      integer::id_inter(ninter),nodeid
      integer::iproc
      real(8)::g0xyz(maxcrd),node_coord(3)
      real(8)::promat(5)
      real(8)::w
      real(8)::c11,c12,c44
      real(8)::c11n,c12n,c44n
      real(8)::DG_FA(9,ninter),DG_FC(9,ninter)
      real(8)::E_A(6,ninter),E_C(6,ninter)
      real(8)::Einter(3*ninter),Ebar(3*ninter)
      real(8)::Sinter(3*ninter),S_A(6,ninter),S_C(6,ninter)
      real(8)::E_shear(3*ninter),S_shear(3*ninter)
      real(8)::LS_A(2,2),LS_B(2),detA
      real(8)::totE_shear, totS_shear

      LS_A(:,:) = 0.d0
      LS_B(:) = 0.d0
      totE_shear = 0.d0
      totS_shear = 0.d0
      c11 = promat(1)
      c12 = promat(2)
      c44 = promat(5)
      c11 = 241.0e3
      c12 = 150.0e3
      c44 = 127.0e3
      do i = 1,ninter
         nodeid = id_inter(i)
         node_coord(1:3) = g0xyz(nodeid*3-2:nodeid*3)
         call DGF2strain(DG_FA(:,i),E_A(:,i))
         call DGF2strain(DG_FC(:,i),E_C(:,i))
         S_C(1,i) = E_C(1,i)*c11+E_C(2,i)*c12+E_C(3,i)*c12
         S_C(2,i) = E_C(1,i)*c12+E_C(2,i)*c11+E_C(3,i)*c12
         S_C(3,i) = E_C(1,i)*c12+E_C(2,i)*c12+E_C(3,i)*c11
         S_C(4,i) = E_C(4,i)*c44  
         S_C(5,i) = E_C(5,i)*c44  
         S_C(6,i) = E_C(6,i)*c44 
         
         w = natoms_IVC(i)

         if(abs(node_coord(3)-0)<1e-2)then        ! F_A for these nodes are not correct
           E_A(3,i) = E_C(3,i)
         endif
           Sinter(3*i-3+1) = S_C(1,i)*w
           Sinter(3*i-3+2) = S_C(2,i)*w
           Sinter(3*i-3+3) = S_C(3,i)*w

           S_shear(3*i-3+1) = S_C(4,i)*w
           S_shear(3*i-3+2) = S_C(5,i)*w
           S_shear(3*i-3+3) = S_C(6,i)*w
           Einter(3*i-3+1) = E_A(1,i)*w
           Einter(3*i-3+2) = E_A(2,i)*w
           Einter(3*i-3+3) = E_A(3,i)*w

           E_shear(3*i-3+1) = E_A(4,i)*w
           E_shear(3*i-3+2) = E_A(5,i)*w
           E_shear(3*i-3+3) = E_A(6,i)*w
           Ebar(3*i-3+1) = (E_A(2,i)+E_A(3,i))*w
           Ebar(3*i-3+2) = (E_A(3,i)+E_A(1,i))*w
           Ebar(3*i-3+3) = (E_A(1,i)+E_A(2,i))*w
      enddo
      do i =1,ninter*3
         LS_A(1,1) = LS_A(1,1)+Einter(i)*Einter(i)
         LS_A(1,2) = LS_A(1,2)+Einter(i)*Ebar(i)
         LS_A(2,1) = LS_A(2,1)+Einter(i)*Ebar(i)
         LS_A(2,2) = LS_A(2,2)+Ebar(i)*Ebar(i)
         LS_B(1)   = LS_B(1) + Sinter(i)*Einter(i)
         LS_B(2)   = LS_B(2) + Sinter(i)*Ebar(i)
         totE_shear = totE_shear + E_shear(i)
         totS_shear = totS_shear + S_shear(i)
      enddo
      detA = LS_A(1,1)*LS_A(2,2)-LS_A(1,2)*LS_A(2,1)
      if(iproc == 0)then
        write(*,*) 'LS_A',LS_A,'detA',detA,'LS_B',LS_B,'stiff0',c11,c12
       endif
      if(abs(detA)<1e-10) then
         if(iproc == 0)then
           write(*,*) 'determint is 0, LS not working, keep old c11 c12'
         endif
         c11n = c11
         c12n = c12
      else
         c11n = (LS_B(1)*LS_A(2,2)-LS_B(2)*LS_A(1,2))/detA
         c12n = (LS_B(2)*LS_A(1,1)-LS_B(1)*LS_A(1,2))/detA
      endif
      if(totE_shear.eq.0) then
         if(iproc == 0)then
           write(*,*) 'totE_shear is 0, LS not working, keep old c44'
         endif
         c44n = c44
      else
         c44n = totS_shear/totE_shear
      endif
      if(c11n > 50)then
          promat(1) = c11n*0.1 + c11*0.9
      else
          promat(1) = c11
      endif
      if(c12n > 20)then
          promat(2) = c12n*0.1 + c12*0.9
      else
          promat(2) = c12
      endif
!      if(c44n>0)then   
!        promat(5) = c44n
!      endif
      if(iproc==0)then
         write(*,*) 'new elastic modulus,c11 c12',
     &   c11n,c12n,promat(1),promat(2)
      endif
      return
      end

      subroutine DGF2strain(DGF,E)
      implicit none
      real(8)::DGF(9),E(6),F(3,3),C(3,3)
      integer::i,j,k

      F(1,1:3) = DGF(1:3)
      F(2,1:3) = DGF(4:6)
      F(3,1:3) = DGF(7:9)
      C(:,:) = 0.d0
      do i=1,3
         do j = 1,3   
            do k = 1,3   
               C(i,j) = C(i,j)+F(k,i)*F(k,j)
            enddo
         enddo
      enddo
      E(1) = 0.5*(C(1,1)-1)
      E(2) = 0.5*(C(2,2)-1)
      E(3) = 0.5*(C(3,3)-1)
      E(4) = 0.5*(C(1,2))
      E(5) = 0.5*(C(1,3))
      E(6) = 0.5*(C(3,1))
      return
      end
