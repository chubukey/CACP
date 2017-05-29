      subroutine uel_disp_t(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     & props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     & kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     & lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period,iflag,
     & ngauss,Fp_dgamma_local,euler,flag_SC,
     & elem_strain,elem_stress,elem_V,elem_energy)
      
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

! ------------ SC elem stress strain energy volume ------
      integer:: flag_SC
      real(8):: elem_stress(6),elem_strain(6),elem_energy,elem_V
      real(8):: SC_FTF(3,3),SC_FT(3,3)
      real(8):: stiff_media(6,6),xcenter(3)
! --------------------------------------------------------    



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


      call shape_tet4(coords,ipt,shp0,dvol0,xc0)
      elem_V = dvol0
!      write(*,*) 'coords and V',coords,dvol0

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
   
           
         call umat(stress,statev,ddsdde,sse,spd,scd,                 
     &        rpl,ddsddt,drplde,drpldt,                                
     &        stran,dstran,time,dtime,temp,dtemp,predef1,dpred,      
     &        cmname,ndi,nshr,ntens,nstatv,props,nprops,coords_ipt,   
     &        drot,pnewdt,celent,dfgrd0,dfgrd1,jelem,ipt,layer,      
     &        kspt,kstep,kinc,euler,gradient_GND)

       ! for SC model calc elem_V elem_strain elem_stress, elem_energy
!       open(unit=401,file='outtest.txt')
!       write (401,*),('check 1')        
       

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

            ! added by jiaxi for SC, elem_stress,strain and V
!            if(flag_SC==1)then
              SC_FT(1,1) = dfgrd1(1,1)
              SC_FT(2,2) = dfgrd1(2,2)
              SC_FT(3,3) = dfgrd1(3,3)
              SC_FT(1,2) = dfgrd1(2,1)
              SC_FT(2,3) = dfgrd1(3,2)
              SC_FT(3,1) = dfgrd1(1,3)
              SC_FT(2,3) = dfgrd1(3,2)
              SC_FT(3,1) = dfgrd1(1,3)
              SC_FT(1,2) = dfgrd1(2,1)
              call mat33(SC_FTF,SC_FT,dfgrd1,3) 
              elem_strain(1) = (SC_FTF(1,1)-1)*elem_V/2
              elem_strain(2) = (SC_FTF(2,2)-1)*elem_V/2
              elem_strain(3) = (SC_FTF(3,3)-1)*elem_V/2
              elem_strain(4) = (SC_FTF(2,3)-0)*elem_V/2
              elem_strain(5) = (SC_FTF(3,1)-0)*elem_V/2
              elem_strain(6) = (SC_FTF(1,2)-0)*elem_V/2

              elem_stress(1:6) =  stress(1:6)*elem_V
!            endif
             elem_energy= 0.5*(
     &         elem_strain(1)*stress(1)+elem_strain(2)*stress(2)
     &       + elem_strain(3)*stress(3)+elem_strain(4)*stress(4)
     &       + elem_strain(5)*stress(5)+elem_strain(6)*stress(6))
  
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
