      subroutine dispRadical(u,coords)
!     user defined displacement field pointing to center
      implicit none
      real(8)::u(3),coords(3)
      real(8)::x0,y0,z0
      real(8)::a

      a = 0.001
      x0 = coords(1)
      y0 = coords(2)
      z0 = coords(3)

      u(1) = -x0/(x0*x0+y0*y0)*a
      u(2) = -y0/(x0*x0+y0*y0)*a
      u(3) = 0.0

      
      return 
      end

!----------------------------------------------
      subroutine dispX(u,time,coords)                                           ! user defined displacment field in x direction
	implicit real*8 (a-h,o-z)
	dimension u(3),time(2), coords(3)
	common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin,uamp,
     &umin
	common/load_mono/sload,st_rate
	common/icyc_flag/icyc

	u(3)=0.0
	u(2)=0.0
      

	if(icyc.eq.1)then

	   umax=uamp+umin
	   
	   if(time(1).le.tramp)then
	      u(1)=umin*time(1)/tramp
	   else
	      tcyc=time(1)-tramp
	      icyc1=tcyc/tperiod
	      t1=tcyc-icyc1*tperiod
	      
	      if(t1.ge.0.d0.and.t1.lt.t_nodwell)then
		 u(1)=umin+t1/t_nodwell*uamp
	      elseif(t1.ge.t_nodwell.and.t1.lt.(t_nodwell+t_dwell))then
		 u(1)=umax
	      else
		 t2=tperiod-t1
		 u1=t2/t_nodwell*uamp
		 u(1)=umin+u1
	      endif
	      
	   endif
	else 
	  
	   x0=coords(1)                              !  g0xyz
	   y0=coords(2)
	   z0=coords(3)
         r=50.0*100.0/time(1)
	!   r=50.0
	!  u(1)=sin(x0/r)*(y0+r)-x0
	!    u(1)=sin(x0/r)*r-x0  
	   u(1)=x0*y0/r
	
	
	!   u(1)=u(1)*time(1)/100.0	   
      !   u(1)=10.d0*st_rate*time(1)

	endif
        write(*,*) 'by calling DISPX u is',u(1)
	return
	end
	
	
	subroutine dispY(u,time,coords)                                           ! user defined displacment field in y direction
	implicit real*8 (a-h,o-z)
	dimension u(3),time(2), coords(3)
	common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin,uamp,
     &umin
	common/load_mono/sload,st_rate
	common/icyc_flag/icyc

	u(3)=0.0
	u(2)=0.0


	if(icyc.eq.1)then

	   umax=uamp+umin
	   
	   if(time(1).le.tramp)then
	      u(1)=umin*time(1)/tramp
	   else
	      tcyc=time(1)-tramp
	      icyc1=tcyc/tperiod
	      t1=tcyc-icyc1*tperiod
	      
	      if(t1.ge.0.d0.and.t1.lt.t_nodwell)then
		 u(1)=umin+t1/t_nodwell*uamp
	      elseif(t1.ge.t_nodwell.and.t1.lt.(t_nodwell+t_dwell))then
		 u(1)=umax
	      else
		 t2=tperiod-t1
		 u1=t2/t_nodwell*uamp
		 u(1)=umin+u1
	      endif
	      
	   endif
	   
	   
	   
	else 
	   x0=coords(1)                              !  g0xyz
	   y0=coords(2)
	   z0=coords(3)
          r=50.0*100.0/time(1)
	 !    r=50.0
      !	   u(1)=(cos(x0/r)-1.0)*(y0+r)	
	!     u(1)=(cos(x0/r)-1.0)*r
	     u(1)=-x0*x0/2.0/r-0.3*y0*y0/2.0/r
	
	!    u(1)=u(1)*time(1)/100.0   
      !    u(1)=10.d0*st_rate*time(1)

	endif

        write(*,*) 'by calling DISPY u is',u(1)

	return
	end
	
!----------------------------------------------------------------------------


          
      subroutine disp2(u,time)
      implicit real*8 (a-h,o-z)
      dimension u(3),time(2)
      common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin,uamp,
     &umin
      common/load_mono/sload,st_rate
      common/icyc_flag/icyc
        
      u(3)=0.0
      u(2)=0.0

      if(icyc.eq.1)then

        umax=uamp+umin
   
        if(time(1).le.tramp)then
           u(1)=umin*time(1)/tramp
        else
          tcyc=time(1)-tramp
          icyc1=tcyc/tperiod
          t1=tcyc-icyc1*tperiod
      
          if(t1.ge.0.d0.and.t1.lt.t_nodwell)then
             u(1)=umin+t1/t_nodwell*uamp
          elseif(t1.ge.t_nodwell.and.t1.lt.(t_nodwell+t_dwell))then
             u(1)=umax
          else
             t2=tperiod-t1
             u1=t2/t_nodwell*uamp
             u(1)=umin+u1
          endif
      
        endif
   
      else
        
         u(1)=-1000.d0*(dexp(st_rate*time(1))-1.d0)
      !   u(1)=10.d0*st_rate*time(1)
        
      endif

!     write(*,*) 'by calling DISP2 u is',u(1)

      return
      end
          
      subroutine disp(u,time)
      implicit real*8 (a-h,o-z)
      dimension u(3),time(2)
      common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin,uamp,
     &umin
      common/load_mono/sload,st_rate
      common/icyc_flag/icyc
        
      u(3)=0.0
      u(2)=0.0

      if(icyc.eq.1)then

	   umax=uamp+umin
	   
	   if(time(1).le.tramp)then
	      u(1)=umin*time(1)/tramp
	   else
	      tcyc=time(1)-tramp
	      icyc1=tcyc/tperiod
	      t1=tcyc-icyc1*tperiod
	      
	      if(t1.ge.0.d0.and.t1.lt.t_nodwell)then
		 u(1)=umin+t1/t_nodwell*uamp
	      elseif(t1.ge.t_nodwell.and.t1.lt.(t_nodwell+t_dwell))then
		 u(1)=umax
	      else
		 t2=tperiod-t1
		 u1=t2/t_nodwell*uamp
		 u(1)=umin+u1
	      endif
	      
	   endif
	   
	else
        
	   u(1)=100.d0*(dexp(st_rate*time(1))-1.d0)
      !    u(1)=10.d0*st_rate*time(1)
!        write(*,*) 'by calling DISP u is',u(1)
	endif

        
	return
	end
	
!-----------------------------------------------
	subroutine dload(jelem,pr,time)
	implicit real*8 (a-h,o-z)
	dimension time(2)
	common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin,uamp,
     &umin
	common/load_mono/sload,st_rate
	common/icyc_flag/icyc



	if(icyc.eq.1)then
	
	   smax=samp+smin
	   
	   if(time(1).le.tramp)then
	      pr=smin*time(1)/tramp
	   else
	      tcyc=time(1)-tramp
	      icyc1=tcyc/tperiod
	      t1=tcyc-icyc1*tperiod
	      
	      if(t1.ge.0.d0.and.t1.lt.t_nodwell)then
		 pr=smin+t1/t_nodwell*samp
	      elseif(t1.ge.t_nodwell.and.t1.lt.(t_nodwell+t_dwell))then
		 pr=smax
	      else
		 t2=tperiod-t1
		 s1=t2/t_nodwell*samp
		 pr=smin+s1
	      endif
	      
	   endif
	   
	else

	   if(time(1).le.tramp)then

	      pr=sload*time(1)/tramp

	   else

	      pr=sload

	   endif
	      
	endif

	pr=-pr

	
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
       subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     &     rpl,ddsddt,drplde,drpldt,
     &     stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     &     cmname,ndi,nshr,ntens,nstatv,props_umat,nprops,coords,
     &     drot,pnewdt,celent,dfgrd0,dfgrd1,noel,npt,layer,
     &     kspt,kstep,kinc,euler,gradient_GND)
     
!      INCLUDE 'ABA_PARAM.INC'
      
        implicit double precision(a-h,o-z)
      
      
        character*80 cmname
C       islip = max no of slip systems allowed             
  
        parameter (itotslips=8,nslip=30)
 
	 ! itotslip is types of dislocation
	 ! 1: 3 basal
	 ! 2: 3 prismatic
	 ! 3: 6 pyramidal (a)
	 ! 4: 12 1st-order pyramidal (c+a)
	 ! 5: 6 2nd-order pyramidal (c+a)
	 ! 6:	
	 ! 7:
	 ! 8:
	
              
        dimension stress(ntens),statev(nstatv),re_tau_t(3,3),
     &     ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     &     stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     &     props_umat(nprops),coords(3,8),drot(3,3),dfgrd0(3,3),
     &     dfgrd1(3,3),
     &     euler(3),g_alpha_t(nslip),g_alpha_tau(nslip),fp_t(3,3),
     &     dsrot(3,nslip),dmrot(3,nslip),s_tan(3,3,nslip),
     &     ddsdde_4d(3,3,3,3),de_gg(3,3,3,3),de_gg_2d(6,6),
     &     gst(3,nslip),gmt(3,nslip),delta_gamma(nslip),fp_tau(3,3),
     &     tr(3,3),tpk_1d_t(6),tpk_1d_tau(6),tot_gamma(nslip),
     &     c_alpha_1d(6,nslip),tau(nslip),fe_tau(3,3),
     &     re_tau(3,3),ue_tau(3,3),tlg(3,3),ics(itotslips),
     &     ga(itotslips),xh0a(itotslips),xs_bara(itotslips),
     &     xs_bar(nslip),xh0(nslip),
     &     nslips(itotslips),
     &     cst(3,nslip),cmt(3,nslip),
     &     xn(nslip),xr(nslip),xri(itotslips),xni(itotslips),cst_1d(6),
     &     q_bcc_hcp(3,3),tr1(3,3),
     &     xh_0(5),xh_inf(5),xg_0(5),xg_inf(5),a_dot_bcc(5),
     &     rate_exp_bcc(5),g0_bcc(5),creep(6),
     &     stress1(6),ddsdde1(6,6),wt(2),q_bcc_hcpt(3,3),
     &     xlast_rot(3,3),tr3(3,3),xn_trb_comp(itotslips),
     &     xr_trb_comp(itotslips),xh0_trb_comp(itotslips),
     &     xs_trb_comp(itotslips),g0_trb_comp(itotslips),
     &     xn_pra_comp(itotslips),xr_pra_comp(itotslips),
     &     xh0_pra_comp(itotslips),xs_pra_comp(itotslips),
     &     g0_pra_comp(itotslips),ps(3),cauchy_2d(3,3),
     &     fe_tau_inv_t(3,3),cbasal(3),bas_tract(3),props(nprops),
     &     chi_t(nslip),chi_tau(nslip),s_gg_2d(6,6),g1_bcc(5),g2_bcc(5)

 
      dimension bn_stress(1)

      dimension ghcp(4),gbcc(1),gpyr(2)
     
     
      dimension ddsdde2(6,6)
     

!----------------------------------------------------------------------
! GND
      real(8):: gradient_GND(3,12),Fp_gradient(3,3,3)
!      real(8):: dgamma_gradient(3,30), rho_GND(30)
!      real(8):: rho_GND_1(3), tmp_GND_g(3), rho_GND_vector(3)
!      real(8):: rho_GND_s_rate(30), rho_GND_et_rate(30)
!      real(8):: rho_GND_en_rate(30)
!      real(8):: type1(3,3),type2(3,3),type3(3,3)
      real(8):: rho_mobile(30), rho_mobile_tt
      real(8):: rho_et_GND(30),rho_en_GND(30)
      real(8):: Curl_Fp(3,3), Fp_t_n(3), rho_s_GND(30)
      real(8):: tmp_m_GND(3), tmp_n_GND(3), tmp_t_GND(3)
      real(8):: ctt(3,nslip), gtt(3,nslip)
      real(8):: rho_F_GND(30), rho_P_GND(30),x_a_b(30,30)
      real(8):: a1,a2,a3,a4,a5,a6 
      real(8):: tau_Pass_GND(30),tau_Cut_GND(30),Nye_rate(3,3) 
      real(8):: a_dot_m(30), f_twin(6), Nye_tensor(3,3)
      real(8):: b_park(3,36), t_park(3,36), A_park(9,36)
      real(8):: AA_T(9,9), A_T(36,9), AtAAt_inv(36,9)
      real(8):: AAt_inv(9,9), rho36(36),Nye_vector(9)
      real(8):: EE, mumu, zo, bur_hcp,sample_size
      real(8):: rho_GND_tt(12),S_critical
      real(8):: screw_cos(12,12),screw_sin(12,12)
      real(8):: edge_cos(12,12),edge_sin(12,12)
      

!twin
      real(8):: twin_R(3,3),twin_X(3,3), twin_X_inv(3,3) 
      real(8):: twin_U_temp(3,3), twin_U(3,3),twin_S(3,3)
      real(8):: twin_Q(3,3), twin_Q_inv(3,3)
      real(8):: cst_tmp(3),cmt_tmp(3),cst_tw(3),cmt_tw(3)
      real(8):: det_twin_X ,s_tan_tw(3,3,6) 
      integer:: grainID,itag

!non_linear
      real(8):: dfgrd0_t(3,3),C_el(3,3),e_el(3,3)
      real(8):: dfgrd_back(3,3),dfgrd_back_t(3,3),tr_t(3,3)    
      real(8):: dfgrd_rot1(3,3),C_el_mid(3,3),C_el_prime(3,3)
      integer(8):: iflag_method
! F = R^T F' R, dfgrd_rot1 = F'*R, dfgrd_back = R^T * dfgrd_rot1
!----------------------------------------------------------------------
!Masoud
!        real(8):: cfp(3,3)
!        real(8):: g_alpha0(nelx,2,nslip),hgnd(30,30)
!        real(8):: ktag,aaaa,bbbb ,aaaa1,bbbb1
        
!-------------------------------------------------------------
      

      common/slip1/tot_gamma

      common/twin_kalidandi/f_twin,f_twin_max, twin_nu      
!      common/GND_part_m/a_dot_m
      common/crys_const/a_dot,rate_exp
      common/time_step/pnewdt1
      common/iter_check/ifirst_check
      common/slips/islip
      common/crys2/ics,ga,xri,xni
      common/crys3/xh0a,xs_bara,nslips
      common/crys4/itag
      common/crys5/xh_0, xh_inf, xg_0,xg_inf, itag_el
      common/hcphard/xn,xr,xh0,xs_bar
      common/crys6/a_dot_bcc,rate_exp_bcc
      common/crys7/a_dot_hcp,rate_exp_hcp
      common/crys8/g0_bcc,g1_bcc,g2_bcc
      common/write_flags/ikinc,ngauss
      common/crys_switch/ihcp,ibcc
      common/poly_inicrys/ipoly_switch
      common/tencomp/isign
      common/kinematic/xc,xd,chi_t,chi_tau


      character*200 fname

       
!      open(unit=4001,file='outtest.txt')
!      open(unit=4002,file='outtest1.txt')
       

      

      ikinc=kinc                                                                        ! increment number
      ngauss=npt                                                                        ! integration point number
      pnewdt1=pnewdt                                                                    ! ratio of time increment

          



!------------added by jiahao   05/02/2013--------------      
      props=0.0
      do i=1,nprops
      props(i)=props_umat(i)
      enddo
      
      nprops1=209
      
      
      
      G_hcp_a=props(nprops1+9)
      xnu_hcp_a=props(nprops1+10)   
      bur_hcp=props(nprops1+11)
      son=props(nprops1+12)
      nGNDswitch=props(nprops1+13)
      sample_size=props(nprops1+14)
      grainID=props(nprops1+16)
      
      call hallpetch(props,nprops1,nprops,grainsize,ghcp,gbcc,gpyr) 
!------------------------------------------------------


c**************     MODIFIED G VALUES FOR HALL-PETCH   ***************
 

c     props for prim alpha 

      if(kinc.le.1.and.kstep.eq.1) then
          call twin_nucleat(twin_nu,dfgrd1,tpk_1d_t,fp_t,s_tan_tw,
     & grainsize,s_gg_2d,S_critical,grainID,1)
      endif
        
 
      do isign=1,2                                                                    !   isign-1----tension,2---compression
         if(isign.eq.1)then
            index=0                                                                         ! all 3 are index for reading props
            index1=0
            index2=0
         endif
     
!     basal
         props(128+index1)=ghcp(1)                                       !  hall-petch results
!     prism     
         props(129+index1)=ghcp(2)
!     <c+a> pyramidal	     
         props(130+index1)=gpyr(1)
!     extension twin
!     props(131+index1)=ghcp(3)
         props(131+index1)=S_critical
!     basal     
         props(132+index1)=ghcp(1)
!     prism   
         props(133+index1)=ghcp(2)
!     <c+a> pyramidal	 
         props(134+index1)=gpyr(1)
!     compression twin
         props(135+index1)=ghcp(4)
         index=6
         index1=8
         index2=2
      enddo
       
       
      
      do i=1,6
         creep(i)=0.d0
         stress1(i)=0.d0
         do j=1,6
            ddsdde1(i,j)=0.d0
         enddo
      enddo





!******************   additional rotation   ************************

      do i=1,3                                                                        !  additional rotation
         do j=1,3
            xlast_rot(i,j)=0.d0
         enddo
      enddo

! -----------for basal loading----------------
!	xlast_rot(1,1)=1.d0
!	xlast_rot(2,2)=cos(0.89)
!	xlast_rot(3,3)=cos(0.89)
!	xlast_rot(2,3)=-sin(0.89)
!	xlast_rot(3,2)=sin(0.89)

!-----------default, no rotation--------------
        xlast_rot(1,1)=1.0d0
        xlast_rot(2,2)=1.0d0
        xlast_rot(3,3)=1.0d0

!----------for rotation in xy plane--------------------
      
!	xlast_rot(1,1)=cos(1.570796326794897*0.0)
!	xlast_rot(2,2)=cos(1.570796326794897*0.0)
!	xlast_rot(1,2)=-sin(1.570796326794897*0.0)
!	xlast_rot(2,1)=sin(1.570796326794897*0.0)
!	xlast_rot(3,3)=1.0d0
!------------------------------------------------




!******************   Start reading properities    *********************

!       now check the phase of the current material point
c.......ipoly=1......Primary Alpha................                             
c.......ipoly=2......Transformed Beta.............


        ipoly_switch=1

!	if(ipoly_switch.eq.1)then
        wt(1)=1.d0                                                                    ! weight for hcp and bcc
        wt(2)=0.d0
!	elseif(ipoly_switch.eq.2)then                                               
!	   wt(1)=0.93d0
!	   wt(2)=0.07d0
!	endif


        ikin_inp=208                                                                     ! both these two are index for reading props
!	ikin_state=220
!	xc=props(ikin_inp)                                                              ! 208, 209: xc,xd: bck stress  related
!	xd=props(ikin_inp+1)

        xc=0.d0
        xd=0.d0

        do itag_el=1,1
 
!	   if(itag_el.eq.1)then
        
           c11=props(1)                                                              ! elastic modulus for ddsdde
           c12=props(2)
           c13=props(3)
           c33=props(4)
           c55=props(5)

           c11 = props(1)
           c12 = props(2)
           c44 = props(5)
    
           c11 = 241.4e3
           c12 = 149.8e3
           c44 = 127.1e3
           
!          if(ipoly(noel).eq.1) ind1=125                                  if phase 1, ind1= 1, if phase 2, ind1=0             
           ind1=125
           itot=itotslips
       
           ind1=ind1+1
           a_dot_hcp=props(ind1)                                                      !  126,  6, gamma 0 dot

           ind1=ind1+1
           rate_exp_hcp=props(ind1)                                                   !  127,  7, 1/m for power law 



           do i=1,itot
              ind1=ind1+1
              ga(i)=props(ind1)                                                              !  128-135, 8-15,   g_0  
           enddo



                                                                                 !  !  136-143, 16-23 g_0 compression
           do i=1,itot
              ind1=ind1+1
              g0_pra_comp(i)=props(ind1)                                       
           enddo                                                              

           do i=1,itot
              ind1=ind1+1
              xni(i)=props(ind1)                                                             ! 144-151, 24-31    n
           enddo




           do i=1,itot
              ind1=ind1+1
              xn_pra_comp(i)=props(ind1)                                                 ! 152-159, 32-39
           enddo



           do i=1,itot
              ind1=ind1+1
              xri(i)=props(ind1)                                                              !   160-167, 40-47      r
           enddo


           do i=1,itot
              ind1=ind1+1
              xr_pra_comp(i)=props(ind1)                                                 !     168-175, 48-55
           enddo




           do i=1,itot
              ind1=ind1+1
              xh0a(i)=props(ind1)                                                             !   56-63, 176- 183    h0
           enddo

           do i=1,itot
              ind1=ind1+1
              xh0_pra_comp(i)=props(ind1)                                                  !  184-191, 64-71
           enddo


           do i=1,itot
              ind1=ind1+1
              xs_bara(i)=props(ind1)                                                          ! 192-199, 72-79    saturation stress
           enddo


           do i=1,itot
              ind1=ind1+1
              xs_pra_comp(i)=props(ind1)                                                   ! 200-207
           enddo



           ihcp=1
           ibcc=0

           call inicrys(icrys,ics,cst,cmt,nslips)
       
            ! nslips(1)=3                             3 basal
            ! nslips(2)=3                             3 prismatic
            ! nslips(3)=6                             6 <c+a> pyramidal
            ! nslips(4)=6                             6 twin 
            ! nslips(5)=3                             3 basal in twin
            ! nslips(6)=3                             3 prismatic in twin
            ! nslips(7)=3                             3  <c+a> pyramidal in twin
            ! nslips(8)=3                             3  <c+a> pyramidal in twin
  

            



!**********************************************************************************	      
!--------------------------       Twin initiation    ----------------------------
        
           if(kinc.le.1.and.kstep.eq.1) then
            
             twin_nu=0.0
             i_twin_active=0 
            
             do isys=1,6
               f_twin(isys)=0.0d0
             enddo
           endif
                        
            
            
            
!************************************************************************************            
!---------------------------       GND    part 1  -----------------------------------
!           Jiahao : I introduced t as orthognal to m and n	      
           do n=1,nslip
              do i=1,3
                 tmp_m_GND(i)=cst(i,n)
                 tmp_n_GND(i)=cmt(i,n)
              enddo
              call crossproduct(tmp_m_GND,tmp_n_GND,tmp_t_GND)
              do i=1,3
                 ctt(i,n)=tmp_t_GND(i)
              enddo
           enddo
!------------------------------------------------------------------------------------	      
!************************************************************************************	      
            ! nslips(1)=3                             3 basal
            ! nslips(2)=3                             3 prismatic
            ! nslips(3)=6                             6 <c+a> pyramidal
            ! nslips(4)=6                             6 twin 
            ! nslips(5)=3                             3 basal in twin
            ! nslips(6)=3                             3 prismatic in twin
            ! nslips(7)=3                             3  <c+a> pyramidal in twin
            ! nslips(8)=3                             3  <c+a> pyramidal in twin
            
            

           ic=0
           do i=1,itot
              do isys=1,nslips(i)
                 xh0(isys+ic)=xh0a(i)
                 xs_bar(isys+ic)=xs_bara(i)
                 xn(isys+ic)=xni(i)
                 xr(isys+ic)=xri(i)
              enddo
              ic=ic+nslips(i)
           enddo

!       Forming the elasticity matrix

           call rclear66(ddsdde)
!              ddsdde(1,1)=c11
!              ddsdde(2,2)=c11
!              ddsdde(3,3)=c33
!              ddsdde(4,4)=0.5d0*(c11-c12)*2
!              ddsdde(5,5)=c55*2
!              ddsdde(6,6)=c55*2
!              ddsdde(1,2)=c12
!              ddsdde(2,1)=c12
!              ddsdde(1,3)=c13
!              ddsdde(3,1)=c13
!              ddsdde(2,3)=c13
!              ddsdde(3,2)=c13
    
     
!  for fcc material parameters are different
!              write(*,*) 'c11 c13 is',c11,c13
           call rclear66(ddsdde)
             
           !! jiaxi modify here
           ! original code, C = de_gg_2d constant, C*E = S, 
           ! new non-linear formulation, de_gg_2d is function of e_el
           call trans(dfgrd0,dfgrd0_t)
           call euler_slip(euler,tr1)            !from euler angle got  rotation matrix
           call mat33(tr,tr1,xlast_rot,3)
           call trans(tr,tr_t)

           iflag_method = 0
           if(iflag_method == 1) then
           ! F = R^T F' R, dfgrd_rot1 = F'*R, dfgrd_back = R^T * dfgrd_rot1
             call mat33(dfgrd_rot1,dfgrd0,tr)
             call mat33(dfgrd_back,tr_t,dfgrd_rot1)
             call trans(dfgrd_back,dfgrd_back_t)
             call mat33(C_el,dfgrd_back_t,dfgrd_back,3)
           elseif(iflag_method == 0)then
           ! C = F^T F = R^T F'^T F' R  = R^T C' R
             call mat33(C_el_prime,dfgrd0,dfgrd0_t,3)
             call mat33(C_el_mid,C_el_prime,tr,3)
             call mat33(C_el,tr_t,C_el_mid,3)
           endif


           do i = 1,3
              do j = 1,3
                 if(i==j) then
                   e_el(i,j) = (C_el(i,j)-1)/2.0
                 else
                   e_el(i,j) =  C_el(i,j)/2.0
                 endif
              enddo
           enddo
           iflag_nonlinear = 1
           if(iflag_nonlinear == 1)then      ! use c111,c112,c123 and e1,e2,e3 to modify current stiffness,de_gg_2d
              
              !12569640
              c11 = 245.0e3
              c12 = 151.8e3
              c44 = 127.1e3*2
              c111 = -880.0e3 
              c112 = -620.0e3
              c123 = -255.0e3
              c144 = -140.0e3
              c166 = -460.0e3         
              c456 = - 36.0e3

              !12569650
              c11 = 245.0e3
              c12 = 151.8e3
              c44 = 127.1e3*2
              c111 = -880.0e3 
              c112 = -620.0e3
              c123 = -255.0e3
              c144 = -200.0e3
              c166 = -320.0e3         
              c456 = - 36.0e3
              ddsdde(1,1)=c11+c111*e_el(1,1)+
     &            c112*(e_el(2,2)+e_el(3,3))
              ddsdde(2,2)=c11+c111*e_el(2,2)+
     &            c112*(e_el(3,3)+e_el(1,1))
              ddsdde(3,3)=c11+c111*e_el(3,3)+
     &            c112*(e_el(1,1)+e_el(2,2))
              ddsdde(1,2)=c12+c123*e_el(3,3)+
     &            c112*(e_el(1,1)+e_el(2,2))
              ddsdde(2,1)=c12+c123*e_el(3,3)+
     &            c112*(e_el(1,1)+e_el(2,2))
              ddsdde(1,3)=c12+c123*e_el(2,2)+
     &            c112*(e_el(1,1)+e_el(3,3))
              ddsdde(3,1)=c12+c123*e_el(2,2)+
     &            c112*(e_el(1,1)+e_el(3,3))
              ddsdde(2,3)=c12+c123*e_el(1,1)+
     &            c112*(e_el(2,2)+e_el(3,3))
              ddsdde(3,2)=c12+c123*e_el(1,1)+
     &            c112*(e_el(2,2)+e_el(3,3))

              ddsdde(4,4)=c44+c144*e_el(1,1)+
     &            c166*(e_el(2,2)+e_el(3,3))
              ddsdde(4,5)=c456*e_el(1,2)
              ddsdde(4,6)=c456*e_el(3,1)
              ddsdde(5,4)=c456*e_el(1,2)
              ddsdde(5,5)=c44+c144*e_el(2,2)+
     &            c166*(e_el(3,3)+e_el(1,1))
              ddsdde(5,6)=c456*e_el(2,3)
              ddsdde(6,4)=c456*e_el(3,1)
              ddsdde(6,5)=c456*e_el(2,3)
              ddsdde(6,6)=c44+c144*e_el(3,3)+
     &            c166*(e_el(1,1)+e_el(2,2))

              ddsdde(1,4)=c144*e_el(2,3)
              ddsdde(1,5)=c166*e_el(3,1)
              ddsdde(1,6)=c166*e_el(1,2)
              ddsdde(2,4)=c166*e_el(2,3)
              ddsdde(2,5)=c144*e_el(3,1)
              ddsdde(2,6)=c166*e_el(1,2)
              ddsdde(3,4)=c166*e_el(2,3)
              ddsdde(3,5)=c166*e_el(3,1)
              ddsdde(3,6)=c144*e_el(1,2)
              ddsdde(4,1)=c144*e_el(2,3)
              ddsdde(5,1)=c166*e_el(3,1)
              ddsdde(6,1)=c166*e_el(1,2)
              ddsdde(4,2)=c166*e_el(2,3)
              ddsdde(5,2)=c144*e_el(3,1)
              ddsdde(6,2)=c166*e_el(1,2)
              ddsdde(4,3)=c166*e_el(2,3)
              ddsdde(5,3)=c166*e_el(3,1)
              ddsdde(6,3)=c144*e_el(1,2)

!              call tr2to4(ddsdde,ddsdde_4d)
!              call rotate4d(ddsdde_4d,tr,de_gg)
!              call tr4to2(de_gg,de_gg_2d)
           else
              c11 = 245.0e3
              c12 = 151.8e3
              c44 = 127.1e3*2
              ddsdde(1:6,1:6) = 0.0
              do ii = 1,6
                 do jj = 1,6 
                    if(ii<=3)then
                       if(jj==ii)then
                          ddsdde(ii,jj) = c11
                       elseif(jj<=3)then
                          ddsdde(ii,jj) = c12
                       endif
                    else
                        if(jj==ii)then
                           ddsdde(ii,jj) = c44
                        endif
                    endif
                 enddo
              enddo
           endif


           call tr2to4(ddsdde,ddsdde_4d)
!           write(*,*) 'old c_4d:',ddsdde_4d(1,1,1,1),ddsdde_4d(1,2,1,2),
!     & ddsdde_4d(2,2,1,1),ddsdde_4d(2,2,2,2),ddsdde_4d(2,2,1,3)
!	      call mat33(tr,xlast_rot,tr1,3)
           call rotate_crys_vector(cst,tr,gst)
           call rotate_crys_vector(cmt,tr,gmt)
           call rotate_crys_vector(ctt,tr,gtt)
!          write(*,*) 'tr 11 22 33 12 23 31 21 32 13 is',
!     & tr(1,1),tr(2,2),tr(3,3),tr(1,2),tr(2,3),
!     & tr(3,1),tr(2,1),tr(3,2),tr(1,3)
           call rotate4d(ddsdde_4d,tr,de_gg)
!           write(*,*) 'new c_4d:',de_gg(1,1,1,1),de_gg(1,1,2,2),
!     & de_gg(2,2,1,1),de_gg(2,2,2,2),de_gg(2,2,1,3)
           call tr4to2(de_gg,de_gg_2d)
!--------------------------------------------
!           write(*,*) 'de_gg is',de_gg_2d
!     for feTfe
!              
           do j=1,6
              do i=1,6
                 s_gg_2d(i,j)=de_gg_2d(i,j)                                        
              enddo
           enddo
              
           call matinv(s_gg_2d,6,6)                                               ! inverse of stiffness
!----------------------------------------------    
!	   endif
           ihcp=1
           ibcc=0

!------------------------------------------------------------------------
           do isys=1,islip                                                                ! islip = 30 for hcp
              do j=1,3
                 do i=1,3
                    s_tan(i,j,isys)=gst(i,isys)*gmt(j,isys)                                      ! schmid tensor
                 enddo
              enddo
           enddo

! ***************************************************************************************************
!     Initiallization 
! ***************************************************************************************************
   
          if(kinc.le.1.and.kstep.eq.1) then                                             !  first step first increment

             do isys=1,islip
                chi_t(isys)=0.d0                                                               !   for kinematic hardening   
             enddo

             if(itag_el.eq.1)then
                ic=0
                do i=1,itot
                   do isys=1,nslips(i)
                      g_alpha_t(isys+ic)=ga(i)
                   enddo
                 ic=ic+nslips(i)
                enddo
             endif

             do isys=1,islip
                tot_gamma(isys)=0.d0
             enddo 
             call rclear(tpk_1d_t,6)
             call cal_ckrone_2d(fp_t) 
             w_p=0.0d0        !Masoud_crack

!-----------------initialize GND---------------------------------	     

             do isys=1,islip
                rho_s_GND(isys)=0.0d0
                rho_et_GND(isys)=0.0d0
                rho_en_GND(isys)=0.0d0
                tau_Cut_GND(isys)=0.0d0
                tau_Pass_GND(isys)=0.0d0
             enddo
           
!---------------initialize twin volume fraction-----------------           
             do isys=1,6 
                f_twin(isys)=0.0d0
             enddo
!-----------------------------------------------------------------
          endif

!*************************************************************************************************
!      After the first step                             
!*************************************************************************************************

         if(kinc.gt.1.or.kstep.gt.1)then
!-------------------------------------------------------------	        
	      if(itag_el.eq.1)index2=0
	      do isys=1,islip
		 index2=index2+1
		 g_alpha_t(isys)=statev(index2)                                                 ! statev 1-30 :   resistance
	      enddo
	      do isys=1,islip
		 index2=index2+1
		 tot_gamma(isys)=statev(index2)                                                 ! statev 31-60:  tot_gamma
	      enddo
	      do i=1,6
		 index2=index2+1
		 tpk_1d_t(i)=statev(index2)                                                     ! statev 61-66:   1d pk stress 
	      enddo
	      do i=1,3
		 do j=1,3
		    index2=index2+1
		    fp_t(i,j)=statev(index2)                                                    ! statev 67-75: fp 
		 enddo
	      enddo
	      ikin_init_state=0
	      ikin_state=75
	      do isys=1,islip
              chi_t(isys)=statev(ikin_init_state+ikin_state+isys)                           !  statev 76-105: back stress related
           enddo
           w_p=statev(299)      !Masoud_crack

!*************************** twin nucleation *********************************

           twin_nu=statev(289)
           i_twin_active=statev(287)
           gsize=grainsize

           if (twin_nu.eq.0.0) then
           
           do isys=1,6
           do i=1,3
           do j=1,3
           s_tan_tw(i,j,isys)=s_tan(i,j,12+isys)
           enddo
           enddo
           enddo
           
           call twin_nucleat(twin_nu,dfgrd1,tpk_1d_t,fp_t,s_tan_tw,
     &     gsize,s_gg_2d,S_critical,grainID,2)
           endif 
           

           do isys=1,6
           f_twin(isys)=statev(280+isys)
           enddo
           
           

           if ((twin_nu.eq.1.0).and.(i_twin_active.eq.0)) then
           
           f_twin_max= 1.0d-15                                                        !  set here to switch on/off slip in twin

           do isys=1,6 
           if (f_twin(isys).gt.f_twin_max) then
           i_twin_active=isys
           f_twin_max=f_twin(isys)
           endif
           enddo           
           
           endif
       
!   if no twin system has large twin volume fraction,
!  then set max_twin volume fraction back to 0 

           if(i_twin_active.eq.0)f_twin_max=0.0d0
           
           
!-------------------------  Twin rotation    --------------------------------           
           
!   if max_twin_volume fraction is large enough, then active slip inside twin region
        if ((twin_nu.eq.1.0).and.(i_twin_active.ne.0)) then
           
           twin_shear=0.1089
           
           twin_R=reshape((/-1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,
     & 0.0d0,0.0d0,1.0d0/),(/3,3/))
     
          kk=i_twin_active+12
        twin_X(1,1)=cst(1,kk)
        twin_X(1,2)=cst(2,kk)
        twin_X(1,3)=cst(3,kk)
        twin_X(3,1)=cmt(1,kk)
        twin_X(3,2)=cmt(2,kk)
        twin_X(3,3)=cmt(3,kk)
        twin_X(2,1)=cst(3,kk)*cmt(2,kk)-cst(2,kk)*cmt(3,kk)
        twin_X(2,2)=cst(1,kk)*cmt(3,kk)-cst(3,kk)*cmt(1,kk)
        twin_X(2,3)=cst(2,kk)*cmt(1,kk)-cst(1,kk)*cmt(2,kk)
        call matinv3(twin_X, twin_X_inv,det_twin_X)
        call mat33(twin_U_temp,twin_X_inv, twin_R,3)
        call mat33(twin_U,twin_U_temp,twin_X,3)
        
        do i=1,3
        do j=1,3
        if (i==j) then
        akd=1.0d0
        else
        akd=0.0d0
        endif
        twin_S(i,j)=akd+twin_shear*cst(i,kk)*cmt(j,kk)
        enddo
        enddo
        
        
        call mat33(twin_Q, twin_U, twin_S, 3)

        call matinv3(twin_Q, twin_Q_inv, det_twin_Q)
        
       ! reorintate slip plane and slip system 
        
        do isys=1,12
        
        do i=1,3
        cst_tmp(i)=cst(i,isys)
        cmt_tmp(i)=cmt(i,isys)
        enddo
        
        
        do i=1,3
        cst_tw(i)=0.0d0
        cmt_tw(i)=0.0d0
        do j=1,3
            cst_tw(i)=cst_tw(i)+twin_Q(i,j)*cst_tmp(j)
            cmt_tw(i)=cmt_tw(i)+twin_Q_inv(j,i)*cmt_tmp(j)
        enddo
        enddo
        
        
        do i=1,3
        cst(i,isys+18)=cst_tw(i)
        cmt(i,isys+18)=cmt_tw(i)
        
        enddo
        
        enddo                                          ! isys =1,12
           
        endif                                          ! if f_twin_max.gt.1.0d-5 
                
        
        

!********************************   GND   *************************************      

      
      
      
      do i=1,3
      do j=1,3
      n=3+(i-1)*3+j                                                 
      do k=1,3
      Fp_gradient(j,i,k)=gradient_GND(k,n)/sample_size                             !  Fp_transpose
      enddo
      enddo
      enddo
      
      
   !   do i=1,30  
   !   do k=1,3                                                          !  30 = nslip
   !   dgamma_gradient(i,k)=gradient_GND(k,i)/sample_size
   !   enddo
   !   enddo
      
      
      call curlmat(Fp_gradient,Curl_Fp)
     
       rho_GND_total=0.0d0
  
      do i=1,3
      do j=1,3
      
      Nye_tensor(j,i)=-1.0/bur_hcp*Curl_Fp(i,j)                       ! Nye_ transpose
       
      enddo
      enddo
   
  
  
 !     if (noel.eq.100) then
 !     write(4001,*) 'sample_size is', sample_size
 !     write(4001,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
 !     write(4001,*) 'time is', time
 !     write(4001,*) '***************twin is*************'
 !     do i=1,6
 !     write(4001,*) f_twin(i)
 !     enddo      
 !     endif

  
  
  
      
!      do i=1,30                                                                   !               1 
      
!      call matmul(Curl_Fp,gmt(:,i),rho_GND_vector,3,3,1)
!      do j=1,3
!      rho_GND_vector(j)=statev(150+i)*rho_GND_vector(j)
!      enddo
      
!      do k=1,3
!      tmp_GND_g(k)=dgamma_gradient(i,k)
!      enddo
      
!      call matmul(fp_t,gmt(:,i),Fp_t_n,3,3,1)
      
      
      !  rho_GND_1 is the second part of curl vector
!      call crossproduct(tmp_GND_g, Fp_t_n, rho_GND_1)
      
!      do j=1,3
!      rho_GND_vector(j)=rho_GND_vector(j)+rho_GND_1(j)
!      enddo
      
      ! now we have rho_GND_vector, and will calculate Nye_density_tensor
      
!      do k=1,3
!      do j=1,3
!      Nye_rate(k,j)=-1.0/1.0E-10*gst(k,i)*rho_GND_vector(j)
!      enddo
!      enddo
      
      ! now we calculate the dislocation density for 3 types
      
!      do j=1,3
!      tmp_m_GND(j)=gst(j,i)
!      tmp_n_GND(j)=gmt(j,i)
!      tmp_t_GND(j)=gtt(j,i)
!      enddo
      
!      call tensorproduct(tmp_m_GND,tmp_m_GND,type1)
!      call tensorproduct(tmp_m_GND,tmp_t_GND,type2)
!      call tensorproduct(tmp_m_GND,tmp_n_GND,type3)
      
!      call doubledot(Nye_rate,type1,tmp_rho)
!      rho_GND_s_rate(i)=dabs(tmp_rho)
!      call doubledot(Nye_rate,type2,tmp_rho)
!      rho_GND_et_rate(i)=dabs(tmp_rho)
!      call doubledot(Nye_rate,type3,tmp_rho)
!      rho_GND_en_rate(i)=dabs(tmp_rho)
       

      
      ! next is the 3 types of dislocation density
      
!      rho_s_GND(i)=statev(180+i)+ rho_GND_s_rate(i)*dtime
!      rho_et_GND(i)=statev(220+i)+ rho_GND_et_rate(i)*dtime
!      rho_en_GND(i)=statev(250+i)+ rho_GND_en_rate(i)*dtime
      
!      rho_GND_total=rho_GND_total+sqrt(rho_s_GND(i)**2.0
!     &+rho_et_GND(i)**2.0+rho_en_GND(i)**2.0) 
      
!      enddo                                                                       !  i=1,30       1






      
!%%%%%%%%%%%%%%%%%%%%%%    Now use Park's method to find rho_GND   %%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      Nye_vector(1)=Nye_tensor(1,1)
      Nye_vector(2)=Nye_tensor(1,2)
      Nye_vector(3)=Nye_tensor(1,3)
      Nye_vector(4)=Nye_tensor(2,1)
      Nye_vector(5)=Nye_tensor(2,2)
      Nye_vector(6)=Nye_tensor(2,3)
      Nye_vector(7)=Nye_tensor(3,1)
      Nye_vector(8)=Nye_tensor(3,2)
      Nye_vector(9)=Nye_tensor(3,3)
          
      
      do i=1,3
      do k=1,12                              ! only consider 12 slip systems for now
      
      b_park(i,k)=gst(i,k)
      t_park(i,k)=gst(i,k)
      b_park(i,k+12)=gst(i,k)
      t_park(i,k+12)=gtt(i,k)
      b_park(i,k+24)=gst(i,k)
      t_park(i,k+24)=gmt(i,k)
      
      enddo
      enddo
      
      do j=1,36
        A_park(1,j)=b_park(1,j)*t_park(1,j);
        A_park(2,j)=b_park(1,j)*t_park(2,j);
        A_park(3,j)=b_park(1,j)*t_park(3,j);
        A_park(4,j)=b_park(2,j)*t_park(1,j);
        A_park(5,j)=b_park(2,j)*t_park(2,j);
        A_park(6,j)=b_park(2,j)*t_park(3,j);
        A_park(7,j)=b_park(3,j)*t_park(1,j);
        A_park(8,j)=b_park(3,j)*t_park(2,j);
        A_park(9,j)=b_park(3,j)*t_park(3,j);
       enddo
       
       do i=1,9
       do j=1,36
       A_T(j,i)=A_park(i,j)
       enddo
       enddo
       
       call matmul(A_park,A_T,AA_T,9,36,9)
       
       do j=1,9
       do i=1,9
       AAt_inv(i,j)=AA_T(i,j)
       enddo
       enddo
       
       call matinv(AAt_inv,9,9)
       call matmul(A_T,AAt_inv,AtAAt_inv,36,9,9)      
       
       call matmul(AtAAt_inv,Nye_vector,rho36,36,9,1)
         
       
       
! now seperate the rho36 into 36 types GND
      do i=1,30
      rho_s_GND(i)=0.0d0
      rho_et_GND(i)=0.0d0
      rho_en_GND(i)=0.0d0
      enddo

      do i=1,12                    !  so far only consider 12 vslip systems 
      rho_s_GND(i) =  rho36(i)
      rho_et_GND(i)=  rho36(i+12)  
      rho_en_GND(i)=  rho36(i+24)
      enddo



      
      do i=1,30
      rho_GND_total=rho_GND_total+sqrt(rho_s_GND(i)**2.0
     &+rho_et_GND(i)**2.0+rho_en_GND(i)**2.0) 
      enddo
      
      !---------------------------------------------------------------
!      if (noel.eq.100) then
!      write(4001,*) 'rho_3_types_are:'
!      do i=1,30
!      write(4001,*) rho_s_GND(i),rho_et_GND(i),rho_en_GND(i)
!      enddo
      !write(4001,*) '--------------------------------------------'
!      endif
      !---------------------------------------------------------------
                                                                                  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      call screw_edge_cos_sin(screw_cos,screw_sin,edge_cos,edge_sin)
      

      
      ! now define x_alpha_beta
      do i=1,30
      do j=1,30
      if (i.eq.j) x_a_b(i,j)=1.0d0
      if (i.ne.j) x_a_b(i,j)=1.4d0
      enddo
      enddo
      
      ! now calculate rho_Forest and rho_Pass
     
      
      
      do i=1,30                                                                   !   1
      rho_F_GND(i)=0.0d0
      rho_P_GND(i)=0.0d0
      do j=1,30
      
      a1=dabs(gmt(1,i)*gst(1,j)+gmt(2,i)*gst(2,j)+gmt(3,i)*gst(3,j))
      a2=dabs(gmt(1,i)*gtt(1,j)+gmt(2,i)*gtt(2,j)+gmt(3,i)*gtt(3,j))
      a3=dabs(gmt(1,i)*gmt(1,j)+gmt(2,i)*gmt(2,j)+gmt(3,i)*gmt(3,j))
      
      if (a1.gt.1.0) a1=1.0d0
      if (a2.gt.1.0) a2=1.0d0
      if (a3.gt.1.0) a3=1.0d0
      
      a4=sqrt(1.0-a1**2.0)
      a5=sqrt(1.0-a2**2.0)
      a6=sqrt(1.0-a3**2.0)
      

      
      
      rho_F_GND(i)= rho_F_GND(i)+ 
     & x_a_b(i,j)*(rho_s_GND(j)*a1+rho_et_GND(j)*a2+rho_en_GND(j)*a3)
     
      rho_P_GND(i)= rho_P_GND(i)+ 
     & x_a_b(i,j)*(rho_s_GND(j)*a4+rho_et_GND(j)*a5+rho_en_GND(j)*a6)
     
      enddo
      
      tau_Pass_GND(i)=dabs(0.033*17.0d3*bur_hcp*sqrt(abs(rho_P_GND(i))))
      ! 0.04 is c1
      !   Mg shear modulus 17 GPa
      !   1E-10 is burger's vector
      
      tau_Cut_GND(i)=dabs(1.0d-6 *1.602d-19*0.33 /2.0 /(bur_hcp**2.0d0)
     & *sqrt(abs(rho_F_GND(i))))
      ! c2=2.0
      ! c3=0.4  thus is not shown in above equation
      ! activaton energy for Mg alloy AZ63: 1ev =1.602E-19 J
      ! 1.0d-6 is for Mpa
      enddo                                                                       !   2
        
            

      if (nGNDswitch.eq.0) then      
      do isys=1,30
      tau_Pass_GND(isys)=0.0d0
      tau_Cut_GND(isys)=0.0d0
      enddo
      endif
      
        
      !------------next is a back stress based on total GND, it's not used, just ignore it----  
      !tau_r= 1.0d-6 *0.33 * 17.0E9 * 1.0E-10 * sqrt(total_rho_GND)
      
      !   1.0d-6 is for Mpa
      !   assume alpha is 0.33  for now
      !   Mg shear modulus 17 GPa
      !   1E-10 is burger's vector
      
      !if (tau_r.gt.0.0) tau_r =0.0
      !--------------------------------------------------------------------------------------
      	  

      	rho_mobile_tt=0.0d0
      	      
	  do isys=1,30
	  
	  if ((rho_F_GND(isys).gt.0.0).and.
     &	(rho_P_GND(isys).gt.0.0)) then  
	  rho_mobile(isys)=sqrt(rho_F_GND(isys)* rho_P_GND(isys))
	  else
	  rho_mobile(isys)=0.0d0
	  endif 
!	  a_dot_m(isys)=4.0E-11*rho_mobile(isys)*0.0d0 + a_dot_hcp
	  a_dot_m(isys)=a_dot_hcp
	  rho_mobile_tt=rho_mobile_tt + rho_mobile(isys)
	  enddo          
!---------------------end of reading last time step -----------------------------	   
 
	  endif      



	  
!-------------   twin hall petch strengthening   ----------------------------------    
            xksoft=1.0d0       
           
           if (f_twin_max.gt.1.0d-6) then           
              xlen_hcp=xksoft/sqrt(f_twin_max*1.0)
           else
              xlen_hcp=xksoft/sqrt(1.0d-6*1.0)
           endif                                       !if (i_twin_active.ne.0)
           
           
           do isys=19,30
              g_alpha_t(isys)=g_alpha_t(isys)+ xlen_hcp
           enddo   
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
           if(itag_el.eq.1)then
              a_dot=a_dot_hcp
              rate_exp=rate_exp_hcp
           endif

           isign=1                                                                      !   tension

           call ssmatc3(dfgrd1,g_alpha_t,tpk_1d_t,s_tan,dtime,
     &     fp_t,de_gg_2d,
     &     c_alpha_1d,delta_gamma,tpk_1d_tau,g_alpha_tau,tau,xn,xr,
     &     xs_bar,xh0,s_gg_2d,cst,cmt,ctt,noel,hgnd,Curl_Fp,
     &     tau_Cut_GND,tau_Pass_GND,a_dot_m)

           if(ifirst_check.eq.1)then
              pnewdt=0.5d0
              return
           endif

           call dstress_fp_fe(dfgrd1,fp_t,s_tan,delta_gamma,tpk_1d_tau,
     &       gst,gmt,fp_tau,fe_tau,stress,dsrot,dmrot,ntens,e_plastic,
     &       cst_1d,tr,cauchy_2d,fe_tau_inv_t,w_p)      

!        Calculating Stresses normal to the slip plane
   
!***************************************************************
   
           if(itag_el.eq.1)then
             do k1=1,1
               do i11=1,3
                 cbasal(i11)=0.d0
                 do j11=1,3
                   cbasal(i11)=cbasal(i11)+fe_tau_inv_t(i11,j11)*
     &             gmt(j11,k1)
                 enddo
               enddo      
       
               do i11=1,3
                 bas_tract(i11)=0.d0
                 do j11=1,3
                   bas_tract(i11)=bas_tract(i11)+cauchy_2d(i11,j11)*
     &             cbasal(j11)
                 enddo
               enddo
      
               bn_stress(k1)=0.d0
               do i11=1,3
                bn_stress(k1)=bn_stress(k1)+bas_tract(i11)*cbasal(i11)
               enddo
             enddo
       
           endif

!******************************************************************
!       Resolved shear stress
!           statev(217)=tau(1)

!      Normal stress
!           statev(218)=bn_stress(1)

!******************************************************************
            
           statev(299)=w_p    !Masoud_crack

	   call rudcmp(fe_tau,re_tau,ue_tau)
	   call trans(re_tau,re_tau_t)
	   call mat33(tlg,re_tau,tr,3)
	   call stereograph(tlg(3,1),tlg(3,2),tlg(3,3),xpole,ypole)
	
	   call jacobian(dfgrd0,dfgrd1,fp_t,fp_tau,s_tan,c_alpha_1d,
     &     tau,delta_gamma,dtime,rate_exp,a_dot,de_gg,tpk_1d_tau,
     &     g_alpha_tau,ddsdde,ntens,tau_Cut_GND,tau_Pass_GND,a_dot_m)




	   if(itag_el.eq.1)index1=0
            index1=0
 !--------------------------------------------------------------------------------------------           
 !        if(kinc.gt.1.or.kstep.gt.1)then
         
           do isys=19,30
           g_alpha_tau(isys)=g_alpha_tau(isys) - xlen_hcp
           enddo
           
  !       endif
!-----------------------------------------------------------------------------------------------            
            
            
	   do isys=1,islip
	      index1=index1+1
	      statev(index1)=g_alpha_tau(isys)                                                !     1:60   g and tot_gamma
	   enddo

	   do isys=1,islip
	      index1=index1+1
	      statev(index1)=tot_gamma(isys)+dabs(delta_gamma(isys))
!-------------------------------------------------------------------------                 ! Jiahao- 151-180: delta_gamma
!	      statev(150+isys)=dabs(delta_gamma(isys))
!-------------------------------------------------------------------------	      
	   enddo
	
	   do i=1,6
	      index1=index1+1
	      statev(index1)=tpk_1d_tau(i)                                                      61-66    1st PK stress
	   enddo
	   do i=1,3
	      do j=1,3
		 index1=index1+1
		 statev(index1)=fp_tau(i,j)                                                          67-75      Fp stress
	      enddo
	   enddo
	   
	   do i=1,6
	      creep(i)=creep(i)+wt(itag_el)*cst_1d(i)
	      stress1(i)=stress1(i)+wt(itag_el)*stress(i)
	      do j=1,6
		 ddsdde1(i,j)=ddsdde1(i,j)+wt(itag_el)*ddsdde(i,j)
	      enddo
	   enddo
c           if(itag_el.eq.1)statev(218)=bn_stress(1)
	   
	   ikin_init_state=0
         ikin_state=75
	   do isys=1,islip
	      statev(ikin_init_state+ikin_state+isys)=chi_tau(isys)                         !   76-105   back stress
	   enddo


!************************Creep*****************************
!	   do i=1,6
!	      statev(210+i)=stress(i)

!	   enddo
!**********************************************************

!	   call calc_fpdot(delta_gamma,fp_tau,s_tan,islip,dtime,nslip,itag_el)


!	   if(pnewdt1.lt.pnewdt)pnewdt=pnewdt1
	   pnewdt=pnewdt1

	enddo


	do i=1,6
	   stress(i)=stress1(i)
	   do j=1,6
	      ddsdde(i,j)=ddsdde1(i,j)
	   enddo
	enddo

	do i=1,6
	   statev(210+i)=stress(i)                                                          211-216    cauchy stress 
	   
	enddo
      
      
      do i=1,30
      statev(180+i)=tau_Cut_GND(i)                                                    !   181-210    Tau_CUT_GND
      statev(220+i)=tau_Pass_GND(i)                                                   !   221-250    Tau_PASS_GND
      enddo
      
      
!      if (noel.eq.100) then
!       write(4001,*) 'tot_Gamma is'
!       do i=1,3
!       write(4001,*)  tot_gamma(i)
!       enddo  
!       endif
       
            
      do i=1,6
      f_twin(i)=f_twin(i)+delta_gamma(i+12)/0.1089
      statev(280+i)=f_twin(i)                                                          !    281-286   f_twin_volume_fraction
      enddo
      


      statev(287)=i_twin_active                                                        ! 287   itwin_active
      statev(288)=rho_GND_total                                                        ! 288  total_GND
	 statev(289)=twin_nu	                                                           !  289 twin_nu
      
      
      do i=1,12
      rho_GND_tt(i)=sqrt(rho_s_GND(i)**2.0
     &+rho_et_GND(i)**2.0+rho_en_GND(i)**2.0)                                          ! 251-262   rho_GND_tt
      statev(250+i)=rho_GND_tt(i)  
      enddo
      
      
      index1=0
      do i=1,3
      do j=1,3
      index1=index1+1
      statev(289+index1)=Curl_Fp(i,j)                                               !290-298     curl_fp
      enddo
      enddo

      do i=1,12
      statev(110+i)=rho_mobile(i)
	enddo
	statev(123)=rho_mobile_tt
!	 statev(301)=bur_hcp
!      statev(302)=sample_size
!      statev(303)=nGNDswitch
	
	
	statev(124)=S_critical
	statev(125)=gsize
	
	
	      
	return

	end
                
                
                
                
                 
                
                
                
                
                
                
                
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
C  CALCULATES THE INVERSE OF A 3*3 MATRIX

	subroutine matinv3(a,ai,det)
c	implicit double precision*8(a-h,o-z)
       implicit double precision (a-h,o-z)

c      include 'aba_param.inc'

      dimension a(3,3), ai(3,3)
c
      det=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(2,1)*a(1,2)
     &     *a(3,3)+a(2,1)*a(1,3)*a(3,2)+a(3,1)*a(1,2)*a(2,3)-a(3,1)
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


C-------------------------------------------------------------------------

C    THIS SUBROUTINE INITIALISES THE CRYSTAL SLIP SYSTEMS, 
C    CALCULATES THE TRANSFORMATION MATRIX AND CONTAINS 
C    STATE VARIABLES

c----------------------------------------------------------
c     cmt - slip plane in crystal system
c     cst - slip directions in crystal system
c     tlg0 - transformation matrix (kalidindis paper)
c------------------------------------------------------------------------


      subroutine inicrys(icrys,ics,cst,cmt,nslips)
     
      implicit real*8 (a-h,o-z)
c     include 'aba_param.inc'    

	parameter(nslip=30)
	parameter(itotslips=8)
	common/slips/islip
	dimension cmt(3,nslip),cst(3,nslip),ics(itotslips),
     &  nslips(itotslips)
  	
	common/crys_switch/ihcp
	common/poly_inicrys/ipoly_switch
	r2=1.d0/2.d0
	r3=dsqrt(3.d0)
	ic=0


	if(ihcp.eq.1)then
	   nslips(1)=3
	   nslips(2)=3
	   nslips(3)=6
	   nslips(4)=6
	   nslips(5)=3
	   nslips(6)=3
	   nslips(7)=3
	   nslips(8)=3
	endif

		

c---hcp slip system 1
	   cmt(1,1)=  0.0d0
	   cmt(2,1)=  0.0d0
	   cmt(3,1)=  1.0d0

	   cst(1,1)=  0.866025403784439
	   cst(2,1)=  0.5d0
	   cst(3,1)=  0.0d0

c---hcp slip system 2
	   cmt(1,2)=  0.0d0
	   cmt(2,2)=  0.0d0
	   cmt(3,2)=  1.0d0
	   
	   cst(1,2)= -0.866025403784439
	   cst(2,2)= 0.5d0
	   cst(3,2)= 0.0d0

c---hcp slip system 3
	   cmt(1,3)=  0.0d0
	   cmt(2,3)=  0.0d0
	   cmt(3,3)=  1.0d0
	  
	   cst(1,3)=  0.0d0
	   cst(2,3)=  -1.0d0
	   cst(3,3)=  0.0d0

!----3 prismatic <a> systems

c---hcp slip system 4
	   cmt(1,4)=  1.0d0
	   cmt(2,4)=  0.0d0
	   cmt(3,4)=  0.0d0
	  
	   cst(1,4)=  0.0d0
	   cst(2,4)=  1.0d0
	   cst(3,4)=  0.0d0

c---hcp slip system 5
	   cmt(1,5)=  0.5d0
	   cmt(2,5)=  0.866025403784438
	   cmt(3,5)=  0.0d0
	  
	   cst(1,5)=  -0.866025403784438
	   cst(2,5)=  0.5d0
	   cst(3,5)=  0.0d0

c---hcp slip system 6
	   cmt(1,6)=  -0.5d0
	   cmt(2,6)=  0.866025403784438
	   cmt(3,6)=  0.0d0
	  
	   cst(1,6)=  -0.866025403784438
	   cst(2,6)=  -0.5d0
	   cst(3,6)=  0.0d0

! ----6 2nd order pyramidal <c+a> systems

c---hcp slip system 7


	   cmt(1,7)=  0.737376788497700
	   cmt(2,7)=  0.425724687333328
	   cmt(3,7)=  0.524436804943893
	  
	   cst(1,7)=  -0.454175595760956
	   cst(2,7)=  -0.262218402471946
	   cst(3,7)=  0.851449374666657

c---hcp slip system 8
	   cmt(1,8)=  0.0d0
	   cmt(2,8)=  0.851449374666657
	   cmt(3,8)=  0.524436804943893
	  
	   cst(1,8)=  0.0d0
	   cst(2,8)=  -0.524436804943893
	   cst(3,8)=  0.851449374666657

c---hcp slip system 9
	   cmt(1,9)=  -0.737376788497700
	   cmt(2,9)=  0.425724687333328 
	   cmt(3,9)=  0.524436804943893
	  
	   cst(1,9)=  0.454175595760956
	   cst(2,9)=  -0.262218402471946
	   cst(3,9)=  0.851449374666657

c---hcp slip system 10
	   cmt(1,10)=  -0.737376788497700
	   cmt(2,10)=  -0.425724687333328
	   cmt(3,10)=  0.524436804943893
	  
	   cst(1,10)=  0.454175595760956
	   cst(2,10)=  0.262218402471946
	   cst(3,10)=  0.851449374666657

c---hcp slip system 11
	   cmt(1,11)=  0.0d0
	   cmt(2,11)=  -0.851449374666657
	   cmt(3,11)=  0.524436804943893
	  
	   cst(1,11)=  0.0d0
	   cst(2,11)=  0.524436804943893 
	   cst(3,11)=  0.851449374666657

c---hcp slip system 12

	   cmt(1,12)=  0.737376788497700
	   cmt(2,12)=  -0.425724687333328 
	   cmt(3,12)=  0.524436804943893
	  
	   cst(1,12)=  -0.454175595760956
	   cst(2,12)=  0.262218402471946
	   cst(3,12)=  0.851449374666657

c---hcp tension twin system 1, total 13
         

         cmt(1,13)=-0.6840d0
         cmt(2,13)=0.0d0
         cmt(3,13)=0.7295d0 
	          
         cst(1,13)=0.7295d0
         cst(2,13)=0.0d0
         cst(3,13)=0.6840d0
         
c---hcp tension twin system 1, total 14

         cmt(1,14)= -0.3420d0
         cmt(2,14)= -0.5923d0
         cmt(3,14)= 0.7295d0
        
         cst(1,14)= 0.3647d0
         cst(2,14)= 0.6318d0
         cst(3,14)= 0.6840d0
         
c---hcp tension twin system 1, total 15
  

         cmt(1,15)= 0.3420d0
         cmt(2,15)= -0.5923d0
         cmt(3,15)= 0.7295d0
         
         cst(1,15)= -0.3647d0 
         cst(2,15)=  0.6318d0
         cst(3,15)=  0.6840d0
         
c---hcp tension twin system 1, total 16


         cmt(1,16)= 0.6840d0
         cmt(2,16)= 0.0d0
         cmt(3,16)= 0.7295d0 
         
         cst(1,16)= -0.7295d0
         cst(2,16)= 0.0d0
         cst(3,16)= 0.6840d0
         
c---hcp tension twin system 1, total 17

         cmt(1,17)= 0.3420d0
         cmt(2,17)= 0.5923d0
         cmt(3,17)= 0.7295d0
         
         cst(1,17)= -0.3647d0
         cst(2,17)= -0.6318d0
         cst(3,17)=  0.6840d0      


c---hcp tension twin system 1, total 18

         cmt(1,18)= -0.3420d0
         cmt(2,18)= 0.5923d0
         cmt(3,18)= 0.7295d0  
         
         cst(1,18)= 0.3647d0
         cst(2,18)= -0.6318d0
         cst(3,18)= 0.6840d0





       do iicc= 19,30	 
	 cmt(1,iicc)=  0.0d0
	 cmt(2,iicc)=  0.0d0 
	 cmt(3,iicc)=  0.0d0
	 
	 cst(1,iicc)=  0.0d0
	 cst(2,iicc)=  0.0d0 
	 cst(3,iicc)=  0.0d0
      enddo


!	if(ihcp.eq.1)islip=30
      islip=30

	return
	end



C-----------------------------------------------------

C    THIS SUBROUTINE CALCULATES THE TRANSFORMATION MATRIX
c----------------------------------------------------------
c     phi   - euler(1)
c     theta - euler(2)
c     omega - euler(3)
c---------------------------------------------------
      subroutine euler_slip(euler,tlgt)
      implicit double precision (a-h,o-z)  

c      include 'aba_param.inc'

      dimension euler(3),tlg(3,3),tlgt(3,3)

      pi=4.d0*datan(1.d0)

      phi=euler(1)
      theta =euler(2)
      omega  =euler(3)
  
      sp=dsin(phi)                      
      cp=dcos(phi)                     
      st=dsin(theta)                     
      ct=dcos(theta)                    
      so=dsin(omega)                    
      co=dcos(omega)   
      tlg(1,1)=co*cp-so*sp*ct
      tlg(1,2)=co*sp+so*ct*cp   
      tlg(1,3)=so*st   
      tlg(2,1)=-so*cp-sp*co*ct 
      tlg(2,2)=-so*sp+ct*co*cp
      tlg(2,3)=co*st
      tlg(3,1)=sp*st       
      tlg(3,2)=-st*cp       
      tlg(3,3)=ct

	call trans(tlg,tlgt)

      return
      end   

C---------------------------------------------------------

C   REFERING TO THE THESIS THIS SUBROUTINE CALCULATES A,T*tr,B,C
C   AND IT CALLS THE UPDATE IMPLICIT SUBROUTINE WHICH USES 
C   NEWTON RAPHSON FOR CALCULATING T* 
c-----------------------------------------------------------------
c
c         fp_t_inv - inverse of theplastic deformation gradient
c         a_tan    - A (from thesis Pg 175)
c         t_tr_1d  - T*tr
c         b_alpha  - B (Pg 175)
c         c_alpha_1d - C (Pg 175)
c----------------------------------------------------------------

      subroutine ssmatc3(dfgrd1,g_alpha_t,tpk_1d_t,s_tan,dtime,
     &fp_t,de_gg_2d,
     &c_alpha_1d,delta_gamma,tpk_1d_tau,g_alpha_tau,tau,xn,xr,
     &xs_bar,xh0,s_gg_2d,gst,gmt,gtt,noel,hgnd,Curl_Fp,
     &tau_Cut_GND,tau_Pass_GND,a_dot_m)

       implicit real*8(a-h,o-z)
      
c      include 'aba_param.inc'

       parameter (nslip=30)
       common/slips/islip
       dimension g_alpha_t(nslip),g_alpha_tau(nslip),
     &      delta_gamma(nslip),hgnd(30,30),Curl_Fp(3,3)

       dimension tpk_1d_t(6),tpk_1d_tau(6),s_tan(3,3,nslip),
     &      fp_t(3,3),fp_t_inv(3,3),fp_t_inv_t(3,3),
     &      dfgrd1(3,3),dfgrd1_t(3,3),xs_bar(nslip),
     &      cc(3,3),a1_tan(3,3),a_tan(3,3),a_tan_1d(6),
     &      b_alpha(3,3,nslip),b_alpha_1d(6,nslip),c_alpha_1d(6,nslip),
     &      t_tr_1d(6),ckrone_1d(6),de_gg_2d(6,6),tau(nslip),
     &      xh0(nslip),xn(nslip),xr(nslip),s_gg_2d(6,6)
      dimension gst(3,nslip),gmt(3,nslip),gtt(3,nslip)
      dimension a_dot_m(30)
      dimension tau_Cut_GND(30), tau_Pass_GND(30)
!      common/GND_part_m/a_dot_m 
      common/crys_const/a_dot,rate_exp
       

      call matinv3(fp_t,fp_t_inv,det_fp_t) 

      if(det_fp_t.eq.0d0) then
        write(*,*) '0 divide'
        stop
      endif    

      call trans(fp_t_inv,fp_t_inv_t)
      call trans(dfgrd1,dfgrd1_t)
      call mat33(cc,dfgrd1_t,dfgrd1,3)
      call mat33(a1_tan,fp_t_inv_t,cc,3)
      call mat33(a_tan,a1_tan,fp_t_inv,3)
      call tr2to1(a_tan,a_tan_1d)
      call cal_ckrone_1d(ckrone_1d)

      do i=1,6
        t_tr_1d(i)=0d0
        do j=1,6
          t_tr_1d(i)=t_tr_1d(i)+0.5d0*de_gg_2d(i,j)
     &   *(a_tan_1d(j)-ckrone_1d(j))
        enddo
      enddo


      do isys=1,islip
        do i=1,3
          do j=1,3
            b_alpha(i,j,isys)=0d0
            do k=1,3
              b_alpha(i,j,isys)=b_alpha(i,j,isys)
     &        +a_tan(i,k)*s_tan(k,j,isys)+s_tan(k,i,isys)*a_tan(k,j)
            enddo
          enddo
        enddo
      enddo


      call tr2to1crys(b_alpha,b_alpha_1d)



      do isys=1,islip
        do i=1,6
          c_alpha_1d(i,isys)=0d0
          do j=1,6
            c_alpha_1d(i,isys)=c_alpha_1d(i,isys)
     &      +0.5d0*de_gg_2d(i,j)*b_alpha_1d(j,isys)
          enddo
        enddo
      enddo

      call update_implicit(dtime,s_tan,delta_gamma,t_tr_1d,c_alpha_1d,
     &g_alpha_t,tpk_1d_t,tpk_1d_tau,g_alpha_tau,tau,xn,xr,xs_bar,xh0,
     &s_gg_2d,gst,gmt,gtt,noel,de_gg_2d,dfgrd1,fp_t,
     &hgnd,Curl_Fp,tau_Cut_GND,tau_Pass_GND,a_dot_m)

      return
      end




c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine sigmat(x,xx)    
c                                                                     
c  arranges into a symmetric array xx, the six commponents of a vector x
c                                                                     
       implicit double precision (a-h,o-z)
c      include 'aba_param.inc'
      dimension x(6),xx(3,3)
c
      xx(1,1)=x(1)                                                    
      xx(2,2)=x(2)                                                    
      xx(3,3)=x(3)                                                    
      xx(2,3)=x(6)                                                    
      xx(3,2)=x(6)                                                    
      xx(3,1)=x(5)                                                    
      xx(1,3)=x(5)                                                    
      xx(1,2)=x(4)                                                    
      xx(2,1)=x(4)                                                    
      return
      end




C---------------------------------------------------------

C   THIS USES NEWTON RAPHSON TO CALCULATE T*
c ---------------------------------------------
c      xjn_2d - Jacobian of Newton Raphson for T*
c      g_n    - Gn (Residual) (Pg 175)
c      xjn_gn - Jn_inv*[Gn] (Pg 175)
c--------------------------------------------

      subroutine update_implicit(dtime,s_tan,delta_gamma,t_tr_1d,
     &c_alpha_1d,
     &g_alpha_t,tpk_1d_t,tpk_1d_tau,g_alpha_tau,tau,xn,xr,xs_bar,xh0,
     &s_gg_2d,gst,gmt,gtt,noel,de_gg_2d,dfgrd1,fp_t,
     &hgnd,Curl_Fp,tau_Cut_GND,tau_Pass_GND,a_dot_m)
      
       implicit real*8(a-h,o-z)
c      include 'aba_param.inc'
     
      parameter (nslip=30,g0=8.0d0)
      common/slips/islip

      dimension tpk_1d_tau(6),tpk_1d_t(6),g_alpha_tau(nslip),
     &xs_bar(nslip),g_alpha_t(nslip),s_tan(3,3,nslip),t_tr_1d(6),
     &g_n(6),c_alpha_1d(6,nslip),delta_gamma(nslip),tau(nslip),
     &xh0(nslip),
     &xjn_2d(6,6),xjn_2d_inv(6,6),xjn_gn(6),xhard(nslip,nslip),
     &xn(nslip),xr(nslip),chi_t(nslip),chi_tau(nslip),s_gg_2d(6,6)
      dimension gst(3,nslip),gmt(3,nslip),gtt(3,nslip)
      dimension a_dot_m(30), hgnd(30,30),Curl_Fp(3,3)
      dimension tau_Cut_GND(30), tau_Pass_GND(30)
      dimension de_gg_2d(6,6),dfgrd1(3,3),fp_t(3,3)      
      real(8):: f_twin(6), f_twin_max
      real(8):: screw_cos(12,12),screw_sin(12,12)
      real(8):: edge_cos(12,12),edge_sin(12,12)
      

      common/twin_kalidandi/f_twin, f_twin_max
!      common/GND_part_m/a_dot_m 
      common/crys_const/a_dot,rate_exp
      common/iter_check/ifirst_check
      common/time_step/pnewdt1

      common/kinematic/xc,xd,chi_t,chi_tau
      common/iatg_check/iatg

      iatg=0


      do i=1,6
       tpk_1d_tau(i)=tpk_1d_t(i)
      enddo

      do isys=1,islip
        g_alpha_tau(isys)=g_alpha_t(isys)
        chi_tau(isys)=chi_t(isys)
      enddo

	g_tol=1.d0
	xtol0=1.0e-10
	chi_tol=1.d0
	
	do i=1,6
	   if(dabs(tpk_1d_tau(i)).gt.tpk_1d_max)then
	      tpk_1d_max=dabs(tpk_1d_tau(i))
	   endif
	enddo
        
	xtol=xtol0*tpk_1d_max
	if(xtol.lt.xtol0)xtol=xtol0

	g_max=2.0*g_tol
	chi_max=2.0*chi_tol
      
      iter1=0
      iter_max=200

C     rnorm, g_max defined to eneter the iteration loops                                           
!******************************************************************************
!   this is level two iteration, to update resistence and back stress       ***
!******************************************************************************

      call screw_edge_cos_sin(screw_cos,screw_sin,
     &  edge_cos,edge_sin)
     
      do while((g_max.gt.g_tol.or.chi_max.gt.chi_tol).and.(iter1.lt.           ! g_max is the max of (g_alpha_tau(isys) - g_prev(isys) )      
     &iter_max))                                                                ! that is the increment of resistence     
                                                                                !  this is level two iteration, to update resistence and back stress
         iter=0
         rnorm=2*xtol
         gamma_tol=0.02d0

         do while(rnorm.gt.xtol.and.iter.le.iter_max)                            !  rnorm is the norm of g_n (residual)
                                                                                  !  This is level one iteration, to update stress
C     Calculating the residual

            xeta=1.d0                                                             !  coefficient of stress increment

            call stress_residue(tpk_1d_tau,s_tan,a_dot,dtime,rate_exp,
     & g_alpha_tau,t_tr_1d,c_alpha_1d,tau,delta_gamma,g_n,xjn_2d,1,
     & s_gg_2d,de_gg_2d,dfgrd1,fp_t,noel,tau_Cut_GND,tau_Pass_GND,
     & a_dot_m)

	    if(iatg.eq.1)then
	       ifirst_check=1
	       write(*,*)'tau/g exceed 2'                                                !  stress increment is to large
	       return
	    endif


            do isys=1,islip
               if(dabs(delta_gamma(isys)).gt.gamma_tol) xeta=0.25d0               ! xeta: thesis 176 A.14     
            enddo

            
            do i=1,6
               do j=1,6 
                  xjn_2d_inv(i,j)=xjn_2d(i,j)
               enddo
            enddo
            
            call matinv(xjn_2d_inv,6,6)
            do i=1,6   
               xjn_gn(i)=0d0
               do k=1,6
                  xjn_gn(i)=xjn_gn(i)+xjn_2d_inv(i,k)*g_n(k)
               enddo
            enddo
            
            do i7=1,6
               tpk_1d_tau(i7)=tpk_1d_tau(i7)-xeta*xjn_gn(i7)     
            enddo   
            
            call stress_residue(tpk_1d_tau,s_tan,a_dot,dtime,rate_exp,             ! this subroutine is to calculate stress residual and jacobine      
     & g_alpha_tau,t_tr_1d,c_alpha_1d,tau,delta_gamma,g_n,xjn_2d,0, 
     & s_gg_2d,de_gg_2d,dfgrd1,fp_t,noel,tau_Cut_GND,tau_Pass_GND,
     & a_dot_m)
           

	    if(iatg.eq.1)then
	       ifirst_check=1
	       write(*,*)'tau/g exceed 2'
	       return
	    endif


            call magn_vect(g_n,rnorm,6)                                           ! g_n is the residual
            
            iter=iter+1
         enddo                                                                  !   end of level one iteration         
         
!----------------------end of level one iteration -----------------------------
	 ifirst_check=0
         if(rnorm.gt.xtol)then                                                   !  residual still too large after iteration limits
            write(6,*)'The First Level Iteration did not converge'
	    write(6,*)'rnorm=',rnorm
            ifirst_check=1
            return
         endif
         
    !----------------added Jiahao 04/26/2013------------------------     
!         iatg=0
!         call stress_residue(tpk_1d_tau,s_tan,a_dot,dtime,rate_exp,     
!     & g_alpha_tau,t_tr_1d,c_alpha_1d,tau,delta_gamma,g_n,xjn_2d,2, 
!     & s_gg_2d,de_gg_2d,dfgrd1,fp_t,noel,tau_Cut_GND,tau_Pass_GND,
!      & a_dot_m)
!         if(iatg.eq.1)then
!	       ifirst_check=1
!	       write(*,*)'tau/g exceed 2'
!	       return
!	    endif
    !----------------------------------------------------------------     
      
      
C     Hardness Iteration                                                            
      call make_hard(xh0,g_alpha_tau,xs_bar,delta_gamma,
     &   dtime,a_dot,xr,xn,xhard,tau,gst,gmt,gtt,chi_tau,a_dot_m)
         g_max=0.d0

 

!------------------- Modified by Jiahao  ------------------  	 
!------------------SLip and twin in parent------------------------------

!      call hard_GND(nslip,Curl_Fp,gmt,g_alpha_tau,noel,nelx,   
!     &      islip,hgnd,ktag)


         do isys1=1,12
            g_prev=g_alpha_tau(isys1)
            g_alpha_tau(isys1)=g_alpha_t(isys1)
            do isys2=1,12  
     
           g_alpha_tau(isys1)=g_alpha_tau(isys1)+                                   !  important, here is the increase of resistence
     &     xhard(isys1,isys2)*dabs(delta_gamma(isys2)* !0.5)
     & 1.0*(0.0*edge_cos(isys1,isys2)+1.0*screw_cos(isys1,isys2)))
 
        
!       	 g_alpha_tau(isys1)=g_alpha_tau(isys1)+      !Masoud-kh
!     &     dabs(hgnd(isys1,isys2))*dabs(delta_gamma(isys2))
    
    
            enddo
            if(dabs(g_alpha_tau(isys1)-g_prev).gt.g_max)g_max=
     &      dabs(g_alpha_tau(isys1)-g_prev)
         enddo
         do isys1=13,18
            g_prev=g_alpha_tau(isys1)
            g_alpha_tau(isys1)=g_alpha_t(isys1)
            do isys2=13,18  
               g_alpha_tau(isys1)=g_alpha_tau(isys1)+                                   !  important, here is the increase of resistence
     &     xhard(isys1,isys2)*dabs(delta_gamma(isys2)*1.0)
       
           enddo
      enddo
 
!-------------------------Below slip in twin -----------------------------	 

!-------------------------------------------------------------------------
!       Back Stress Iteration

           chi_max=0.d0
           do isys1=1,12
              chi_prev=chi_tau(isys1)
    
!	    chi_tau(isys)=chi_t(isys)+xc*delta_gamma(isys)
!	    chi_tau(isys)=chi_tau(isys)/(1.d0+xd*
!     &      dabs(delta_gamma(isys)))
        
             chi_tau(isys1)=chi_t(isys1)
        
             do isys2=1,12
        
               chi_tau(isys1)=chi_tau(isys1)+  
     &         xhard(isys1,isys2)*dabs(delta_gamma(isys2)* !0.5)     
     & 1.0*(0.0*edge_sin(isys1,isys2)+1.0*screw_sin(isys1,isys2)))                     !   back stress is set to tau_pass_increment
      
             enddo
             if(dabs(chi_tau(isys1)-chi_prev).gt.chi_max)chi_max=
     &        dabs(chi_tau(isys1)-chi_prev)
           enddo

          iter1=iter1+1
        enddo
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!**************************************************************************
!***      	        end of level two iteration                          ***
!**************************************************************************
       call stress_residue(tpk_1d_tau,s_tan,a_dot,dtime,rate_exp,
     &     g_alpha_tau,t_tr_1d,c_alpha_1d,tau,delta_gamma,g_n,xjn_2d,0,
     & s_gg_2d,de_gg_2d,dfgrd1,fp_t,noel,tau_Cut_GND,tau_Pass_GND,
     & a_dot_m)

           
      
	if(iatg.eq.1)then
	   ifirst_check=1
	   write(*,*)'tau/g exceed 2'
	   return
	endif

	gamma_tol=0.02d0
	ni=0
        do isys=1,islip
	   if(dabs(delta_gamma(isys)).gt.gamma_tol)pnewdt1=0.5d0
	   if(dabs(delta_gamma(isys)).le.0.01d0*gamma_tol)ni=ni+1
	enddo
	if(ni.eq.islip)pnewdt1=2.0d0 
	
	if(g_max.gt.g_tol)then
	   write(6,*)'Second Level Iteration did not converge'
	   write(6,*)'g_max=',g_max
	   ifirst_check=1
	   return 
	endif


	if(chi_max.gt.chi_tol)then
	   write(6,*)'Back Stress iteration did not converge'
	   write(6,*)'chi_max=',chi_max
	   ifirst_check=1
	   return
	endif
	   
	
	
	return
	end  

C-----------------------------------------------------------------------


C-----------------------------------------------------------------------

      subroutine make_hard(xh0,g_alpha_tau,xs_bar,delta_gamma,
     &  dtime,a_dot,xr,xn,xhard,tau,gst,gmt,gtt,chi_tau,a_dot_m)
      implicit real*8(a-h,o-z)
c       include 'aba_param.inc'

      parameter(nslip=30)
      parameter (itotslips=8)
      common/slips/islip 
      dimension g_alpha_tau(nslip),delta_gamma(nslip),
     &  xhard(nslip,nslip)
      dimension qab(nslip,nslip),gamma_dot(nslip),xh0(nslip),
     & xs_bar(nslip),
     & xn(nslip),xr(nslip),tot_gamma(nslip),xh_0(5),xh_inf(5),
     & xg_0(5),xg_inf(5),tau(nslip),xn_used(nslip),xr_used(nslip),
     & xh0_used(nslip),xs_bar_used(nslip),xn_trb_comp(itotslips),
     & xr_trb_comp(itotslips),xh0_trb_comp(itotslips),
     & xs_trb_comp(itotslips),
     & xn_pra_comp(itotslips),xr_pra_comp(itotslips),
     & xh0_pra_comp(itotslips),
     & xs_pra_comp(itotslips)
      dimension gst(3,nslip),gmt(3,nslip),gtt(3,nslip)

      
      dimension a_dot_m(30),chi_tau(nslip)
      
      common/slip1/tot_gamma
!      common/GND_part_m/a_dot_m 
      common/crys5/xh_0, xh_inf, xg_0,
     &	     xg_inf, itag_el
	common/crys7/a_dot_hcp,rate_exp_hcp
	common/poly_inicrys/ipoly_switch
	common/write_flags/ikinc,ngauss	
	common/tencomp/isign


	do isys_a=1,islip
	   do isys_b=1,islip
	      xhard(isys_a,isys_b)=0.d0
	   enddo
	enddo

	if(dabs(dtime).gt.1d-10)then                                       	           !                                             0000000 
	do isys_a=1,islip                                                               ! q_alpha_beta
	   do isys_b=1,islip
	      if(isys_a.eq.isys_b)then
		 qab(isys_a,isys_b)=1.0d0
	      else
		 qab(isys_a,isys_b)=1.0d0
	      endif
	   enddo
	enddo

!****************************************************************************
!**-------------------HCP--------------------------------------------------**
!****************************************************************************

	
	if(itag_el.eq.1)then                                                             ! hcp phase, need                             11111
	   a_dot=a_dot_hcp
	   do isys=1,islip
	      gamma_dot(isys)=delta_gamma(isys)/dtime
	   enddo
       
	   do index=1,islip                                                             !    index
	      xs_bar_used(index)=xs_bar(index)                                           ! saturation stress 
	      xh0_used(index)=xh0(index)
	      xn_used(index)=xn(index)
	      xr_used(index)=xr(index)
	      
	      
!**************************************************************************************
	      
	   enddo                                                                        !                    end of index
	   

!--------------------------------------------------------------------
!----------what really matters is below------------------------------	   
	   do isys_a=1,islip                                          !
	      do isys_b=1,islip                                       !
		 if(dabs(gamma_dot(isys_b)).gt.0)then    
		 
!		 if (dabs(gamma_dot(isys_b)).le.(a_dot_m(isys_b))) then
                     !  increase of saturation stress 
!		     xs=xs_bar_used(isys_b) *(dabs(gamma_dot(isys_b)/                          !  increase of saturation stress
!     &               a_dot_m(isys_b)))**xn_used(isys_b)                         
!	  	 else
	      	      xs=xs_bar_used(isys_b)
!	  	 endif
	    	 
		 xtemp=1.d0-(dabs(g_alpha_tau(isys_b))+
     &  dabs(chi_tau(isys_b)))/xs                       !  
		    xhard(isys_a,isys_b)=qab(isys_a,isys_b)*                !                     !  increase of h 
     &              xh0_used(isys_b)*                               !   
     &	    ((dabs(xtemp))**xr_used(isys_b))*dsignf(xtemp)          !
     
		 else                                                      !
		    xxxx=1                
		    xhard(isys_a,isys_b)=0.0                                  !   
		 endif                                                     !  
	      enddo                                                  !
	   enddo                                                     !
!---------------------------------------------------------------------	
!---------------------twin-------------------------------------------


	   do isys_a=13,18                                          !
	      do isys_b=13,18                                       !
		 if(dabs(gamma_dot(isys_b)).gt.0)then    
		 
!		 if (dabs(gamma_dot(isys_b)).le.(a_dot_m(isys_b))) then
                     !  increase of saturation stress 
!		     xs=xs_bar_used(isys_b) *(dabs(gamma_dot(isys_b)/                          !  increase of saturation stress
!     &               a_dot_m(isys_b)))**xn_used(isys_b)                         
!	  	 else
	      	      xs=xs_bar_used(isys_b)
!	  	 endif
	    	 
		 xtemp=1.d0-(dabs(g_alpha_tau(isys_b))+
     &  dabs(chi_tau(isys_b)))/xs                       !  
		    xhard(isys_a,isys_b)=qab(isys_a,isys_b)*                !                     !  increase of h 
     &              xh0_used(isys_b)*                               !   
     &	    ((dabs(xtemp))**1.1)*dsignf(xtemp)          !
     
		 else                                                      !
		    xxxx=1                
		    xhard(isys_a,isys_b)=0.0                                  !   
		 endif                                                     !  
	      enddo                                                  !
	   enddo  


!----------------------------------------------------------------------   
	endif                                                                            !itag_el=1         (HCP phase)                              11111


!**********************************************************************
!**-------------------------end of HCP ------------------------------**	
!**********************************************************************
	endif                                                                            ! if(dabs(dtime).gt.1d-10)then                  0000000
	return
	end






      subroutine rclear (a,max)                                       
c                                                                     
c  fills a real array a with zeros
c                                                                     
       implicit double precision (a-h,o-z)                             
c      include 'aba_param.inc'
      dimension a(max)
c
      do i=1,max                                                 
      a(i) = 0.0                                                     
      enddo
      return
      end         
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine prod(a,b,c)                                          
c                                                                     
c  computes the matrix product c=a*b  all matricies are 3x3
c                                                                     
       implicit double precision (a-h,o-z)                             
c      include 'aba_param.inc'
      dimension a(3,3),b(3,3),c(3,3)                                  
c                                                                     
      do 200 j=1,3                                                    
      do 200 i=1,3                                                    
        s=0.0                                                           
      do 100 k=1,3                                                    
        s=s+a(i,k)*b(k,j)                                             
 100  continue
        c(i,j)=s                                                        
 200  continue
      return
      end   


C---------------------------------------------------------------------
C      THIS SUBROUTINE CALCULATES THE VALUE OF Fp,NORMALIZES IT AND
C      CALCULATES THE CAUCHY STRESS AND THE NEW TEXTURE
c-----------------------------------------------------------------
c      gamma_s_2d - delta_gamma*So (used in calculation of Fp_tau)
c      Fe_tau     - elastic part of the deformation gradient
c      dsrot      - new texture
c      dmrot      - new texture
c----------------------------------------------------------------

      subroutine dstress_fp_fe(dfgrd1,fp_t,s_tan,delta_gamma,tpk_1d,
     &gst,gmt,fp_tau,fe_tau,stress,dsrot,dmrot,ntens,e_plastic,cst_1d,
     &tr,cauchy_2d,fe_tau_inv_t,w_p)

      implicit real*8 (a-h,o-z)   
c       include 'aba_param.inc'
      parameter (nslip=30)
      common/slips/islip
      
        dimension dfgrd1(3,3),fp_t(3,3),s_tan(3,3,nslip),
     &   fp_tau(3,3),gst(3,nslip),gmt(3,nslip),delta_gamma(nslip),
     &   tpk_1d(6),gamma_s_2d(3,3),gamma_s_fp(3,3),fe_tau_inv_t(3,3),
     &   fp_tau_inv(3,3),fe_tau(3,3),fe_tau_t(3,3),fe_tau_inv(3,3),
     &   fe_tpk_2d(3,3),fe_tpk_2d_fe_t(3,3),cauchy_2d(3,3),
     &   dmrot(3,nslip),dsrot(3,nslip),tpk_2d_g(3,3),stress(ntens),
     &   up_tau(3,3),rp_tau(3,3),cstrain(3,3),cst_1d(6),ckrone_2d(3,3),
     &   fp_tau_t(3,3),fp_t_fp(3,3),tr(3,3),fe_tr(3,3),vect(3),
     &   vect_0001(3),vect_0001n(3)
     
      
      real(8):: f_twin(6), f_twin_max
      common/twin_kalidandi/f_twin, f_twin_max 
      
      do i=1,3
         do j=1,3
            up_tau(i,j)=0.0d0
            rp_tau(i,j)=0.0d0
         enddo
      enddo



      do i=1,3
        do j=1,3
          gamma_s_2d(i,j)=0d0
  
          do isys=1,12
            gamma_s_2d(i,j)=gamma_s_2d(i,j)
!     &    +(1.0-f_twin_max)*delta_gamma(isys)*s_tan(i,j,isys)
     &    + delta_gamma(isys)*s_tan(i,j,isys) 
          enddo
          do isys=13,18
            gamma_s_2d(i,j)=gamma_s_2d(i,j)
     &    +delta_gamma(isys)*s_tan(i,j,isys)
          enddo
          do isys=19,islip
            gamma_s_2d(i,j)=gamma_s_2d(i,j)
     &    +f_twin_max*delta_gamma(isys)*s_tan(i,j,isys)
          enddo

!---------------------------          
          
        enddo 
      enddo

  
      call mat33(gamma_s_fp,gamma_s_2d,fp_t,3)


      
      do i=1,3
        do j=1,3
          fp_tau(i,j)=fp_t(i,j)+gamma_s_fp(i,j)
        enddo
      enddo

      call matinv3(fp_tau,fp_tau_inv,det_fp_tau)
      if(det_fp_tau.eq.0d0) then
        write(*,*) '0 divide'        
        stop
      endif


      if(det_fp_tau.ne.1.) then
        do i=1,3
          do j=1,3
            a=-1./3.
            fp_tau(i,j)=det_fp_tau**a*fp_tau(i,j)
          enddo
        enddo
            fp_tau(1,1) = 1
            fp_tau(2,2) = 1
            fp_tau(3,3) = 1
         call matinv3(fp_tau,fp_tau_inv,det_fp_tau)
      endif
     
     
      call sigmat(tpk_1d,tpk_2d_g)
      call mat33(fe_tau,dfgrd1,fp_tau_inv,3)
      call trans(fe_tau,fe_tau_t)      
      call matinv3(fe_tau,fe_tau_inv,det_fe_tau)
      call trans(fe_tau_inv,fe_tau_inv_t)

      if(det_fe_tau.eq.0d0) then
        write(*,*) '0 divide'
        stop
      endif

      call mat33(fe_tpk_2d,fe_tau,tpk_2d_g,3)
      call mat33(fe_tpk_2d_fe_t,fe_tpk_2d,fe_tau_t,3)

      do i=1,3
        do j=1,3 
          cauchy_2d(i,j)=(1/det_fe_tau)*fe_tpk_2d_fe_t(i,j)
!         cauchy_2d(i,j)=tpk_2d_g(i,j)         
          
        enddo
      enddo   

      stress(1)=cauchy_2d(1,1)
      stress(2)=cauchy_2d(2,2)
      stress(3)=cauchy_2d(3,3)
      stress(4)=0.5d0*(cauchy_2d(1,2)+cauchy_2d(2,1))
      stress(5)=0.5d0*(cauchy_2d(1,3)+cauchy_2d(3,1))
      stress(6)=0.5d0*(cauchy_2d(2,3)+cauchy_2d(3,2))


      do i=1,3
        do isys=1,islip
          dmrot(i,isys)=0d0
          dsrot(i,isys)=0d0
          do j=1,3
            dsrot(i,isys)=dsrot(i,isys)+fe_tau(i,j)*gst(j,isys)
            dmrot(i,isys)=dmrot(i,isys)+fe_tau_inv_t(i,j)*gmt(j,isys)
          enddo
        enddo
      enddo


c---------------creep strain-------------------------
	call trans(fp_tau,fp_tau_t)
	call mat33(fp_t_fp,fp_tau_t,fp_tau,3)
	call cal_ckrone_2d(ckrone_2d)
	do i=1,3
	   do j=1,3
	      cstrain(i,j)=0.5*(fp_t_fp(i,j)-ckrone_2d(i,j))
	   enddo
	enddo
	cst_1d(1)=cstrain(1,1)
	cst_1d(2)=cstrain(2,2)
	cst_1d(3)=cstrain(3,3)
	cst_1d(4)=cstrain(1,2)
	cst_1d(5)=cstrain(1,3)
	cst_1d(6)=cstrain(2,3)
	
c----------------------------------------------------
c--------Plastic work-------------------------------
c$$$        w_p=0.d0
c$$$        do i=1,3
c$$$           do j=1,3
c$$$              w_p=w_p+cauchy_2d(i,j)*cstrain(i,j)
c$$$           enddo
c$$$        enddo
c$$$        w_p=0.5*w_p
        do i=1,3
           do j=1,3
              w_p=w_p+0.5d0*cauchy_2d(i,j)*
     $                (gamma_s_2d(i,j)+gamma_s_2d(j,i))
           enddo
        enddo    !Masoud_crack(End)
c---------------------------------------------------
C-------------Final Pole Figure---------------


       call mat33(fe_tr,fe_tau,tr,3)
       vect(1)=0.0d0
       vect(2)=0.0d0
       vect(3)=1.0d0
       do i=1,3
           vect_0001(i)=0.0d0
           do j=1,3
              vect_0001(i)=vect_0001(i)+fe_tr(i,j)*vect(j)
           enddo
       enddo

       xnormfac=sqrt((vect_0001(1))**2+(vect_0001(2))**
     &  2+(vect_0001(3))**2)
       do i=1,3
          vect_0001n(i)=vect_0001(i)/xnormfac
       enddo

       call stereograph(vect_0001n(1),vect_0001n(2),
     &  vect_0001n(3),xpole1,ypole1)
      

C----------------------------------------------
      
      return
      end  



c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
C    MULTIPLICATION OF TWO 3*3 MATRICES

      subroutine mat33(x,y,z,l)
c       
       implicit double precision(a-h,o-z)
c       include 'aba_param.inc'
       dimension x(3,l),y(3,3),z(3,l)
c
        do i=1,3
          do j=1,l
             x(i,j)=0d0
             do  k=1,3
                x(i,j)=x(i,j)+y(i,k)*z(k,j)
             enddo
          enddo
        enddo
        return
        end
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--               
      subroutine rotate(tch,r,t)                                      
c                                                                     
c  computes the matrix t=r*tch*(rtrans)  (all matrices are 3x3)
c                                                                     
       implicit double precision (a-h,o-z)
c      include 'aba_param.inc'
       dimension t(3,3),r(3,3),tch(3,3)                                
c                                                                     
      do 200 j=1,3                                                    
      do 200 i=1,3                                                    
      s=0.0                                                           
      do 100 l=1,3                                                    
        rjl=r(j,l)                                                    
        do 100 k=1,3                                                  
        s=s+r(i,k)*rjl*tch(k,l)                                       
 100  continue                                                      
      t(i,j)=s                                                        
 200  continue                                                        
      return
      end 


c------------------------------------------------------------------------

c    CALCULATES THE JACOBIAN FOR UMAT
C-------------------------------------------------------------------------
C      REFERING TO THESIS Pg 185
c      cl_4d - Lijkl
c      de_g  - Cijkl
c      d_tang_4d -Dijkl
c      g_tan_4d -Gmnkl
c      cj_tan_4d - Jijkl
c      b_tan_2d  - B
c      ck_tan_4d - Kijkl
c      q_tan     -Qijkl
c      r_tan     -Rij
c      Sijkl_tan - Sijkl
c      w_tan     - Wijkl 
C------------------------------------------------------------------------
      subroutine jacobian(dfgrd0,dfgrd1,fp_t,fp_tau,s_tan,c_alpha_1d,
     &tau,delta_gamma,dtime,rate_exp,a_dot,de_g,tpk_1d_tau,g_alpha_tau,
     &ddsdde,ntens,tau_Cut_GND,tau_Pass_GND,a_dot_m)
      
      
      implicit real*8(a-h,o-z)
c     include 'aba_param.inc'
      
      parameter (nslip=30)
      common/slips/islip
      dimension fe_t_t(3,3),fe_tau_inv(3,3),
     &     cl1(3,3),cl2(3,3),fe_gamma_s(3,3),de_g(3,3,3,3),
     &     cl_4d(3,3,3,3),d_tang1_4d(3,3,3,3),defe_1(nslip),
     &     d_tang_4d(3,3,3,3),g_tan_4d(3,3,3,3,nslip),tpk_2d_g(3,3),
     &     cj_tan_4d(3,3,3,3,nslip),cj1_tan_4d(3,3,3,3,nslip), 
     &     b_tan_2d(3,3,nslip),ck_tan_4d(3,3,3,3),ck_tan_2d(6,6),
     &     ckrone(3,3,3,3),cb_tan(3,3,3,3),c_alpha(3,3,nslip),
     &     gamma_j(3,3,3,3),q_tan(3,3,3,3),r_tan(3,3,nslip),
     &     ck_tan_2d_inv(6,6),ck_tan_4d_inv(3,3,3,3),q1_tan(3,3,3,3),
     &     sijkl_tan_1st(3,3,3,3),sijkl_tan_2nd(3,3,3,3),r_s(3,3,3,3),
     &   t_fe_t(3,3),w_tan_1st(3,3,3,3),fe_q(3,3,3,3),dfgrd0_inv(3,3),
     &     w_tan_2nd(3,3,3,3),w_tan(3,3,3,3),w_tan_3rd(3,3,3,3),
     &     w_tan_4th(3,3,3,3),fe_stress_2d(3,3),fe_stress_2d_fe_t(3,3),
     &   sijkl_tan_3rd(3,3,3,3),sijkl_tan(3,3,3,3),g_alpha_tau(nslip),
     &     sijkl_tan_fe_inv(3,3),rt_u_fe_t(3,3),dfgrd0(3,3),
     &     w_tan_2d(6,6),ddsdde(ntens,ntens),tau(nslip),dfgrd1(3,3)

      dimension q2_tan(3,3,3,3) 
      dimension c_alpha_1d(6,nslip),fe_t(3,3),fe_tau(3,3),
     &     fe_tau_t(3,3),
     &     s_tan(3,3,nslip),
     &     fp_t(3,3),fp_t_inv(3,3),
     &     delta_gamma(islip),fp_tau(3,3),fp_tau_inv(3,3),
     &     gamma_s_2d(3,3),ckrone_2d(3,3),
     &     tpk_1d_tau(6),
     &     dfgrdt(3,3),rt_tau(3,3),ut_tau(3,3)    

	dimension xh_0(5),xh_inf(5),xg_0(5),xg_inf(5),tot_gamma(nslip)

	dimension a_dot_bcc(5),rate_exp_bcc(5)
      dimension a_dot_m(30)
      dimension tau_Cut_GND(30), tau_Pass_GND(30)
      
      real(8):: f_twin(6), f_twin_max
      common/twin_kalidandi/f_twin, f_twin_max
      common/slip1/tot_gamma
!      common/GND_part_m/a_dot_m 
      common/crys5/xh_0, xh_inf, xg_0,
     &     xg_inf, itag_el
      common/crys6/a_dot_bcc,rate_exp_bcc
      common/crys7/a_dot_hcp,rate_exp_hcp

      isign=1



      call cal_ckrone_2d(ckrone_2d)


      call matinv3(dfgrd0,dfgrd0_inv,det)      


      call mat33(dfgrdt,dfgrd1,dfgrd0_inv,3)

      call rudcmp(dfgrdt,rt_tau,ut_tau)


      call matinv3(fp_t,fp_t_inv,det)




      call mat33(fe_t,dfgrd0,fp_t_inv,3)
      call mat33(cl1,ut_tau,fe_t,3)
      call trans(fe_t,fe_t_t)
      call mat33(cl2,fe_t_t,ut_tau,3)
      call sigmat(tpk_1d_tau,tpk_2d_g)


      call matinv3(fp_tau,fp_tau_inv,det)



      call mat33(fe_tau,dfgrd1,fp_tau_inv,3)
      call trans(fe_tau,fe_tau_t)
      call matinv3(fe_tau,fe_tau_inv,det_fe_tau)





      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              cl_4d(i,j,k,l)=fe_t_t(i,k)*cl1(l,j)+cl2(i,k)*fe_t(l,j)
            enddo
          enddo
        enddo
      enddo

    

      call mat44(d_tang1_4d,de_g,cl_4d) 




      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              d_tang_4d(i,j,k,l)=0.5d0*d_tang1_4d(i,j,k,l)
            enddo
          enddo
        enddo
      enddo
      

      do isys=1,islip
        do m=1,3
          do n=1,3
            do k=1,3
              do l=1,3
                g_tan_4d(m,n,k,l,isys)=0d0
                do ip=1,3
                  g_tan_4d(m,n,k,l,isys)=g_tan_4d(m,n,k,l,isys)
     &    +cl_4d(m,ip,k,l)*s_tan(ip,n,isys)
     &    +s_tan(ip,m,isys)*cl_4d(ip,n,k,l)
                enddo
              enddo
            enddo
          enddo
        enddo 
      enddo



      call mat44crys(cj1_tan_4d,de_g,g_tan_4d)
      do isys=1,islip
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                cj_tan_4d(i,j,k,l,isys)=0.5d0
     & *cj1_tan_4d(i,j,k,l,isys)
              enddo
            enddo
          enddo
        enddo
      enddo
 
 

 
 

      do isys=1,islip
        g_inv=1d0/(dabs(g_alpha_tau(isys))+dabs(tau_Cut_GND(isys)))
        rate_inv=1d0/rate_exp
        sc=tau(isys)/(g_alpha_tau(isys)+tau_Cut_GND(isys))
        atg=dabs(sc)

        

	if(itag_el.eq.1)then
	   a_dot=a_dot_hcp
	   rate_exp=rate_exp_hcp
	endif

	if(itag_el.eq.2)then
	   if (isys.le.12)then
	      a_dot=a_dot_bcc(1)
	      rate_exp=rate_exp_bcc(1)
	   endif
	   if(isys.ge.13.and.isys.le.48)then
	      if(isign.eq.1)then
		 if(isys.le.24)then
		    a_dot=a_dot_bcc(2)
		    rate_exp=rate_exp_bcc(2)
		 else
		    a_dot=a_dot_bcc(4)
		    rate_exp=rate_exp_bcc(4)
		 endif
	      else
		 if(isys.le.24)then
		    a_dot=a_dot_bcc(3)
		    rate_exp=rate_exp_bcc(3)
		 else
		    a_dot=a_dot_bcc(5)
		    rate_exp=rate_exp_bcc(5)
		 endif
	      endif
	   endif
	    endif




        defe_1(isys)=a_dot_m(isys)*dtime*g_inv*atg**(rate_inv-1)
     &  *rate_inv

        do i=1,3
          do j=1,3
            b_tan_2d(i,j,isys)=0.5d0*defe_1(isys)
     &        *(s_tan(i,j,isys)+s_tan(j,i,isys))
          enddo
        enddo
      enddo




      call sigmat_crys(c_alpha_1d,c_alpha,islip)  
      call mat33crys(cb_tan,c_alpha,b_tan_2d)
      call cal_ckrone(ckrone)
c
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              ck_tan_4d(i,j,k,l)=ckrone(i,j,k,l)+cb_tan(i,j,k,l)
            enddo
          enddo
        enddo
      enddo          

c--- (b) qijkl=
      call tr4to2(ck_tan_4d,ck_tan_2d)

      do i=1,6
        do j=1,6 
          ck_tan_2d_inv(i,j)=ck_tan_2d(i,j)
        enddo
      enddo

      call matinv(ck_tan_2d_inv,6,6)
      call tr2to4(ck_tan_2d_inv,ck_tan_4d_inv)

      do m=1,3
        do n=1,3
          do k=1,3
            do l=1,3
              gamma_j(m,n,k,l)=0d0
                            
              do isys=1,12 
                gamma_j(m,n,k,l)=gamma_j(m,n,k,l)
!     & +(1.0-f_twin_max)*delta_gamma(isys)*cj_tan_4d(m,n,k,l,isys)
     & + delta_gamma(isys)*cj_tan_4d(m,n,k,l,isys)
              enddo
              do isys=13,18 
                gamma_j(m,n,k,l)=gamma_j(m,n,k,l)
     & + delta_gamma(isys)*cj_tan_4d(m,n,k,l,isys)
              enddo
              do isys=19,islip 
                gamma_j(m,n,k,l)=gamma_j(m,n,k,l)
     & +f_twin_max*delta_gamma(isys)*cj_tan_4d(m,n,k,l,isys)
              enddo
!--------------------------------
              
            enddo
          enddo
        enddo
      enddo

      call mat44(q1_tan,ck_tan_4d_inv,d_tang_4d)
      call mat44(q2_tan,ck_tan_4d_inv,gamma_j)
      
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              q_tan(i,j,k,l)=q1_tan(i,j,k,l)-q2_tan(i,j,k,l)
            enddo
          enddo
        enddo
      enddo


      do isys=1,islip
        do i=1,3
          do j=1,3
            r_tan(i,j,isys)=0d0
            do k=1,3 
              do l=1,3
                r_tan(i,j,isys)=r_tan(i,j,isys)
     &    +b_tan_2d(k,l,isys)*q_tan(k,l,i,j)      
              enddo
            enddo
          enddo
        enddo
      enddo 


      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              sijkl_tan_1st(i,j,k,l)=rt_tau(i,k)*fe_t(l,j)
            enddo
          enddo
        enddo
      enddo
      
      
      do i=1,3
         do j=1,3
            gamma_s_2d(i,j)=0.d0
            
      do isys=1,12
               gamma_s_2d(i,j)=gamma_s_2d(i,j)+
!     & (1.0-f_twin_max)*delta_gamma(isys)*s_tan(i,j,isys)
     &  delta_gamma(isys)*s_tan(i,j,isys) 
      enddo
      
       do isys=13,18         
               gamma_s_2d(i,j)=gamma_s_2d(i,j)+
     &  delta_gamma(isys)*s_tan(i,j,isys)
      enddo
      
      do isys=19,30
               gamma_s_2d(i,j)=gamma_s_2d(i,j)+
     &  f_twin_max*delta_gamma(isys)*s_tan(i,j,isys)
      enddo
!-------------------------            
            
         enddo
      enddo

      call mat33(fe_gamma_s,fe_t,gamma_s_2d,3)


      do i=1,3 
        do j=1,3
          do k=1,3
            do l=1,3
              sijkl_tan_2nd(i,j,k,l)=rt_tau(i,k)*fe_gamma_s(l,j)
            enddo
          enddo
        enddo
      enddo

      call mat33crys(r_s,r_tan,s_tan)
      call mat33(rt_u_fe_t,dfgrdt,fe_t,3)

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              sijkl_tan_3rd(i,j,k,l)=0d0
              do ip=1,3
                sijkl_tan_3rd(i,j,k,l)=sijkl_tan_3rd(i,j,k,l)
     &    +rt_u_fe_t(i,ip)*r_s(k,l,ip,j) 
              enddo
            enddo
          enddo
        enddo
      enddo

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              sijkl_tan(i,j,k,l)=sijkl_tan_1st(i,j,k,l)
     &  -sijkl_tan_2nd(i,j,k,l)-sijkl_tan_3rd(i,j,k,l)
            enddo
          enddo
        enddo
      enddo

      call mat33(t_fe_t,tpk_2d_g,fe_tau_t,3) 
      call mat42(w_tan_1st,sijkl_tan,t_fe_t)


      do i=1,3
        do n=1,3
          do k=1,3
            do l=1,3
              fe_q(i,n,k,l)=0d0
              do m=1,3
                fe_q(i,n,k,l)=fe_q(i,n,k,l)+fe_tau(i,m)
     &    *q_tan(m,n,k,l)
              enddo
            enddo
          enddo
        enddo
      enddo

      call mat42(w_tan_2nd,fe_q,fe_tau_t)


      call mat33(fe_stress_2d,fe_tau,tpk_2d_g,3)

      

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              w_tan_3rd(i,j,k,l)=0d0
              do n=1,3
                w_tan_3rd(i,j,k,l)=w_tan_3rd(i,j,k,l)
     &  +fe_stress_2d(i,n)*sijkl_tan(j,n,k,l)
              enddo
            enddo
          enddo
        enddo
      enddo


      call mat33(fe_stress_2d_fe_t,fe_stress_2d,fe_tau_t,3)
      do k=1,3
        do l=1,3
          sijkl_tan_fe_inv(k,l)=0d0 
          do ip=1,3
            do iq=1,3
              sijkl_tan_fe_inv(k,l)=sijkl_tan_fe_inv(k,l)
     &  +sijkl_tan(ip,iq,k,l)*fe_tau_inv(iq,ip)   
            enddo
          enddo
        enddo
      enddo

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              w_tan_4th(i,j,k,l)=fe_stress_2d_fe_t(i,j)
     &  *sijkl_tan_fe_inv(k,l)
            enddo
          enddo
        enddo
      enddo

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              w_tan(i,j,k,l)=1/det_fe_tau*(w_tan_1st(i,j,k,l)
     &  +w_tan_2nd(i,j,k,l)
     &  +w_tan_3rd(i,j,k,l)-w_tan_4th(i,j,k,l))
            enddo
          enddo
        enddo
      enddo

      call tr4to2(w_tan,w_tan_2d)

       do i=1,6
	 do j=1,3
	    ddsdde(i,j)=w_tan_2d(i,j)
	 enddo
         do j=4,6
            ddsdde(i,j)=w_tan_2d(i,j)/2
         enddo
      enddo



    
      return
      end



c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine mat44crys(x,y,z)
c
       implicit double precision(a-h,o-z)
c      include 'aba_param.inc'
      parameter (nslip=30)
      common/slips/islip	
      dimension x(3,3,3,3,nslip),y(3,3,3,3),z(3,3,3,3,nslip) 
c
      do isys=1,islip
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                x(i,j,k,l,isys)=0d0
                do m=1,3
                  do n=1,3
                    x(i,j,k,l,isys)=x(i,j,k,l,isys)
     &                +y(i,j,m,n)*z(m,n,k,l,isys)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
      return
      end
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
C    CONVERTS FOURTH ORDER TENSOR TO EQUIVALENT SECOND ORDER TENSOR

      subroutine tr4to2(a_4d,a_2d)
      implicit real*8(a-h,o-z)
c     include 'aba_param.inc'
      dimension a_4d(3,3,3,3),a_2d(6,6)

      call rclear66(a_2d)

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

      return
      end


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
C     CONVERTS A SECOND ORDER TENSOR TO A EQUIVALENT FOURTH ORDER TENSOR

      subroutine tr2to4(a_2d,a_4d)
      implicit real*8(a-h,o-z)
c     include 'aba_param.inc'
      parameter (nslip=30)
      common/slips/islip	
      dimension a_4d(3,3,3,3),a_2d(6,6)

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a_4d(i,j,k,l)=0d0
               enddo
            enddo
         enddo
      enddo
c
c---1-3
      do i=1,3
         a_4d(i,i,1,1)=a_2d(i,1)
         a_4d(i,i,2,2)=a_2d(i,2)
         a_4d(i,i,3,3)=a_2d(i,3)
         a_4d(i,i,1,2)=a_2d(i,4)
         a_4d(i,i,1,3)=a_2d(i,5)
         a_4d(i,i,2,3)=a_2d(i,6)
      enddo

      a_4d(1,2,1,1)=a_2d(4,1)
      a_4d(1,2,2,2)=a_2d(4,2)
      a_4d(1,2,3,3)=a_2d(4,3)
      a_4d(1,2,1,2)=a_2d(4,4)
      a_4d(1,2,1,3)=a_2d(4,5)
      a_4d(1,2,2,3)=a_2d(4,6)

      a_4d(1,3,1,1)=a_2d(5,1)
      a_4d(1,3,2,2)=a_2d(5,2)
      a_4d(1,3,3,3)=a_2d(5,3)
      a_4d(1,3,1,2)=a_2d(5,4)
      a_4d(1,3,1,3)=a_2d(5,5)
      a_4d(1,3,2,3)=a_2d(5,6)

      a_4d(2,3,1,1)=a_2d(6,1)
      a_4d(2,3,2,2)=a_2d(6,2)
      a_4d(2,3,3,3)=a_2d(6,3)
      a_4d(2,3,1,2)=a_2d(6,4)
      a_4d(2,3,1,3)=a_2d(6,5)
      a_4d(2,3,2,3)=a_2d(6,6)

      a_4d(2,1,1,1)=a_2d(4,1)
      a_4d(2,1,2,2)=a_2d(4,2)
      a_4d(2,1,3,3)=a_2d(4,3)
      a_4d(2,1,1,2)=a_2d(4,4)
      a_4d(2,1,1,3)=a_2d(4,5)
      a_4d(2,1,2,3)=a_2d(4,6)

      a_4d(3,1,1,1)=a_2d(5,1)
      a_4d(3,1,2,2)=a_2d(5,2)
      a_4d(3,1,3,3)=a_2d(5,3)
      a_4d(3,1,1,2)=a_2d(5,4)
      a_4d(3,1,1,3)=a_2d(5,5)
      a_4d(3,1,2,3)=a_2d(5,6)

      a_4d(3,2,1,1)=a_2d(6,1)
      a_4d(3,2,2,2)=a_2d(6,2)
      a_4d(3,2,3,3)=a_2d(6,3)
      a_4d(3,2,1,2)=a_2d(6,4)
      a_4d(3,2,1,3)=a_2d(6,5)
      a_4d(3,2,2,3)=a_2d(6,6)
      
      return
      end

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--       
      subroutine matinv(a,n,id)                                         
c     ******************************************************************
c     this subroutine computes the inverse and determinant of matrix a *
c     of order n,by the gauss-jordan method, a-inverse replaces a ,and *
c     the determinant of a is placed in determ. if m=1 the vector b    *
c     contans the constant vector when matinv is called,and this is    *
c     replaced with the solution vector if m=0,no simultaneous         *
c     equation solution is called for, and b is not pertinent. n is not*
c     to exceed 100.                                                   *
c      a--is the matrix of coefficients or the matrix to be inverted.  *
c      a contans a-inverse after execution.                            *
c      n-- is the order of the square matrix, a.                       *
c      b--is the matrix contaning column vectors of constants (each    *
c         column vector is associated with a in the following          *
c         manner--ax=b.).                                              *
c      m--is the number of criterion vectors (column vectors of        *
c         simultaneous solutions) to be calculated and stored in b.    *
c      m=0 results in the computation of only the inverse and          *
c          determinant of a.                                           *
c      determ--contans the value of the determinant after execution.   *
c     ******************************************************************
c     
       implicit double precision(a-h,o-z)
c      include 'aba_param.inc'
c
      dimension ipivot(1000),index(1000,2),pivot(1000)                  
      dimension a(id,1),b(500,1)                                     

c     initialization                                                    
c      determ=1.0e0                                                      
c     search for pivot element                                          
      do 30 j=1,n                                                       
   30 ipivot (j)=0                                                      
      do 470 i=1,n                                                      
      amax=0.0e0                                                        
      do 150 j=1,n                                                      
      if(ipivot(j).eq.1) go to 150                                      
      do 140 k=1,n                                                      
      if(ipivot(k).gt.1) go to 590                                      
      if(ipivot(k).eq.1) go to 140                                      
      if(dabs(amax).ge.dabs(a(j,k))) go to 140                          
      irow=j                                                            
      icolum=k                                                          
      amax=a(j,k)                                                       
  140 continue                                                          
  150 continue                                                          
      ipivot(icolum)=ipivot(icolum)+1                                   
c     interchange rows to put pivot element on diagonal                 
      if(irow.eq.icolum) go to 280                                      
c      determ=-determ                                                    
      do 220 l=1,n                                                      
      swap=a(irow,l)                                                    
      a(irow,l)=a(icolum,l)                                             
  220 a(icolum,l)=swap                                                  
c      do 270 l=1,m                                                      
c      swap=b(irow,l)                                                    
c      b(irow,l)=b(icolum,l)                                             
c  270 b(icolum,l)=swap                                                  
  280 index(i,1)=irow                                                   
      index(i,2)=icolum                                                 
      pivot(i)=a(icolum,icolum)                                         
c      determ=determ*pivot(i)                                            
c     divide pivot row by pivot element                                 
      a(icolum,icolum)=1.0e0                                            
      do 340 l=1,n                                                      
  340 a(icolum,l)=a(icolum,l)/pivot(i)                                  
c      if(m.le.0) go to 380                                              
c      do 370 l=1,m                                                      
c  370 b(icolum,l)=b(icolum,l)/pivot(i)                                  
c     reduce non-pivot rows                                             
  380 do 470 l1=1,n                                                     
      if(l1.eq.icolum) go to 470                                        
      t=a(l1,icolum)                                                    
      a(l1,icolum)=0.0e0                                                
      do 430 l=1,n                                                      
  430 a(l1,l)=a(l1,l)-a(icolum,l)*t                                     
c      if(m.le.0) go to 470                                              
c      do 460 l=1,m                                                      
c  460 b(l1,l)=b(l1,l)-b(icolum,l)*t                                     
  470 continue                                                          
c     interchange columns                                               
      do 580 i=1,n                                                      
      l=n+1-i                                                           
      if(index(l,1).eq.index(l,2)) go to 580                            
      jrow=index(l,1)                                                   
      jcolum=index(l,2)                                                 
      do 570 k=1,n                                                      
      swap=a(k,jrow)                                                    
      a(k,jrow)=a(k,jcolum)                                             
      a(k,jcolum)=swap   
c       write(115,*) k,jcolum,a(k,jcolum)                                               
  570 continue                                                          
  580 continue                                                          
  590 return                                                            
      end 

c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine mat44(x,y,z)
c
       implicit double precision(a-h,o-z)
c      include 'aba_param.inc'
      dimension x(3,3,3,3),y(3,3,3,3),z(3,3,3,3) 
c
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              x(i,j,k,l)=0d0
              do m=1,3
                do n=1,3
                  x(i,j,k,l)=x(i,j,k,l)+y(i,j,m,n)*z(m,n,k,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
      return
      end           
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
C    CALCULATES THE TRANSPOSE OF A MATRIX

       subroutine trans(x,x_t)
        implicit double precision(a-h,o-z)
c       include 'aba_param.inc'
       dimension x(3,3),x_t(3,3) 
c
       do i=1,3
         do j=1,3
           x_t(j,i)=x(i,j)
         enddo
       enddo
       return
       end
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine mat33crys(x,y,z)

       implicit double precision (a-h,o-z)
c      include 'aba_param.inc'
      parameter (nslip=30)
      common/slips/islip	
      dimension x(3,3,3,3),y(3,3,islip),z(3,3,nslip)
c
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              x(i,j,k,l)=0d0
              do isys=1,islip
                x(i,j,k,l)=x(i,j,k,l)
     &    +y(i,j,isys)*z(k,l,isys)
              enddo
            enddo
          enddo
        enddo
      enddo
c
      return
      end





      subroutine mat42(x,y,z)
       implicit real*8 (a-h,o-z)
c      include 'aba_param.inc'
      dimension x(3,3,3,3),y(3,3,3,3),z(3,3)
c
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              x(i,j,k,l)=0d0
                do n=1,3
                  x(i,j,k,l)=x(i,j,k,l)
     &          +y(i,n,k,l)*z(n,j)
                enddo
              enddo
            enddo
          enddo
        enddo
c
        return
        end


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-
C     FORMS A FOURTH ORDER IDENTITY TENSOR

      subroutine cal_ckrone(ckrone)

       implicit double precision (a-h,o-z)
c      include 'aba_param.inc'
      dimension ckrone(3,3,3,3)

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              ckrone(i,j,k,l)=0d0
            enddo
           enddo
         enddo
      enddo
c  
      ckrone(1,1,1,1)=1d0
      ckrone(2,2,2,2)=1d0
      ckrone(3,3,3,3)=1d0
      ckrone(1,2,1,2)=1d0
      ckrone(1,3,1,3)=1d0
      ckrone(2,3,2,3)=1d0
      ckrone(2,1,2,1)=1.d0
      ckrone(3,1,3,1)=1.d0
      ckrone(3,2,3,2)=1.d0
c
c
      return
      end
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
C      FORMS A 1-D IDENTITY VECTOR

      subroutine cal_ckrone_1d(ckrone_1d)
       implicit double precision (a-h,o-z)
c      include 'aba_param.inc'    
      dimension ckrone_1d(6)
c
      do i=1,6
        if(i.eq.1.or.i.eq.2.or.i.eq.3) then
          ckrone_1d(i)=1d0
        else
          ckrone_1d(i)=0d0
        endif
      enddo

      return
      end
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
C    FORMS A SECOND ORDER IDENTITY TENSOR

      subroutine cal_ckrone_2d(ckrone_2d)

       implicit double precision (a-h,o-z)
c      include 'aba_param.inc'
       dimension ckrone_2d(3,3)
c
       do i=1,3
          do j=1,3
             if(i.eq.j) then
                ckrone_2d(i,j)=1.0d0
             else
                ckrone_2d(i,j)=0.0d0
             endif
          enddo
       enddo
c     
       return
       end
c---  +----1----+----2----+----3----+----4----+----5----+----6----+----7--
C     SETS ALL THE COMPONENTS OF A 6*6 MATRIX TO ZERO

      subroutine rclear66(x)

       implicit double precision (a-h,o-z)
c      include 'aba_param.inc'
      dimension x(6,6)
c
      do i=1,6
        do j=1,6
          x(i,j)=0d0
        enddo
      enddo
c
      return
      end
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine sigmat_crys(x,xx,islip)    
c                                                                     
c  arranges into a symmetric array xx, the six commponents of a vector x
c                                                                     
       implicit double precision (a-h,o-z)
c      include 'aba_param.inc'
      dimension x(6,islip),xx(3,3,islip)
c
      do isys=1,islip
        xx(1,1,isys)=x(1,isys)
        xx(2,2,isys)=x(2,isys)
        xx(3,3,isys)=x(3,isys)
        xx(2,3,isys)=x(6,isys)
        xx(3,2,isys)=x(6,isys)
        xx(3,1,isys)=x(5,isys)
        xx(1,3,isys)=x(5,isys)
        xx(1,2,isys)=x(4,isys)
        xx(2,1,isys)=x(4,isys)
      enddo

      return
      end

c******************************************************

        
      subroutine rotate_crys(tch,r,t)                                  
c                                                                     
c  computes the matrix t=r*tch*(rtrans)  (all matrices are 3x3)
c                                                                     
      implicit double precision (a-h,o-z)
      parameter(nslip=30)
      common/slips/islip	
c      include 'aba_param.inc'
      dimension t(3,3,nslip),r(3,3),tch(3,3,nslip)                     
      
      do isys=1,islip
        do j=1,3                                                    
          do i=1,3                                                    
            s=0.0                                                      
            do l=1,3                                                    
              rjl=r(j,l)                                               
              do k=1,3                                                  
                s=s+r(i,k)*rjl*tch(k,l,isys)                           
              enddo
            enddo                                                      
            t(i,j,isys)=s
          enddo
        enddo
      enddo      
                                    
      return
      end 
	
c----------------------------------------------------
C     CALCULATES THE SIGN 
             
      double precision function dsignf(x)
      double precision x
      if(x.ge.0)then
         dsignf=+1.d0
      else
         dsignf=-1.d0
      endif
      return
      end
      
c-----------------------------------------------------
      subroutine rotate_crys_vector(a,q,a_prime)
      implicit real*8(a-h,o-z)
      parameter(nslip=30)
      common/slips/islip	
      dimension a(3,nslip),a_prime(3,nslip),q(3,3)
      do isys=1,islip
         do i=1,3
            a_prime(i,isys)=0.d0
            do j=1,3
               a_prime(i,isys)=a_prime(i,isys)+q(i,j)*a(j,isys)
            enddo
         enddo
      enddo
      return
      end


c------------------------------------------------------------
C   TRANSFORMS THE FOURTH ORDER ELASTICITY TENSOR FROM 
C   CRYSTAL SYSTEM TO GLOBAL 

      subroutine rotate4d(ddsdde_4d,tr,de_gg)
      implicit real*8(a-h,o-z)
C     transforms a fourth order tensor
      dimension ddsdde_4d(3,3,3,3),de_gg(3,3,3,3),tr(3,3)
      dimension de_g11(3,3,3,3),de_g22(3,3,3,3),de_g33(3,3,3,3)

      do i=1,3
         do n=1,3
            do io=1,3
               do ip=1,3
                  de_g11(i,n,io,ip)=0d0
                  do m=1,3
                     de_g11(i,n,io,ip)=de_g11(i,n,io,ip)
     &                    +tr(i,m)*ddsdde_4d(m,n,io,ip)
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      do i=1,3
         do j=1,3
            do io=1,3
               do ip=1,3
              de_g22(i,j,io,ip)=0d0
              do n=1,3
                 de_g22(i,j,io,ip)=de_g22(i,j,io,ip)
     &                +tr(j,n)*de_g11(i,n,io,ip)
              enddo
           enddo
          enddo
       enddo
      enddo
      
      do i=1,3
         do j=1,3
            do k=1,3
               do ip=1,3
                  de_g33(i,j,k,ip)=0d0
                  do io=1,3
                     de_g33(i,j,k,ip)=de_g33(i,j,k,ip)
     &                    +tr(k,io)*de_g22(i,j,io,ip)
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      do i=1,3
         do j=1,3
          do k=1,3
             do l=1,3
              de_gg(i,j,k,l)=0d0
              do ip=1,3
                 de_gg(i,j,k,l)=de_gg(i,j,k,l)
     &                +tr(l,ip)*de_g33(i,j,k,ip)
              enddo
           enddo
        enddo
      enddo
      enddo
      
      return
      end


      subroutine tr2to1(a_2d,a_1d)
      implicit real*8(a-h,o-z)
      dimension a_2d(3,3),a_1d(6)
      a_1d(1)=a_2d(1,1)
      a_1d(2)=a_2d(2,2)
      a_1d(3)=a_2d(3,3)
      a_1d(4)=0.5d0*(a_2d(1,2)+a_2d(2,1))
      a_1d(5)=0.5d0*(a_2d(1,3)+a_2d(3,1))
      a_1d(6)=0.5d0*(a_2d(2,3)+a_2d(3,2))
      return
      end

      subroutine tr2to1crys(a_2d,a_1d)
      implicit real*8(a-h,o-z)
      parameter(nslip=30)
      common/slips/islip	
      dimension a_2d(3,3,nslip),a_1d(6,nslip)
      do isys=1,islip
         a_1d(1,isys)=a_2d(1,1,isys)
         a_1d(2,isys)=a_2d(2,2,isys)
         a_1d(3,isys)=a_2d(3,3,isys)
         a_1d(4,isys)=0.5d0*(a_2d(1,2,isys)+a_2d(2,1,isys))
         a_1d(5,isys)=0.5d0*(a_2d(1,3,isys)+a_2d(3,1,isys))
         a_1d(6,isys)=0.5d0*(a_2d(2,3,isys)+a_2d(3,2,isys))
      enddo
      return
      end
c----------------------------------------------------------------
C       CALCULATES THE RESIDUE FOR NEWTON RAPHSON FOR T*
C------------------------------------------------------
C       g_n    - Residual
c       tau    - resolved shear stress
c       xjn_2d - Jacobain for Newton Raphson
c----------------------------------------------------------------
      subroutine stress_residue(tpk_1d_tau,s_tan,a_dot,dtime,rate_exp,
     & g_alpha_tau,t_tr_1d,c_alpha_1d,tau,delta_gamma,g_n,xjn_2d,iflag,
     & s_gg_2d,de_gg_2d,dfgrd1,fp_t,noel,tau_Cut_GND,tau_Pass_GND,
     & a_dot_m)

      implicit real*8(a-h,o-z)
      parameter(nslip=30)
      common/slips/islip
      dimension tpk_1d_tau(6),s_tan(3,3,nslip),g_alpha_tau(nslip)
      dimension c_alpha_1d(6,nslip),tau(nslip),delta_gamma(nslip)
      dimension tpk_2d(3,3),gamma_c_alpha(6),t_tr_1d(6),g_n(6)
      dimension defe_11(nslip),b_alpha_der(3,3,nslip),xjn(3,3,3,3)
      dimension xjn_2d(6,6),c_alpha_2dd(3,3,nslip)
      dimension a_dot_bcc(5),rate_exp_bcc(5),chi_t(nslip),
     &  chi_tau(nslip),s_gg_2d(6,6),e_strain(6,1),e_strain_2d(3,3),
     &     fe_t_fe_tpk(3,3),fe_t_fe(3,3)
      dimension xh_0(5),xh_inf(5),xg_0(5),xg_inf(5),tot_gamma(nslip)
      dimension a_dotj(nslip),rate_expj(nslip)
      dimension a_dot_m(30),f_twin(6)
      dimension tau_Cut_GND(30), tau_Pass_GND(30)
      dimension de_gg_2d(6,6),dfgrd1(3,3),fp_t(3,3)
      dimension fp_tau(3,3),fe_tau(3,3),fp_tau_inv(3,3)
      dimension fe_tau_t(3,3),fe_tau_t_fe_tau(3,3),ckrone_2d(3,3)
      dimension tpk_1d_new(6)
      
      !added by jiaxi for nonlinear constitutive elastic 
      dimension ddsdde(6,6)
     
      common/twin_kalidandi/f_twin, f_twin_max, twin_nu    
!      common/GND_part_m/a_dot_m
      common/slip1/tot_gamma 
      common/crys5/xh_0, xh_inf, xg_0,xg_inf, itag_el
      common/crys6/a_dot_bcc,rate_exp_bcc
      common/crys7/a_dot_hcp,rate_exp_hcp
      common/write_flags/ikinc,ngauss
      common/tencomp/isign
      common/kinematic/xc,xd,chi_t,chi_tau
      common/forjacob1/a_dotj,rate_expj
      common/iatg_check/iatg

      if((iflag.eq.0).or.(iflag.eq.1).or.(iflag.eq.2))then

c--------------------------------------------------------------------------------------------------
c
c       for feTfe
c-------------------------------------------------------------------------------------------------

         call matmul(s_gg_2d,tpk_1d_tau,e_strain,6,6,1)
         call sigmat(e_strain,e_strain_2d)
         do i=1,3
            do j=1,3
               if(i.eq.j)then
                  fe_t_fe(i,j)=1.d0+(2.0*e_strain_2d(i,j))
               else
                  fe_t_fe(i,j)=2.0*e_strain_2d(i,j)
               endif
               enddo
               enddo        
         call sigmat(tpk_1d_tau,tpk_2d)
         call mat33(fe_t_fe_tpk,fe_t_fe,tpk_2d,3)
c
c-------------------------------------------------------------
c
c       call sigmat(tpk_1d_tau,tpk_2d)
      ! modified jiaxi, forcing tau(isys) = 0 for all isys
        do isys=1,islip
           tau(isys)=0.0d0
        enddo
       
        do isys = 1,islip
           delta_gamma(isys) = 0.d0
        enddo     
        
        iflag_elastic = 1        
        if(iflag_elastic == 0)then
          do isys=1,islip
            sc=tau(isys)/(g_alpha_tau(isys)+tau_Cut_GND(isys))
            atg=dabs(sc)
            if(itag_el.eq.1)then
               a_dot=a_dot_hcp
               rate_exp=rate_exp_hcp
            endif
!-------------------------itag =2   bcc-------------------------------------------	    
            if(itag_el.eq.2)then
               if (isys.le.12)then
                  a_dot=a_dot_bcc(1)
                  rate_exp=rate_exp_bcc(1)
               endif
               if(isys.ge.13.and.isys.le.48)then
                  if(isign.eq.1)then
                     if(isys.le.24)then
                       a_dot=a_dot_bcc(2)
                       rate_exp=rate_exp_bcc(2)
                     else 
                       a_dot=a_dot_bcc(4)
                       rate_exp=rate_exp_bcc(4)
                     endif
                  else
                     if(isys.le.24)then
                       a_dot=a_dot_bcc(3)
                       rate_exp=rate_exp_bcc(3)
                     else
                       a_dot=a_dot_bcc(5)
                       rate_exp=rate_exp_bcc(5)
                    endif
                 endif
              endif
           endif
!----------------------------------------------------------------------------------
    
         a_dotj(isys)=a_dot_m(isys)
         rate_expj(isys)=rate_exp

         if (iflag.eq.2) then       
            if ((isys.lt.13).or.(isys.gt.18)) then
               if(atg.gt.2.d0)then
                  iatg=1
                  write(*,*)'tau/g >>2,isys and noel is',isys,noel
                  return
               endif
            endif
         endif
         
         delta_gamma(isys)=a_dot_m(isys)*dtime*atg**(1./rate_exp)*
     &   dsignf(tau(isys))
     
 !--------------------twin-------------------------------------------------------           
         if ((isys.gt.12).and.(isys.lt.19)) then                   
            if (dabs(delta_gamma(isys)).ge.0.01) then
                 delta_gamma(isys)=0.01*dsignf(delta_gamma(isys))
            endif
         endif
 !-------------------------------------------------------------------------------           
         
         enddo
       endif

!----------------------------Jiahao--Twin_kalidindi----------------------------     
            
       f_twin_tot=0.0
       if(iflag_elastic.eq.0)then
          do ii=1,6
          f_twin_tot=f_twin_tot + f_twin(ii)
          enddo
          do ii=1,6
             if ((f_twin(ii)+delta_gamma(ii+12)/0.1089).lt.0.0) then
               delta_gamma(ii+12)=0.0d0
               tau(ii+12)=0.0d0
             endif   
             if (f_twin_tot.gt.1.0) then
               delta_gamma(ii+12)=0.0d0
               tau(ii+12)=0.0d0
             endif

             if (twin_nu.ne.1.0)then
               delta_gamma(ii+12)=0.0d0
               tau(ii+12)=0.0d0
             endif
          enddo

        endif

!------------------------------------------------------------------------------
!************************ Jiahao - I modified the following part to make it consistent with AVE model****

        do i=1,3
           do j=1,3
              tmp11=0.d0     
              do isys=1,12
                do k=1,3    
                 tmp11=tmp11+delta_gamma(isys)*s_tan(i,k,isys)*fp_t(k,j)     
                enddo
              enddo
              fp_tau(i,j)=fp_t(i,j)+tmp11
           enddo
        enddo
 
        call matinv3(fp_tau,fp_tau_inv,det_fp_tau)

!	 do i=1,3
!	    do j=1,3
!	       fp_tau(i,j)=det_fp_tau**(-1.d0/3.d0)*fp_tau(i,j)
!	    enddo
!	 enddo
  
        call matinv3(fp_tau,fp_tau_inv,det_fp_tau)
 
        call mat33(fe_tau,dfgrd1,fp_tau_inv,3)
 
        call trans(fe_tau,fe_tau_t)
 
        call mat33(fe_tau_t_fe_tau,fe_tau_t,fe_tau,3)
 
        call cal_ckrone_2d(ckrone_2d)

        do i=1,3
           do j=1,3
           e_strain_2d(i,j)=0.5d0*(fe_tau_t_fe_tau(i,j)-ckrone_2d(i,j))
           enddo
        enddo
 
 
        call tr2to1(e_strain_2d,e_strain)

        call matmul(de_gg_2d,e_strain,tpk_1d_new,6,6,1) 

         
         
!****************************************************************************

         
        do it=1,6
           gamma_c_alpha(it)=0.0d0
           
           do isys=1,12
              gamma_c_alpha(it)=gamma_c_alpha(it)
!     &  +(1.0-f_twin_max)*delta_gamma(isys)*c_alpha_1d(it,isys)
     &        +delta_gamma(isys)*c_alpha_1d(it,isys)  
           enddo
            
           do isys=13,18
              gamma_c_alpha(it)=gamma_c_alpha(it)
     &        +delta_gamma(isys)*c_alpha_1d(it,isys)
           enddo
           do isys=19,islip
              gamma_c_alpha(it)=gamma_c_alpha(it)
     &        +f_twin_max*delta_gamma(isys)*c_alpha_1d(it,isys)
           enddo
           
!--------------------------------            
        enddo
        
        do i2=1,6
           g_n(i2)=tpk_1d_tau(i2)-t_tr_1d(i2)+gamma_c_alpha(i2)
           
!           g_n(i2)=tpk_1d_tau(i2)-tpk_1d_new(i2)
         
        enddo
         
      endif
!-------------------------------------------------------------------------------------------------     
!************************ Jiahao - I modified the above  part to make it consistent with AVE model****



    
         
C     Calculating the jacobian for nr iteration
	if(iflag.eq.1)then
!----------------------------------------------------------------------	
	   do isys=1,islip
	      
	      if(itag_el.eq.1)then
		 a_dot=a_dot_hcp
		 rate_exp=rate_exp_hcp
	      endif
	      
!------------------------------------------------------	      
	      if(itag_el.eq.2)then
		 if (isys.le.12)then
		    a_dot=a_dot_bcc(1)
		    rate_exp=rate_exp_bcc(1)
		 endif
		 if(isys.ge.13.and.isys.le.48)then
		    
		    if(isign.eq.1)then
		    
		       if(isys.le.24)then
			  a_dot=a_dot_bcc(2)
			  rate_exp=rate_exp_bcc(2)
		       else
			  a_dot=a_dot_bcc(4)
			  rate_exp=rate_exp_bcc(4)
		       endif
		    
		    else
		    
		       if(isys.le.24)then
			  a_dot=a_dot_bcc(3)
			  rate_exp=rate_exp_bcc(3)
		       else
			  a_dot=a_dot_bcc(5)
			  rate_exp=rate_exp_bcc(5)
		       endif
		       
		    endif
		    
		 endif
	      endif
!-------------------------------------------------------	      
	      
	      
	   g_inv=1.d0/(dabs(g_alpha_tau(isys))+dabs(tau_Cut_GND(isys)))
	      rate_inv=1.d0/rate_exp
	       sc=tau(isys)/(g_alpha_tau(isys)+tau_Cut_GND(isys))

	      atg=dabs(sc)
	     
	      

	      defe_11(isys)=a_dot_m(isys)*dtime*g_inv*atg**(rate_inv-1)
     &           *rate_inv
	      do i=1,3
		 do j=1,3
		    b_alpha_der(i,j,isys)=0.5d0*defe_11(isys)
     &                 *(s_tan(i,j,isys)+s_tan(j,i,isys))
		 enddo
	      enddo
	      
	   enddo
!-----------------------------------------------------------------------------	   
	   
	   call sigmat_crys(c_alpha_1d,c_alpha_2dd,islip)
	   
	   do i=1,3
	      do j=1,3
		 do k=1,3
		    do l=1,3
		       xjn(i,j,k,l)=0.0d0

               do isys=1,12
			  xjn(i,j,k,l)=xjn(i,j,k,l)
!     &   + (1.0-f_twin_max)*c_alpha_2dd(i,j,isys)*b_alpha_der(k,l,isys)
     &    + c_alpha_2dd(i,j,isys)*b_alpha_der(k,l,isys)     
		       enddo
		       do isys=13,18
			  xjn(i,j,k,l)=xjn(i,j,k,l)
     &    +c_alpha_2dd(i,j,isys)*b_alpha_der(k,l,isys)
		       enddo
		       do isys=19,islip
			  xjn(i,j,k,l)=xjn(i,j,k,l)
     &   + f_twin_max*c_alpha_2dd(i,j,isys)*b_alpha_der(k,l,isys)
		       enddo
!----------------------------------		       
		       
		       
		    enddo
		 enddo
	      enddo
	   enddo
	   
	   
	   
	   call tr4to2(xjn,xjn_2d)      
	   
	   do i=1,6
	      do j=1,6
		 if(i.eq.j)xjn_2d(i,j)=xjn_2d(i,j)+1.d0
	      enddo
	   enddo
	   
	endif
	
C	write(6,*)'exit sr'
	return
	end



c------------------------------------------
C    CALCULATES THE NORM OF THE RESIDUAL
      subroutine magn_vect(a,r,n)
      implicit real*8(a-h,o-z)
      dimension a(n)
      r=0.d0
      do i=1,n
         r=r+a(i)**2
      enddo
      r=dsqrt(r)
      return
      end
c-------------------------------------------
C    RU DECOMPOSITON          

	subroutine rudcmp(f,rd,ud) 
      implicit real*8 (a-h,o-z)                                        
C      implicit integer*8 (i-n)                                          
C      real *8 li1,li2,li3,lamda1,lamda2,lamda3
c                 
      dimension rd(3,3),f(3,3),ud(3,3) 

c    
C	write(6,*)'entering ru'
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
      ui3=c11*s11+c12*(c2313-c12*c33)
     1      +c13*(c1223-c22*c13)
      ui1s=ui1*ui1
      q    =dsqrt(-dmin1(o3*ui2-ui1s,0.d0))                     
      r    =0.5*(ui3-ui1*ui2)+ui1*ui1s
      xmod =q*q*q
      scl1 =.5d0+dsign(.5d0,xmod-1.d-30)                          
      scl2 =.5d0+dsign(.5d0,xmod-dabs(r))                   
      scl0 =dmin1(scl1,scl2)                               
      scl1 =1.-scl0

      if(scl1.eq.0)then
        xmodscl1=xmod
      else
        xmodscl1=xmod+scl1
      endif

      sdetm=dacos(r/(xmodscl1))*o3

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
C	write(6,*)'exiting ru'
      return
      end


	subroutine stereograph(xl,xm,xn,x,y)
	implicit real*8 (a-h,o-z)
	x=dsignf(xn)*xl/(xn*dsignf(xn)+1.d0)
	y=dsignf(xn)*xm/(xn*dsignf(xn)+1.d0)
	return
	end

	integer function ngrain(nel,m,n)
	implicit real*8 (a-h,o-z)
	n3=(nel-1)/(m*n)**2+1
	nz=(n3-1)/m+1
	nel_2d=nel-(n3-1)*(m*n)**2
	nt1=nel_2d-1
	nt2=m*n
	n1=nt1-(nt1/nt2)*nt2+1
	nx=(n1-1)/m+1
	n2=(nel_2d-n1)/(m*n)+1
	ny=(n2-1)/m+1
	ngrain=n**2*(nz-1)+n*(ny-1)+nx
	return
	end

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
C    MATRIX MULTIPLICATION SUBROUTINE 

      subroutine matmul(aa,bb,cc,n,m,l)
      
c      include 'aba_param.inc'
       implicit double precision(a-h,o-z)

c   this subroutine returns martix [cc] = matrix [aa] * matrix [bb]
c   n=no. of rows of [aa]      =no. of rows of [cc]
c   m=no. of columns of [aa]   =no. of rows of [bb]
c   l=no. of columns of [bb]   =no. of columns of [cc]

       dimension aa(n,m),bb(m,l),cc(n,l)

       do 10 i=1,n
       do 10 j=1,l
         cc(i,j)=0.0
       do 10 k=1,m
   10    cc(i,j)=cc(i,j)+aa(i,k)*bb(k,j)
       return 
       end  

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7--

      subroutine hallpetch(props,nprops1,nprops,grainsize,
     &  ghcp,gbcc,gpyr)
      implicit real*8(A-H,O-Z)
	parameter(intel=1)   
!      include 'fea3d.h' 
      dimension g0_bcc(5),gi_b1(5),gi_b2(5),gi_b3(5)
      dimension ghcp(4),gbcc(1)
      dimension gpyr(2), props(nprops)
      
     

      
      
      pi=4.d0*datan(1.d0)

c     initialize values
         do jj=1,4
            ghcp(jj)=0.d0
         enddo
         do jj=1,2
            gpyr(jj)=0.d0
         enddo
 



      
  
c     material constants
!	    open(unit=303,file='mat_singlet.dat')


c     G_hcp------shear modulus 
c     g0_ab-----hardness of basal
c     g0_ap-----hardness of prism
c     g0_bcc---hardness of bcc
c     gi_ab----internal g of basal
c     gi_ap----internal g of prism
c     gi_b1,gi_b2,gi_b3-----internal g of three families of bcc

c     Material constants and hardness
c     primary alpha
!         call readstr(303)
!         read(303,*)g0_ab_a,g0_ap_a,g0_cpyr,g0_twin_ex,
!     &g0_ab_a2,g0_ap_a2,g0_cpyr2,g0_twin_cm
c     shear modulus and poissons ratio
!         call readstr(303)
!         read(303,*)G_hcp_a,xnu_hcp_a
c     burger vector
!         call readstr(303)
!         read(303,*)bur_hcp
!         call readstr(303)
!         read(303,*)son   ! single or not
!         call readstr(303)
!         read(303,*)nGNDswitch, sample_size
!         close(303)
         
         
      G_hcp_a=props(nprops1+9)
      xnu_hcp_a=props(nprops1+10)   
      bur_hcp=props(nprops1+11)
      son=props(nprops1+12)
      nGNDswitch=props(nprops1+13)
      sample_size=props(nprops1+14)  
      grainsize=props(nprops1+15)                  !     Grain and Lath Size
         
      g0_ab_a=props(nprops1+1)
      g0_ap_a=props(nprops1+2)
      g0_cpyr=props(nprops1+3)
      g0_twin_ex=props(nprops1+4)
      g0_ab_a2=props(nprops1+5)
      g0_ap_a2=props(nprops1+6)
      g0_cpyr2=props(nprops1+7)
      g0_twin_cm=props(nprops1+8)
         
         
         
c     alpha
         alpha_hcp_a=(1.0+(1.0-xnu_hcp_a))*0.5d0
        




            index=0
            index1=0
            index2=0

            do isys=1,1
c     primary alpha hcp
               bur=bur_hcp
               if(1.eq.1)then
                  Gmod=G_hcp_a
                  alpha=alpha_hcp_a
                  xnu=xnu_hcp_a
                  
                  
                     g0_ab=g0_ab_a
                     g0_ap=g0_ap_a
                     xbar=0.24*Gmod
                     
               endif




                     xksoft=(sqrt((2.0-xnu)*pi*Gmod*xbar*
     &                    bur/(2.0*(1.0-xnu)))/(10.d0**6))*(10.d0**3)
                     xkhard=(sqrt(4.0*Gmod*xbar*bur/
     &                    (alpha*pi))/(10.d0**6))*(10.d0**3)  
                        xlen_hcp=xksoft/sqrt(grainsize)


               
               
c     primary alpha
                 
                     gi_ab=g0_ab
                     gi_ap=g0_ap
                  
c     basal & prism                 
             
                     ghcp(1)=g0_ab_a+son*xlen_hcp
                     ghcp(2)=g0_ap_a+son*xlen_hcp
                     ghcp(3)=g0_twin_ex+son*xlen_hcp
                     ghcp(4)=g0_twin_cm+son*xlen_hcp

c     <a> pyr and <c+a> pyramidal
                  gpyr(1)=g0_cpyr+son*xlen_hcp
                  gpyr(2)=g0_cpyr+son*xlen_hcp      
                  
            enddo          
      

      return

      end




!-----------------------------------------------------------------------
      subroutine curlmat(T,C)
      implicit none
    !   C: curl tensor
    !   T: input gradient tensor (3*3*3) 
    
      real(8):: T(3,3,3), C(3,3), alter(3,3,3)
      integer::  i,j,k,r,s
      
      do i=1,3
      do j=1,3
      do k=1,3
      alter(i,j,k)=0.0d0
      enddo
      enddo
      enddo
      
      alter(1,2,3)=1.0d0
      alter(2,3,1)=1.0d0
      alter(3,1,2)=1.0d0
      alter(2,1,3)=-1.0d0
      alter(1,3,2)=-1.0d0
      alter(3,2,1)=-1.0d0
      
      do i=1,3
      do j=1,3
      C(i,j)=0.0d0
      enddo
      enddo
      
      do i=1,3
      do j=1,3
      do s=1,3
      do r=1,3
      C(i,j)=C(i,j)+alter(i,r,s)*T(j,s,r)
      enddo
      enddo
      enddo
      enddo
    
      return
      end
      
      
      
      
      subroutine crossproduct(v1,v2,C)
      
      real(8):: v1(3), v2(3), C(3),alter(3,3,3)
      
      integer:: i,j,k
      
      do i=1,3
      do j=1,3
      do k=1,3
      alter(i,j,k)=0.0d0
      enddo
      enddo
      enddo
      
      alter(1,2,3)=1.0d0
      alter(2,3,1)=1.0d0
      alter(3,1,2)=1.0d0
      alter(2,1,3)=-1.0d0
      alter(1,3,2)=-1.0d0
      alter(3,2,1)=-1.0d0
      
      do i=1,3
      C(i)=0.0d0
      do j=1,3
      do k=1,3
      C(i)=C(i)+alter(i,j,k)*v1(j)*v2(k)
      enddo
      enddo
      enddo
      
      return
      end


      subroutine doubledot(v1,v2,C)
      
      real(8):: v1(3,3), v2(3,3), C
      
      integer:: i,j,k
      
      C=0.0d0
      
      do i=1,3
      do j=1,3
      
      C=C+v1(i,j)*v2(i,j)
      enddo
      enddo
      
      return
      end
      
      
      subroutine tensorproduct(v1,v2,C)
      
      real(8):: v1(3), v2(3), C(3,3)
      
      integer:: i,j,k
      
      do i=1,3
      do j=1,3
      C(i,j)=0.0d0
      enddo
      enddo
      
      do i=1,3
      do j=1,3
      C(i,j)=v1(i)*v2(j)
      enddo
      enddo
      
      return
      end


















!----------------Masoud----------------------------------------------

!      subroutine hard_GND(nslip,cfp,gmt,g_alpha_tau,noel,nelx,    !Masoud-kh
!     &      islip,hgnd,ktag)     
      
!       implicit double precision(a-h,o-z)

!      dimension cfp(3,3),gmt(3,nslip),g_alpha_tau(nslip),
!     &           a(3),hgnd(30,30)
              


        

!      xk0=5.0d0
!      alpha=0.33333333333333d0
!      xmiu=48.0d3   !in Mpa
!      bur=1.0d-4    !in micron

!      do i=1,30         !islip
!        do j=1,30       !islip
!          hgnd(i,j)=0.0d0
!        end do
!      end do          






!      do i=1,3      !islip
!        do j=1,3    !islip

c     calculate lamda_alpha
     
!          do ii=1,3
!            a(ii)=0.0d0
!            do jj=1,3
!              a(ii)=a(ii)+cfp(jj,ii)*gmt(jj,j)
!            end do
!          end do

!          xlamda=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

c                     

 
!          c1=g_alpha_tau(i)-20.0
     
          
!          if (abs(c1).gt.1.0d0) then    !test
!            c2=(0.000026d0*xk0*bur*(alpha*xmiu)**2.0d0)/c1
!          else
!            c1=1.0d0
!            c2=(0.000026d0*xk0*bur*(alpha*xmiu)**2.0d0)/c1
!          end if

!          hgnd(i,j)=c2*xlamda
            
!            hgnd(i,j)=1.0d0
                     

c        if (noel.eq.1.and.ktag.eq.1.and.i.eq.1.and.j.eq.1)    !test
c     &     write (45,*) c1,c2,xlamda,(cfp(1,kkk),kkk=1,3),
c     &       (cfp(2,kkk),kkk=1,3),(cfp(3,kkk),kkk=1,3)
c        if (noel.eq.4707.and.ktag.eq.1.and.i.eq.1.and.j.eq.1)    !test
c     &     write (46,*) c1,c2,xlamda,cfp(1,1),gmt(1,1)  

!        end do
!      end do
            

            
            
!       return 
!       end  


!----------------------Jiahao SSD_GND_cos_sin------------------------------

        subroutine screw_edge_cos_sin(screw_cos,screw_sin,
     &  edge_cos,edge_sin)
     
      real(8):: screw_cos(12,12), screw_sin(12,12)
      real(8):: edge_cos(12,12),  edge_sin(12,12)
      integer:: i
      
      i=1
      edge_cos(i,1)=0.0d0
      edge_cos(i,2)=0.0d0
      edge_cos(i,3)=0.0d0
      edge_cos(i,4)=1.0d0
      edge_cos(i,5)=1.0d0
      edge_cos(i,6)=1.0d0
      edge_cos(i,7)=0.0d0
      edge_cos(i,8)=0.0d0
      edge_cos(i,9)=0.0d0
      edge_cos(i,10)=0.0d0
      edge_cos(i,11)=0.0d0
      edge_cos(i,12)=0.0d0
      
      i=2     
      edge_cos(i,1)=0.0d0
      edge_cos(i,2)=0.0d0
      edge_cos(i,3)=0.0d0
      edge_cos(i,4)=1.0d0
      edge_cos(i,5)=1.0d0
      edge_cos(i,6)=1.0d0
      edge_cos(i,7)=0.0d0
      edge_cos(i,8)=0.0d0
      edge_cos(i,9)=0.0d0
      edge_cos(i,10)=0.0d0
      edge_cos(i,11)=0.0d0
      edge_cos(i,12)=0.0d0
      
      i=3
      edge_cos(i,1)=0.0d0
      edge_cos(i,2)=0.0d0
      edge_cos(i,3)=0.0d0
      edge_cos(i,4)=1.0d0
      edge_cos(i,5)=1.0d0
      edge_cos(i,6)=1.0d0
      edge_cos(i,7)=0.0d0
      edge_cos(i,8)=0.0d0
      edge_cos(i,9)=0.0d0
      edge_cos(i,10)=0.0d0
      edge_cos(i,11)=0.0d0
      edge_cos(i,12)=0.0d0
      
            i=4
      edge_cos(i,1)=0.5d0
      edge_cos(i,2)=0.5d0
      edge_cos(i,3)=1.0d0
      edge_cos(i,4)=0.0d0
      edge_cos(i,5)=0.0d0
      edge_cos(i,6)=0.0d0
      edge_cos(i,7)=0.5d0
      edge_cos(i,8)=1.0d0
      edge_cos(i,9)=0.5d0
      edge_cos(i,10)=0.5d0
      edge_cos(i,11)=1.0d0
      edge_cos(i,12)=0.5d0
      
            i=5
      edge_cos(i,1)=0.5d0
      edge_cos(i,2)=1.0d0
      edge_cos(i,3)=0.5d0
      edge_cos(i,4)=0.0d0
      edge_cos(i,5)=0.0d0
      edge_cos(i,6)=0.0d0
      edge_cos(i,7)=0.5d0
      edge_cos(i,8)=0.5d0
      edge_cos(i,9)=1.0d0
      edge_cos(i,10)=0.5d0
      edge_cos(i,11)=0.5d0
      edge_cos(i,12)=1.0d0
      
            i=6
      edge_cos(i,1)=1.0d0
      edge_cos(i,2)=0.5d0
      edge_cos(i,3)=0.5d0
      edge_cos(i,4)=0.0d0
      edge_cos(i,5)=0.0d0
      edge_cos(i,6)=0.0d0
      edge_cos(i,7)=1.0d0
      edge_cos(i,8)=0.5d0
      edge_cos(i,9)=0.5d0
      edge_cos(i,10)=1.0d0
      edge_cos(i,11)=0.5d0
      edge_cos(i,12)=0.5d0
      
            i=7
      edge_cos(i,1)=0.0d0
      edge_cos(i,2)=0.7374
      edge_cos(i,3)=0.7374
      edge_cos(i,4)=0.5243
      edge_cos(i,5)=0.5243
      edge_cos(i,6)=0.5243
      edge_cos(i,7)=0.0d0
      edge_cos(i,8)=0.7374
      edge_cos(i,9)=0.7374
      edge_cos(i,10)=0.0d0
      edge_cos(i,11)=0.7374
      edge_cos(i,12)=0.7374
      
            i=8
      edge_cos(i,1)=0.7374
      edge_cos(i,2)=0.7374
      edge_cos(i,3)=0.0d0
      edge_cos(i,4)=0.5243
      edge_cos(i,5)=0.5243
      edge_cos(i,6)=0.5243
      edge_cos(i,7)=0.7374
      edge_cos(i,8)=0.0d0
      edge_cos(i,9)=0.7374
      edge_cos(i,10)=0.7374
      edge_cos(i,11)=0.0d0
      edge_cos(i,12)=0.7374
      
            i=9
      edge_cos(i,1)=0.7374
      edge_cos(i,2)=0.0d0
      edge_cos(i,3)=0.7374
      edge_cos(i,4)=0.5243
      edge_cos(i,5)=0.5243
      edge_cos(i,6)=0.5243
      edge_cos(i,7)=0.7374
      edge_cos(i,8)=0.7374
      edge_cos(i,9)=0.0d0
      edge_cos(i,10)=0.7374
      edge_cos(i,11)=0.7374
      edge_cos(i,12)=0.0d0
      
            i=10
      edge_cos(i,1)=0.0d0
      edge_cos(i,2)=0.7374
      edge_cos(i,3)=0.7374
      edge_cos(i,4)=0.5243
      edge_cos(i,5)=0.5243
      edge_cos(i,6)=0.5243
      edge_cos(i,7)=0.0d0
      edge_cos(i,8)=0.7374
      edge_cos(i,9)=0.7374
      edge_cos(i,10)=0.0d0
      edge_cos(i,11)=0.7374
      edge_cos(i,12)=0.7374
      
           i=11
      edge_cos(i,1)=0.7374
      edge_cos(i,2)=0.7374
      edge_cos(i,3)=0.0d0
      edge_cos(i,4)=0.5243
      edge_cos(i,5)=0.5243
      edge_cos(i,6)=0.5243
      edge_cos(i,7)=0.7374
      edge_cos(i,8)=0.0d0
      edge_cos(i,9)=0.7374
      edge_cos(i,10)=0.7374
      edge_cos(i,11)=0.0d0
      edge_cos(i,12)=0.7374
      
           i=12
      edge_cos(i,1)=0.7374
      edge_cos(i,2)=0.0d0
      edge_cos(i,3)=0.7374
      edge_cos(i,4)=0.5243
      edge_cos(i,5)=0.5243
      edge_cos(i,6)=0.5243
      edge_cos(i,7)=0.7374
      edge_cos(i,8)=0.7374
      edge_cos(i,9)=0.0d0
      edge_cos(i,10)=0.7374
      edge_cos(i,11)=0.7374
      edge_cos(i,12)=0.0d0
      
      i=1
      edge_sin(i,1)=1.0d0
      edge_sin(i,2)=1.0d0
      edge_sin(i,3)=1.0d0
      edge_sin(i,4)=0.0d0
      edge_sin(i,5)=0.0d0
      edge_sin(i,6)=0.0d0
      edge_sin(i,7)=1.0d0
      edge_sin(i,8)=1.0d0
      edge_sin(i,9)=1.0d0
      edge_sin(i,10)=1.0d0
      edge_sin(i,11)=1.0d0
      edge_sin(i,12)=1.0d0
      
      i=2
      edge_sin(i,1)=1.0d0
      edge_sin(i,2)=1.0d0
      edge_sin(i,3)=1.0d0
      edge_sin(i,4)=0.0d0
      edge_sin(i,5)=0.0d0
      edge_sin(i,6)=0.0d0
      edge_sin(i,7)=1.0d0
      edge_sin(i,8)=1.0d0
      edge_sin(i,9)=1.0d0
      edge_sin(i,10)=1.0d0
      edge_sin(i,11)=1.0d0
      edge_sin(i,12)=1.0d0
      
            i=3
      edge_sin(i,1)=1.0d0
      edge_sin(i,2)=1.0d0
      edge_sin(i,3)=1.0d0
      edge_sin(i,4)=0.0d0
      edge_sin(i,5)=0.0d0
      edge_sin(i,6)=0.0d0
      edge_sin(i,7)=1.0d0
      edge_sin(i,8)=1.0d0
      edge_sin(i,9)=1.0d0
      edge_sin(i,10)=1.0d0
      edge_sin(i,11)=1.0d0
      edge_sin(i,12)=1.0d0
      
            i=4
      edge_sin(i,1)=0.8660
      edge_sin(i,2)=0.8660
      edge_sin(i,3)=0.0d0
      edge_sin(i,4)=1.0d0
      edge_sin(i,5)=1.0d0
      edge_sin(i,6)=1.0d0
      edge_sin(i,7)=0.8660
      edge_sin(i,8)=0.0d0
      edge_sin(i,9)=0.8660
      edge_sin(i,10)=0.8660
      edge_sin(i,11)=0.0d0
      edge_sin(i,12)=0.8660
      
            i=5
      edge_sin(i,1)=0.8660
      edge_sin(i,2)=0.8660
      edge_sin(i,3)=0.0d0
      edge_sin(i,4)=1.0d0
      edge_sin(i,5)=1.0d0
      edge_sin(i,6)=1.0d0
      edge_sin(i,7)=0.8660
      edge_sin(i,8)=0.8660
      edge_sin(i,9)=0.0d0
      edge_sin(i,10)=0.8660
      edge_sin(i,11)=0.8660
      edge_sin(i,12)=0.0d0
      
            i=6
      edge_sin(i,1)=0.0d0
      edge_sin(i,2)=0.8660
      edge_sin(i,3)=0.8660
      edge_sin(i,4)=1.0d0
      edge_sin(i,5)=1.0d0
      edge_sin(i,6)=1.0d0
      edge_sin(i,7)=0.0d0
      edge_sin(i,8)=0.8660
      edge_sin(i,9)=0.8660
      edge_sin(i,10)=0.0d0
      edge_sin(i,11)=0.8660
      edge_sin(i,12)=0.8660
      
            i=7
      edge_sin(i,1)=1.0d0
      edge_sin(i,2)=0.6754
      edge_sin(i,3)=0.6754
      edge_sin(i,4)=0.8515
      edge_sin(i,5)=0.8515
      edge_sin(i,6)=0.8515
      edge_sin(i,7)=1.0d0
      edge_sin(i,8)=0.6754
      edge_sin(i,9)=0.6754
      edge_sin(i,10)=1.0d0
      edge_sin(i,11)=0.6754
      edge_sin(i,12)=0.6754
      
            i=8
      edge_sin(i,1)=0.6754
      edge_sin(i,2)=0.6754
      edge_sin(i,3)=1.0d0
      edge_sin(i,4)=0.8515
      edge_sin(i,5)=0.8515
      edge_sin(i,6)=0.8515
      edge_sin(i,7)=0.6754
      edge_sin(i,8)=1.0d0
      edge_sin(i,9)=0.6754
      edge_sin(i,10)=0.6754
      edge_sin(i,11)=1.0d0
      edge_sin(i,12)=0.6754
      
            i=9
      edge_sin(i,1)=0.6754
      edge_sin(i,2)=1.0d0
      edge_sin(i,3)=0.6754
      edge_sin(i,4)=0.8515
      edge_sin(i,5)=0.8515
      edge_sin(i,6)=0.8515
      edge_sin(i,7)=0.6754
      edge_sin(i,8)=0.6754
      edge_sin(i,9)=1.0d0
      edge_sin(i,10)=0.6754
      edge_sin(i,11)=0.6754
      edge_sin(i,12)=1.0d0
      
            i=10
      edge_sin(i,1)=1.0d0
      edge_sin(i,2)=0.6754
      edge_sin(i,3)=0.6754
      edge_sin(i,4)=0.8515
      edge_sin(i,5)=0.8515
      edge_sin(i,6)=0.8515
      edge_sin(i,7)=1.0d0
      edge_sin(i,8)=0.6754
      edge_sin(i,9)=0.6754
      edge_sin(i,10)=1.0d0
      edge_sin(i,11)=0.6754
      edge_sin(i,12)=0.6754
      
            i=11
      edge_sin(i,1)=0.6754
      edge_sin(i,2)=0.6754
      edge_sin(i,3)=1.0d0
      edge_sin(i,4)=0.8515
      edge_sin(i,5)=0.8515
      edge_sin(i,6)=0.8515
      edge_sin(i,7)=0.6754
      edge_sin(i,8)=1.0d0
      edge_sin(i,9)=0.6754
      edge_sin(i,10)=0.6754
      edge_sin(i,11)=1.0d0
      edge_sin(i,12)=0.6754
      
            i=12
      edge_sin(i,1)=0.6754
      edge_sin(i,2)=1.0d0
      edge_sin(i,3)=0.6754
      edge_sin(i,4)=0.8515
      edge_sin(i,5)=0.8515
      edge_sin(i,6)=0.8515
      edge_sin(i,7)=0.6754
      edge_sin(i,8)=0.6754
      edge_sin(i,9)=1.0d0
      edge_sin(i,10)=0.6754
      edge_sin(i,11)=0.6754
      edge_sin(i,12)=1.0d0
      
      
      
      
           i=1
      screw_cos(i,1)=0.0d0
      screw_cos(i,2)=0.0d0
      screw_cos(i,3)=0.0d0
      screw_cos(i,4)=0.0d0
      screw_cos(i,5)=0.0d0
      screw_cos(i,6)=0.0d0
      screw_cos(i,7)=0.8514
      screw_cos(i,8)=0.8514
      screw_cos(i,9)=0.8514
      screw_cos(i,10)=0.8514
      screw_cos(i,11)=0.8514
      screw_cos(i,12)=0.8514
      
      i=2     
      screw_cos(i,1)=0.0d0
      screw_cos(i,2)=0.0d0
      screw_cos(i,3)=0.0d0
      screw_cos(i,4)=0.0d0
      screw_cos(i,5)=0.0d0
      screw_cos(i,6)=0.0d0
      screw_cos(i,7)=0.8514
      screw_cos(i,8)=0.8514
      screw_cos(i,9)=0.8514
      screw_cos(i,10)=0.8514
      screw_cos(i,11)=0.8514
      screw_cos(i,12)=0.8514
      
      i=3
      screw_cos(i,1)=0.0d0
      screw_cos(i,2)=0.0d0
      screw_cos(i,3)=0.0d0
      screw_cos(i,4)=0.0d0
      screw_cos(i,5)=0.0d0
      screw_cos(i,6)=0.0d0
      screw_cos(i,7)=0.8514
      screw_cos(i,8)=0.8514
      screw_cos(i,9)=0.8514
      screw_cos(i,10)=0.8514
      screw_cos(i,11)=0.8514
      screw_cos(i,12)=0.8514
      
            i=4
      screw_cos(i,1)=0.8660
      screw_cos(i,2)=0.8660
      screw_cos(i,3)=0.0d0
      screw_cos(i,4)=0.0d0
      screw_cos(i,5)=0.8660
      screw_cos(i,6)=0.8660
      screw_cos(i,7)=0.4542
      screw_cos(i,8)=0.0d0
      screw_cos(i,9)=0.4542
      screw_cos(i,10)=0.4542
      screw_cos(i,11)=0.0d0
      screw_cos(i,12)=0.4542
      
            i=5
      screw_cos(i,1)=0.8660
      screw_cos(i,2)=0.0d0
      screw_cos(i,3)=0.8660
      screw_cos(i,4)=0.8660
      screw_cos(i,5)=0.0d0
      screw_cos(i,6)=0.8660
      screw_cos(i,7)=0.4542
      screw_cos(i,8)=0.4542
      screw_cos(i,9)=0.0d0
      screw_cos(i,10)=0.4542
      screw_cos(i,11)=0.4542
      screw_cos(i,12)=0.0d0
      
            i=6
      screw_cos(i,1)=0.0d0
      screw_cos(i,2)=0.8660
      screw_cos(i,3)=0.8660
      screw_cos(i,4)=0.8660
      screw_cos(i,5)=0.8660
      screw_cos(i,6)=0.0d0
      screw_cos(i,7)=0.0d0
      screw_cos(i,8)=0.4542
      screw_cos(i,9)=0.4542
      screw_cos(i,10)=0.0d0
      screw_cos(i,11)=0.4542
      screw_cos(i,12)=0.4542
      
            i=7
      screw_cos(i,1)=0.8514
      screw_cos(i,2)=0.4257
      screw_cos(i,3)=0.4257
      screw_cos(i,4)=0.4257
      screw_cos(i,5)=0.4257
      screw_cos(i,6)=0.8514
      screw_cos(i,7)=0.0d0
      screw_cos(i,8)=0.2233
      screw_cos(i,9)=0.6698
      screw_cos(i,10)=0.8931
      screw_cos(i,11)=0.6698
      screw_cos(i,12)=0.2233
      
            i=8
      screw_cos(i,1)=0.4257
      screw_cos(i,2)=0.4257
      screw_cos(i,3)=0.8514
      screw_cos(i,4)=0.8514
      screw_cos(i,5)=0.4257
      screw_cos(i,6)=0.4257
      screw_cos(i,7)=0.2233
      screw_cos(i,8)=0.0d0
      screw_cos(i,9)=0.2233
      screw_cos(i,10)=0.6698
      screw_cos(i,11)=0.8931
      screw_cos(i,12)=0.6698
      
            i=9
      screw_cos(i,1)=0.4257
      screw_cos(i,2)=0.8514
      screw_cos(i,3)=0.4257
      screw_cos(i,4)=0.4257
      screw_cos(i,5)=0.8514
      screw_cos(i,6)=0.4257
      screw_cos(i,7)=0.6698
      screw_cos(i,8)=0.2233
      screw_cos(i,9)=0.0d0
      screw_cos(i,10)=0.2233
      screw_cos(i,11)=0.6698
      screw_cos(i,12)=0.8931
      
            i=10
      screw_cos(i,1)=0.8514
      screw_cos(i,2)=0.4257
      screw_cos(i,3)=0.4257
      screw_cos(i,4)=0.4257
      screw_cos(i,5)=0.4257
      screw_cos(i,6)=0.8514
      screw_cos(i,7)=0.8931
      screw_cos(i,8)=0.6698
      screw_cos(i,9)=0.2233
      screw_cos(i,10)=0.0d0
      screw_cos(i,11)=0.2233
      screw_cos(i,12)=0.6698
      
           i=11
      screw_cos(i,1)=0.4257
      screw_cos(i,2)=0.4257
      screw_cos(i,3)=0.8514
      screw_cos(i,4)=0.8514
      screw_cos(i,5)=0.4257
      screw_cos(i,6)=0.4257
      screw_cos(i,7)=0.6698
      screw_cos(i,8)=0.8931
      screw_cos(i,9)=0.6698
      screw_cos(i,10)=0.2233
      screw_cos(i,11)=0.0d0
      screw_cos(i,12)=0.2233
      
           i=12
      screw_cos(i,1)=0.4257
      screw_cos(i,2)=0.8514
      screw_cos(i,3)=0.4257
      screw_cos(i,4)=0.4257
      screw_cos(i,5)=0.8514
      screw_cos(i,6)=0.4257
      screw_cos(i,7)=0.2233
      screw_cos(i,8)=0.6698
      screw_cos(i,9)=0.8931
      screw_cos(i,10)=0.6698
      screw_cos(i,11)=0.2233
      screw_cos(i,12)=0.0d0
      
      i=1
      screw_sin(i,1)=1.0d0
      screw_sin(i,2)=1.0d0
      screw_sin(i,3)=1.0d0
      screw_sin(i,4)=1.0d0
      screw_sin(i,5)=1.0d0
      screw_sin(i,6)=1.0d0
      screw_sin(i,7)=0.5244
      screw_sin(i,8)=0.5244
      screw_sin(i,9)=0.5244
      screw_sin(i,10)=0.5244
      screw_sin(i,11)=0.5244
      screw_sin(i,12)=0.5244
      
      i=2
      screw_sin(i,1)=1.0d0
      screw_sin(i,2)=1.0d0
      screw_sin(i,3)=1.0d0
      screw_sin(i,4)=1.0d0
      screw_sin(i,5)=1.0d0
      screw_sin(i,6)=1.0d0
      screw_sin(i,7)=0.5244
      screw_sin(i,8)=0.5244
      screw_sin(i,9)=0.5244
      screw_sin(i,10)=0.5244
      screw_sin(i,11)=0.5244
      screw_sin(i,12)=0.5244      
            i=3
      screw_sin(i,1)=1.0d0
      screw_sin(i,2)=1.0d0
      screw_sin(i,3)=1.0d0
      screw_sin(i,4)=1.0d0
      screw_sin(i,5)=1.0d0
      screw_sin(i,6)=1.0d0
      screw_sin(i,7)=0.5244
      screw_sin(i,8)=0.5244
      screw_sin(i,9)=0.5244
      screw_sin(i,10)=0.5244
      screw_sin(i,11)=0.5244
      screw_sin(i,12)=0.5244
      
            i=4
      screw_sin(i,1)=0.5d0
      screw_sin(i,2)=0.5d0
      screw_sin(i,3)=1.0d0
      screw_sin(i,4)=1.0d0
      screw_sin(i,5)=0.5d0
      screw_sin(i,6)=0.5d0
      screw_sin(i,7)=0.8909
      screw_sin(i,8)=1.0d0
      screw_sin(i,9)=0.8909
      screw_sin(i,10)=0.8909
      screw_sin(i,11)=1.0d0
      screw_sin(i,12)=0.8909
      
            i=5
      screw_sin(i,1)=0.5d0
      screw_sin(i,2)=1.0d0
      screw_sin(i,3)=0.5d0
      screw_sin(i,4)=0.5d0
      screw_sin(i,5)=1.0d0
      screw_sin(i,6)=0.5d0
      screw_sin(i,7)=0.8909
      screw_sin(i,8)=0.8909
      screw_sin(i,9)=1.0d0
      screw_sin(i,10)=0.8909
      screw_sin(i,11)=0.8909
      screw_sin(i,12)=1.0d0
      
            i=6
      screw_sin(i,1)=1.0d0
      screw_sin(i,2)=0.5d0
      screw_sin(i,3)=0.5d0
      screw_sin(i,4)=0.5d0
      screw_sin(i,5)=0.5d0
      screw_sin(i,6)=1.0d0
      screw_sin(i,7)=1.0d0
      screw_sin(i,8)=0.8909
      screw_sin(i,9)=0.8909
      screw_sin(i,10)=1.0d0
      screw_sin(i,11)=0.8909
      screw_sin(i,12)=0.8909
      
            i=7
      screw_sin(i,1)=0.5244
      screw_sin(i,2)=0.9049
      screw_sin(i,3)=0.9049
      screw_sin(i,4)=0.9049
      screw_sin(i,5)=0.9049
      screw_sin(i,6)=0.5244
      screw_sin(i,7)=1.0d0
      screw_sin(i,8)=0.9748
      screw_sin(i,9)=0.7425
      screw_sin(i,10)=0.4499
      screw_sin(i,11)=0.7425
      screw_sin(i,12)=0.9748
      
            i=8
      screw_sin(i,1)=0.9049
      screw_sin(i,2)=0.9049
      screw_sin(i,3)=0.5244
      screw_sin(i,4)=0.5244
      screw_sin(i,5)=0.9049
      screw_sin(i,6)=0.9049
      screw_sin(i,7)=0.9748
      screw_sin(i,8)=1.0d0
      screw_sin(i,9)=0.9748
      screw_sin(i,10)=0.7425
      screw_sin(i,11)=0.4499
      screw_sin(i,12)=0.7425
      
            i=9
      screw_sin(i,1)=0.9049
      screw_sin(i,2)=0.5244
      screw_sin(i,3)=0.9049
      screw_sin(i,4)=0.9049
      screw_sin(i,5)=0.5244
      screw_sin(i,6)=0.9049
      screw_sin(i,7)=0.7425
      screw_sin(i,8)=0.9748
      screw_sin(i,9)=1.0d0
      screw_sin(i,10)=0.9748
      screw_sin(i,11)=0.7425
      screw_sin(i,12)=0.4499
      
            i=10
      screw_sin(i,1)=0.5244
      screw_sin(i,2)=0.9049
      screw_sin(i,3)=0.9049
      screw_sin(i,4)=0.9049
      screw_sin(i,5)=0.9049
      screw_sin(i,6)=0.5244
      screw_sin(i,7)=0.4499
      screw_sin(i,8)=0.7425
      screw_sin(i,9)=0.9748
      screw_sin(i,10)=1.0d0
      screw_sin(i,11)=0.9748
      screw_sin(i,12)=0.7425
      
            i=11
      screw_sin(i,1)=0.9049
      screw_sin(i,2)=0.9049
      screw_sin(i,3)=0.5244
      screw_sin(i,4)=0.5244
      screw_sin(i,5)=0.9049
      screw_sin(i,6)=0.9049
      screw_sin(i,7)=0.7425
      screw_sin(i,8)=0.4499
      screw_sin(i,9)=0.7425
      screw_sin(i,10)=0.9748
      screw_sin(i,11)=1.0d0
      screw_sin(i,12)=0.9748
     
            i=12
      screw_sin(i,1)=0.9049
      screw_sin(i,2)=0.5244
      screw_sin(i,3)=0.9049
      screw_sin(i,4)=0.9049
      screw_sin(i,5)=0.5244
      screw_sin(i,6)=0.9049
      screw_sin(i,7)=0.9748
      screw_sin(i,8)=0.7425
      screw_sin(i,9)=0.4499
      screw_sin(i,10)=0.7425
      screw_sin(i,11)=0.9748
      screw_sin(i,12)=1.0d0
      
      return
      end
      
      
      subroutine twin_nucleat(twin_nu,dfgrd1,tpk_1d,fp,s_tan_tw,
     &  gsize, s_gg_2d,S_critical,grainID,iflag)
      
      implicit none
      
      integer:: nslip
      parameter (nslip=30)
      real(8):: twin_nu, dfgrd1(3,3), tpk_1d(6)
      real(8):: fp(3,3),e_strain(6),e_strain_2d(3,3)
      real(8):: fe_t_fe(3,3), tpk_2d(3,3), fe_t_fe_tpk(3,3)
      real(8):: tau_twin(6), s_tan_tw(3,3,6), gsize
      real(8):: tau_critical, alpha, S_min, S_max
      real(8):: n_star,Y,S_critical, s_gg_2d(6,6)
      real(8):: tau_twin_1(6), rand_chance(46)
      
      integer:: i,j,k,n,isys,iflag
      integer:: grainID
      
      alpha = 2.0
      tau_critical = 100.0     !Mpa
      n_star=3.1415926*gsize/2.0d0    !micron
      
!      CALL RANDOM_SEED()
!      CALL RANDOM_NUMBER(Y)

      rand_chance(1)= 0.445586200710899
      rand_chance(2)=0.646313010111265
      rand_chance(3)=0.709364830858073
      rand_chance(4)=0.754686681982361
      rand_chance(5)=0.276025076998578
      rand_chance(6)=0.679702676853675
      rand_chance(7)=0.278498218867048
      rand_chance(8)=0.546881519204984
      rand_chance(9)=0.957506835434298
      rand_chance(10)=0.964888535199277
      rand_chance(11)=0.157613081677548
      rand_chance(12)=0.970592781760616
      rand_chance(13)=0.957166948242946
      rand_chance(14)=0.485375648722841
      rand_chance(15)= 0.800280468888800
      rand_chance(16)=0.141886338627215
      rand_chance(17)= 0.421761282626275
      rand_chance(18)=0.915735525189067
       rand_chance(19)=0.792207329559554
       rand_chance(20)=0.959492426392903
       rand_chance(21)=0.655740699156587
       rand_chance(22)=0.035711678574190
       rand_chance(23)=0.849129305868777
       rand_chance(24)=0.933993247757551
       rand_chance(25)=0.678735154857773
       rand_chance(26)=0.757740130578333
       rand_chance(27)=0.743132468124916
       rand_chance(28)= 0.392227019534168
       rand_chance(29)=0.655477890177557
       rand_chance(30)=0.171186687811562
      rand_chance(31)= 0.706046088019609
      rand_chance(32)= 0.031832846377421
       rand_chance(33)= 0.276922984960890
       rand_chance(34)=0.046171390631154
       rand_chance(35)= 0.097131781235848
       rand_chance(36)= 0.823457828327293
       rand_chance(37)= 0.694828622975817
       rand_chance(38)=0.317099480060861
       rand_chance(39)=0.950222048838355
       rand_chance(40)=0.034446080502909
       rand_chance(41)=0.438744359656398
       rand_chance(42)=0.381558457093008
       rand_chance(43)=0.765516788149002
       rand_chance(44)=0.795199901137063
       rand_chance(45)=0.186872604554379
       rand_chance(46)=0.489764395788231

!      Y=rand_chance(grainID)
      Y=0.2
      
      S_min=tau_critical*(-log(1.0-Y)/n_star)**(1.0/alpha)
      S_max=tau_critical*(-log(1.0-Y**(1.0/n_star)))**(1.0/alpha)
      S_critical=S_min+((S_max-S_min)/3.0)
      

!----------------------------------------------------------      

       if (iflag.eq.2) then        
      
         call matmul(s_gg_2d,tpk_1d,e_strain,6,6,1)
         call sigmat(e_strain,e_strain_2d)
         do i=1,3
            do j=1,3
               if(i.eq.j)then
                  fe_t_fe(i,j)=1.d0+(2.0*e_strain_2d(i,j))
               else
                  fe_t_fe(i,j)=2.0*e_strain_2d(i,j)
               endif
               enddo
               enddo        
         call sigmat(tpk_1d,tpk_2d)
         call mat33(fe_t_fe_tpk,fe_t_fe,tpk_2d,3)
c
c-------------------------------------------------------------
c
	  call sigmat(tpk_1d,tpk_2d)
	 do k=1,6
	      tau_twin(k)=0.0d0
	      do i=1,3
		 do j=1,3
		    tau_twin(k)=tau_twin(k)+tpk_2d(i,j)*s_tan_tw(i,j,k)
!      tau_twin_1(k)=tau_twin(k)+fe_t_fe_tpk(i,j)*s_tan_tw(i,j,k)
		 enddo
	      enddo
	      
	  
!	  write(4002,*) 'tau_twin is', tau_twin(k) 
!       write(4002,*) 'S_critical is', S_critical  
	      
	  if (abs(tau_twin(k)).gt.S_critical) then
        twin_nu=0.0
        endif
        
       enddo

      endif

       return
       end
