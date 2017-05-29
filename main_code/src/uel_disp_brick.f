
	subroutine uel_disp_b(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period,iflag,
     4 ngauss,Fp_dgamma_local,euler)
 

	include 'aba_param.inc'
	include 'pardis.h'
	parameter (ntens=6)
	dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),
     1     svars(msv),energy(8),coords(mcrd,nnode),u(ndofel),
     2     du(mlvarx,*),v(ndofel),a(ndofel),time(2),params(*),
     3     jdltyp(mdload,*),adlmag(mdload,*),ddlmag(mdload,*),
     4     predef(2,npredf,nnode),Lflags(*),jprops(*)

	dimension u_j(ndofel),f1_j(24),a11matrx(nnode*mcrd,nnode*mcrd)
     $     ,rhs1(mlvarx,1),amatrx1(ndofel,ndofel)
     $     ,svars1(msv),xii(8),etai(8),chii(8),xi(8),eta(8),chi(8)
     $     ,f1(24),du_j(mlvarx,1),ske(mcrd*nnode,mcrd*nnode)
     &     , ske_j(mcrd*nnode,mcrd*nnode)
     
      real(8)::  Fp_dgamma_local(4,12),euler(3)
	character*80 cmname
	dimension stress(6),
     1     ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),drot(3,3),
     2     stran(ntens),dstran(ntens),predef1(1),dpred(1),coords_ipt(3)

	common/wr/iwr
	common/gauss/xii,etai,chii,xi,eta,chi
	common/inits/igauss_init
	common/conv/iconv1
	common/numj/inumj
!	data igauss_init/0/
!	data iinp/0/

	if(igauss_init.eq.0)then
	   igauss_init=1
	   call init_gauss(mcrd)
	endif
	if(lflags(3).eq.4)then
	   amatrx(1:ndofel,1:ndofel)=0.d0
	   return
	endif
	inumj=0
	call calc_resid_disp(rhs,amatrx,svars,energy,ndofel,nrhs
     $     ,nsvars,props,nprops,coords,mcrd,nnode,u,du,v,a,jtype
     $     ,time,dtime,kstep,kinc,jelem,params,ndload,jdltyp
     $     ,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload
     $     ,pnewdt,jprops,njprop,period,f1,ske,ngauss,euler)
	
	if(pnewdt.lt.1.d0)return
	
	rhs(1:ndofel,1)=-f1(1:ndofel)
	amatrx(1:ndofel,1:ndofel)=ske(1:ndofel,1:ndofel)

	
	return
	end



      
	subroutine calc_resid_disp(rhs,amatrx,svars,energy,ndofel,nrhs,
     1 nsvars,props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period,f1,ske,
     4 ngauss,euler)

	include 'aba_param.inc'
	include 'pardis.h'
	parameter (ntens=6)

	dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),
     1     svars(msv),energy(8),coords(mcrd,nnode),u(ndofel),
     2     du(mlvarx,*),v(ndofel),a(ndofel),time(2),params(*),
     3     jdltyp(mdload,*),adlmag(mdload,*),ddlmag(mdload,*),
     4     predef(2,npredf,nnode),lflags(*),jprops(*)
      
	dimension stress(6),ddsdde(ntens,ntens),ddsddt(ntens),
     &   drplde(ntens), drot(3,3), stran(ntens),dstran(ntens),
     &   predef1(1),dpred(1) 
      real(8):: coords_ipt(3)
      
	dimension ue(mcrd,nnode),dfgrd1(3,3),dfgrd0(3
     $     ,3),ue_t(mcrd,nnode),xii(8),etai(8),chii(8),xi(8),eta(8)
     $     ,chi(8),dfgrd1_inv(3,3),pk1_2d(3,3),pk1_1d(9),xb(9,24),f1(24)
     $     ,shp(3,8),shpf(8),sigma(3,3)

	dimension xe_t(mcrd,nnode),xe_tau(mcrd,nnode),xe_mid(mcrd,nnode),
     &     shp0(3,8),shpf0(8),shptau(3,8),shpftau(8),shpmid(3,8),
     &     shpfmid(8),shp0_st(3,8,8),shptau_st(3,8,8),shpmid_st(3,8,8),
     &     shpf0_st(8,8),shpftau_st(8,8),shpfmid_st(8,8),dvol0_st(8),
     &     dvoltau_st(8),dvolmid_st(8),xb_av(mcrd*nnode)
     &     ,xb_tr(mcrd*nnode)
	dimension dfgrd1_st(3,3,8),dfgrd0_st(3,3,8),xcmid(mcrd),
     &    xbar_st(2*mcrd,mcrd*nnode,8),xc0(mcrd),xctau(mcrd),
     &    xjac_st(nnode),xbtr_st(mcrd*nnode,8),
     &    ske2(mcrd*nnode,mcrd*nnode),ske3(mcrd*nnode,mcrd*nnode),
     &    ske(mcrd*nnode,mcrd*nnode),stmat(mcrd**2,mcrd**2),
     &    stmat1(2*mcrd,2*mcrd),xmtemp1(mcrd**2,mcrd*nnode),
     &    rt_tau(3,3),ut_tau(3,3),dfgrdt(3,3),crot(2*mcrd,2*mcrd),
     &    xtm(2*mcrd,mcrd*nnode),drotm(2*mcrd,2*mcrd),ftemp(mcrd*nnode),
     &    xg_st(mcrd**2,mcrd*nnode,8),dfgrd0_inv(3,3),
     &    skt(ndofel,ndofel),ct(ndofel,ndofel),statev(maxstatev),
     &    ske1(mcrd*nnode,mcrd*nnode),  xbar(ndofel,ndofel)

	dimension ibelem(maxel),pressb(maxel),ilf(maxel),iface_tet(12)
	dimension xcg(3),uv1(3),uv2(3),uv3(1:3),uv4(1:3),ibnode(4)
	dimension shp_surf(3,12),fload(12),shp_surf_xi(3,12),
     &     shp_surf_eta(3,12),dx_xi(3),dx_eta(3),ske_load(12,12),
     &     xel_surf(12),da(3),tempv(3),tempv1(3),tempv2(3),
     &     iface_brick(24)

      real(8):: euler(3), gradient_GND(3,12)

	dimension shp_surfj(3,12),floadj(12),shp_surf_xij(3,12),
     &     shp_surf_etaj(3,12),ske_loadj(12,12),
     &     xel_surfj(12),daj(3),skdiff(12,12)
	

	common/press/pressb,ibelem,ilf,iface_tet,iface_brick


	character*80 cmname      

	common/gauss/xii,etai,chii,xi,eta,chi
	common/inits/igauss_init
	common/conv/iconv1
	common/numj/inumj


	ue(1:mcrd,1:nnode)=reshape(u,(/mcrd,nnode/))
	ue_t(1:mcrd,1:nnode)=ue(1:mcrd,1:nnode)-
     &  reshape(du(1:mcrd*nnode,1),(/mcrd,nnode/))
	xe_tau(1:mcrd,1:nnode)=coords(1:mcrd,1:nnode)+ue(1:mcrd,1:nnode)
	xe_t(1:mcrd,1:nnode)=coords(1:mcrd,1:nnode)+ue_t(1:mcrd,1:nnode)
	xe_mid(1:mcrd,1:nnode)=0.5d0*(xe_t(1:mcrd,1:nnode)
     & +xe_tau(1:mcrd,1:nnode))	
	
	
	tvol0=0.d0
	tvoltau=0.d0
	tvolmid=0.d0
	if(mcrd.eq.2)nstr=3
	if(mcrd.eq.3)nstr=6

	ske(1:ndofel,1:ndofel)=0.d0
	f1(1:ndofel)=0.d0

	
	skt(1:ndofel,1:ndofel)=0.d0
	do ipt=1,ngauss
	   if(mcrd.eq.2)then
	      call shape2(coords,ipt,shp0,shpf0,dvol0,xc0)
	      call shape2(xe_tau,ipt,shptau,shpftau,dvoltau,xctau)
	      call shape2(xe_mid,ipt,shpmid,shpfmid,dvolmid,xcmid)
	   endif
	   if(mcrd.eq.3)then
	      call shape3(coords,ipt,shp0,shpf0,dvol0,xc0)
	      call shape3(xe_tau,ipt,shptau,shpftau,dvoltau,xctau)
	      call shape3(xe_mid,ipt,shpmid,shpfmid,dvolmid,xcmid)
	   endif
	   dvol0_st(ipt)=dvol0
	   dvoltau_st(ipt)=dvoltau
	   dvolmid_st(ipt)=dvolmid
	   tvol0=tvol0+dvol0
	   tvoltau=tvoltau+dvoltau
	   tvolmid=tvolmid+dvolmid
	   shp0_st(1:mcrd,1:nnode,ipt)=shp0(1:mcrd,1:nnode)
	   shptau_st(1:mcrd,1:nnode,ipt)=shptau(1:mcrd,1:nnode)
	   shpmid_st(1:mcrd,1:nnode,ipt)=shpmid(1:mcrd,1:nnode)
	   shpf0_st(1:nnode,ipt)=shpf0(1:nnode)
	   shpftau_st(1:nnode,ipt)=shpftau(1:nnode)
	   shpfmid_st(1:nnode,ipt)=shpfmid(1:nnode)
	enddo


	call calc_dfg(ue_t,nnode,ndofel,mcrd,dfgrd0_st,dvol0_st,xjac_st,
     &     shp0_st,xjbar,0,ngauss)

	call calc_dfg(ue,nnode,ndofel,mcrd,dfgrd1_st,dvol0_st,xjac_st,
     &     shp0_st,xjbar,0,ngauss)

	call makegrad(xbar_st,xbtr_st,xg_st,shptau_st,dvoltau_st,
     & nnode,mcrd,ndofel,ngauss) 	

	f1(1:ndofel)=0.d0
	
	nstatv=nsvars
	do ipt=1,ngauss

	   dfgrd1(1:3,1:3)=dfgrd1_st(1:3,1:3,ipt)
	   dfgrd0(1:3,1:3)=dfgrd0_st(1:3,1:3,ipt)
	   
	   kstep=1
	   
	   istv_start=nsvars*(ipt-1)
	   do ista=1,nstatv
	      statev(ista)=svars(istv_start+ista)
	   enddo

	   call umat(stress,statev,ddsdde,sse,spd,scd,
     & 		 rpl,ddsddt,drplde,drpldt,
     &           stran,dstran,time,dtime,temp,dtemp,predef1,dpred,
     &           cmname,ndi,nshr,ntens,nstatv,props,nprops,coords_ipt,
     &           drot,pnewdt,celent,dfgrd0,dfgrd1,jelem,ipt,layer,
     &           kspt,kstep,kinc,euler,gradient_GND)

	   do ista=1,nstatv
	      svars(istv_start+ista)=statev(ista)
	   enddo


	   if(iconv1.eq.1)return

	   
	   xbar(1:nstr,1:ndofel)=xbar_st(1:nstr,1:ndofel,ipt)*
     &  dvoltau_st(ipt)
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

	      xmtemp1(1:mcrd**2,1:ndofel)=matmul(stmat(1:mcrd**2,
     & 1:mcrd**2) ,xg_st(1:mcrd**2,1:ndofel,ipt))*dvoltau_st(ipt)
	      ske2(1:ndofel,1:ndofel)=matmul(transpose(
     &  xg_st(1:mcrd**2,1:ndofel,ipt)),xmtemp1(1:mcrd**2,1:ndofel))	    

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
		 
	      do i=1,ndofel
		 do j=1,ndofel
		    do k=1,nstr
		       skt(i,j)=skt(i,j)+xbar_st(k,i,ipt)*ct(k,j)*
     &	       dvoltau_st(ipt)
		    enddo
		 enddo
	      enddo

	     

	      crot(1:nstr,1:nstr)=matmul(ddsdde(1:nstr,1:nstr)
     &  ,drotm(1:nstr,1:nstr))    
	      crot(1:nstr,1:nstr)=(crot(1:nstr,1:nstr)-
     &	 stmat1(1:nstr,1:nstr))*dvoltau_st(ipt)           

	      xtm(1:nstr,1:ndofel)=matmul(crot(1:nstr,1:nstr),
     &   xbar_st(1:nstr,1:ndofel,ipt))	      
	      ske3(1:ndofel,1:ndofel)=matmul(transpose(
     &	      xbar_st(1:nstr,1:ndofel,ipt)),xtm(1:nstr,1:ndofel))
	      ske(1:ndofel,1:ndofel)=ske(1:ndofel,1:ndofel)+
     &     ske1(1:ndofel,1:ndofel)+ske3(1:ndofel,1:ndofel)
c	      ske(1:ndofel,1:ndofel)=ske(1:ndofel,1:ndofel)+
c     &     ske1(1:ndofel,1:ndofel)+ske3(1:ndofel,1:ndofel)+ske2(1:ndofel,1:ndofel)
	   endif
	enddo
c	ske(1:ndofel,1:ndofel)=skt(1:ndofel,1:ndofel)


	if(ibelem(jelem).ne.0)then
	   ifs=ibelem(jelem)
           ist=4*(ifs-1)
	   ibnode(1:4)=iface_brick(ist+1:ist+4)

	   do i=1,4
	      xel_surf(3*i-2:3*i)=xe_tau(1:3,ibnode(i))

	   enddo
	   
	   xcg(1:3)=1.d0/8.d0*sum(xe_tau(1:3,:),2)
	   uv1(1:3)=xcg-xe_tau(1:3,ibnode(1))
	   uv2(1:3)=xe_tau(1:3,ibnode(4))-xe_tau(1:3,ibnode(1))
	   uv3(1:3)=xe_tau(1:3,ibnode(2))-xe_tau(1:3,ibnode(1))
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
	   fload=0
	   ske_load=0.d0
	   do ipt1=1,4
	      call surf_prop(reshape(xel_surf,(/3,4/)),ipt1,da,
     &   shp_surf,shp_surf_xi,shp_surf_eta)
	      fload=fload+matmul(transpose(shp_surf),da)
	      dx_xi=matmul(shp_surf_xi,xel_surf)
	      dx_eta=matmul(shp_surf_eta,xel_surf)
	      do i=1,12
		 call cross(tempv1,dx_xi,shp_surf_eta(1:3,i))
		 call cross(tempv2,dx_eta,shp_surf_xi(1:3,i))
		 tempv=tempv1-tempv2
		 ske_load(:,i)=ske_load(:,i)+matmul(transpose(shp_surf),
     &		 tempv)
	      enddo
	   enddo

	      

	   do i=1,4
	      in=ibnode(i)
	      f1(3*in-2:3*in)=f1(3*in-2:3*in)-pr*fload(3*i-2:3*i)
	      do j=1,4
		 jn=ibnode(j)
		 ske(3*in-2:3*in,3*jn-2:3*jn)=ske(3*in-2:3*in,3*jn-2:3*jn)
     &                       -ske_load(3*i-2:3*i,3*j-2:3*j)*pr
	      enddo
	   enddo
        endif

	return
	end


	subroutine makegrad(xb_st,xbtr_st,xg_st,shptau_st,
     &  dvoltau_st,nnode,mcrd,ndofel,ngauss)	
	implicit real*8(a-h,o-z)
	dimension xb_st(2*mcrd,mcrd*nnode,8),dvoltau_st(8),
     &    nshr(3),xb_av(mcrd*nnode),xbtr_st(mcrd*nnode,8),
     &    shptau_st(3,8,8), xg_st(mcrd**2,mcrd*nnode,8),nshr0(3)

	
	nshr(1)=2
	nshr(2)=3
	nshr(3)=3

	nshr0(1)=1
	nshr0(2)=1
	nshr0(3)=2
	   	
	if(mcrd.eq.3)then
	   nstr=6
	endif
	if(mcrd.eq.2)then
	   nstr=3
	endif
	xb_av(1:ndofel)=0.d0

	do ipt=1,ngauss
	   do i=1,nstr
	      do j=1,ndofel
		 xb_st(i,j,ipt)=0.d0
		 jtemp1=(j-1)/mcrd
		 jtemp2=j-jtemp1*mcrd
		 if(i.eq.jtemp2.and.i.le.mcrd)xb_st(i,j,ipt)=
     &	 shptau_st(i,jtemp1+1,ipt)
		 if(i.gt.mcrd)then
		    is=i-mcrd
		    ishr=nshr(is)
		    is=nshr0(is)
	 if(jtemp2.eq.is)xb_st(i,j,ipt)=shptau_st(ishr,jtemp1+1,ipt)
	 if(jtemp2.eq.ishr)xb_st(i,j,ipt)=shptau_st(is,jtemp1+1,ipt)
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
		 if(jtemp2.eq.itemp1)xg_st(i,j,ipt)=    
     &		 shptau_st(itemp2,jtemp1,ipt)
		 if(itemp2.eq.1)xg_st(i,j,ipt)=xg_st(i,j,ipt)-1.d0/3.d0*
     &		 xbtr_st(j,ipt)
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


	subroutine calc_dfgrd(ipt,u,mcrd,nnode,shp_st,dfgrd)
	implicit real*8 (a-h,o-z)
	dimension u(mcrd,nnode),dfgrd(3,3),shp_st(3,8,nnode)

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


	subroutine init_gauss(mcrd)
	implicit real*8 (a-h,o-z)
	dimension xii(8),etai(8),chii(8),xi(8),eta(8),chi(8)
	common/gauss/xii,etai,chii,xi,eta,chi
	os3=1.d0/dsqrt(3.d0)
	xii(1)=-1.d0
	xii(2)=1.d0
	xii(3)=1.d0
	xii(4)=-1.d0
	if(mcrd.eq.3)then
	   xii(5)=-1.d0
	   xii(6)=1.d0
	   xii(7)=1.d0
	   xii(8)=-1.d0
	endif
	
	etai(1)=-1.d0
	etai(2)=-1.d0
	etai(3)=1.d0
	etai(4)=1.d0
	if(mcrd.eq.3)then
	   etai(5)=-1.d0
	   etai(6)=-1.d0
	   etai(7)=1.d0
	   etai(8)=1.d0
	endif
	if(mcrd.eq.3)then
	   chii(1)=-1.d0
	   chii(2)=-1.d0
	   chii(3)=-1.d0
	   chii(4)=-1.d0
	   chii(5)=1.d0
	   chii(6)=1.d0
	   chii(7)=1.d0
	   chii(8)=1.d0
	endif
	xi(1)=-os3
	xi(2)=os3
	xi(3)=os3
	xi(4)=-os3
	if(mcrd.eq.3)then
	   xi(5)=-os3
	   xi(6)=os3
	   xi(7)=os3
	   xi(8)=-os3
	endif
	
	eta(1)=-os3
	eta(2)=-os3
	eta(3)=os3
	eta(4)=os3
	if(mcrd.eq.3)then
	   eta(5)=-os3
	   eta(6)=-os3
	   eta(7)=os3
	   eta(8)=os3
	endif
	
	if(mcrd.eq.3)then
	   chi(1)=-os3
	   chi(2)=-os3
	   chi(3)=-os3
	   chi(4)=-os3
	   chi(5)=os3
	   chi(6)=os3
	   chi(7)=os3
	   chi(8)=os3
	endif	
	return
	end


	subroutine shape2(xel,ipt,shp,shpf,dvol,xc)
	implicit real*8 (a-h,o-z)
	dimension xii(8),etai(8),chii(8),xi(8),eta(8),chi(8),
     &   shp1(3,8),shp(3,8),xjac(2,2),xjacinv(2,2),xel(2,4),xc(2),
     &   shpf(8)
	common/gauss/xii,etai,chii,xi,eta,chi
	oe=1.d0/4.d0
	os3=1.d0/dsqrt(3.d0)

	xip=xi(ipt)
	etap=eta(ipt)

	do i=1,2
	   do j=1,2
	      xjac(i,j)=0.d0
	   enddo
	enddo
	
	do inode=1,4
	   shp1(1,inode)=oe*xii(inode)*(1.d0+etai(inode)*etap)
	   shp1(2,inode)=oe*etai(inode)*(1.d0+xii(inode)*xip)
	 shpf(inode)=oe*(1.d0+xii(inode)*xip)*(1.d0+etai(inode)*etap)
	enddo
	
	do i=1,2
	   xc(i)=0.d0
	enddo

	do inode=1,4
	   do j=1,2
	      xc(j)=xc(j)+shpf(inode)*xel(j,inode)
	   enddo
	enddo

	
	do i=1,2
	   do j=1,2
	      xjac(i,j)=0.d0
	      do inode=1,4
		 xjac(i,j)=xjac(i,j)+shp1(i,inode)*xel(j,inode)
	      enddo
	   enddo
	enddo


	call matinv2_uel(xjac,xjacinv,dvol)

	do inode=1,4
	   do i=1,2
	      shp(i,inode)=0.d0
	      do j=1,2
		 shp(i,inode)=shp(i,inode)+xjacinv(i,j)*shp1(j,inode)
	      enddo
	   enddo
	enddo
	return
	end






      

      subroutine matinv2_uel(a,ai,det)
      implicit real*8(a-h,o-z)
      dimension a(2,2),ai(2,2)
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)

      ai(1,1)=a(2,2)/det
      ai(2,1)=-a(2,1)/det
      ai(1,2)=-a(1,2)/det
      ai(2,2)=a(1,1)/det
      end


      subroutine shape3(xel,ipt,shp,shpf,dvol,xc)
      implicit real*8 (a-h,o-z)
      dimension xii(8),etai(8),chii(8),xi(8),eta(8),chi(8),
     &     shp1(3,8),shp(3,8),xjac(3,3),xjacinv(3,3),xel(3,8),xc(3),
     &     shpf(8)
      common/gauss/xii,etai,chii,xi,eta,chi

      oe=1.d0/8.d0
      os3=1.d0/dsqrt(3.d0)

      xip=xi(ipt)
      etap=eta(ipt)
      chip=chi(ipt)

      do inode=1,8
         do i=1,3
            do j=1,3
               xjac(i,j)=0.d0
            enddo
         enddo
         
         shp1(1,inode)=oe*xii(inode)*(1.d0+etai(inode)*etap)*(1.d0
     $        +chii(inode)*chip)
         shp1(2,inode)=oe*etai(inode)*(1.d0+xii(inode)*xip)*(1.d0
     $        +chii(inode)*chip)
         shp1(3,inode)=oe*chii(inode)*(1.d0+xii(inode)*xip)*(1.d0
     $        +etai(inode)*etap)
       shpf(inode)=oe*(1.d0+xii(inode)*xip)*(1.d0+etai(inode)*etap)*
     &        (1.d0+chii(inode)*chip)
      enddo
      
      do i=1,3
         xc(i)=0.d0
      enddo

      do inode=1,8
         do j=1,3
            xc(j)=xc(j)+shpf(inode)*xel(j,inode)
         enddo
      enddo
      do i=1,3
         do j=1,3
            xjac(i,j)=0.d0
            do inode=1,8
               xjac(i,j)=xjac(i,j)+shp1(i,inode)*xel(j,inode)
            enddo
         enddo
      enddo
      call matinv3_uel(xjac,xjacinv,dvol)
      do inode=1,8
         do i=1,3
            shp(i,inode)=0.d0
            do j=1,3
               shp(i,inode)=shp(i,inode)+xjacinv(i,j)*shp1(j,inode)
            enddo
         enddo
      enddo
      return
      end









	subroutine calc_dfg(ue,nnode,ndofel,mcrd,dfgrd_st,dvol_st,
     &	xjac_st,shp_st,xjbar,ibar,ngauss)
	implicit real*8 (a-h,o-z)
	dimension dfgrd0(3,3),dfgrd1(3,3),ue(mcrd,nnode),
     &   xjac_st(nnode),dfgrd_st(3,3,8),dfgrd_inv(3,3),
     &   dvol_st(8),shp_st(3,8,8),dfgrd(3,3)

	xjbar=0.d0
	tv=0.d0
	do ipt=1,ngauss
	   call calc_dfgrd(ipt,ue,mcrd,nnode,shp_st,dfgrd)
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





	
	
	
	

		subroutine surf_prop(xel,ipt,da,shp_surf,shp_surf_xi,
     &		shp_surf_eta)
	implicit real*8 (a-h,o-z)

	dimension xii(8),etai(8),chii(8),xi(8),eta(8),chi(8),
     &     shp1(3,8),shp(3,8),shpf(4),xel(3,4),da(3),xjac(2,3),
     &     dx_xi(3),dx_eta(3),shp_surf(3,12),shp_surf_xi(3,12),
     &     shp_surf_eta(3,12)

	common/gauss/xii,etai,chii,xi,eta,chi
	oe=1.d0/4.d0
	os3=1.d0/dsqrt(3.d0)

	xip=xi(ipt)
	etap=eta(ipt)

	xjac=0.d0
	
	do inode=1,4
	   shp1(1,inode)=oe*xii(inode)*(1.d0+etai(inode)*etap)
	   shp1(2,inode)=oe*etai(inode)*(1.d0+xii(inode)*xip)
	shpf(inode)=oe*(1.d0+xii(inode)*xip)*(1.d0+etai(inode)*etap)
	enddo
	
	shp_surf=0
	shp_surf(1,1:12:3)=shpf(1:4)
	shp_surf(2,2:12:3)=shpf(1:4)
	shp_surf(3,3:12:3)=shpf(1:4)


	shp_surf_xi=0
	shp_surf_xi(1,1:12:3)=shp1(1,1:4)
	shp_surf_xi(2,2:12:3)=shp1(1,1:4)
	shp_surf_xi(3,3:12:3)=shp1(1,1:4)
	
	shp_surf_eta=0
	shp_surf_eta(1,1:12:3)=shp1(2,1:4)
	shp_surf_eta(2,2:12:3)=shp1(2,1:4)
	shp_surf_eta(3,3:12:3)=shp1(2,1:4)

	
	do i=1,2
	   do j=1,3
	      xjac(i,j)=0.d0
	      do inode=1,4
		 xjac(i,j)=xjac(i,j)+shp1(i,inode)*xel(j,inode)
	      enddo
	   enddo
	enddo

	dx_xi=xjac(1,1:3)
	dx_eta=xjac(2,1:3)

	call cross(da,dx_xi,dx_eta)
	

	return
	end
