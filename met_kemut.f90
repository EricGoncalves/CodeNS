      subroutine met_kemut( &
                 l, &
                 sn,vol,t, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 dist,v,mu,mut, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT
!        Calculs de la viscosite turbulente mut modele k-eps
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
      use proprieteflu
      use chainecarac
implicit none
integer :: inc
integer :: ind
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: sn
double precision :: vol
double precision :: t
double precision :: dvxx
double precision :: dvxy
double precision :: dvxz
double precision :: dvyx
double precision :: dvyy
double precision :: dvyz
double precision :: dvzx
double precision :: dvzy
double precision :: dvzz
double precision :: dist
double precision :: v
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
integer :: i1
integer :: i2
integer :: i2m1
integer :: j1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k2
integer :: k2m1
integer :: m
integer :: n
integer :: n0
integer :: nci
integer :: ncj
integer :: nck
integer :: nid
integer :: nijd
integer :: njd
!
!-----------------------------------------------------------------------
!
      real mu,mut,mut0,retur,fmu,a1,coef1,coef2,rota,zeta,exp2x,f2
      dimension mu(ip12),mut(ip12)
      dimension v(ip11,ip60)
      dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
                dvyx(ip00),dvyy(ip00),dvyz(ip00), &
                dvzx(ip00),dvzy(ip00),dvzz(ip00),t(ip00)
      dimension sn(ip31*ndir),vol(ip11),dist(ip12)
      dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
      ind(i,j,k)=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
!
      n0=npc(l)
      i1=ii1(l)
      i2=ii2(l)
      j1=jj1(l)
      j2=jj2(l)
      k1=kk1(l)
      k2=kk2(l)
!
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd= nid*njd
!
      nci = inc(1,0,0)
      ncj = inc(0,1,0)
      nck = inc(0,0,1)
!
      if(equatt(1:3).eq.'2JL') then
       if(kesst.eq.0) then  ! modele de base
        do k=k1,k2m1
         do j=j1,j2m1
          n=ind(i1-1,j,k)
          do i=i1,i2m1
           n=n+nci
           retur=(v(n,6)**2)/(v(n,7)*mu(n))
           mut0=cmu*retur*mu(n)
!          lois de paroi
           if((lparoi.eq.1).or.(lparoi.eq.2)) then
            fmu=1.
           else
            fmu=exp(-2.5/(1.+retur/50.))
           endif
           mut(n)=fmu*mut0
          enddo
         enddo
        enddo
!
!     modele avec correction SST
       elseif(kesst.eq.1) then
        a1=sqrt(cmu)
!       Calcul du gradient de la vitesse. Les tableaux ont ete reutilises
!       dans la phase implicite.
        call teq_gradv( &
             l, &
             sn, &
             vol,v, &
             t, &
             dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
             cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
        do k=k1,k2m1
         do j=j1,j2m1
          n=ind(i1-1,j,k)
          do i=i1,i2m1
           n=n+nci
           m=n-n0
           retur=(v(n,6)**2)/(v(n,7)*mu(n))
           mut0=cmu*retur*mu(n)
!          lois de paroi
           if((lparoi.eq.1).or.(lparoi.eq.2)) then
             fmu=1.
            else
             fmu=exp(-2.5/(1.+retur/50.))
           endif
           coef1=2.*v(n,6)**1.5/(v(n,7)*sqrt(v(n,1))*dist(n))
           coef2=500.*mu(n)*cmu*v(n,6)/(v(n,7)*dist(n)**2*v(n,1))
           rota = sqrt((dvzy(m)-dvyz(m))**2 &
                      +(dvxz(m)-dvzx(m))**2 &
                      +(dvyx(m)-dvxy(m))**2)
           zeta=max(coef1,coef2)
           exp2x=exp(min(2.*zeta**2,25.))
           f2=(exp2x-1.)/(exp2x+1.)
!           mut(n)=v(n,6)/max(v(n,7)/(cmu*fmu*v(n,6)),rota*f2/a1)
           mut(n)=fmu*mut0/max(1.,rota*f2*v(n,6)*a1*fmu/v(n,7))
          enddo
         enddo
        enddo
      endif
!
      elseif(equatt(1:3).eq.'2LS') then
       do k=k1,k2m1
        do j=j1,j2m1
         n=ind(i1-1,j,k)
         do i=i1,i2m1
          n=n+nci
          retur=(v(n,6)**2)/(v(n,7)*mu(n))
          mut0=cmu*retur*mu(n)
          if((lparoi.eq.1).or.(lparoi.eq.2)) then
            fmu=1.
           else
            fmu=exp(-3.4/((1.+retur/50.)**2) )
          endif
          mut(n)=fmu*mut0
         enddo
        enddo
       enddo
      endif
!
      return
      end
