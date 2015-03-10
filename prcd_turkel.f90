      subroutine prcd_turkel( &
                 lm,u,v, &
                 sn,lgsnlt, &
                 temp,cson)
!
!***********************************************************************
!
!_DA  DATE_C : avril 2003 - Eric Goncalves / LEGI
!
!     ACT
!       Preconditionnement basse vitesse de Turkel.
!       Calcul totalement explicite.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use proprieteflu
      use schemanum
      use definition
implicit none
integer :: indc
integer :: i
integer :: j
integer :: k
integer :: lm
double precision :: u
double precision :: v
double precision :: sn
integer :: lgsnlt
double precision :: temp
double precision :: cson
double precision :: a2
double precision :: beta2
double precision :: cndsi
double precision :: cndsj
double precision :: cndsk
integer :: i1
integer :: i1m1
integer :: i1p1
integer :: i2
integer :: i2m1
integer :: ind1
integer :: ind2
integer :: j1
integer :: j1m1
integer :: j1p1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k1m1
integer :: k1p1
integer :: k2
integer :: k2m1
integer :: m
integer :: n
integer :: n0c
integer :: nid
integer :: nijd
integer :: njd
double precision :: q2
double precision :: qinf
double precision :: uu
double precision :: vv
double precision :: ww
!
!-----------------------------------------------------------------------
!
      real get,ge,gd
      real p11,p12,p13,p14,p15, &
           p21,p22,p23,p24,p25, &
           p31,p32,p33,p34,p35, &
           p41,p42,p43,p44,p45, &
           p51,p52,p53,p54,p55
      dimension u(ip11,ip60),v(ip11,ip60)
      dimension sn(lgsnlt,nind,ndir)
      dimension temp(ip11),cson(ip11)
!
      indc(i,j,k)=n0c+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
!
      n0c=npc(lm)
      i1=ii1(lm)
      i2=ii2(lm)
      j1=jj1(lm)
      j2=jj2(lm)
      k1=kk1(lm)
      k2=kk2(lm)
!
      nid = id2(lm)-id1(lm)+1
      njd = jd2(lm)-jd1(lm)+1
      nijd = nid*njd
!
      i1p1=i1+1
      j1p1=j1+1
      k1p1=k1+1
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
!
      qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
!
      ind1 = indc(i1  ,j1,  k1  )
      ind2 = indc(i2m1,j2m1,k2m1)
      do n=ind1,ind2
       m=n-n0c
       cndsi=sqrt(sn(m,1,1)*sn(m,1,1)+ &
                  sn(m,1,2)*sn(m,1,2)+ &
                  sn(m,1,3)*sn(m,1,3))
       cndsj=sqrt(sn(m,2,1)*sn(m,2,1)+ &
                  sn(m,2,2)*sn(m,2,2)+ &
                  sn(m,2,3)*sn(m,2,3))
       cndsk=sqrt(sn(m,3,1)*sn(m,3,1)+ &
                  sn(m,3,2)*sn(m,3,2)+ &
                  sn(m,3,3)*sn(m,3,3))
       uu=(v(n,2)*sn(m,1,1)+v(n,3)*sn(m,1,2)+v(n,4)*sn(m,1,3))/(v(n,1)*cndsi)
       vv=(v(n,2)*sn(m,2,1)+v(n,3)*sn(m,2,2)+v(n,4)*sn(m,2,3))/(v(n,1)*cndsj)
       ww=(v(n,2)*sn(m,3,1)+v(n,3)*sn(m,3,2)+v(n,4)*sn(m,3,3))/(v(n,1)*cndsk)
       q2=uu**2+vv**2+ww**2
       a2=cson(n)**2
       beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
       get=v(n,5)/v(n,1)                              !energie totale
       ge=get-0.5*q2
       gd=(beta2-1.)/ge
!      matrice de preconditionnement Pc
       p11 = 1. + gd*(0.5*q2)
       p12 = -gd*uu
       p13 = -gd*vv
       p14 = -gd*ww
       p15 = gd
       p21 = gd*(0.5*q2)*uu
       p22 = 1. - gd*uu**2
       p23 = -gd*uu*vv
       p24 = -gd*uu*ww
       p25 = gd*uu
       p31 = gd*(0.5*q2)*vv
       p32 = -gd*uu*vv
       p33 = 1. - gd*vv**2
       p34 = -gd*vv*ww
       p35 = gd*vv
       p41 = gd*(0.5*q2)*ww
       p42 = -gd*uu*ww
       p43 = -gd*vv*ww
       p44 = 1. - gd*ww**2
       p45 = gd*ww
       p51 = gd*(0.5*q2)*get
       p52 = -gd*uu*get
       p53 = -gd*vv*get
       p54 = -gd*ww*get
       p55 = 1. + gd*get
!
       u(n,1)=p11*u(n,1)+p12*u(n,2)+p13*u(n,3)+p14*u(n,4)+p15*u(n,5)
       u(n,2)=p21*u(n,1)+p22*u(n,2)+p23*u(n,3)+p24*u(n,4)+p25*u(n,5)
       u(n,3)=p31*u(n,1)+p32*u(n,2)+p33*u(n,3)+p34*u(n,4)+p35*u(n,5)
       u(n,4)=p41*u(n,1)+p42*u(n,2)+p43*u(n,3)+p44*u(n,4)+p45*u(n,5)
       u(n,5)=p51*u(n,1)+p52*u(n,2)+p53*u(n,3)+p54*u(n,4)+p55*u(n,5)
      enddo
!
      return
      end
