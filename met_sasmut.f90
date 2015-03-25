module mod_met_sasmut
  implicit none
contains
  subroutine met_sasmut( &
       l, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       v,mu,mut,dist, &
       sn,vol,t, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_EG         decembre 2001 - Eric Goncalves / SINUMEF
!
!     ACT
!_A   modele de Spalart Allmaras
!_A   calcul du mu_t avec correction SST de Menter
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!
!     OUT
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use mod_teq_gradv
    implicit none
  integer          ::    i,  i1,  i2,i2m1,ind1
  integer          :: ind2,   j,  j1,  j2,j2m1
  integer          ::    k,  k1,  k2,k2m1,   l
  integer          ::    m,   n,  n0, nid,nijd
  integer          ::  njd
  double precision ::            a1,  cmui1(ip21),  cmui2(ip21),  cmuj1(ip21),  cmuj2(ip21)
  double precision ::   cmuk1(ip21),  cmuk2(ip21),         cmut,        coef1,        coef2
  double precision ::            ct,         cv13,   dist(ip12),          ds3,   dvxx(ip00)
  double precision ::    dvxy(ip00),   dvxz(ip00),   dvyx(ip00),   dvyy(ip00),   dvyz(ip00)
  double precision ::    dvzx(ip00),   dvzy(ip00),   dvzz(ip00),        exp2x,           f2
  double precision ::           fv1,     mu(ip12),    mut(ip12),           rk,         rota
  double precision :: sn(ip31*ndir),      t(ip00), v(ip11,ip60),    vol(ip11),         xkhi
  double precision ::         xkhi3,         zeta
!
!-----------------------------------------------------------------------
!
!

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
!c    correction SST
    ds3=2./3.
    cmu=0.09
    a1=sqrt(cmu)
    cv13=cv1**3
    ct=1./(0.41*cmu**0.25)
!     Calcul du gradient de la vitesse. Les tableaux ont ete reutilises
!     dans la phase implicite.
    call teq_gradv( &
         l, &
         sn, &
         vol,v, &
         t , &
         dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=ind(i1,j,k)
          ind2=ind(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0
             xkhi=v(n,6)/mu(n)
             xkhi3=xkhi**3
             fv1=xkhi3/(xkhi3+cv13)
             cmut=v(n,6)*fv1
!           toxx=ds3*cmut*(2*dvxx(m)-dvyy(m)-dvzz(m))
!           toyy=ds3*cmut*(2*dvyy(m)-dvxx(m)-dvzz(m))
!           tozz=ds3*cmut*(2*dvzz(m)-dvxx(m)-dvyy(m))
!           rk=0.5*v(n,1)*(toxx**2+toyy**2+tozz**2)
             rk=(cmut*ct/dist(n))**2/v(n,1)
             rk=max(rk,rokinf)
!
!           if(rk.lt.1.e-06) then
!             mut(n)=cmut
!           else
             coef1=2.*v(n,6)/(cmu*sqrt(rk*v(n,1))*dist(n))
             coef2=500.*mu(n)*v(n,6)/(rk*v(n,1)*dist(n)**2)
             rota = sqrt((dvzy(m)-dvyz(m))**2 &
                  +(dvxz(m)-dvzx(m))**2 &
                  +(dvyx(m)-dvxy(m))**2)
             zeta=max(coef1,coef2)
             exp2x=exp(min(2.*zeta**2,25.))
             f2=(exp2x-1.)/(exp2x+1.)
!            mut(n)=min(cmut,rk*a1/rota)
             mut(n)=cmut/max(1.,f2*rota*cmut/rk*a1)
!           endif
          enddo
       enddo
    enddo
!
    return
  contains
    function    ind(i,j,k)
      implicit none
  integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine met_sasmut
end module mod_met_sasmut
