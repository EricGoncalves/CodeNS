module mod_met_kemut
  implicit none
contains
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
    use mod_teq_gradv
    implicit none
    integer          ::    i,  i1,  i2,i2m1,  id
    integer          ::    j,  j1,  j2,j2m1,  jd
    integer          ::    k,  k1,  k2,k2m1,  kd
    integer          ::    l,   m,   n,  n0, nci
    integer          ::  ncj, nck, nid,nijd, njd
    double precision ::            a1,  cmui1(ip21),  cmui2(ip21),  cmuj1(ip21),  cmuj2(ip21)
    double precision ::   cmuk1(ip21),  cmuk2(ip21),        coef1,        coef2,   dist(ip12)
    double precision ::    dvxx(ip00),   dvxy(ip00),   dvxz(ip00),   dvyx(ip00),   dvyy(ip00)
    double precision ::    dvyz(ip00),   dvzx(ip00),   dvzy(ip00),   dvzz(ip00),        exp2x
    double precision ::            f2,          fmu,     mu(ip12),    mut(ip12),         mut0
    double precision ::         retur,         rota,sn(ip31*ndir),      t(ip00), v(ip11,ip60)
    double precision ::     vol(ip11),         zeta
!$OMP MASTER
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
    nci = inc(1,0,0)
    ncj = inc(0,1,0)
    nck = inc(0,0,1)
!
    if(equatt(1:3).eq.'2JL') then
       if(kesst.eq.0) then  ! modele de base
          do k=k1,k2m1
             do j=j1,j2m1
                n=ind(i1-1,j,k)
!$OMP SIMD
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
!$OMP SIMD
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
!$OMP SIMD
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
!$OMP END MASTER
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine met_kemut
end module mod_met_kemut
