module mod_met_komut
  use mod_teq_gradv
  implicit none
contains
  subroutine met_komut( &
       l, &
       sn,vol,t, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       dist,v,mu,mut, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT
!_A   Calcul de la viscosite turbulente a partir de k et omega
!_A   Modele de Wilcox et Menter
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!
!     OUT
!_O    mut        : arg real(ip12      ) ; viscosite turbulente
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use proprieteflu
    implicit none
    integer          ::    i,  i1,  i2,i2m1,  id
    integer          ::    j,  j1,  j2,j2m1,  jd
    integer          ::    k,  k1,  k2,k2m1,  kd
    integer          ::    l,   m,   n,  n0, nci
    integer          ::  nid,nijd, njd
    double precision ::            a1,       betae2,  cmui1(ip21),  cmui2(ip21),  cmuj1(ip21)
    double precision ::   cmuj2(ip21),  cmuk1(ip21),  cmuk2(ip21),        coef1,        coef2
    double precision ::    dist(ip12),   dvxx(ip00),   dvxy(ip00),   dvxz(ip00),   dvyx(ip00)
    double precision ::    dvyy(ip00),   dvyz(ip00),   dvzx(ip00),   dvzy(ip00),   dvzz(ip00)
    double precision ::         exp2x,           f2,     mu(ip12),    mut(ip12),         omeg
    double precision ::          rota,sn(ip31*ndir),      t(ip00), v(ip11,ip60),    vol(ip11)
    double precision ::          zeta
!
!-----------------------------------------------------------------------
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
!
    if(komsst.eq.0) then
!-----------------------------------------------------
!       modele de base Wilcox ou de Menter
!----------------------------------------------------
!
       do k=k1,k2m1
          do j=j1,j2m1
             n=ind(i1-1,j,k)
             do i=i1,i2m1
                n=n+nci
                mut(n)=v(n,1)*v(n,6)/v(n,7)
             enddo
          enddo
       enddo
!
    elseif(komsst.eq.1) then
!---------------------------------------------------
!       Modele Menter SST
!---------------------------------------------------
       betae2=betae/2.
       a1    =sqrt(betae)
!       Calcul du gradient de la vitesse. Les tableaux ont ete reutilises
!       dans la phase implicite.
!
       call teq_gradv( &
            l, &
            sn, &
            vol,v, &
            t , &
            dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
            cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)   !cEG-TOM
!
       do k=k1,k2m1
          do j=j1,j2m1
             n=ind(i1-1,j,k)
             do i=i1,i2m1
                n=n+nci
                m=n-n0
                omeg=v(n,7)/v(n,1)
                coef1=sqrt(v(n,6)/v(n,1))/(betae2*omeg*dist(n))
                coef2=500.*mu(n)/(v(n,7)*dist(n)**2)
                rota = sqrt( (dvzy(m)-dvyz(m))**2 &
                     +(dvxz(m)-dvzx(m))**2 &
                     +(dvyx(m)-dvxy(m))**2)
                zeta=max(coef1,coef2)
                exp2x=exp(min(2.*zeta**2,25.))
                f2=(exp2x-1.)/(exp2x+1)
                mut(n)=v(n,6)/max(omeg,rota*f2/a1)
             enddo
          enddo
       enddo
    endif
!
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
  end subroutine met_komut
end module mod_met_komut
