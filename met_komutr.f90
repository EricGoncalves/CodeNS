module mod_met_komutr
  implicit none
contains
  subroutine met_komutr( &
       l, &
       sn,vol,t, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       dist,v,mu,mut, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_H   DATE_C : septembre 2002  - AUTEUR: Eric Goncalves - LEGI
!_H
!     ACT
!_A   Calcul de la viscosite turbulente a partir de k et omega
!_A   Modele de Menter avec conditions de realisabilite de Menter
!_A   avec ou sans correction SST.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!
!     OUT
!_O    mut        : arg real(ip12      ) ; viscosite turbulente
!
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
    double precision ::            a1,        alpha,       betae2,  cmui1(ip21),  cmui2(ip21)
    double precision ::   cmuj1(ip21),  cmuj2(ip21),  cmuk1(ip21),  cmuk2(ip21),        coef1
    double precision ::         coef2,   dist(ip12),   dvxx(ip00),   dvxy(ip00),   dvxz(ip00)
    double precision ::    dvyx(ip00),   dvyy(ip00),   dvyz(ip00),   dvzx(ip00),   dvzy(ip00)
    double precision ::    dvzz(ip00),        exp2x,           f2,     mu(ip12),    mut(ip12)
    double precision ::          omeg,         rcmu,         rmut,         rota,         smut
    double precision :: sn(ip31*ndir),           ss,      t(ip00), v(ip11,ip60),    vol(ip11)
    double precision ::          zeta
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
    rcmu=1./sqrt(0.09)         !correspond a c=0.52
    betae2=betae/2.
    a1=sqrt(betae)
!
!       Calcul du gradient de la vitesse
    call teq_gradv( &
         l, &
         sn, &
         vol,v, &
         t , &
         dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    if(komsst.eq.0) then
!       Modele sans correction SST
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = ind(i1  ,j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do n=ind1,ind2
                m=n-n0
                ss=v(n,1)*sqrt(4.*(dvxx(m)**2+dvyy(m)**2+ &
                     dvzz(m)**2)/3. + (dvzy(m)+dvyz(m))**2 &
                     + (dvxz(m)+dvzx(m))**2 &
                     + (dvyx(m)+dvxy(m))**2)/v(n,7)
                alpha=min(1.,rcmu/ss)
                mut(n)=v(n,1)*alpha*v(n,6)/v(n,7)
             enddo
          enddo
       enddo
!
    else if(komsst.eq.1) then
!       Modele avec correction SST
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = ind(i1  ,j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do n=ind1,ind2
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
                smut=v(n,6)/max(omeg,rota*f2/a1)
                ss=v(n,1)*sqrt(4.*(dvxx(m)**2+dvyy(m)**2+ &
                     dvzz(m)**2)/3. + (dvzy(m)+dvyz(m))**2 &
                     + (dvxz(m)+dvzx(m))**2 &
                     + (dvyx(m)+dvxy(m))**2)/v(n,7)
                alpha=min(1.,rcmu/ss)
                rmut=v(n,1)*alpha*v(n,6)/v(n,7)
                mut(n)=min(smut,rmut)
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
  end subroutine met_komutr
end module mod_met_komutr
