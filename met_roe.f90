module mod_met_roe
  implicit none
contains
  subroutine met_roe( &
       l,t,d, &
       equat, &
       sn,lgsnlt, &
       vol)
!
!***********************************************************************
!
!_A   Auteur:       Eric GONCALVES / SINUMEF
!
!     Schema decentre de Roe - ordre 1.
!     Stockage de la dissipation dans le tableau d.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::      i,    i1,  i1m1,  i1p1,    i2
    integer          ::   i2m1,    id,  ind1,  ind2,     j
    integer          ::     j1,  j1m1,  j1p1,    j2,  j2m1
    integer          ::     jd,     k,    k1,  k1m1,  k1p1
    integer          ::     k2,  k2m1,    kd,  kdir,     l
    integer          :: lgsnlt,     m,     n,   n0c,    n1
    integer          ::    nci,   ncj,   nck,   nid,  nijd
    integer          ::   ninc,   njd
    double precision ::                 cnds,        d(ip11,ip60),                 di6,                 di7,                 dj6
    double precision ::                  dj7,                 dk6,                 dk7,sn(lgsnlt,nind,ndir),        t(ip11,ip60)
    double precision ::                   um,                  vm,                  vn,           vol(ip11),                  wm
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
!


!
    n0c=npc(l)
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    nid  = id2(l)-id1(l)+1
    njd  = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    i1p1 = i1+1
    j1p1 = j1+1
    k1p1 = k1+1
    i2m1 = i2-1
    j2m1 = j2-1
    k2m1 = k2-1
    i1m1 = i1-1
    j1m1 = j1-1
    k1m1 = k1-1
!
    nci   = inc(1,0,0)
    ncj   = inc(0,1,0)
    nck   = inc(0,0,1)
!
! ------------------------------------------------------------
!
    if(equat(3:5).eq.'2dk') then
       ind1=indc(i1m1,j1m1,k1  )
       ind2=indc(i2  ,j2  ,k2m1)
    elseif(equat(3:4).eq.'3d') then
       ind1=indc(i1m1,j1m1,k1m1)
       ind2=indc(i2  ,j2  ,k2  )
    endif
    do n=ind1,ind2
       d(n,6)=0.
       d(n,7)=0.
    enddo
!
!*****************************************************************************
!    calcul du flux numerique par direction
!*****************************************************************************
!
!-----direction i-------------------------------------------------------
!
    kdir=1
    ninc=nci
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2,j,k)
!!!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             n1=n-ninc
!        vecteur normal unitaire a la face consideree
             cnds=sqrt(sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3))
!        calcul d'un etat moyen
             um=0.5*(t(n,2)/t(n,1)+t(n1,2)/t(n1,1))
             vm=0.5*(t(n,3)/t(n,1)+t(n1,3)/t(n1,1))
             wm=0.5*(t(n,4)/t(n,1)+t(n1,4)/t(n1,1))
!        valeur propre
             vn=abs(um*sn(m,kdir,1)+vm*sn(m,kdir,2)+wm*sn(m,kdir,3))/cnds
             di6=vn*(t(n,6)-t(n1,6))
             di7=vn*(t(n,7)-t(n1,7))
             d(n,6)=d(n,6)-di6
             d(n,7)=d(n,7)-di7
             d(n-ninc,6)=d(n-ninc,6)+di6
             d(n-ninc,7)=d(n-ninc,7)+di7
          enddo
       enddo
    enddo
!
!------direction j----------------------------------------------
!
    kdir=2
    ninc=ncj
!
    do k=k1,k2m1
       do j=j1,j2
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
!!!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             n1=n-ninc
!        vecteur normal unitaire a la face consideree
             cnds=sqrt(sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3))
!        calcul d'un etat moyen
             um=0.5*(t(n,2)/t(n,1)+t(n1,2)/t(n1,1))
             vm=0.5*(t(n,3)/t(n,1)+t(n1,3)/t(n1,1))
             wm=0.5*(t(n,4)/t(n,1)+t(n1,4)/t(n1,1))
!        valeur propre
             vn=abs(um*sn(m,kdir,1)+vm*sn(m,kdir,2)+wm*sn(m,kdir,3))/cnds
             dj6=vn*(t(n,6)-t(n1,6))
             dj7=vn*(t(n,7)-t(n1,7))
             d(n,6)=d(n,6)-dj6
             d(n,7)=d(n,7)-dj7
             d(n-ninc,6)=d(n-ninc,6)+dj6
             d(n-ninc,7)=d(n-ninc,7)+dj7
          enddo
       enddo
    enddo
!
!------direction k-------------------------------------------------------
!
    if(equat(3:4).eq.'3d') then
       kdir=3
       ninc=nck
!
       do k=k1,k2
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
!!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                n1=n-ninc
!        vecteur normal unitaire a la face consideree
                cnds=sqrt(sn(m,kdir,1)*sn(m,kdir,1)+ &
                     sn(m,kdir,2)*sn(m,kdir,2)+ &
                     sn(m,kdir,3)*sn(m,kdir,3))
!        calcul d'un etat moyen
                um=0.5*(t(n,2)/t(n,1)+t(n1,2)/t(n1,1))
                vm=0.5*(t(n,3)/t(n,1)+t(n1,3)/t(n1,1))
                wm=0.5*(t(n,4)/t(n,1)+t(n1,4)/t(n1,1))
!        valeur propre
                vn=abs(um*sn(m,kdir,1)+vm*sn(m,kdir,2)+wm*sn(m,kdir,3))/cnds
                dk6=vn*(t(n,6)-t(n1,6))
                dk7=vn*(t(n,7)-t(n1,7))
                d(n,6)=d(n,6)-dk6
                d(n,7)=d(n,7)-dk7
                d(n-ninc,6)=d(n-ninc,6)+dk6
                d(n-ninc,7)=d(n-ninc,7)+dk7
             enddo
          enddo
       enddo
!
    endif
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine met_roe
end module mod_met_roe
