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
integer :: inc
integer :: indc
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: t
double precision :: d
double precision :: sn
integer :: lgsnlt
double precision :: vol
double precision :: cnds
double precision :: di6
double precision :: di7
double precision :: dj6
double precision :: dj7
double precision :: dk6
double precision :: dk7
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
integer :: kdir
integer :: m
integer :: n
integer :: n0c
integer :: n1
integer :: nci
integer :: ncj
integer :: nck
integer :: nid
integer :: nijd
integer :: ninc
integer :: njd
double precision :: um
double precision :: vm
double precision :: vn
double precision :: wm
!
!-----------------------------------------------------------------------
!
      character(len=7 ) :: equat
      dimension t(ip11,ip60),d(ip11,ip60)
      dimension vol(ip11)
      dimension sn(lgsnlt,nind,ndir)
!
      indc(i,j,k)=n0c+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
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
!!DEC$ IVDEP
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
!!DEC$ IVDEP
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
!!DEC$ IVDEP
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
      end subroutine
end module
