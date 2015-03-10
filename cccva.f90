      subroutine cccva( &
                 l,equat, &
                 v,mut, &
                 vv1,vv2,vv3,vv4,vv5,vv6,vv7,mmut)
!
!***********************************************************************
!
!     ACT
!_A    Transfert de variables aux centres a variables aux noeuds.
!
!_I    l          : arg int              ; numero de domaine
!_I    equat      : arg char             ; type d'equations modelisant l'ecoulement
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!
!     OUT
!_O    vv1        : arg real(ip00      ) ; variable de calcul numero 1
!_O    vv2        : arg real(ip00      ) ; variable de calcul numero 2
!_O    vv3        : arg real(ip00      ) ; variable de calcul numero 3
!_O    vv4        : arg real(ip00      ) ; variable de calcul numero 4
!_O    vv5        : arg real(ip00      ) ; variable de calcul numero 5
!_O    mmut       : arg real(ip00      ) ; viscosite turbulente
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
   use maillage
implicit none
integer :: inc
integer :: ind
integer :: l
double precision :: v
double precision :: vv1
double precision :: vv2
double precision :: vv3
double precision :: vv4
double precision :: vv5
double precision :: vv6
double precision :: vv7
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: i1
integer :: i1m1
integer :: i2
integer :: i2m1
integer :: i2p1
integer :: imax
integer :: imin
integer :: j1
integer :: j1m1
integer :: j2
integer :: j2m1
integer :: j2p1
integer :: jmax
integer :: jmin
integer :: k1
integer :: k1m1
integer :: k2
integer :: k2m1
integer :: k2p1
integer :: kmax
integer :: kmin
integer :: m
integer :: n
integer :: n0
integer :: n1
integer :: n2
integer :: nci
integer :: ncij
integer :: ncijk
integer :: ncik
integer :: ncj
integer :: ncjk
integer :: nck
integer :: nid
integer :: nijd
integer :: njd
!
!-----------------------------------------------------------------------
!
      character(len=7 ) :: equat
      real mmut,mut
      dimension v(ip11,ip60),mut(ip12)
      dimension vv1(ip00),vv2(ip00),vv3(ip00),vv4(ip00),vv5(ip00), &
                vv6(ip00),vv7(ip00),mmut(ip00)
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
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
      i2p1=i2+1
      j2p1=j2+1
      k2p1=k2+1
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd = nid*njd
!
      nci = inc(1,0,0)
      ncj = inc(0,1,0)
      nck = inc(0,0,1)
      ncij = inc(1,1,0)
      ncik = inc(1,0,1)
      ncjk = inc(0,1,1)
      ncijk= inc(1,1,1)
!
      n1=ind(i1m1,j1m1,k1m1)-n0
      n2=ind(i2p1,j2p1,k2p1)-n0
      do m=n1,n2
      vv1(m)=0.
      vv2(m)=0.
      vv3(m)=0.
      vv4(m)=0.
      vv5(m)=0.
      mmut(m)=0.
      enddo
!
      if(equat(6:7).eq.'ke') then
      do m=n1,n2
      vv6(m)=0.
      vv7(m)=0.
      enddo
      end if
!
      imin=i1m1
      imax=i2
      jmin=j1m1
      jmax=j2
      kmin=k1m1
      kmax=k2
!
      if (equat(3:5).eq.'2di') then
      imin=i1
      imax=i2m1
      end if
      if (equat(3:5).eq.'2dj') then
      jmin=j1
      jmax=j2m1
      end if
      if (equat(3:5).eq.'2dk') then
      kmin=k1
      kmax=k2m1
      end if
      if (equat(3:5).eq.'2xk') then
      kmin=k1
      kmax=k2m1
      end if
!
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
      n=ind(i,j,k)
      m=n-n0
!
      vv1(m      )=vv1(m      )+.125*v(n,1)
      vv1(m+nci  )=vv1(m+nci  )+.125*v(n,1)
      vv1(m+ncj  )=vv1(m+ncj  )+.125*v(n,1)
      vv1(m+nck  )=vv1(m+nck  )+.125*v(n,1)
      vv1(m+ncij )=vv1(m+ncij )+.125*v(n,1)
      vv1(m+ncik )=vv1(m+ncik )+.125*v(n,1)
      vv1(m+ncjk )=vv1(m+ncjk )+.125*v(n,1)
      vv1(m+ncijk)=vv1(m+ncijk)+.125*v(n,1)
!
      vv2(m      )=vv2(m      )+.125*v(n,2)
      vv2(m+nci  )=vv2(m+nci  )+.125*v(n,2)
      vv2(m+ncj  )=vv2(m+ncj  )+.125*v(n,2)
      vv2(m+nck  )=vv2(m+nck  )+.125*v(n,2)
      vv2(m+ncij )=vv2(m+ncij )+.125*v(n,2)
      vv2(m+ncik )=vv2(m+ncik )+.125*v(n,2)
      vv2(m+ncjk )=vv2(m+ncjk )+.125*v(n,2)
      vv2(m+ncijk)=vv2(m+ncijk)+.125*v(n,2)
!
      vv3(m      )=vv3(m      )+.125*v(n,3)
      vv3(m+nci  )=vv3(m+nci  )+.125*v(n,3)
      vv3(m+ncj  )=vv3(m+ncj  )+.125*v(n,3)
      vv3(m+nck  )=vv3(m+nck  )+.125*v(n,3)
      vv3(m+ncij )=vv3(m+ncij )+.125*v(n,3)
      vv3(m+ncik )=vv3(m+ncik )+.125*v(n,3)
      vv3(m+ncjk )=vv3(m+ncjk )+.125*v(n,3)
      vv3(m+ncijk)=vv3(m+ncijk)+.125*v(n,3)
!
      vv4(m      )=vv4(m      )+.125*v(n,4)
      vv4(m+nci  )=vv4(m+nci  )+.125*v(n,4)
      vv4(m+ncj  )=vv4(m+ncj  )+.125*v(n,4)
      vv4(m+nck  )=vv4(m+nck  )+.125*v(n,4)
      vv4(m+ncij )=vv4(m+ncij )+.125*v(n,4)
      vv4(m+ncik )=vv4(m+ncik )+.125*v(n,4)
      vv4(m+ncjk )=vv4(m+ncjk )+.125*v(n,4)
      vv4(m+ncijk)=vv4(m+ncijk)+.125*v(n,4)
!
      vv5(m      )=vv5(m      )+.125*v(n,5)
      vv5(m+nci  )=vv5(m+nci  )+.125*v(n,5)
      vv5(m+ncj  )=vv5(m+ncj  )+.125*v(n,5)
      vv5(m+nck  )=vv5(m+nck  )+.125*v(n,5)
      vv5(m+ncij )=vv5(m+ncij )+.125*v(n,5)
      vv5(m+ncik )=vv5(m+ncik )+.125*v(n,5)
      vv5(m+ncjk )=vv5(m+ncjk )+.125*v(n,5)
      vv5(m+ncijk)=vv5(m+ncijk)+.125*v(n,5)
      enddo
      enddo
      enddo
!
      if (equat(1:2).eq.'ns') then
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
      n=ind(i,j,k)
      m=n-n0
      mmut(m      )=mmut(m      )+.125*mut(n)
      mmut(m+nci  )=mmut(m+nci  )+.125*mut(n)
      mmut(m+ncj  )=mmut(m+ncj  )+.125*mut(n)
      mmut(m+nck  )=mmut(m+nck  )+.125*mut(n)
      mmut(m+ncij )=mmut(m+ncij )+.125*mut(n)
      mmut(m+ncik )=mmut(m+ncik )+.125*mut(n)
      mmut(m+ncjk )=mmut(m+ncjk )+.125*mut(n)
      mmut(m+ncijk)=mmut(m+ncijk)+.125*mut(n)
      enddo
      enddo
      enddo
      endif
!
      if(equat(6:7).eq.'ke') then
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
      n=ind(i,j,k)
      m=n-n0
      vv6(m      )=vv6(m      )+.125*v(n,6)
      vv6(m+nci  )=vv6(m+nci  )+.125*v(n,6)
      vv6(m+ncj  )=vv6(m+ncj  )+.125*v(n,6)
      vv6(m+nck  )=vv6(m+nck  )+.125*v(n,6)
      vv6(m+ncij )=vv6(m+ncij )+.125*v(n,6)
      vv6(m+ncik )=vv6(m+ncik )+.125*v(n,6)
      vv6(m+ncjk )=vv6(m+ncjk )+.125*v(n,6)
      vv6(m+ncijk)=vv6(m+ncijk)+.125*v(n,6)
!
      vv7(m      )=vv7(m      )+.125*v(n,7)
      vv7(m+nci  )=vv7(m+nci  )+.125*v(n,7)
      vv7(m+ncj  )=vv7(m+ncj  )+.125*v(n,7)
      vv7(m+nck  )=vv7(m+nck  )+.125*v(n,7)
      vv7(m+ncij )=vv7(m+ncij )+.125*v(n,7)
      vv7(m+ncik )=vv7(m+ncik )+.125*v(n,7)
      vv7(m+ncjk )=vv7(m+ncjk )+.125*v(n,7)
      vv7(m+ncijk)=vv7(m+ncijk)+.125*v(n,7)
      enddo
      enddo
      enddo
      end if
!
      if (equat(3:4).eq.'2d') then
!
      do k=k1,k2
      do j=j1,j2
      do i=i1,i2
      n=ind(i,j,k)
      m=n-n0
      vv1(m)=2*vv1(m)
      vv2(m)=2*vv2(m)
      vv3(m)=2*vv3(m)
      vv4(m)=2*vv4(m)
      vv5(m)=2*vv5(m)
      enddo
      enddo
      enddo
!
      if (equat(1:2).eq.'ns') then
      do k=k1,k2
      do j=j1,j2
      do i=i1,i2
      n=ind(i,j,k)
      m=n-n0
      mmut(m)=2*mmut(m)
      enddo
      enddo
      enddo
      end if
!
      if(equat(6:7).eq.'ke') then
      do k=k1,k2
      do j=j1,j2
      do i=i1,i2
      n=ind(i,j,k)
      m=n-n0
      vv6(m)=2*vv6(m)
      vv7(m)=2*vv7(m)
      enddo
      enddo
      enddo
      end if
!
      end if
!
      return
      end
