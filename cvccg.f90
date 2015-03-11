module mod_cvccg
implicit none
contains
      subroutine cvccg( &
                 l, &
                 x,y,z, &
                 xx,yy,zz)
!
!***********************************************************************
!
!     ACT
!_A    Transfert de coordonnees aux noeuds a coordonnees aux centres,
!_A    y compris pour les mailles fictives.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
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
!_O    xx         : arg real(ip00      ) ; coordonnee sur l'axe x du centre de
!_O                                        la maille
!_O    yy         : arg real(ip00      ) ; coordonnee sur l'axe y du centre de
!_O                                        la maille
!_O    zz         : arg real(ip00      ) ; coordonnee sur l'axe z du centre de
!_O                                        la maille
!
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
double precision :: x
double precision :: y
double precision :: z
double precision :: xx
double precision :: yy
double precision :: zz
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: i1
integer :: i1m1
integer :: i2
integer :: inc1
integer :: inc2
integer :: inc3
integer :: j1
integer :: j1m1
integer :: j2
integer :: k1
integer :: k1m1
integer :: k2
integer :: m
integer :: n
integer :: n0
integer :: nid
integer :: nijd
integer :: njd
!
!-----------------------------------------------------------------------
!
      dimension x(ip21),y(ip21),z(ip21)
      dimension xx(ip00),yy(ip00),zz(ip00)
!
      ind(i,j,k)=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
!
      n0=npn(l)
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
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd = nid*njd
!
      inc1= inc(1,0,0)
      inc2= inc(0,1,0)
      inc3= inc(0,0,1)
!
      do k=k1m1,k2
      do j=j1m1,j2
      do i=i1m1,i2
      n=ind(i,j,k)
      m=n-n0
!
      xx(m)=.125*(x(n          )+x(n+          inc3) &
                 +x(n+inc1     )+x(n+inc1     +inc3) &
                 +x(n     +inc2)+x(n     +inc2+inc3) &
                 +x(n+inc1+inc2)+x(n+inc1+inc2+inc3))
      yy(m)=.125*(y(n          )+y(n+          inc3) &
                 +y(n+inc1     )+y(n+inc1     +inc3) &
                 +y(n     +inc2)+y(n     +inc2+inc3) &
                 +y(n+inc1+inc2)+y(n+inc1+inc2+inc3))
      zz(m)=.125*(z(n          )+z(n+          inc3) &
                 +z(n+inc1     )+z(n+inc1     +inc3) &
                 +z(n     +inc2)+z(n     +inc2+inc3) &
                 +z(n+inc1+inc2)+z(n+inc1+inc2+inc3))
      enddo
      enddo
      enddo
!
      return
      end subroutine
end module
