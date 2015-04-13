module mod_cvcvg
  implicit none
contains
  subroutine cvcvg( &
       l, &
       x,y,z, &
       xx,yy,zz)
!
!***********************************************************************
!
!     ACT
!_A    Transfert de coordonnees aux noeuds a coordonnees aux noeuds.
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
!_O    xx         : arg real(ip00      ) ; coordonnee sur l'axe x
!_O    yy         : arg real(ip00      ) ; coordonnee sur l'axe y
!_O    zz         : arg real(ip00      ) ; coordonnee sur l'axe z
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::    i,  i1,  i2,   j,  j1
    integer          ::   j2,   k,  k1,  k2,   l
    integer          ::    m,   n,  n0, nid,nijd
    integer          ::  njd
    double precision ::  x(ip21),xx(ip00), y(ip21),yy(ip00), z(ip21)
    double precision :: zz(ip00)
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!

!
    n0=npn(l)
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    do k=k1,k2
       do j=j1,j2
!$OMP SIMD
          do i=i1,i2
             n=ind(i,j,k)
             m=n-n0
!
             xx(m)=x(n)
             yy(m)=y(n)
             zz(m)=z(n)
!
          enddo
       enddo
    enddo
!
!$OMP END MASTER
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine cvcvg
end module mod_cvcvg
