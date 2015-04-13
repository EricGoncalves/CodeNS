module mod_extmhg
  implicit none
contains
  subroutine extmhg(l,x,y,z,ex1,ex2)
!
!***********************************************************************
!
!     ACT
!_A    Extrapolation des coordonnees du maillage aux points fictifs.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    ex1        : arg real             ; premier coef d'interpolation
!_I                                        du couple de coef ex1, ex2
!_I    ex2        : arg real             ; deuxieme coef d'interpolation
!_I                                        du couple de coef ex1, ex2
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
!     I/O
!_/    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_/    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_/    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::    i,  i1,i1m1,  i2,i2p1
    integer          ::   id,iinc,  is,   j,  j1
    integer          :: j1m1,  j2,j2p1,  jd,jinc
    integer          ::   js,   k,  k1,k1m1,  k2
    integer          :: k2p1,  kd,kinc,  ks,   l
    integer          ::    n,  n0, nci,nci2, ncj
    integer          :: ncj2, nck,nck2, nid,nijd
    integer          ::  njd
    double precision ::     ex1,    ex2,x(ip21),y(ip21),z(ip21)
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!
!----- extrapolation d' ordre 0 ou 1
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
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
    i2p1=i2+1
    j2p1=j2+1
    k2p1=k2+1
!
    nci = inc(1,0,0)
    ncj = inc(0,1,0)
    nck = inc(0,0,1)
    nci2= inc(2,0,0)
    ncj2= inc(0,2,0)
    nck2= inc(0,0,2)
!
!     extrapolation des coordonnees
!
    iinc = i2p1-i1m1
    jinc = j2p1-j1m1
    kinc = k2p1-k1m1
!
    ks = -1
    do k = k1m1,k2p1,kinc
       ks = -ks
       do j = j1,j2
!$OMP SIMD
          do i = i1,i2
             n = ind(i,j,k)
             x(n) = ex1*x(n+ks*nck)+ex2*x(n+ks*nck2)
             y(n) = ex1*y(n+ks*nck)+ex2*y(n+ks*nck2)
             z(n) = ex1*z(n+ks*nck)+ex2*z(n+ks*nck2)
          enddo
       enddo
    enddo
!
    do k = k1m1,k2p1
       do j = j1,j2
          is = -1
          do i = i1m1,i2p1,iinc
             is = -is
             n = ind(i,j,k)
             x(n) = ex1*x(n+is*nci)+ex2*x(n+is*nci2)
             y(n) = ex1*y(n+is*nci)+ex2*y(n+is*nci2)
             z(n) = ex1*z(n+is*nci)+ex2*z(n+is*nci2)
          enddo
       enddo
    enddo
!
    do k = k1m1,k2p1
       js = -1
       do j = j1m1,j2p1,jinc
          js = -js
!$OMP SIMD
          do i = i1m1,i2p1
             n = ind(i,j,k)
             x(n) = ex1*x(n+js*ncj)+ex2*x(n+js*ncj2)
             y(n) = ex1*y(n+js*ncj)+ex2*y(n+js*ncj2)
             z(n) = ex1*z(n+js*ncj)+ex2*z(n+js*ncj2)
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
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine extmhg
end module mod_extmhg
