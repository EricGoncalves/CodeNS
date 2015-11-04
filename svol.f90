module mod_svol
  implicit none
contains
  subroutine svol( &
       l,x,y,z, &
       sn,lgsnlt, &
       vol, &
       siv,sjv,skv)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des volumes de toutes les cellules d'un domaine structure.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    sn         : arg real(lgsnlt,
!_I                          nind,ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    lgsnlt     : arg int              ; nombre de noeuds du dom (dont fic.)
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
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
!_O    vol        : arg real(ip11      ) ; volume d'une cellule
!
!     LOC
!_L    siv        : arg real(ip00      ) ; flux vect(x,y,z) a travers facette i
!_L    sjv        : arg real(ip00      ) ; flux vect(x,y,z) a travers facette j
!_L    skv        : arg real(ip00      ) ; flux vect(x,y,z) a travers facette k
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::      i,    i1,  i1m1,    i2,  i2p1
    integer          ::      j,    j1,  j1m1,    j2
    integer          ::   j2p1,     k,    k1,  k1m1
    integer          ::     k2,  k2p1,     l,lgsnlt
    integer          ::      m,    m1,    m2,    m4,    m5
    integer          ::      n,   n0c,   n0n,   nid,  nijd
    integer          ::    njd
    double precision ::            siv(ip00),           sjv(ip00),           skv(ip00),sn(lgsnlt,nind,ndir),           vol(ip11)
    double precision ::              x(ip21),             y(ip21),             z(ip21)
!
!-----------------------------------------------------------------------
!
!



!
    n0n=npn(l)
    n0c=npc(l)
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
!     calcul de volume des cellules
!
    do k=k1m1,k2p1
       do j=j1m1,j2
          do i=i1m1,i2
             n=indn(i,j,k)
             m=n-n0n
             skv(m)=(x(n)+x(n+inc(1,0,0))+x(n+inc(0,1,0))+x(n+inc(1,1,0))) &
                  *sn(m,3,1)+ &
                  (y(n)+y(n+inc(1,0,0))+y(n+inc(0,1,0))+y(n+inc(1,1,0))) &
                  *sn(m,3,2)+ &
                  (z(n)+z(n+inc(1,0,0))+z(n+inc(0,1,0))+z(n+inc(1,1,0))) &
                  *sn(m,3,3)
          enddo
       enddo
    enddo
!
    do j=j1m1,j2p1
       do k=k1m1,k2
          do i=i1m1,i2
             n=indn(i,j,k)
             m=n-n0n
             sjv(m)=(x(n)+x(n+inc(1,0,0))+x(n+inc(0,0,1))+x(n+inc(1,0,1))) &
                  *sn(m,2,1)+ &
                  (y(n)+y(n+inc(1,0,0))+y(n+inc(0,0,1))+y(n+inc(1,0,1))) &
                  *sn(m,2,2)+ &
                  (z(n)+z(n+inc(1,0,0))+z(n+inc(0,0,1))+z(n+inc(1,0,1))) &
                  *sn(m,2,3)
          enddo
       enddo
    enddo
!
    do i=i1m1,i2p1
       do k=k1m1,k2
          do j=j1m1,j2
             n=indn(i,j,k)
             m=n-n0n
             siv(m)=(x(n)+x(n+inc(0,1,0))+x(n+inc(0,0,1))+x(n+inc(0,1,1))) &
                  *sn(m,1,1)+ &
                  (y(n)+y(n+inc(0,1,0))+y(n+inc(0,0,1))+y(n+inc(0,1,1))) &
                  *sn(m,1,2)+ &
                  (z(n)+z(n+inc(0,1,0))+z(n+inc(0,0,1))+z(n+inc(0,1,1))) &
                  *sn(m,1,3)
          enddo
       enddo
    enddo
!
    do k=k1m1,k2
       do j=j1m1,j2
          do i=i1m1,i2
             n = indc(i,j,k)
             m1 = n-n0c
             m2 = m1+inc(1,0,0)
             m4 = m1+inc(0,1,0)
             m5 = m1+inc(0,0,1)
             vol(n) = 0.25*(siv(m2)-siv(m1)+sjv(m4)-sjv(m1) &
                  +skv(m5)-skv(m1))/3
          enddo
       enddo
    enddo
!
    return
  contains
    function    indn(i,j,k)
      implicit none
      integer          ::    i,indn,   j,   k
      indn=n0n+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indn
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine svol
end module mod_svol
