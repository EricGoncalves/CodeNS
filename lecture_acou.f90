module mod_lecture_acou
  implicit none
contains
  subroutine lecture_acou(l,v)
!
!***********************************************************************
!
!_DA  DATE_C : mai 2002 -- AUTEUR : Eric GONCALVES / SINUMEF
!
!     ACT
!_A    Lecture des variables (v) pour calculs acoustiques
!
!
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
!_O    v          : arg real(ip11,ip60 ) ; variables de calcul
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::    i,  i1,  i2,i2m1,   j
    integer          ::   j1,  j2,j2m1,   k,  k1
    integer          ::   k2,k2m1,   l,   n,  n0
    integer          ::  nid,nijd, njd
    double precision :: v(ip11,ip60)
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
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
!
    open(unit=77,file='sortieacou',form='formatted',status='unknown')
!
    do j=j1,j2m1
       do i=i1,i2m1
          n=ind(i,j,0)
          read(77,*) v(n,1),v(n,2),v(n,3),v(n,4),v(n,5)
       enddo
    enddo
!
!      read(98)(( v(ind(i,0,k),1),i=i1,i2m1),k=k1,k2m1)
!      read(98)(( v(ind(i,0,k),2),i=i1,i2m1),k=k1,k2m1)
!      read(98)(( v(ind(i,0,k),3),i=i1,i2m1),k=k1,k2m1)
!      read(98)(( v(ind(i,0,k),5),i=i1,i2m1),k=k1,k2m1)
!
!$OMP END MASTER
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine lecture_acou
end module mod_lecture_acou
