module mod_initis
  implicit none
contains
  subroutine initis( &
       l,imin,imax,jmin,jmax,kmin,kmax, &
       indfl,ncbd, &
       mt,m0)
!
!***********************************************************************
!
!     ACT
!_A    Determination du nombre de facettes d'une frontiere d'un
!_A    domaine structure et des indices dans les tableaux tous
!_A    domaines des cellules fictives formant cette frontiere.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    imin       : arg int              ; indice min en i
!_I    imax       : arg int              ; indice max en i
!_I    jmin       : arg int              ; indice min en j
!_I    jmax       : arg int              ; indice max en j
!_I    kmin       : arg int              ; indice min en k
!_I    kmax       : arg int              ; indice max en k
!_I    indfl      : arg char             ; type de plan de la frontiere
!_I    m0         : arg int              ; pointeur fin de front precedente
!_I                                        dans tab toutes front
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!
!     OUT
!_O    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_O                                        cellule frontiere fictive
!_O    mt         : arg int              ; nombre de points de la frontiere
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::          i,        ii,      imax,      imin,         j
    integer          ::         jj,      jmax,      jmin,         k,        kk
    integer          ::       kmax,      kmin,         l,         m,        m0
    integer          ::         mt,        n0,ncbd(ip41),       nid,      nijd
    integer          ::        njd
!
!-----------------------------------------------------------------------
!
    character(len=2 ) :: indfl
!$OMP MASTER
!

!
    n0=npc(l)
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    m=m0
!
    if (indfl.eq.'i1') then
       ii=imin-1
       do k=kmin,kmax-1
!$OMP SIMD
          do j=jmin,jmax-1
             m=m+1
             ncbd(m)=ind(ii,j,k)
          enddo
       enddo
       mt=(kmax-kmin)*(jmax-jmin)
!
    else if (indfl.eq.'i2') then
       ii=imax
       do k=kmin,kmax-1
!$OMP SIMD
          do j=jmin,jmax-1
             m=m+1
             ncbd(m)=ind(ii,j,k)
          enddo
       enddo
       mt=(kmax-kmin)*(jmax-jmin)
!
    elseif (indfl.eq.'j1') then
       jj=jmin-1
       do k=kmin,kmax-1
!$OMP SIMD
          do i=imin,imax-1
             m=m+1
             ncbd(m)=ind(i,jj,k)
          enddo
       enddo
       mt=(kmax-kmin)*(imax-imin)
!
    elseif (indfl.eq.'j2') then
       jj=jmax
       do k=kmin,kmax-1
!$OMP SIMD
          do i=imin,imax-1
             m=m+1
             ncbd(m)=ind(i,jj,k)
          enddo
       enddo
       mt=(kmax-kmin)*(imax-imin)
!
    elseif (indfl.eq.'k1') then
       kk=kmin-1
       do j=jmin,jmax-1
!$OMP SIMD
          do i=imin,imax-1
             m=m+1
             ncbd(m)=ind(i,j,kk)
          enddo
       enddo
       mt=(jmax-jmin)*(imax-imin)
!
    elseif (indfl.eq.'k2') then
       kk=kmax
       do j=jmin,jmax-1
!$OMP SIMD
          do i=imin,imax-1
             m=m+1
             ncbd(m)=ind(i,j,kk)
          enddo
       enddo
       mt=(jmax-jmin)*(imax-imin)
!
    end if
!
!$OMP END MASTER
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine initis
end module mod_initis
