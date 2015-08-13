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
    integer          ::         mt,        n0,       nid,      nijd
    integer          ::        njd
    integer,allocatable :: ncbd(:),ncbd1(:)
!
!-----------------------------------------------------------------------
!
    character(len=2 ) :: indfl
!
   allocate(ncbd1(ip41))
   ncbd1(1:size(ncbd))=ncbd
   deallocate(ncbd)
   
!
    n0=npc(l)
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    m=m0
!
    select case(indfl)
    case('i1')
       ii=imin-1
       mt=(kmax-kmin)*(jmax-jmin)
       ip41=ip41+mt
       allocate(ncbd(ip41))
       ncbd(:size(ncbd1))=ncbd1
       do k=kmin,kmax-1
          do j=jmin,jmax-1
             m=m+1
             ncbd(m)=ind(ii,j,k)
          enddo
       enddo

!
    case('i2')
       ii=imax
       mt=(kmax-kmin)*(jmax-jmin)
       ip41=ip41+mt
       allocate(ncbd(ip41))
       ncbd(:size(ncbd1))=ncbd1
       do k=kmin,kmax-1
          do j=jmin,jmax-1
             m=m+1
             ncbd(m)=ind(ii,j,k)
          enddo
       enddo

!
    case('j1')
       jj=jmin-1
       mt=(kmax-kmin)*(imax-imin)
       ip41=ip41+mt
       allocate(ncbd(ip41))
       ncbd(:size(ncbd1))=ncbd1
       do k=kmin,kmax-1
          do i=imin,imax-1
             m=m+1
             ncbd(m)=ind(i,jj,k)
          enddo
       enddo
!
    case('j2')
       jj=jmax
       mt=(kmax-kmin)*(imax-imin)
       ip41=ip41+mt
       allocate(ncbd(ip41))
       ncbd(:size(ncbd1))=ncbd1
       do k=kmin,kmax-1
          do i=imin,imax-1
             m=m+1
             ncbd(m)=ind(i,jj,k)
          enddo
       enddo
!
    case('k1')
       kk=kmin-1
       mt=(jmax-jmin)*(imax-imin)
       ip41=ip41+mt
       allocate(ncbd(ip41))
       ncbd(:size(ncbd1))=ncbd1
       do j=jmin,jmax-1
          do i=imin,imax-1
             m=m+1
             ncbd(m)=ind(i,j,kk)
          enddo
       enddo
!
    case('k2')
       kk=kmax
       mt=(jmax-jmin)*(imax-imin)
       ip41=ip41+mt
       allocate(ncbd(ip41))
       ncbd(:size(ncbd1))=ncbd1
       do j=jmin,jmax-1
          do i=imin,imax-1
             m=m+1
             ncbd(m)=ind(i,j,kk)
          enddo
       enddo
!
    end select
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine initis
end module mod_initis
