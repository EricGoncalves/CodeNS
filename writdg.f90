module mod_writdg
  implicit none
contains
  subroutine writdg( &
       l,kdg, &
       imin,imax,jmin,jmax,kmin,kmax, &
       x,y,z)
!
!***********************************************************************
!
!     ACT
!_A    Ecriture sur l'unite logique kdg des coordonnees (x, y, z)
!_A    de points d'un domaine structure.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    kdg        : arg int              ; unite logique, maillage
!_I    imin       : arg int              ; indice min en i
!_I    imax       : arg int              ; indice max en i
!_I    jmin       : arg int              ; indice min en j
!_I    jmax       : arg int              ; indice max en j
!_I    kmin       : arg int              ; indice min en k
!_I    kmax       : arg int              ; indice max en k
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use mod_mpi
    implicit none
    integer          ::    i,imax,imin,   j,jmax
    integer          :: jmin,   k, kdg,kmax,kmin
    integer          ::    l, nid,nijd, njd,ll
    integer(8)       :: pos
    double precision :: x(ip00),y(ip00),z(ip00)
!
!-----------------------------------------------------------------------
!
!

!
    ll=bl_to_bg(l)
    if(ll.eq.1) rewind kdg
    pos=FTELL(kdg)
    call start_keep_order(ll,bg_to_proc,pos)
    CALL my_FSEEK(kdg, pos)
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    write(kdg)(((x(ind(i,j,k)),i=imin,imax),j=jmin,jmax),k=kmin,kmax)
    write(kdg)(((y(ind(i,j,k)),i=imin,imax),j=jmin,jmax),k=kmin,kmax)
    write(kdg)(((z(ind(i,j,k)),i=imin,imax),j=jmin,jmax),k=kmin,kmax)
!
    pos=FTELL(kdg)
    call END_KEEP_ORDER(ll,bg_to_proc,pos)
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine writdg
end module mod_writdg
