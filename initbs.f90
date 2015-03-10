module mod_initbs
implicit none
contains
      subroutine initbs( &
                 mfbi,l,indf, &
                 ncin,ncbd, &
                 mt,m0)
!
!***********************************************************************
!
!     ACT
!_A    Determination des indices dans les tableaux tous domaines
!_A    des cellules interieures adjacentes a une frontiere de
!_A    domaine structure.
!
!     INP
!_I    mfbi       : arg int              ; numero interne de frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mt         : arg int              ; nombre de points de la frontiere
!_I    m0         : arg int              ; pointeur fin de front precedente
!_I                                        dans tab toutes front
!_I    indfl      : com char(mtb       ) ; type de plan de la frontiere
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_I    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_I    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!
!     OUT
!_O    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_O                                        cell. interieure adjacente a la front
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use boundary
implicit none
integer :: mfbi
integer :: l
integer :: ncin
integer :: ncbd
integer :: mt
integer :: m0
integer :: m
integer :: mdnmfb
integer :: ml
integer :: nid
integer :: njd
!
!-----------------------------------------------------------------------
!
      character(len=2 ) :: indf
      dimension ncin(ip41),ncbd(ip41)
!
      nid=id2(l)-id1(l)+1
      njd=jd2(l)-jd1(l)+1
!
      if (indf(1:1).eq.'i') &
              mdnmfb=isign(1,-iabs(ii1(l)-iminb(mfbi)))*1
      if (indf(1:1).eq.'j') &
              mdnmfb=isign(1,-iabs(jj1(l)-jminb(mfbi)))*nid
      if (indf(1:1).eq.'k') &
              mdnmfb=isign(1,-iabs(kk1(l)-kminb(mfbi)))*nid*njd
!
      do m=1,mt
      ml=m0+m
      ncin(ml)=ncbd(ml)+mdnmfb
      enddo
!
      return
      end
end module
