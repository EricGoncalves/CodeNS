module mod_inbdb
  implicit none
contains
  subroutine inbdb( &
       ncbd,ncin, &
       mfbe,clmf,kibdb, &
       ibdcst,ibdcfl,ibddim,nvbc,vbc,bceqt)
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des donnees de base pour une frontiere, soit
!_A    le type de condition physique dont releve cette frontiere
!_A    et l'indice du centre de la maille qui s'appuie sur la frontiere.
!
!     INP
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    clmf       : arg char             ; type de cond lim a appliquer
!_I    kibdb      : arg int              ; cle initialisation de base front
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    kfb        : com int              ; unite logiq, tableaux de base front
!
!     OUT
!_O    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_O    cl         : com char(mtb       ) ; type de cond lim a appliquer
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use maillage
    use sortiefichier
    use mod_initbs
    use mod_inbdbfl
    use mod_inbdbdf
    use mod_inbdbst
    implicit none
  integer          ::     ibdcfl,    ibdcst,    ibddim,       img,     kibdb
  integer          ::          l,        lm,         m,        m0,      mfbe
  integer          ::       mfbi,     mfbim,      mflm,        ml,        mt
  integer          :: ncbd(ip41),ncin(ip41),        nv,      nvbc
  double precision :: bceqt(ip41,neqt),  vbc(ista*lsta)
!
!-----------------------------------------------------------------------
!
    character(len=4 ) :: clmf
!
    mfbi=nfei(mfbe)
    cl(mfbi)=clmf
    l=ndlb(mfbi)
!
!  Utilisation d'un etat thermodynamique pour appliquer une condition limite
!
    if(ibdcst.ne.0) then
       call inbdbst(clmf,ibdcst,nvbc,vbc)
    endif
!
!  Utilisation des informations contenues dans la grille de donnees
!
    if(nvbc.ne.0) then
       call inbdbdf(clmf,ibddim,nvbc,vbc)
    endif
!
    nbdc(mfbi)=nvbc
    do nv=1,nvbc
       bc(mfbi,nv)=vbc(nv)
    enddo
!
!  Preparation du tableau de conditions aux limites
!
    do nv=1,nvbc
       do img=1,lgx
          mflm = mfbi + (img-1)*mtb
          mt = mmb(mflm)
          do m=1,mt
             ml=mpb(mflm)+m
             bceqt(ml,nv)=vbc(nv)
          enddo
       enddo
    enddo
!
    if(ibdcfl.ne.0) then
       call inbdbfl(ibdcfl,mfbi,bceqt,mpb)
    endif
!
    do img=1,lgx
!
       lm=l+(img-1)*lz
       mfbim=mfbi+(img-1)*mtb
       mt=mmb(mfbim)
       m0=mpb(mfbim)
!
       if(kibdb.eq.1) then
          call initbs( &
               mfbim,lm,indfl(mfbi), &
               ncin,ncbd, &
               mt,m0)
       elseif(kibdb.eq.0) then
!            call readfb( &
!                 kfb,ncin, &
!                 mt,m0)
       endif
    enddo
!
    return
  end subroutine inbdb
end module mod_inbdb
