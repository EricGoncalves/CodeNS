module mod_dpbdn
  implicit none
contains
  subroutine dpbdn( &
       mfbe,img, &
       ncbd,nxn,nyn,nzn)
!
!***********************************************************************
!
!     ACT
!_A    Ecriture formatee sur le fichier fimpdes donnees specifiques
!_A    des points des frontieres a normales stockees.
!
!     INP
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!
!     COM
!_C   Impression des normales interieures de composantes nxn ,nyn ,nzn.
!_C   Le tableau ncbd doit avoir ete rempli au prealable.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use maillage
    use sortiefichier
    implicit none
    integer          ::          i,       img,         j,         k,         l
    integer          ::         lm,       m0b,       m0n,        mb,      mfbe
    integer          ::       mfbi,     mfbim,        mm,        mn,        mt
    integer          ::          n,       n0c,ncbd(ip41),       nid,      nijd
    integer          ::        njd
    double precision :: nxn(ip42),nyn(ip42),nzn(ip42)
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!
!
    mfbi=nfei(mfbe)
!
    mfbim=mfbi+(img-1)*mtb
!
    mt=mmb(mfbim)
    m0b=mpb(mfbim)
    m0n=mpn(mfbim)
!
    l=ndlb(mfbi)
    lm=l+(img-1)*lz
!
    n0c =npc(lm)
    nid =id2(lm)-id1(lm)+1
    njd =jd2(lm)-jd1(lm)+1
    nijd=nid*njd
!
    do mm=1,mt
       mb=m0b+mm
       mn=m0n+mm
       n=ncbd(mb)
       i=iin(n)
       j=jjn(n)
       k=kkn(n)
       if (mod(mm,50).eq.1) write(imp,1900)
       write(imp,1910) mm,l,img,n,i,j,k,nxn(mn),nyn(mn),nzn(mn)
    enddo
!
1900 format('1     m',4x,'  l',2x,'img',2x,'     n',2x,'  i', &
         2x,'  j',2x,'  k',2x,'    nxn   ',2x,'    nyn   ', &
         2x,'    nzn   '/)
1910 format(1x,i5,4x,i3,2x,i3,2x,i6,3(2x,i3),3(2x,f10.6))
!
!$OMP END MASTER
    return
  contains
    function    iin(n)
      implicit none
      integer          :: iin,  n
      iin=id1(lm)+mod(n-n0c-1,nid)
    end function iin
    function    jjn(n)
      implicit none
      integer          :: jjn,  n
      jjn=jd1(lm)+mod((n-n0c-1-(iin(n)-id1(lm)))/nid,njd)
    end function jjn
    function    kkn(n)
      implicit none
      integer          :: kkn,  n
      kkn=kd1(lm)+ &
           (n-n0c-1-(iin(n)-id1(lm))-(jjn(n)-jd1(lm))*nid)/nijd
    end function kkn
  end subroutine dpbdn
end module mod_dpbdn
