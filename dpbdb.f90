module mod_dpbdb
  implicit none
contains
  subroutine dpbdb( &
       mfbe,img, &
       ncbd,ncin)
!
!***********************************************************************
!
!     ACT
!_A    Ecriture formatee sur le fichier fimp des donnees de base
!_A     des points des frontieres des domaines structures.
!
!     INP
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
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
!
!     COM
!_C   Impression des indices immediatement interieurs : ncin.
!_C   Le tableau ncbd doit avoir ete rempli au prealable.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use sortiefichier
    implicit none
    integer          ::          i,       img,         j,         k,         l
    integer          ::         lm,         m,        m0,      mfbe,      mfbi
    integer          ::      mfbim,        mm,        mt,         n,       n0c
    integer          :: ncbd(ip41),ncin(ip41),       nid,      nijd,       njd
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!
!
    mfbi=nfei(mfbe)
!
    l=ndlb(mfbi)
    lm=l+(img-1)*lz
!
    n0c =npc(lm)
    nid =id2(lm)-id1(lm)+1
    njd =jd2(lm)-jd1(lm)+1
    nijd=nid*njd
!
    mfbim=mfbi+(img-1)*mtb
!
    m0=mpb(mfbim)
    mt=mmb(mfbim)
!
    do mm=1,mt
       m=m0+mm
       n=ncbd(m)
       i=iin(n)
       j=jjn(n)
       k=kkn(n)
       if (mod(mm,50).eq.1) write(imp,1900)
       write(imp,1910) mm,l,img,i,j,k,ncin(m)
    enddo
!
1900 format('1     m',4x,'  l',2x,'img',2x,'  i',2x,'  j', 2x,'  k',2x,'      ncin'/)
1910 format(1x,i5,2x,5(2x,i3),5x,i7)
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
  end subroutine dpbdb
end module mod_dpbdb
