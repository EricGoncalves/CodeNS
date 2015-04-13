module mod_dpbdc
  implicit none
contains
  subroutine dpbdc( &
       mfbe,img, &
       mnc, &
       pi,pj,pk)
!
!***********************************************************************
!
!     ACT
!_A    Ecriture formatee sur le fichier fimp des donnees specifiques
!_A     des points des frontieres coincidentes avec des domaines structures.
!
!     INP
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule coincidente
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    mpc        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux front coinc
!_I    ndcc       : com int (mtb       ) ; numero du dom coicident
!
!     LOC
!_L    pi         : arg real(ip40      ) ; indice i d'une maille recouvrante
!_L    pj         : arg real(ip40      ) ; indice j d'une maille recouvrante
!_L    pk         : arg real(ip40      ) ; indice k d'une maille recouvrante
!
!     COM
!_C   Impression des indices coincidents pour fr. coincidentes.
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
    integer          ::        ii,      img,       jj,       kk,       lc
    integer          ::       lcm,       ll,        m,      m0c,     mfbe
    integer          ::      mfbi,    mfbim,mnc(ip43),       mt,      n0c
    integer          ::        nc,      nid,      njd,       nn
    double precision :: pi(ip40),pj(ip40),pk(ip40)
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!
    mfbi=nfei(mfbe)
    mfbim=mfbi+(img-1)*mtb
!
    lc=ndcc(mfbi)
    lcm=lc+(img-1)*lz
!
    n0c=npc(lcm)
    nid = id2(lcm)-id1(lcm)+1
    njd = jd2(lcm)-jd1(lcm)+1
!
    ll=ndlb(mfbi)
!
    m0c=mpc(mfbim)
    mt =mmb(mfbim)
!
    write(imp,1900) mfbe,ll,lc,img
!
    do m=1,mt
       nc=mnc(m0c+m)
       nn=nc-1-n0c
       ii=id1(lcm)+mod(nn,nid)
       nn=(nn-(ii-id1(lcm)))/nid
       jj=jd1(lcm)+mod(nn,njd)
       kk=kd1(lcm)+(nn-(jj-jd1(lcm)))/njd
       pi(m)=ii
       pj(m)=jj
       pk(m)=kk
       if (mod(m,50).eq.1) write(imp,1910)
       write(imp,1920) m,nc,nint(pi(m)),nint(pj(m)),nint(pk(m))
    enddo
!
1900 format('1 Frontiere ',i5,' appartenant au domaine ',i5, &
         ' coincidente avec le domaine ',i5, &
         ' grille numero ',i3//)
1910 format(1x,'    m', &
         4x,'  nc',2x,' ic',2x,' jc', 2x,' kc'/)
1920 format(1x,i5,2x,i6,3(2x,i3))
!
!$OMP END MASTER
    return
  end subroutine dpbdc
end module mod_dpbdc
