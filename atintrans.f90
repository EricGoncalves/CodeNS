module mod_atintrans
  implicit none
contains
  subroutine atintrans(ncin,fgam)
!
!***********************************************************************
!
!     ACT
!_A   Lecture des donnees pour fixer la transition
!_A   initialisation du tableau de la fonction d'intermittence "fgam"
!
!_A   - Recherche du mot cle TRANSITION en colonne 1
!_A   - lecture du nombre de zones laminaires
!_A   - lecture des indices extremes dans les 3 directions pour chaque zone
!_A     les indices sont les numeros des noeuds du maillage
!
!     INP
!_I    mdimtnx    : com int              ; nbr de pts de ttes les front
!_I                                        a normales stockees
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
!_I    kd2        : com int (lt        ) ; indice max en k fictif
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    iminb      : com int (mtt       ) ; indice min en i d'une front
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une front
!_I    jminb      : com int (mtt       ) ; indice min en j d'une front
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une front
!_I    kminb      : com int (mtt       ) ; indice min en k d'une front
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une front
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!
!     LOC
!_L    l          : arg int              ; numero de domaine
!_L    ilmin       : arg int              ; indice min en i
!_L    ilmax       : arg int              ; indice max en i
!_L    jlmin       : arg int              ; indice min en j
!_L    jlmax       : arg int              ; indice max en j
!_L    klmin       : arg int              ; indice min en k
!_L    klmax       : arg int              ; indice max en k
!
!     OUT
!_O    fgam       : arg real(ip42      ) ; fonction d'intermittence pour
!_O                                        transition
!_O    ktransi    : com int              ; cle transition 1->transition fixee
!_O                                                       0->pas de transition
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use sortiefichier
    use modeleturb
    implicit none
    integer          ::         i1,        i2,       idm,     ilmax,     ilmin
    integer          ::      imaxf,     iminf,        j1,        j2,     jlmax
    integer          ::      jlmin,     jmaxf,     jminf,        k1,        k2
    integer          ::      klmax,     klmin,     kmaxf,     kminf,         l
    integer          ::        lig,       m0b,       m0n,        m1,     m1deb
    integer          ::      m1fin,     m1max,   m1maxm1,     m1min,        m2
    integer          ::      m2deb,     m2fin,     m2max,   m2maxm1,     m2min
    integer          ::      mfacn,       mfe,       mfl,         n,       nci
    integer          :: ncin(ip41),       ncj,       nck,       nfr,     nfrmx
    integer          ::        nid,      nijd,       njd
    double precision :: fgam(ip42)
!
!-----------------------------------------------------------------------

    character(len=80) ::  ligne
!
    open(99,file='fatdon',status='old',err=100)
!
!     recherche de la ligne "TRANSITION"
!
    lig=0
    do
       lig=lig+1
       read(99,300,err=304,end=301)ligne
       if(ligne(1:6).EQ.'TRANSI') exit
    enddo
!
!     transition fixee
    ktransi=1
!
!     initialisation de "fgam" a 1
    do n=1,mdimtnx
       fgam(n)=1.
    enddo
!
!     lecture des frontieres laminaires
!
    lig=lig+1
    read(99,*,err=302)nfrmx
    if(nfrmx.le.0) then
!       le nombre de zones laminaires doit etre positif
       write(imp,'(/,"!!!atintrans: domaine laminaire mal defini",14x,"nfrmx=",i5)')nfrmx
       stop
    endif
!
    do nfr=1,nfrmx
       read(99,*,err=303)mfe,ilmin,ilmax,jlmin,jlmax,klmin,klmax
!
!       calcul numero facette de la frontiere
!
       mfl=nfei(mfe)
       l =ndlb(mfl)
!
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
       nci=1
       ncj=nid
       nck=nijd
!
       m0b=mpb(mfl)
       m0n=mpn(mfl)
       iminf=iminb(mfl)
       imaxf=imaxb(mfl)
       jminf=jminb(mfl)
       jmaxf=jmaxb(mfl)
       kminf=kminb(mfl)
       kmaxf=kmaxb(mfl)
       m1min=1
       m2min=1
!
!           m1deb m1fin : limites de variation des facettes parois laminaires
!           m2deb m2fin : dans numerotation des frontieres a normale stockee
!

!                     ilmin           ilmax
!                       |               |
!                       |   laminaire   |
!                       |===============|
!       ..........x-----------------------------x.............
!               iminf   frontiere "mfl"       imaxf  indice en "i" domaine l
!                |
!            m0n=mpn(mfl)                            pointeur normales stockees
!
!
       if (iminf.eq.imaxf) then
!         frontiere i=cste
          m1max=jmaxf-jminf+1
          m1maxm1=m1max-1
          m2max=kmaxf-kminf+1
          m2maxm1=m2max-1
          m1deb=jlmin-jminf+1
          m1fin=jlmax-jminf
          m2deb=klmin-kminf+1
          m2fin=klmax-kminf
!
       elseif (jminf.eq.jmaxf) then
!         frontiere j=cste
          m1max=imaxf-iminf+1
          m1maxm1=m1max-1
          m2max=kmaxf-kminf+1
          m2maxm1=m2max-1
          m1deb=ilmin-iminf+1
          m1fin=ilmax-iminf
          m2deb=klmin-kminf+1
          m2fin=klmax-kminf
!
       elseif (kminf.eq.kmaxf) then
!         frontiere k=cste
          m1max=imaxf-iminf+1
          m1maxm1=m1max-1
          m2max=jmaxf-jminf+1
          m2maxm1=m2max-1
          m1deb=ilmin-iminf+1
          m1fin=ilmax-iminf
          m2deb=jlmin-jminf+1
          m2fin=jlmax-jminf
       endif
!
!
       idm=m1max-m1min
       do m2=m2deb,m2fin
!         boucle sur les bandes
          do m1=m1deb,m1fin
!           boucle sur les cellules de la bande
             mfacn=m0n+m1+(m2-1)*idm
             fgam(mfacn)=0.
          enddo
       enddo
!     fin boucle sur frontieres laminaires
    enddo
    close(99)
!
300 format(a80)
!
    return
!
301 continue
!     fin de fichier. Pas de transition
    ktransi=0
    write(imp,'(/,"===>atintrans: ktransi=0. Calcul tout laminaire ou tout turbulent")')
    close(99)
    return
!
100 continue
    write(imp,'(/,"!!!atintrans: erreur ouverture fichier fatdon")')
    stop
!
302 continue
    write(imp,'(/,"!!!atintrans: erreur nombre de frontieres laminaires ")')
    stop
!
303 continue
    write(imp,'(/,"!!!atintrans: erreur lecture frontiere ligne",i4)')lig
    stop
!
304 continue
    write(imp,'(/,"!!!atintrans: erreur positionnement dans fatdon")')
    stop
!
    return
  end subroutine atintrans
end module mod_atintrans
