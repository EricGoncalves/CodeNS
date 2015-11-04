module mod_atcaldis3
  use mod_at_lecdist
  use mod_at_ecrdist
  implicit none
contains
  subroutine atcaldis3( &
       x,y,z,nxn,nyn,nzn, &
       xpar,ypar,zpar,xcc,ycc,zcc,dist2, &
       dist,mnpar,fgam, &
       ncin,mnc,ncbd, &
       m1tb,m2tb,nfrtb)
!
!***********************************************************************
!
!     ACT
!_A    Calcul de la distance caracteristique de chaque cellule a une
!_A    paroi et rattachement de chaque cellule a une facette de paroi.
!_A    Optimisation par utilisation d'un niveau de maillage grossier.
!
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    reelmx     : com real             ; nombre reel grand
!_I    lzx        : com int              ; nbr total de domaines
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    mtnx       : com int              ; nbr total de frontieres  a normales stockes
!_I    mtcx       : com int              ; nbr total de frontieres coincidentes
!_I    mtrx       : com int              ; nbr total de frontieres recouvertes
!_I    mtax       : com int              ; nbr total de frontieres autres
!_L    ndimctbx   : com int              ; nbr de cellules de tts les domaines
!_L                                        (pts fictifs inclus)
!_L    mdimtnx    : com int              ; nbr de pts de ttes les front
!_L                                        a normales stockees
!_I    nnn        : com int (lt        ) ; nombre de noeuds du dom (dont fic.)
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    nnc        : com int (lt        ) ; nombre de cellules du dom (dont fic.)
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
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_L    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_L    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_L                                        en fct du numero externe
!_I    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une frontiere
!_I    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une frontiere
!_I    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une frontiere
!_I    nfbr       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere recouverte
!_I    nfbc       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere coincidente
!_I    nfbn       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une front a normales stockees
!_I    cl         : com char(mtb       ) ; type de cond lim a appliquer
!_I    indfl      : arg char             ; type de plan de la frontiere
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    kcaldis    : com int              ; cle calcul distance aux parois
!_I    reelmx     : com real             ; nombre reel grand
!
!     OUT
!_O    dist       : arg real(ip12      ) ; distance a la paroi
!_O    mnpar      : arg real(ip12      ) ; pointeur dans tableaux front normales
!_O                                        stockees du point de rattach normale
!_O    nbdko      : com int              ; nombre de parois a traiter
!_O    lbdko      : com int (mtt       ) ; numero interne des parois
!_O                                        (necessaire pour k-omega)
!_O    nbdrat     : com int (lz        ) ; nb parois rattachees au domaine l
!_O    lbdrat     : com int (mtb       ) ; no interne des parois rattachees
!_O                                        au domaine l
!_O    npbrat     : com int (lz        ) ; pointeur fin liste frontieres a
!_O                                        traiter pour domaine precedent
!
!     LOC
!_L    tn1        : arg real(ip00      ) ; tableau de travail -> xpar
!_L    tn2        : arg real(ip00      ) ; tableau de travail -> ypar
!_L    tn3        : arg real(ip00      ) ; tableau de travail -> zpar
!_L    tn4        : arg real(ip00      ) ; tableau de travail -> xcc
!_L    tn5        : arg real(ip00      ) ; tableau de travail -> ycc
!_L    tn6        : arg real(ip00      ) ; tableau de travail -> zcc
!_L    tn7        : arg real(ip00      ) ; tableau de travail ->dist2
!_L    xpar       : arg real(ip00      ) ; coordonnee x centre facette paroi
!_L    ypar       : arg real(ip00      ) ; coordonnee x centre facette paroi
!_L    zpar       : arg real(ip00      ) ; coordonnee y centre facette paroi
!_L    xcc        : arg real(ip00      ) ; coordonnee z centre cellule
!_L    ycc        : arg real(ip00      ) ; coordonnee y centre cellule
!_L    zcc        : arg real(ip00      ) ; coordonnee z centre cellule
!_L    dist2      : arg real(ip00      ) ; travail
!_L    nbd        : com int              ; nombre de frontieres a traiter
!_L    lbd        : com int (mtt       ) ; numero interne de front a traiter
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use maillage
    use constantes
    use boundary
    use mod_atccfp
    use mod_atindnor
    use mod_at_grdist
    use mod_at_lecopt
    use mod_atccc
    use mod_at_dlist
    use mod_atdist_3
    use mod_at_fidist
    use mod_mpi
    implicit none
    integer          ::         dm1,        dm2,        dm3,    idefaut,       ierr
    integer          ::         if0,        if1,    igr(lz),     isens3,    jgr(lz)
    integer          ::     kgr(lz),          l,     ldismx,          m,        m10
    integer          ::       m1max,      m1min, m1tb(ip00),        m20,      m2max
    integer          ::       m2min, m2tb(ip00),        m30,      m3max,      m3min
    integer          ::          mf,  mnc(ip43),mnpar(ip12),          n, ncbd(ip41)
    integer          ::  ncin(ip41),         nf,       nfbe,       nfbi,nfrtb(ip00)
    integer          ::          no,proc,nbdkog
    integer,allocatable :: nbdko_proc(:)
    double precision ::  dist(ip12),dist2(ip00), fgam(ip42),  nxn(ip42),  nyn(ip42)
    double precision ::   nzn(ip42),    x(ip21),  xcc(ip00), xpar(ip00),    y(ip21)
    double precision ::   ycc(ip00), ypar(ip00),    z(ip21),  zcc(ip00), zpar(ip00)
    double precision,allocatable :: raptat(:)
    allocate(raptat(mtb))
!
!-----------------------------------------------------------------------
!
!
!     verification compatibilite du dimensionnement
!     Les centres des facettes parois a normales stockees sont mis
!     dans xpar(ip00), ypar(ip00),zpar(ip00). Il faut donc :
!     ip00 > mdimtnx
    if(ip00 .le. mdimtnx) then
       write(imp,'(/,''!!!atcaldis3: ip00='',i6,4x,''inferieur a mdimtnx='',i6)')ip00,mdimtnx
       stop
    end if
!
!     -----------------------------------------------------------------
!     Liste des parois - numero interne
!                        nba : rang de traitement
!                        nbd : nombre de frontieres a traiter
!
    ierr =0
    nbd  =0
    nbdko=0
    do no=1,mtbx
       nfbe=nba(no)
       nfbi=nfei(nfbe)
!
       if((cl(nfbi)(1:2).eq.'pa').or.(cl(nfbi)(1:2).eq.'lp')) then
!
!
!         la frontiere est une paroi
!
          nbd         =nbd+1
          lbd(nbd)    =nfbi
          nbdko       =nbdko+1
          lbdko(nbdko)=nfbi
!
!         verification qu'elle est definie a normale stockee
!
          do m=1,mtnx
             if(nfbn(m).eq.nfbi) then
                ierr=-1
                exit
             endif
          enddo
!         la frontiere n'est pas a normale stockee
          ierr=ierr+1
       endif
    enddo
    if(ierr.ne.0) then
       write(imp,'(/,"!!!atcaldis3: pb declaration parois - arret")')
       stop
    end if
    if(nbd.eq.0) then
       write(imp,'(/,"!!!atcaldis3: pas de paroi pour calcul des distances aux parois - arret")')
       stop
    end if

    allocate(nbdko_proc(nprocs))
    call gather([nbdko],nbdko_proc,1)
    nbdkog=sum(nbdko_proc)
    allocate(lbdko_to_lbdkog(nbdko))
    allocate(lbdkog_to_proc(nbdkog))
    proc=0
    do mf=1,nbdkog
      if(mf>sum(nbdko_proc(1:proc+1))) then
        do l=proc,nprocs
           proc=proc+1
           if(mf<=sum(nbdko_proc(1:proc+1))) exit
        enddo
      endif
      lbdkog_to_proc(mf)=proc
    enddo
    do mf=1,nbdko
      lbdko_to_lbdkog(mf)=mf+sum(nbdko_proc(1:rank))
    enddo
    write(stderr,*) rank, nbdko_proc
    write(stderr,*) rank, lbdko_to_lbdkog
    write(stderr,*) rank, lbdkog_to_proc

!
!     initialisation pour chaque domaine de la liste des frontieres
!     qui seront explorees
!
!     Verifications
!
    if(lzx.gt.lz) then
       write(imp,'(/,"!!!atcaldis3: trop de domaines  lzx=",i3,3x,"lz=",i3)')lzx,lz
       stop
    end if
    if(nbd.ge.mtb) then
       write(imp,'(/,"!!!atcaldis3: trop de frontieres  nbd=",i3,3x,"mtb=",i3)')nbd,mtb
       stop
    end if
!
!     lecture des informations pour calcul des distances
    call at_lecopt( &
         igr,jgr,kgr,raptat,idefaut)
!
!     initialisation de la liste des frontieres a explorer pour chaque domaine
!
    if(idefaut.eq.1) then
!       toutes les parois pour tous les domaine
!
       npbrat(1)=0
       do l=1,lzx
          nbdrat(l)=nbd
          if(l.lt.lz) npbrat(l+1)=npbrat(l)+nbdrat(l)
          do nf=1,nbd
             lbdrat(npbrat(l)+nf)=lbd(nf)
          end do
       end do
    end if
!
    write(imp,'(/,"==>atcaldis3: rattachement des domaines aux parois")')
    do l=1,lzx
       if0=npbrat(l)+1
       if1=npbrat(l)+nbdrat(l)
    end do
!
!     ----------------------------------------------------------
!
!       initialisation des distances a une valeur infinie
!
    do n=1,ndimctbx
       dist(n) =reelmx
       mnpar(n)=0
    enddo
!
    if(klecdis.eq.1) then
!       lecture des distances
!
       call at_lecdist( &
            ldismx, &
            dist,mnpar)
!
    else
       ldismx=0
    end if
!     ----------------------------------------------------------
!
!     boucle sur les frontieres a traiter (parois)
!
    do mf=1,nbd
       nfbi=lbd(mf)
!
!       calcul des indices pour deplacement "parallelement" a la paroi
       call atindnor( &
            nfbi, &
            m10,m20,m30, &
            m1min,m1max,m2min,m2max,m3min,m3max, &
            dm1,dm2,dm3,isens3)
!
!       calcul des centres des facettes formant les parois
       call atccfp( &
            x,y,z, &
            xpar,ypar,zpar, &
            ncin, &
            nfbi,dm1,dm2,dm3,isens3)
!
!       relation entre les indices "m1 m2" des facettes paroi et leur numero
!
       call at_dlist( &
            x,y,z, &
            xpar,ypar,zpar, &
            ncin, &
            m1tb,m2tb,nfrtb, &
            nfbi)
!
    enddo
!
!     -------------------------------------------------------------
!
!     Calcul de la distance de toutes les cellules a la facette
!     la plus proche. Toutes les facettes de toutes les parois
!     sont esssayees.
!     boucle sur les domaines
!
    do l=1,lzx
       if(l.gt.ldismx) then
!
!         calcul des centres des cellules
          call atccc( &
               x,y,z, &
               xcc,ycc,zcc, &
               l)
!
!         calcul distance des cellules d'un domaine aux parois dans le
!         maillage grossier par exploration systematique des frontieres
!
!         igr=5
!         jgr=3
!         kgr=1
          if(igr(l).ne.1 .or. jgr(l).ne.1 .or. kgr(l).ne.1) then
             call at_grdist( &
                  igr(l),jgr(l),kgr(l), &
                  x,y,z, &
                  xpar,ypar,zpar, &
                  xcc,ycc,zcc,dist2, &
                  dist,mnpar, &
                  m1tb,m2tb,nfrtb, &
                  l)
!
!           retour au maillage fin en se limitant a l'exploration de
!           portions reduites de frontieres quand une cellule du maillage
!           grossier est completement rattachee a un domaine reduit de
!           frontiere
!
!           fraction la plus grande autorisee sur les indices "m1,m2"
!           pour decider qu'une portion de frontiere est contigue.
!           On peut mettre 0 pour les frontieres qui ne se referment pas.
!           Une frontiere se referme si les indices min et max correspondent
!           physiquement a un meme noeud. Exemple d'un bord de fuite pour
!           un profil defini en une seule frontiere, du BF au BF en passant
!           par le bord d'attaque.
!           raptat=0.5
!           ==========

             call at_fidist( &
                  igr(l),jgr(l),kgr(l),raptat(l), &
                  x,y,z, &
                  xpar,ypar,zpar, &
                  xcc,ycc,zcc,dist2, &
                  dist,mnpar, &
                  m1tb,m2tb,nfrtb, &
                  l)
          else
!           recherche systematique avec liste
!
             call atdist_3( &
                  x,y,z, &
                  xpar,ypar,zpar, &
                  xcc,ycc,zcc,dist2, &
                  dist,mnpar, &
                  l)
          end if
          if(kecrdis.eq.1) then
!           ecriture disque des distances (fichiers separes "fdist_l")
             call at_ecrdist( &
                  l,          &
                  dist,mnpar)
          end if
       end if
!
    enddo
!
    deallocate(raptat)
    return
  end subroutine atcaldis3
end module mod_atcaldis3
