module mod_at_lecopt
  implicit none
contains
  subroutine at_lecopt( &
       igr,jgr,kgr,raptat,idefaut)
!
!***********************************************************************
!
!     ACT
!_A   Lectures des donnees pour calcul optimise des distances
!
!     INP
!_I    imp        : com int              ; unite logique fichier fimp
!
!     OUT
!_O    nbdrat     : com int (lz        ) ; nb parois rattachees au domaine l
!_O    lbdrat     : com int (mtb       ) ; no interne des parois rattachees
!_O                                        au domaine l
!_O    npbrat     : com int (lz        ) ; pointeur fin liste frontieres a
!_O                                        traiter pour domaine precedent
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use maillage
    implicit none
    integer          :: idefaut,   ierr,igr(lz),   igrl,jgr(lz)
    integer          ::    jgrl,kgr(lz),   kgrl,      l,    lig
    integer          ::     llu,     nf
    double precision :: raptat(mtb),    raptatl
!
!-----------------------------------------------------------------------
!
    character(len=80) ligne
!
    write(imp,'(/,"==>at_lecopt: lecture relations domaines-paroi pour calcul optimise de la distance",/)')
!
    open(99,file='fatdon',form='formatted',err=100)
!
    lig=0
    do
!       boucle sur les lignes de "fatdon"
!
       read(99,  300,err=20)ligne
       write(imp,300)ligne
       lig=lig+1
       if(ligne(1:6).eq.'optimi' .or. ligne(1:6).eq.'OPTIMI') then
!
!         lecture des informations obligatoires pour l'optimisation
!
          read(99,*,err=11)igrl
          lig=lig+1
          read(99,*,err=11)jgrl
          lig=lig+1
          read(99,*,err=11)kgrl
          lig=lig+1
          read(99,*,err=11)raptatl
          lig=lig+1
!
!         valeurs par defaut
!$OMP SIMD
          do l=1,lzx
             igr(l)=igrl
             jgr(l)=jgrl
             kgr(l)=kgrl
             raptat(l)=raptatl
          enddo
!
!         modification eventuelle pour certains domaines
          read(99,300,err=20)ligne
          lig=lig+1
          if(ligne(1:6).ne.'modifi' .and. ligne(1:6).ne.'MODIFI') then
             write(imp,'("!!!at_lecopt: manque MODIFICATION")')
             stop
          endif
          do
             read(99,300,err=21)ligne
             lig=lig+1
             if(ligne(1:3).ne.'fin' .and. ligne(1:3).ne.'FIN') then
                backspace 99
                read(99,*,err=22)l,igr(l),jgr(l),kgr(l),raptat(l)
             else
                exit
             endif
          enddo
!
          read(99,300,err=24,end=31)ligne
          lig=lig+1
          if(ligne(1:5).eq.'liste' .or. ligne(1:5).eq.'LISTE') then
!         lecture facultative des listes de parois rattachees a chaque domaine
             idefaut=0
             ierr=0
             npbrat(1)=0
             do l=1,lzx
                read(99,*,err=25)llu,nbdrat(l)
                if(abs(llu).ne.l) then
                   write(imp,'("!!!at_lecopt: erreur sur numero bloc l=",i4,3x,"l lu=",i4)')
                   ierr=ierr+1
                endif
                if(l.lt.lz) then
                   npbrat(l+1)=npbrat(l)+nbdrat(l)
                endif
                read(99,*,err=25)(lbdrat(npbrat(l)+nf),nf=1,nbdrat(l))
             enddo
             if(ierr.ne.0) then
!             detection erreurs lecture. Arret.
                write(imp,'("!!!at_lecopt: detection erruers. ARRET")')
                stop
             endif
             exit  !fin des lectures
          else
             idefaut=1
          endif ! fin lecture facultative listes parois attachees a chaq
!
          exit
       endif ! pas de mot cle OPTIMI trouve, lire ligne suivante
    enddo
!
31  continue
    close(99)
    return
!
20  continue
    write(imp,'("!!!at_lecopt: fin de fichier prematuree")')
    stop
21  continue
    write(imp,'("!!!at_lecopt: pas de FIN sequence MODIFICATION")')
    stop
22  continue
    write(imp,'("!!!at_lecopt: erreur lecture donnees MODIFICATION")')
    stop
23  continue
    write(imp,'("!!!at_lecopt: pas de FIN sequence MODIFICATION")')
    stop
24  continue
    write(imp,'("!!!at_lecopt: erreur recherche LISTE")')
    stop
25  continue
    write(imp,'("!!!at_lecopt: erreur lecture donnees LISTE")')
    stop
10  continue
    write(imp,'("!!!at_lecopt: erreur lecture mot-cle OPTIMI")')
    stop
11  continue
    write(imp,'("!!!at_lecopt: erreur lecture igr,jgr,kgr,raptat")')
    stop
12  continue
    write(imp,'("!!!at_lecopt: erreur lecture mot-cle LISTE")')
    stop
100 continue
!     erreur ouverture fichier
    write(imp,'(/,''!!!at_lecopt: probleme ouverture fichier '',''atlecdon'')')
    stop
!
300 format(a80)
!
    return
  end subroutine at_lecopt
end module mod_at_lecopt
