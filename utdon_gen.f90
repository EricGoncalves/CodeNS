module mod_utdon_gen
  implicit none
contains
  subroutine utdon_gen( &
       config,cl,x,y,z,omg1, &
       ncbd,v, &
       nxn,nyn,nzn)
!
!***********************************************************************
!
!     ACT
!_A    Lecture de donnees specifiques du cas de calcul et stockees dans
!_A    des common pour etre accessibles pendant le calcul.
!_A
!_A    Lecture de donnees concernant le calcul de grandeurs globales
!_A      au cours des iterations (utit) et en fin de calcul (utsor), soit
!_A      point (xref, yref, zref), surface (sref) et longueur (xlref)
!_A      de reference,
!_A      nombre de frontieres sur lesquelles se feront les integrations,
!_A      numeros des frontieres concernees.
!
!_I    ip11       : arg int              ; dim, nbr max de cellules de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip41       : arg int              ; dim, nbr max de pts de ttes les front
!_I    ip60       : arg int              ; dim, nbr max d'equations
!_I    kimp       : arg int              ; niveau de sortie sur unite logi imp
!_I    don1       : com int              ; unite logiq, donnees utilisateur
!_I    ta1        : com real             ; etat de reference utilisateur
!_I                                        adimensionne, temperature d'arret
!_I    pa1        : com real             ; pression d'arret de l'etat
!_I                                        de reference utilisateur adimensionne
!
!     OUT
!_O    tpar       : com real             ; temperature paroi
!_O    xref       : com real             ; x du pt de ref pour calcul grandeurs
!_O                                        globales
!_O    yref       : com real             ; y du pt de ref pour calcul grandeurs
!_O                                        globales
!_O    zref       : com real             ; z du pt de ref pour calcul grandeurs
!_O                                        globales
!_O    sref       : com real             ; surface de ref pour calcul grandeurs
!_O                                        globales
!_O    xlref      : com real             ; longueur de ref pour calcul grandeurs
!_O                                        globales
!_O    kvglo      : com int              ; cle calcul des grandeurs globales
!_O    nbfll      : com int              ; nb de front d'integration des
!_O                                        grandeurs globales
!_O    nmfint     : com int (mtb       ) ; no des front d'integration des
!_O                                        grandeurs globales
!     LOC
!_L    mdimtb     : com int              ; dim, nbr max de pts de ttes les front
!
!***********************************************************************
!
    use para_var
    use para_fige
    use boundary,only :new2old_f
    use sortiefichier
    use maillage
    use constantes
    use modeleturb
    use schemanum
    use definition
    implicit none
    integer          ::   idefconf,  idefxref,      ierr,     ligne,      mflu
    integer          ::         nb,ncbd(ip41),l1
    double precision ::    nxn(ip42),   nyn(ip42),   nzn(ip42),        omg1,          p2
    double precision ::          rpi,         rti,        tpar,v(ip11,ip60),     x(ip21)
    double precision ::      y(ip21),     z(ip21)
!
!**********************************************************************
!
    character(len=4 ) :: config,cl(mtb)
    character(len=80) ::  cmtlec
!
!
    open(don1,file='fdon1',err=101)
!
!     initialisation des cles globales de controles
!
    kvglo=0
    nbfll=0
!
!     initialisation des cles locales de controles de la lecture de FDON1
!
    ierr=0
    idefxref=0
    idefconf=0
    ligne=0
!
!     initialisation des grandeurs obilgatoires pour verification en
!     fin de lecture
!
    kditur=intmx
!
    cmtlec=""
    do while(cmtlec(1:4 ).ne.'stop' .or.cmtlec(1:4 ).ne.'STOP')
!
       read(don1,7000,err=99,end=98) cmtlec
       ligne=ligne+1
       if(cmtlec(1:1) .eq.'#') then
          continue
!
       elseif(cmtlec(1:10).eq.'schema num' .or. &
            cmtlec(1:10).eq.'SCHEMA NUM') then
!
          write(imp,'(/,"SCHEMA NUMerique:")')
!
!         donnees pour schema numerique d'integration de k-epsilon
!         schema centre de Jameson ou decentre de Roe
!
          read(don1,7000) cmtlec
          ligne=ligne+1
          read(don1,*,err=99) kditur
          ligne=ligne+1
          write(imp,6803) kditur
          read(don1,7000) cmtlec
          ligne=ligne+1
          read(don1,*,err=99) klroe
          ligne=ligne+1
          write(imp,6802) klroe
          read(don1,7000) cmtlec
          ligne=ligne+1
          read(don1,*,err=99) epsroe
          ligne=ligne+1
          write(imp,6800) epsroe
!
       elseif(cmtlec(1:10).eq.'configurat' .or. &
            cmtlec(1:10).eq.'CONFIGURAT') then
!
          if(kimp.gt.1) then
             write(imp,'(/,"CONFIGURATion >>>")')
          endif
!
!         Donnees propres a la configuration : exemple
!         0.72092786       : P2/PI0
!         1.00             : TPAR/TI0
!         1.00  1.00       : RPI   RTI
!         0.70  0.00  0.00 : RMACH ALPHA BETA
!
          read(don1,*,err=99) p2
          ligne=ligne+1
          read(don1,*,err=99) tpar
          ligne=ligne+1
          read(don1,*,err=99) rpi,rti
          ligne=ligne+1
          read(don1,*,err=99) rm0,al0,be0
          ligne=ligne+1
          idefconf=1
!
!         passage en coordonnees aerodynamiques dans:
!         utsor.f utsorfr.f utit.f utitfr.f
          alpha0=al0
          beta0 =be0
!
          if(kimp.gt.1) then
             write(imp,'(9x," p2=",1pe11.3,3x,"tpar=",1pe11.3,' &
                  //'3x,"  rpi=",1pe11.3,3x,"rti=",1pe11.3)') p2,tpar,rpi,rti
             write(imp,'(9x,"rm0=",1pe11.3,3x," al0=",1pe11.3,' &
                  //'3x,"beta0=",1pe11.3)')rm0,al0,be0
          endif
!
          p2=p2*pa1
          tpar=tpar*ta1
!
       elseif(cmtlec(1:10).eq.'surfaces i' .or. &
            cmtlec(1:10).eq.'SURFACES I') then
!
          if(kimp.gt.1) then
             write(imp,'(/,"SURFACES I >>> ")')
          end if
!
          kvglo=1
          idefxref=1
          read(don1,*,err=99) xref,yref,zref,sref,xlref
          ligne=ligne+1
          read(don1,*,err=99) p0spi0,q0spi0,v0
          ligne=ligne+1
!
          if(kimp.gt.1) then
             write(imp,'("integration pression et frottement ",' &
                  //'"sur les parois. kvglo=",i3)') kvglo
             write(imp,'("grandeurs de reference :",/,' &
                  //'5x,"Xref=",1pe12.4,"   Yref=",1pe12.4,"  Zref=",1pe12.4,/,' &
                  //'5x,"Sref=",1pe12.4,"  XLref=",1pe12.4)') &
                  xref,yref,zref, sref,xlref
             write(imp,'(3x,"Po/Pio=",1pe12.4,1x,"Qo/Pio=",1pe12.4,' &
                  //'4x,"Vo=",1pe12.4)')p0spi0,q0spi0,v0
          endif
!
       elseif(cmtlec(1:10).eq.'liste surf' .or. &
            cmtlec(1:10).eq.'LISTE SURF') then
!         liste d'un nombre mon fixe de numeros de parois
!
          if(kimp.gt.1) then
             write(imp,'(/,"LISTE SURFaces >>> pour calcul epaisseurs et Cx_f")')
          endif
!
          if(kvglo.eq.0) then
             write(imp,'(/,"!!!utdon_gen: il faut mettre : SURFACES Integration avant LISTE SURFaces")')
             ierr=ierr+1
          endif
!
          do
             read(don1,7000,err=99,end=98) cmtlec
             ligne=ligne+1
             if(cmtlec(1:3).ne.'fin' .and. cmtlec(1:3).ne.'FIN') then
                backspace(don1)
                read(don1,*,err=99) mflu
               do l1=1,mtb
                if (new2old_f(l1)==mflu) then
                  if(l1.ge.1 .and. l1.le.mtbx) then
                     if(nbfll.lt.mtb) then
                        nbfll=nbfll+1
                        nmfint(nbfll)=l1
                     else
  !                 suite de la liste non prise en compte
                        write(imp,'(/,"!!!utdon_gen: trop de surfaces ")')
                     endif
                  else
  !               mauvais numero de surface
                     write(imp,'(/,"!!!utdon_gen: numero de surface incorrect: mf=",i4)')mflu
                  endif
               endif
               enddo
             else !fin sequence, continuer la lecture generale des mots-c
                if(kimp.gt.1) then
                   write(imp,'( (20i4) )')(nmfint(nb),nb=1,nbfll)
                endif
                exit
             endif
          enddo
!
       elseif(cmtlec(1:3).eq.'vrt' .or. cmtlec(1:3).eq.'VRT') then
          read(don1,*,err=99) vrtcz
          ligne=ligne+1
          if(kimp.gt.1) then
             write(imp,'(/,"VRT  >>> correction de vorticite vrtcz=",1pe11.4)')vrtcz
             write(imp,'(/,9x,"le profil doit etre dans le plan x-z")')
             if(idefxref.ne.1 .or. idefconf.ne.1) then
                write(imp,'(/,"!!!utdon_gen: il faut VRT apres ",' &
                     //'"SURFACES Integration pour definition de ",' &
                     //'"xlref, xref et zref",/,32x,' &
                     //'"CONFIGURATion pour rmach, alpha et beta=0")')
                stop
             endif
          endif
!
          vrtmac=rm0
          vrtalp=al0
          vrtlre=xlref
          vrtxre=xref
          vrtzre=zref
!
!         pour un profil de corde 1: xref=.25, zref=0., xlref=1. place le
!         vortex a 25% de la corde
!
!         EN SOMME, il faut fournir a clvrt l'origine du vortex, ce qui est fait
!         par vrtlre, vrtxre, vrtzre et le cz pour en calculer l'intensite.
!         Le cz passe par utit au cours des iterations et par utdon en cas
!         de reprise.
!
       endif
    enddo
!
98  continue
!
!     VERIFICATIONS
!
    if(kditur.eq.intmx) then
       ierr=ierr+1
       write(imp,'(/,"!!!utdon_gen: schema numerique non defini mot-cle SCHEMA NUM")')
    end if
!
    if(ierr.ne.0) then
       write(imp,'(/,"!!!utdon_gen: ierr=",i3," STOP")')ierr
       stop
    end if
    if (kimp.gt.1) then
       write(imp,'(/,70("-"),/)')
    end if
!
!
7000 format(a)
7100 format(5f10.0)
7200 format(16i5)
!
6800 format(7x,'epsroe = ' ,f10.4)
6802 format(7x,'klroe  = ' ,i5)
6803 format(7x,'kditur = ' ,i5)
!
    return
!
99  continue
    write(imp,'(/,"!!!utdon_gen: erreur lecture fdon1 ligne :",i4)') &
         ligne
    stop
101 continue
    write(imp,'(/,"!!!utdon_gen: erreur ouverture fdon1")')
    stop
!
    return
  end subroutine utdon_gen
end module mod_utdon_gen
