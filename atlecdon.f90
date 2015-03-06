      subroutine atlecdon
!
!     ACT
!_A   Lecture des donnees pour les modeles de turbulence
!
!
!     OUT
!_O    kcaldis    : com int              ; cle calcul distance aux parois
!_O    kecrdis    : com int              ; cle ecriture disque distance parois
!_O    klecdis    : com int              ; cle lecture  disque distance parois
!_O    equatt     : com char*7           ; code modele de turbulence
!_O    kcutke     : com int              ; cle limiteurs variables turb (k-e)
!_O    kfludis    : com int              ; cle calcul flux dissipatifs
!_O    ksecmb     : com int              ; cle calcul second membre
!_O    kcmut      : com int              ; cle calcul mut
!_O    kclkep     : com int              ; cle comdition paroi
!_O    kinke      : com int              ; cle initialisation k-e
!_O    ncytuke0   : com int              ; nombre de cycles initiaux d'integ.
!_O                                        de k-e uniquement (u fixe)
!_O    ncycke     : com int              ; frequence d'integration des
!_O                                        equations de transport
!_O    imxclko    : com int              ; nb pts pres paroi pour ponderation
!_O                                        sur omega (k-omega)
!_cc    kutau : cle pour modeles de turbulence qui utilisent la vitesse de
!cc             frottement U_tau (Chien, Wilcox bas Reynolds, Wilcox rugueux)
!cc     kparoi : cle pour modeles de turbulence qui presentent une condition a
!cc              la paroi particuliere (Wilcox, Menter, Wilcox bas reynolds)
!_O    epspid     : com real             ; precision sur Pi pour calcul delta
!_O    epstaud    : com real             ; precision sur tau pour calcul delta
!_O    epsvord    : com real             ; precision sur vort pour calcul delta
!
!***********************************************************************
!
        use maillage
        use modeleturb
        use chainecarac
        use sortiefichier
        use schemanum
!
!-----------------------------------------------------------------------
!
      character *80 ligne
!
!      write(imp,'(/,''==>atlecdon: lecture atlecdon et '',
!     &''initialisation modeles de turbulence'')')
!
      open(99,file='fatdon',form='formatted',err=100)
!
      lig=1
      read(99,  300,err=200)ligne
!      write(imp,300)ligne
      lig=lig+1
      read(99,  300,err=200)ligne
!      write(imp,300)ligne
      lig=lig+1
!
      read(99,  300,err=200)ligne
!      write(imp,300)ligne
      backspace(99)
      read(99,*,err=200)kcaldis
!
      lig=lig+1
      read(99,  300,err=200)ligne
!      write(imp,300)ligne
      backspace(99)
      read(99,*,err=200)kecrdis
!
      lig=lig+1
      read(99,  300,err=200)ligne
!      write(imp,300)ligne
      backspace(99)
      read(99,*,err=200)klecdis
!
      lig=lig+1
      read(99,  300,err=200)ligne
!      write(imp,300)ligne
      equatt=ligne(1:7)
!
      lig=lig+1
      read(99,  300,err=200)ligne
!      write(imp,300)ligne
      backspace(99)
      read(99,*,err=200)ncytuke0
!
      lig=lig+1
      read(99,  300,err=200)ligne
!      write(imp,300)ligne
      backspace(99)
      read(99,*,err=200)ncycke
!
      ierrdis=0
      keasm=0
!
!---------------------------------------------------------------------------
!       Modele k-eps de Jones et Launder (1972)
!       2JLS -> SST    ; 2JLR -> realisable   ; 2JLM -> Cmu variable
!       2JLC -> compressible   ;  2JLL -> SAS 
!---------------------------------------------------------------------------
!
      if(equatt(1:3).eq.'2JL') then
!
        kcutke =1
        kfludis=1
        ksecmb =1
        kcmut  =1
        kclkep =1
        kinke  =1
        kparoi =0
        kutau  =0
        kesst  =0
        cke1   =1.55
        cke2   =2.
        alfak  =1.
        alfae  =1.3
        cmu    =0.09
!       modele Jones Launder avec correction SST
        if(equatt(4:4).eq.'S') then
          kesst=1
        endif
        if(equatt(1:4).eq.'2JL2') then
!
!         Utilisation de la fonction de raccord de Menter pour passer
!         a la fonction f_mu de Smith dans les regions externes et les sillages
          kcmut  =6
          cmukl  =0.09
          cklb1  =18.
          xkappa =0.41
        endif
!
        if(equatt(4:4).ne.' ' .and. equatt(4:4).ne.'R' .and. &
           equatt(4:4).ne.'S' .and. equatt(4:4).ne.'M' .and. &
           equatt(4:4).ne.'C' .and. equatt(4:4).ne.'L') then
          write(imp,'(/,''!!!atlecdon: modele'',a7, &
          '' non prevu'')')equatt
          stop
        endif
!
!---------------------------------------------------------------------------
!       Modele k-eps de  Launder Sharma
!---------------------------------------------------------------------------
!
      elseif(equatt(1:3).eq.'2LS') then
!
        kcutke =1
!       decouplage limiteur sur k et epsilon pres paroi
!       kcutke =3
        kfludis=1
        ksecmb =1
        kcmut  =1
        kclkep =1
        kinke  =1
        kparoi=0
        kutau=0
!
        if(equatt(1:4).eq.'2LS2') then
!
!         Utilisation de la fonction de raccord de Menter pour passer
!         a la fonction f_mu de Smith dans les regions externes et les sillages
!
          kcmut  =6
          cmukl  =0.09
          cklb1  =18.
          xkappa =0.41
        endif
!
!---------------------------------------------------------------------------
!       Modele de Spalart-Allmaras
!       1SAD -> DES Spalart  ; 1SAL -> SAS
!---------------------------------------------------------------------------
!
      elseif(equatt(1:3).eq.'1SA') then
!
        kcutke =2
        kfludis=2
        ksecmb =2
        kcmut  =2
        kclkep =1
        kinke  =2
        kparoi=0
        kutau=0
        cb1 = 0.1355
        sigma=2./3.
        cb2=0.622
        kappa=0.41
        cw1=cb1/kappa**2+(1+cb2)/sigma
        cw2=0.3
        cw3=2.
        cv1=7.1
        ct1=1.
        ct2=2.
        ct3=1.1
        ct4=2.
!
        if(equatt(4:4).ne.' ' .and. equatt(4:4).ne.'D' .and. equatt(4:4).ne.'L') then
          write(imp,'(/,''!!!atlecdon: modele'',a7, &
          '' non prevu'')')equatt
          stop
        endif
        if(kcaldis.eq.0 .and. klecdis.eq.0) ierrdis=1
!
!------------------------------------------------------------------
!       Modele k-l de Smith
!       2SmS -> SST  ; 2SmR -> realisable ; 2SmL -> SAS
!--------------------------------------------------------------------
!
      elseif(equatt(1:3).eq.'2Sm') then
!
        kinke  =3
        kcutke =3
        kfludis=1
        ksecmb =3
        kcmut  =3
        kclkep =3
        kparoi=0
        kutau=0
        cmukl  =0.09
        cklb1  =18.
        ckle2  =1.2
        sigmak =1.43
        sigmal =1.43
        xkappa =0.41
!
        if(equatt(4:4).ne.' ' .and. equatt(4:4).ne.'L' .and. &
           equatt(4:4).ne.'S' .and. equatt(4:4).ne.'R') then
          write(imp,'(/,''!!!atlecdon: modele'',a7, &
          '' non prevu'')')equatt
          stop
        endif
!
        if(kcaldis.eq.0 .and. klecdis.eq.0) ierrdis=1
!
!------------------------------------------------------------------
!       Modele k-w de Wilcox et Menter (BSL et SST)
!--------------------------------------------------------------------
!
      elseif(equatt(1:3).eq.'2MT' .or. equatt(1:3).eq.'2WL') then
!
        kinke  =4
        kcutke =3
        kfludis=4
        ksecmb =4
        kcmut  =4
        kparoi=1
        kutau=0
        betae  =0.09
        okappa =0.41
        sigme1 =0.5
        sigma1 =0.5
        beta1  =0.075
        wsig1  =0.
        sigme2 =1.
        sigma2 =0.856
        beta2  =0.0828
        wsig2  =0.856
!       komsst=0->Menter de base  komsst=1->Menter SST
        komsst =0
!
!       imxclko : nombre de cellules pres fontiere pour
!                               ponderation entre omega calcule et impose
!                               Doit etre superieur ou egal 1
!       imxclko=1
        imxclko=4
        if(kcaldis.eq.0 .and. klecdis.eq.0) ierrdis=1
        if(equatt(1:3).eq.'2WL') then
          kfracom=0
        elseif(equatt(1:3).eq.'2MT') then
          kfracom=1
          sigme2=0.85
          if(equatt(4:4).eq.'S') then
           komsst=1
          endif
        endif
        if(equatt(4:4).eq.'G') keasm=1
!
!------------------------------------------------------------------
!       Modele k-eps de Chieng
!--------------------------------------------------------------------
!
      else if(equatt(1:3).eq.'2CH') then
!
        kcutke =1
        kfludis=1
        ksecmb =5
        kcmut  =5
        kclkep =1
        kinke  =1
        kparoi=0
        kutau=1
!       utaumin : valeur mini de u_tau
!       utau / a_i = u_0/ a_i sqrt( Cf0/2 )
!
        utaumin=2.e-4
!
        if(kcaldis.eq.0 .and. klecdis.eq.0) ierrdis=1
!
!--------------------------------------------------------------------------
!       Modele k-eps RNG de Yakhot, Orzag, Thangam, Gatski et Speziale
!----------------------------------------------------------------------------
!
      else if(equatt(1:3).eq.'2RN') then
!
!       kinke  =6         ! pour interdire l'initialisation
        kinke  =3         ! initialisation comme Smith
        kcutke =1
        kfludis=6
        ksecmb =6
        kcmut  =7
        kclkep =3
        kparoi=0
        kutau=0
        cmukl  =0.09
        cklb1  =18.
        ckle2  =1.2
        sigmak =1.43
        sigmal =1.43
        xkappa =0.41
!         rtrac  : valeur centrale de R_t pour raccord des modeles
!         drtrac : amplitude du raccord en R_t
!                  F1 varie de         1       a      0      pour
!                  R_t variant de rtrac-drtrac a rtrac+drtrac
!         ncycrac: frequence de recalcul de la fonction F_1
!         naprng : application correction RNG sur C_epsilon_1
!                  0->non  1->oui
        rtrac  =400.
        drtrac =50.
        ncycrac=100
        naprng =1
!
!--------------------------------------------------------------------------
!       Modele longueur de melange de Michel
!----------------------------------------------------------------------------
!
      elseif(equatt(1:3).eq.'0LM') then
!
!      rapvisq: valeur de tau/tau_max pour definir la frontiere de la couche limite
!      rapvisq=0.02
        rapvisq=0.005
!
!--------------------------------------------------------------------------
!       Modele k-w de Wilcox compressible
!----------------------------------------------------------------------------
!
      elseif(equatt(1:3).eq.'2KO') then
!
        kinke  =4
        kcutke =3
        kfludis=5
        ksecmb =8
        kcmut  =8
        kparoi =1
        kutau  =0
        betas  =0.09
        beta   =0.075
        sigmk  =0.5
        sigmw  =0.5
        sigmd  =0.65
        kwsst  =0

        if(equatt(4:4).eq.'S') then
!       version SST
         kwsst=1
        endif
!
        if(equatt(4:4).ne.' ' .and. equatt(4:4).ne.'R' .and. &
           equatt(4:4).ne.'S' ) then
          write(imp,'(/,''!!!atlecdon: modele'',a7, &
          '' non prevu'')')equatt
          stop
        endif

!
!--------------------------------------------------------------------------
!       Modele k-e-v2 de Durbin
!----------------------------------------------------------------------------
!
      elseif(equatt(1:3).eq.'3KEV') then
!
        kinke  =4
        kcutke =3
        kfludis=5
        ksecmb =8
        kcmut  =8
        kparoi =1
        kutau  =0
        ccmu   =0.09
        cc1    =1.4
        cc2    =0.3
        ce2    =1.9
        cgl    =0.25
        ceta   =85.
        sigk   =2./3.
        sige   =0.5
!
!--------------------------------------------------------------------------
!       CALCUL EULER
!----------------------------------------------------------------------------
!
      else if(equatt(1:3).eq.'EU') then
!
      else
         write(imp,'(/,''!!!atlecdon: modele'',a7, &
         '' non prevu'')')equatt
         stop
      endif
!
      if(klecdis.eq.1 .and. (kcaldis.eq.1 .or. kcaldis.eq.2) ) then
        write(imp,'(/,''!!!atlecdon: lecture et calcul distance '', &
        ''incompatible sauf si kcaldis=3'')')
        if(ierrdis.eq.0) stop
      endif
      if(ierrdis.ne.0) then
        write(imp,'(/,''!!!atlecdon: il faut un calcul de '', &
        ''distance avec ce modele'')')
        stop
      endif
!
!     --------------------------------------------------------
!     lecture des precisions pour calcul de delta
!
      muttrav=0
      epspid =-1.e+10
      epstaud=-1.e+10
      epsvord=-1.e+10
      rap_h0=0.
      kamort=0
      do
        read(99,300,err=200,end=51)ligne
        lig=lig+1
        if(ligne(1:5).eq.'preci') then
!         lecture des precisions pour calcul de delta
          read(99,300,err=200,end=201)ligne
          backspace(99)
          read(99,*,err=200)epspid
          lig=lig+1
          read(99,300,err=200,end=201)ligne
          backspace(99)
          read(99,*,err=200)epstaud
          lig=lig+1
          read(99,300,err=200,end=201)ligne
          backspace(99)
          read(99,*,err=200)epsvord
          lig=lig+1
        elseif(ligne(1:5).eq.'travs') then
          read(99,  300,err=200,end=201)ligne
          backspace(99)
          read(99,*,err=200)muttrav
          lig=lig+1
        elseif(ligne(1:10).eq.'ACOUSTIQUE') then
          read(99,300,err=200,end=201)ligne
          backspace(99)
          read(99,*,err=200)lacou
          lig=lig+1
          read(99,*,err=200)x0,y0,z0,freq,cga
          lig=lig+1
        elseif(ligne(1:10).eq.'SORTIEPLOT') then
          read(99,300,err=200,end=201)ligne
          backspace(99)
          read(99,*,err=200)lsortie
          lig=lig+1
          read(99,*,err=200)nfreq
          lig=lig+1
        endif
       enddo
!
   51 continue
!
      close(99)
!
      return
!
  100 continue
!     erreur ouverture fichier
      write(imp,'(/,''!!!atlecdon: probleme ouverture fichier '', &
      ''atlecdon'')')
      stop
  200 continue
!     erreur lecture fichier
      write(imp,'(/,''!!!atlecdon: probleme lecture fichier '', &
      ''fatdon ligne'',i2)')lig
      stop
  201 continue
!     fin prematuree lecture fichier
      write(imp,'(/,''!!!atlecdon: fin fichier fatdon. '', &
      ''manquer sequence precision pour delta'',i2)')lig
      stop
!
  300 format(a80)
!
      return
      end
