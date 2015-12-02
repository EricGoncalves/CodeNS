module mod_inbdc
  implicit none
contains
  subroutine inbdc( &
       exs1,exs2, &
       x,y,z, &
       ncbd,ncin,mnc, &
       krr,mfbea,mfbeb,kibdc,epsmsh, &
       iba,jba,kba,tvi,tvj,tvk)
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des donnees necessaires aux frontieres
!_A    coincidentes, principalement pour chaque point l'indice de la maille
!_A    coincidente.
!
!     INP
!_I    ip00       : arg int              ; dim, nbr max de pts d'un dom
!_I                                        (pts fictifs inclus)
!_I    ip11       : arg int              ; dim, nbr max de cellules de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip12       : arg int              ; dim, nbr max de cellules de tous les
!_I                                        dom (pts fictifs inclus)
!_I                                        pour Navier-Stokes
!_I    ip13       : arg int              ; dim, nbr max de cellules de tous les
!_I                                        dom (pts fictifs inclus)
!_I                                        pour k-eps
!_I    ip21       : arg int              ; dim, nbr max de noeuds de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip31       : arg int              ; dim, nbr max de facettes de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip40       : arg int              ; dim, nbr max de pts d'une front
!_I    ip41       : arg int              ; dim, nbr max de pts de ttes les front
!_I    ip42       : arg int              ; dim, nbr max de pts de ttes les front
!_I                                        a normales stockees
!_I    ip43       : arg int              ; dim, nbr max de pts de ttes les front
!_I                                        coincidentes
!_I    ip44       : arg int              ; dim, nbr max de pts de ttes les front
!_I                                        recouvertes
!_I    ip60       : arg int              ; dim, nbr max d'equations
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    krr        : arg int              ; cle info sur front coinc
!_I    epsmsh     : arg real             ; dist max entre deux pts confondus
!_I    kibdc      : arg int              ; cle initialisation front coinc
!_I    exs1       : arg real             ; premier coef d'interpolation a l'
!_I                                        l'ordre 0 du couple de coef exs1,exs2
!_I    exs2       : arg real             ; deuxieme coef d'interpolation a l'
!_I                                        l'ordre 0 du couple de coef exs1,exs2
!_I    la         : arg int              ; numero de domaine
!_I    mfbea      : arg int              ; numero externe de front
!_I    lb         : arg int              ; numero du domaine coincident
!_I    mfbeb      : arg int              ; numero externe de front coincidente
!_I    iba        : arg int              ; ind i du pt coinc au pt d'indices
!_I                                        min de la frontiere
!_I    jba        : arg int              ; ind j du pt coinc au pt d'indices
!_I                                        min de la frontiere
!_I    kba        : arg int              ; ind k du pt coinc au pt d'indices
!_I                                        min de la frontiere
!_I    tvi        : arg char             ; sens de variation sur la front coinc
!_I                                        pour une var de l'ind i sur la front
!_I    tvj        : arg char             ; sens de variation sur la front coinc
!_I                                        pour une var de l'ind j sur la front
!_I    tvk        : arg char             ; sens de variation sur la front coinc
!_I                                        pour une var de l'ind k sur la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    indfl      : com char(mtb       ) ; type de plan de la frontiere
!_I    cl         : com char(mtb       ) ; type de cond lim a appliquer
!_I    equat      : com char             ; type d'equations modelisant l'ecoule-
!_I                                        ment
!_I    kfc        : com int              ; unite logiq, tableaux front coinc
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une frontiere
!_I    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une frontiere
!_I    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une frontiere
!
!     OUT
!_O    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_O                                        cellule coincidente
!_O    mpc        : com int (mtt       ) ; pointeur fin de front precedente
!_O                                        dans tableaux front coinc
!_O    nfbc       : com int (mtb       ) ; numero dans numerotation interne
!_O                                        d'une frontiere coincidente
!_O    ndcc       : com int (mtb       ) ; numero du dom coicident
!_O    mdnc       : com int (mtt       ) ; saut d'ind entre pt front fictif et
!_O                                        pt interieur pour la front coinc
!
!     I/O
!_/    mtcx       : com int              ; nbr total de frontieres coincidentes
!_/    mdimtcx    : com int              ; nbr de pts de ttes les front
!_/                                        coincidentes
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use sortiefichier
    use chainecarac
    use mod_initcs
    use mod_extmhg
    use mod_mpi
    implicit none
    integer          ::        ia1,       ia2,       ib1,       ib2,       iba
    integer          ::       ibam,       img,      imgi,      imgj,      imgk
    integer          ::        ja1,       ja2,       jb1,       jb2,       jba
    integer          ::       jbam,       ka1,       ka2,       kb1,       kb2
    integer          ::        kba,      kbam,     kibdc,    kinitc,       krr
    integer          ::         la,       lam,        lb,       lbm,       m0c
    integer          ::      mdncb,     mfbea,     mfbeb,     mfbia,    mfbiam
    integer          ::      mfbib,    mfbibm,       mfc,       mlb, mnc(ip43)
    integer          ::         mt,ncbd(ip41),ncin(ip41),buff(16),lag,lbg,mfbea1,mfbeb1
    integer          ::  id1b,id2b,jd1b,jd2b,kd1b,kd2b,npnb,npcb,nb,nbi,nbj,nbk,nidlb,njdlb,nijdlb
    double precision ::  epsmsh,   exs1,   exs2,x(ip21),y(ip21)
    double precision :: z(ip21),xb,xbi,xbj,xbk,yb,ybi,ybj,ybk,zb,zbi,zbj,zbk,BUFF2(12)
!
!-----------------------------------------------------------------------
!
    character(len=2 ) :: tvi,tvj,tvk
    character(len=7 ) :: eqt
    character(len=6 ) :: typa,typb
!
!
    mfbea1=bcg_to_bcl(mfbea)
    mfbeb1=bcg_to_bcl(mfbeb)
    lag=bcg_to_proc(mfbea)
    lbg=bcg_to_proc(mfbeb)

    if(rank==lbg) then

    mfbib=nfei(mfbeb1)
!
    lb=ndlb(mfbib)
!
    do img=1,lgx
!
       mfbibm=mfbib+(img-1)*mtb
!
       mlb=mpb(mfbibm)+1
       mdncb=ncin(mlb)-ncbd(mlb)
!
       if(kibdc.eq.1) then
!
          typb=indfl(mfbib)
          ib1=iminb(mfbibm)
          ib2=imaxb(mfbibm)
          jb1=jminb(mfbibm)
          jb2=jmaxb(mfbibm)
          kb1=kminb(mfbibm)
          kb2=kmaxb(mfbibm)
!
          lbm=lb+(img-1)*lz
!
          if(krr.eq.1) then
             call extmhg(lbm,x,y,z,exs1,exs2)
          endif
!

      nidlb = id2(lb)-id1(lb)+1
      njdlb = jd2(lb)-jd1(lb)+1
      nijdlb = nidlb*njdlb

      nb = npn(lb)+1+(ib1-id1(lb)) &
                    +(jb1-jd1(lb))*nidlb &
                    +(kb1-kd1(lb))*nijdlb
!
      nbi = nb+1
      nbj = nb+nidlb
      nbk = nb+nijdlb

          xb    =x(nb)
          xbi   =x(nbi)
          xbj   =x(nbj)
          xbk   =x(nbk)
          yb    =y(nb)
          ybi   =y(nbi)
          ybj   =y(nbj)
          ybk   =y(nbk)
          zb    =z(nb)
          zbi   =z(nbi)
          zbj   =z(nbj)
          zbk   =z(nbk)



      call MPI_TRANS(BUFF,[mdncb,lbm,ib1,ib2,jb1,jb2,kb1,kb2,&
          id1(lbm),id2(lbm),jd1(lbm),jd2(lbm),kd1(lbm),kd2(lbm),npn(lbm),npc(lbm)],lbg,lag)
      call MPI_TRANS(BUFF2,[xb,xbi,xbj,xbk,yb,ybi,ybj,ybk,zb,zbi,zbj,zbk],lbg,lag)
      call MPI_TRANS(typb,typb,lbg,lag)
      endif
      enddo
    endif
    if(rank==lag) then

    mfbia=nfei(mfbea1)
!
    la=ndlb(mfbia)
!
    if (cl(mfbia)(1:2).ne.'rh') then
       mtcx=mtcx+1
       mfc =mtcx
       nfbc(mfc)=mfbia
    endif
!
    do img=1,lgx
!
       mfbiam=mfbia+(img-1)*mtb
!
!
       mt=mmb(mfbiam)
!
       mpc(mfbiam)=mdimtcx
       m0c=mpc(mfbiam)
!
       mdimtcx=mdimtcx+mt
!
       ndcc(mfbia)=mfbeb
!
       kinitc=0
!
       if(kibdc.eq.1) then
!
          typa=indfl(mfbia)
          ia1=iminb(mfbiam)
          ia2=imaxb(mfbiam)
          ja1=jminb(mfbiam)
          ja2=jmaxb(mfbiam)
          ka1=kminb(mfbiam)
          ka2=kmaxb(mfbiam)
          eqt=equat
!
          lam=la+(img-1)*lz
!
          if(krr.eq.1) then
             call extmhg(lam,x,y,z,exs1,exs2)
          endif
!
          imgi=img
          imgj=img
          imgk=img
          if (equat(3:5).eq.'2di') imgi = 1
          if (equat(3:5).eq.'2dj') imgj = 1
          if (equat(3:5).eq.'2dk') imgk = 1
          if (equat(3:5).eq.'2xk') imgk = 1
!
          ibam=(iba-ii1(lam))/2**(imgi-1)+ii1(lam)
          jbam=(jba-jj1(lam))/2**(imgj-1)+jj1(lam)
          kbam=(kba-kk1(lam))/2**(imgk-1)+kk1(lam)
!
          call MPI_TRANS(BUFF,BUFF,lbg,lag)
          call MPI_TRANS(BUFF2,BUFF2,lbg,lag)
          call MPI_TRANS(typb,typb,lbg,lag)
          mdncb=BUFF(1)
          lbm  =BUFF(2)
          ib1  =BUFF(3)
          ib2  =BUFF(4)
          jb1  =BUFF(5)
          jb2  =BUFF(6)
          kb1  =BUFF(7)
          kb2  =BUFF(8)
          id1b  =BUFF(9)
          id2b  =BUFF(10)
          jd1b  =BUFF(11)
          jd2b  =BUFF(12)
          kd1b  =BUFF(13)
          kd2b  =BUFF(14)
          npnb  =BUFF(15)
          npcb  =BUFF(16)
          xb    =BUFF2(1)
          xbi   =BUFF2(2)
          xbj   =BUFF2(3)
          xbk   =BUFF2(4)
          yb    =BUFF2(5)
          ybi   =BUFF2(6)
          ybj   =BUFF2(7)
          ybk   =BUFF2(8)
          zb    =BUFF2(9)
          zbi   =BUFF2(10)
          zbj   =BUFF2(11)
          zbk   =BUFF2(12)

          mdnc(mfbiam)=mdncb
          call initcs( &
               x,y,z,krr,epsmsh,kinitc, &
               mfbiam, &
               lam,typa,ia1,ia2,ja1,ja2,ka1,ka2, &
               lbm,typb,ib1,ib2,jb1,jb2,kb1,kb2, &
               xb,xbi,xbj,xbk,yb,ybi,ybj,ybk,zb,zbi,zbj,zbk,&
               ibam,jbam,kbam,tvi,tvj,tvk, &
               eqt, &
               mnc,id1b,id2b,jd1b,jd2b,kd1b,kd2b,npnb,npcb)
!
!     arret en cas de probleme dans les raccords coincidents
          if (kinitc.eq.1) stop 'initc'

       elseif(kibdc.eq.0) then
!           call readfc( &
!                kfc,mnc, &
!                mt,m0c)
       end if
    enddo
    end if
    call barrier
!
    return
  end subroutine inbdc
end module mod_inbdc
