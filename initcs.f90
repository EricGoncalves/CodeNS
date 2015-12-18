module mod_initcs
  implicit none
contains
  subroutine initcs( &
       x,y,z,krr,eps,kinitc, &
       mfba, &
       la,typa,ia1,ia2,ja1,ja2,ka1,ka2, &
       lb,typb,ib1,ib2,jb1,jb2,kb1,kb2, &
       xb,xbi,xbj,xbk,yb,ybi,ybj,ybk,zb,zbi,zbj,zbk,&
       iba,jba,kba,tvi,tvj,tvk, &
       equat, &
       mnc,id1b,id2b,jd1b,jd2b,kd1b,kd2b,npnb,npcb)
!
!***********************************************************************
!
!     ACT
!_A    Determination des points coincidents a une frontiere
!_A    d'un domaine structure.
!
!     INP
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    krr        : arg int              ; cle info sur front coinc
!_I    eps        : arg real             ; distance max entre 2 pts confondus
!_I    mfba       : arg int              ; numero de frontiere
!_I    la         : arg int              ; numero de domaine
!_I    typa       : arg char             ; type de plan de la frontiere a
!_I    ia1        : arg int              ; indice min en i de frontiere
!_I    ia2        : arg int              ; indice max en i de frontiere
!_I    ja1        : arg int              ; indice min en j de frontiere
!_I    ja2        : arg int              ; indice max en j de frontiere
!_I    ka1        : arg int              ; indice min en k de frontiere
!_I    ka2        : arg int              ; indice max en k de frontiere
!_I    lb         : arg int              ; numero du domaine coincident
!_I    typb       : arg char             ; type de plan de la frontiere b
!_I    ib1        : arg int              ; indice min en i de front coinc
!_I    ib2        : arg int              ; indice max en i de front coinc
!_I    jb1        : arg int              ; indice min en j de front coinc
!_I    jb2        : arg int              ; indice max en j de front coinc
!_I    kb1        : arg int              ; indice min en k de front coinc
!_I    kb2        : arg int              ; indice max en k de front coinc
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
!_I    equat      : arg char             ; type d'equations modelisant l'ecoulement
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
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
!_I    mpc        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux front coinc
!
!     OUT
!_O    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_O                                        cellule coincidente
!
!     I/O
!_/    kinitc     : arg int              ; cle controle initialisation des
!_/                                        raccords coincidents
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use definition
    use sortiefichier
    use mod_mpi
    implicit none
    integer          ::      i,   ia1,   ia2,   ib1,   ib2
    integer          ::    iba,    ii,  imax,   img, inpcb
    integer          ::      j,   ja1,   ja2,   jb1
    integer          ::    jb2,   jba,    jj,  jmax
    integer          ::      k,   ka1,   ka2,   kb1,   kb2
    integer          ::    kba,kinitc,    kk,  kmax
    integer          ::    krr,     l,    la,    lb
    integer          ::      m,   m0c,    mc,  mfba,   mnc(ip43)
    integer          ::     na,    nb
    integer          ::    nbv, ncilb, ncjlb, ncklb, nidla
    integer          ::  nidlb,nijdla,nijdlb, njdla, njdlb
    integer          ::    nvi,   nvj,   nvk
    integer          ::  id1b,id2b,jd1b,jd2b,kd1b,kd2b,npnb,npcb
    double precision ::  dist,   eps,x(ip21),y(ip21),z(ip21)
    double precision ::  xb,xbi,xbj,xbk,yb,ybi,ybj,ybk,zb,zbi,zbj,zbk
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=6 ) :: typa,typb
    character(len=7 ) :: equat
    character(len=2 ) :: tvi,tvj,tvk
!
!---------------------------------------------------------------------
!
      nidla = id2(la)-id1(la)+1
      njdla = jd2(la)-jd1(la)+1
      nijdla = nidla*njdla
!
!
!
      nidlb = id2b-id1b+1
      njdlb = jd2b-jd1b+1
      nijdlb = nidlb*njdlb

      ncilb = 1
      ncjlb = nidlb
      ncklb = nijdlb
!
      if (krr.eq.1) then
!---------------------------------------------------------------------
! look for the first point
      do k=ka1,ka2
       do j=ja1,ja2
        do i=ia1,ia2
         iba=i
         jba=j
         kba=k
         na = npn(la)+1+(i-id1(la)) &
                       +(j-jd1(la))*nidla &
                       +(k-kd1(la))*nijdla
!
         dist=sqrt( (x(na)-xb)**2+(y(na)-yb)**2+(z(na)-zb)**2 )
         if(dist.lt.eps) go to 51
        enddo
       enddo
      enddo 
      kinitc=1
   51 continue
!
!-------
! look for directions
      do k=ka1,ka2
       do j=ja1,ja2
        do i=ia1,ia2
         ii=i
         jj=j
         kk=k
         na = npn(la)+1+(i-id1(la)) &
                       +(j-jd1(la))*nidla &
                       +(k-kd1(la))*nijdla
                       
         dist=sqrt( (x(na)-xbi)**2+(y(na)-ybi)**2+(z(na)-zbi)**2 )
         if(dist.lt.eps) go to 61
        enddo
       enddo
      enddo 
   61 continue
!
      tvi = 'fa'
!
      if(equat(3:5).ne.'2di') then
      if(ii-iba.eq.+1) tvi = '+i'
      if(ii-iba.eq.-1) tvi = '-i'
      endif
!
      if(equat(3:5).ne.'2dj') then
      if(jj-jba.eq.+1) tvi = '+j'
      if(jj-jba.eq.-1) tvi = '-j'
      endif
!
!VC   if(equat(3:5).ne.'2dk') then
      if(equat(3:5).ne.'2dk' .and. equat(3:5).ne.'2xk' ) then
      if(kk-kba.eq.+1) tvi = '+k'
      if(kk-kba.eq.-1) tvi = '-k'
      endif
!
!-------
!
      do k=ka1,ka2
       do j=ja1,ja2
        do i=ia1,ia2
         ii=i
         jj=j
         kk=k
         na = npn(la)+1+(i-id1(la)) &
                       +(j-jd1(la))*nidla &
                       +(k-kd1(la))*nijdla
         dist=sqrt( (x(na)-xbj)**2+(y(na)-ybj)**2+(z(na)-zbj)**2 )
         if(dist.lt.eps) go to 71
        enddo
       enddo
      enddo 
   71 continue
!
      tvj = 'fa'
!
      if(equat(3:5).ne.'2di') then
      if(ii-iba.eq.+1) tvj = '+i'
      if(ii-iba.eq.-1) tvj = '-i'
      endif
!
      if(equat(3:5).ne.'2dj') then
      if(jj-jba.eq.+1) tvj = '+j'
      if(jj-jba.eq.-1) tvj = '-j'
      endif
!
!VC   if(equat(3:5).ne.'2dk') then
      if(equat(3:5).ne.'2dk' .and. equat(3:5).ne.'2xk' ) then
      if(kk-kba.eq.+1) tvj = '+k'
      if(kk-kba.eq.-1) tvj = '-k'
      endif
! 
!-------
!
      do k=ka1,ka2
       do j=ja1,ja2
        do i=ia1,ia2
         ii=i
         jj=j
         kk=k
         na = npn(la)+1+(i-id1(la)) &
                       +(j-jd1(la))*nidla &
                       +(k-kd1(la))*nijdla
         dist=sqrt( (x(na)-xbk)**2+(y(na)-ybk)**2+(z(na)-zbk)**2 )
         if(dist.lt.eps) go to 81
        enddo
       enddo
      enddo
   81 continue
!
      tvk = 'fa'
!
      if(equat(3:5).ne.'2di') then
      if(ii-iba.eq.+1) tvk = '+i'
      if(ii-iba.eq.-1) tvk = '-i'
      endif
!
      if(equat(3:5).ne.'2dj') then
      if(jj-jba.eq.+1) tvk = '+j'
      if(jj-jba.eq.-1) tvk = '-j'
      endif
!
!VC   if(equat(3:5).ne.'2dk') then
!      if(equat(3:5).ne.'2dk' .and. equat(3:5).ne.'2xk' ) then
      if(kk-kba.eq.+1) tvk = '+k'
      if(kk-kba.eq.-1) tvk = '-k'
!      endif
!
!---------------------------------------------------------------------

      if(tvi.eq.'fa') iba = ib1
      if(tvi.eq.'+i') iba = ib1
      if(tvi.eq.'-i') iba = ib2
      if(tvi.eq.'+j') iba = jb1
      if(tvi.eq.'-j') iba = jb2
      if(tvi.eq.'+k') iba = kb1
      if(tvi.eq.'-k') iba = kb2
      
      if(tvj.eq.'fa') jba = jb1
      if(tvj.eq.'+i') jba = ib1
      if(tvj.eq.'-i') jba = ib2
      if(tvj.eq.'+j') jba = jb1
      if(tvj.eq.'-j') jba = jb2
      if(tvj.eq.'+k') jba = kb1
      if(tvj.eq.'-k') jba = kb2
      
      if(tvk.eq.'fa') kba = kb1
      if(tvk.eq.'+i') kba = ib1
      if(tvk.eq.'-i') kba = ib2
      if(tvk.eq.'+j') kba = jb1
      if(tvk.eq.'-j') kba = jb2
      if(tvk.eq.'+k') kba = kb1
      if(tvk.eq.'-k') kba = kb2


      endif
!
      if (kimp.ge.2) then
!
       l=mod(la,lz)+1
       img=(la-l+1)/lz+1
!
       form='(/3x,''numero de grille  : '',i3,/,' &
            //'3x,''indices du point du domaine b coincident '',' &
            //'''avec le point ia1,ja1,ka1 du domaine a'',/,' &
            //'19x,''iba = '',i5,''    jba = '',i5,''    kba = '',i5,/' &
            //'3x,''sens de variation des indices dans le domaine b '',' &
            //'''en fonction de la variation des indices dans '',' &
            //'''le domaine a'',/,' &
            //'19x,''tvi ='',4x,a,4x,''tvj ='',4x,a,4x,''tvk ='',4x,a,/' &
            //'19x,''rotation de frt b (nb de pas)  ='',i3)'
       write(imp,form) bl_to_bg(l),iba,jba,kba,tvi,tvj,tvk,0
      endif
!
!-----------------------------------------------------------------------------
!
!---- variation des indices du domaine b sur la frontiere commune
!     en fonction de la variation des indices du domaine a
!
      if(tvi.eq.'fa') nvi = 0
      if(tvi.eq.'+i') nvi = ncilb
      if(tvi.eq.'-i') nvi =-ncilb
      if(tvi.eq.'+j') nvi = ncjlb
      if(tvi.eq.'-j') nvi =-ncjlb
      if(tvi.eq.'+k') nvi = ncklb
      if(tvi.eq.'-k') nvi =-ncklb
!
      if(tvj.eq.'fa') nvj = 0
      if(tvj.eq.'+i') nvj = ncilb
      if(tvj.eq.'-i') nvj =-ncilb
      if(tvj.eq.'+j') nvj = ncjlb
      if(tvj.eq.'-j') nvj =-ncjlb
      if(tvj.eq.'+k') nvj = ncklb
      if(tvj.eq.'-k') nvj =-ncklb
!
      if(tvk.eq.'fa') nvk = 0
      if(tvk.eq.'+i') nvk = ncilb
      if(tvk.eq.'-i') nvk =-ncilb
      if(tvk.eq.'+j') nvk = ncjlb
      if(tvk.eq.'-j') nvk =-ncjlb
      if(tvk.eq.'+k') nvk = ncklb
      if(tvk.eq.'-k') nvk =-ncklb
!
!-----------------------------------------------------------------------------
!
!---- type de frontiere du domaine b :
!     calcul de l'increment vers le point fictif
!
!---- faces
      if(typb.eq.'i1    ')   inpcb =  0
      if(typb.eq.'i2    ')   inpcb = -ncilb
      if(typb.eq.'j1    ')   inpcb =  0
      if(typb.eq.'j2    ')   inpcb = -ncjlb
      if(typb.eq.'k1    ')   inpcb  = 0
      if(typb.eq.'k2    ')   inpcb = -ncklb
!
!--   adaptation cell centered
      imax=ia2-1
      jmax=ja2-1
      kmax=ka2-1
!
      if(typa.eq.'i1    ')   imax=ia2
      if(typa.eq.'i2    ')   imax=ia2
      if(typa.eq.'j1    ')   jmax=ja2
      if(typa.eq.'j2    ')   jmax=ja2
      if(typa.eq.'k1    ')   kmax=ka2
      if(typa.eq.'k2    ')   kmax=ka2
!
      m0c=mpc(mfba)
      m=0
      do k=ka1,kmax
       do j=ja1,jmax
        do i=ia1,imax
         m=m+1
         mc=m0c+m
!
         nb = npcb+1+(iba-id1b) &
                       +(jba-jd1b)*nidlb &
                       +(kba-kd1b)*nijdlb &
                       +(i-ia1)*nvi &
                       +(j-ja1)*nvj &
                       +(k-ka1)*nvk 
!
!--   modification du stockage du point du domaine b
!     prenant en compte le stockage en bas a gauche
!     si le pas est negatif on decroit suivant l'indice
!om   if(nvi.lt.0) nb=nb - pas suivant = nb-abs(nvi) = nb+nvi
!om   if(nvj.lt.0) nb=nb - pas suivant = nb-abs(nvj) = nb+nvj
!om   if(nvk.lt.0) nb=nb - pas suivant = nb-abs(nvk) = nb+nvk
         if(nvi.lt.0) nb=nb+nvi
         if(nvj.lt.0) nb=nb+nvj
         if(nvk.lt.0) nb=nb+nvk
!
!---- tableau de coincidence des indices
!
         nbv=nb+inpcb
         mnc(mc)=nbv
        enddo
       enddo
      enddo
!
      return
  end subroutine initcs
end module mod_initcs

