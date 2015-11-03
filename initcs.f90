module mod_initcs
  implicit none
contains
  subroutine initcs( &
       x,y,z,krr,eps,kinitc, &
       mfba, &
       la,typa,ia1,ia2,ja1,ja2,ka1,ka2, &
       lb,typb,ib1,ib2,jb1,jb2,kb1,kb2, &
       iba,jba,kba,tvi,tvj,tvk, &
       equat, &
       mnc)
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
    implicit none
    integer          ::      i,   ia1,   ia2,   ib1,   ib2
    integer          ::    iba,    ii,  imax,   img, inpcb
    integer          ::    ipb,     j,   ja1,   ja2,   jb1
    integer          ::    jb2,   jba,    jj,  jmax,   jpb
    integer          ::      k,   ka1,   ka2,   kb1,   kb2
    integer          ::    kba,kinitc,    kk,  kmax,   kpb
    integer          ::    krr,     l,    la,    lb
    integer          ::      m,   m0c,    mc,  mfba,   mnc
    integer          ::     na,   nai,   naj,   nak,    nb
    integer          ::    nbv, ncilb, ncjlb, ncklb, nidla
    integer          ::  nidlb,nijdla,nijdlb, njdla, njdlb
    integer          ::    nvi,   nvj,   nvk
    double precision :: cnrota,  dist,   eps,snrota,     x
    double precision ::      y,   ynb,     z,   znb
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=6 ) :: typa,typb
    character(len=7 ) :: equat
    character(len=2 ) :: tvi,tvj,tvk
    dimension x(ip21),y(ip21),z(ip21)
    dimension mnc(ip43)
!
!---------------------------------------------------------------------
!
      nidla = id2(la)-id1(la)+1
      njdla = jd2(la)-jd1(la)+1
      nijdla = nidla*njdla
      nidlb = id2(lb)-id1(lb)+1
      njdlb = jd2(lb)-jd1(lb)+1
      nijdlb = nidlb*njdlb
!
      ncilb = 1
      ncjlb = nidlb
      ncklb = nijdlb
!
      if (krr.eq.1) then
!---------------------------------------------------------------------
      na = npn(la)+1+(ia1-id1(la)) &
                    +(ja1-jd1(la))*nidla &
                    +(ka1-kd1(la))*nijdla
!
      nai = na+1
      naj = na+nidla
      nak = na+nijdla
!
      do k=kb1,kb2
       do j=jb1,jb2
        do i=ib1,ib2
         iba=i
         jba=j
         kba=k
         nb = npn(lb)+1+(i-id1(lb)) &
                       +(j-jd1(lb))*nidlb &
                       +(k-kd1(lb))*nijdlb
!
!print*,x(na),x(nb),i
!print*,y(na),y(nb),j
!print*,z(na),z(nb),k
!print*,la   ,lb     
         dist=sqrt( (x(nb)-x(na))**2+(y(nb)-y(na))**2+(z(nb)-z(na))**2 )
         if(dist.lt.eps) go to 51
!         if(dist.lt.eps) exit
        enddo
       enddo
      enddo 
      kinitc=1
   51 continue
!
!-------
!
      do k=kb1,kb2
       do j=jb1,jb2
        do i=ib1,ib2
         ii=i
         jj=j
         kk=k
         nb = npn(lb)+1+(i-id1(lb)) &
                       +(j-jd1(lb))*nidlb &
                       +(k-kd1(lb))*nijdlb
         dist=sqrt( (x(nb)-x(nai))**2+(y(nb)-y(nai))**2+(z(nb)-z(nai))**2 )
!         if(dist.lt.eps) exit
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
      do k=kb1,kb2
       do j=jb1,jb2
        do i=ib1,ib2
         ii=i
         jj=j
         kk=k
         nb = npn(lb)+1+(i-id1(lb)) &
                       +(j-jd1(lb))*nidlb &
                       +(k-kd1(lb))*nijdlb
         dist=sqrt( (x(nb)-x(naj))**2+(y(nb)-y(naj))**2+(z(nb)-z(naj))**2 )
         if(dist.lt.eps) go to 71
!         if(dist.lt.eps) exit 
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
      do k=kb1,kb2
       do j=jb1,jb2
        do i=ib1,ib2
         ii=i
         jj=j
         kk=k
         nb = npn(lb)+1+(i-id1(lb)) &
                       +(j-jd1(lb))*nidlb &
                       +(k-kd1(lb))*nijdlb
         dist=sqrt( (x(nb)-x(nak))**2+(y(nb)-y(nak))**2+(z(nb)-z(nak))**2 )
!         if(dist.lt.eps) exit 
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
      if(equat(3:5).ne.'2dk' .and. equat(3:5).ne.'2xk' ) then
      if(kk-kba.eq.+1) tvk = '+k'
      if(kk-kba.eq.-1) tvk = '-k'
      endif
!
!---------------------------------------------------------------------
      endif
!
      if (kimp.ge.2) then
!
      l=mod(la,lz)
      img=(la-l)/lz+1
!
      form='(/3x,''numero de grille  : '',i3,/, &
             3x,''indices du point du domaine b coincident '', &
                ''avec le point ia1,ja1,ka1 du domaine a'',/, &
             19x,''iba = '',i5,''    jba = '',i5,''    kba = '',i5,/ &
             3x,''sens de variation des indices dans le domaine b '', &
                ''en fonction de la variation des indices dans '', &
                ''le domaine a'',/, &
              19x,''tvi ='',4x,a,4x,''tvj ='',4x,a,4x,''tvk ='',4x,a,/ &
              19x,''rotation de frt b (nb de pas)  ='',i3)'
      write(imp,form) img,iba,jba,kba,tvi,tvj,tvk,0
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
         nb = npc(lb)+1+(iba-id1(lb)) &
                       +(jba-jd1(lb))*nidlb &
                       +(kba-kd1(lb))*nijdlb &
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

