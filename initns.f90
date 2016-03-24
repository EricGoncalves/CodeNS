module mod_initns
  implicit none
contains
  subroutine initns( &
       mfb,l,indfb, &
       imin,imax,jmin,jmax,kmin,kmax, &
       eqt,x,y,z, &
       ncbd,nxn,nyn,nzn, &
       tnix,tniy,tniz, &
       tnjx,tnjy,tnjz, &
       tnkx,tnky,tnkz)
!
!***********************************************************************
!
!     ACT
!_A    Determination des normales interieures aux facettes
!_A    d'une frontiere d'un domaine structure.
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
!_I    mfb        : arg int              ; numero de frontiere
!_I    l          : arg int              ; numero de domaine
!_I    indfb      : arg char             ; type de plan de la frontiere
!_I    imin       : arg int              ; indice min en i
!_I    imax       : arg int              ; indice max en i
!_I    jmin       : arg int              ; indice min en j
!_I    jmax       : arg int              ; indice max en j
!_I    kmin       : arg int              ; indice min en k
!_I    kmax       : arg int              ; indice max en k
!_I    eqt        : arg char             ; type d'equations modelisant l'ecoulement
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
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
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!
!     OUT
!_O    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_O                                        normal a une facette frontiere
!_O    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_O                                        normal a une facette frontiere
!_O    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_O                                        normal a une facette frontiere
!
!     LOC
!_L    tnix       : arg real(ip00      ) ; composante en x du vect normal a une
!_L                                        facette i et de norme egale a la
!_L                                        surface de celle-ci
!_L    tniy       : arg real(ip00      ) ; composante en y du vect normal a une
!_L                                        facette i et de norme egale a la
!_L                                        surface de celle-ci
!_L    tniz       : arg real(ip00      ) ; composante en z du vect normal a une
!_L                                        facette i et de norme egale a la
!_L                                        surface de celle-ci
!_L    tnjx       : arg real(ip00      ) ; composante en x du vect normal a une
!_L                                        facette j et de norme egale a la
!_L                                        surface de celle-ci
!_L    tnjy       : arg real(ip00      ) ; composante en y du vect normal a une
!_L                                        facette j et de norme egale a la
!_L                                        surface de celle-ci
!_L    tnjz       : arg real(ip00      ) ; composante en z du vect normal a une
!_L                                        facette j et de norme egale a la
!_L                                        surface de celle-ci
!_L    tnkx       : arg real(ip00      ) ; composante en x du vect normal a une
!_L                                        facette k et de norme egale a la
!_L                                        surface de celle-ci
!_L    tnky       : arg real(ip00      ) ; composante en y du vect normal a une
!_L                                        facette k et de norme egale a la
!_L                                        surface de celle-ci
!_L    tnkz       : arg real(ip00      ) ; composante en z du vect normal a une
!_L                                        facette k et de norme egale a la
!_L                                        surface de celle-ci
!
!     COM
!_C    Le tableau ncbd doit avoir ete rempli au prealable.
!_C
!_C    La frontiere est une partie structuree du plan de maillage
!_C    d'indice indfb avec indfb (variable caractere) qui vaut
!_C    'i1', 'j1', 'k1', 'i2', 'j2' ou 'k2'.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use mod_norm
    implicit none
    integer          ::         dm,         i,        i1,        i2,        id
    integer          ::       imax,      imin,        is,         j,        j1
    integer          ::         j2,        jd,      jmax,      jmin,        js
    integer          ::          k,        k1,        k2,        kd,      kmax
    integer          ::       kmin,        ks,         l,       m0b,       m0n
    integer          ::         m1,     m1max,     m1min,        m2,     m2max
    integer          ::      m2min,        mb,       mfb,        mn,        n0
    integer          :: ncbd(ip41),       nci,       ncj,       nck,       nid
    integer          ::       nijd,      nip0,       njd
    double precision ::        eps,        nn, nxn(ip42), nyn(ip42), nzn(ip42)
    double precision :: tnix(ip00),tniy(ip00),tniz(ip00),tnjx(ip00),tnjy(ip00)
    double precision :: tnjz(ip00),tnkx(ip00),tnky(ip00),tnkz(ip00),   x(ip21)
    double precision ::    y(ip21),   z(ip21)
!
!-----------------------------------------------------------------------
!
    character(len=2 ) :: indfb
    character(len=7 ) :: eqt
!
    eps=1.d-10
!
    n0=npc(l)
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
    m0n=mpn(mfb)
    m0b=mpb(mfb)
!
    nip0=ip00
!
!     calcul des normales interieures
!
    if( (indfb.eq.'i1').or.(indfb.eq.'i2') ) then
       m1min=1
       m1max=jmax-jmin+1-1
       m2min=1
       m2max=kmax-kmin+1-1
       dm=m1max-m1min+1
!
!     test sur l'extremite pour correspondance avec les tableaux de normales
!     (test sur imin=imax)
       if(imin.eq.i1) is= 1
       if(imin.eq.i2) is=-1
       if(imin.eq.i1) id= nci
       if(imin.eq.i2) id= 0
!
       call norm( &
            l,x,y,z, &
            eqt,nip0, &
            imin,imax,jmin,jmax,kmin,kmax, &
            tnix,tniy,tniz, &
            tnjx,tnjy,tnjz, &
            tnkx,tnky,tnkz)
!
       mn=m0n
       mb=m0b
       do k=kmin,kmax-1
          do j=jmin,jmax-1
             do i=imin,imax
                mn=mn+1
                mb=mb+1
                nxn(mn)=is*tnix(ncbd(mb)+id-n0)
                nyn(mn)=is*tniy(ncbd(mb)+id-n0)
                nzn(mn)=is*tniz(ncbd(mb)+id-n0)
             enddo
          enddo
       enddo
!
    elseif ( (indfb.eq.'j1').or.(indfb.eq.'j2') ) then
       m1min=1
       m1max=imax-imin+1-1
       m2min=1
       m2max=kmax-kmin+1-1
       dm   =m1max-m1min+1
!
!     test sur l'extremite pour correspondance avec les tableaux de normales
!     (test sur jmin=jmax)
       if(jmin.eq.j1) js= 1
       if(jmin.eq.j2) js=-1
       if(jmin.eq.j1) jd= ncj
       if(jmin.eq.j2) jd= 0
!
       call norm( &
            l,x,y,z, &
            eqt,nip0, &
            imin,imax,jmin,jmax,kmin,kmax, &
            tnix,tniy,tniz, &
            tnjx,tnjy,tnjz, &
            tnkx,tnky,tnkz)
!
       mn=m0n
       mb=m0b
       do k=kmin,kmax-1
          do j=jmin,jmax
             do i=imin,imax-1
                mn=mn+1
                mb=mb+1
                nxn(mn)=js*tnjx(ncbd(mb)+jd-n0)
                nyn(mn)=js*tnjy(ncbd(mb)+jd-n0)
                nzn(mn)=js*tnjz(ncbd(mb)+jd-n0)
             enddo
          enddo
       enddo
!
    elseif ( (indfb.eq.'k1').or.(indfb.eq.'k2') ) then
       m1min=1
       m1max=imax-imin+1-1
       m2min=1
       m2max=jmax-jmin+1-1
       dm=m1max-m1min+1
!
!     test sur l'extremite pour correspondance avec les tableaux de normales
!     (test sur kmin=kmax)
       if(kmin.eq.k1) ks= 1
       if(kmin.eq.k2) ks=-1
       if(kmin.eq.k1) kd= nck
       if(kmin.eq.k2) kd= 0
!
       call norm( &
            l,x,y,z, &
            eqt,nip0, &
            imin,imax,jmin,jmax,kmin,kmax, &
            tnix,tniy,tniz, &
            tnjx,tnjy,tnjz, &
            tnkx,tnky,tnkz)
!
       mn=m0n
       mb=m0b
       do k=kmin,kmax
          do j=jmin,jmax-1
             do i=imin,imax-1
                mn=mn+1
                mb=mb+1
                nxn(mn)=ks*tnkx(ncbd(mb)+kd-n0)
                nyn(mn)=ks*tnky(ncbd(mb)+kd-n0)
                nzn(mn)=ks*tnkz(ncbd(mb)+kd-n0)
             enddo
          enddo
       enddo
!
    end if
!
!     m1 : m1min a m1max ,  m2 : m2min a m2max
!
    do m2=m2min,m2max
       do m1=m1min,m1max
          mn=m0n+m1+(m2-1)*dm
!
          nn=sqrt(nxn(mn)*nxn(mn)+nyn(mn)*nyn(mn)+nzn(mn)*nzn(mn))
          nn=nn+(sign(1.D0,-abs(nn))+1.D0)*eps
          nxn(mn)=nxn(mn)/nn
          nyn(mn)=nyn(mn)/nn
          nzn(mn)=nzn(mn)/nn
       enddo
    enddo
!
    return
  end subroutine initns
end module mod_initns
