module mod_clglis
implicit none
contains
      subroutine clglis( &
                 mfb,l,indf, &
                 ncin,ncbd, &
                 x,y,z, &
                 sn,lgsnlt,vol, &
                 v,pression,ztemp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement de la condition de glissement (qn=0.)
!_A    Extrapolation d'ordre superieur:
!_A      varext = varint + dvarint . dx  .
!_A    Choix sur les relation utilisee:
!_A      1) ro = const ,  pt=const , qt//qtint , p=pext  ,
!_A      2) h  = const , ent=const , qt//qtint , p=pext  .
!_A    Normales interieures.
!
!     INP
!_I    mfb        : arg int              ; numero de frontiere
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    sn         : arg real(lgsnlt,
!_I                          nind,ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    lgsnlt     : arg int              ; nombre de noeuds du dom (dont fic.)
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    indfl      : com char(mtb       ) ; type de plan de la frontiere
!_I    gam        : com real             ; rapport des chaleurs specifiques
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    reelmn     : com real             ; nombre reel petit
!
!     COM
!_C    Notation s      : valeurs issues du schema.
!_C    Notation   ou 1 : valeurs a l' issue de l' application de la condition
!_C                      a la limite.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use proprieteflu
      use constantes
      use boundary
use mod_idirch
implicit none
integer :: inc
integer :: mfb
integer :: l
integer :: ncin
integer :: ncbd
double precision :: x
double precision :: y
double precision :: z
double precision :: sn
integer :: lgsnlt
double precision :: vol
double precision :: v
double precision :: pression
double precision :: ztemp
double precision :: cson
integer :: id
integer :: jd
integer :: kd
double precision :: a2s
double precision :: coef
double precision :: dp
double precision :: ds
double precision :: dsn
double precision :: dx
double precision :: dy
double precision :: dz
double precision :: ec
double precision :: ecs
double precision :: ent
double precision :: ents
double precision :: gm
double precision :: gm1
double precision :: gm1rg
double precision :: grgm1
double precision :: h
double precision :: hs
integer :: incf
integer :: incm
integer :: isens
integer :: m
integer :: mb
integer :: mn
integer :: mnode
integer :: mt
integer :: n0c
integer :: n0n
integer :: nc
integer :: nc12f
integer :: nc1f
integer :: nc2f
integer :: nci
integer :: ncij
integer :: ncijk
integer :: ncj
integer :: ncjk
integer :: nck
integer :: ncki
integer :: nd
integer :: nid
integer :: nijd
integer :: nindf
integer :: nindfa
integer :: njd
integer :: nlc
integer :: nln
integer :: nlnf
integer :: nn
double precision :: p
double precision :: pext
double precision :: ps
double precision :: pts
double precision :: q2
double precision :: q2s
double precision :: qns
double precision :: qt2s
double precision :: qtx
double precision :: qtxs
double precision :: qty
double precision :: qtys
double precision :: qtz
double precision :: qtzs
double precision :: r2vol
double precision :: rgm
double precision :: rgm1
double precision :: ro
double precision :: roe
double precision :: roec
double precision :: roes
double precision :: ros
double precision :: rou
double precision :: rous
double precision :: rov
double precision :: rovs
double precision :: row
double precision :: rows
double precision :: rros
double precision :: rvol
double precision :: sx
double precision :: sxn
double precision :: sy
double precision :: syn
double precision :: sz
double precision :: szn
double precision :: tworgm1
double precision :: us
double precision :: vs
double precision :: ws
double precision :: xcp
double precision :: xcs
double precision :: ycp
double precision :: ycs
double precision :: zcp
double precision :: zcs
!
!-----------------------------------------------------------------------
!
      character(len=2 ) :: indf
      dimension x(ip21),y(ip21),z(ip21)
      dimension vol(ip11),pression(ip11),ztemp(ip11),cson(ip11)
      dimension v(ip11,ip60)
      dimension sn(lgsnlt,nind,ndir)
      dimension ncbd(ip41),ncin(ip41)
!
!.....calcul des indices
!
      inc(id,jd,kd)=id+jd*nid+kd*nijd
!
!.....info frontiere: nb. pts., n0 dom.
!
      mt     = mmb(mfb)
!
!.....indice de face
!
      nindf  = idirch(indf)
      nindfa = abs(nindf)
      isens  = sign(1,nindf)
!
!.....calcul des indices
!
      n0n    = npn(l)
      n0c    = npc(l)
      nid    = id2(l)-id1(l)+1
      njd    = jd2(l)-jd1(l)+1
      nijd   = nid*njd
!
      nci    = inc(1,0,0)
      ncj    = inc(0,1,0)
      nck    = inc(0,0,1)
      ncij   = inc(1,1,0)
      ncjk   = inc(0,1,1)
      ncki   = inc(1,0,1)
      ncijk  = inc(1,1,1)
!
      if      (nindfa.eq.1) then
       nc1f  = inc(0,1,0)
       nc2f  = inc(0,0,1)
       nc12f = inc(0,1,1)
      else if (nindfa.eq.2) then
       nc1f  = inc(0,0,1)
       nc2f  = inc(1,0,0)
       nc12f = inc(1,0,1)
      else if (nindfa.eq.3) then
       nc1f  = inc(1,0,0)
       nc2f  = inc(0,1,0)
       nc12f = inc(1,1,0)
      end if
!
      gm     = gam
      gm1    = gm-1.0
      rgm    = 1.0/gm
      rgm1   = 1.0/gm1
      grgm1  = gm*rgm1
      tworgm1= 2.0*rgm1
      gm1rg  = 1.0/grgm1
!
      do m=1,mt
!
!.......initialisation
!
        dp   = 0.0
!
!........indice de la paroi (p)
!........ m  = indice relatif de frontiere
!........ mb = indice global de frontiere de base
!........ nd = indice global dans domaine de la cellule
!
        mb   = mpb(mfb) + m
        nd   = ncbd(mb)
        incm = ncin(mb)-nd
        incf = incm*((isens+1)/2)
        mnode= nd-n0c + incf
!
        nlc  = nd
        nln  = nd-n0c+n0n
        nlnf = nd-n0c+n0n+incf
!
!.......indice de la cellule (s) et du noeud en bas a gauche
!
        nc   = nlc+incm
        nn   = nln+incm
!
!.......indice du noeud dans le tableau sn
!
        mn   = nn - n0n
!
!........normale frontiere interieure
!
        sx   = isens*sn(mnode,nindfa,1)
        sy   = isens*sn(mnode,nindfa,2)
        sz   = isens*sn(mnode,nindfa,3)
        ds   = sqrt(sx*sx + sy*sy + sz*sz)
        ds   = max(ds,reelmn)
        sxn  = sx/ds
        syn  = sy/ds
        szn  = sz/ds
!
!.......valeur du champ en s
!
        ros  = v(nc,1)
        rous = v(nc,2)
        rovs = v(nc,3)
        rows = v(nc,4)
        roes = v(nc,5)
        rros = 1.0/ros
        us   = rous*rros
        vs   = rovs*rros
        ws   = rows*rros
        q2s  = us*us + vs*vs + ws*ws
        ecs  = 0.5*q2s
        qns  = us*sxn + vs*syn + ws*szn
        qtxs = us - qns*sxn
        qtys = vs - qns*syn
        qtzs = ws - qns*szn
        qt2s = qtxs*qtxs + qtys*qtys + qtzs*qtzs
        ps   = gm1*(roes - pinfl - ros*ecs)
        a2s  = gm*ps/ros
        pts  = ps*(1.0 + gm1*ecs/a2s)**grgm1
        hs   = (roes + ps)*rros
        ents = ps*rros**gm
!
!.......calcul distance paroi p - centre maille s
!
        xcs  = ( x(nn    )+x(nn+nci )+x(nn+ncj )+x(nn+ncij ) &
               + x(nn+nck)+x(nn+ncki)+x(nn+ncjk)+x(nn+ncijk) ) &
              *0.125
        ycs  = ( y(nn    )+y(nn+nci )+y(nn+ncj )+y(nn+ncij ) &
               + y(nn+nck)+y(nn+ncki)+y(nn+ncjk)+y(nn+ncijk) ) &
              *0.125
        zcs  = ( z(nn    )+z(nn+nci )+z(nn+ncj )+z(nn+ncij ) &
               + z(nn+nck)+z(nn+ncki)+z(nn+ncjk)+z(nn+ncijk) ) &
              *0.125
!
        xcp  = (x(nlnf)+x(nlnf+nc12f)+x(nlnf+nc1f)+x(nlnf+nc2f))*0.25
        ycp  = (y(nlnf)+y(nlnf+nc12f)+y(nlnf+nc1f)+y(nlnf+nc2f))*0.25
        zcp  = (z(nlnf)+z(nlnf+nc12f)+z(nlnf+nc1f)+z(nlnf+nc2f))*0.25
!
        dx   = xcs - xcp
        dy   = ycs - ycp
        dz   = zcs - zcp
!
!.......calcul des gradients
!
!.......    b u
!.......     \|
!.......   l--s--r
!.......      |\
!.......      d f
!
!.......flux l
!
        roec = 0.5*(v(nc-nci,2)*v(nc-nci,2)+v(nc-nci,3)*v(nc-nci,3) &
                   +v(nc-nci,4)*v(nc-nci,4))/v(nc-nci,1)
        p    = gm1*(v(nc-nci,5) -pinfl - roec)
        dsn  =-(sn(mn,1,1)*dx+sn(mn,1,2)*dy+sn(mn,1,3)*dz)
        dp   = dp   + (ps   + p        )*dsn
!
!.......flux r
!
        roec = 0.5*(v(nc+nci,2)*v(nc+nci,2)+v(nc+nci,3)*v(nc+nci,3) &
                   +v(nc+nci,4)*v(nc+nci,4))/v(nc+nci,1)
        p    = gm1*(v(nc+nci,5) -pinfl - roec)
        dsn  = (sn(mn+nci,1,1)*dx+sn(mn+nci,1,2)*dy+sn(mn+nci,1,3)*dz)
        dp   = dp   + (ps   + p        )*dsn
!
!.......flux u
!
        roec = 0.5*(v(nc-ncj,2)*v(nc-ncj,2)+v(nc-ncj,3)*v(nc-ncj,3) &
                   +v(nc-ncj,4)*v(nc-ncj,4))/v(nc-ncj,1)
        p    = gm1*(v(nc-ncj,5) -pinfl - roec)
        dsn  =-(sn(mn,2,1)*dx+sn(mn,2,2)*dy+sn(mn,2,3)*dz)
        dp   = dp   + (ps   + p        )*dsn
!
!.......flux d
!
        roec = 0.5*(v(nc+ncj,2)*v(nc+ncj,2)+v(nc+ncj,3)*v(nc+ncj,3) &
                   +v(nc+ncj,4)*v(nc+ncj,4))/v(nc+ncj,1)
        p    = gm1*(v(nc+ncj,5) -pinfl - roec)
        dsn  = (sn(mn+ncj,2,1)*dx+sn(mn+ncj,2,2)*dy+sn(mn+ncj,2,3)*dz)
        dp   = dp   + (ps   + p        )*dsn
!
!.......flux b
!
        roec = 0.5*(v(nc-nck,2)*v(nc-nck,2)+v(nc-nck,3)*v(nc-nck,3) &
                   +v(nc-nck,4)*v(nc-nck,4))/v(nc-nck,1)
        p    = gm1*(v(nc-nck,5) -pinfl - roec)
        dsn  =-(sn(mn,3,1)*dx+sn(mn,3,2)*dy+sn(mn,3,3)*dz)
        dp   = dp   + (ps   + p        )*dsn
!
!.......flux f
!
        roec = 0.5*(v(nc+nck,2)*v(nc+nck,2)+v(nc+nck,3)*v(nc+nck,3) &
                   +v(nc+nck,4)*v(nc+nck,4))/v(nc+nck,1)
        p    = gm1*(v(nc+nck,5) -pinfl - roec)
        dsn  = (sn(mn+nck,3,1)*dx+sn(mn+nck,3,2)*dy+sn(mn+nck,3,3)*dz)
        dp   = dp   + (ps   + p        )*dsn
!
!.......retrait de la contribution paroi
!
!....... -------
!....... |  s  |
!....... -------
!....... |  p  |
!....... -------
!
!.......flux p
!
        roec = 0.5*(v(nlc,2)*v(nlc,2)+v(nlc,3)*v(nlc,3) &
                   +v(nlc,4)*v(nlc,4))/v(nlc,1)
        p    = gm1*(v(nlc,5) -pinfl - roec)
        dsn  =-(sx*dx + sy*dy + sz*dz)
        dp   = dp   - (ps   + p        )*dsn
!
!.......prise en compte du volume
!
        rvol = 1.0/vol(nc)
        r2vol= 0.5 *rvol
        dp   = dp  *r2vol
        ds   = dsn *rvol
!
!.......valeur a la paroi
!
        pext = (ps   - dp  )/(1.0 + ds)
!
!.......1) ro = const  , pt=const , qt//qtint , p=pext
!
!       ro   = ros
!       pt   = pts
!       p    = pext
!
!       a2   = max(0.,gm*p/ro)
!       xm2  = tworgm1*((pt/p)**gm1rg - 1.0)
!       q2   = xm2*a2
!       ec   = 0.5*q2
!       coef = sqrt(q2/qt2s)
!       qtx  = coef*qtxs
!       qty  = coef*qtys
!       qtz  = coef*qtzs
!       rou  = ro*qtx
!       rov  = ro*qty
!       row  = ro*qtz
!       roe  = rgm1*p+ro*ec
!
!.......2) h = const , ent=const , qt//qtint , p=pext
!
        h    = hs
        ent  = ents
        p    = pext
!
        ro   = (p/ent)**rgm
        q2   = 2.0*(h - grgm1*p/ro)
        q2   = max(0.,q2)
        ec   = 0.5*q2
        coef = sqrt(q2/qt2s)
        qtx  = coef*qtxs
        qty  = coef*qtys
        qtz  = coef*qtzs
        rou  = ro*qtx
        rov  = ro*qty
        row  = ro*qtz
        roe  = p*rgm1 + pinfl + ro*ec
!
!.......valeur paroi
!
        v(nlc,1)    = ro
        v(nlc,2)    = rou
        v(nlc,3)    = rov
        v(nlc,4)    = row
        v(nlc,5)    = roe
!
        pression(nlc)=p
        ztemp(nlc)=gam*p/ro
        cson(nlc)=sqrt(ztemp(nlc))
        enddo
!
      return
      end
end module
