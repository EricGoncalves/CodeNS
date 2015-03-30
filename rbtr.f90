module mod_rbtr
  implicit none
contains
  subroutine rbtr( &
       ncbd,mnr,xnr,ynr,znr, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       ncin)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des composantes du tenseur visqueux et des flux de chaleur
!_A    sur des facettes frontieres recouvertes (continuite par interpolations et
!_A    moyennes).
!
!     INP
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mnr        : arg int (ip44      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule recouvrante
!_I    xnr        : arg real(ip44      ) ; coefficient d'interpolation en x
!_I                                        dans une cellule recouvrante
!_I    ynr        : arg real(ip44      ) ; coefficient d'interpolation en y
!_I                                        dans une cellule recouvrante
!_I    znr        : arg real(ip44      ) ; coefficient d'interpolation en z
!_I                                        dans une cellule recouvrante
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpr        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front recouvertes
!_I    ndrr       : com int (mtb       ) ; numero de domaine recouvrant
!_I    srotr      : com real(mtb       ) ; rotation amenant une front periodique
!_I                                        dans sa position recouverte, sinus
!_I    crotr      : com real(mtb       ) ; rotation amenant une front periodique
!_I                                        dans sa position recouverte, cosinus
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!     OUT
!
!     I/O
!_/    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_/                                        contraintes visqueuses
!_/    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_/                                        contraintes visqueuses
!_/    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_/                                        contraintes visqueuses
!_/    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_/                                        contraintes visqueuses
!_/    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_/                                        contraintes visqueuses
!_/    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_/                                        contraintes visqueuses
!_/    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_/    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_/    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use maillage
    implicit none
    integer          ::          l,         m,        mb,        mf,       mfb
    integer          ::       mfbm, mnr(ip44),        mr,        mt,ncbd(ip41)
    integer          :: ncin(ip41),        nd,       nid,       njd,        nr
    double precision ::         cr, qcx(ip12),      qcxr, qcy(ip12),      qcyr
    double precision ::  qcz(ip12),      qczr,        sr,toxx(ip12),toxy(ip12)
    double precision :: toxz(ip12),toyy(ip12),toyz(ip12),tozz(ip12),      txxr
    double precision ::       txyr,      txzr,      tyyr,      tyzr,      tzzr
    double precision ::  xnr(ip44), ynr(ip44), znr(ip44)
!
!-----------------------------------------------------------------------
!
!
!
    do mf=1,nbd
!
       mfbm=lbd(mf)
       mt=mmb(mfbm)
!
       mfb =mod(mfbm,mtb)
       if (mfb.eq.0) mfb=mtb
       sr =-srotr(mfb)
       cr = crotr(mfb)
       l  = ndrr(mfb)
!
       nid=id2(l)-id1(l)+1
       njd=jd2(l)-jd1(l)+1
!
!!$OMP SIMD
!$OMP SIMD
       do m=1,mt
          mr=mpr(mfbm)+m
          nr=mnr(mr)
          mb=mpb(mfbm)+m
          nd=ncbd(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
          toxx(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*toxx(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *toxx(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*toxx(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *toxx(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*toxx(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *toxx(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*toxx(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *toxx(nr+1+nid+nid*njd)
          toyy(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*toyy(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *toyy(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*toyy(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *toyy(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*toyy(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *toyy(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*toyy(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *toyy(nr+1+nid+nid*njd)
          tozz(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*tozz(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *tozz(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*tozz(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *tozz(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*tozz(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *tozz(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*tozz(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *tozz(nr+1+nid+nid*njd)
          toxy(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*toxy(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *toxy(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*toxy(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *toxy(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*toxy(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *toxy(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*toxy(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *toxy(nr+1+nid+nid*njd)
!
          toxz(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*toxz(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *toxz(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*toxz(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *toxz(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*toxz(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *toxz(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*toxz(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *toxz(nr+1+nid+nid*njd)
!
          toyz(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*toyz(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *toyz(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*toyz(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *toyz(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*toyz(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *toyz(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*toyz(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *toyz(nr+1+nid+nid*njd)
!
          qcx(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))* qcx(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) * qcx(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))* qcx(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) * qcx(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))* qcx(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) * qcx(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))* qcx(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) * qcx(nr+1+nid+nid*njd)
          qcy(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))* qcy(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) * qcy(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))* qcy(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) * qcy(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))* qcy(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) * qcy(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))* qcy(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) * qcy(nr+1+nid+nid*njd)
          qcz(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))* qcz(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) * qcz(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))* qcz(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) * qcz(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))* qcz(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) * qcz(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))* qcz(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) * qcz(nr+1+nid+nid*njd)
!
!  Attention : le sens de rotation est le sens direct.
!  ----------
!
!                  | 1      0      0 |
!          P^-1  = | 0    cos   -sin |
!                  | 0    sin    cos |
!
!      | x |              | x |
!      | y |     = P^-1 * | y |
!      | z |new           | z |old
!
!
!        A       = P^-1 * A * P
!
          txxr=toxx(nd)
          txyr=cr*toxy(nd)-sr*toxz(nd)
          txzr=cr*toxz(nd)+sr*toxy(nd)
          tyyr=cr*cr*toyy(nd)-2.*cr*sr*toyz(nd)+sr*sr*tozz(nd)
          tyzr=-cr*sr*(tozz(nd)-toyy(nd))+(2.*cr*cr-1.)*toyz(nd)
          tzzr=sr*sr*toyy(nd)+2.*cr*sr*toyz(nd)+cr*cr*tozz(nd)
          qcxr=qcx(nd)
          qcyr=cr*qcy(nd)-sr*qcz(nd)
          qczr=cr*qcz(nd)+sr*qcy(nd)
!
          toxx(nd)=(txxr+toxx(ncin(mb)))*.5
          toxy(nd)=(txyr+toxy(ncin(mb)))*.5
          toxz(nd)=(txzr+toxz(ncin(mb)))*.5
          toyy(nd)=(tyyr+toyy(ncin(mb)))*.5
          toyz(nd)=(tyzr+toyz(ncin(mb)))*.5
          tozz(nd)=(tzzr+tozz(ncin(mb)))*.5
          qcx(nd)=(qcxr+ qcx(ncin(mb)))*.5
          qcy(nd)=(qcyr+ qcy(ncin(mb)))*.5
          qcz(nd)=(qczr+ qcz(ncin(mb)))*.5
!
       enddo
    enddo
!
    return
  end subroutine rbtr
end module mod_rbtr
