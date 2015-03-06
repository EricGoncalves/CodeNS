      subroutine rbvr( &
                 ncbd,mnr,xnr,ynr,znr, &
                 t, &
                 ncin)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables (t) sur des facettes frontieres recouvertes
!_A    (continuite par interpolations et moyennes).
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
!     I/O
!_/    t          : arg real(ip11,ip60 ) ; variables de calcul
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use boundary
      use maillage
!
!-----------------------------------------------------------------------
!
      dimension t(ip11,ip60)
      dimension ncbd(ip41),ncin(ip41)
      dimension mnr(ip44),xnr(ip44),ynr(ip44),znr(ip44)
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
!DEC$ IVDEP
      do m=1,mt
      mr=mpr(mfbm)+m
      nr=mnr(mr)
      ml=mpb(mfbm)+m
      nd=ncbd(ml)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
      t(nd,1)= &
           ((1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*t(nr              ,1)+ &
            (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *t(nr      +nid*njd,1)+ &
            (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*t(nr  +nid        ,1)+ &
            (1-xnr(mr))*   ynr(mr) *   znr(mr) *t(nr  +nid+nid*njd,1)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*t(nr+1            ,1)+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *t(nr+1    +nid*njd,1)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*t(nr+1+nid        ,1)+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *t(nr+1+nid+nid*njd,1) &
           +t(ncin(ml),1))*.5
      t(nd,2)= &
           ((1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*t(nr              ,2)+ &
            (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *t(nr      +nid*njd,2)+ &
            (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*t(nr  +nid        ,2)+ &
            (1-xnr(mr))*   ynr(mr) *   znr(mr) *t(nr  +nid+nid*njd,2)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*t(nr+1            ,2)+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *t(nr+1    +nid*njd,2)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*t(nr+1+nid        ,2)+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *t(nr+1+nid+nid*njd,2) &
           +t(ncin(ml),2))*.5
      t(nd,3)= &
            (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*t(nr              ,3)+ &
            (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *t(nr      +nid*njd,3)+ &
            (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*t(nr  +nid        ,3)+ &
            (1-xnr(mr))*   ynr(mr) *   znr(mr) *t(nr  +nid+nid*njd,3)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*t(nr+1            ,3)+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *t(nr+1    +nid*njd,3)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*t(nr+1+nid        ,3)+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *t(nr+1+nid+nid*njd,3)
      t(nd,4)= &
            (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*t(nr              ,4)+ &
            (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *t(nr      +nid*njd,4)+ &
            (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*t(nr  +nid        ,4)+ &
            (1-xnr(mr))*   ynr(mr) *   znr(mr) *t(nr  +nid+nid*njd,4)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*t(nr+1            ,4)+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *t(nr+1    +nid*njd,4)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*t(nr+1+nid        ,4)+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *t(nr+1+nid+nid*njd,4)
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
      aux =t(nd,3)
      t(nd,3) =aux*cr-t(nd,4)*sr
      t(nd,4) =t(nd,4)*cr+aux*sr
!
      t(nd,3)=(t(nd,3)+t(ncin(ml),3))*.5
      t(nd,4)=(t(nd,4)+t(ncin(ml),4))*.5
!
      t(nd,5)= &
           ((1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*t(nr              ,5)+ &
            (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *t(nr      +nid*njd,5)+ &
            (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*t(nr  +nid        ,5)+ &
            (1-xnr(mr))*   ynr(mr) *   znr(mr) *t(nr  +nid+nid*njd,5)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*t(nr+1            ,5)+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *t(nr+1    +nid*njd,5)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*t(nr+1+nid        ,5)+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *t(nr+1+nid+nid*njd,5) &
           +t(ncin(ml),5))*.5
!
      enddo
      enddo
!
      return
      end
