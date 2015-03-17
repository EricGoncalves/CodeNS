module mod_utdebit
  implicit none
contains
  subroutine utdebit( &
       mfb,bceqt,pres, &
       nxn,nyn,nzn,ncbd,v, &
       mmb,mpb,mpn,lm, &
       sn,lgsnlt)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 2006 -- AUTEUR : E. Goncalves / LEGI
!
!     ACT
!_A    Sous-programme utilisateur de preparation des donnees pour le
!_A    sous-programme cldebit.
!_A    Calcul du debit global correspondant l'etat schema par
!_A    integration (methode des trapezes)
!_A    ATTENTION: le debit est stocke dans le tableau 'pres'.
!
!_I    mfb        : arg int              ; numero de la frontiere
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!
!     OUT
!_O    pres       : arg real(ip40      ) ; debit local rho*vitesse normale
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::      i,    i1,    i2,    id,   inc
    integer          ::   ind1,  ind2,  indc,     j,    j1
    integer          ::     j2,    jd,     k,    k1,    k2
    integer          ::     kd,  kdir,lgsnlt,    lm,     m
    integer          ::     mb,   mfb,   mmb,    mn,   mpb
    integer          ::    mpn,    mt,     n,   n0c,  ncbd
    integer          ::    nci,   ncj,   nck,   nid,  nijd
    integer          ::   ninc,   njd,    nl
    double precision :: bceqt,  dqm,  nxn,  nyn,  nzn
    double precision ::  pres,  qms, rqns,   sn,    v
!
!-----------------------------------------------------------------------
!
!
    dimension v(ip11,ip60)
    dimension bceqt(ip41,neqt)
    dimension pres(ip40),ncbd(ip41),nxn(ip42),nyn(ip42),nzn(ip42)
    dimension mmb(mtt),mpb(mtt),mpn(mtt)
    dimension sn(lgsnlt,nind,ndir)
!
    indc(i,j,k)=n0c+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
    inc(id,jd,kd)=id+jd*nid+kd*nijd
!
    n0c=npc(lm)
    i1=ii1(lm)
    i2=ii2(lm)
    j1=jj1(lm)
    j2=jj2(lm)
    k1=kk1(lm)
    k2=kk2(lm)
!
    nid = id2(lm)-id1(lm)+1
    njd = jd2(lm)-jd1(lm)+1
    nijd= nid*njd
!
    nci=inc(1,0,0)
    ncj=inc(0,1,0)
    nck=inc(0,0,1)
!
    qms=0.
    kdir=1
    ninc=ncj
!
!     calcul du debit global etat schema
    do k=k1,k2-1
       ind1 = indc(i2,j1,k)
       ind2 = indc(i2,j2,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          dqm= (v(n,2)+v(n-ninc,2))*sn(m,kdir,1) &
               +(v(n,3)+v(n-ninc,3))*sn(m,kdir,2) &
               +(v(n,4)+v(n-ninc,4))*sn(m,kdir,3)
          qms=qms+dqm
       enddo
    enddo
!
    mt=mmb(mfb)
    do m=1,mt
       mb=mpb(mfb)+m
       mn=mpn(mfb)+m
       nl=ncbd(mb)
       rqns=-v(nl,2)*nxn(mn)-v(nl,3)*nyn(mn)-v(nl,4)*nzn(mn)   !attentio
       pres(m)=(bceqt(mb,1)/qms)*rqns
    enddo
!
    return
  end subroutine utdebit
end module mod_utdebit
