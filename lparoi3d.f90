      subroutine lparoi3d(  &
                 ncyc, &
                 nxn,nyn,nzn, &
                 ncin,ncbd,mfb, &
                 toxx,toxy,toxz,toyy,toyz,tozz, &
                 qcx,qcy,qcz, &
                 v,utau,temp,topz)
!
!***********************************************************************
!
!_DA  DATE_C : Eric Goncalves
!
!     ACT
!_A    Lois de paroi - Approche de Smith - CAS 3D
!_A    - parois adiabatiques
!_A    - integration du profil de vitesse entre la paroi et 
!_A      le premier point du maillage
!_A    - contrainte de frottement imposee a la paroi
!_A    - flux de chaleur imposee a la paroi
!_A
!
!     INP
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    npfb       : com int (lt        ) ; pointeur fin de domaine precedent
!_I                                        dans tableau toutes facettes
!_I    nnn        : com int (lt        ) ; nombre de noeuds du domaine (dont fictif)
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_I                                        contraintes visqueuses
!_I    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_I                                        contraintes visqueuses
!_I    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_I                                        contraintes visqueuses
!_I    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_I                                        contraintes visqueuses
!_I    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_I                                        contraintes visqueuses
!_I    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_I                                        contraintes visqueuses
!_I    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_I    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_I    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
!_I    img        : com int              ; niveau de grille (multigrille)
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    mmb        : com int (mtt       ) ; nombre de facettes d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    npfb       : com int (lt        ) ; pointeur fin de domaine precedent
!_I                                        dans tableau toutes facettes
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!_I    utau       : arg real (ip 42    ) : sert a stocker topx
!
!
!***********************************************************************
!
      use para_var
      use para_fige
      use maillage
      use boundary
      use proprieteflu
      use modeleturb
      use definition
!
!-----------------------------------------------------------------------
!
      real nxn,nyn,nzn,n1,n2,n3,e1,e2,e3
!
      dimension v(ip11,ip60)
      dimension toxx(ip12),toxy(ip12),toxz(ip12), &
                toyy(ip12),toyz(ip12),tozz(ip12), &
                qcx(ip12),qcy(ip12),qcz(ip12)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
      dimension ncin(ip41),ncbd(ip41)
      dimension utau(ip42),temp(ip11),topz(ip11)
!
      mt=mmb(mfb)
      m0ns=mpn(mfb)
!
      do m=1,mt
        mb=mpb(mfb)+m
        ni=ncin(mb)
        nc=ncbd(mb)
        nfacns=m0ns+m
!       vitesse cellule 1
        v1x=v(ni,2)/v(ni,1)
        v1y=v(ni,3)/v(ni,1)
        v1z=v(ni,4)/v(ni,1)
!       normale a la paroi
        n1=nxn(nfacns)
        n2=nyn(nfacns)
        n3=nzn(nfacns)
!       tangente normee a la paroi
        tn=v1x*n1 + v1y*n2 + v1z*n3
        t1=v1x-tn*n1
        t2=v1y-tn*n2
        t3=v1z-tn*n3
        tt=sqrt(t1**2+t2**2+t3**2)
        t1=t1/tt
        t2=t2/tt
        t3=t3/tt
!       vecteur direction profondeur : e = t vectoriel n
        e1=n2*t3-n3*t2
        e2=n3*t1-n1*t3
        e3=n1*t2-n2*t1
!       masse volumique a la paroi
        rop=v(ni,1)*temp(ni)/temp(nc)
!       calcul du topx
        topx=rop*utau(nfacns)*abs(utau(nfacns))
!       tenseur des contraintes a la paroi dans repere general
        toxx(nc)=2*t1*n1*topx + 2.*e1*n1*topz(nc)
        toyy(nc)=2*t2*n2*topx + 2.*e2*n2*topz(nc)
        tozz(nc)=2*t3*n3*topx + 2.*e3*n3*topz(nc)
        toxy(nc)=(t2*n1+t1*n2)*topx+(e1*n2+n1*e2)*topz(nc)
        toxz(nc)=(t1*n3+t3*n1)*topx+(e1*n3+n1*e3)*topz(nc)
        toyz(nc)=(t2*n3+t3*n2)*topx+(e2*n3+n2*e3)*topz(nc)
!       tenseur des contraintes en 1 dans repere general
        toxx(ni)=toxx(nc)
        toyy(ni)=toyy(nc)
        tozz(ni)=tozz(nc)
        toxy(ni)=toxy(nc)
        toxz(ni)=toxz(nc)
        toyz(ni)=toyz(nc)
!       flux de chaleur dans cellule 1 dans repere general
!       ATTENTION! le signe 'moins' provient de la convention utilisee
        qc1=topx*(v1x*t1+v1y*t2+v1z*t3)+topz(nc)*(v1x*e1+v1y*e2+v1z*e3)
        qcx(ni)=-n1*qc1
        qcy(ni)=-n2*qc1
        qcz(ni)=-n3*qc1
!       flux de chaleur a la paroi dans repere general
        qcx(nc)=0.
        qcy(nc)=0.
        qcy(nc)=0.
!     fin boucle sur facettes d'une frontiere paroi
      end do
!
      return
      end
