module mod_lpsa1
implicit none
contains
      subroutine lpsa1( &
                 v,mu,mut,dist, &
                 nxn,nyn,nzn, &
                 ncin,ncbd,mfb,l, &
                 mnpar,fgam,ncyc, &
                 temp)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 1999-- AUTEUR :  Eric Goncalves
!
!     ACT
!_A    Lois de paroi : modele de Spalart&Allmaras
!_A    - parois adiabatiques
!_A    - nut imposee a la paroi
!_A
!
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    u          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    cl         : com char(mtb       ) ; type de cond lim a appliquer
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    mmb        : com int (mtt       ) ; nombre de facettes d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!_I    mnpar      : arg real(ip12      ) ; pointeur dans tableaux front normales
!_I                                        stockees du point de rattach normale
!_I    fgam       : arg real(ip42      ) ; fonction d'intermittence pour
!_I                                        transition
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use boundary
      use proprieteflu
      use definition
      use modeleturb
implicit none
double precision :: v
double precision :: dist
integer :: ncin
integer :: ncbd
integer :: mfb
integer :: l
integer :: mnpar
double precision :: fgam
integer :: ncyc
double precision :: temp
double precision :: c1
double precision :: c13
double precision :: c18
double precision :: c2
double precision :: c3
double precision :: c4
double precision :: ca
double precision :: cb
double precision :: co
double precision :: dudy
integer :: iter
integer :: m
integer :: m0ns
integer :: mb
integer :: mpar
integer :: mt
integer :: n0c
integer :: nc
integer :: nfacns
integer :: ni
integer :: nii
double precision :: pka
double precision :: pkb
double precision :: rc4
double precision :: rnutilde
double precision :: rop
double precision :: t1
double precision :: t2
double precision :: t3
double precision :: temp1
double precision :: tn
double precision :: top
double precision :: tparoi
double precision :: tt
double precision :: upyp1
double precision :: v1t
double precision :: v1x
double precision :: v1y
double precision :: v1z
double precision :: yp02
double precision :: yplus1
!
!-----------------------------------------------------------------------
!
      logical lamin
      real mu,mut,mup,sv
      real nxn,nyn,nzn,n1,n2,n3
!
      dimension mu(ip12),mut(ip12)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
      dimension ncin(ip41),ncbd(ip41)
      dimension v(ip11,ip60),dist(ip12)
      dimension mnpar(ip12),fgam(ip42)
      dimension temp(ip11)
!
!     lois viscosite : sutherland pour gaz et loi exponentielle pour liquide
!     loi sutherland mu=mu0*sqrt(T/T0)*(1+S/T0)/(1+S/T)
!     pour vapeur d'eau S=548K, mu0=9.73e-6 Pa.s et T0=293K
!     loi exponentielle mu=A*exp(B/T)
!     pour l'eau A=1.24e-6 Pa.s et B=1968K
!      if(iflu.eq.1) then  !eau
!        sv=548./tnz
!        s0=tnz/293.
!        mu0=9.73E-6
!        al=1.214E-6
!        bl=1968./tnz
!      elseif(iflu.eq.2) then !freon R114
!     pour vapeur R114 : S=260K, mu0=10.527e-6 Pa.s et T0=293K
!     pour le R114: A=10.336e-6 Pa.s et B=976.9738K
!        sv=260./tnz
!        s0=tnz/293.
!        mu0=10.527E-6
!        al=10.336E-6
!        bl=976.9738/tnz
!      elseif(iflu.eq.3) then !air
!        sv=110.4/tnz !air
!        STOP 
!      endif
!
      sv=110.4/tnz !air
      c13=1./3.
      c18=1./18.
!
      mt=mmb(mfb)
      m0ns=mpn(mfb)
      n0c=npc(l)
!
!       boucle sur les facettes d'une frontiere paroi
        do m=1,mt
          mb=mpb(mfb)+m
          ni=ncin(mb)
          nc=ncbd(mb)
          nfacns=m0ns+m
          mpar=mnpar(ni)
          nii=ni-n0c
!         test sur transition et regime d'ecoulement
          if((fgam(mpar).lt.1.e-3).and.(ktransi.gt.0)) then
!           laminaire
            lamin=.true.
          else
!           turbulent
            lamin=.false.
          end if
!         vitesse cellule adjacente a la paroi (cellule 1)
          v1x=v(ni,2)/v(ni,1)
          v1y=v(ni,3)/v(ni,1)
          v1z=v(ni,4)/v(ni,1)
!         normale a la paroi
          n1=nxn(nfacns)
          n2=nyn(nfacns)
          n3=nzn(nfacns)
!         tangente normee a la paroi
          tn=v1x*n1+v1y*n2+v1z*n3
          t1=v1x-tn*n1
          t2=v1y-tn*n2
          t3=v1z-tn*n3
          tt=sqrt(t1**2+t2**2+t3**2)
          t1=t1/tt
          t2=t2/tt
          t3=t3/tt
!         composante tangentielle de la vitesse dans repere paroi : v1t
          v1t=v1x*t1+v1y*t2+v1z*t3
!         temperature cellule 1 : temp1
          temp1=temp(ni)
!         temperature a la paroi : tparoi
          pka=cp*(mu(ni)/pr+mut(ni)/prt)
          pkb=mu(ni)+mut(ni)
          tparoi=temp1+0.5*pkb*v1t**2/pka
!          tparoi=temp(nc)
!         viscosite moleculaire a la paroi
          mup=mu(ni)*sqrt(tparoi/temp1)*(1.+sv/temp1)/(1.+sv/tparoi)
!         masse volumique a la paroi
          rop=v(ni,1)*temp1/tparoi
!         correction de compressibilite (loi de Van Driest)
          co=sqrt(2.*pka*tparoi/pkb)
          v1t=co*asin(v1t/co)
!         contrainte de frottement a la paroi : top
          upyp1=rop*v1t*dist(ni)/mup
          yp02=yp0**2
!         loi standard
          if(upyp1.le.yp02 .or. lamin) then
!           loi lineaire
            top=mup*v1t/dist(ni)
            dudy=top/mup
          else
!           loi logarithmique
            top=mup*v1t/dist(ni)
            do iter=1,10
              top=rop*v1t**2/(log(dist(ni)*sqrt(rop*top)/mup)/vkar+cllog)**2
            enddo
            dudy=sqrt(top/rop)/(vkar*dist(ni))
          endif
!
!        loi raffinee avec interpolation
!         if(upyp1.le.9.) then
!          loi lineaire
!           top=max(1.e-10,mup*v1t/dist(ni))
!       else
!         yplus0=3.
!         do ii=1,10
!           yplus0=upyp1/(log(yplus0)/vkar+cllog)
!          end do
!           if(yplus0.gt.40.) then
!            loi logarithmique
!           top=max(1.e-10,mup*v1t/dist(ni))
!           do iter=1,10
!             top=max(1.e-10,rop)*v1t**2/(log(dist(ni)*sqrt(rop &
!                     *top)/mup)/vkar+cllog)**2)
!            end do
!           else
!            region tampon : interpolation par polynome de degre 4
!           top=max(1.e-10,mup*v1t/dist(ni))
!           do jj=1,10
!             yy=log(dist(ni)*sqrt(rop*top)/mup)
!              top=max(1.e-10,rop*v1t** 2/(0.17962*yy**4- & 
!                       2.2117*yy**3+9.2052*yy**2-10.804*yy &
!                   +6.4424)**2)
!            end do
!           end if
!        end if
!
!          if(lamin) then
!          loi lineaire
!           top=mup*v1t/dist(ni)
!           dudy=top/mup
!          else
!          loi de Spalding
!           top=mup*v1t/dist(ni)
!           upl=v1t*sqrt(rop/top)
!           do jj=1,15
!            fu=rop*dist(ni)*v1t/(mup*upl) - upl - exp(-vkar*cllog)*
!     &         (exp(vkar*upl) -1. -vkar*upl - 0.5*(vkar*upl)**2 -
!     &          ((vkar*upl)**3)/6.)
!            dfu=-rop*dist(ni)*v1t/(mup*upl**2)-1.-exp(-vkar*cllog)*
!     &          vkar*(exp(vkar*upl) -1. -vkar*upl - 0.5*(vkar*upl)**2)
!            upl=upl-fu/dfu
!           enddo
!           top=rop*(v1t/upl)**2
!           dudy=sqrt(top/rop)/(vkar*dist(ni))
!          endif
!         yplus cellule 1
          yplus1=dist(ni)*sqrt(top*rop)/mup
!         calcul de nu_tilde en cellule 1
          mut(ni)=v(ni,1)*vkar**2*dist(ni)**2*dudy* &
                    (1.-exp(-yplus1/26.))**2
          ca=mut(ni)
          cb=mu(ni)**3*mut(ni)*cv1**3
          c1=cb*sqrt(3.)*sqrt(256.*cb+27.*ca**4)
          c2=0.5*cb*ca**2+c18*c1
          c3=-0.5*cb*ca**2+c18*c1
          c4=4.*abs(c3**c13-c2**c13)
          rc4=sqrt(ca**2+c4)
          rnutilde=0.25*(ca+rc4)+0.25*sqrt(2.)*sqrt(ca**2*rc4-2.*c3**c13*rc4 &
             +2.*c2**c13*rc4+ca**3)/rc4**0.5
!         ATTENTION! si probleme en maillage fin, utiliser nu=kappa*dist*utau
!          rnutilde=v(ni,1)*kappa*dist(ni)*sqrt(top/rop)
          v(ni,6)=max(rnutilde,epsk)
          v(nc,6)=0.
!       fin boucle sur facettes paroi
        enddo
!
      return
      end subroutine

end module
