module mod_lp2kl3d
implicit none
contains
      subroutine lp2kl3d( &
                v,mu,mut,dist, &
                nxn,nyn,nzn, &
                ncin,ncbd,mfb,l, &
                vol,sn,ncyc, &
                mnpar,fgam, &
                tprod,utau,topz,  &
                dpdx,dpdy,dpdz, &
                ps,temp)
!
!***********************************************************************
!
!_DA  DATE_C : octobre 2000-- AUTEUR : Eric Goncalves
!
!     ACT
!_A    Lois de paroi - Approche de Smith - CAS 3D
!_A    discretisation sur deux points
!_A    - parois adiabatiques
!_A    - integration du profil de vitesse et de l'equation de l'energie
!_A      entre la paroi et le premier point du maillage
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
!
!***********************************************************************
!
      use para_var
      use para_fige
      use maillage
      use boundary
      use proprieteflu
      use definition
      use modeleturb
use mod_pgrad
implicit none
double precision :: v
double precision :: dist
integer :: ncin
integer :: ncbd
integer :: mfb
integer :: l
double precision :: vol
double precision :: sn
integer :: ncyc
integer :: mnpar
double precision :: fgam
double precision :: tprod
double precision :: utau
double precision :: topz
double precision :: dpdx
double precision :: dpdy
double precision :: dpdz
double precision :: ps
double precision :: temp
double precision :: ca
double precision :: cc3
double precision :: cklb2
double precision :: cklb3
double precision :: coefa
double precision :: coefb
double precision :: coefc
double precision :: coefc1
double precision :: coefp
double precision :: conv
double precision :: ctk
double precision :: ctmu
double precision :: ctmu2
double precision :: dconv
double precision :: dconvx
double precision :: dconvz
double precision :: dist2
double precision :: dpc
double precision :: dpde
double precision :: dpdt
double precision :: dy
double precision :: e1
double precision :: e2
double precision :: e3
double precision :: f1i
double precision :: fmui
integer :: ij
integer :: in
integer :: iter
integer :: j
integer :: jk
integer :: kk
integer :: lgsnlt
integer :: m
integer :: m0ns
integer :: mb
integer :: mpar
integer :: mt
integer :: n0c
integer :: nc
integer :: ndis
integer :: nfacns
integer :: ni
integer :: nii
integer :: nn
integer :: npsn
integer :: ntab
double precision :: rhoi
double precision :: rhol
double precision :: rk2
double precision :: rok2
double precision :: rop
double precision :: som1
double precision :: som2
double precision :: som3
double precision :: t1
double precision :: t2
double precision :: t3
double precision :: tn
double precision :: top
double precision :: toparx
double precision :: toparz
double precision :: topinix
double precision :: topiniz
double precision :: topx0
double precision :: topz0
double precision :: tt
double precision :: upyp1
double precision :: v1t
double precision :: v1x
double precision :: v1y
double precision :: v1z
double precision :: xi
double precision :: yi
double precision :: yp02
!
!-----------------------------------------------------------------------
!
      parameter( ntab=50  )
      logical lamin
      real mu,mut,mup,sv,ki,li,nxn,nyn,nzn
      real dtopx,dtopz,seuil
!
      dimension mu(ip12),mut(ip12)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
      dimension ncin(ip41),ncbd(ip41)
      dimension v(ip11,ip60),dist(ip12)
      dimension mnpar(ip12),fgam(ip42),utau(ip42)
      dimension temp(ip11),vol(ip11),ps(ip11),topz(ip11)
      dimension sn(ip31*ndir)
      dimension tprod(ip00),dpdx(ip00),dpdy(ip00),dpdz(ip00)
!
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: alfaa,betaa,ff,vitx,vitz,mui,muti,tempi,topcx,topcz
      ALLOCATE(alfaa(ntab),betaa(ntab),ff(ntab),topcx(2),topcz(2), &
            vitx(ntab),vitz(ntab),mui(ntab),muti(ntab),tempi(ntab))

      ndis=30
      dtopx=0.001
      dtopz=0.001
      seuil=1.e-6
      sv=110.4/tnz !air
!
!     constantes du calcul
      cklb2=cklb1**(4./3.)
      cklb3=cklb1**(1./3.)
      coefp=4.*sigmal/(vkar**2*cklb1)
!
      mt=mmb(mfb)
      m0ns=mpn(mfb)
      n0c=npc(l)
!
!    mise a zero des tableaux
      do in=1,ntab
        alfaa(in)=0.
        betaa(in)=0.
        tempi(in)=0.
        vitx(in)=0.
        vitz(in)=0.
      end do
!
!--------------------------------------------------------------
!-----initialisation de utau-----------------------------------
      if(ncyc.eq.icytur0) then
       do m=1,mt
        mb=mpb(mfb)+m
        ni=ncin(mb)
        nc=ncbd(mb)
        nfacns=m0ns+m
        mpar=mnpar(ni)
!       test sur transition et regime d'ecoulement
        if((fgam(mpar).lt.1.e-3).and.(ktransi.gt.0)) then
!         laminaire
          lamin=.true.
        else
!         turbulent
          lamin=.false.
        endif
!       vitesse cellule 1
        v1x=v(ni,2)/v(ni,1)
        v1y=v(ni,3)/v(ni,1)
        v1z=v(ni,4)/v(ni,1)
!       tangente normee a la paroi     
        tn=v1x*nxn(nfacns)+v1y*nyn(nfacns)+v1z*nzn(nfacns)
        t1=v1x-tn*nxn(nfacns)
        t2=v1y-tn*nyn(nfacns)
        t3=v1z-tn*nzn(nfacns)
        tt=sqrt(t1**2+t2**2+t3**2)
        t1=t1/tt
        t2=t2/tt
        t3=t3/tt
!       composante tangentielle de la vitesse dans repere paroi : v1t
        v1t=v1x*t1+v1y*t2+v1z*t3
!       masse volumique a la paroi
        rop=v(ni,1)*temp(ni)/temp(nc)
!       viscosite moleculaire a la paroi
        mup=mu(ni)*sqrt(temp(nc)/temp(ni))*(1.+sv/temp(ni))/(1.+sv/temp(nc))
!       contrainte de frottement a la paroi : top
        upyp1=rop*v1t*dist(ni)/mup
        yp02=yp0**2
        if(upyp1.le.yp02 .or. lamin) then
!        loi lineaire
         top=mup*v1t/dist(ni)
        else
!        loi logarithmique
         top=mup*v1t/dist(ni)
         do iter=1,10
          top=rop*v1t**2/(log(dist(ni)*sqrt(rop*top)/mup)/vkar+cllog)**2
         enddo
        endif
        utau(nfacns)=sign(1.,top)*sqrt(abs(top)/rop) 
       enddo !fin boucle sur facettes paroi
      endif
!-----fin initialisation de utau--------------------------------
!---------------------------------------------------------------
!
!     calcul du terme de derivee de la pression dp/dx a la paroi
      npsn=ndir*npfb(l)+1
      lgsnlt=nnn(l)
       call pgrad( &
          sn(npsn),vol,lgsnlt,l, &
          ps,dpdx,dpdy,dpdz)
!
!      boucle sur les facettes d'une frontiere paroi
       do m=1,mt
        mb=mpb(mfb)+m
        ni=ncin(mb)
        nc=ncbd(mb)
        nii=ni-n0c
        nfacns=m0ns+m
        mpar=mnpar(ni)
!       vitesse cellule 1
        v1x=v(ni,2)/v(ni,1)
        v1y=v(ni,3)/v(ni,1)
        v1z=v(ni,4)/v(ni,1)
!       tangente normee a la paroi
        tn=v1x*nxn(nfacns)+v1y*nyn(nfacns)+v1z*nzn(nfacns)
        t1=v1x-tn*nxn(nfacns)
        t2=v1y-tn*nyn(nfacns)
        t3=v1z-tn*nzn(nfacns)
        tt=sqrt(t1**2+t2**2+t3**2)
        t1=t1/tt
        t2=t2/tt
        t3=t3/tt
!       vecteur direction profondeur : e
        e1=nyn(nfacns)*t3-nzn(nfacns)*t2
        e2=nzn(nfacns)*t1-nxn(nfacns)*t3
        e3=nxn(nfacns)*t2-nyn(nfacns)*t1
!       composante tangentielle de la vitesse dans repere paroi : v1t
        v1t=v1x*t1+v1y*t2+v1z*t3
!       masse volumique a la paroi
        rop=v(ni,1)
!       terme de derivee de pression (dp/dt) a la paroi
        dpdt=dpdx(nii)*t1 + dpdy(nii)*t2 + dpdz(nii)*t3
!       terme de derivee de pression (dp/de) a la paroi
        dpde=dpdx(nii)*e1 + dpdy(nii)*e2 + dpdz(nii)*e3
!       calcul coefficients intermediaires
        coefc1=-cklb2*coefp*(dpdt**2+dpde**2)
        rk2=v(ni,6)**2
        dist2=dist(ni)**2
        dpc=dist(ni)**(sqrt(coefp))
!       calcul du topx a partir de utau
        toparx=rop*utau(nfacns)*abs(utau(nfacns))
!       initialisation topz
        toparz=topz(nc)
!       distance
        dy=dist(ni)/(ndis-1)
!c**************************************************************************
!c--------boucle Newton sur top---------------------------------------------
!c***************************************************************************
        dconv=1.
        do while (dconv.gt.seuil)
          topx0=toparx
          topz0=toparz
!         boucle pour calculer la derivee de la fonction pour le Newton
          do kk=1,2
           topinix=toparx
           topiniz=toparz
           conv=1.
           coefa=-cklb2*coefp*(toparx**2+toparz**2)
           coefb=-2.*cklb2*coefp*(toparx*dpdt+toparz*dpde)
           coefc=coefc1+4.*rk2/dist2
           ca=(4.*rk2+coefa/coefp-coefb*dist(ni)/(1.-coefp) &
                 -coefc*dist2/(4.-coefp)-4.*rk2/(4.-coefp))/dpc
!          viscosite a la paroi
           mup=mu(ni)
!          balayage pour calculet muti
           do nn=2,ndis
            yi=(nn-1)*dy
            mui(nn)=mu(ni)
            if((fgam(mpar).lt.1.e-3).and.(ktransi.gt.0)) then
             muti(nn)=0.
            else
             rhoi=v(ni,1)
             rok2=ca*yi**(sqrt(coefp))-coefa/coefp+coefb*yi/ &
                  (1.-coefp)+coefc*yi**2/(4.-coefp) + &
                   4.*rk2*yi**2/(dist2*(4.-coefp))
             ki=sqrt(rok2)/(2.*rhoi)
             li=vkar*yi+(v(ni,7)/v(ni,1)-vkar*dist(ni))*(yi**2/dist2)
             xi=rhoi*sqrt(2.*ki)*li/(mui(nn)*cklb3)
             f1i=exp(-50.*(li/(vkar*yi))**2)
             fmui=((25.5**4*f1i+2.**2*xi**2+xi**4)/ &
                   (25.5**4+2.**2*xi**2+xi**4))**0.25
             muti(nn)=mui(nn)*xi*fmui
            end if
            ff(nn)=muti(nn)/(mui(nn)+muti(nn))**2
           enddo
!---------------------------------------------------------------
!          integration de la vitesse U
!---------------------------------------------------------------
!          initialisation du calcul
           vitx(1)=0.
           vitx(ndis)=v1t
           ctk=-0.5*(mup+mui(2)+muti(2))/dy
           ctmu2=0.5*(mui(2)+muti(2)+mui(3)+muti(3))/dy
           betaa(ndis)=v1t
!          calcul des alfa, beta
           do ij=ndis-1,3,-1
             ctmu=0.5*(mui(ij)+muti(ij)+mui(ij+1)+muti(ij+1))/dy
             alfaa(ij)=(ctk+ctmu*alfaa(ij+1))/ctmu
             betaa(ij)=-(dpdt*(ij-1)*dy-ctmu*betaa(ij+1))/ctmu
           enddo
           vitx(2)=(dy*dpdt-ctmu2*betaa(3))/(-ctmu2+ctk+ctmu2*alfaa(3))
           do j=3,ndis-1
             vitx(j)=alfaa(j)*vitx(2)+betaa(j)
           enddo
!---------------------------------------------------------------
!          integration de la vitesse W
!--------------------------------------------------------------
!          initialisation du calcul
           vitz(1)=0.
           vitz(ndis)=0.
           do jk=1,ndis
             alfaa(jk)=0.
             betaa(jk)=0.
           enddo
!          calcul des alfa, beta
           do ij=ndis-1,3,-1
             ctmu=0.5*(mui(ij)+muti(ij)+mui(ij+1)+muti(ij+1))/dy
             alfaa(ij)=(ctk+ctmu*alfaa(ij+1))/ctmu
             betaa(ij)=-(dpde*(ij-1)*dy-ctmu*betaa(ij+1))/ctmu
           enddo
           vitz(2)=(dy*dpde-ctmu2*betaa(3))/(-ctmu2+ctk+ctmu2*alfaa(3))
           do j=3,ndis-1
             vitz(j)=alfaa(j)*vitz(2)+betaa(j)
           enddo
           mup=mu(ni)
           topcx(kk)=-ctk*vitx(2)-dy*0.5*dpdt
           topcz(kk)=-ctk*vitz(2)-dy*0.5*dpde
           toparx=topinix+dtopx
           toparz=topiniz+dtopz
          enddo
          toparx=topx0-(topcx(1)-topx0)*dtopx/(topcx(2)-topcx(1)-dtopx)
          toparz=topz0-(topcz(1)-topz0)*dtopz/(topcz(2)-topcz(1)-dtopz)
          dconvx=abs(topx0-toparx)
          dconvz=abs(topz0-toparz)
          dconv=max(dconvx,dconvz)
          dtopx=dconvx
          dtopz=dconvz
!************************************************************************
!------fin boucle du Newton----------------------------------------------
!************************************************************************
        enddo
!       test sur la transition
        if((fgam(mpar).lt.1.e-3).and.(ktransi.gt.0)) then
!        tprod(nii)=0.
         v(ni,7)=roeinf
         v(ni,6)=rokinf
        else
!        calcul de la production de k -> valeur moyenne sur la cellule adjacente
         som1=0.
         som2=0.
         som3=0.
         ff(1)=0.
         do kk=1,ndis-1
          som1=som1+(ff(kk)+ff(kk+1))*0.5*dy
          som2=som2+(ff(kk)*(kk-1)+ff(kk+1)*kk)*0.5*dy**2
          som3=som3+(ff(kk)*(kk-1)**2+ff(kk+1)*kk**2)*0.5*dy**3
         enddo
         cc1=toparx**2 + toparz**2
         cc2=2.*toparx*dpdt + 2.*toparz*dpde
         cc3=dpdt**2 + dpde**2
         tprod(nii)=(cc1*som1 + cc2*som2 + cc3*som3)/dist(ni)
!        calcul de l en cellule 1
         rhol=vkar*dist(ni)*v(ni,1)
         v(ni,7)=max(rhol,epse)
         v(nc,7)=0.
        endif
!       vitesse de frottement utau
        utau(nfacns)=sign(1.,toparx)*sqrt(abs(toparx)/rop)
!       topz
        topz(nc)=toparz
!       fin boucle sur facettes d'une frontiere paroi
       enddo

DEALLOCATE(alfaa,betaa,ff,vitx,vitz,mui,muti,tempi,topcx,topcz)

      return
      end
end module
