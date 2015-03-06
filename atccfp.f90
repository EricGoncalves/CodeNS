        subroutine atccfp( &
                 x,y,z, &
                 xpar,ypar,zpar, &
                 ncin, &
                 nfbi,dm1,dm2,dm3,isens3)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des coordonnees des centres des facettes formant les parois
!
!     INP
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontie
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    nfbi       : arg int              ; numero interne de frontiere
!
!     OUT
!_O    xpar       : arg real(ip00      ) ; coordonnee x centre facette paroi
!_O    ypar       : arg real(ip00      ) ; coordonnee y centre facette paroi
!_O    zpar       : arg real(ip00      ) ; coordonnee z centre facette paroi
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
	  use boundary
!
!-----------------------------------------------------------------------
!
      integer dm1,dm2,dm3
      dimension x(ip21),y(ip21),z(ip21)
      dimension xpar(ip00),ypar(ip00),zpar(ip00)
      dimension ncin(ip41)
!
      m0 =mpb(nfbi)
      m0n=mpn(nfbi)
      mt=mmb(nfbi)
      if(isens3.eq.-1) then
!       paroi en m3max
        ndm3=dm3
      else
!       paroi en m3min
        ndm3=0
      end if
!
!     boucle sur tous les points de la paroi
      do m=1,mt
        mb =m0 +m
        mbn=m0n+m
        nci=ncin(mb)+ndm3
        xpar(mbn)=(x(nci)+x(nci+dm1)+x(nci+dm2)+x(nci+dm1+dm2))*.25
        ypar(mbn)=(y(nci)+y(nci+dm1)+y(nci+dm2)+y(nci+dm1+dm2))*.25
        zpar(mbn)=(z(nci)+z(nci+dm1)+z(nci+dm2)+z(nci+dm1+dm2))*.25
       enddo
!
      return
      end
