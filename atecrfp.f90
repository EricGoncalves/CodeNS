        subroutine atecrfp( &
                 x,y,z, &
                 xpar,ypar,zpar, &
                 nxn,nyn,nzn, &
                 ncin, &
                 nfbi)
!
!***********************************************************************
!
!     ACT
!_A    ecriture des coordonnees des centres des facettes formant les parois
!
!     INP
!_I    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une frontiere
!_I    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une frontiere
!_I    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontie
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    nfbi       : arg int              ; numero interne de frontiere
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nfbn       : com int (mtb       ) ; numero dans numerotation interne
!
!     OUT
!_O    xpar       : arg real(ip00      ) ; coordonnee x centre facette paroi
!_O                                        dans tab normales stockees
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
      real nxn,nyn,nzn
!
      dimension x(ip21),y(ip21),z(ip21)
      dimension xpar(ip00),ypar(ip00),zpar(ip00)
      dimension ncin(ip41)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
!
      l=ndlb(nfbi)
      m0b=mpb(nfbi)
      m0n=mpn(nfbi)
!
      iminf=iminb(nfbi)
      imaxf=imaxb(nfbi)
      jminf=jminb(nfbi)
      jmaxf=jmaxb(nfbi)
      kminf=kminb(nfbi)
      kmaxf=kmaxb(nfbi)
!
      if (iminf.eq.imaxf) then
        m1max=jmaxf-jminf+1
        m2max=kmaxf-kminf+1
       elseif (jminf.eq.jmaxf) then
        m1max=imaxf-iminf+1
        m2max=kmaxf-kminf+1
      elseif (kminf.eq.kmaxf) then
        m1max=imaxf-iminf+1
        m2max=jmaxf-jminf+1
      end if
      idm=m1max-1
!
!     boucle sur toutes les facettes de la paroi
      do m2=1,m2max-1
        do m1=1,m1max-1
          mfac =m0b+m1+(m2-1)*idm
          mfacn=m0n+m1+(m2-1)*idm
        enddo
      enddo
!
      return
      end
