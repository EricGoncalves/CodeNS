          subroutine at_fidist( &
                 igr,jgr,kgr,raptat, &
                 x,y,z, &
                 xpar,ypar,zpar, &
                 xcc,ycc,zcc,dist2, &
                 dist,mnpar, &
                 m1tb,m2tb,nfrtb, &
                 l)
!
!***********************************************************************
!
!     ACT
!_A    Distance des cellules aux noeuds du maillage contenu dans
!_A    un maillage grossier par exploration de la portion de frontiere
!_A    limitee par les valeurs extremes des indices parois rattaches
!_A    aux noeuds maillage grossier. Traitement d'un domaine.
!
!     INP
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    xpar       : arg real(ip00      ) ; coordonnee x centre facette paroi
!_I    ypar       : arg real(ip00      ) ; coordonnee y centre facette paroi
!_I    zpar       : arg real(ip00      ) ; coordonnee z centre facette paroi
!_I    xcc        : arg real(ip00      ) ; coordonnee x centre cellule
!_I    ycc        : arg real(ip00      ) ; coordonnee y centre cellule
!_I    zcc        : arg real(ip00      ) ; coordonnee z centre cellule
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero interne de front a traiter
!_I    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une frontiere
!_I    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une frontiere
!_I    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une frontiere
!_I    igr        : arg int              ; pas en i pour definition maillage
!_I                                        grossier
!_I    jgr        : arg int              ; pas en j pour definition maillage
!_I                                        grossier
!_I    kgr        : arg int              ; pas en k pour definition maillage
!_I                                        grossier
!_I    raptat     : arg int              ; parametre definition portion
!_I                                        contigue de frontiere
!_I    nbdrat     : com int (lz        ) ; nb parois rattachees au domaine l
!_I    lbdrat     : com int (mtb       ) ; no interne des parois rattachees
!_I                                        au domaine l
!_I    npbrat     : com int (lz        ) ; pointeur fin liste frontieres a
!_I                                        traiter pour domaine precedent
!
!     OUT
!_O    dist       : arg real(ip12      ) ; distance a la paroi
!_O    mnpar      : arg real(ip12      ) ; pointeur tablaux front normales
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
   use sortiefichier
   use maillage
   use boundary
!
!-----------------------------------------------------------------------
!
      logical nondeg,degi,degj,degk,contig
!
      dimension x(ip21),y(ip21),z(ip21)
      dimension xpar(ip00),ypar(ip00),zpar(ip00)
      dimension xcc (ip00),ycc (ip00),zcc (ip00),dist2(ip00)
      dimension dist(ip12),mnpar(ip12)
      dimension m1tb(ip00),m2tb(ip00),nfrtb(ip00)
      dimension nc(8)
!
      ind(i,j,k)=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
!
      n0=npn(l)
      i1=ii1(l)
      j1=jj1(l)
      k1=kk1(l)
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
!
      nijd = nid*njd
      i2=ii2(l)
      j2=jj2(l)
      k2=kk2(l)
!
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
      degi  = i2-1.eq.i1
      degj  = j2-1.eq.j1
      degk  = k2-1.eq.k1
      nondeg= (.not. degi) .and. (.not. degj) .and. (.not. degk)
      if(nondeg)then
        nsom=8
      else
        nsom=4
      end if
!
      i2m2=i2-2
      j2m2=j2-2
      k2m2=k2-2
      if(degj) j2m2=j1
      if(degk) k2m2=k1
!
!     ncelgr : compteur cellules maillage grossier
!     ncelat : compteur cellules maillage grossier rattachees
!
      ncelgr=0
      ncelat=0
!
!               nfpar0 : numero de la paroi de rattachement. Evite de
!                        recalculer pour chaque maille maillage grossier
!                        la correspondance "m1,m2"<->pointeur
!                 igrp : indice i noeud suivant en i dans maillage grossier
!                 jgrp : indice j noeud suivant en j dans maillage grossier
!                 kgrp : indice k noeud suivant en k dans maillage grossier
!
      nfpar0=0
!
!     boucle sur les cellules interieures au maillage grossier,
!     premieres et dernieres lignes comprises
!
      do kg=k1,k2m2,kgr
        kgrp=kg+kgr
        if(kgrp.gt.k2m2) kgrp=k2m2-kg
        do jg=j1,j2m2,jgr
          jgrp=jgr
          if(jg+jgrp.gt.j2m2) jgrp=j2m1-jg
          do ig=i1,i2m2,igr
!           boucles dans le maillage grossier
!
            igrp=igr
            if(ig+igrp.gt.i2m2) igrp=i2m1-ig
!
            if(nondeg) then
!
!             maillage non degenere
!
              idind=igrp
              jdind=jgrp*nid
              kdind=kgrp*nijd
              nc(1)=ind(ig,jg,kg)
              nc(2)=nc(1)+idind
              nc(3)=nc(2)+jdind
              nc(4)=nc(1)+jdind
              nc(5)=nc(1)+kdind
              nc(6)=nc(2)+kdind
              nc(7)=nc(3)+kdind
              nc(8)=nc(4)+kdind
            else
              if(degj) then

                idind=igrp
                kdind=kgrp*nijd
                nc(1)=ind(ig,jg,kg)
                nc(2)=nc(1)+idind
                nc(3)=nc(2)+kdind
                nc(4)=nc(1)+kdind
!
              else if(degk) then
!
                idind=igrp
                jdind=jgrp*nid
                nc(1)=ind(ig,jg,kg)
                nc(2)=nc(1)+idind
                nc(3)=nc(2)+jdind
                nc(4)=nc(1)+jdind
!
              end if
            end if
!
!           les cellules aux sommets de la maille du maillage grossier
!           sont-elles toutes rattachees a une meme paroi pour former
!           un element de surface contigue?
!
            idif=0
            mfacp=mnpar(nc(1))
            nfpar=nfrtb(mfacp)
            m1mi =m1tb (mfacp)
            m2mi =m2tb (mfacp)
            m1mx =m1tb (mfacp)
            m2mx =m2tb (mfacp)
            do ns=2,nsom
              mfacp=mnpar(nc(ns))
              if(nfpar.eq.nfrtb(mfacp) ) then
                m1mi =min( m1mi,m1tb(mfacp) )
                m2mi =min( m2mi,m2tb(mfacp) )
                m1mx =max( m1mx,m1tb(mfacp) )
                m2mx =max( m2mx,m2tb(mfacp) )
              else
                idif=idif+1
              end if
            end do
!
            if(idif.eq.0) then
!             les facettes de rattachement forment-elles un element contigue ?
!
              if(nfpar.ne.nfpar0) then
!
!               changement de paroi de rattachement.
!               Recalcul des parametres
!
                nfpar0=nfpar
                m0b=mpb(nfpar)
                m0n=mpn(nfpar)
!
                iminf=iminb(nfpar)
                imaxf=imaxb(nfpar)
                jminf=jminb(nfpar)
                jmaxf=jmaxb(nfpar)
                kminf=kminb(nfpar)
                kmaxf=kmaxb(nfpar)
!
                if(iminf.eq.imaxf) then
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
                iespacem=int( float(m1max)*raptat)
                jespacem=int( float(m2max)*raptat)
              end if
!
              if(raptat.gt.0.) then
!
                contig= m1mx-m1mi.le.iespacem .and. &
                        m2mx-m2mi.le.jespacem
              else
                contig=.true.
              end if
            else
              contig=.false.
            end if
!
            ncelgr=ncelgr+1
            igmx=min(ig+igr,i2m1)
            jgmx=min(jg+jgr,j2m1)
            kgmx=min(kg+kgr,k2m1)
!
            if(contig) then
!
!             la maille grossiere est completement rattachee a un element
!             reduit de paroi qu'il faut explorer
!
              ncelat=ncelat+1
!
              do k=kg,kgmx
               do j=jg,jgmx
                ndeb=ind(0   ,j,k)
                do i=ig,igmx
                  n=ndeb+i
                  if(mnpar(n).eq.0) then
!                   calcul de la distance mini
!
                    mc=n-n0
                    dmini=1.e+20
!                   boucle coupee en 2 pour vectorisation
                    do m1=m1mi,m1mx
                     do m2=m2mi,m2mx
                       mfac0=m1+(m2-1)*idm
                       dist2(mfacn)=(xcc(mc)-xpar(mfacn))**2+ &
                                   (ycc(mc)-ypar(mfacn))**2+ &
                                   (zcc(mc)-zpar(mfacn))**2
                     end do
                    end do
                    do m1=m1mi,m1mx
                     do m2=m2mi,m2mx
                       mfac0=m1+(m2-1)*idm
                       mfacn=m0n+mfac0
                       if(dist2(mfacn).lt.dmini) then
                         imini=mfacn
                         dmini=dist2(mfacn)
                       end if
                     end do
                    end do
!
                    dist(n) =sqrt(dmini)
                    mnpar(n)=imini
!
                  end if
!                 fin des boucles d'exploration du maillage fin dans
!                 une maille du maillage grossier
                enddo
               enddo
              enddo
!
            else
!
!             la maille grossiere n'est pas completement rattachee a un
!             element bien defini de paroi. On revient a l'exploration
!             complete des parois.
!
              do k=kg,kgmx
               do j=jg,jgmx
                do i=ig,igmx
                  n=ind(i,j,k)
                  if(mnpar(n).eq.0) then
!                   calcul de la distance mini
!                   Boucle sur les facettes des parois. Recherche du
!                   minimum de la distance paroi par paroi.
!                      nfbi : numero interne paroi
!                      mpn  : pointeur fin frontiere precedente norm. stockees
!                      mmb  : nombre de facettes sur une paroi
!                      n    : pointeur cellule tab toutes cellules
!                     imin  : pointeur facette paroi la plus proche dans
!                             vecteur des facettes de toutes les parois
!
                    mc=n-n0
                    dmini=1.e+20
                    mpar0=npbrat(l)+1
                    mpar1=npbrat(l)+nbdrat(l)
                    do mp=mpar0,mpar1
                      nfbi=lbdrat(mp)
                      m0b=mpn(nfbi)
                      mbmx=mmb(nfbi)
!
                      do mb=1,mbmx
                        mbb=m0b+mb
                        dist2(mb)=(xcc(mc)-xpar(mbb))**2+ &
                                  (ycc(mc)-ypar(mbb))**2+ &
                                  (zcc(mc)-zpar(mbb))**2
                      enddo
!
!                     existe une fonction pour le minimum
                      do mb=1,mbmx
                        if(dist2(mb).lt.dmini) then
                          imini=m0b+mb
                          dmini=dist2(mb)
                        end if
                      enddo
                    enddo
                    dist(n) =sqrt(dmini)
                    mnpar(n)=imini
!
                  end if
!
                  enddo
                enddo
              enddo
!             fin teste sur rattachement d'une cellule du maillage grossier
            end if
!
            enddo
           enddo
          enddo
!
      return
      end
