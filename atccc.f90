          subroutine atccc( &
                 x,y,z, &
                 xcc,ycc,zcc, &
                 l)
!
!***********************************************************************
!
!     ACT
!_A    calcul des coordonnees des centres des cellules
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
!_O    xcc        : arg real(ip00      ) ; coordonnee x centre facette paroi
!_O    ycc        : arg real(ip00      ) ; coordonnee y centre facette paroi
!_O    zcc        : arg real(ip00      ) ; coordonnee z centre facette paroi
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
	  use maillage
!
!-----------------------------------------------------------------------
!
      dimension x(ip21),y(ip21),z(ip21)
      dimension xcc(ip00),ycc(ip00),zcc(ip00)
!
      ind(  j,k)=n0+1+(j-jd1(l))*nid+(k-kd1(l))*nijd
!
      n0=npc(l)
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
!     calcul des centres des facettes pour les plans k=constante (facteur 4)
      do k=k1,k2
        do j=j1,j2-1
          indjk=ind(j,k)
          do i=i1,i2-1
            n=indjk+(i-id1(l))
            m=n-n0
            xcc(m)=x(n)+x(n+1)+x(n+nid)+x(n+1+nid)
            ycc(m)=y(n)+y(n+1)+y(n+nid)+y(n+1+nid)
            zcc(m)=z(n)+z(n+1)+z(n+nid)+z(n+1+nid)
          enddo
         enddo
        enddo
!
!     moyenne sur k pour coordonnees au centre des cellules
      do k=k1,k2-1
        do j=j1,j2-1
          indjk=ind(j,k)
          do i=i1,i2-1
            n=indjk+(i-id1(l))
            m=n-n0
            xcc(m)=(xcc(m)+xcc(m+nijd))*0.125
            ycc(m)=(ycc(m)+ycc(m+nijd))*0.125
            zcc(m)=(zcc(m)+zcc(m+nijd))*0.125
           enddo
         enddo
        enddo
!
      return
      end
