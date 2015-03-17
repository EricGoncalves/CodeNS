module mod_atdist_3
  implicit none
contains
  subroutine atdist_3( &
       x,y,z, &
       xpar,ypar,zpar, &
       xcc,ycc,zcc,dist2, &
       dist,mnpar, &
       l)
!
!***********************************************************************
!
!     ACT
!_A    calcul de la facette de paroi la plus proche d'une cellule
!_A    et rattachement de la cellule au pointeur de la facette dans
!_A    tableaux des frontieres a normales stockees.
!-A    Pour chaque domaine, une liste de paroi doit etre donnee.
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
    implicit none
    integer          ::     i,   i1,   i2,imini,  ind
    integer          ::     j,   j1,   j2,    k,   k1
    integer          ::    k2,    l,  m0b,   mb,  mbb
    integer          ::  mbmx,   mc,mnpar,   mp,mpar0
    integer          :: mpar1,   n0,   nc, nfbi,  nid
    integer          ::  nijd,  njd
    double precision ::  dist,dist2,dmini,    x,  xcc
    double precision ::  xpar,    y,  ycc, ypar,    z
    double precision ::   zcc, zpar
!
!-----------------------------------------------------------------------
!
    dimension x(ip21),y(ip21),z(ip21)
    dimension xpar(ip00),ypar(ip00),zpar(ip00)
    dimension xcc (ip00),ycc (ip00),zcc (ip00),dist2(ip00)
    dimension dist(ip12),mnpar(ip12)
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
!     boucle sur les cellules
    do k=k1,k2-1
       do j=j1,j2-1
          do i=i1,i2-1
             nc=ind(i,j,k)
             mc=nc-n0
!
!           Boucle sur les facettes des parois. Recherche du minimum de la
!           distance paroi par paroi de la liste.
!                    nfbi : numero interne paroi
!                    mpn  : pointeur fin frontiere precedente norm. stockees
!                    mmb  : nombre de facettes sur une paroi
!                    nc   : pointeur cellule tab toutes cellules
!                   imin  : pointeur facette paroi la plus proche dans
!                           vecteur des facettes de toutes les parois
!
             dmini=1.e+20
             mpar0=npbrat(l)+1
             mpar1=npbrat(l)+nbdrat(l)
             do mp=mpar0,mpar1
                nfbi=lbdrat(mp)
                m0b=mpn(nfbi)
                mbmx=mmb(nfbi)
                do mb=1,mbmx
                   mbb=m0b+mb
                   dist2(mb)=(xcc(mc)-xpar(mbb))**2+(ycc(mc)-ypar(mbb))**2+ &
                        (zcc(mc)-zpar(mbb))**2
                enddo
!
!             existe une fonction pour le minimum
                do mb=1,mbmx
                   if(dist2(mb).lt.dmini) then
                      imini=m0b+mb
                      dmini=dist2(mb)
                   end if
                enddo
             enddo
!
             dist(nc) =sqrt(dmini)
             mnpar(nc)=imini
          enddo
       enddo
    enddo
!
    return
  end subroutine atdist_3
end module mod_atdist_3
