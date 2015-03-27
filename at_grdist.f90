module mod_at_grdist
  implicit none
contains
  subroutine at_grdist( &
       igr,jgr,kgr, &
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
!_A    Distance des cellules aux noeuds d'un maillage grossier par
!_A    calcul de la facette de paroi la plus proche de la cellule
!_A    et rattachement de la cellule au pointeur de la facette dans
!_A    tableaux des frontieres a normales stockees.
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
!_I    igr        : arg int              ; pas en i pour definition maillage
!_I                                        grossier
!_I    jgr        : arg int              ; pas en j pour definition maillage
!_I                                        grossier
!_I    kgr        : arg int              ; pas en k pour definition maillage
!_I                                        grossier
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
    implicit none
    integer          ::           i,         i1,         i2,         ib,        igr
    integer          ::        igr0,      imini,          j,         j1,         j2
    integer          ::          jb,        jgr,       jgr0,          k,         k1
    integer          ::          k2,         kb,        kgr,       kgr0,          l
    integer          ::         m0b, m1tb(ip00), m2tb(ip00),         mb,        mbb
    integer          ::        mbmx,         mc,mnpar(ip12),         mp,      mpar0
    integer          ::       mpar1,         n0,         nc,       nfbi,nfrtb(ip00)
    integer          ::         nid,       nijd,        njd
    double precision ::  dist(ip12),dist2(ip00),      dmini,    x(ip21),  xcc(ip00)
    double precision ::  xpar(ip00),    y(ip21),  ycc(ip00), ypar(ip00),    z(ip21)
    double precision ::   zcc(ip00), zpar(ip00)
!
!-----------------------------------------------------------------------
!
!

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
!     boucle sur les cellules maillage grossier, premieres a dernieres lignes
!
    kgr0=kgr
    jgr0=jgr
    igr0=igr
    if(mod(k2-1,kgr).eq.1 .or. kgr.eq.1) kgr0=0
    if(mod(j2-1,jgr).eq.1 .or. jgr.eq.1) jgr0=0
    if(mod(i2-1,igr).eq.1 .or. igr.eq.1) igr0=0

    do kb=k1-kgr0,k2-1,kgr
       k=kb
       if(kb.le.0) k=k2-1
       do jb=j1-jgr0,j2-1,jgr
          j=jb
          if(jb.le.0) j=j2-1
          do ib=i1-igr0,i2-1,igr
             i=ib
             if(ib.le.0) i=i2-1
!
             nc=ind(i,j,k)
             mc=nc-n0
!
!           Boucle sur les facettes des parois. Recherche du minimum de la
!           distance paroi par paroi.
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
!           boucle sur les parois de la liste
!
                nfbi=lbdrat(mp)
                m0b=mpn(nfbi)
                mbmx=mmb(nfbi)
!
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
!
             enddo
             dist(nc) =sqrt(dmini)
             mnpar(nc)=imini
!
!         boucle cellule suivante
          enddo
       enddo
    enddo
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine at_grdist
end module mod_at_grdist
