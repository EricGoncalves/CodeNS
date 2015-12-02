module mod_atdist_2
  implicit none
contains
  subroutine atdist_2(img, &
       x,y,z, &
       xpar,ypar,zpar, &
       dist2,dist,mnpar)
!
!***********************************************************************
!
!     ACT
!_A    calcul de la facette de paroi la plus proche d'une cellule
!_A    et rattachement de la cellule au pointeur de la facette dans
!_A    tableaux des frontieres a normales stockees
!_A    VERSION POUR LE PARALLELE
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
    use mod_mpi
    use tools
    use mod_atccc
    use mod_at_ecrdist
    implicit none
    integer          ::           i,         i1,         i2,          j
    integer          ::          j1,         j2,          k,         k1,         k2
    integer          ::           l,        m0b,         mb,        mbb,       mbmx
    integer          ::          mc,mnpar(ip12),         n0,         nc
    integer          ::         nid,       nijd,        njd,num,numt
    integer          ::         bcg,        bcl,        img,         lm
    logical          ::     isparoi
    double precision ::  dist(ip12),dist2(ip00),    x(ip21),  xcc(ip00)
    double precision ::  xpar(ip00),    y(ip21),  ycc(ip00), ypar(ip00),    z(ip21)
    double precision ::   zcc(ip00), zpar(ip00)
    double precision,allocatable :: buff(:,:)
!
!-----------------------------------------------------------------------
!
!

   numt=0
   do bcl=1,mtb
   if ((cl(bcl)(1:2).eq.'pa').or. &
       (cl(bcl)(1:2).eq.'lp').or. &
       (cl(bcl)(1:2).eq.'gl'))      numt=numt+mmb(bcl)
   enddo
   call sum_mpi(numt)
   num=0
   dist=Huge(1.)
!  boucle sur toutes les conditions limites
   do bcg=1,num_bcg

!     test si paroi
      isparoi=.False.
      if(rank==bcg_to_proc(bcg)) then
        bcl=bcg_to_bcl(bcg)
        isparoi=((cl(bcl)(1:2).eq.'pa').or.(cl(bcl)(1:2).eq.'lp').or. &
              (cl(bcl)(1:2).eq.'gl'))
      endif
      call lor_mpi(isparoi)
      if(isparoi) then

!       broadcast the coordinates of points on the boundary
        mbmx=0
        if(rank==bcg_to_proc(bcg)) mbmx=mmb(bcl)
        call bcast(mbmx,bcg_to_proc(bcg))
        num=num+mbmx
        if (rank==0) write(*,*) "computing dist with the boundary",bcg," : ",100*num/numt,"%"
        call reallocate(buff,3,mbmx)
        if(rank==bcg_to_proc(bcg)) then
           m0b=mpn(bcl)
           do mb=1,mbmx
              mbb=m0b+mb
              buff(1,mb)=xpar(mbb)
              buff(2,mb)=ypar(mbb)
              buff(3,mb)=zpar(mbb)
           enddo
        endif
        call bcast(buff,bcg_to_proc(bcg))

!       boucle sur les domaines locaux
       do l=1,lzx
          lm=l+(img-1)*lz

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

!         calcul des centres des cellules
          call atccc( &
               x,y,z, &
               xcc,ycc,zcc, &
               lm)

!         boucle sur les cellules du domaine
          do k=k1,k2-1
             do j=j1,j2-1
                do i=i1,i2-1
                   nc=ind(i,j,k)
                   mc=nc-n0

!                  boucle sur les points de la condition limite
                   do mb=1,mbmx
!                    calcul des distances
                     dist2(mc)=(xcc(mc)-buff(1,mb))**2+(ycc(mc)-buff(2,mb))**2+ &
                          (zcc(mc)-buff(3,mb))**2
                     if(dist2(mc).lt.dist(nc)) then
                        mnpar(nc)=mb
!                        mnpar2(nc)=bcg
                        dist(nc)=dist2(mc)
                     end if
                   enddo
                enddo
             enddo
          enddo
       enddo
     endif
   enddo

  dist=sqrt(dist)

    if(kecrdis.eq.1) then
      do l=1,lzx
        lm=l+(img-1)*lz
!
!           ecriture disque des distances (fichiers separes "fdist_l")
         call at_ecrdist( &
              lm,         &
              dist,mnpar)
      enddo
    endif
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine atdist_2
end module mod_atdist_2
