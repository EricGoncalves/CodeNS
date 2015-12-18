module mod_atintrans2
  implicit none
contains
  subroutine atintrans2(x,y,z,ncin,fgam)
!
!***********************************************************************
!
!     ACT
!_A   Lecture des donnees pour fixer la transition
!_A   initialisation du tableau de la fonction d'intermittence "fgam"
!
!_A   - Recherche du mot cle TRANSITION en colonne 1
!_A   - lecture du nombre de zones laminaires
!_A   - lecture des indices extremes dans les 3 directions pour chaque zone
!_A     les indices sont les numeros des noeuds du maillage
!
!     INP
!_I    mdimtnx    : com int              ; nbr de pts de ttes les front
!_I                                        a normales stockees
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!_I    kd2        : com int (lt        ) ; indice max en k fictif
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    iminb      : com int (mtt       ) ; indice min en i d'une front
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une front
!_I    jminb      : com int (mtt       ) ; indice min en j d'une front
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une front
!_I    kminb      : com int (mtt       ) ; indice min en k d'une front
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une front
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!
!     LOC
!_L    l          : arg int              ; numero de domaine
!_L    ilmin       : arg int              ; indice min en i
!_L    ilmax       : arg int              ; indice max en i
!_L    jlmin       : arg int              ; indice min en j
!_L    jlmax       : arg int              ; indice max en j
!_L    klmin       : arg int              ; indice min en k
!_L    klmax       : arg int              ; indice max en k
!
!     OUT
!_O    fgam       : arg real(ip42      ) ; fonction d'intermittence pour
!_O                                        transition
!_O    ktransi    : com int              ; cle transition 1->transition fixee
!_O                                                       0->pas de transition
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use sortiefichier
    use modeleturb
    use mod_mpi
    implicit none
    integer          ::          l,       lig,       mfe,       mfl,         n
    integer          :: ncin(ip41),       nfr,     nfrmx
    integer          ::        nid,      nijd,       njd,mfbe,m,mb,nc,i,j,k,xyz,mn
    double precision ::    x(ip21), y(ip21), z(ip21)
    double precision :: fgam(ip12),xmin,xmax,ymin,ymax,zmin,zmax
    double precision :: xcc,ycc,zcc
    double precision,parameter :: eps=1e-10
!
!-----------------------------------------------------------------------

    character(len=80) ::  ligne
!
    open(99,file='fatdon',status='old',err=100)
!
!     recherche de la ligne "TRANSITION"
!
    lig=0
    do
       lig=lig+1
       read(99,300,err=304,end=301)ligne
       if(ligne(1:6).EQ.'TRANSI') exit
    enddo
!
!     transition fixee
    ktransi=1
!
!     initialisation de "fgam" a 1
    do n=1,mdimtnx
       fgam(n)=1.
    enddo
!
!     lecture des frontieres laminaires
!
    lig=lig+1
    read(99,*,err=302)nfrmx
    if(nfrmx.le.0) then
!       le nombre de zones laminaires doit etre positif
        if (rank==0) write(imp,'(/,"!!!atintrans2: domaine laminaire mal defini",14x,"nfrmx=",i5)')nfrmx
       stop
    endif
!
    do nfr=1,nfrmx
      ! lecture du numero de la frontiere concernée
      ! et des coordonées de deux points extremes définissant une zone
      ! Tout point de cette frontiere inclus dans cette zone
       read(99,*,err=303)mfe,xmin,ymin,zmin,xmax,ymax,zmax
!
!       calcul numero facette de la frontiere
!
       do mfbe=1,mtb
       if(bcg_to_bci(bcl_to_bcg(mfbe))==mfe) then
         mfl=nfei(mfbe)
         l =ndlb(mfl)
  !
         nid = id2(l)-id1(l)+1
         njd = jd2(l)-jd1(l)+1
         nijd = nid*njd
  !
       do m=1,mmb(mfl)
          mb=mpb(mfl)+m
          mn=mpn(mfl)+m
          nc=ncin(mb)
          
          !calcul de la coordonnée du centre de la cellule
          xcc=0.
          ycc=0.
          zcc=0.
          do k=0,1
            do j=0,1
              do i=0,1
                xyz =nc + i + j*nid + k*nijd
                xcc=xcc+x(xyz)  
                ycc=ycc+y(xyz)
                zcc=zcc+z(xyz)
              enddo
            enddo
          enddo
          xcc=xcc*0.125
          ycc=ycc*0.125
          zcc=zcc*0.125
          
          if((abs(abs(xcc-xmin)+abs(xcc-xmax)-abs(xmin-xmax))<=eps).and. &
             (abs(abs(ycc-ymin)+abs(ycc-ymax)-abs(ymin-ymax))<=eps).and. &
             (abs(abs(zcc-zmin)+abs(zcc-zmax)-abs(zmin-zmax))<=eps)) then
             fgam(mn)=0.
          endif
                  
        enddo
        endif
!     fin boucle sur frontieres laminaires
      enddo
    enddo
    close(99)
!
300 format(a80)
!
    return
!
301 continue
!     fin de fichier. Pas de transition
    ktransi=0
    if (rank==0) write(imp,'(/,"===>atintrans2: ktransi=0. Calcul tout laminaire ou tout turbulent")')
    close(99)
    return
!
100 continue
    if (rank==0) write(imp,'(/,"!!!atintrans2: erreur ouverture fichier fatdon")')
    stop
!
302 continue
    if (rank==0) write(imp,'(/,"!!!atintrans2: erreur nombre de frontieres laminaires ")')
    stop
!
303 continue
    if (rank==0) write(imp,'(/,"!!!atintrans2: erreur lecture frontiere ligne",i4)')lig
    stop
!
304 continue
    if (rank==0) write(imp,'(/,"!!!atintrans2: erreur positionnement dans fatdon")')
    stop
!
    return
  end subroutine atintrans2
end module mod_atintrans2
