module mod_ingr
  implicit none
contains
  subroutine ingr(l,x,y,z,king)
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des coordonnees des noeuds du maillage,
!_A    c'est-a-dire des tableaux x, y et z.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    king       : arg int              ; cle initialisation maillage
!_I    td         : com char(lz        ) ; type de domaine (struct./non struct.)
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    ptrans     : com real             ; distance pour periodicite
!_I    protat     : com real             ; angle(rad) pour periodicite
!_I    kdgv       : com int              ; unite logiq, coordonnees des noeuds
!
!     OUT
!_O    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_O    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_O    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!     I/O
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage 
    use sortiefichier
    use definition
    use mod_readdg
    use mod_smg_fcm
    implicit none
    integer          ::  img,king,   l
    double precision :: x,y,z
!
!-----------------------------------------------------------------------
!
    dimension x(ip21),y(ip21),z(ip21)
!
    if(king.ne.0) then
!     remplissage des tableaux x , y , z dans un sous-programme d'initialisation
!        call utinig(l,x,y,z)
    elseif(king.eq.0) then
!     remplissage des tableaux x , y , z par lecture d' un fichier de reprise
       call readdg(l,kdgv,x,y,z)
    endif
!
    if(lgx.gt.1) then
       do img = 2,lgx
          call smg_fcm( &
               img-1,img, &
               l, &
               x,y,z)
       enddo
    end if
!
    return
  end subroutine ingr
end module mod_ingr
