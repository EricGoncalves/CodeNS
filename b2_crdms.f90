module mod_b2_crdms
  implicit none
contains
  subroutine b2_crdms(l)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'crdm'.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    nnn        : arg int              ; nombre de noeuds du dom (dont fic.)
!_I    nnc        : arg int              ; nombre de cellules du dom (dont fic.)
!_I    nnfb       : arg int              ; nombre de facettes du dom (dont fic.)
!_I    id1        : arg int              ; indice min en i fictif
!_I    ii1        : arg int              ; indice min en i reel
!_I    ii2        : arg int              ; indice max en i reel
!_I    id2        : arg int              ; indice max en i fictif
!_I    jd1        : arg int              ; indice min en j fictif
!_I    jj1        : arg int              ; indice min en j reel
!_I    jj2        : arg int              ; indice max en j reel
!_I    jd2        : arg int              ; indice max en j fictif
!_I    kd1        : arg int              ; indice min en k fictif
!_I    kk1        : arg int              ; indice min en k reel
!_I    kk2        : arg int              ; indice max en k reel
!_I    kd2        : arg int              ; indice max en k fictif
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use maillage
    implicit none
  integer          :: img,  l, lm
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
!
    form='(/ 2x,''numero de la grille      : '',11x,i5/' &
         //'2x,''nb de noeuds (dont fic.) : '',10x,i6/' &
         //'2x,''tot. noeuds dom preced   : '', 9x,i7/' &
         //'2x,''nb de cellules(dont fic.): '',10x,i6/' &
         //'2x,''tot. cellules dom preced : '', 9x,i7/' &
         //'2x,''nb de facettes(dont fic.): '', 9x,i7/' &
         //'2x,''tot. facettes dom preced : '', 9x,i7/' &
         //'2x,''domaine d''''indicage       : '',/' &
         //'2x,5h id1=,i5,6h  ii1=,i5,6h  ii2=,i5,6h  id2=,i5,/' &
         //'2x,5h jd1=,i5,6h  jj1=,i5,6h  jj2=,i5,6h  jd2=,i5,/' &
         //'2x,5h kd1=,i5,6h  kk1=,i5,6h  kk2=,i5,6h  kd2=,i5)'
!
    do img=1,lgx
       lm=l+(img-1)*lz
       write(imp,form) img, &
            nnn(lm),npn(lm),nnc(lm),npc(lm),nnfb(lm),npfb(lm), &
            id1(lm),ii1(lm),ii2(lm),id2(lm), &
            jd1(lm),jj1(lm),jj2(lm),jd2(lm), &
            kd1(lm),kk1(lm),kk2(lm),kd2(lm)
    enddo
!
    return
  end subroutine b2_crdms
end module mod_b2_crdms
