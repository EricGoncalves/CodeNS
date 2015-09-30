module mod_crdms
  implicit none
contains
  subroutine crdms( &
       l,ni,nj,nk)
!
!***********************************************************************
!
!     ACT
!_A    Creation d'un domaine structure essentiellement par la
!_A    determination de son encombrement en nombre de points.
!
!     INP
!_I    ni         : arg int              ; nbr de plans indices i
!_I    nj         : arg int              ; nbr de plans indices j
!_I    nk         : arg int              ; nbr de plans indices k
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    nfi        : com int              ; nbr de rangees de pts fictifs
!
!     OUT
!_O    id1        : com int (lz,lg     ) ; indice min en i fictif
!_O    ii1        : com int (lz,lg     ) ; indice min en i reel
!_O    ii2        : com int (lz,lg     ) ; indice max en i reel
!_O    id2        : com int (lz,lg     ) ; indice max en i fictif
!_O    jd1        : com int (lz,lg     ) ; indice min en j fictif
!_O    jj1        : com int (lz,lg     ) ; indice min en j reel
!_O    jj2        : com int (lz,lg     ) ; indice max en j reel
!_O    jd2        : com int (lz,lg     ) ; indice max en j fictif
!_O    kd1        : com int (lz,lg     ) ; indice min en k fictif
!_O    kk1        : com int (lz,lg     ) ; indice min en k reel
!_O    kk2        : com int (lz,lg     ) ; indice max en k reel
!_O    kd2        : com int (lz,lg     ) ; indice max en k fictif
!
!     I/O
!_/    lzx        : com int              ; nbr total de domaines
!_/    ndimubx    : com int              ; nbr de cellules du plus grd domaine
!_/                                        (pts fictifs inclus)
!_/    ndimctbx   : com int              ; nbr de cellules de tts les domaines
!_/                                        (pts fictifs inclus)
!_/    ndimntbx   : com int              ; nbr de noeuds de tts les domaines
!_/                                        (pts fictifs inclus)
!_/    npn        : com int (lz,lg     ) ; pointeur fin de dom precedent
!_/                                        dans tab tous noeuds
!_/    nnn        : com int (lz,lg     ) ; nombre de noeuds du dom (dont fic.)
!_/    npc        : com int (lz,lg     ) ; pointeur fin de dom precedent
!_/                                        dans tab toutes cellules
!_/    nnc        : com int (lz,lg     ) ; nombre de cellules du dom (dont fic.)
!_/    npfb       : com int (lz,lg     ) ; pointeur fin de dom precedent
!_/                                        dans tab toutes facettes
!_/    nnfb       : com int (lz,lg     ) ; nombre de facettes du dom (dont fic.)
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use maillage
    use kcle
    use chainecarac
    use schemanum
    use tools
    use mod_mpi
    implicit none
    integer          ::   img, imgi, imgj, imgk,    l
    integer          ::    lm,   ni,  nid,   nj,  njd
    integer          ::    nk,  nkd,nptfs

    num_bg=num_bg+1
    num_bi=num_bi+1
    call reallocate_s(bg_to_proc,num_bg)
    call reallocate_s(bg_to_bl,num_bg)
!    call reallocate_s(bg_to_bi,num_bg)
    bg_to_proc(l)=modulo((l-1),nprocs)
    bg_to_bl(l)=0
!    bg_to_bi(l)=l
!
    if(bg_to_proc(l)==rank) then
      num_bl=num_bl+1
  !
  !-----------------------------------------------------------------------
  !
      lzx=lzx+1
      klzx=2

  !    lz=lzx
      lt=lgx*lzx
      call reallocate_s(ii1,lt)
      call reallocate_s(jj1,lt)
      call reallocate_s(kk1,lt)
      call reallocate_s(ii2,lt)
      call reallocate_s(jj2,lt)
      call reallocate_s(kk2,lt)
      call reallocate_s(id1,lt)
      call reallocate_s(jd1,lt)
      call reallocate_s(kd1,lt)
      call reallocate_s(id2,lt)
      call reallocate_s(jd2,lt)
      call reallocate_s(kd2,lt)
      call reallocate_s(nnn,lt)
      call reallocate_s(nnc,lt)
      call reallocate_s(nnfb,lt)
      call reallocate_s(npn,lt)
      call reallocate_s(npc,lt)
      call reallocate_s(npfb,lt)
      call reallocate_s(bl_to_bg,lt)

      bl_to_bg(lzx)=l
      bg_to_bl(l)=lzx
  !
      do img=1,lgx
  !
         lm=lzx+(img-1)*lz
  !
         imgi = img
         imgj = img
         imgk = img
         if (equat(3:5).eq.'2di') imgi = 1
         if (equat(3:5).eq.'2dj') imgj = 1
         if (equat(3:5).eq.'2dk') imgk = 1
  !
         ii1(lm)=1
         jj1(lm)=1
         kk1(lm)=1
  !
         ii2(lm)= (ni-1)/2**(imgi-1)+1
         jj2(lm)= (nj-1)/2**(imgj-1)+1
         kk2(lm)= (nk-1)/2**(imgk-1)+1
  !
         id1(lm)=ii1(lm)-nfi
         jd1(lm)=jj1(lm)-nfi
         kd1(lm)=kk1(lm)-nfi
         id2(lm)=ii2(lm)+nfi
         jd2(lm)=jj2(lm)+nfi
         kd2(lm)=kk2(lm)+nfi
  !
         nid=id2(lm)-id1(lm)+1
         njd=jd2(lm)-jd1(lm)+1
         nkd=kd2(lm)-kd1(lm)+1
         nptfs =nid*njd*nkd
         nnn(lm)=nptfs
!         nptfs =(nid-1)*(njd-1)*(nkd-1)
         nnc(lm)=nptfs
         nnfb(lm)=nind*nptfs
  !
         npn(lm)=ndimntbx
         npc(lm)=ndimctbx
         npfb(lm)=nind*ndimntbx
  !
         ndimubx =max(ndimubx,nnn(lm))
!         ndimubx =max(ndimubx,nnc(lm))
         ndimctbx=ndimctbx+nnc(lm)
         ndimntbx=ndimntbx+nnn(lm)
      enddo
    endif
!
    return
  end subroutine crdms
end module mod_crdms
