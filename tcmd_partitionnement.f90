module mod_tcmd_partitionnement
  implicit none
contains
  subroutine tcmd_partitionnement(mot,imot,nmot,nblocks,nsplit,nsplit_dir)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action partitionnement.
!
!-----parameters figes--------------------------------------------------
!
    use mod_mpi,only : num_bg
    use chainecarac,only : ci
    use para_fige,only : nmx
    use mod_valenti
    implicit none
    integer          ::      icmt,imot(nmx),       nm,     nmot,kval,nblocks,lgi
    integer          :: nsplit(num_bg),nsplit_dir(3,num_bg),dir
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
!
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
    kval=0
!
    nm=1
    do while(nm.lt.nmot)
      nm=nm+1
      select case(mot(nm)(1:imot(nm)))
      case('nblocks')! read number of minimum block the user want
        nm=nm+1
        if(nmot.lt.nm) then
           comment=ci
           call synterr(mot,imot,nmot,comment)
        else
           call valenti(mot,imot,nm,nblocks,kval)
        endif
      case('nsplit')! read number of number of sub_block the user want in a block
        nm=nm+1
        if(nmot.lt.nm) then
           comment=ci
           call synterr(mot,imot,nmot,comment)
        else
           call valenti(mot,imot,nm,lgi,kval)
        endif
        nm=nm+1
        if(nmot.lt.nm) then ! firt the block number
           comment=ci
           call synterr(mot,imot,nmot,comment)
        else
           call valenti(mot,imot,nm,nsplit(lgi),kval)
        endif
      case('force') ! read the splitting the user want in a block
        nm=nm+1
        if(nmot.lt.nm) then ! firt the block number
           comment=ci
           call synterr(mot,imot,nmot,comment)
        else
           call valenti(mot,imot,nm,lgi,kval)
        endif
        do dir=1,3
          nm=nm+1
          if(nmot.lt.nm) then ! then the number of block per direction
             comment=ci
             call synterr(mot,imot,nmot,comment)
          else
             call valenti(mot,imot,nm,nsplit_dir(dir,lgi),kval)
          endif
        enddo
      case default
        comment=ci
        call synterr(mot,imot,nmot,comment)
      end select
    enddo

!
    return
  end subroutine tcmd_partitionnement
end module mod_tcmd_partitionnement
