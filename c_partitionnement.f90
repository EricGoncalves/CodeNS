module mod_c_partitionnement
  implicit none
contains
  subroutine c_partitionnement(mot,imot,nmot,x,y,z,ncbd,mnc,ncin,bceqt,exs1,exs2)
    use mod_partitionnement
    use mod_tcmd_partitionnement
    use para_fige,only : nmx
    use mod_mpi,only : num_bg
    implicit none
    double precision,allocatable,intent(inout) :: bceqt(:,:)
    integer         ,allocatable,intent(inout) :: ncin(:),mnc(:),ncbd(:)
    double precision,allocatable,intent(inout) :: x(:),y(:),z(:)
    double precision            ,intent(inout) :: exs1,exs2
    integer,intent(in)          ::  imot(nmx),nmot
    character(len=32),intent(in) ::  mot(nmx)
    integer          :: nblocks
    integer          :: nsplit(num_bg),nsplit_dir(3,num_bg)

    nblocks=1 ! initial values
    nsplit=1
    nsplit_dir=0 ! 0 mean automatic

    call tcmd_partitionnement(mot,imot,nmot,nblocks,nsplit,nsplit_dir)

    call partitionnement(x,y,z,ncbd,mnc,ncin,bceqt,exs1,exs2,nblocks,nsplit,nsplit_dir)

  end subroutine c_partitionnement
end module mod_c_partitionnement
