module tools

  interface reallocate
     module procedure reallocate_1r,reallocate_1i, &
          reallocate_2r,reallocate_2i, &
          reallocate_3r,reallocate_3i, &
          reallocate_4r,reallocate_4i, &
          reallocate_1c
  end interface reallocate
  interface reallocate_s
     module procedure reallocate_s_1r,reallocate_s_4i,reallocate_s_1i,reallocate_1c_s,reallocate_s_2i
  end interface reallocate_s

contains

subroutine reallocate_1r(in,newsize)
 implicit none
 double precision,allocatable,intent(inout) :: in(:)
 integer,intent(in)             :: newsize

! if(newsize/=size(in)) then
   if(allocated(in)) deallocate(in)
   allocate(in(newsize))
   in=0.
! endif

end subroutine reallocate_1r

subroutine reallocate_2r(in,newsize1,newsize2)
 implicit none
 double precision,allocatable,intent(inout) :: in(:,:)
 integer,intent(in)             :: newsize1,newsize2

! if(newsize1*newsize2/=size(in)) then
   if(allocated(in)) deallocate(in)
   allocate(in(newsize1,newsize2))
   in=0.
! endif

end subroutine reallocate_2r

subroutine reallocate_3r(in,newsize1,newsize2,newsize3)
 implicit none
 double precision,allocatable,intent(inout) :: in(:,:,:)
 integer,intent(in)             :: newsize1,newsize2,newsize3

! if(newsize1*newsize2*newsize3/=size(in)) then
   if(allocated(in)) deallocate(in)
   allocate(in(newsize1,newsize2,newsize3))
   in=0.
! endif

end subroutine reallocate_3r

subroutine reallocate_4r(in,newsize1,newsize2,newsize3,newsize4)
 implicit none
 double precision,allocatable,intent(inout) :: in(:,:,:,:)
 integer,intent(in)                :: newsize1,newsize2,newsize3,newsize4

! if(newsize1*newsize2*newsize3*newsize4/=size(in)) then
   if(allocated(in)) deallocate(in)
   allocate(in(newsize1,newsize2,newsize3,newsize4))
  in=0.
! endif

end subroutine reallocate_4r

subroutine reallocate_1i(in,newsize)
 implicit none
 integer,allocatable,intent(inout) :: in(:)
 integer,intent(in)                :: newsize

! if(newsize/=size(in)) then
   if(allocated(in)) deallocate(in)
   allocate(in(newsize))
   in=0
! endif

end subroutine reallocate_1i

subroutine reallocate_2i(in,newsize1,newsize2)
 implicit none
 integer,allocatable,intent(inout) :: in(:,:)
 integer,intent(in)                :: newsize1,newsize2

! if(newsize1*newsize2/=size(in)) then
   if(allocated(in)) deallocate(in)
   allocate(in(newsize1,newsize2))
   in=0
! endif

end subroutine reallocate_2i

subroutine reallocate_3i(in,newsize1,newsize2,newsize3)
 implicit none
 integer,allocatable,intent(inout) :: in(:,:,:)
 integer,intent(in)                :: newsize1,newsize2,newsize3

! if(newsize1*newsize2*newsize3/=size(in)) then
   if(allocated(in)) deallocate(in)
   allocate(in(newsize1,newsize2,newsize3))
   in=0
! endif

end subroutine reallocate_3i

subroutine reallocate_4i(in,newsize1,newsize2,newsize3,newsize4)
 implicit none
 integer,allocatable,intent(inout) :: in(:,:,:,:)
 integer,intent(in)                :: newsize1,newsize2,newsize3,newsize4

! if(newsize1*newsize2*newsize3*newsize4/=size(in)) then
   if(allocated(in)) deallocate(in)
   allocate(in(newsize1,newsize2,newsize3,newsize4))
   in=0
! endif

end subroutine reallocate_4i

subroutine reallocate_1c(in,newsize)
 implicit none
 character(*),allocatable,intent(inout) :: in(:)
 integer,intent(in)                :: newsize

! if(newsize/=size(in)) then
   if(allocated(in)) deallocate(in)
   allocate(in(newsize))
   in=""
! endif

end subroutine reallocate_1c

subroutine reallocate_1c_s(tab,newsize)
  implicit none
  character(len=2),allocatable::tab(:),tab1(:)
  integer :: newsize
  allocate(tab1(size(tab)))
  tab1=tab
  deallocate(tab)
  allocate(tab(newsize))
  tab(:size(tab1))=tab1
end subroutine reallocate_1c_s

subroutine reallocate_s_1i(in,newsize)
 implicit none
 integer,allocatable::in(:),tmp(:)
 integer :: newsize
 if(newsize>size(in)) then
   allocate(tmp(size(in)))
   tmp=in
   call reallocate(in,newsize)
   in=0
   in(:size(tmp))=tmp
   deallocate(tmp)
 elseif(newsize<size(in)) then
  write(*,*) "Probleme dans les tailles : ",newsize,size(in)
!  call backtrace
  call abort
 endif
end subroutine reallocate_s_1i

subroutine reallocate_s_2i(in,newsize1,newsize2)
 implicit none
 integer,allocatable::in(:,:),tmp(:,:)
 integer :: newsize1,newsize2
 if(newsize1*newsize2>size(in)) then
   allocate(tmp(size(in,1),size(in,2)))
   tmp=in
   call reallocate(in,newsize1,newsize2)
   in=0
   in(:size(tmp,1),:size(tmp,2))=tmp
   deallocate(tmp)
 elseif(newsize1*newsize2<size(in)) then
  write(*,*) "Probleme dans les tailles : ",newsize1*newsize2,size(in)
!  call backtrace
  call abort
 endif
end subroutine reallocate_s_2i

subroutine reallocate_s_4i(in,newsize1,newsize2,newsize3,newsize4)
 implicit none
 integer,allocatable::in(:,:,:,:),tmp(:,:,:,:)
 integer :: newsize1,newsize2,newsize3,newsize4
 if(newsize1*newsize2*newsize3*newsize4>size(in)) then
   allocate(tmp(size(in,1),size(in,2),size(in,3),size(in,4)))
   tmp=in
   call reallocate(in,newsize1,newsize2,newsize3,newsize4)
   in=0
   in(:size(tmp,1),:size(tmp,2),:size(tmp,3),:size(tmp,4))=tmp
   deallocate(tmp)
 elseif(newsize1*newsize2*newsize3*newsize4<size(in)) then
  write(*,*) "Probleme dans les tailles : ",newsize1*newsize2*newsize3*newsize4,size(in)
!  call backtrace
  call abort
 endif
end subroutine reallocate_s_4i

subroutine reallocate_s_1r(in,newsize)
 implicit none
 double precision,allocatable::in(:),tmp(:)
 integer :: newsize
 if(newsize>size(in)) then
   allocate(tmp(size(in)))
   tmp=in
   call reallocate(in,newsize)
   in=0.
   in(:size(tmp))=tmp
   deallocate(tmp)
 elseif(newsize<size(in)) then
  write(*,*) "Probleme dans les tailles : ",newsize,size(in)
!  call backtrace
  call abort
 endif
end subroutine reallocate_s_1r
end module tools
