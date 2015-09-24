#define WITH_MPI

module mod_mpi
#if defined(WITH_MPI)
  USE MPI
#endif
#if defined(__INTEL_COMPILER)
  use ifport,only : fseek,ftell,getpid
#endif
  implicit none
  integer             :: rank,NPROCS
  integer             :: num_bcg=0,num_bci=0,num_bcl=0
  integer             :: num_bg=0,num_bi=0,num_bl=0
  integer,allocatable :: bcg_to_proc(:),bcg_to_bcl(:),bcl_to_bcg(:),bcg_to_bci(:)
  integer,allocatable :: bg_to_proc(:),bg_to_bl(:),bl_to_bg(:),bg_to_bi(:)
  !ab_to_ab where a is :
  ! bc : boundary condition
  ! b  : block
  !and b is :
  ! i : initial numerotation
  ! g : global numerotation
  ! l : local numerotation
  !ab can be proc to have the rank of the process that own the object

  INTERFACE SUM_MPI
    !      SUM_MPI(A,B)
    ! COMPUTE THE SUM OF A IN B
    MODULE PROCEDURE SUM_MPI_0R,SUM_MPI_0I,&
          SUM_MPI_1R,SUM_MPI_1I
  END INTERFACE SUM_MPI

  INTERFACE MAX_MPI
    !      MAX_MPI(A,B)
    ! COMPUTE THE MAX OF A IN B
    MODULE PROCEDURE MAX_MPI_0I
  END INTERFACE MAX_MPI

  INTERFACE BCAST
    !        BCAST(A,ORIG)
    ! BROADCAST THE MESSAGE A FROM ORIG FOR EVERY PROC
    ! RETURN WHEN EVERYTHING IS DONE
    MODULE PROCEDURE BCAST_0R,BCAST_0I,&
          BCAST_1R,BCAST_1I,&
          BCAST_2R,BCAST_2I
  END INTERFACE BCAST

  INTERFACE GATHER
    ! THIS ROUTINE IS GATHERING DATA FOR EVERY PROC FROM EVERY PROC
    ! EXEMPLE
    ! PROC 1 : IN=[A,B,C,D] , SIZE=3 , OUT=[A,B,C,E,F]
    ! PROC 2 : IN=[E,F]     , SIZE=2 , OUT=[A,B,C,E,F]
    ! PROC 3 : IN=[]        , SIZE=0 , OUT=[A,B,C,E,F]
    ! MPI_ALLGATHERV(IN, SIZE, MPI_TYPE(IN),OUT, SIZE(OUT), SHIFT, MPI_TYPE(OUT), MPI_COMM_WORLD,IERR)
    MODULE PROCEDURE GATHER_R,GATHER_I
  END INTERFACE GATHER

  INTERFACE MPI_TRANS
    !        SEND A MESSAGE WITH MPI
    MODULE PROCEDURE MPI_TRANS_R1,MPI_TRANS_R2,MPI_TRANS_I,MPI_TRANS_I0,MPI_TRANS_C0
  END INTERFACE MPI_TRANS

  INTERFACE MPI_ITRANS2
    !        SEND A MESSAGE WITH MPI
    MODULE PROCEDURE MPI_ITRANS2_R1,MPI_ITRANS2_R2!,MPI_ITRANS2_I,MPI_ITRANS2_I0,MPI_ITRANS2_C0
  END INTERFACE MPI_ITRANS2

contains

  !************************************
  SUBROUTINE  INIMPI
      !************************************
      IMPLICIT NONE
      integer :: ierr

#if defined(WITH_MPI)
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
#else
      RANK=0
      NPROCS=1
#endif

      allocate(bcg_to_proc(0),bcg_to_bcl(0),bcl_to_bcg(0),bcg_to_bci(0))
      allocate(bg_to_proc(0),bg_to_bl(0),bl_to_bg(0),bg_to_bi(0))
  END SUBROUTINE  INIMPI

  subroutine endmpi
      implicit none
      integer :: ierr
#if defined(WITH_MPI)
      CALL MPI_FINALIZE(IERR)
#endif

  end subroutine endmpi




  SUBROUTINE MPI_TRANS_R1(A,B,ORIG,DEST)
      !ORIG SEND THE MESSAGE A TO DEST
      !DEST RECV THE MESSAGE B FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      double precision   ,INTENT(INOUT) :: A(:)
      double precision   ,INTENT(IN)    :: B(:)
      integer,INTENT(IN)    :: ORIG,DEST
#ifdef WITH_MPI
      integer :: STATUS(MPI_STATUS_SIZE),TAG
      integer :: ierr
#endif

      IF (RANK==ORIG.AND.RANK==DEST) THEN ! SENDING A MESSAGE TO MYSELF
        A=B
#ifdef WITH_MPI
      ELSE
        TAG=ORIG*NPROCS+DEST

        IF(RANK==ORIG)  & ! I'M ORIG, I SEND THE MESSAGE B TO DEST
            CALL MPI_SEND(B(1),SIZE(B),MPI_REAL8,DEST, &
            TAG,MPI_COMM_WORLD,IERR)

        IF(RANK==DEST) &  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
            CALL MPI_RECV(A(1),SIZE(A),MPI_REAL8,ORIG, &
            TAG,MPI_COMM_WORLD,STATUS,IERR)
#endif
      ENDIF

  END SUBROUTINE MPI_TRANS_R1

  SUBROUTINE MPI_TRANS_R2(A,B,ORIG,DEST)
      !ORIG SEND THE MESSAGE A TO DEST
      !DEST RECV THE MESSAGE B FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      double precision   ,INTENT(INOUT) :: A(:,:)
      double precision   ,INTENT(IN)    :: B(:,:)
      integer,INTENT(IN)    :: ORIG,DEST
#ifdef WITH_MPI
      integer :: STATUS(MPI_STATUS_SIZE),TAG
      integer :: ierr
#endif

      IF (RANK==ORIG.AND.RANK==DEST) THEN ! SENDING A MESSAGE TO MYSELF
        A=B
#ifdef WITH_MPI
      ELSE
        TAG=ORIG*NPROCS+DEST

        IF(RANK==ORIG)  & ! I'M ORIG, I SEND THE MESSAGE B TO DEST
            CALL MPI_SEND(B(1,1),SIZE(B),MPI_REAL8,DEST, &
            TAG,MPI_COMM_WORLD,IERR)

        IF(RANK==DEST) &  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
            CALL MPI_RECV(A(1,1),SIZE(A),MPI_REAL8,ORIG, &
            TAG,MPI_COMM_WORLD,STATUS,IERR)
#endif
      ENDIF

  END SUBROUTINE MPI_TRANS_R2

  SUBROUTINE MPI_TRANS_I(A,B,ORIG,DEST)
      !ORIG SEND THE MESSAGE A TO DEST
      !DEST RECV THE MESSAGE B FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      integer,INTENT(INOUT) :: A(:)
      integer,INTENT(IN)    :: B(:)
      integer,INTENT(IN)    :: ORIG,DEST
#ifdef WITH_MPI
      integer :: STATUS(MPI_STATUS_SIZE),TAG
      integer :: ierr
#endif

      IF (RANK==ORIG.AND.RANK==DEST) THEN ! SENDING A MESSAGE TO MYSELF
        A=B
#ifdef WITH_MPI
      ELSE
        TAG=ORIG*NPROCS+DEST

        IF(RANK==ORIG)  & ! I'M ORIG, I SEND THE MESSAGE B TO DEST
            CALL MPI_SEND(B(1),SIZE(B),MPI_INTEGER,DEST, &
            TAG,MPI_COMM_WORLD,IERR)

        IF(RANK==DEST) &  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
            CALL MPI_RECV(A(1),SIZE(A),MPI_INTEGER,ORIG, &
            TAG,MPI_COMM_WORLD,STATUS,IERR)
#endif
      ENDIF

  END SUBROUTINE MPI_TRANS_I


  SUBROUTINE MPI_TRANS_I0(A,B,ORIG,DEST)
      !ORIG SEND THE MESSAGE A TO DEST
      !DEST RECV THE MESSAGE B FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      integer,INTENT(INOUT) :: A
      integer,INTENT(IN)    :: B
      integer,INTENT(IN)    :: ORIG,DEST
#ifdef WITH_MPI
      integer :: STATUS(MPI_STATUS_SIZE),TAG
      integer :: ierr
#endif

      IF (RANK==ORIG.AND.RANK==DEST) THEN ! SENDING A MESSAGE TO MYSELF
        A=B
#ifdef WITH_MPI
      ELSE
        TAG=ORIG*NPROCS+DEST

        IF(RANK==ORIG)  & ! I'M ORIG, I SEND THE MESSAGE B TO DEST
            CALL MPI_SEND(B,1,MPI_INTEGER,DEST, &
            TAG,MPI_COMM_WORLD,IERR)

        IF(RANK==DEST) &  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
            CALL MPI_RECV(A,1,MPI_INTEGER,ORIG, &
            TAG,MPI_COMM_WORLD,STATUS,IERR)
#endif
      ENDIF

  END SUBROUTINE MPI_TRANS_I0

  SUBROUTINE MPI_TRANS_C0(A,B,ORIG,DEST)
      !ORIG SEND THE MESSAGE A TO DEST
      !DEST RECV THE MESSAGE B FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      character(* ),INTENT(INOUT) :: A
      character(* ),INTENT(IN)    :: B
      integer,INTENT(IN)    :: ORIG,DEST
#ifdef WITH_MPI
      integer :: STATUS(MPI_STATUS_SIZE),TAG
      integer :: ierr
#endif

      IF (RANK==ORIG.AND.RANK==DEST) THEN ! SENDING A MESSAGE TO MYSELF
        A=B
#ifdef WITH_MPI
      ELSE
        TAG=ORIG*NPROCS+DEST

        IF(RANK==ORIG)  & ! I'M ORIG, I SEND THE MESSAGE B TO DEST
            CALL MPI_SEND(B,len(b),MPI_CHARACTER,DEST, &
            TAG,MPI_COMM_WORLD,IERR)

        IF(RANK==DEST) &  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
            CALL MPI_RECV(A,len(a),MPI_CHARACTER,ORIG, &
            TAG,MPI_COMM_WORLD,STATUS,IERR)
#endif
      ENDIF

  END SUBROUTINE MPI_TRANS_C0

  SUBROUTINE MPI_ITRANS_BEGIN(A,ORIG,DEST,REQS,REQR)
      !PREPARE THE REQUESTS ARRAY
      IMPLICIT NONE
      double precision   ,INTENT(INOUT) :: A(:)
      integer,INTENT(IN)    :: ORIG,DEST
      integer,INTENT(INOUT) :: REQS,REQR
      integer :: ierr
      integer :: TAG

#ifdef WITH_MPI
      TAG = (ORIG+1)*NPROCS*2+DEST+1

      IF(DEST/=ORIG) THEN
        CALL MPI_SEND_INIT(A(1), SIZE(A), MPI_REAL8, DEST, &
            TAG, MPI_COMM_WORLD,REQS,IERR)! SEND THE BUFFER
        CALL MPI_RECV_INIT(A(1), SIZE(A), MPI_REAL8, ORIG, &
            TAG, MPI_COMM_WORLD,REQR,IERR)! RECV THE BUFFER
      ENDIF
#endif

  END SUBROUTINE MPI_ITRANS_BEGIN

  SUBROUTINE MPI_ITRANS(REQ)
      !ORIG SEND THE MESSAGE A TO DEST
      !DEST RECV THE MESSAGE B FROM ORIG
      !RETURN WHEN EVERYTHING IS STARTED
      !(THE MESSAGES HAVEN'T BEEN DELIVERED YET)
      IMPLICIT NONE
      integer,INTENT(INOUT) :: REQ
      integer :: ierr

#ifdef WITH_MPI
      IF(REQ/=MPI_REQUEST_NULL)  CALL MPI_START(REQ,IERR)
#endif

  END SUBROUTINE MPI_ITRANS

  SUBROUTINE MPI_ITRANS2_R1(A,ORIG,DEST,REQ)
      !ORIG SEND THE MESSAGE A TO DEST
      !DEST RECV THE MESSAGE B FROM ORIG
      !RETURN WHEN EVERYTHING IS STARTED
      !(THE MESSAGES HAVEN'T BEEN DELIVERED YET)
      IMPLICIT NONE
      double precision   ,INTENT(INOUT) :: A(:)
      integer,INTENT(IN)    :: ORIG,DEST
      integer,INTENT(INOUT) :: REQ
      integer :: ierr
      integer :: TAG

#ifdef WITH_MPI
      TAG = (ORIG+1)*NPROCS*2+DEST+1
      if (rank==ORIG) &
          CALL MPI_ISEND(A(1), SIZE(A), MPI_REAL8, DEST, &
          TAG, MPI_COMM_WORLD,REQ,IERR)! SEND THE BUFFER

      if (rank==dest) &
          CALL MPI_IRECV(A(1), SIZE(A), MPI_REAL8, ORIG, &
          TAG, MPI_COMM_WORLD,REQ,IERR)! RECV THE BUFFER
#endif

  END SUBROUTINE MPI_ITRANS2_R1


  SUBROUTINE MPI_ITRANS2_R2(A,ORIG,DEST,REQ)
      !ORIG SEND THE MESSAGE A TO DEST
      !DEST RECV THE MESSAGE B FROM ORIG
      !RETURN WHEN EVERYTHING IS STARTED
      !(THE MESSAGES HAVEN'T BEEN DELIVERED YET)
      IMPLICIT NONE
      double precision   ,INTENT(INOUT) :: A(:,:)
      integer,INTENT(IN)    :: ORIG,DEST
      integer,INTENT(INOUT) :: REQ
      integer :: ierr
      integer :: TAG

#ifdef WITH_MPI
      TAG = (ORIG+1)*NPROCS*2+DEST+1
      if (rank==ORIG) &
          CALL MPI_ISEND(A(1,1), SIZE(A), MPI_REAL8, DEST, &
          TAG, MPI_COMM_WORLD,REQ,IERR)! SEND THE BUFFER

      if (rank==dest) &
          CALL MPI_IRECV(A(1,1), SIZE(A), MPI_REAL8, ORIG, &
          TAG, MPI_COMM_WORLD,REQ,IERR)! RECV THE BUFFER
#endif

  END SUBROUTINE MPI_ITRANS2_R2

  SUBROUTINE WAIT_MPI(REQ)
      !WAIT FOR A COMMUNICATION TO END
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      integer :: REQ
#ifdef WITH_MPI
      integer :: STATUS(MPI_STATUS_SIZE)
      integer :: ierr

      IF(REQ/=MPI_REQUEST_NULL)   CALL MPI_WAIT(REQ, STATUS, IERR)
#endif

  END SUBROUTINE WAIT_MPI

  SUBROUTINE BARRIER()
      ! SYNC POINT FOR ALL PROCESS
      use sortiefichier
      IMPLICIT NONE

      interface
        !function fsync (fd) bind(c,name="fsync")
            !use iso_c_binding, only: c_int
            !integer(c_int), value :: fd
            !integer(c_int) :: fsync
        !end function fsync
      end interface

      integer :: ierr,i
      !do i =1,10000
#ifdef WITH_MPI
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR )
#endif

        flush(imp)
        !ierr = fsync(fnum(imp)) ! not supported on ifort

#ifdef WITH_MPI
        CALL MPI_Barrier(MPI_COMM_WORLD,IERR )
#endif
      !enddo
  END SUBROUTINE BARRIER

  SUBROUTINE GATHER_P(IN,OUT)
      ! THIS ROUTINE IS GATHERING INTEGER FOR EVERY PROC FROM EVERY PROC
      ! EXEMPLE
      ! PROC 1 : IN=A , OUT=[A , B , C]
      ! PROC 2 : IN=B , OUT=[A , B , C]
      ! PROC 3 : IN=C , OUT=[A , B , C]
      IMPLICIT NONE
      integer,INTENT(IN)  :: IN
      integer,INTENT(OUT) :: OUT(:)
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_ALLGATHER(IN, 1, MPI_INTEGER, OUT(1), 1, MPI_INTEGER, MPI_COMM_WORLD,IERR)
#else
      OUT(:)=IN
#endif

  END SUBROUTINE GATHER_P

  SUBROUTINE MAXLOC_MPI(IN,OUT)
      IMPLICIT NONE
      ! THIS ROUTINE IS SEARCHING FOR THE MAX OF IN(1,J)
      ! IN(2,J) IS USED TO IDENTIFY THE OWNER OF THE MAX VALUE
      ! EXEMPLE
      ! PROC 1 : IN=[A,1,B,1] , OUT=[A,1,B,2]
      ! PROC 2 : IN=[C,2,D,2] , OUT=[A,1,B,2]

      double precision,INTENT(IN)    ::  IN(:,:)
      double precision,INTENT(INOUT) :: OUT(:,:)
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_ALLREDUCE(IN(1,1), OUT(1,1), SIZE(IN,2), MPI_2DOUBLE_PRECISION,MPI_MAXLOC, MPI_COMM_WORLD,IERR)
#endif

  END SUBROUTINE MAXLOC_MPI

  SUBROUTINE LOR_MPI(IN,OUT)
      IMPLICIT NONE
      ! THIS ROUTINE COMPUTE AN "OR" OPERATION BETWEEN IN

      LOGICAL,INTENT(IN)  ::  IN
      LOGICAL,INTENT(OUT) :: OUT
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_ALLREDUCE(IN, OUT, 1, MPI_LOGICAL,MPI_LOR, MPI_COMM_WORLD,IERR)
#else
      OUT=IN
#endif

  END SUBROUTINE LOR_MPI

  SUBROUTINE FREE_MPI_REQ(REQ)
      IMPLICIT NONE
      ! THIS ROUTINE COMPUTE AN "OR" OPERATION BETWEEN IN

      integer,INTENT(INOUT)  ::  REQ
      integer :: ierr
#ifdef WITH_MPI
      if(req/=MPI_REQUEST_NULL) CALL MPI_REQUEST_FREE(REQ,IERR)
#endif

  END SUBROUTINE FREE_MPI_REQ



  SUBROUTINE REALLOCATE_1R(IN,SIZE)
      IMPLICIT NONE
      double precision,ALLOCATABLE,INTENT(INOUT) :: IN(:)
      integer,INTENT(IN)             :: SIZE

      IF(ALLOCATED(IN)) DEALLOCATE(IN)
      ALLOCATE(IN(SIZE))

  END SUBROUTINE REALLOCATE_1R

  SUBROUTINE REALLOCATE_2R(IN,SIZE1,SIZE2)
      IMPLICIT NONE
      double precision,ALLOCATABLE,INTENT(INOUT) :: IN(:,:)
      integer,INTENT(IN)             :: SIZE1,SIZE2

      IF(ALLOCATED(IN)) DEALLOCATE(IN)
      ALLOCATE(IN(SIZE1,SIZE2))

  END SUBROUTINE REALLOCATE_2R

  SUBROUTINE REALLOCATE_3R(IN,SIZE1,SIZE2,SIZE3)
      IMPLICIT NONE
      double precision,ALLOCATABLE,INTENT(INOUT) :: IN(:,:,:)
      integer,INTENT(IN)             :: SIZE1,SIZE2,SIZE3

      IF(ALLOCATED(IN)) DEALLOCATE(IN)
      ALLOCATE(IN(SIZE1,SIZE2,SIZE3))

  END SUBROUTINE REALLOCATE_3R


  SUBROUTINE REALLOCATE_1I(IN,SIZE)
      IMPLICIT NONE
      integer,ALLOCATABLE,INTENT(INOUT) :: IN(:)
      integer,INTENT(IN)                :: SIZE

      IF(ALLOCATED(IN)) DEALLOCATE(IN)
      ALLOCATE(IN(SIZE))

  END SUBROUTINE REALLOCATE_1I

  SUBROUTINE REALLOCATE_2I(IN,SIZE1,SIZE2)
      IMPLICIT NONE
      integer,ALLOCATABLE,INTENT(INOUT) :: IN(:,:)
      integer,INTENT(IN)                :: SIZE1,SIZE2

      IF(ALLOCATED(IN)) DEALLOCATE(IN)
      ALLOCATE(IN(SIZE1,SIZE2))

  END SUBROUTINE REALLOCATE_2I

  SUBROUTINE REALLOCATE_3I(IN,SIZE1,SIZE2,SIZE3)
      IMPLICIT NONE
      integer,ALLOCATABLE,INTENT(INOUT) :: IN(:,:,:)
      integer,INTENT(IN)                :: SIZE1,SIZE2,SIZE3

      IF(ALLOCATED(IN)) DEALLOCATE(IN)
      ALLOCATE(IN(SIZE1,SIZE2,SIZE3))

  END SUBROUTINE REALLOCATE_3I

  SUBROUTINE SUM_MPI_0R(A,B)
      !COMPUTE THE SUM OF A IN B
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      double precision   ,INTENT(INOUT)    :: A
      double precision   ,INTENT(OUT),optional :: B
#ifdef WITH_MPI
      double precision ::C
      integer :: ierr
      CALL MPI_ALLREDUCE(A, C, 1, MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,IERR)
      if (present(B)) then
        B=C
      else
        A=C
      endif
#else
      if (present(B))  B=A
#endif

  END SUBROUTINE SUM_MPI_0R

  SUBROUTINE SUM_MPI_1R(A,B)
      !COMPUTE THE SUM OF A IN B
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      double precision   ,INTENT(INOUT)    :: A(:)
      double precision   ,INTENT(OUT),optional :: B(:)
#ifdef WITH_MPI
      double precision,allocatable ::C(:)
      integer :: ierr
      allocate(c(size(A)))
      CALL MPI_ALLREDUCE(A(1), C(1), SIZE(A), MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,IERR)
      if (present(B)) then
        B=C
      else
        A=C
      endif
#else
      if (present(B))  B=A
#endif

  END SUBROUTINE SUM_MPI_1R


  SUBROUTINE SUM_MPI_0I(A,B)
      !COMPUTE THE SUM OF A IN B
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      integer   ,INTENT(INOUT)    :: A
      integer   ,INTENT(OUT),optional :: B
#ifdef WITH_MPI
      integer ::C
      integer :: ierr
      CALL MPI_ALLREDUCE(A, C, 1, MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD,IERR)
      if (present(B)) then
        B=C
      else
        A=C
      endif
#else
      if (present(B))  B=A
#endif

  END SUBROUTINE SUM_MPI_0I

  SUBROUTINE MAX_MPI_0I(A,B)
      !COMPUTE THE SUM OF A IN B
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      integer   ,INTENT(INOUT)    :: A
      integer   ,INTENT(OUT),optional :: B
#ifdef WITH_MPI
      integer ::C
      integer :: ierr
      CALL MPI_ALLREDUCE(A, C, 1, MPI_INTEGER,MPI_MAX, MPI_COMM_WORLD,IERR)
      if (present(B)) then
        B=C
      else
        A=C
      endif
#else
      if (present(B))  B=A
#endif
  END SUBROUTINE MAX_MPI_0I

  SUBROUTINE SUM_MPI_1I(A,B)
      !COMPUTE THE SUM OF A IN B
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      integer   ,INTENT(INOUT)    :: A(:)
      integer   ,INTENT(OUT),optional :: B(:)
#ifdef WITH_MPI
      integer,allocatable ::C(:)
      integer :: ierr
      allocate(c(size(a)))
      CALL MPI_ALLREDUCE(A(1), C(1), SIZE(A), MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD,IERR)
      if (present(B)) then
        B=C
      else
        A=C
      endif
#else
      if (present(B))  B=A
#endif
  END SUBROUTINE SUM_MPI_1I

  SUBROUTINE BCAST_0I(IN,ORIG)
      !BROADCAST THE MESSAGE IN FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      integer   ,INTENT(INOUT) :: IN
      integer,INTENT(IN)    :: ORIG
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_BCAST( IN,1, MPI_INTEGER, ORIG, MPI_COMM_WORLD,IERR)
#endif

  END SUBROUTINE BCAST_0I

  SUBROUTINE BCAST_0R(IN,ORIG)
      !BROADCAST THE MESSAGE IN FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      double precision   ,INTENT(INOUT) :: IN
      integer,INTENT(IN)    :: ORIG
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_BCAST( IN,1, MPI_REAL8, ORIG, MPI_COMM_WORLD,IERR)
#endif

  END SUBROUTINE BCAST_0R

  SUBROUTINE BCAST_1R(IN,ORIG)
      !BROADCAST THE MESSAGE IN FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      double precision   ,INTENT(INOUT) :: IN(:)
      integer,INTENT(IN)    :: ORIG
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_BCAST( IN(1),SIZE(IN), MPI_REAL8, ORIG, MPI_COMM_WORLD,IERR)
#endif

  END SUBROUTINE BCAST_1R

  SUBROUTINE BCAST_2R(IN,ORIG)
      !BROADCAST THE MESSAGE IN FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      double precision   ,INTENT(INOUT) :: IN(:,:)
      integer,INTENT(IN)    :: ORIG
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_BCAST( IN(1,1),SIZE(IN), MPI_REAL8, ORIG, MPI_COMM_WORLD,IERR)
#endif

  END SUBROUTINE BCAST_2R

  SUBROUTINE BCAST_1I(IN,ORIG)
      !BROADCAST THE MESSAGE IN FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      integer,INTENT(INOUT) :: IN(:)
      integer,INTENT(IN)    :: ORIG
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_BCAST( IN(1),SIZE(IN), MPI_INTEGER, ORIG, MPI_COMM_WORLD,IERR)
#endif

  END SUBROUTINE BCAST_1I

  SUBROUTINE BCAST_2I(IN,ORIG)
      !BROADCAST THE MESSAGE IN FROM ORIG
      !RETURN WHEN EVERYTHING IS DONE
      IMPLICIT NONE
      integer,INTENT(INOUT) :: IN(:,:)
      integer,INTENT(IN)    :: ORIG
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_BCAST( IN(1,1),SIZE(IN), MPI_INTEGER, ORIG, MPI_COMM_WORLD,IERR)
#endif

  END SUBROUTINE BCAST_2I

  SUBROUTINE GATHER_I(IN,OUT,SIZE)
      ! THIS ROUTINE IS GATHERING INTEGER FOR EVERY PROC FROM EVERY PROC
      ! EXEMPLE
      ! PROC 1 : IN=[A,B,C,D] , SIZE=3 , OUT=[A,B,C,E,F]
      ! PROC 2 : IN=[E,F]     , SIZE=2 , OUT=[A,B,C,E,F]
      ! PROC 3 : IN=[]        , SIZE=0 , OUT=[A,B,C,E,F]
      IMPLICIT NONE
      integer,INTENT(IN)  :: IN(:),SIZE
      integer,INTENT(OUT) :: OUT(:)
      integer             :: I,DISP(NPROCS),SIZE_ALL(NPROCS)
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_ALLGATHER(SIZE, 1, MPI_INTEGER, &
          SIZE_ALL, 1, MPI_INTEGER, MPI_COMM_WORLD,IERR)
      DISP(1)=0
      DO I=2,NPROCS
      DISP(I)=DISP(I-1) + SIZE_ALL(I-1)
      ENDDO

      CALL MPI_ALLGATHERV(IN(1), SIZE, MPI_INTEGER, &
          OUT(1), SIZE_ALL, DISP, MPI_INTEGER, MPI_COMM_WORLD,IERR)
#else
      OUT=IN
#endif

  END SUBROUTINE GATHER_I

  SUBROUTINE GATHER_R(IN,OUT,SIZE)
      IMPLICIT NONE
      ! THIS ROUTINE IS GATHERING REAL FOR EVERY PROC FROM EVERY PROC
      ! EXEMPLE
      ! PROC 1 : IN=[A,B,C,D] , SIZE=3 , OUT=[A,B,C,E,F]
      ! PROC 2 : IN=[E,F]     , SIZE=2 , OUT=[A,B,C,E,F]
      ! PROC 3 : IN=[]        , SIZE=0 , OUT=[A,B,C,E,F]

      double precision,INTENT(IN)     :: IN(:)
      integer,INTENT(IN)  :: SIZE
      double precision,INTENT(OUT)    :: OUT(:)
      integer             :: I,DISP(NPROCS),SIZE_ALL(NPROCS)
      integer :: ierr
#ifdef WITH_MPI
      CALL MPI_ALLGATHER(SIZE, 1, MPI_INTEGER, &
          SIZE_ALL, 1, MPI_INTEGER, MPI_COMM_WORLD,IERR)
      DO I=1,NPROCS
      DISP(I)=SUM(SIZE_ALL(1:I-1))
      ENDDO

      CALL MPI_ALLGATHERV(IN(1), SIZE, MPI_REAL8, &
          OUT(1), SIZE_ALL, DISP, MPI_REAL8, MPI_COMM_WORLD,IERR)
#else
      OUT=IN
#endif
  END SUBROUTINE GATHER_R


  SUBROUTINE START_KEEP_ORDER(relais_in)
      IMPLICIT NONE
      integer,intent(inout),optional :: relais_in
      integer :: relais
      relais=0
      if (present(relais_in)) relais=relais_in
#ifdef WITH_MPI
      if(rank>0) call mpi_trans(relais,relais,rank-1,rank)
#endif
      if (present(relais_in)) relais_in=relais
  END SUBROUTINE START_KEEP_ORDER


  SUBROUTINE END_KEEP_ORDER(relais_in)
      IMPLICIT NONE
      integer,intent(inout),optional :: relais_in
      integer :: relais
      if (present(relais_in)) relais=relais_in
#ifdef WITH_MPI
      if(rank<nprocs-1) call mpi_trans(relais,relais,rank,rank+1)
#endif
      if (present(relais_in)) relais_in=relais
  END SUBROUTINE END_KEEP_ORDER

  subroutine my_fseek(unit,pos)
      implicit none
      integer :: unit,pos
#if defined(__INTEL_COMPILER)
      integer :: i
      i=fseek(unit,pos,0)
#else
      call fseek(unit,pos,0)
#endif
  end subroutine my_fseek
end module mod_mpi

