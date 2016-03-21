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
  integer,allocatable :: bcg_to_proc(:),bcg_to_bcl(:),bcl_to_bcg(:)
  integer,allocatable :: bg_to_proc(:),bg_to_bl(:),bl_to_bg(:),bcg_to_bg(:)
  integer,allocatable :: bcint_to_bcintg(:),bcintg_to_proc(:)
  integer,allocatable :: lbdko_to_lbdkog(:),lbdkog_to_proc(:)
  !ab_to_ab where a is :
  ! bc : boundary condition
  ! b  : block
  !and b is :
  ! i : initial numerotation
  ! g : global numerotation
  ! l : local numerotation
  !ab can be proc to have the rank of the process that own the object

#ifndef WITH_MPI
  integer,parameter :: mpi_request_null = 0
  integer,parameter :: mpi_comm_world = 0
#endif
#if defined(__PGI)
  integer,external :: fseek,ftell,getpid
#endif

  INTERFACE SUM_MPI
    !       SUM_MPI(A,B)
    ! COMPUTE THE SUM OF A IN B
    !       SUM_MPI(A)
    ! COMPUTE THE SUM OF A IN A
    ! RETURN WHEN EVERYTHING IS DONE
    MODULE PROCEDURE SUM_MPI_0R,SUM_MPI_0I,&
          SUM_MPI_1R,SUM_MPI_1I
  END INTERFACE SUM_MPI

  INTERFACE MAX_MPI
    !       MAX_MPI(A,B)
    ! COMPUTE THE MAX OF A IN B
    !       MAX_MPI(A)
    ! COMPUTE THE MAX OF A IN A
    ! RETURN WHEN EVERYTHING IS DONE
    MODULE PROCEDURE MAX_MPI_0I,MAX_MPI_1R
  END INTERFACE MAX_MPI

  INTERFACE MIN_MPI
    !       MIN_MPI(A,B)
    ! COMPUTE THE MAX OF A IN B
    !       MIN_MPI(A)
    ! COMPUTE THE MAX OF A IN A
    ! RETURN WHEN EVERYTHING IS DONE
    MODULE PROCEDURE MIN_MPI_0R
  END INTERFACE MIN_MPI

  INTERFACE BCAST
    !       BCAST(A,ORIG)
    ! BROADCAST THE MESSAGE A FROM ORIG FOR EVERY PROC
    ! RETURN WHEN EVERYTHING IS DONE
    MODULE PROCEDURE BCAST_0R,BCAST_0I,&
          BCAST_1R,BCAST_1I,BCAST_0C,&
          BCAST_2R,BCAST_2I,BCAST_3I
  END INTERFACE BCAST

  INTERFACE GATHER
    ! THIS ROUTINE IS GATHERING DATA FOR EVERY PROC FROM EVERY PROC
    ! IF IN IS A VECTOR
    !       GATHER(IN,OUT,SIZE)
    ! EXEMPLE
    ! PROC 1 : IN=[A,B,C,D] , SIZE=3 , OUT=[A,B,C,E,F]
    ! PROC 2 : IN=[E,F]     , SIZE=2 , OUT=[A,B,C,E,F]
    ! PROC 3 : IN=[]        , SIZE=0 , OUT=[A,B,C,E,F]
    ! IF IN IS A SCALAR
    !       GATHER(IN,OUT)
    ! EXEMPLE
    ! PROC 1 : IN=A , OUT=[A,B,C]
    ! PROC 2 : IN=B , OUT=[A,B,C]
    ! PROC 3 : IN=C , OUT=[A,B,C]
    MODULE PROCEDURE GATHER_R,GATHER_I,GATHER_I0,GATHER_C
    END INTERFACE GATHER

    INTERFACE MPI_TRANS
      !       MPI_TRANS(A,B,ORIG,DEST)
      ! SEND A MESSAGE WITH MPI AND WAIT FOR IT TO BE RECIEVED
      MODULE PROCEDURE MPI_TRANS_R1,MPI_TRANS_R2,MPI_TRANS_R4,&
            MPI_TRANS_I1,MPI_TRANS_I2,MPI_TRANS_I0,&
            MPI_TRANS_C0
    END INTERFACE MPI_TRANS

    INTERFACE MPI_ITRANS2
      !       MPI_TRANS(A,B,ORIG,DEST,REQ)
      ! SEND A MESSAGE WITH MPI AND DOES NOT WAIT FOR IT TO BE RECIEVED
      MODULE PROCEDURE MPI_ITRANS2_R1,MPI_ITRANS2_R2
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

        allocate(bcg_to_proc(1),bcg_to_bcl(1),bcl_to_bcg(1))
        allocate(bg_to_proc(1),bg_to_bl(1),bl_to_bg(1),bcg_to_bg(1))

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
        double precision,allocatable :: buff(:) ! security necessary
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
          STATUS=MPI_STATUS_IGNORE
          IF(RANK==ORIG)  then ! I'M ORIG, I SEND THE MESSAGE B TO DEST
              allocate(buff(size(b)))
              buff=b
              CALL MPI_SEND(Buff(1),SIZE(B),MPI_REAL8,DEST, &
              TAG,MPI_COMM_WORLD,IERR)
              deallocate(buff)
          endif
          IF(RANK==DEST) then  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
              allocate(buff(size(a)))
              CALL MPI_RECV(buff(1),SIZE(A),MPI_REAL8,ORIG, &
              TAG,MPI_COMM_WORLD,STATUS,IERR)
              a=buff
              deallocate(buff)
          endif
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
        double precision,allocatable :: buff(:,:) ! security necessary
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
          STATUS=MPI_STATUS_IGNORE
          IF(RANK==ORIG)  then ! I'M ORIG, I SEND THE MESSAGE B TO DEST
              allocate(buff(size(b,1),size(b,2)))
              buff=b
              CALL MPI_SEND(Buff(1,1),SIZE(B),MPI_REAL8,DEST, &
              TAG,MPI_COMM_WORLD,IERR)
              deallocate(buff)
          endif
          IF(RANK==DEST) then  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
              allocate(buff(size(a,1),size(a,2)))
              CALL MPI_RECV(buff(1,1),SIZE(A),MPI_REAL8,ORIG, &
              TAG,MPI_COMM_WORLD,STATUS,IERR)
              a=buff
              deallocate(buff)
          endif
#endif
        ENDIF

    END SUBROUTINE MPI_TRANS_R2


    SUBROUTINE MPI_TRANS_R4(A,B,ORIG,DEST)
        !ORIG SEND THE MESSAGE B TO DEST
        !DEST RECV THE MESSAGE A FROM ORIG
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        double precision   ,INTENT(INOUT) :: A(:,:,:,:)
        double precision   ,INTENT(IN)    :: B(:,:,:,:)
        double precision,allocatable :: buff(:,:,:,:) ! security necessary
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
          STATUS=MPI_STATUS_IGNORE
          IF(RANK==ORIG)  then ! I'M ORIG, I SEND THE MESSAGE B TO DEST
              allocate(buff(size(b,1),size(b,2),size(b,3),size(b,4)))
              buff=b
              CALL MPI_SEND(Buff(1,1,1,1),SIZE(B),MPI_REAL8,DEST, &
              TAG,MPI_COMM_WORLD,IERR)
              deallocate(buff)
          endif
          IF(RANK==DEST) then  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
              allocate(buff(size(a,1),size(a,2),size(a,3),size(a,4)))
              CALL MPI_RECV(buff(1,1,1,1),SIZE(A),MPI_REAL8,ORIG, &
              TAG,MPI_COMM_WORLD,STATUS,IERR)
              a=buff
              deallocate(buff)
          endif
#endif
        ENDIF

    END SUBROUTINE MPI_TRANS_R4

    SUBROUTINE MPI_TRANS_I1(A,B,ORIG,DEST)
        !ORIG SEND THE MESSAGE A TO DEST
        !DEST RECV THE MESSAGE B FROM ORIG
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        integer,INTENT(INOUT) :: A(:)
        integer,INTENT(IN)    :: B(:)
        integer,allocatable :: buff(:) ! security necessary
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
          STATUS=MPI_STATUS_IGNORE
          IF(RANK==ORIG)  then ! I'M ORIG, I SEND THE MESSAGE B TO DEST
              allocate(buff(size(b)))
              buff=b
              CALL MPI_SEND(Buff(1),SIZE(B),MPI_INTEGER,DEST, &
              TAG,MPI_COMM_WORLD,IERR)
              deallocate(buff)
          endif
          IF(RANK==DEST) then  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
              allocate(buff(size(a)))
              CALL MPI_RECV(buff(1),SIZE(A),MPI_INTEGER,ORIG, &
              TAG,MPI_COMM_WORLD,STATUS,IERR)
              a=buff
              deallocate(buff)
          endif
#endif
        ENDIF

    END SUBROUTINE MPI_TRANS_I1

    SUBROUTINE MPI_TRANS_I2(A,B,ORIG,DEST)
        !ORIG SEND THE MESSAGE A TO DEST
        !DEST RECV THE MESSAGE B FROM ORIG
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        integer,INTENT(INOUT) :: A(:,:)
        integer,INTENT(IN)    :: B(:,:)
        integer,allocatable :: buff(:,:) ! security necessary
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
          STATUS=MPI_STATUS_IGNORE
          IF(RANK==ORIG)  then ! I'M ORIG, I SEND THE MESSAGE B TO DEST
              allocate(buff(size(b,1),size(b,2)))
              buff=b
              CALL MPI_SEND(Buff(1,1),SIZE(B),MPI_INTEGER,DEST, &
              TAG,MPI_COMM_WORLD,IERR)
              deallocate(buff)
          endif
          IF(RANK==DEST) then  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
              allocate(buff(size(a,1),size(a,2)))
              CALL MPI_RECV(buff(1,1),SIZE(A),MPI_INTEGER,ORIG, &
              TAG,MPI_COMM_WORLD,STATUS,IERR)
              a=buff
              deallocate(buff)
          endif
#endif
        ENDIF

    END SUBROUTINE MPI_TRANS_I2


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
          STATUS=MPI_STATUS_IGNORE
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
          STATUS=MPI_STATUS_IGNORE
          IF(RANK==ORIG)  & ! I'M ORIG, I SEND THE MESSAGE B TO DEST
              CALL MPI_SEND(B,len(b),MPI_CHARACTER,DEST, &
              TAG,MPI_COMM_WORLD,IERR)

          IF(RANK==DEST) &  ! I'M DEST, I RECIEVE THE MESSAGE A FORM ORIG
              CALL MPI_RECV(A,len(a),MPI_CHARACTER,ORIG, &
              TAG,MPI_COMM_WORLD,STATUS,IERR)
#endif
        ENDIF

    END SUBROUTINE MPI_TRANS_C0

    SUBROUTINE MPI_ITRANS2_R1(A,ORIG,DEST,REQ,TAG)
        !ORIG SEND THE MESSAGE A TO DEST
        !DEST RECV THE MESSAGE B FROM ORIG
        !RETURN WHEN EVERYTHING IS STARTED
        !(THE MESSAGES HAVEN'T BEEN DELIVERED YET)
        use sortiefichier
        IMPLICIT NONE
        double precision   ,INTENT(INOUT) :: A(:)
        integer,INTENT(IN)    :: ORIG,DEST,TAG
        integer,INTENT(INOUT) :: REQ
        integer :: ierr

#ifdef WITH_MPI
        if (size(A)>0) then
          if (rank==ORIG) &
              CALL MPI_ISEND(A(1), SIZE(A), MPI_REAL8, DEST, &
              TAG, MPI_COMM_WORLD,REQ,IERR)! SEND THE BUFFER

          if (rank==dest) &
              CALL MPI_IRECV(A(1), SIZE(A), MPI_REAL8, ORIG, &
              TAG, MPI_COMM_WORLD,REQ,IERR)! RECV THE BUFFER
        endif
#endif

    END SUBROUTINE MPI_ITRANS2_R1


    SUBROUTINE MPI_ITRANS2_R2(A,ORIG,DEST,REQ,TAG)
        !ORIG SEND THE MESSAGE A TO DEST
        !DEST RECV THE MESSAGE B FROM ORIG
        !RETURN WHEN EVERYTHING IS STARTED
        !(THE MESSAGES HAVEN'T BEEN DELIVERED YET)
        IMPLICIT NONE
        double precision   ,INTENT(INOUT) :: A(:,:)
        integer,INTENT(IN)    :: ORIG,DEST,TAG
        integer,INTENT(INOUT) :: REQ
        integer :: ierr

#ifdef WITH_MPI
        if (size(A)>0) then

          if (rank==ORIG) &
              CALL MPI_ISEND(A(1,1), SIZE(A), MPI_REAL8, DEST, &
              TAG, MPI_COMM_WORLD,REQ,IERR)! SEND THE BUFFER

          if (rank==dest) &
              CALL MPI_IRECV(A(1,1), SIZE(A), MPI_REAL8, ORIG, &
              TAG, MPI_COMM_WORLD,REQ,IERR)! RECV THE BUFFER
        endif
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
        STATUS=MPI_STATUS_IGNORE
        IF(REQ/=MPI_REQUEST_NULL)   CALL MPI_WAIT(REQ, STATUS, IERR)
#endif

    END SUBROUTINE WAIT_MPI

    SUBROUTINE BARRIER
        ! SYNC POINT FOR ALL PROCESS
        use sortiefichier
        IMPLICIT NONE
        integer :: ierr

#ifdef WITH_MPI
        CALL MPI_Barrier(MPI_COMM_WORLD,IERR )
#endif

        flush(imp)
        flush(stderr)

#ifdef WITH_MPI
        CALL MPI_Barrier(MPI_COMM_WORLD,IERR )
#endif

    END SUBROUTINE BARRIER

    SUBROUTINE GATHER_I0(IN,OUT,COMM)
        ! THIS ROUTINE IS GATHERING INTEGER FOR EVERY PROC FROM EVERY PROC
        ! EXEMPLE
        ! PROC 1 : IN=A , OUT=[A , B , C]
        ! PROC 2 : IN=B , OUT=[A , B , C]
        ! PROC 3 : IN=C , OUT=[A , B , C]
        IMPLICIT NONE
        integer,INTENT(IN)  :: IN
        integer,INTENT(OUT) :: OUT(nprocs)
        integer,INTENT(IN),optional  :: COMM
        integer :: ierr
#ifdef WITH_MPI
        if (present(comm)) then
          CALL MPI_ALLGATHER(IN, 1, MPI_INTEGER, OUT(1), 1, MPI_INTEGER, COMM,IERR)
        else
          CALL MPI_ALLGATHER(IN, 1, MPI_INTEGER, OUT(1), 1, MPI_INTEGER, MPI_COMM_WORLD,IERR)
        endif
#else
        OUT(:)=IN
#endif

    END SUBROUTINE GATHER_I0

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

        LOGICAL,INTENT(INOUT)  ::  IN
        LOGICAL,INTENT(OUT),optional :: OUT
        LOGICAL             :: TMP
        integer :: ierr
#ifdef WITH_MPI
        CALL MPI_ALLREDUCE(IN, TMP, 1, MPI_LOGICAL,MPI_LOR, MPI_COMM_WORLD,IERR)
        if (present(out))  then
          OUT=TMP
        else
          IN=TMP
        endif
#else
        if (present(out))  OUT=IN
#endif

    END SUBROUTINE LOR_MPI

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

    SUBROUTINE MIN_MPI_0R(A,B)
        !COMPUTE THE SUM OF A IN B
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        double precision   ,INTENT(INOUT)    :: A
        double precision   ,INTENT(OUT),optional :: B
#ifdef WITH_MPI
        double precision ::C
        integer :: ierr
        CALL MPI_ALLREDUCE(A, C, 1, MPI_REAL8,MPI_MIN, MPI_COMM_WORLD,IERR)
        if (present(B)) then
          B=C
        else
          A=C
        endif
#else
        if (present(B))  B=A
#endif
    END SUBROUTINE MIN_MPI_0R

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

    SUBROUTINE MAX_MPI_1R(A,B)
        !COMPUTE THE SUM OF A IN B
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        double precision   ,INTENT(INOUT)    :: A(:)
        double precision   ,INTENT(OUT),optional :: B(:)
#ifdef WITH_MPI
        double precision,allocatable ::C(:)
        integer :: ierr
        allocate(c(size(a)))
        CALL MPI_ALLREDUCE(A(1), C(1), SIZE(A), MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,IERR)
        if (present(B)) then
          B=C
        else
          A=C
        endif
#else
        if (present(B))  B=A
#endif
    END SUBROUTINE MAX_MPI_1R

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

    SUBROUTINE BCAST_0I(IN,ORIG,COMM)
        !BROADCAST THE MESSAGE IN FROM ORIG
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        integer   ,INTENT(INOUT) :: IN
        integer,INTENT(IN)    :: ORIG
        integer,intent(in),optional :: comm
        integer :: ierr
#ifdef WITH_MPI
        if (present(comm)) then
          CALL MPI_BCAST( IN,1, MPI_INTEGER, ORIG, comm,IERR)
        else
          CALL MPI_BCAST( IN,1, MPI_INTEGER, ORIG, MPI_COMM_WORLD,IERR)
        endif
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

    SUBROUTINE BCAST_1R(IN,ORIG,COMM)
        !BROADCAST THE MESSAGE IN FROM ORIG
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        double precision,INTENT(INOUT) :: IN(:)
        integer,INTENT(IN)    :: ORIG
        integer,intent(in),optional :: comm
        double precision,allocatable :: buff(:) ! securiy necessary
        integer :: ierr
#ifdef WITH_MPI
        allocate(buff(size(in)))
        buff=in
        if (present(comm)) then
          CALL MPI_BCAST( buff(1),SIZE(IN), MPI_REAL8, ORIG, comm,IERR)
        else
          CALL MPI_BCAST( buff(1),SIZE(IN), MPI_REAL8, ORIG, MPI_COMM_WORLD,IERR)
        endif
        in=buff
        deallocate(buff)
#endif

    END SUBROUTINE BCAST_1R

    SUBROUTINE BCAST_2R(IN,ORIG)
        !BROADCAST THE MESSAGE IN FROM ORIG
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        double precision   ,INTENT(INOUT) :: IN(:,:)
        integer,INTENT(IN)    :: ORIG
        double precision,allocatable :: buff(:,:) ! securiy necessary
        integer :: ierr
#ifdef WITH_MPI
        allocate(buff(size(in,1),size(in,2)))
        buff=in
        CALL MPI_BCAST( buff(1,1),SIZE(IN), MPI_REAL8, ORIG, MPI_COMM_WORLD,IERR)
        in=buff
        deallocate(buff)
#endif

    END SUBROUTINE BCAST_2R

    SUBROUTINE BCAST_1I(IN,ORIG,COMM)
        !BROADCAST THE MESSAGE IN FROM ORIG
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        integer,INTENT(INOUT) :: IN(:)
        integer,INTENT(IN)    :: ORIG
        integer,intent(in),optional :: comm
        integer,allocatable :: buff(:) ! securiy necessary
        integer :: ierr
#ifdef WITH_MPI
        allocate(buff(size(in)))
        buff=in
        if (present(comm)) then
          CALL MPI_BCAST( buff(1),SIZE(IN), MPI_INTEGER, ORIG, comm,IERR)
        else
          CALL MPI_BCAST( buff(1),SIZE(IN), MPI_INTEGER, ORIG, MPI_COMM_WORLD,IERR)
        endif
        in=buff
        deallocate(buff)
#endif

    END SUBROUTINE BCAST_1I

    SUBROUTINE BCAST_2I(IN,ORIG)
        !BROADCAST THE MESSAGE IN FROM ORIG
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        integer,INTENT(INOUT) :: IN(:,:)
        integer,INTENT(IN)    :: ORIG
        integer,allocatable :: buff(:,:) ! securiy necessary
        integer :: ierr
#ifdef WITH_MPI
        allocate(buff(size(in,1),size(in,2)))
        buff=in
        CALL MPI_BCAST( buff(1,1),SIZE(IN), MPI_INTEGER, ORIG, MPI_COMM_WORLD,IERR)
        in=buff
        deallocate(buff)
#endif

    END SUBROUTINE BCAST_2I


    SUBROUTINE BCAST_3I(IN,ORIG)
        !BROADCAST THE MESSAGE IN FROM ORIG
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        integer,INTENT(INOUT) :: IN(:,:,:)
        integer,INTENT(IN)    :: ORIG
        integer,allocatable :: buff(:,:,:) ! securiy necessary
        integer :: ierr
#ifdef WITH_MPI
        allocate(buff(size(in,1),size(in,2),size(in,3)))
        buff=in
        CALL MPI_BCAST( buff(1,1,1),SIZE(IN), MPI_INTEGER, ORIG, MPI_COMM_WORLD,IERR)
        in=buff
        deallocate(buff)
#endif

    END SUBROUTINE BCAST_3I

    SUBROUTINE BCAST_0C(IN,ORIG)
        !BROADCAST THE MESSAGE IN FROM ORIG
        !RETURN WHEN EVERYTHING IS DONE
        IMPLICIT NONE
        character(*),INTENT(INOUT) :: IN
        integer,INTENT(IN)    :: ORIG
        integer :: ierr
#ifdef WITH_MPI
        CALL MPI_BCAST( IN,len(in), MPI_CHARACTER, ORIG, MPI_COMM_WORLD,IERR)
#endif

    END SUBROUTINE BCAST_0C

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

    SUBROUTINE GATHER_C(IN,OUT,SIZE)
        ! THIS ROUTINE IS GATHERING CHARCTER FOR EVERY PROC FROM EVERY PROC
        ! EXEMPLE
        ! PROC 1 : IN=[A,B,C,D] , SIZE=3 , OUT=[A,B,C,E,F]
        ! PROC 2 : IN=[E,F]     , SIZE=2 , OUT=[A,B,C,E,F]
        ! PROC 3 : IN=[]        , SIZE=0 , OUT=[A,B,C,E,F]
        IMPLICIT NONE
        integer,INTENT(IN)       :: SIZE
        character(*),INTENT(IN)  :: IN(:)
        character(*),INTENT(OUT) :: OUT(:)
        integer             :: I,DISP(NPROCS),SIZE_ALL(NPROCS)
        integer :: ierr
#ifdef WITH_MPI
        CALL MPI_ALLGATHER(SIZE, 1, MPI_INTEGER, &
            SIZE_ALL, 1, MPI_INTEGER, MPI_COMM_WORLD,IERR)
        SIZE_ALL=SIZE_ALL*len(IN)
        DISP(1)=0
        DO I=2,NPROCS
          DISP(I)=DISP(I-1) + SIZE_ALL(I-1)
        ENDDO

        CALL MPI_ALLGATHERV(IN(1), SIZE*len(IN), MPI_character, &
            OUT(1), SIZE_ALL, DISP, MPI_character, MPI_COMM_WORLD,IERR)
#else
        OUT=IN
#endif

    END SUBROUTINE GATHER_C

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


    SUBROUTINE START_KEEP_ORDER(index,queue,relais_in)
        IMPLICIT NONE
        integer,intent(inout),optional :: relais_in
        integer,intent(in),optional :: index,queue(:)
        integer :: relais
        relais=0
        if (present(relais_in)) relais=relais_in
#ifdef WITH_MPI
        if(present(index)) then
          if(index>1) call mpi_trans(relais,relais,queue(index-1),rank)
        else
          if(rank>0) call mpi_trans(relais,relais,rank-1,rank)
        endif
#endif
        if (present(relais_in)) relais_in=relais
    END SUBROUTINE START_KEEP_ORDER


    SUBROUTINE END_KEEP_ORDER(index,queue,relais_in)
        IMPLICIT NONE
        integer,intent(inout),optional :: relais_in
        integer,intent(in),optional :: index,queue(:)
        integer :: relais
        if (present(relais_in)) relais=relais_in
#ifdef WITH_MPI
        if(present(index)) then
          if(index<size(queue)) call mpi_trans(relais,relais,rank,queue(index+1))
        else
          if(rank<nprocs-1) call mpi_trans(relais,relais,rank,rank+1)
        endif
#endif
        if (present(relais_in)) relais_in=relais
    END SUBROUTINE END_KEEP_ORDER

    subroutine my_fseek(unit,pos)
        implicit none
        integer :: unit,pos
#if defined(__GFORTRAN__)
        call fseek(unit,pos,0)
#else
        integer :: i
        i=fseek(unit,pos,0)
#endif
    end subroutine my_fseek


    subroutine get_time(time)
        implicit none
        double precision,intent(out) :: time
        logical,parameter :: wall_time=.true.
#ifdef WITH_MPI
        if (wall_time) then
          time=MPI_Wtime()
        else
          call CPU_TIME(time)
          call sum_mpi(time)                      ! THIS MIGHT NOT BE RELEVANT
        endif
#else
        integer(kind=8) :: clock_rate,clock
        if (wall_time) then
          call system_clock(clock,clock_rate)
          time=clock/clock_rate
        else
          call CPU_TIME(time)
        endif
#endif
    end subroutine get_time
  end module mod_mpi

