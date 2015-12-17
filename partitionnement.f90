module mod_partitionnement
  use tools
  implicit none
  integer         ,parameter :: verbosity=2  ! from 0 to 3
  double precision,parameter :: eps=1e-10
contains
  subroutine partitionnement(x,y,z,ncbd,mnc,ncin,bceqt,exs1,exs2,nblocks,nsplit,nsplit_dir)
    use boundary
    use para_var
    use kcle,only:klzx,kmtbx
    use sortiefichier
    use mod_inbdc
    use mod_c_inbdc
    use mod_inbdb
    use mod_c_inbdb
    use mod_mpi
    !
    !***********************************************************************
    !
    !     act
    !_a    realisation du partitionnement
    !      ecrit par Alexandre Poux - sept 2015
    !
    !***********************************************************************
    !
    !   for the index, we use the notation abci :
    !   a being :
    !    - l    for a block
    !    - fr   for a boundary
    !    - proc for a proc
    !   b being :
    !    - l if it's local (inside a proc)
    !    - g if it's global (shared by all proc)
    !    - nothing if a is proc
    !   c being :
    !    - i for the initial index (before partitionning)
    !    - f for the final index (after partitionning)
    !   i being an arbirary number
    !
    implicit none
    double precision,allocatable,intent(inout) :: bceqt(:,:)
    integer         ,allocatable,intent(inout) :: ncin(:),mnc(:),ncbd(:)
    double precision,allocatable,intent(inout) :: x(:),y(:),z(:)
    double precision            ,intent(inout) :: exs1,exs2
    integer,intent(inout) :: nsplit(:),nsplit_dir(:,:)

    integer             :: save_lt,save_num_bcg,save_num_bg
    integer             :: nblocks,nxyza
    integer             :: imot(nmx),nmot
    integer             :: nmin,nmax,nmin1,nmax1
    integer             :: imax,imin,jmax,jmin,kmax,kmin
    integer             :: nid,njd,nijd,sblock
    integer             :: nsub
    integer             :: i2,j2,k2,i,j,k
    integer             :: xs,ys,zs,xe,ye,ze
    integer             :: lli1,lgi1,lgi2,lli2
    integer             :: llf1,lgf1,llf2,lgf2,llf3,lgf3
    integer             :: frli1,frli2,frgi1,frgi2,frgi3
    integer             :: frlf1,frgf1,frgf2,frlf3,frgf3,frgf4
    integer             :: proci1,procf1,proci2,procf2,proc,load(0:nprocs-1),nload(0:nprocs-1)
    integer             :: dir(3,3),offset(3),iba,ii,jba,jj,kba,kk,na,nb,nbi,nbj,nbk,i1,j1,k1

    integer,allocatable :: save_ii2(:),save_jj2(:),save_kk2(:)
    integer,allocatable :: save_id1(:),save_jd1(:),save_kd1(:)
    integer,allocatable :: save_id2(:),save_jd2(:),save_kd2(:)
    integer,allocatable :: save_npn(:)
    integer,allocatable :: save_iminb(:),save_jminb(:),save_kminb(:)
    integer,allocatable :: save_imaxb(:),save_jmaxb(:),save_kmaxb(:)
    integer,allocatable :: save_bcg_to_bcl(:),save_bcg_to_bci(:),save_bcg_to_bg(:)
    integer,allocatable :: save_bg_to_proc(:),save_bg_to_bl(:),save_bg_to_bi(:),save_bl_to_bg(:)
    integer,allocatable :: ni(:,:,:,:),nj(:,:,:,:),nk(:,:,:,:),nijk(:,:,:,:)
    integer,allocatable :: ii2g(:),jj2g(:),kk2g(:),sub_bc(:,:),tmp(:)

    double precision    :: vbc(ista*lsta),sub_bc_c(2,6),sub_bc2(4,3)
    double precision    :: dist
    double precision,allocatable :: save_x(:),save_y(:),save_z(:)

    character(len=2),allocatable :: save_indfl(:)
    character(len=32) ::  mot(nmx)
    character(len=2) :: indmf

    logical :: test

    !############################################################################################
    !############################## GET PARAMETERS ##############################################
    !############################################################################################

    nblocks=max(nblocks,nprocs) ! minimum number of bloc that we need
    nblocks=max(nblocks,sum(nsplit)) ! number of bloc that we need
    nblocks=max(nblocks,num_bg) ! number of bloc that we need

    if(.true.) then ! always do the partitionning in order to reequilibrate
    !if(nblocks/=num_bg) then ! there is some partitionning to do

       !############################################################################################
       !############################# SAVE OLD MESH FOR CHECKING PURPOSE ###########################
       !############################################################################################
       if(verbosity>=3) call write_mesh("orig",x,y,z)

       !############################################################################################
       !################# SAVE ALL MESH AND BOUNDARIES VARIABLES ###################################
       !#################            AND CLEAR EVERYTHING        ###################################
       !#################       WE WILL RECREATE EVERYTHING      ###################################
       !############################################################################################

       ! allocate statments are non necessary with newer fortran norm, but ifort don't support it yet
       allocate(save_x(size(x)))
       allocate(save_y(size(y)))
       allocate(save_z(size(z)))
       allocate(save_indfl(size(indfl)))
       allocate(save_iminb(size(iminb)))
       allocate(save_imaxb(size(imaxb)))
       allocate(save_jminb(size(jminb)))
       allocate(save_jmaxb(size(jmaxb)))
       allocate(save_kminb(size(kminb)))
       allocate(save_kmaxb(size(kmaxb)))
       allocate(save_ii2(size(ii2)))
       allocate(save_jj2(size(jj2)))
       allocate(save_kk2(size(kk2)))
       allocate(save_id1(size(id1)))
       allocate(save_jd1(size(jd1)))
       allocate(save_kd1(size(kd1)))
       allocate(save_id2(size(id2)))
       allocate(save_jd2(size(jd2)))
       allocate(save_kd2(size(kd2)))
       allocate(save_npn(size(npn)))
       allocate(save_bcg_to_bcl(size(bcg_to_bcl)))
       allocate(save_bcg_to_bci(size(bcg_to_bci)))
       allocate(save_bg_to_proc(size(bg_to_proc)))
       allocate(save_bg_to_bl(size(bg_to_bl)))
       allocate(save_bg_to_bi(size(bg_to_bi)))
       allocate(save_bl_to_bg(size(bl_to_bg)))
       allocate(save_bcg_to_bg(size(bcg_to_bg)))

       save_lt=lt
       save_x=x
       save_y=y
       save_z=z
       save_ii2 = ii2
       save_jj2 = jj2
       save_kk2 = kk2
       save_id1 = id1
       save_jd1 = jd1
       save_kd1 = kd1
       save_id2 = id2
       save_jd2 = jd2
       save_kd2 = kd2
       save_npn = npn
       save_indfl=indfl
       save_iminb=iminb
       save_imaxb=imaxb
       save_jminb=jminb
       save_jmaxb=jmaxb
       save_kminb=kminb
       save_kmaxb=kmaxb
       save_num_bcg=num_bcg
       save_num_bg=num_bg
       save_bcg_to_bcl=bcg_to_bcl
       save_bcg_to_bci=bcg_to_bci
       save_bg_to_proc=bg_to_proc
       save_bg_to_bl=bg_to_bl
       save_bg_to_bi=bg_to_bi
       save_bl_to_bg=bl_to_bg
       save_bcg_to_bg=bcg_to_bg

       !   reinitialisation

       call reallocate(x,0)
       call reallocate(y,0)
       call reallocate(z,0)
       call reallocate(ndlb,0)
       call reallocate(nfei,0)
       call reallocate(indfl,0)
       call reallocate(iminb,0)
       call reallocate(imaxb,0)
       call reallocate(jminb,0)
       call reallocate(jmaxb,0)
       call reallocate(kminb,0)
       call reallocate(kmaxb,0)
       call reallocate(mpb,0)
       call reallocate(mmb,0)
       call reallocate(ncbd,0)
       call reallocate(ii1,0)
       call reallocate(jj1,0)
       call reallocate(kk1,0)
       call reallocate(ii2,0)
       call reallocate(jj2,0)
       call reallocate(kk2,0)
       call reallocate(id1,0)
       call reallocate(jd1,0)
       call reallocate(kd1,0)
       call reallocate(id2,0)
       call reallocate(jd2,0)
       call reallocate(kd2,0)
       call reallocate(nnn,0)
       call reallocate(nnc,0)
       call reallocate(nnfb,0)
       call reallocate(npn,0)
       call reallocate(npc,0)
       call reallocate(npfb,0)
       call reallocate(bcg_to_proc,0)
       call reallocate(bcg_to_bcl,0)
       call reallocate(bcg_to_bci,0)
       call reallocate(bcl_to_bcg,0)
       call reallocate(bg_to_proc,0)
       call reallocate(bg_to_bl,0)
       call reallocate(bg_to_bi,0)
       call reallocate(bl_to_bg,0)
       call reallocate(bcg_to_bg,0)

       ndimubx = 0
       ndimctbx=0
       ndimntbx=0
       lzx =0
       klzx=0
       lt=0
       mtbx=0
       mtcx=0
       kmtbx=0
       mtb=0
       mtt=0
       mdimubx=0
       mdimtbx=0
       frgf1=0
       ip41=0
       vbc=0.
       num_bcg=0
       num_bci=0
       num_bcl=0
       num_bg=0
       num_bi=0
       num_bl=0

       !############################################################################################
       !####################### COMPUTE SPLITTING ##################################################
       !############################################################################################

       !counting total number of points
       allocate(ii2g(save_num_bg),jj2g(save_num_bg),kk2g(save_num_bg))

       ii2g=0
       jj2g=0
       kk2g=0

       ! Gather all block size
       ! in order to keep the order, a gather call, won't do the trick
       do lli1=1,save_lt
          lgi1=save_bl_to_bg(lli1)
          ii2g(lgi1)=save_ii2(lli1)
          jj2g(lgi1)=save_jj2(lli1)
          kk2g(lgi1)=save_kk2(lli1)
       enddo
       call sum_mpi(ii2g) ! for each block, only one proc has a non-zero value
       call sum_mpi(jj2g) ! thus this will gather all the block size
       call sum_mpi(kk2g)

       nxyza=sum(ii2g*jj2g*kk2g)  ! total number of points

       ! routine calculant combiens de fois splitter chaque block
       ! sortie : nsplit   : nombre de splits de chaque block initial
       !          nsplit_dir  : nombre de splits par direction de chaque block initial
       !          ni,nj,nk  : nombre de points par direction de chaque nouveau block
       ! entrée : tout le reste
       call compute_split(nsplit,nsplit_dir,ni,nj,nk,                  &
            nxyza,save_num_bg,nblocks,save_lt,          &
            save_ii2 ,save_jj2 ,save_kk2 ,save_bl_to_bg,&
            ii2g,jj2g,kk2g,save_bg_to_proc)

       deallocate(ii2g,jj2g,kk2g)

       sblock=sum(nsplit_dir(1,:)*nsplit_dir(2,:)*nsplit_dir(3,:))
       if(sblock/=nblocks) then
          stop 'partitionnement impossible'
       else
          if(verbosity>=1.and.rank==0) then
             print*,''
             print*,'découpage réussis : '
             print*,"Block | Nb block x | Nb block y | Nb block z | unbalance"
             nmax=0
             nmin=nxyza
             do lgi1=1,save_num_bg
                nmax1=maxval(ni(:nsplit_dir(1,lgi1),:nsplit_dir(2,lgi1),:nsplit_dir(3,lgi1),lgi1)* &
                     nj(:nsplit_dir(1,lgi1),:nsplit_dir(2,lgi1),:nsplit_dir(3,lgi1),lgi1)* &
                     nk(:nsplit_dir(1,lgi1),:nsplit_dir(2,lgi1),:nsplit_dir(3,lgi1),lgi1))
                nmin1=minval(ni(:nsplit_dir(1,lgi1),:nsplit_dir(2,lgi1),:nsplit_dir(3,lgi1),lgi1)* &
                     nj(:nsplit_dir(1,lgi1),:nsplit_dir(2,lgi1),:nsplit_dir(3,lgi1),lgi1)* &
                     nk(:nsplit_dir(1,lgi1),:nsplit_dir(2,lgi1),:nsplit_dir(3,lgi1),lgi1))
                write(*,500) lgi1,nsplit_dir(:,lgi1),( nmax1- nmin1 )*100./ nmin1
                nmax=max(nmax,nmax1)
                nmin=min(nmin,nmin1)
             enddo
             write(*,'(A,I4,A,F6.2,A)') " Total : ",sblock," blocks, unbalance : ",( nmax- nmin )*100./ nmin,"%"
          endif
       end if
500    format(4x,i2," | ",6x,i4," | ",6x,i4," | ",6x,i4," | ",F5.2,"%")

       !############################################################################################
       !######################### PLACE BLOCKS ON PROCS ############################################
       !############################################################################################

       ! construct new bg_to_proc
       call reallocate(bg_to_proc,sblock)
       allocate(nijk(size(ni,1),size(ni,2),size(ni,3),size(ni,4)))
       allocate(tmp(4))
       nijk=ni*nj*nk ! size of the new blocks before attribution
       load=0
       nload=0

       do while (maxval(nijk)/=0)                        ! there is still block to be attributed
         proc=minloc(load,1)-1                           ! the process with the minimum load
         tmp=maxloc(nijk)                                ! take the biggest non attributed block
         lgf1=old2new_b(tmp(1),tmp(2),tmp(3),tmp(4))
         nload(proc)=nload(proc)+1
         load(proc)=load(proc)+nijk(tmp(1),tmp(2),tmp(3),tmp(4))
         bg_to_proc(lgf1)=proc
         nijk(tmp(1),tmp(2),tmp(3),tmp(4))=0             ! block is attributed
       enddo
       deallocate(nijk)
       deallocate(tmp)

      if(verbosity>=1.and.rank==0) then
         print*,''
         print*,"Proc | Nb blocks | Load"
         do proc=0,nprocs-1
            write(*,'(1x,i4," | ",5x,i4," | ",i8)') proc,nload(proc),load(proc)
         enddo
         write(*,'(A,I8,A,I4,A,I4,A,F6.2,A)') " Total : ",sum(load)," dispatched on "&
                                    ,sum(nload)," blocks on "&
                                    ,nprocs," process, unbalance : ",( maxval(load)- minval(load) )*100./ minval(load),"%"
      endif

       !############################################################################################
       !######################### RECREATE GRID ####################################################
       !############################################################################################


       !    print*,'create new grid '
       do lgi1=1,save_num_bg     !     create new grid and initialize it
          do k=1,nsplit_dir(3,lgi1)!  tout le monde parcours les nouveaux blocs dans le même ordre
             do j=1,nsplit_dir(2,lgi1)
                do i=1,nsplit_dir(1,lgi1)

                   lgf1=old2new_b(i,j,k,lgi1)

                   call create_grid(lgf1,ni(i,j,k,lgi1),nj(i,j,k,lgi1),nk(i,j,k,lgi1),&
                          save_bg_to_bi(lgi1),bg_to_proc(lgf1))

                   proci1=save_bg_to_proc(lgi1)
                   procf1=bg_to_proc(lgf1)

                   if(rank==proci1.or.rank==procf1) then ! TODO : non blocking communication, removing need of copy_grid
                      if(proci1==procf1) then ! sending to myself
                         lli1=save_bg_to_bl(lgi1)
                         llf1=bg_to_bl(lgf1)

                         call new2old_p(0,0,0,xs,ys,zs,i,j,k,lgi1)

                         call reallocate_s(x,ndimntbx)
                         call reallocate_s(y,ndimntbx)
                         call reallocate_s(z,ndimntbx)

                         call copy_grid(lli1,llf1,xs,ys,zs,x,y,z,save_x,save_y,save_z,              &
                              save_id1,save_id2,save_jd1,save_jd2,save_kd1,save_kd2,save_npn)

                      elseif(rank==proci1) then
                         lli1=save_bg_to_bl(lgi1)

                         call new2old_p(1,1,1,xs,ys,zs,i,j,k,lgi1)
                         call new2old_p(ni(i,j,k,lgi1),nj(i,j,k,lgi1),nk(i,j,k,lgi1),xe,ye,ze,i,j,k,lgi1)

                         nid = save_id2(lli1)-save_id1(lli1)+1
                         njd = save_jd2(lli1)-save_jd1(lli1)+1
                         nijd = nid*njd

                         call send_grid(save_x,save_y,save_z,xs,ys,zs,xe,ye,ze,       &
                              save_id1(lli1),save_jd1(lli1),save_kd1(lli1),nid,nijd,save_npn(lli1),proci1,procf1)

                      elseif(rank==procf1) then
                         llf1=bg_to_bl(lgf1)

                         call recv_grid(x,y,z,llf1,proci1,procf1)
                      endif
                   endif

                enddo
             enddo
          enddo
       enddo

       !############################################################################################
       !######################### RECREATE OLD BOUNDARIES ##########################################
       !############################################################################################
       allocate(tmp(6))
       !    print*,'recreate old boundaries '
       do frgi1=1,save_num_bcg
          frli1=save_bcg_to_bcl(frgi1) ! tout le monde parcours les nouvelles condition limites dans le même ordre
          lgi1=save_bcg_to_bg(frgi1)
          lli1=save_bg_to_bl(lgi1)
          proci1=save_bg_to_proc(lgi1)
          do k=1,nsplit_dir(3,lgi1)
             do j=1,nsplit_dir(2,lgi1)
                do i=1,nsplit_dir(1,lgi1)

                   lgf1=old2new_b(i,j,k,lgi1)
                   llf1=bg_to_bl(lgf1)

                   nsub=0 ! count sub_boundaries
                   call reallocate(sub_bc,6,1)

                   if (rank==proci1) then

                      call new2old_p(1,1,1,xs,ys,zs,i,j,k,lgi1)
                      call new2old_p(ni(i,j,k,lgi1),nj(i,j,k,lgi1),nk(i,j,k,lgi1),xe,ye,ze,i,j,k,lgi1)

                      if ( save_iminb(frli1)<=xe .and. &
                           save_imaxb(frli1)>=xs .and. &
                           save_jminb(frli1)<=ye .and. &
                           save_jmaxb(frli1)>=ys .and. & ! there is a part of the boundary in this block
                           save_kminb(frli1)<=ze .and. &
                           save_kmaxb(frli1)>=zs ) then

                         ! part of the boundary which concern this block

                         tmp(1)=min(xe,max(xs,save_iminb(frli1))) ! coordinate of the new boundary
                         tmp(2)=min(xe,max(xs,save_imaxb(frli1))) ! in the old block ref
                         tmp(3)=min(ye,max(ys,save_jminb(frli1)))
                         tmp(4)=min(ye,max(ys,save_jmaxb(frli1)))
                         tmp(5)=min(ze,max(zs,save_kminb(frli1)))
                         tmp(6)=min(ze,max(zs,save_kmaxb(frli1)))

                         if((tmp(2)-tmp(1)>0.and.tmp(4)-tmp(3)>0).or. &
                            (tmp(6)-tmp(5)>0.and.tmp(4)-tmp(3)>0).or. &
                            (tmp(2)-tmp(1)>0.and.tmp(6)-tmp(5)>0)) then

                           nsub=1
                           sub_bc(:,1)=tmp
                           indmf=save_indfl(frli1)
                         endif
                      endif
                   endif
                   call bcast(nsub,proci1)

                   if (nsub==0) cycle

                   if(tab_raccord(frgi1)/=0) then ! if boundary is shared with another block
                      frgi2=tab_raccord(frgi1)      ! a split on the other side induce split here
                      frli2=save_bcg_to_bcl(frgi2)
                      lgi2=save_bcg_to_bg(frgi2)
                      lli2=save_bg_to_bl(lgi2)
                      proci2=save_bg_to_proc(lgi2)

                      if(rank==proci1) then
                         nid = save_id2(lli1)-save_id1(lli1)+1
                         njd = save_jd2(lli1)-save_jd1(lli1)+1
                         nijd = nid*njd

                          nb = save_npn(lli1)+1+(sub_bc(1,1)-save_id1(lli1)) &
                                               +(sub_bc(3,1)-save_jd1(lli1))*nid &
                                               +(sub_bc(5,1)-save_kd1(lli1))*nijd
!
                          nbi = nb+1
                          nbj = nb+nid
                          nbk = nb+nijd

                          sub_bc2(1,1)   =save_x(nb)
                          sub_bc2(2,1)   =save_x(nbi)
                          sub_bc2(3,1)   =save_x(nbj)
                          sub_bc2(4,1)   =save_x(nbk)
                          sub_bc2(1,2)   =save_y(nb)!   coordinate of the first point of the new boundary
                          sub_bc2(2,2)   =save_y(nbi)
                          sub_bc2(3,2)   =save_y(nbj)
                          sub_bc2(4,2)   =save_y(nbk)
                          sub_bc2(1,3)   =save_z(nb)
                          sub_bc2(2,3)   =save_z(nbi)
                          sub_bc2(3,3)   =save_z(nbj)
                          sub_bc2(4,3)   =save_z(nbk)

                         sub_bc(1,1)=sub_bc(1,1)-save_iminb(frli1)
                         sub_bc(2,1)=sub_bc(2,1)-save_iminb(frli1) ! index of the new boundary
                         sub_bc(3,1)=sub_bc(3,1)-save_jminb(frli1) ! in the old boundary ref
                         sub_bc(4,1)=sub_bc(4,1)-save_jminb(frli1)
                         sub_bc(5,1)=sub_bc(5,1)-save_kminb(frli1)
                         sub_bc(6,1)=sub_bc(6,1)-save_kminb(frli1)
                      endif

                      call MPI_TRANS(sub_bc,sub_bc,proci1,proci2)
                      call MPI_TRANS(sub_bc2,sub_bc2,proci1,proci2)

                      if(rank==proci2) then
                         nsub=0
                         dir=0
                         nid = save_id2(lli2)-save_id1(lli2)+1
                         njd = save_jd2(lli2)-save_jd1(lli2)+1
                         nijd = nid*njd

                         ! look for the first point
                          search1:do k1=save_kminb(frli2),save_kmaxb(frli2)
                           do j1=save_jminb(frli2),save_jmaxb(frli2)
                            do i1=save_iminb(frli2),save_imaxb(frli2)
                             iba=i1
                             jba=j1
                             kba=k1
                            na = save_npn(lli2)+1+(i1-save_id1(lli2)) &
                                                 +(j1-save_jd1(lli2))*nid &
                                                 +(k1-save_kd1(lli2))*nijd
                             dist=sqrt( (save_x(na)-sub_bc2(1,1))**2+(save_y(na)-sub_bc2(1,2))**2+(save_z(na)-sub_bc2(1,3))**2 )
                             if(dist.lt.eps) exit search1
                            enddo
                           enddo
                          enddo search1
                          if(dist.ge.eps) stop "problem in boundary "
                         
                          ! look for directions
                          search2:do k1=save_kminb(frli2),save_kmaxb(frli2)
                           do j1=save_jminb(frli2),save_jmaxb(frli2)
                            do i1=save_iminb(frli2),save_imaxb(frli2)
                             ii=i1
                             jj=j1
                             kk=k1
                            na = save_npn(lli2)+1+(i1-save_id1(lli2)) &
                                                 +(j1-save_jd1(lli2))*nid &
                                                 +(k1-save_kd1(lli2))*nijd
                             dist=sqrt( (save_x(na)-sub_bc2(2,1))**2+(save_y(na)-sub_bc2(2,2))**2+(save_z(na)-sub_bc2(2,3))**2 )
                             if(dist.lt.eps) exit search2
                            enddo
                           enddo
                          enddo search2
                          if(dist.lt.eps) then
                            dir(1,1)=ii-iba ; dir(1,1)=ii-iba
                            dir(1,2)=jj-jba ; dir(2,1)=jj-jba
                            dir(1,3)=kk-kba ; dir(3,1)=kk-kba
                          endif
                          search3:do k1=save_kminb(frli2),save_kmaxb(frli2)
                           do j1=save_jminb(frli2),save_jmaxb(frli2)
                            do i1=save_iminb(frli2),save_imaxb(frli2)
                             ii=i1
                             jj=j1
                             kk=k1
                            na = save_npn(lli2)+1+(i1-save_id1(lli2)) &
                                                 +(j1-save_jd1(lli2))*nid &
                                                 +(k1-save_kd1(lli2))*nijd
                             dist=sqrt( (save_x(na)-sub_bc2(3,1))**2+(save_y(na)-sub_bc2(3,2))**2+(save_z(na)-sub_bc2(3,3))**2 )
                             if(dist.lt.eps) exit search3
                            enddo
                           enddo
                          enddo search3
                          if(dist.lt.eps) then
                            dir(2,1)=ii-iba ; dir(1,2)=ii-iba
                            dir(2,2)=jj-jba ; dir(2,2)=jj-jba
                            dir(2,3)=kk-kba ; dir(3,2)=kk-kba
                          endif
                          search4:do k1=save_kminb(frli2),save_kmaxb(frli2)
                           do j1=save_jminb(frli2),save_jmaxb(frli2)
                            do i1=save_iminb(frli2),save_imaxb(frli2)
                             ii=i1
                             jj=j1
                             kk=k1
                            na = save_npn(lli2)+1+(i1-save_id1(lli2)) &
                                                 +(j1-save_jd1(lli2))*nid &
                                                 +(k1-save_kd1(lli2))*nijd
                             dist=sqrt( (save_x(na)-sub_bc2(4,1))**2+(save_y(na)-sub_bc2(4,2))**2+(save_z(na)-sub_bc2(4,3))**2 )
                             if(dist.lt.eps) exit search4
                            enddo
                           enddo
                          enddo search4
                          if(dist.lt.eps) then
                            dir(3,1)=ii-iba ; dir(1,3)=ii-iba
                            dir(3,2)=jj-jba ; dir(2,3)=jj-jba
                            dir(3,3)=kk-kba ; dir(3,3)=kk-kba
                          endif
                          
                          ! on change de repère

                          tmp(1:3)=matmul(dir,sub_bc(1:5:2,1))
                          tmp(4:6)=matmul(dir,sub_bc(2:6:2,1))
                          tmp(1:4:3)=tmp(1:4:3)+save_iminb(frli2) ! coordinate of the new boundary
                          tmp(2:5:3)=tmp(2:5:3)+save_jminb(frli2) ! in the old block ref
                          tmp(3:6:3)=tmp(3:6:3)+save_kminb(frli2)

                          offset=[iba,jba,kba]-tmp(1:3)
                         
                         imin=min(tmp(1),tmp(4))+offset(1)
                         imax=max(tmp(1),tmp(4))+offset(1)
                         jmin=min(tmp(2),tmp(5))+offset(2)
                         jmax=max(tmp(2),tmp(5))+offset(2)
                         kmin=min(tmp(3),tmp(6))+offset(3)
                         kmax=max(tmp(3),tmp(6))+offset(3)

                         do k2=1,nsplit_dir(3,lgi2)
                            do j2=1,nsplit_dir(2,lgi2)
                               do i2=1,nsplit_dir(1,lgi2)

                                  lgf2=old2new_b(i2,j2,k2,lgi2)
                                  llf2=bg_to_bl(lgf2)
                                  lli2=save_bg_to_bl(lgi2)

                                   call new2old_p(1,1,1,xs,ys,zs,i2,j2,k2,lgi2)
                                   call new2old_p(ni(i2,j2,k2,lgi2),nj(i2,j2,k2,lgi2),nk(i2,j2,k2,lgi2),xe,ye,ze,i2,j2,k2,lgi2)

                                  if ( imin<=xe .and. &
                                       imax>=xs .and. &
                                       jmin<=ye .and. &
                                       jmax>=ys .and. & ! there is a part of the boundary in this block
                                       kmin<=ze .and. &
                                       kmax>=zs ) then

                                     ! part of the boundary which concern this block

                                     tmp(1)=min(xe,max(xs,imin))-save_iminb(frli2) ! coordinate of the new boundary
                                     tmp(2)=min(xe,max(xs,imax))-save_iminb(frli2) ! in the old boundary ref
                                     tmp(3)=min(ye,max(ys,jmin))-save_jminb(frli2)
                                     tmp(4)=min(ye,max(ys,jmax))-save_jminb(frli2)
                                     tmp(5)=min(ze,max(zs,kmin))-save_kminb(frli2)
                                     tmp(6)=min(ze,max(zs,kmax))-save_kminb(frli2)

                                     if((tmp(2)-tmp(1)>0.and.tmp(4)-tmp(3)>0).or. &
                                        (tmp(6)-tmp(5)>0.and.tmp(4)-tmp(3)>0).or. &
                                        (tmp(2)-tmp(1)>0.and.tmp(6)-tmp(5)>0)) then

                                       nsub=nsub+1
                                       call reallocate_s(sub_bc,6,nsub)
                                       sub_bc(:,nsub)=tmp
                                     endif
                                  endif
                               enddo
                            enddo
                         enddo
                         if (nsub==0) then
                            print*,"Problem in the boundary ",frgi1,frgi2
                            call abort
                         endif
                         ! on rechange de repère
                         do i2=1,nsub
                            tmp=sub_bc(:,i2)
                            
                            tmp(1:3)=matmul(dir,sub_bc(1:5:2,i2)-offset)
                            tmp(4:6)=matmul(dir,sub_bc(2:6:2,i2)-offset)
                           
                             imin=min(tmp(1),tmp(4))
                             imax=max(tmp(1),tmp(4))
                             jmin=min(tmp(2),tmp(5))
                             jmax=max(tmp(2),tmp(5))
                             kmin=min(tmp(3),tmp(6))
                             kmax=max(tmp(3),tmp(6))
                         
                            sub_bc(:,i2)=[imin,imax,jmin,jmax,kmin,kmax]
                         enddo
                      endif
                      call MPI_TRANS(nsub,nsub,proci2,proci1)
                      if(rank==proci1) call reallocate_s(sub_bc,6,nsub)
                      call MPI_TRANS(sub_bc,sub_bc,proci2,proci1)
                      if(rank==proci1) then
                         sub_bc(1,:)=sub_bc(1,:)+save_iminb(frli1)
                         sub_bc(2,:)=sub_bc(2,:)+save_iminb(frli1) ! coordinate of the new boundary
                         sub_bc(3,:)=sub_bc(3,:)+save_jminb(frli1) ! in the old block ref
                         sub_bc(4,:)=sub_bc(4,:)+save_jminb(frli1)
                         sub_bc(5,:)=sub_bc(5,:)+save_kminb(frli1)
                         sub_bc(6,:)=sub_bc(6,:)+save_kminb(frli1)
                      endif
                   endif
                   call bcast(nsub,proci1)
                   if(rank/=proci1) call reallocate(sub_bc,6,nsub)
                   call bcast(sub_bc,proci1)
                   call bcast(indmf,proci1)
                   do i2=1,nsub
                      frgf1=frgf1+1
                      ! coordinate of the new boundary in the new block ref
                      call old2new_p(sub_bc(1,i2),sub_bc(3,i2),sub_bc(5,i2),sub_bc(1,i2),sub_bc(3,i2),sub_bc(5,i2),i,j,k,lgi1)
                      call old2new_p(sub_bc(2,i2),sub_bc(4,i2),sub_bc(6,i2),sub_bc(2,i2),sub_bc(4,i2),sub_bc(6,i2),i,j,k,lgi1)

                      call create_boundary(frgf1,lgf1,indmf,ncbd,save_bcg_to_bci(frgi1),&
                           sub_bc(1,i2),sub_bc(2,i2),sub_bc(3,i2),sub_bc(4,i2),sub_bc(5,i2),sub_bc(6,i2))

                   enddo
                enddo
             enddo
          enddo
       enddo
       deallocate(tmp)
       !############################################################################################
       !########################### CREATE NEW BOUNDARIES ##########################################
       !############################################################################################

       !    print*,'create new boundaries '
       call reallocate(sub_bc,6,2)
       do lgi1=1,save_num_bg
          do k=1,nsplit_dir(3,lgi1)
             do j=1,nsplit_dir(2,lgi1)
                do i=1,nsplit_dir(1,lgi1) !           New coincident boundaries
                   lgf1=old2new_b(i,j,k,lgi1)
                   llf1=bg_to_bl(lgf1)
                   procf1=bg_to_proc(lgf1)
                   sub_bc=0

                   if (i>1) then
                      lgf2=old2new_b(i-1,j,k,lgi1)
                      llf2=bg_to_bl(lgf2)
                      procf2=bg_to_proc(lgf2)

                      if (rank==procf1) then
                         sub_bc(1,1)=ii1(llf1)
                         sub_bc(2,1)=ii1(llf1)
                         sub_bc(3,1)=jj1(llf1)
                         sub_bc(4,1)=jj2(llf1)
                         sub_bc(5,1)=kk1(llf1)
                         sub_bc(6,1)=kk2(llf1)
                      endif
                      if (rank==procf2) then
                         sub_bc(1,2)=ii2(llf2)
                         sub_bc(2,2)=ii2(llf2)
                         sub_bc(3,2)=jj1(llf2)
                         sub_bc(4,2)=jj2(llf2)
                         sub_bc(5,2)=kk1(llf2)
                         sub_bc(6,2)=kk2(llf2)
                      endif

                      call bcast(sub_bc(:,1),procf1)
                      call bcast(sub_bc(:,2),procf2)

                      frgf1=frgf1+1
                      call create_boundary(frgf1,lgf1,"i1",ncbd,0,&
                           sub_bc(1,1),sub_bc(2,1),sub_bc(3,1),sub_bc(4,1),sub_bc(5,1),sub_bc(6,1))

                      frgf1=frgf1+1
                      call create_boundary(frgf1,lgf2,"i2",ncbd,0,&
                           sub_bc(1,2),sub_bc(2,2),sub_bc(3,2),sub_bc(4,2),sub_bc(5,2),sub_bc(6,2))
                   endif
                   if (j>1) then
                      lgf2=old2new_b(i,j-1,k,lgi1)
                      llf2=bg_to_bl(lgf2)
                      procf2=bg_to_proc(lgf2)

                      if (rank==procf1) then
                         sub_bc(1,1)=ii1(llf1)
                         sub_bc(2,1)=ii2(llf1)
                         sub_bc(3,1)=jj1(llf1)
                         sub_bc(4,1)=jj1(llf1)
                         sub_bc(5,1)=kk1(llf1)
                         sub_bc(6,1)=kk2(llf1)
                      endif
                      if (rank==procf2) then
                         sub_bc(1,2)=ii1(llf2)
                         sub_bc(2,2)=ii2(llf2)
                         sub_bc(3,2)=jj2(llf2)
                         sub_bc(4,2)=jj2(llf2)
                         sub_bc(5,2)=kk1(llf2)
                         sub_bc(6,2)=kk2(llf2)
                      endif

                      call bcast(sub_bc(:,1),procf1)
                      call bcast(sub_bc(:,2),procf2)

                      frgf1=frgf1+1
                      call create_boundary(frgf1,lgf1,"j1",ncbd,0,&
                           sub_bc(1,1),sub_bc(2,1),sub_bc(3,1),sub_bc(4,1),sub_bc(5,1),sub_bc(6,1))

                      frgf1=frgf1+1
                      call create_boundary(frgf1,lgf2,"j2",ncbd,0,&
                           sub_bc(1,2),sub_bc(2,2),sub_bc(3,2),sub_bc(4,2),sub_bc(5,2),sub_bc(6,2))
                   endif
                   if (k>1) then
                      lgf2=old2new_b(i,j,k-1,lgi1)
                      llf2=bg_to_bl(lgf2)
                      procf2=bg_to_proc(lgf2)

                      if (rank==procf1) then
                         sub_bc(1,1)=ii1(llf1)
                         sub_bc(2,1)=ii2(llf1)
                         sub_bc(3,1)=jj1(llf1)
                         sub_bc(4,1)=jj2(llf1)
                         sub_bc(5,1)=kk1(llf1)
                         sub_bc(6,1)=kk1(llf1)
                      endif
                      if (rank==procf2) then
                         sub_bc(1,2)=ii1(llf2)
                         sub_bc(2,2)=ii2(llf2)
                         sub_bc(3,2)=jj1(llf2)
                         sub_bc(4,2)=jj2(llf2)
                         sub_bc(5,2)=kk2(llf2)
                         sub_bc(6,2)=kk2(llf2)
                      endif

                      call bcast(sub_bc(:,1),procf1)
                      call bcast(sub_bc(:,2),procf2)

                      frgf1=frgf1+1
                      call create_boundary(frgf1,lgf1,"k1",ncbd,0,&
                           sub_bc(1,1),sub_bc(2,1),sub_bc(3,1),sub_bc(4,1),sub_bc(5,1),sub_bc(6,1))

                      frgf1=frgf1+1
                      call create_boundary(frgf1,lgf2,"k2",ncbd,0,&
                           sub_bc(1,2),sub_bc(2,2),sub_bc(3,2),sub_bc(4,2),sub_bc(5,2),sub_bc(6,2))
                   endif
                enddo
             enddo
          enddo
       enddo



       !############################################################################################
       !############################# SAVE NEW MESH FOR CHECKING PURPOSE ###########################
       !############################################################################################

       if(verbosity>=3) call write_mesh("test",x,y,z)

       deallocate(save_x,save_y,save_z,save_indfl,save_iminb,save_jminb,save_kminb)
       deallocate(save_imaxb,save_jmaxb,save_kmaxb,save_bcg_to_bcl,save_bcg_to_bci)
       deallocate(save_bg_to_proc,save_bg_to_bl,save_bg_to_bi,save_bl_to_bg,save_bcg_to_bg)
       deallocate(save_ii2,save_jj2,save_kk2,save_id1,save_jd1,save_kd1,save_id2,save_jd2)
       deallocate(save_kd2,save_npn)

       ip21=ndimntbx
       ip40=mdimubx            ! Nb point frontiere
       ip41=mdimtbx            ! Nb point frontiere
       ip42=mdimtbx            ! Nb point frontiere
       ip43=mdimtbx            ! Nb point frontiere
       ip44=0!mdimubx            !TODO
       test=.false.
       call reallocate(cl,mtb)
       call reallocate(ndcc,mtb)
       call reallocate(nbdc,mtb)
       call reallocate(nfbc,mtb)
       call reallocate(bc,mtb,ista*lsta)
       call reallocate(mdnc,mtt)
       call reallocate(mper,mtt)
       call reallocate(mpc,mtt)
       call reallocate(ncin,ip41)
       call reallocate(bceqt,ip41,neqt)
       call reallocate(mnc,ip43)

       !############################################################################################
       !################### INITIALIZE COINCIDENT BOUNDARIES #######################################
       !############################################################################################

    endif
    call reallocate(tmp,3)

       !    print*,'initialization '
       do frgf1=1,num_bcg
          frlf1=bcg_to_bcl(frgf1)
          lgf1=bcg_to_bg(frgf1)
          llf1=bg_to_bl(lgf1)
          frgi1=bcg_to_bci(frgf1)
          test=.false.
          procf1=bcg_to_proc(frgf1)
          if (frgi1==0) then  ! new boundary , it's easy to find the number of the associated bc
             test=.true.
             if (rank==procf1) then
                if (indfl(frlf1)(2:2)=="1") frgf2=frgf1+1
                if (indfl(frlf1)(2:2)=="2") frgf2=frgf1-1
             endif
             call bcast(frgf2,procf1)
          elseif(tab_raccord(frgi1)/=0) then ! old raccord boundary, it's more difficult
             test=.true.
             if (rank==procf1) &
                  call get_coords_box(sub_bc_c(1,1),sub_bc_c(1,2),sub_bc_c(1,3),sub_bc_c(1,4),sub_bc_c(1,5),sub_bc_c(1,6),  &
                  iminb(frlf1),imaxb(frlf1),jminb(frlf1),jmaxb(frlf1),kminb(frlf1),kmaxb(frlf1), &
                  id1(llf1),id2(llf1),jd1(llf1),jd2(llf1),kd1(llf1),kd2(llf1),npn(llf1),       &
                  x,y,z)
             call bcast(sub_bc_c,procf1)

             frgf4=0
             find_otherblock: do frlf3=1,mtb
                frgf3=bcl_to_bcg(frlf3)
                lgf3=bcg_to_bg(frgf3)
                llf3=bg_to_bl(lgf3)
                frgi3=bcg_to_bci(frgf3)
                if(frgi3/=0)then                       ! old boundary
                   if(tab_raccord(frgi3)==frgi1) then  ! potential new boundary number of the associated boundary

                      call get_coords_box(sub_bc_c(2,1),sub_bc_c(2,2),sub_bc_c(2,3),sub_bc_c(2,4),sub_bc_c(2,5),sub_bc_c(2,6),   &
                           iminb(frlf3),imaxb(frlf3),jminb(frlf3),jmaxb(frlf3),kminb(frlf3),kmaxb(frlf3), &
                           id1(llf3),id2(llf3),jd1(llf3),jd2(llf3),kd1(llf3),kd2(llf3),npn(llf3),       &
                           x,y,z)

                      if (abs(sub_bc_c(1,1)-sub_bc_c(2,1))<=eps .and. &
                           abs(sub_bc_c(1,2)-sub_bc_c(2,2))<=eps .and. &
                           abs(sub_bc_c(1,3)-sub_bc_c(2,3))<=eps .and. & ! It's me !
                           abs(sub_bc_c(1,4)-sub_bc_c(2,4))<=eps .and. &
                           abs(sub_bc_c(1,5)-sub_bc_c(2,5))<=eps .and. &
                           abs(sub_bc_c(1,6)-sub_bc_c(2,6))<=eps) then

                         frgf4=frgf3
                         exit find_otherblock
                      endif
                   endif
                endif
             enddo find_otherblock
             call sum_mpi(frgf4) ! there must be exactly one non zero value in this sum
             if (frgf4==0) then
                print*,"Coincident boundary not found ",frgi1,frgf1
                stop
             endif
             frgf2=frgf4
          endif
          if (test) then   ! raccord boundary

             lgf2=bcg_to_bg(frgf2)
             if (rank==bg_to_proc(lgf2)) then
                tmp(1)=iminb(bcg_to_bcl(frgf2))
                tmp(2)=jminb(bcg_to_bcl(frgf2))
                tmp(3)=kminb(bcg_to_bcl(frgf2))
             endif
             call bcast(tmp,bcg_to_proc(frgf2))

             if(verbosity>=2) then
                mot="" ; nmot=6   ; imot=0
                mot(1)="init"     ; imot(1)=4
                mot(2)="boundary" ; imot(2)=8
                mot(3)="basic"    ; imot(3)=5
                call str(mot,imot,nmx,4 ,frgf1)
                mot(5)="rc"       ; imot(5)=2
                call str(mot,imot,nmx,6 ,1)
                call c_inbdb( mot,imot,nmot,ncbd,ncin,bceqt,partition=.true.)
             else
                call inbdb( &
                     ncbd,ncin, &
                     frgf1,"rc  ",1, &
                     0,0,0,0,vbc,bceqt)
             endif

             if (rank==procf1) indmf=indfl(frlf1)
             call bcast(indmf,procf1)

             if(verbosity>=2) then
                nmot=12
                mot(1)="init"     ; imot(1)=4
                mot(2)="boundary" ; imot(2)=8
                mot(3)="coin"     ; imot(3)=4
                call str(mot,imot,nmx,4 ,frgf1)
                mot(5)="frc"      ; imot(5)=3
                call str(mot,imot,nmx,6 ,frgf2)
                mot(7)="kibdc"    ; imot(7)=5
                call str(mot,imot,nmx,8 ,1)
                mot(9)="krr"      ; imot(9)=3
                call str(mot,imot,nmx,10 ,1)
                mot(11)="epsmsh"  ; imot(11)=6
                write(mot(12),'(e11.3)') eps ;  imot(12)=11
                call c_inbdc(  mot,imot,nmot, exs1,exs2, x,y,z, ncbd,ncin,mnc)
             else
                call inbdc( &
                     exs1,exs2, &
                     x,y,z, &
                     ncbd,ncin,mnc, &
                     1,frgf1,frgf2,1,eps, &
                     tmp(1),tmp(2),tmp(3),"  ","  ","  ")
             endif
          endif
       enddo
      deallocate(tmp)
!      write(stderr,*) "done"
    return


  contains

    integer function old2new_b(i,j,k,lg)
      implicit none
      integer,intent(in) :: i,j,k,lg
      old2new_b=sum(nsplit(1:lg-1))+i+(j-1)*nsplit_dir(1,lg)+(k-1)*nsplit_dir(2,lg)*nsplit_dir(1,lg)
    end function old2new_b

    subroutine old2new_p(old_i,old_j,old_k,new_i,new_j,new_k,i,j,k,l)
      implicit none
      integer,intent(in) :: old_i,old_j,old_k,i,j,k,l
      integer,intent(out)  :: new_i,new_j,new_k
      new_i=old_i-sum(ni(:i-1,j,k,l))+i-1
      new_j=old_j-sum(nj(i,:j-1,k,l))+j-1
      new_k=old_k-sum(nk(i,j,:k-1,l))+k-1
    end subroutine old2new_p

    subroutine new2old_p(new_i,new_j,new_k,old_i,old_j,old_k,i,j,k,l)
      implicit none
      integer,intent(out) :: old_i,old_j,old_k
      integer,intent(in)  :: new_i,new_j,new_k,i,j,k,l
      old_i=new_i+sum(ni(:i-1,j,k,l))-i+1
      old_j=new_j+sum(nj(i,:j-1,k,l))-j+1
      old_k=new_k+sum(nk(i,j,:k-1,l))-k+1
    end subroutine new2old_p

  end subroutine partitionnement

  subroutine write_mesh(prefix,x,y,z)
    use boundary
    use mod_mpi
    use mod_vtk
    use para_var
    implicit none
    double precision,intent(in) :: x(:),y(:),z(:)
    character(*),intent(in) :: prefix

    integer :: l,fr,nid,njd,nijd,xyz,i,j,k,typ
    character(len=50)::fich,format
    character(len=4)::ext

!    typ=0 ! gnuplot format
!    typ=1 ! csv format
    typ=2 ! vtk format

    if(typ==0) format='(3e11.3,i8)'       
    if(typ==1) format='(3(e11.3,","),i8)'
    if(typ==2) format='(3e11.3)'

    if(typ==0) ext=".dat"
    if(typ==1) ext=".csv"
    if(typ==2) ext=".vts"

    ! write grid
    if(typ==2) then 
       call vtk_start_collection(prefix//"mesh.pvd","fields")
    endif

    do l=1,lt
       write(fich,'(A,I0.2,A)') prefix//"mesh_",bl_to_bg(l),ext

       if(typ==2) then 
        call vtk_open(fich,x,y,z,l)
        call vtk_writer(fich,bl_to_bg(l),"block_num",l)
        call vtk_close(fich)
       else
       open(42,file=fich,status="replace")

         if(typ==1)  write(42,*) "x,y,z,block"

       do k=kk1(l),kk2(l)
          do j=jj1(l),jj2(l)
             do i=ii1(l),ii2(l)

                nid = id2(l)-id1(l)+1
                njd = jd2(l)-jd1(l)+1
                nijd = nid*njd

                xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

                  write(42,format) x(xyz),y(xyz),z(xyz),bl_to_bg(l)

             enddo
               if (typ==0) write(42,*) ""
          enddo
       enddo
       close(42)
       endif
    enddo

    if(typ==2) then 
       call vtk_end_collection()
    endif
    
    call barrier

    if(typ==2) then 
       call vtk_start_collection(prefix//"bnd.pvd","fields")
    endif

    ! write boundaries
    do fr=1,mtb

       write(fich,'(A,I0.3,A)') prefix//"bnd_",bcl_to_bcg(fr),ext

       if(typ==2) then
        call vtk_open(fich,x,y,z,l)
        call vtk_writer(fich,bcl_to_bcg(fr),"boundary_num",ndlb(fr), &
                iminb(fr),imaxb(fr),jminb(fr),jmaxb(fr),kminb(fr),kmaxb(fr))
        call vtk_close(fich)
       else
       open(42,file=fich,status="replace")
         if(typ==1)  write(42,*) "x,y,z,boundary"

       do k=kminb(fr),kmaxb(fr)
          do j=jminb(fr),jmaxb(fr)
             do i=iminb(fr),imaxb(fr)
                l=ndlb(fr)
                nid = id2(l)-id1(l)+1
                njd = jd2(l)-jd1(l)+1
                nijd = nid*njd

                xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

                  write(42,format) x(xyz),y(xyz),z(xyz),bcl_to_bcg(fr)

             enddo
               if(typ==0.and.indfl(fr)(1:1)/="i") write(42,*) ""
          enddo
            if(typ==0.and.indfl(fr)(1:1)=="i") write(42,*) ""
       enddo
       close(42)
       endif
    enddo
    if(typ==2) then 
       call vtk_end_collection()
    endif
    call barrier
  end subroutine write_mesh

  subroutine compute_split(nsplit,nsplit_dir,nigf,njgf,nkgf,  &
       npoints,num_bgi,num_bgf,num_bli,   &
       nili,njli,nkli,bli_to_bgi,         &
       nigi,njgi,nkgi,bgi_to_proc)
    use mod_mpi
    use sortiefichier,only:stderr
    implicit none
    integer,intent(in)    :: npoints,num_bgi,num_bli
    integer,intent(inout) :: num_bgf,nsplit(:),nsplit_dir(:,:)
    integer,dimension(num_bli),intent(in) :: nili,njli,nkli,bli_to_bgi
    integer,dimension(num_bgi),intent(in) :: nigi,njgi,nkgi,bgi_to_proc
    integer,allocatable,intent(out) :: nigf(:,:,:,:),njgf(:,:,:,:),nkgf(:,:,:,:)

    integer,allocatable :: nilf(:,:,:),njlf(:,:,:),nklf(:,:,:)
    integer :: lg,ll,proc

    allocate(nilf(1,1,1),njlf(1,1,1),nklf(1,1,1))
    allocate(nigf(1,1,1,num_bgi),njgf(1,1,1,num_bgi),nkgf(1,1,1,num_bgi))

    ! routine calculant combiens de fois splitter chaque block initial
    ! sortie : nsplit
    ! entrée : tout le reste
    call num_split(nsplit,num_bgi,npoints,num_bgf,nigi,njgi,nkgi)

    do ll=1,num_bli
       lg=bli_to_bgi(ll)

       ! calcule le split pour un block, c'est à dire
       ! le nombre de blocs par direction              (nsplit_dir)
       ! le nombre de points dans chaque nouveau block (nilf,njlf,nklf)
       call triv_split(nsplit(lg),npoints,nili(ll) ,njli(ll) ,nkli(ll), &
            nsplit_dir(:,lg),nilf,njlf,nklf)

       ! ajoute le split avec les splits des autres blocks localement
       call reallocate_s(nigf,maxval(nsplit_dir(1,:)),maxval(nsplit_dir(2,:)),maxval(nsplit_dir(3,:)),num_bgi)
       nigf(:size(nilf,1),:size(nilf,2),:size(nilf,3),lg)=nilf

       call reallocate_s(njgf,maxval(nsplit_dir(1,:)),maxval(nsplit_dir(2,:)),maxval(nsplit_dir(3,:)),num_bgi)
       njgf(:size(njlf,1),:size(njlf,2),:size(njlf,3),lg)=njlf

       call reallocate_s(nkgf,maxval(nsplit_dir(1,:)),maxval(nsplit_dir(2,:)),maxval(nsplit_dir(3,:)),num_bgi)
       nkgf(:size(nklf,1),:size(nklf,2),:size(nklf,3),lg)=nklf
    enddo

    do lg=1,num_bgi
       ! propage à tout le monde les informations sur le découpage
       proc=bgi_to_proc(lg)
       call bcast(nsplit_dir(:,lg),proc)
    enddo
    do lg=1,num_bgi
       ! propage à tout le monde les informations sur le découpage
       proc=bgi_to_proc(lg)
       call reallocate_s(nigf,maxval(nsplit_dir(1,:)),maxval(nsplit_dir(2,:)),maxval(nsplit_dir(3,:)),num_bgi)
       call reallocate_s(njgf,maxval(nsplit_dir(1,:)),maxval(nsplit_dir(2,:)),maxval(nsplit_dir(3,:)),num_bgi)
       call reallocate_s(nkgf,maxval(nsplit_dir(1,:)),maxval(nsplit_dir(2,:)),maxval(nsplit_dir(3,:)),num_bgi)
       call bcast(nigf(:,:,:,lg),proc)
       call bcast(njgf(:,:,:,lg),proc)
       call bcast(nkgf(:,:,:,lg),proc)
    enddo

  end subroutine compute_split

  subroutine num_split(nblock2,lt,nxyza,nblocks,ii2,jj2,kk2)
    use mod_mpi,only:nprocs
    implicit none
    integer,intent(in)    :: lt,nxyza
    integer,intent(in)    :: ii2(lt),jj2(lt),kk2(lt)
    integer,intent(inout) :: nblocks
    integer,intent(out)   :: nblock2(lt)
    integer               :: rsize,sblock(lt),i,j,k,nblock(lt),li
    double precision      :: unbalance,unbalance1

    ! TODO
    ! switch to the alternative version which permit to have ideal blocks size
    ! need a criteria to avoid too small block, and need to manage the residual block
    ! todo
    sblock=ii2*jj2*kk2
    rsize=minval(sblock) ! smallest block
    nblocks=max(nblocks,nint(nxyza*1./rsize))   ! have all the block to be around the size of the smallest one
    nblocks=nblocks+mod(nblocks,nprocs)         ! have a multiple of the number of process
!nblocks=7
    !   compute number of spliting of each blocks with the best equilibrium
    ! do i=lt,nblocks-1                       ! split until lt>=nblocks
    !    unbalance=10000 ! a lot
    !    do j=1,lt
    !      nblock=nblock2                    ! try every split
    !      nblock(j)=nblock(j)+1
    !      sblock=ceiling(ii2*jj2*kk2*1./nblock) ! compute the current size of blocks (approx.)
    !      unbalance1=(maxval(sblock)-minval(sblock))*1./maxval(sblock) ! and unbalance
    !      if (unbalance1< unbalance) then ! if better remember it
    !        unbalance = unbalance1
    !        k=j
    !      endif
    !    end do
    !    nblock2(k)=nblock2(k)+1
    ! end do
    where(nblock2>0) nblock2=-nblock2 ! block user values
    where(nblock2==0) nblock2=1 
    li=sum(abs(nblock2))

    !   compute number of spliting of each blocks with the best equilibrium
    do i=li,nblocks-1                       ! split until we have nblocks blocks
       sblock=ceiling(ii2*jj2*kk2*1./nblock2) ! compute the current size of blocks
       j=maxloc(sblock,1)                        ! split the first bigest block
       if(sblock(j)<0) stop "partitionning impossible with user parameters"
       do k=j+1,lt
          if(sblock(k)==sblock(j) &          ! if more than one bigest block
               .and.nblock2(k)>nblock2(j)) &      ! split the most splitted
               j=k
       end do
       nblock2(j)=nblock2(j)+1
    end do
    where(nblock2<0) nblock2=-nblock2 ! restore user values

    !    compute number of spliting of each blocks with the ideal equilibrium
    !     rsize=nint(nxyza*1./nblocks)              ! ideal size of a block
    !     nblock2=ceiling(ii2*jj2*kk2*1./rsize) ! number of split needed
    !!     sblock=mod(ii2*jj2*kk2,rsize)         ! size of the smallest block
    !!     do i=1,lt
    !!       if (sblock(i) <= 3*3*3) &          ! allow for small imbalance in order to avoid too small blocks
    !!            nblock2(i)=nblock2(i)-1
    !!     end do
    ! do i=sum(nblock2),nblocks-1,-1
    !    sblock=ceiling(ii2*jj2*kk2*1./nblock2) ! compute the current size of blocks (approx.)
    !    unbalance=(maxval(sblock)-minval(sblock))*1./maxval(sblock) ! and unbalance
    !    k=0
    !    do j=1,lt
    !      nblock=nblock2                    ! try every split
    !      nblock(j)=nblock(j)-1
    !      sblock=ceiling(ii2*jj2*kk2*1./nblock) ! compute the current size of blocks (approx.)
    !      unbalance1=(maxval(sblock)-minval(sblock))*1./maxval(sblock) ! and unbalance
    !      if (unbalance1< unbalance) then ! if better remember it
    !        unbalance = unbalance1
    !        k=j
    !      endif
    !    end do
    !    print*,sum(nblock2),unbalance
    !    if (k==0) exit
    !    nblock2(k)=nblock2(k)-1
    ! end do

  end subroutine num_split

  subroutine triv_split(nblock2,nxyza,ii2,jj2,kk2, &
       nblockd,new_ii2,new_jj2,new_kk2)
    implicit none
    integer,intent(in)  :: ii2,jj2,kk2,nxyza,nblock2
    integer,allocatable,intent(out) :: new_ii2(:,:,:),new_jj2(:,:,:),new_kk2(:,:,:)
    integer,intent(inout)             :: nblockd(3)

    integer             :: i,j,k,i1,j1,k1,i2,j2,k2
    integer,allocatable :: tmp_ii2(:,:,:),tmp_jj2(:,:,:),tmp_kk2(:,:,:),num_cft(:,:,:), num_cf2(:,:,:)

    !   trivial spliting : divide my block in nblock2 subblock
    !                      test all possiblities constisting in dividing
    !                      i times in the x direction, j times in the y direction and k times in the z direction
    call reallocate(num_cf2,1,1,1)
    num_cf2=nxyza*nblock2*6 ! useless initial big value
    i1=nblockd(1)
    i2=nblockd(1)
    j1=nblockd(2)
    j2=nblockd(2)
    k1=nblockd(3)
    k2=nblockd(3)

    if(i1==0) i1=1
    if(i2==0) i2=nblock2
    if(j1==0) j1=1
    if(j2==0) j2=nblock2
    if(k1==0) k1=1
    if(k2==0) k2=nblock2


    do k=k1,k2
       do j=j1,j2
          do i=i1,i2
             if(i*j*k==nblock2) then !           if we get the right number of blocks
                allocate(num_cft(i,j,k))
                !       compute sizes of sub-blocks
                call compute_size(i,j,k,ii2,jj2,kk2,tmp_ii2,tmp_jj2,tmp_kk2)
!
                if (min(minval(tmp_ii2),minval(tmp_jj2),minval(tmp_kk2))>1) then ! if the splitting is acceptable
                   !             compute sizes of new communication, must be over evaluated (including boundary condition)
                   num_cft=2*(tmp_jj2*tmp_ii2 + tmp_ii2*tmp_kk2 + tmp_jj2*tmp_kk2)
                   !             choose the best splitting (less comm)
                   if(sum(num_cft)<sum(num_cf2)) then  !  TODO : is sum better than maxval ?
                      nblockd=[i,j,k]
                      call reallocate(  new_jj2,i,j,k) ;   new_jj2=tmp_jj2
                      call reallocate(  new_ii2,i,j,k) ;   new_ii2=tmp_ii2
                      call reallocate(  new_kk2,i,j,k) ;   new_kk2=tmp_kk2
                      call reallocate(num_cf2,i,j,k) ; num_cf2=num_cft
                   end if
                end if
                deallocate(tmp_jj2,tmp_ii2,tmp_kk2,num_cft)
             end if
          end do
       end do
    end do
  end subroutine triv_split

subroutine compute_size(i,j,k,ii2,jj2,kk2,new_ii2,new_jj2,new_kk2)
    !       compute sizes of sub-blocks
    implicit none
    integer,intent(in)  :: ii2,jj2,kk2,i,j,k
    integer,allocatable,intent(out) :: new_ii2(:,:,:),new_jj2(:,:,:),new_kk2(:,:,:)
    integer             :: i1,j1,k1

    allocate(new_ii2(i,j,k),new_jj2(i,j,k),new_kk2(i,j,k))
    do k1=1,k
       do j1=1,j
          do i1=1,i
             new_ii2(i1,j1,k1)=nint(i1*(ii2+i-1)*1./i) - nint((i1-1.)*(ii2+i-1)*1./i)
             new_jj2(i1,j1,k1)=nint(j1*(jj2+j-1)*1./j) - nint((j1-1.)*(jj2+j-1)*1./j)    ! count interface twice
             new_kk2(i1,j1,k1)=nint(k1*(kk2+k-1)*1./k) - nint((k1-1.)*(kk2+k-1)*1./k)
          end do
       end do
    end do
  end subroutine compute_size

  subroutine create_grid(l,ni,nj,nk,bi,proc)
    use para_fige,only:nmx
    use mod_c_crdms
    use mod_crdms
    implicit none
    integer,intent(in) :: l,ni,nj,nk,bi,proc
    integer            :: imot(nmx),nmot
    character(len=32)  :: mot(nmx)

    if(verbosity>=2) then
       mot="" ; nmot=7 ; imot=0
       mot(1)="create" ; imot(1)=6
       mot(2)="dom"    ; imot(2)=3
       mot(3)="st"     ; imot(3)=2
       call str(mot,imot,nmx,4 ,l)
       call str(mot,imot,nmx,5 ,ni)
       call str(mot,imot,nmx,6 ,nj)
       call str(mot,imot,nmx,7 ,nk)
       call c_crdms( mot,imot,nmot,bi,proc)
    else
       call crdms(l,ni,nj,nk,bi,proc)
    endif
  end subroutine create_grid

  subroutine create_boundary(fr,l,indmf,ncbd,bci,&
       imin,imax,jmin,jmax,kmin,kmax)
    use para_fige,only:nmx
    use mod_c_crbds
    use mod_crbds
    implicit none
    integer            ,intent(in)    :: fr,l,bci,imin,imax,jmin,jmax,kmin,kmax
    integer,allocatable,intent(inout) :: ncbd(:)
    character(len=2)   ,intent(in)    :: indmf
    integer            :: imot(nmx),nmot
    character(len=32)  :: mot(nmx)

    if(verbosity>=2) then
       mot="" ; nmot=13  ; imot=0
       mot(1)="create"   ; imot(1)=6
       mot(2)="boundary" ; imot(2)=8
       mot(3)="st"       ; imot(3)=2
       call str(mot,imot,nmx,4 ,fr)
       call str(mot,imot,nmx,5 ,1)
       call str(mot,imot,nmx,6 ,l)
       call str(mot,imot,nmx,7 ,imin)
       call str(mot,imot,nmx,8 ,imax)
       call str(mot,imot,nmx,9 ,jmin)
       call str(mot,imot,nmx,10,jmax)
       call str(mot,imot,nmx,11,kmin)
       call str(mot,imot,nmx,12,kmax)
       mot(13)=indmf ; imot(13)=2
       call c_crbds( mot,imot,nmot, ncbd,bci)
    else
       call crbds(fr,1,l, &
            imin,imax,jmin,jmax,kmin,kmax, &
            indmf,ncbd,bci)
    endif
  end subroutine create_boundary

  subroutine copy_grid(l,l2,xs,ys,zs,x,y,z,save_x,save_y,save_z,       &
       save_id1,save_id2,save_jd1,save_jd2,save_kd1,save_kd2,save_npn)
    use maillage,only : kk1,kk2,jj1,jj2,ii1,ii2,kd1,kd2,jd1,jd2,id1,id2,npn
    implicit none
    integer,intent(in)             :: l,l2,xs,ys,zs,save_npn(:)
    integer,intent(in)             :: save_id1(:),save_id2(:),save_jd1(:),save_jd2(:),save_kd1(:),save_kd2(:)
    double precision,intent(in)    :: save_x(:),save_y(:),save_z(:)
    double precision,intent(inout) :: x(:),y(:),z(:)

    integer :: xi,yi,zi,nid,njd,nijd,xyz
    integer :: save_xi,save_yi,save_zi,save_nid,save_njd,save_nijd,save_xyz

    nid = id2(l2)-id1(l2)+1
    njd = jd2(l2)-jd1(l2)+1
    nijd = nid*njd

    save_nid = save_id2(l)-save_id1(l)+1
    save_njd = save_jd2(l)-save_jd1(l)+1
    save_nijd = save_nid*save_njd

    do zi=kk1(l2),kk2(l2)
       do yi=jj1(l2),jj2(l2)
          do xi=ii1(l2),ii2(l2)

             save_xi=xi+xs
             save_yi=yi+ys
             save_zi=zi+zs

             xyz      =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
             save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd

             ! fill grid
             x(xyz)=save_x(save_xyz)
             y(xyz)=save_y(save_xyz)
             z(xyz)=save_z(save_xyz)
          enddo
       enddo
    enddo
  end subroutine copy_grid


  subroutine send_grid(x,y,z,xs,ys,zs,xe,ye,ze,       &
       id1,jd1,kd1,nid,nijd,npn,orig,dest)
    use mod_mpi
    implicit none
    integer,intent(in)             :: xs,ys,zs,xe,ye,ze,id1,jd1,kd1,nid,nijd,npn,orig,dest
    double precision,intent(in)    :: x(:),y(:),z(:)
    double precision,allocatable   :: buff(:,:,:,:)
    integer :: xi,yi,zi,xyz

    allocate(buff(3,xe-xs+1,ye-ys+1,ze-zs+1))
    do zi=zs,ze
       do yi=ys,ye
          do xi=xs,xe

             xyz =npn+1+(xi-id1)+(yi-jd1)*nid+(zi-kd1)*nijd

             ! fill buffer
             buff(1,xi-xs+1,yi-ys+1,zi-zs+1)=x(xyz)
             buff(2,xi-xs+1,yi-ys+1,zi-zs+1)=y(xyz)
             buff(3,xi-xs+1,yi-ys+1,zi-zs+1)=z(xyz)
          enddo
       enddo
    enddo

    call MPI_TRANS(buff,buff,orig,DEST)

    deallocate(buff)
  end subroutine send_grid

  subroutine recv_grid(x,y,z,l,orig,dest)
    use mod_mpi
    use maillage
    implicit none
    integer                     ,intent(in)    :: l,orig,dest
    double precision,allocatable,intent(inout) :: x(:),y(:),z(:)

    double precision,allocatable :: buff(:,:,:,:)
    integer                      :: xi,yi,zi,xyz,nid,njd,nijd

    allocate(buff(3,ii2(l)-ii1(l)+1,jj2(l)-jj1(l)+1,kk2(l)-kk1(l)+1))

    call reallocate_s(x,ndimntbx)
    call reallocate_s(y,ndimntbx)
    call reallocate_s(z,ndimntbx)

    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd = nid*njd

    call MPI_TRANS(buff,buff,orig,dest)

    do zi=kk1(l),kk2(l)
       do yi=jj1(l),jj2(l)
          do xi=ii1(l),ii2(l)

             xyz = npn(l)+1+(xi-id1(l)) + (yi-jd1(l))*nid + (zi-kd1(l))*nijd

             ! fill grid
             x(xyz)=buff(1,xi-ii1(l)+1,yi-jj1(l)+1,zi-kk1(l)+1)
             y(xyz)=buff(2,xi-ii1(l)+1,yi-jj1(l)+1,zi-kk1(l)+1)
             z(xyz)=buff(3,xi-ii1(l)+1,yi-jj1(l)+1,zi-kk1(l)+1)
          enddo
       enddo
    enddo
    deallocate(buff)
  end subroutine recv_grid

  subroutine get_coords_box(xmin,xmax,ymin,ymax,zmin,zmax,  &
       a,b,c,d,e,f,id1,id2,jd1,jd2,kd1,kd2,npn,x,y,z)
    implicit none
    double precision,intent(out) :: xmin,xmax,ymin,ymax,zmin,zmax
    integer,intent(in)           :: a,b,c,d,e,f,id1,id2,jd1,jd2,kd1,kd2,npn
    double precision,intent(in)  :: x(:),y(:),z(:)

    integer :: nid,njd,nijd,p1,p2,p3,p

    nid = id2-id1+1
    njd = jd2-jd1+1
    nijd = nid*njd

    xmin= Huge(1.d0)
    xmax=-Huge(1.d0)
    ymin= Huge(1.d0)
    ymax=-Huge(1.d0)
    zmin= Huge(1.d0)
    zmax=-Huge(1.d0)

    do p1=a,b,max(b-a,1)
        do p2=c,d,max(d-c,1)
            do p3=e,f,max(f-e,1)
              p=npn+1+(p1 - id1)+(p2 - jd1)*nid+(p3 - kd1)*nijd
              xmin=min(xmin,x(p))
              xmax=max(xmax,x(p))
              ymin=min(ymin,y(p))
              ymax=max(ymax,y(p))
              zmin=min(zmin,z(p))
              zmax=max(zmax,z(p))
            enddo
        enddo
    enddo

  end subroutine get_coords_box

  subroutine iniraccord(mot,imot,nmot)
    use chainecarac,only : ci
    use para_fige,only : nmx
    use boundary,only : tab_raccord
    use mod_valenti
    implicit none
    character(len=32) ::  mot(nmx),comment
    integer             :: imot(nmx),nmot,fr1,fr2,nm,kval,icmt

    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
    kval=0
    !
    nm=3
    if(nmot.lt.nm) then ! read number first index of coincident boundary
       comment=ci
       call synterr(mot,imot,nmot,comment)
    else
       call valenti(mot,imot,nm,fr1,kval)
    endif
    nm=nm+1
    if(nmot.lt.nm) then ! read number second index of coincident boundary
       comment=ci
       call synterr(mot,imot,nmot,comment)
    else
       call valenti(mot,imot,nm,fr2,kval)
    endif

    ! fr1 and fr2 are coincident boundary
    tab_raccord(fr1)=fr2
    tab_raccord(fr2)=fr1
  end subroutine iniraccord

  subroutine str(mot,imot,nmx,lmot,val)
    implicit none
    integer          ,intent(in)    :: nmx,lmot,val
    integer          ,intent(inout) :: imot(nmx)
    character(len=32),intent(inout) ::  mot(nmx)

    write(mot(lmot) ,*) val
    mot(lmot)=adjustl(mot(lmot))
    imot(lmot) =len_trim(adjustl(mot(lmot)))
  end subroutine str

end module mod_partitionnement
