!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file SourcePosition.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<
!! Locate physical sources, then expand CG moments into owned element entries.
subroutine SourcePosition (Tdomain)
    use sdomain
    use ssources
    use mpi
    use mshape8
    use mshape27
    implicit none
    type(domain), intent(inout) :: Tdomain
    type(source), allocatable :: local_sources(:), candidates(:)
    type(source) :: item
    real(fpp), allocatable :: coord(:,:)
    real(fpp) :: point(0:2), ref(0:2), jac(0:2,0:2), invjac(0:2,0:2)
    real(fpp) :: eps_x, scale, det, local_weight, total_weight
    real(fpp), parameter :: EPS_REF=1e-10_fpp
    integer :: ns, el, j, axis, dom, nlocal, ncand, i, rg, ierr, owner, proposed
    integer :: local_counts(2), global_counts(2)
    logical :: ok, fluid_source, solid_element, inside

    rg = Tdomain%rank
    nlocal = 0
    allocate(local_sources(0:max(1,Tdomain%n_source)-1), candidates(0:7))
    allocate(coord(0:2,0:Tdomain%n_nodes-1))
    do ns=0,Tdomain%n_source-1
        point = [Tdomain%sSource(ns)%Xsource,Tdomain%sSource(ns)%Ysource,Tdomain%sSource(ns)%Zsource]
        fluid_source = Tdomain%sSource(ns)%i_type_source==3 .or. Tdomain%sSource(ns)%i_type_source==7
        ncand = 0
        ! Scan owned elements, with a tolerant physical bounding-box prefilter.
        ! Unlike the old fixed candidate buffer this has no valence-dependent cap.
        do el=0,Tdomain%n_elem-1
            dom = Tdomain%specel(el)%domain
            solid_element = dom==DM_SOLID_CG .or. dom==DM_SOLID_DG .or. dom==DM_SOLID_CG_PML
            if (fluid_source .eqv. solid_element) cycle
            do j=0,Tdomain%n_nodes-1
                coord(:,j) = Tdomain%Coord_Nodes(:,Tdomain%specel(el)%Control_Nodes(j))
            enddo
            scale = maxval(maxval(coord,dim=2)-minval(coord,dim=2))
            eps_x = max(EPS_REF*scale,64*epsilon(1._fpp)*max(1._fpp,maxval(abs(coord))))
            inside = .true.
            do axis=0,2
                if (point(axis)<minval(coord(axis,:))-eps_x .or. &
                    point(axis)>maxval(coord(axis,:))+eps_x) inside=.false.
            enddo
            if (.not.inside) cycle
            if (Tdomain%n_nodes==8) then
                call shape8_global2local(coord,point(0),point(1),point(2),ref(0),ref(1),ref(2),ok)
            else
                call shape27_global2local(coord,point(0),point(1),point(2),ref(0),ref(1),ref(2),ok)
            endif
            if (.not.ok) cycle
            if (any(abs(ref)>1._fpp+EPS_REF)) cycle
            do axis=0,2
                if (abs(abs(ref(axis))-1._fpp)<=EPS_REF) ref(axis)=sign(1._fpp,ref(axis))
            enddo
            if (Tdomain%n_nodes==8) then
                call shape8_local2jacob(coord,ref(0),ref(1),ref(2),jac)
            else
                call shape27_local2jacob(coord,ref(0),ref(1),ref(2),jac)
            endif
            det = dot_product(jac(0,:),cross(jac(1,:),jac(2,:)))
            if (det<=0._fpp) then
                write(*,*) 'ERROR: nonpositive source-element Jacobian',ns,el,det
                call MPI_Abort(Tdomain%communicateur,1,ierr)
                stop 1
            endif
            invjac=jac
            call invert_3d(invjac,det)
            item=Tdomain%sSource(ns)
            item%elem=el
            item%proc=rg
            item%refcoord=ref
            item%InvGrad=invjac
            item%spatial_weight=1._fpp
            if (item%i_type_source==2 .and. dom==DM_SOLID_CG) &
                item%spatial_weight=cone_fraction(jac,invjac,ref)
            call append(candidates,ncand,item)
        enddo

        local_counts=0
        do i=0,ncand-1
            if (Tdomain%specel(candidates(i)%elem)%domain==DM_SOLID_CG) then
                local_counts(1)=local_counts(1)+1
            else
                local_counts(2)=local_counts(2)+1
            endif
        enddo
        call MPI_Allreduce(local_counts,global_counts,2,MPI_INTEGER,MPI_SUM,Tdomain%communicateur,ierr)
        if (sum(global_counts)==0) then
            if (rg==0) write(*,*) 'ERROR: source not found in a compatible element',ns,point
            call MPI_Abort(Tdomain%communicateur,1,ierr)
            stop 1
        endif

        if (Tdomain%sSource(ns)%i_type_source==2 .and. global_counts(1)>0) then
            if (global_counts(2)>0) then
                if (rg==0) write(*,*) 'ERROR: CG moment on mixed CG/DG/PML interface is unsupported',ns
                call MPI_Abort(Tdomain%communicateur,1,ierr)
                stop 1
            endif
            local_weight=0._fpp
            do i=0,ncand-1
                local_weight=local_weight+candidates(i)%spatial_weight
            enddo
            call MPI_Allreduce(local_weight,total_weight,1,MPI_DOUBLE_PRECISION,MPI_SUM,Tdomain%communicateur,ierr)
            if (total_weight<=0._fpp) then
                call MPI_Abort(Tdomain%communicateur,1,ierr)
                stop 1
            endif
            do i=0,ncand-1
                item=candidates(i)
                item%spatial_weight=item%spatial_weight/total_weight
                call append(local_sources,nlocal,item)
                if (Tdomain%config%verbose_level>1) &
                    write(*,*) 'Source contribution',ns,'rank',rg,'elem',item%elem,'ref',item%refcoord, &
                               'weight',item%spatial_weight
            enddo
            if (rg==0 .and. Tdomain%config%verbose_level>1) &
                write(*,*) 'CG moment incident elements',global_counts(1),'raw cone fraction',total_weight
        else
            ! Point forces and acoustic monopoles retain exactly one global owner.
            ! DG retains its existing one-sided source convention.
            proposed=-1
            if (ncand>0) proposed=rg
            call MPI_Allreduce(proposed,owner,1,MPI_INTEGER,MPI_MAX,Tdomain%communicateur,ierr)
            if (rg==owner) then
                item=candidates(0)
                item%spatial_weight=1._fpp
                call append(local_sources,nlocal,item)
                if (Tdomain%config%verbose_level>1) then
                    write(*,*) 'Found source',ns,'on proc',owner
                    write(*,*) 'Source local coordinates:',item%refcoord
                endif
            endif
        endif
    enddo
    deallocate(Tdomain%sSource)
    Tdomain%n_source=nlocal
    Tdomain%logicD%any_source=nlocal>0
    if (nlocal>0) then
        allocate(Tdomain%sSource(0:nlocal-1))
        Tdomain%sSource=local_sources(0:nlocal-1)
    endif
contains
    subroutine append(buffer,n,value)
        type(source), allocatable, intent(inout) :: buffer(:)
        integer, intent(inout) :: n
        type(source), intent(in) :: value
        type(source), allocatable :: grown(:)
        if (n==size(buffer)) then
            allocate(grown(0:2*size(buffer)-1))
            grown(0:n-1)=buffer
            call move_alloc(grown,buffer)
        endif
        buffer(n)=value
        n=n+1
    end subroutine append

    function cross(a,b) result(c)
        real(fpp), intent(in) :: a(0:2),b(0:2)
        real(fpp) :: c(0:2)
        c=[a(1)*b(2)-a(2)*b(1),a(2)*b(0)-a(0)*b(2),a(0)*b(1)-a(1)*b(0)]
    end function cross

    function cone_fraction(jac,invjac,ref) result(w)
        real(fpp), intent(in) :: jac(0:2,0:2),invjac(0:2,0:2),ref(0:2)
        real(fpp) :: w,u(0:2),v(0:2),z(0:2),nu,nv,nz,denom,cosang
        integer :: active(3),n,a
        n=0
        do a=0,2
            if (abs(ref(a))==1._fpp) then
                n=n+1
                active(n)=a
            endif
        enddo
        select case(n)
        case(0)
            w=1._fpp
        case(1)
            w=0.5_fpp
        case(2)
            ! jac rows are physical tangent rays; inverse columns are normals.
            u=-sign(1._fpp,ref(active(1)))*invjac(:,active(1))
            v=-sign(1._fpp,ref(active(2)))*invjac(:,active(2))
            cosang=dot_product(u,v)/sqrt(dot_product(u,u)*dot_product(v,v))
            w=(M_PI-acos(max(-1._fpp,min(1._fpp,cosang))))/(2*M_PI)
        case(3)
            u=-sign(1._fpp,ref(0))*jac(0,:)
            v=-sign(1._fpp,ref(1))*jac(1,:)
            z=-sign(1._fpp,ref(2))*jac(2,:)
            nu=sqrt(dot_product(u,u)); nv=sqrt(dot_product(v,v)); nz=sqrt(dot_product(z,z))
            denom=nu*nv*nz+dot_product(u,v)*nz+dot_product(v,z)*nu+dot_product(z,u)*nv
            w=2*atan2(abs(dot_product(u,cross(v,z))),denom)/(4*M_PI)
        end select
    end function cone_fraction
end subroutine SourcePosition

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
