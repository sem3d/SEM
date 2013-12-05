program create_model_2D
    !- creation of a 2D mesh; definition of some elastic quantities
    !- data read in files different from the usual ones in the "spectral element"-like
    !     code
    use gid2spec
    implicit none

    integer                        :: n_elem_x,n_elem_z,ngllx,ngllz,n_mat,     &
        model_type,i,j,k,nelem,nnodes,nn, flag,  &
        p_ngllx,p_ngllz,p_up,logic_pml,nm, npow, &
        n_elem_xinit,n_elem_zinit
    integer, parameter             :: RANDOM = 1, NOEUD = 4, n_dime = 2
    real                           :: Vp,Vs,rho,x_zone,z_zone, fact_deep,   &
        Apow,omegac,kc,sm_len,dtmin,x_start,z_start,x_step,z_step
    real, parameter                :: xzmin = 0d0
    real, dimension(:,:), allocatable :: coord_nodes
    integer, dimension(:,:), allocatable           :: specel
    logical                        :: p_filter,log_created
    logical, parameter             :: VRAI = .true. , FAUX = .false.
    character(len=50)              :: file_out,mater_out

    n_mat = 1
    logic_pml = 0
    p_up = 0

    write(*,*) "  -> Mesh creation - 2D"
    !- file reading
    open(11,file="input.create",action="read",status="old")
    read(11,*) log_created   ! if .true. : geometrical file already exists, not to be created
    if(log_created)then
        read(11,*) file_out
        close(11)
    else   ! info to create a model
        read(11,*) file_out
        read(11,*) x_zone,z_zone
        read(11,*) n_elem_x,n_elem_z
        read(11,*) ngllx,ngllz
        read(11,*) Vp,Vs,rho
        read(11,*) dtmin
        read(11,*) logic_pml
        if(logic_pml == 1)then
            read(11,*) p_up
            read(11,*) fact_deep
            read(11,*) p_ngllx,p_ngllz
            read(11,*) p_filter,npow,apow,omegac,kc
        else
            read(11,*) ; read(11,*) ; read(11,*) ; read(11,*)
            fact_deep = 0.0
            p_ngllx = 0
            p_ngllz = 0
            p_filter = 0.0
            npow=0
            apow=0.0
            omegac=0.0
            kc=0.0
        end if
        read(11,*) mater_out
        close(11)

        n_elem_xinit = n_elem_x
        n_elem_zinit = n_elem_z

        x_step = x_zone/n_elem_x
        z_step = z_zone/n_elem_z

        if(logic_pml == 1)then   ! changing in the starting x and z
            x_start = -1d0*x_zone/n_elem_x
            z_start = -1d0*z_zone/n_elem_z
        else
            x_start = 0.0
            z_start = 0.0
        end if
        n_elem_x = n_elem_x + 2*logic_pml
        n_elem_z = n_elem_z + 2*logic_pml
        if(logic_pml == 1 .and. p_up == 0) n_elem_z = n_elem_z - 1

        !- nodes
        nnodes = -1
        allocate(coord_nodes(0:1,0:(n_elem_x+1)*(n_elem_z+1)-1))
        do j = 0, n_elem_z
            do i = 0, n_elem_x
                nnodes = nnodes+1
                Coord_Nodes(0,nnodes) = x_start + i*x_step
                Coord_Nodes(1,nnodes) = z_start + j*z_step
                ! change in PMLs boundary co-ordinates
                if(logic_pml == 1 .and. i == 0) Coord_Nodes(0,nnodes) = Coord_Nodes(0,nnodes)-fact_deep*x_step
                if(logic_pml == 1 .and. i == n_elem_x) Coord_Nodes(0,nnodes) = Coord_Nodes(0,nnodes)+fact_deep*x_step
                if(logic_pml == 1 .and. j == 0) Coord_Nodes(1,nnodes) = Coord_Nodes(1,nnodes)-fact_deep*z_step
                if(logic_pml == 1 .and. p_up == 1 .and. j == n_elem_z) Coord_Nodes(1,nnodes) = Coord_Nodes(1,nnodes)+fact_deep*z_step

            end do
        end do

        !- elements
        nelem = -1
        allocate(specel(0:3,0:n_elem_x*n_elem_z-1))
        do j = 0, n_elem_z-1
            do i = 0, n_elem_x-1
                nelem = nelem + 1
                specel(0,nelem) = 1 + i + j*(n_elem_x+1)
                specel(1,nelem) = 2 + i + j*(n_elem_x+1)
                specel(3,nelem) = 2 + i + n_elem_x + j*(n_elem_x+1)
                specel(2,nelem) = 3 + i + n_elem_x + j*(n_elem_x+1)
            end do
        end do


        !- temporary file creation
        open(10,file="temp.txt",action="write",status="replace")
        write(10,*) "MESH: tmp"
        write(10,*) n_dime,"    Dimension"
        write(10,*) "Quadrilateral  ","   ElemType"
        write(10,*) NOEUD, "  Nnode"
        write(10,*) (n_elem_x+1)*(n_elem_z+1), "  Number of points"
        write(10,*) n_elem_x*n_elem_z, "  Number of elements"
        write(10,*) n_mat + 5*logic_pml + 3*p_up,"   Number of materials"
        write(10,*)
        write(10,*) "Coordinates"
        do i = 0, (n_elem_x+1)*(n_elem_z+1)-1
            write(10,*) i+1, Coord_Nodes(0,i),Coord_Nodes(1,i)
        end do
        write(10,*) "end coordinates"
        write(10,*)
        write(10,*) "Elements"
        nelem = -1
        do j = 0, n_elem_z -1
            do i = 0, n_elem_x -1
                nelem = nelem + 1
                !-  PMLs -------
                if(logic_pml == 1 .and. j == 0 .and. i == 0)then
                    nm = n_mat+1
                    write(10,*) nelem+1,(specel(k,nelem),k=0,3),nm
                    cycle
                end if
                if(logic_pml == 1 .and. j == 0 .and. i == n_elem_x-1)then
                    nm = n_mat+2
                    write(10,*) nelem+1,(specel(k,nelem),k=0,3),nm
                    cycle
                end if
                if(logic_pml == 1 .and. j == n_elem_z-1 .and. i == 0 .and. p_up == 1)then
                    nm = n_mat+3
                    write(10,*) nelem+1,(specel(k,nelem),k=0,3),nm
                    cycle
                end if
                if(logic_pml == 1 .and. j == n_elem_z-1 .and. i == n_elem_x-1 .and. p_up == 1)then
                    nm = n_mat+4
                    write(10,*) nelem+1,(specel(k,nelem),k=0,3),nm
                    cycle
                end if
                if(logic_pml == 1 .and. j == 0)then
                    nm = n_mat+3+2*p_up
                    write(10,*) nelem+1,(specel(k,nelem),k=0,3),nm
                    cycle
                end if
                if(logic_pml == 1 .and. i == 0)then
                    nm = n_mat+4+2*p_up
                    write(10,*) nelem+1,(specel(k,nelem),k=0,3),nm
                    cycle
                end if
                if(logic_pml == 1 .and. i == n_elem_x-1)then
                    nm = n_mat+5+2*p_up
                    write(10,*) nelem+1,(specel(k,nelem),k=0,3),nm
                    cycle
                end if
                if(logic_pml == 1 .and. j == n_elem_z-1 .and. p_up == 1)then
                    nm = n_mat+6+2*p_up
                    write(10,*) nelem+1,(specel(k,nelem),k=0,3),nm
                    cycle
                end if

                write(10,*) nelem+1,(specel(k,nelem),k=0,3),n_mat

                !--------
            end do
        end do

        write(10,*) "end elements"
        close(10)

        deallocate(specel,coord_nodes)


        !- material file
        open(10,file=mater_out,status="replace",action="write")
        write(10,*) "# File in which material properties are stored"
        write(10,*) "# For any material, put type, Vp, Vs, rho, dt, ngllx, ngllz"
        write(10,"(a1,2x,2(f0.2,2x),f5.0,2x,f0.9,2x,i2,2x,i2)") "S",Vp, Vs,rho,  &
            dtmin,ngllx,ngllz
        if(logic_pml == 1)then
            write(10,"(a1,2x,2(f0.2,2x),f5.0,2x,f0.9,2x,i2,2x,i2)") "P",Vp,Vs,rho,  &
                dtmin,p_ngllx,p_ngllz
            write(10,"(a1,2x,2(f0.2,2x),f5.0,2x,f0.9,2x,i2,2x,i2)") "P",Vp,Vs,rho,  &
                dtmin,p_ngllx,p_ngllz
            if(p_up == 1)then   !- PML on top
                write(10,"(a1,2x,2(f0.2,2x),f5.0,2x,f0.9,2x,i2,2x,i2)") "P",Vp,Vs,rho,  &
                    dtmin,p_ngllx,p_ngllz
                write(10,"(a1,2x,2(f0.2,2x),f5.0,2x,f0.9,2x,i2,2x,i2)") "P",Vp,Vs,rho,  &
                    dtmin,p_ngllx,p_ngllz
            end if
            write(10,"(a1,2x,2(f0.2,2x),f5.0,2x,f0.9,2x,i2,2x,i2)") "P",Vp,Vs,rho,  &
                dtmin,ngllx,p_ngllz
            write(10,"(a1,2x,2(f0.2,2x),f5.0,2x,f0.9,2x,i2,2x,i2)") "P",Vp,Vs,rho,  &
                dtmin,p_ngllx,ngllz
            write(10,"(a1,2x,2(f0.2,2x),f5.0,2x,f0.9,2x,i2,2x,i2)") "P",Vp,Vs,rho,  &
                dtmin,p_ngllx,ngllz
            if(p_up == 1)then   !- PML on top
                write(10,"(a1,2x,2(f0.2,2x),f5.0,2x,f0.9,2x,i2,2x,i2)") "P",Vp,Vs,rho,  &
                    dtmin,ngllx,p_ngllz
            end if
            write(10,*) "# Specifications for PMLs"
            write(10,*) "# Filtering, npow, Apow, Px, Left, Pz, Down, omegac, kc"
            write(10,"(l1,2x,i2,2x,f0.1,4(2x,l1),2x,f0.1,2x,f0.1)") p_filter, npow, apow,   &
                VRAI,VRAI,VRAI,VRAI,omegac,kc
            write(10,"(l1,2x,i2,2x,f0.1,4(2x,l1),2x,f0.1,2x,f0.1)") p_filter, npow, apow,   &
                VRAI,FAUX,VRAI,VRAI,omegac,kc
            if(p_up == 1)then ! PML on top
                write(10,"(l1,2x,i2,2x,f0.1,4(2x,l1),2x,f0.1,2x,f0.1)") p_filter, npow, apow,   &
                    VRAI,VRAI,VRAI,FAUX,omegac,kc
                write(10,"(l1,2x,i2,2x,f0.1,4(2x,l1),2x,f0.1,2x,f0.1)") p_filter, npow, apow,   &
                    VRAI,FAUX,VRAI,FAUX,omegac,kc
            end if
            write(10,"(l1,2x,i2,2x,f0.1,4(2x,l1),2x,f0.1,2x,f0.1)") p_filter, npow, apow,   &
                FAUX,FAUX,VRAI,VRAI,omegac,kc
            write(10,"(l1,2x,i2,2x,f0.1,4(2x,l1),2x,f0.1,2x,f0.1)") p_filter, npow, apow,   &
                VRAI,VRAI,FAUX,FAUX,omegac,kc
            write(10,"(l1,2x,i2,2x,f0.1,4(2x,l1),2x,f0.1,2x,f0.1)") p_filter, npow, apow,   &
                VRAI,FAUX,FAUX,FAUX,omegac,kc
            if(p_up == 1)then ! PML on top
                write(10,"(l1,2x,i2,2x,f0.1,4(2x,l1),2x,f0.1,2x,f0.1)") p_filter, npow, apow,   &
                    FAUX,FAUX,VRAI,FAUX,omegac,kc
            end if
        end if
        close(10)

    end if
    !---------------------
    call gid2spec_seq(file_out,log_created)
    write(*,*) "     Mesh: done"

    !------
end program create_model_2D
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

