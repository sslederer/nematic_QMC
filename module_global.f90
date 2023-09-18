module global
    implicit none

    integer, parameter :: LEFT=-1, RIGHT=1, CENTER=0
    integer, parameter :: HORIZONTAL=-1, VERTICAL=1, INV=-1, NO_INV=1
    real(8) PI
    integer, parameter :: $is_complex = 1$ ! ##@@##NEVER CHANGE ANYTHING IN THIS LINE. NOT EVEN SPACES.
    integer, parameter :: $is_local_int = 0$ ! ##@@##NEVER CHANGE ANYTHING IN THIS LINE. NOT EVEN SPACES.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!! Site indexing scheme: ind = is + 2*(ip-1) + 4*(ix-1) + 4*Lx*(iy-1) !!!!!
    !!!!!!!!                            spin, orbital, x, y coordinates
    !!!!!!!! Type Definitions !!!!!!!!!

    !!!! Complex matrix
    type CMatrix
        integer(4) d
        integer:: is_complex = 1
        complex(8), allocatable :: m(:,:)
    end type

	!!!! Real matrix
    type Rmatrix
		integer(4) d
        integer:: is_complex = 0
		real(8), allocatable :: m(:,:)
	end type

    !!!! Complex matrix, non-square
    type CMatrix_nonsq
        integer(4) d1, d2
        complex(8), allocatable :: m(:,:)
    end type

    !!!! Kinetic energy matrix
    type KinMat ! stored as a sparse matrix
        integer(4)  num_pairs ! number of bonds
        integer(4), allocatable :: ind(:,:) !numpairs x2: stores indices of sites connected by each bond
        $thop_type$, allocatable :: thop(:) ! numpairs
        integer(4):: dir, shift
        integer(4), allocatable :: pair_ind(:,:)
    END TYPE

    !!!! order parameter matrix NEWSL
    type OP ! stored as a sparse matrix, there will be 4 of these
        integer(4)  num_pairs, Ntau !num_pairs in this case will equal Lx * Ly / 2 ( we split up odd, even, horiz, vertical)
        integer(4), allocatable :: ind(:,:) !numpairs x2: stores indices of sites connected by each bond
		$op_type$, allocatable :: eta(:,:) ! numpairs x nt: this actually stores the spin values
		integer(4), allocatable :: Q(:,:,:) !numpairs x 4 x 2: stores nearest neighbors Q(j,n,1) = matrix number (1-4) of nth n.n.
					    !(measured clockwise from top right)
						!Q(j,n,2) = link number of nth n.n.
		integer(4), allocatable :: J(:,:) !Lx x Ly; J(ix,iy) is the link number of bond whose first site is (ix,iy)
	end type

    !!!!!! SVD storage
    type SVD
        $G_type$ U, Vd  !!!! Intermediate SVDs
        real(8), allocatable :: S(:)
        integer(4) is_alloc, n1, n2
    end type

    type McParam
        real(8) thop, mu, alpha, V, h, hz,Vt,ht,dht,dVt !Vt and ht are fictional parameters which we may alter to control convergence of global updates
        real(8):: u1, u2 !local interaction parameters: u1 is the nearest-bond coupling, u2 is the on-bond coupling.
        real(8) target_rho
    end type

    type McOpt
        integer(4) num_trial, num_accept, num_trial_global, num_accept_global, num_check_global, size_global_tot
        integer(4) size_global(4), min_cluster_size, num_check, num_adjust_mu, num_add_rho
        real(8) rho
    end type

    type McProf
        real(8) t(10)
    end type

    type McMeasureFlags
		integer(4) :: is_measure_G, is_measure_Gk, is_measure_ETCF, is_measure_ETPC
		integer(4) :: is_measure_tdgf, is_measure_tdgf_K, is_measure_tdcf
		integer(4) :: is_measure_cccf, is_measure_tdpc, is_measure_A, is_measure_E
    	integer(4) :: is_measure_etgf, is_measure_fermions, is_measure_eta_corr
        integer(4) :: is_measure_chi, is_measure_vertex
    	integer(4) :: is_overwrite, n_write_meas
     end type

    !!!! MC configuration
    type McConfig
        real(8) :: dtau            ! MC time step
        integer(4) :: Ntau         ! number of MC time steps
        integer(4) :: Ntau_mul     ! number of MC steps that can be multiplied safely
        integer(4) :: lat_dim      ! Lattice matrix dimension (number of sites, inc. spin and orbital)
        integer(4) :: Lx, Ly       ! Lattice dimensions
        integer(4) :: max_mc_it
        integer(4), allocatable :: ix(:), iy(:), i(:,:) ! conversion tables from spin/orbital/x/y to lattice index and back
        integer(4) :: num_kin_mat
        integer(4) ::n_eq
        type(KinMat), allocatable :: kin_mat(:)  ! kinetic energy matrices
		integer(4) :: Nsvd
		integer(4) :: n, dir, ns
        type(SVD), allocatable :: USV(:)  !!!! Intermediate SVDs
        $G_type$ G            ! Green's function
        real(8) det, det_change, S_eta !NEWSL S_eta
        real(8), allocatable :: S(:)  ! Singular values of 1/G
        type(McParam) par
        integer(4) seed
        integer(4) is_recalc_G, n_write_config, n_write_out, num_global_update
        type(McOpt) opt
        character(100) run_name, purpose
        real(8) :: Fx, Fy		!Flux for x,y directions. I've defined this in a funny way so that Fx is the flux obtained by circling along x (i.e it's actually the magnetic field in the y direction), and similarly for y.
        integer(4) Fz			!Flux along z direction.
        real(8) gauge_factor ! 1>=gauge_factor>=0 determines the gauge for the flux (from Landau along x to symmetric to Landau along y).
        type(McProf) prof, prof_s !prof_s measures real (system) time rather than CPU time
        type(McMeasureFlags) measFlags
        type(SVD) USV_G !G in USV form
        integer(4) is_verbose
        real(8), allocatable::  S0(:,:,:,:)  !the bosonic bare action matrix
        real(8) delta !initial step size for updating a continous field

    	!! Ising parameters:
    	real(8) :: K_cl,K_clt,ps,pt				!classical ising coupling strength in the 't' direction
        real(8) :: weights(5,3,2)

        type(OP), allocatable :: o_p(:) ! ising nematic order parameter


		!! For testing purposes:
		integer(4) is_run_fermions
    end type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Measurements:
    type meas_tdgf
		$G_type$, allocatable:: Gt0(:),G0t(:)  !Gt0=G(tau,0),G0t=G(0,tau) for beta>tau>=0
    end type

    type MCMeasure
		type(mcMeasureFlags) flags
		$G_type$  meanG
		$thop_type$,allocatable :: ETCF(:,:,:,:,:,:,:,:), ETPC(:,:,:,:,:,:,:,:), TDPC(:,:,:,:,:,:,:,:), TDCF(:,:,:,:,:,:,:,:)
		real(8), allocatable :: CCCF(:,:,:), CCCF2(:,:,:,:)
		real(8) E(7)
		type(meas_tdgf) tdgf, mean_tdgf
		$G_type$, allocatable :: G(:), Gh(:)
		integer(4) n_meas !number of measurements
		integer(4) n_write_meas
		integer(4) d(5,2)
		integer(4) minusd(5)
		real(8) Ebonds(2)
        real(8), allocatable :: chi(:,:,:,:) !fermionic nematic correlator and cross term.
        type(cmatrix), allocatable :: Gk1k2t0(:), Gk1k20t(:)
        integer(4) cccf2_N_q ! number of momenta to store in cccf2

        real(8), allocatable :: P_vertex(:,:,:,:)
	end type



    Interface myconj
        module procedure myconj_real, myconj_complex
    end Interface

    contains
        elemental real(8) function myconj_real(x)
            real(8), intent(in) :: x
            myconj_real = x
        end function

        elemental complex(8) function myconj_complex(x)
            complex(8), intent(in) :: x
            myconj_complex = dconjg(x)
        end function
end module global
