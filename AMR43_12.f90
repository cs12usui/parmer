!***********************************************************************
!                               Time-stamp: <2009-02-19 20:54:58 nunami>
!
!       PPPPP        A       RRRRR       MM  MM       EEEEEE     RRRRR
!      PP   PP     AAA      RR   RR     MMM MMM      EE         RR   RR
!     PPPPPP     AA AA     RRRRRR      MM MM MM     EEEEEE     RRRRRR
!    PP        AAAAAAA    RR  RR      MM  M  MM    EE         RR  RR
!   PP        AA    AA   RR    RR    MM      MM   EEEEEEE    RR    RR
!
!
!                          P A R M E R - D D D
!
!     --- full PARticle simulation code 
!                           with adaptic MEsh Refinement technique 
!                        & Dynamic Domain Decomposition capability  ---
!
!   'parmer.F90'   
!
!       by Nunami, Masanori
!          Laboratory of Computer Simulation for Humanospheric Science,
!          Reserch Institute for Sustainable Humanosphere,
!          Kyoto University
!          Gokasho, Uji, Kyoto 611-0011, Japan
!                                  E-mail: nunami@rish.kyoto-u.ac.jp
!                                  PHONE : +81-774-38-4610
!                                  FAX   : +81-774-31-8463
!    'ParmerDDD'
!
!       by Tatsuki Matsui
!          Graduate School of System Informatics
!          Kobe University
!          1-1 Rokkodai, Nada, Kobe, Japan
!                                  E-mail: matsui.tatsuki@gmail.com
!
!  +-----------------------------------------------------------------+
!  |                                                                 |
!  |  Dimensionless units:                                           |
!  |    [coordinate]       x -> x/ \lambda  [coordinate]             |
!  |    [time]             t -> t/T                                  |
!  |    [electric fields]  E -> e E/(m_e \omega c)                   |
!  |    [magnetic fields]  B -> e B/(m_e \omega c)                   |
!  |    [current density]  J -> 8\pi^2 e/(m_e c \omega^2) J          |
!  |    [energy]           expressed in terms of (m_e c^2).          |
!  |                                                                 |
!  |  Units:                                                         |
!  |    T                  period     of incident radiation          |
!  |    \omega             frequency  of incident radiation          |
!  |    \lambda            wavelength of incident radiation          |
!  |    e                  elementary electric charge                |
!  |    m_e                mass of electron                          |
!  |                                                                 |
!  +-----------------------------------------------------------------+
!
! Version History
! AMR4.f90
! Particle movement works for LvMax=1 (but w/o refine_oct).
! AMR5.f90
! Particle movement works with a fixed refined region using refine_oct
! AMR8.f90
! LvMax=1 + no DDD
! AMR9.f90
! LvMax=0 + DDD or LvMax=1 + no DDD
!***********************************************************************
      module param
! +-------------------------------------------------------------------+
! | Module for several parameters                                     |
! +-------------------------------------------------------------------+
!***********************************************************************
        character(LEN=5),parameter :: version = 'AMR31'

        integer(kind=4),parameter :: Nstep0= 0            ! of calculation step (start)
        integer(kind=4),parameter :: Nstep = 3000        ! of calculation step (end)
        integer(kind=4)           :: istep0
        integer(kind=4)           :: istep                ! variable of current step
        integer(kind=4),parameter :: DataIntvl = 10        ! Interval for data output
        integer(kind=4),parameter :: DataIntvl2= 500      ! Interval for data output2
        integer(kind=4),parameter :: PtclIntvl = 20       ! Interval for debug
                                                          ! particle data interval

        real(kind=8),parameter    :: dt    =  0.0025d0     ! time/oct
        real(kind=8),parameter    :: dx(3) =  0.005d0      ! distance/oct
                                                          ! dx(1):dx, dx(2):dy, dx(3):dz
        real(kind=8),parameter    :: pbuff =  0.015d0*7.d0
        integer(kind=4),parameter :: iomp0 = 16           ! thread parallel number
        real(kind=8)              :: dtdx(3), dt2dx(3)    ! dt/dx, 0.5*dt/dx
        real(kind=8)              :: cdx(3), odx(3)       ! variable dx, 1/dx
        real(kind=8)              :: center(3)            ! block center
        real(kind=8)              :: R_lim(0:1,1:3)       ! simulation range
        real(kind=8)              :: R_lim_local(0:1,1:3) ! variable to memory process range

        real(kind=8)              :: charge               ! charge/particle
        real(kind=8),parameter    :: omega   = 1.0d0      ! omgeace/omegape
        real(kind=8),parameter    :: omegace = 0.0d0
        real(kind=8),parameter    :: QM      = 1.d0       ! q/m      (qs/m)/(e/m_e)??? 
        integer(kind=4),parameter :: LvMax   = 2          ! Max refinement level (no change?)
        integer(kind=4),parameter :: LvMax2  = 2**(LvMax+1)  ! LvMax2 = 2**(LvMax+1). This is Nint for Lv=-1
        integer(kind=4),parameter :: intLv   = 2
        integer(kind=4),parameter :: initAMR = 1          ! initial create AMR flag
        integer(kind=4)           :: MaxID(4,-1:LvMax)    ! Max. of cell index in each level. MaxID(2,*) are for the new portion added from the end of Max(1,*) array  
        integer(kind=4)           :: MinID(4,-1:LvMax)    ! Min. of cell index in each level
                                                          ! Min,MaxID(3,iLv):GMesh MinMax
                                                          ! Min,MaxID(4,iLv):new GMesh MinMax

        integer(kind=4) :: MaxIP(-1:LvMax)                ! Max. of particles index in each level
        integer(kind=4) :: irec(0:Lvmax)                  ! flag for copy particle
        integer(kind=4),parameter :: maxN = 1000000       ! Maximal number of particles
        integer(kind=4),parameter :: IonSorts = 2         ! number of sorts of ions (1: ele)
        integer(kind=4) :: PICnumber(IonSorts)            ! number of super particle per cell
        real(kind=8)    :: Rmass(0:IonSorts)              ! ion mass/electron mass
        real(kind=8)    :: Rcharge(0:IonSorts)            ! ion charge/electron charge
        
        real(kind=8)    :: crtr0(0:LvMax), crtr1(0:LvMax) ! criterion
                                                          ! crtr0 is for refine_octB which is required in make_GMesh 
                                                          ! crtr1 is for refine_oct
        real(kind=8),parameter :: TVELE = 0.00225d0
        real(kind=8):: TVION
!        real(kind=8),parameter :: RTEMP=TVELP/TVELE

!Note from Tatsuki: NXB, NYB, and NZB have to be 2^n in order to run properly (at least for now)
        integer(kind=4), parameter:: NXR = 4, NYR = 4, NZR = 4 ! NXR*NYR*NZR Must be equal to the total process number
        integer(kind=4), parameter:: NX = 32, NY=32, NZ=32        ! NX, NY, and NZ have to be equal to or larger than 8.
        integer(kind=4), parameter:: npart_per_cell=32         ! number of particle per cell

        real(kind=8)  :: BD=25.d0
        real(kind=8)  :: fourE(nstep)   ! fourier change
        integer(kind=4) :: MnDisp       ! index = Mn - MnDisp
!yagi added
        !set output type
        integer(kind=4), parameter::AnimeOutput   = 0
        integer(kind=4), parameter::CompOutput    = 0
        character(LEN=7),parameter::comp_version  = 'AMR20'

        !set output timing
        integer(kind=4), parameter::InitialOutput = 1
        integer(kind=4), parameter::LastOutput    = 0

        !choose parameter
        integer(kind=4), parameter::GoctOut       = 0 !Gocts are included into output files.
        integer(kind=4), parameter::Eoutput       = 1 !Efield
        integer(kind=4), parameter::Boutput       = 1 !Bfield
        integer(kind=4), parameter::Joutput       = 0 !charge
        integer(kind=4), parameter::Zoutput       = 1 !density
        integer(kind=4), parameter::Poutput       = 1 !particle
        integer(kind=4), parameter::Routput       = 0 !output file for tecplot
        integer(kind=4), parameter::Doutput       = 0 !domain
        integer(kind=4), parameter::BackBoutput   = 0 !dipole
        integer(kind=4), parameter::eleJoutput    = 0 !only electron
        integer(kind=4), parameter::ionJoutput    = 0 !only ion
        integer(kind=4), parameter::Voutput       = 1 !particle velocity
        integer(kind=4), parameter::ExBoutput     = 0 !ExB velocity
        integer(kind=4), parameter::gradBoutput   = 0 !gradB velocity
        integer(kind=4), parameter::Gridoutput    = 0 !hierarchy grid status !!Need Poutput
        integer(kind=4), parameter::GMeshoutput   = 0 !ghost oct

        !choose trace type parameter
        integer(kind=4)::Ltrace_num=0,Ntrace_num=0
        integer(kind=4), parameter::PTraceOutput  = 0!one particle traced data
        integer(kind=4), parameter::NtraceOutput  = 1!particle number 
        integer(kind=4), parameter::LtraceOutput  = 0!particle loops 
        integer(kind=4), parameter::XtraceOutput  = 1!energy
        integer(kind=4), parameter::RtraceOutput  = 0!output hierarchy cell number
        integer(kind=4), parameter::EtraceOutput  = 0!time evolution of Efield(very heavy)
        integer(kind=4), parameter::TtraceOutput  = 0

        !set histogram parameter
        real(kind=8)   , parameter::HistMax       = 1.0d0 
        real(kind=8)   , parameter::HistMin       = -1.0d0
        integer(kind=4), parameter::HistDatanum   = 100
        integer(kind=4), parameter::VHistOutput   = 0!output particle velocity histogram(heavy)
        integer(kind=4), parameter::PHistOutput   = 0!output particle position histogram(heavy)

        !set model type
        integer(kind=4), parameter::Model         = 5 !watch out PoutIntvl
                                                      !Model=0 electromagnetic wave
                                                      !Model=1 Uniform particles
                                                      !Model=2 moving particle cloud (imbalance model)
                                                      !Model=3 Weibel instability
                                                      !Model=4 Dipole
                                                      !Model=5 inject particle
                                                      !Model=6 beam
                                                      !Model=7 Lunar Gamma with grand
        integer(kind=4), parameter::AMRcondition  = 1 !AMRcondition=1 cut by plasma density
                                                      !AMRcondition=2 cut by magnetic field
                                                      !AMRcondition=3 cut by position
        integer(kind=4), parameter::PoutIntvl     = 1 !interval of outputted particles
        integer(kind=4), parameter::OoctIntvl     = 0   !if this=1 outputted oct is thined down
                                                        !if this=2 outputted oct is 1D
        
        integer(kind=4), parameter::outputOption=0 !change DrawingType for gnuplot
                                                   !0 = symbol
                                                   !1 = lines
        integer(kind=4), parameter::outputQuality=1!change quality of output
                                                   !0 = maximum
                                                   !1 = minimum

        integer(kind=4),parameter::DDD=0
        integer(kind=4),parameter::initDDD=0
        !DDD=1:DDD is using DDD=0:DDD is not using
        !DDD=-1:DDD is debuging. DDD is executed every step.
        integer(kind=4)          ::DDDflag=0
        integer(kind=4)          ::DDDcount=0
        integer(kind=4),parameter::DDDMaxSize=3 !for Lv-1

        integer(kind=4)::debugMode=0 !choose debug Level
        !debugMode <=0 : no debug output
        !debugMode ==1 : watching oct and ptcl connection + routine starting-ending message
        !debugMode ==2 : 1 + detail output
        !debugMode ==3 : 2 + more detail output
        integer(kind=4)::debugChange=-1 !130!59 !debugMode is set to 3 in (debugChange) step
        integer(kind=4)::detailCheck=0
        integer(kind=4)::iFLGDetail=0
        integer(kind=4)::ptclDetail=0
        integer(kind=4)::IDDetail=0
        integer(kind=4)::detailChange=0

        !const for neibouring iPOS
        integer(kind=4),dimension(8,3,0:LvMax)::intvLv
        !const for converting Hierarchy iPOS and Normal iPOS
        integer(kind=4),dimension(-1:LvMax)::sConst
        integer(kind=4),dimension(-1:LvMax)::intNxt

        integer(kind=4),parameter::Nintv=2**LvMax !For Lv0
        integer(kind=4)::MrtNisUsable

        integer(kind=4), parameter:: mempat=1
!add end
        integer(kind=4),parameter:: dipoleFlag = 0
        real(kind=8) :: Dpos(1:3)
        real(kind=8) ::mx,my,mz
        real(kind=8)::flow(3),FIMF(1:6)
        integer(kind=4),parameter::boundaryFlag(1:3)      =(/1,1,1/)
        integer(kind=4),parameter::ParticleDeleteFlag(1:3)=(/1,1,1/)
        integer(kind=4),parameter::ParticleInjectFlag(1:3)=(/1,0,0/)
        integer(kind=4),parameter::wall=9  !<-boundary region oct
        integer(kind=4),parameter::BackgroundParticle=0
        integer(kind=4),parameter::BeamParticle=1
        integer(kind=4),parameter::satellite=0
        integer(kind=4),dimension(1:NY*NZ)::InjectionTable
        integer(kind=4)::tableSize=0,arraySize=0,satelliteArraySize=0
        integer(kind=4)::CIP = 0
        real(kind=8)   ,parameter::injectionFlow = 0.13d0

      end module param
!
!***********************************************************************
      module particle_set
! +-------------------------------------------------------------------+
! | Module for structure 'prtcl'                                      |
! +-------------------------------------------------------------------+
!***********************************************************************
        type prtcl
! -- variables of the particle --
          real(kind=8)         :: R(9)   ! 6-dim. particle coordinate
                                         !  R(1:3) = (x,y,z) at t = dt * n 
                                         !  R(4:6) = (vx, vy, vz) at t = dt (n+1/2)
                                         !  R(7:9) = (vx, vy, vz) at t = dt (n-1/2)
          integer(kind=4)      :: Isort  ! particle species index
                                         !-1 representative particle
                                         ! 0 not use
                                         ! 1 ele.
                                         ! 2 ion1.
                                         ! 3 ion2.
                                         ! 4 ion3...
          integer(kind=4)      :: ioct   ! oct that the particle belongs to(in Mesh)
          integer(kind=4)      :: inum   ! index of particle (in Pesh)
! -- pointers --
          type(prtcl), pointer :: prtNxt ! next particle
		  integer(kind=4)	   :: Vflg
        end type prtcl
      end module particle_set

!***********************************************************************
      module MOD_GMap
! +-------------------------------------------------------------------+
! | Module for structure 'prtcl'                                      |
! +-------------------------------------------------------------------+
!***********************************************************************
        type GMap_type
! -- Auxiliary structure to compute the size Ghost cell
          integer(kind=4)      :: iPOSN(3) !Position index with Normal cordinates
          integer(kind=4)      :: MrtN
          integer(kind=4)      :: proc
          integer(kind=4)      :: octN
       end type GMap_type
      end module MOD_GMap
!
!***********************************************************************
      module oct_set
! +-------------------------------------------------------------------+
! | Module for structure 'oct'                                        |
! +-------------------------------------------------------------------+
!***********************************************************************
        use param
        use particle_set
        type oct
! -- variables of the self oct --
          integer(kind=4)  :: octN        ! oct number in mesh
          integer(kind=4)  :: octLv       ! oct level
          integer(kind=4)  :: octP        ! particle number in the oct
          integer(kind=4)  :: Csort       ! child number for parent oct
          integer(kind=4)  :: iFLG(3)     ! flag for refinement 

!--------------------------------------------------------------------------------------------------
! NOTE
! This flag is the most important information in AMR. In this description, Lower means "Finer", 
! while "Upper" means Coaser resolution layer.
!
! The "actual" physical computation is carried out B/W iFLG(1) = -2 ~ 5
! However, the "REAL" cells are iFLG(1) = 0~3 so the leftover margine iFLG(1) =-2~-1 & 4~5 are discarded (see reconnect_particleP).
!
! iFLG(1) >= 4 : Exising a "REAL" (i.e., Not Overlap) finer level BELOW.
!                      So the field data are pulled out of the finer Level by averaging.
! 
! iFLG(1) = 1~3 : Upper Overlap Region (These are "REAL" cells, and there is No "REAL" finer cells BELOW, but Lower "Overlap" Region exists)
!           iFLG(1) = 1 : The data are handed over to Lower Overlap Region of iFLG(1) = -3 by interpolation 
!                         This is located at the edge of Upper Overlap Region 
!           iFLG(1) = 2 : Data is brought up to Lower Level of iFLG(1) = 4 by averaging
!           iFLG(1) = 3 : Data is brought up to Lower Level of iFLG(1) = 4 by averaging
!
! iFLG(1) = -4~-1 : Lower Overlap Region (So there is no finer cells BELOW.)
!           iFLG(1) = -4 : Doesn't exist so it is to be deleted in sort_oct
!           iFLG(1) = -3 : Interpolated from Upper Overlap Region of iFLG(1) = 1
!           iFLG(1) = -2, -1 : Corresponding to iFLG(1) = 2 & 3 in Upper Overlap Region
!
! iFLG(1) = 0 : Lowest Possible Overlap Region
!                      So there is no finer Level BELOW. The data are brought up to an Upper Level of iFLG(1) = 4
! 
! iFLG(2) = 1:  Indicating existance of children cell 
!               Creation of Children and their pointer connection in add_oct is completed in add_oct
! iFLG(2) = 2 : Call for creation of 8 Children and their pointer connection in add_oct (Correct me if I am wrong!)
! iFLG(3) = 4 : Call for Refinement at set_flag_block
!---------------------------------------------------------------------------------------

          integer(kind=4)  :: iPOS(3)     ! index for grids (This is Hieralcheal Ordering)
          real(kind=8)     :: rPOS(3)     ! position for grids
          integer(kind=4)  :: MrtN        ! Morton number
          integer(kind=4)  :: iC(iomp0)   ! work region for rep-particle
          integer(kind=4) :: prc_bndry ! indicate process boundary buffer (See connect_oct_for_make_GMesh for detail)
!         Process Boundary Buffer is the region of Mesh that has corresponding GMesh cells in other processes.
!
!                     0: This cell is not in the process boundary buffer
!                     1: This cell is in the process boundary buffer
!                     2: This cell is at the edge of process boundary buffer (This flag is not used in this version)
          integer(kind=4) :: octType !indicate the type of oct (-1=BMesh,0=Mesh,1=GMesh)
! -- physical variables --
          real(kind=8)     :: F(18)       ! electromagneitc + current field
                                          ! F(1:3)=E, F(4:6)=B, F(7:9)= J(t+dt/2), F(10:12) = J(t-dt/2)
          real(kind=8)     :: C(6*iomp0)       ! current density (working array) 3*16(num of threads)
          real(kind=8)     :: Z(IonSorts) ! number density, species(Ionsorts)
          real(kind=8)     :: G(1:18)
          real(kind=8)     :: D(1:3)
          real(kind=8)     :: O(1:16)
          integer(kind=4) :: ptcl_loops   ! number of particle_loops
! -- pointers --
          type(oct), pointer :: octPrt    ! parent oct
          type(oct), pointer :: octNb1    ! neighbor oct for -x 
          type(oct), pointer :: octNb2    ! neighbor oct for +x 
          type(oct), pointer :: octNb3    ! neighbor oct for -y 
          type(oct), pointer :: octNb4    ! neighbor oct for +y 
          type(oct), pointer :: octNb5    ! neighbor oct for -z 
          type(oct), pointer :: octNb6    ! neighbor oct for +z
          type(oct), pointer :: octCh1    ! 1st child oct
          type(oct), pointer :: octCh2    ! 2nd child oct
          type(oct), pointer :: octCh3    ! 3rd child oct
          type(oct), pointer :: octCh4    ! 4th child oct
          type(oct), pointer :: octCh5    ! 5th child oct
          type(oct), pointer :: octCh6    ! 6th child oct
          type(oct), pointer :: octCh7    ! 7th child oct
          type(oct), pointer :: octCh8    ! 8th child oct
          type(oct), pointer :: Psort     ! pointer for sorting
! -- pointer for particle --
          type(prtcl), pointer :: ptcl    ! representative particle (0th)
          type(prtcl), pointer :: ptclA   ! representative particle (0th) ! We don't really use this pointer
       end type oct
      end module oct_set

!***********************************************************************
! +-------------------------------------------------------------------+
! | Module for particle injection                                      |
! +-------------------------------------------------------------------+
!***********************************************************************
      module injectionOct
        use param
        type injectionStructure
          integer(kind=4)::octNum, iSort
          real(kind=8)::phi, theta
          real(kind=8)::ratio
        end type injectionStructure
        type(injectionStructure),dimension(1:NY*NZ)::InjectionTableArray
        integer(kind=4),dimension(1:NY*NZ*3)::satelliteTableArray
        
      end module injectionOct

!***********************************************************************
      module time_evolution
! +-------------------------------------------------------------------+
! | Module for initial mesh size                                      |
! +-------------------------------------------------------------------+
!***********************************************************************
        use param
!
        real(kind = 8) :: time
        integer(kind=4),dimension(0:2**(LvMax+2)) :: jstep1,jstep2
        integer(kind=4),dimension(0:LvMax+1)      :: jtemp
!
      end module time_evolution
!
!***********************************************************************
      module init_mesh_size
! +-------------------------------------------------------------------+
! | Module for initial mesh size                                      |
! +-------------------------------------------------------------------+
!***********************************************************************
        use param
        use oct_set
        use particle_set
!===== Parameters for Region Size =====
        integer(kind=4), parameter:: NXYR=NXR*NYR
        integer(kind=4), parameter:: NXB=NX, NYB=NY, NZB=NZ
        integer(kind=4), parameter:: NXB_2=NXB/2, NYB_2=NYB/2, NZB_2=NZB/2! NXB_2=NXB/2 for Lv -1 operation
        integer(kind=4), parameter:: NXYB_2=NXB_2*NYB_2, NYZB_2=NYB_2*NZB_2, NZXB_2=NZB_2*NXB_2
        integer(kind=4), parameter:: NXYZB_2=NXB_2*NYB_2*NZB_2! NXYB,NXYZB for Lv -1 operation
        integer(kind=4), parameter:: NX_2t = NXR*NXB_2, NY_2t=NYR*NYB_2, NZ_2t=NZR*NZB_2! NXB_2=NXB/2 for Lv -1 operation
        integer(kind=4), parameter:: Nall_ini=NXB*NYB*NZB !Nall_ini is the initial total number of cells in a simulation space 
        integer(kind=4)           :: Nall ! Remember, Nall is no longer a constant in DDD version 
        integer(kind=4), parameter :: Niall = Nall_ini*npart_per_cell !This is the number of ions / process
                                                                      !So there are 2*Niall actual particles / process in the simulation
!===== Parameters for Mesh =====
        !for AMR parameters
        real(kind=8),parameter:: refineRatio(-1:10) = (/ 1.d0, 0.5d0, 0.5d0, 0.0d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
        !                                                Lv-1,   Lv0,   Lv1,  Lv2,  Lv3,  Lv4,  Lv5,  Lv6,  Lv7,  Lv8,  Lv9, Lv10
        integer(kind=4)::LevelSize(-1:10),GLevelSize(-1:10)
        integer(kind=4),parameter::BMeshSize=Nall_ini/8*20
        integer(kind=4)::GMeshSizeIni
        integer(kind=4)::MeshSizeIni      
        integer(kind=4)::MeshSize        !MeshSize=MeshSizeIni+BMeshSize+GMeshSizeIni
        type(oct),dimension(:),allocatable,target:: Mesh
        type(oct),dimension(:),allocatable,target:: Mesh2
        integer(kind=4)::BMeshBound
        integer(kind=4)::MeshBound !limitation of Mesh array
                                   !Mesh [1~BMeshBound=BMesh] [BMeshBound+1~MeshBound=Mesh] [MeshBound+1~ =GMesh]     

        integer(kind=4):: GMeshSize_whole !The length of the whole GMesh from Lv=-1 to Lv=Lvmax
        integer(kind=4):: GMeshSize       !The size of GMesh of Lv=-1

!===== Parameters for Pesh =====       
        integer(kind=4)::PeshSize
        type(prtcl),dimension(:,:),target,allocatable:: Pesh
        type(prtcl),dimension(:),target,allocatable:: Pesh2 
        integer(kind=4),dimension(:,:),target,allocatable:: iCp!,iCp2

!===== Parameter for Index ===== 
        integer(kind=4), allocatable:: Mn2CPU(:),Mn2CPU_old(:) ! Records the ending of the Mtn number for each CPU
        integer(kind=4) ::MaxMn, MinMn

      end module init_mesh_size
!
!***********************************************************************
      module const
! +-------------------------------------------------------------------+
! | Module for numerical & physical constants                         |
! +-------------------------------------------------------------------+
!***********************************************************************
        real(kind=8), parameter :: ZERO   = 0.0d0
        real(kind=8), parameter :: EIGHTH = 0.125d0
        real(kind=8), parameter :: QUARTER= 0.25d0
        real(kind=8), parameter :: THIRD  = 0.333333333333333333333d0
        real(kind=8), parameter :: HALF   = 0.5d0
        real(kind=8), parameter :: THR_FOURTH = 0.75d0
        real(kind=8), parameter :: ONE    = 1.0d0
        real(kind=8), parameter :: TWO    = 2.0d0
        real(kind=8), parameter :: PI     = 3.141592653589793238462d0
        real(kind=8), parameter :: PI2    = PI + PI

        real(kind=8), parameter :: emass = 0.51099892D0 ! electron mass in [MeV]

      end module const
!***********************************************************************
      module init_condition
! +-------------------------------------------------------------------+
! | Module for initial condition                                      |
! +-------------------------------------------------------------------+
!***********************************************************************
        real(kind=8)      :: BiniN(3),EiniN(3) ! amplitude
        real(kind=8)      :: kiniN(3)          ! wave number
        real(kind=8)      :: phi               ! inital phase
      end module init_condition

!******************ADDED by TATSUKI******************************
      module message_passing_interface
        !use mpi
        use oct_set
        use init_mesh_size
        implicit none
        include "mpif.h"
        real(kind=8)::wtime !elapse time

        integer(kind=4) ierr, nprocs, rank
        integer(kind=4), parameter:: buf_width=2
        integer(kind=4), parameter:: buf_width_J=0, buf_width_P=0! width of ghost cells. Should be multiples of 2 for Lv -1 operation. 

        !GMesh related variables
        integer(kind=4), dimension(:), allocatable:: n_rcells_proc, n_gcells_proc
        integer(kind=4), dimension(:), allocatable:: n_rcells_sum, n_gcells_sum
        integer(kind=4), dimension(:), allocatable:: iBMesh_arr, iGMesh_arr
        integer(kind=4), dimension(:), allocatable:: GMn2octN

        !variables for refresh routines
        integer(kind=4)::MPIbufsize
        integer(kind=4)::send_onum(2),recv_onum(2),send_pnum(2),recv_pnum(2)
        integer(kind=4),dimension(NXR*NYR*NZR)::send_glist,recv_glist,send_rlist,recv_rlist
        integer(kind=4)::sglistnum,rglistnum,srlistnum,rrlistnum
        real(kind=8),dimension(:,:),allocatable::sbuf_BFR, rbuf_GFR    !F(1:12)=F(1:12)
        real(kind=8),dimension(:,:),allocatable::sbuf_G, rbuf_G    !G(1:18)=G(1:18)
        integer(kind=4),dimension(:,:),allocatable::sbuf_BFI, rbuf_GFI !FI(1)=iFLG FI(2)=OctP
        real(kind=8),dimension(:,:),allocatable::sbuf_GPR, rbuf_BPR    !PR(1:9)=R(1:9)
        integer(kind=4),dimension(:,:),allocatable::sBuf_GPI, rbuf_BPI !PI(1)=isort
        integer(kind=4),dimension(:,:),allocatable:: rbuf_BPO, sbuf_GPO

        !variables for exchange_moved_oct
        integer(kind=4)::DDDbufsize
        real(kind=8),dimension(:),allocatable::sbufL_FR,sbufR_FR,rbufL_FR,rbufR_FR    !FR(1:12)=F(1:12)
        integer(kind=4),dimension(:),allocatable::sbufL_FI,sbufR_FI,rbufL_FI,rbufR_FI !FI(1)=iFLG FI(2)=OctP
        real(kind=8),dimension(:),allocatable::sbufL_PR,sbufR_PR,rbufL_PR,rbufR_PR    !PR(1:9)=R(1:9)
        integer(kind=4),dimension(:),allocatable::sBufL_PI,sbufR_PI,rbufL_PI,rbufR_PI !PI(1)=isort

        integer(kind=4),dimension(4)::sendSize,recvSize !sendSize(1) is the # of CELLs to be sent to the LEFT
                                                        !sendSize(2) is the # of PARTICLESs to be sent to the LEFT
                                                        !sendSize(3) is the # of CELLs to be sent to the RIGHT
                                                        !sendSize(4) is the # of PARTICLESs to be sent to the RIGHT
                                                        !1=Lsize_F,2=Lnum_P,3=Rsize_F,4=Rnum_P
        integer(kind=4)::LeftSendMin,LeftSendMax,RightSendMin,RightSendMax,SavedMin,SavedMax
        integer(kind=4)::LeftRecvMin,LeftRecvMax,RightRecvMin,RightRecvMax,ExistMin,ExistMax
        integer(kind=4)::next,prev,leftDisp,rightDisp
        integer(kind=4),dimension(-1:LvMax)::LeftNum,RightNum

        !for DDD
        integer(kind=4)::nloop

        integer(kind=4)::temptime,steptime
        
        integer(kind=4)::px,py,pz

      end module message_passing_interface

! ---------------- Module for Mersenne twister: A random number generator ------------------
      module mtmod
        implicit none
        ! Default seed
        integer,parameter::defaultsd=4357
        ! Period parameters
        integer,parameter::N=624, N1=N+1  
        ! the array for the state vector
        integer,save,dimension(0:N-1)::mt
        integer,save::mti=N1
      contains
        !Initialization subroutine
        subroutine sgrnd(seed)
          implicit none
          ! setting initial seeds to mt[N] using
          ! the generator Line 25 of Table 1 in
          ! [KNUTH 1981, The Art of Computer Programming
          ! Vol. 2 (2nd Ed.), pp102]
          integer,intent(in)::seed
          mt(0)=iand(seed,-1)
          do mti=1,N-1
             mt(mti)=iand(69069*mt(mti-1),-1)
          enddo
          return
        end subroutine sgrnd
        !Random number generator
        real(8) function grnd()
          implicit integer(a-z)
          ! Period parameters
          integer,parameter::M=397, MATA=-1727483681
          ! constant vector a
          integer,parameter::LMASK=2147483647
          ! least significant r bits
          integer,parameter::UMASK=-LMASK-1
          ! most significant w-r bits
          ! Tempering parameters
          integer,parameter::TMASKB=-1658038656, TMASKC=-272236544
          dimension mag01(0:1)
          data mag01/0, MATA/
          save mag01
          ! mag01(x) = x * MATA for x=0,1
          TSHFTU(y)=ishft(y,-11)
          TSHFTS(y)=ishft(y,7)
          TSHFTT(y)=ishft(y,15)
          TSHFTL(y)=ishft(y,-18)
          if(mti.ge.N) then
             ! generate N words at one time
             if(mti.eq.N+1) then
                ! if sgrnd() has not been called,
                call sgrnd( defaultsd )
                ! a default initial seed is used
             endif
             do kk=0,N-M-1
                y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
                mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
             enddo
             do kk=N-M,N-2
                y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
                mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
             enddo
             y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
             mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
             mti = 0
          endif
          y=mt(mti)
          mti=mti+1 
          y=ieor(y,TSHFTU(y))
          y=ieor(y,iand(TSHFTS(y),TMASKB))
          y=ieor(y,iand(TSHFTT(y),TMASKC))
          y=ieor(y,TSHFTL(y))
          if(y.lt.0) then
             grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
          else
             grnd=dble(y)/(2.0d0**32-1.0d0)
          endif
          return
        end function grnd
      end module mtmod

module Vparam
	use param
	use oct_set
	implicit none
	integer(kind=4),parameter::VNX=1,VNY=NY,VNZ=NZ
	integer(kind=4)::VMeshSize=VNX*VNY*VNZ
	type(oct),allocatable,target,dimension(:) :: VMesh
end module Vparam
		
	  
	  
!***********************************************************************
!******************//                               //******************
!*****************//     START of main program     //*******************
!****************//                               //********************
!
program PARMER
  use oct_set
  use init_mesh_size
  use const
  use param
  use time_evolution
  use message_passing_interface

  implicit none
  integer(kind=4)    :: i
  
  if(Nstep0+1.gt.Nstep) then 
     print *,"error Nstep0+1.gt.Nstep"
     stop
  end if
  
  time = 0.d0
  irec = 0

!---MPI Initialization--- 
  call MPI_init(ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, rank, ierr)
  call nprocs_check
  
  
  call set_BufferLength
  call allocate_Memory
!initial equation domain setup 
  if(DipoleFlag==1)then
    if(Model==7)then
        Dpos(1) = dx(1) * NXR * NX - dble(wall)*dx(1)
    else
        Dpos(1) = dx(1) * NXR * NX * 0.5d0
    end if
    Dpos(2) = dx(2) * NYR * NY * 0.5d0
    Dpos(3) = dx(3) * NZR * NZ * 0.5d0
  end if
  call initial_setting(0)
  call setup_model   
! Make injection region table
  if(BackgroundParticle==1 .and. px==0)then
     call makeInjectionOctTable(tableSize)
  end if
  if(beamParticle==1)then
     call makeInjectionOctTableArray(arraySize)
  end if
  if(satellite==1)then  
	 call makeSatelliteTableArray(satelliteArraySize)
  end if
  if(initAMR/=0)call initial_AMR
  call output_traces
  if(initDDD/=0)call initial_DDD

!-----------------------------------------------       
  Center(1) = dx(1) * NXR * NX * 0.5d0
  Center(2) = dx(1) * NYR * NY * 0.5d0
  Center(3) = dx(1) * NZR * NZ * 0.5d0
  
  wtime=MPI_Wtime()
  steptime=MPI_Wtime()
  
  if(InitialOutput/=0)call output_params(1)
!  if(InitialOutput/=0) call output_traces
! --- calc ptcl density ---
  do i=Lvmax,0,-1
     call density(i)
  end do
! -------------------------------------
  call mpi_barrier(MPI_COMM_WORLD,ierr)
!--------------------------------------
  do istep = Nstep0+1, Nstep
     time = time + dt
     
     if(rank==0)print *, '******* time =',time,'step=', istep,' *********'
     if(istep==debugChange)debugMode=3
     if(istep==detailChange)detailCheck=1
     
     istep0 = istep   

! -------- AMR routines --------
     if( LvMax>0) then
        call refine_oct 
        call sort_oct   
     end if
! ------------------------------
! -- Time Evolution - 
     irec(:)=0
     irec(0)=99
     do i=1,jtemp(Lvmax+1)
        call advance_fieldT(jstep1(i),jstep2(i)) !iLv is entered in jstep1
     end do
	 
! -- Sorting (particle)
     call sort_particle                   
!--- DDD ---
     if(DDD/=0)then
        !call fipp_start()
        call DDD_check
        if(DDDflag/=0)then
           call reset_Gst!
           call load_balance
           call DDD_sendrecv
           call make_GMesh
           if(LvMax>0)then
              call refine_oct_DDD
              call sort_oct
           end if
        end if
        !call fipp_stop()
     end if
!------ 
! --- calc ptcl density ---
     do i=Lvmax,0,-1
        call density(i)
     enddo
!--------- output routines ------
     if(mod(istep,DataIntvl)==0)then
!        call output_params(PoutIntvl)
        call output_traces
     end if
     if(mod(istep,DataIntvl2)==0 )then
        call output_params(PoutIntvl)
     end if
     if(detailCheck/=0)then
        call output_details
     end if
!---------------------------------
  end do
  
  !call checkVParticle
  
  if(rank==0)  print*,"rank=",rank,"CPUtime=",MPI_Wtime()-wtime,"DDDcount=",DDDcount
  
  if(LastOutput/=0)call output_params(npart_per_cell)     
  call output_traces
  
  call deallocate_Memory
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if(rank==0)print *,"parmer completed"
  call MPI_finalize(ierr)
  
  
!
!******************//                               //******************
!*****************//     END of main program       //*******************
!****************//                               //********************
!***********************************************************************
!       
!  contains  
end program PARMER

!***********************************************************************
      subroutine initial_setting(icon)
! +-------------------------------------------------------------------+
! +    set up initial parameters and initial conditions               +
! +-------------------------------------------------------------------+
!***********************************************************************
        use param
        use init_condition
        use const
        use message_passing_interface
        use oct_set
        implicit none       
        integer(kind=4),intent(in)::icon
        
        !if(debugMode>=1)print *,"***** initial setting icon=",icon,"*****",rank
! -- Preparation --
        if(icon==0)then
           call MeshSize_check       
           call preparation
           call CFL_check
           call set_conversionConst
! -- initial setting for time evolution
           call set_evolution
           call set_param
           call setup_parameter
           call set_AMR_Criterion
! -- Initial setting for octs --
           call set_MnMinMaxID
        else
           call set_MnMinMaxID_DDD
        endif
        call set_BMesh(icon)
        call set_RlimLocal
        call set_Mesh

        call initial_connection(icon)

        call make_GMesh

        call full_oct_check   

!        call periodic_boundary      
!        call define_sim_bounds
        return

!************************************************************************************
!routines related to initial_setting
      contains

!**************************************************************
!=============== Parameter Checking Routines ==================
!**************************************************************

        !yagi 2012/01/10
        subroutine MeshSize_check
          implicit none

          MrtNisUsable=0

          if(NXB==0 .or. NYB==0 .or. NZB==0)then
             print *,"Error : NXB, NYB, NZB are out of bounds",NXB,NYB,NZB
             stop
          endif

          if(NXB==NYB .and. NYB==NZB)then
             !check NXB,NYB,NZB  If Size is exponentiation of 2.
             if(rank==0)then
                print *,""
                print *,"IAND(NXB,NXB-1)==",IAND(NXB,NXB-1)
                print *,"IAND(NYB,NYB-1)==",IAND(NYB,NYB-1)
                print *,"IAND(NZB,NZB-1)==",IAND(NZB,NZB-1)
                print *,""
             endif
             if(IAND(NXB,NXB-1)==0 .and. IAND(NYB,NYB-1)==0 .and. IAND(NZB,NZB-1)==0)then
                MrtNisUsable=1
                if(rank==0)then
                   print *,""
                   print *,"*************************************************"
                   print *," NXB,NYB and NZB is exponentiation of 2"
                   print *," Morton Number is chosen to order Octs"
                   print *,"*************************************************"
                   print *,""
                endif
             endif
          endif
          if(MrtNisUsable==0)then
             if(rank==0)then
                print *,""
                print *,"?????????????????????????????????????????????????"
                print *," NXB,NYB and NZB is NOT exponentiation of 2"
                print *," Raster scan is chosen to order Octs"
                print *,"?????????????????????????????????????????????????"
                print *,""
             endif
!             if(NYR/=1 .or. NZR/=1)then
!                print *,"NYR NZR must be 1 NYR=",NYR,"NZR=",NZR
!                stop
!             endif
          endif
          
        end subroutine MeshSize_check
        !
!***********************************************************************************
        subroutine CFL_check
          ! +-------------------------------------------------------------------+
          ! |                                                                   |
          ! | Check the Courant-Friedrichs-Levy condition                       |
          ! |                                                                   |
          ! +-------------------------------------------------------------------+
          !***********************************************************************
          ! ---------------------------------
          implicit none
          real(kind=8) :: dt_max
          ! ---------------------------------
          !
          dt_max=ONE/SQRT(ONE/dx(1)**2 + ONE/dx(2)**2+ ONE/dx(3)**2)
          !

          if (dt > dt_max) then
             write(*,*) ' => dt is too large, Maxwell solver will be unstable!'
             stop
          endif
          if (dt < HALF*dt_max) then
             write(*,*) ' => dt is too small, phase error at small k will be significant!'
          else
          endif
          !
          !print *, ' '
          !
          return

        end subroutine CFL_check

!-------------------------------------------------------------
!================ Parameter Setting Routines =================
!-------------------------------------------------------------

        subroutine preparation
          ! +-------------------------------------------------------------------+
          ! |                                                                   |
          ! | Preparation routine                                               |
          ! |                                                                   |
          ! +-------------------------------------------------------------------+
          !***********************************************************************
          implicit none
          integer(kind=4)   :: k
          integer(kind=4),dimension(1:3)::rankPOS
          !
          ! == BEGIN =============================================================
          !
          ! -- Starting Massage ---
!!$        print *, '                                                              '
!!$        print *, '=============================================================='
!!$        print *, '       PPPPP        A       MM  MM       EEEEEE     RRRRR     '
!!$        print *, '      PP   PP     AAA      MMM MMM      EE         RR   RR    '
!!$        print *, '     PPPPPP     AA AA     MM MM MM     EEEEEE     RRRRRR      '
!!$        print *, '    PP        AAAAAAA    MM  M  MM    EE         RR  RR       '
!!$        print *, '   PP        AA    AA   MM      MM   EEEEEEE    RR    RR      '
!!$        print *, '                                                              '
!!$        print *, '                                       Masanori Nunami        '
!!$        print *, '                                       RISH, Kyoto University '
!!$        print *, '=============================================================='
!!$        print *, '                                                              '
!!$!
          ! -- Prepare general mesh variables --
          !
          do k=1,3
             dtdx(k)  = 0.d0
             dt2dx(k) = 0.d0
             odx(k)   = ONE/dx(k)
          end do

          ! --- Boundary of calculation box ---

          R_lim(0,1) = ZERO
          R_lim(1,1)= NXB*NXR*dx(1)
          R_lim(0,2) = ZERO
          R_lim(1,2) = NYB*NYR*dx(2)
          R_lim(0,3) = ZERO
          R_lim(1,3) = NZB*NZR*dx(3)

!!$          if(debugMode>2)then
!!$             print *,"R_lim(0,1)=",R_lim(0,1),rank
!!$             print *,"R_lim(1,1)=",R_lim(1,1),rank
!!$             print *,"R_lim(0,2)=",R_lim(0,2),rank
!!$             print *,"R_lim(1,2)=",R_lim(1,2),rank
!!$             print *,"R_lim(0,3)=",R_lim(0,3),rank
!!$             print *,"R_lim(1,3)=",R_lim(1,3),rank
!!$          endif

          ! --- Center of calculation box ---
          center(1) = (R_lim(0,1) + R_lim(1,1)) * HALF
          center(2) = (R_lim(0,2) + R_lim(1,2)) * HALF
          center(3) = (R_lim(0,3) + R_lim(1,3)) * HALF
          
          !oct process position
          if(MrtNisUsable/=0)then
            call Inverse_MortonN(rank,rankPOS)
            px = rankPOS(1)-1
            py = rankPOS(2)-1
            pz = rankPOS(3)-1
          else 
            px=mod(rank,NXR)
            py=mod(rank/NXR,NYR)
            pz=mod(rank/(NXR*NYR),NZR)
          end if
          return

        end subroutine preparation

!****************************************************************************************
        !2011/04/19 -> 2012/02/13
        !this subroutine is called in initial_setting
        subroutine set_conversionConst
          implicit none
          integer(kind=4)::intLv0,iLv

          do iLv=0,LvMax
             intLv0= 2**(LvMax-iLv)

             intvLv(1,1,iLv)=intLv0*(-1)
             intvLv(1,2,iLv)=intLv0*(-1)
             intvLv(1,3,iLv)=intLv0*(-1)
             intvLv(2,1,iLv)=intLv0*( 1)
             intvLv(2,2,iLv)=intLv0*(-1)
             intvLv(2,3,iLv)=intLv0*(-1)
             intvLv(3,1,iLv)=intLv0*(-1)
             intvLv(3,2,iLv)=intLv0*( 1)
             intvLv(3,3,iLv)=intLv0*(-1)
             intvLv(4,1,iLv)=intLv0*( 1)
             intvLv(4,2,iLv)=intLv0*( 1)
             intvLv(4,3,iLv)=intLv0*(-1)
             intvLv(5,1,iLv)=intLv0*(-1)
             intvLv(5,2,iLv)=intLv0*(-1)
             intvLv(5,3,iLv)=intLv0*( 1)
             intvLv(6,1,iLv)=intLv0*( 1)
             intvLv(6,2,iLv)=intLv0*(-1)
             intvLv(6,3,iLv)=intLv0*( 1)
             intvLv(7,1,iLv)=intLv0*(-1)
             intvLv(7,2,iLv)=intLv0*( 1)
             intvLv(7,3,iLv)=intLv0*( 1)
             intvLv(8,1,iLv)=intLv0*( 1)
             intvLv(8,2,iLv)=intLv0*( 1)
             intvLv(8,3,iLv)=intLv0*( 1)

          enddo

          do iLv=-1,LvMax
             sConst(iLv)=2**(LvMax-iLv)
             intNxt(iLv)=2**(LvMax-iLv+1)
          enddo

        end subroutine set_conversionConst

!***************************************************************************************
        !yagi 2011/10/27
        subroutine set_MnMinMaxID
          implicit none
          integer(kind=4)::index,i
          integer(kind=4) :: Nall_8

          Nall_8=NXB_2*NYB_2*NZB_2 

          !set MnDisp
          MnDisp=rank*Nall_8-1

          !set Mn2CPU
          do index=1, nprocs
             Mn2CPU(index)=Nall_8 *index-1
          enddo

          !--- Specify the Min and Max Morton on each process.
          if(rank==0)then
             MinMn=0
          else
             MinMn=Mn2CPU(rank)+1
          endif
          MaxMn=Mn2CPU(rank+1)

          !----------------------------
          ! Level -1
          MinID(1,-1) = 1
          MaxID(1,-1) = Nall_ini/8
          MinID(2,-1) = Nall_ini/8
          MaxID(2,-1) = Nall_ini/8

          !----------------------------
          ! Level 0

          MinID(1,0)  = BMeshBound+1
          MaxID(1,0)  = BMeshBound+Nall_ini
          MinID(2,0)  = BMeshBound+Nall_ini
          MaxID(2,0)  = BMeshBound+Nall_ini

          !----------------------------
          ! Level > 1   We do not create Mesh for Lv > 1 in this routine

          if(LvMax>0) then 
             do i=1,LvMax
                MinID(1,i) = BMeshBound+Nall_ini
                MaxID(1,i) = BMeshBound+Nall_ini
                MinID(2,i) = BMeshBound+Nall_ini
                MaxID(2,i) = BMeshBound+Nall_ini
             enddo
          end if
        end subroutine set_MnMinMaxID

        subroutine set_MnMinMaxID_DDD
          implicit none
          integer(kind=4)::Lv0Size,i

          !new MnMin and MnMax have already set at load_balance
          MnDisp=MinMn-1
          !if(debugMode>=3)print *,"MnDisp is set to ",MnDisp,rank
          !-------------------------
          ! Level -1
          MinID(1,-1)=1
          MaxID(1,-1)=MaxMn-MinMn+1
          MinID(2,-1)=MaxMn-MinMn+1
          MaxID(2,-1)=MaxMn-MinMn+1

          Lv0Size=(MaxID(1,-1)-MinID(1,-1)+1)*8
          !Level 0
          MinID(1,0)=BMeshBound+1
          MaxID(1,0)=BMeshBound+Lv0Size
          MinID(2,0)=BMeshBound+Lv0Size
          MaxID(2,0)=BMeshBound+Lv0Size

          !Level > 1
          if(LvMax>0)then
             do i=1,LvMax
                MinID(1,i) = BMeshBound+Lv0Size
                MaxID(1,i) = BMeshBound+Lv0Size
                MinID(2,i) = BMeshBound+Lv0Size
                MaxID(2,i) = BMeshBound+Lv0Size
             enddo
          endif

        end subroutine set_MnMinMaxID_DDD
!******************************************************************************************
        !yagi 2011/10/27
        subroutine set_RlimLocal
          implicit none
          integer(kind=4)::sID,eID

          sID=MinID(1,-1)
          eID=MaxID(1,-1)

          R_lim_local(0,1) = minval(Mesh(sID:eID) % rPOS(1))
          R_lim_local(1,1) = maxval(Mesh(sID:eID) % rPOS(1))
          R_lim_local(0,2) = minval(Mesh(sID:eID) % rPOS(2))
          R_lim_local(1,2) = maxval(Mesh(sID:eID) % rPOS(2))
          R_lim_local(0,3) = minval(Mesh(sID:eID) % rPOS(3))
          R_lim_local(1,3) = maxval(Mesh(sID:eID) % rPOS(3))

          R_lim_local(0,:) = R_lim_local(0,:) - dx(:)
          R_lim_local(1,:) = R_lim_local(1,:) + dx(:)

!!$          if(debugMode>2)then
!!$             print *,"R_lim_local(0,1)=",R_lim_local(0,1),rank
!!$             print *,"R_lim_local(1,1)=",R_lim_local(1,1),rank
!!$             print *,"R_lim_local(0,2)=",R_lim_local(0,2),rank
!!$             print *,"R_lim_local(1,2)=",R_lim_local(1,2),rank
!!$             print *,"R_lim_local(0,3)=",R_lim_local(0,3),rank
!!$             print *,"R_lim_local(1,3)=",R_lim_local(1,3),rank
!!$          endif

        end subroutine set_RlimLocal

!******************************************************************************************
        !yagi 2012/01/10
        subroutine set_param
          implicit none
          integer(kind=4)::i

          do i=1,3
             dtdx(i)  = dt/dx(i)
             dt2dx(i) = HALF*dt/dx(i)
             odx(i)   = ONE/dx(i)
          end do
          
        end subroutine set_param
        !**********************************************************************
        !yagi 2011/12/26
        subroutine set_AMR_Criterion
          implicit none

          crtr0(:) = 9999.d0
          crtr1(:) = 9999.d0

          ! -- criterion setting
          if(AMRcondition==1)then
            if(LvMax.eq.1) then
                crtr0(0)=9999.d0 ! 0.40d0
                crtr1(0)=dble(npart_per_cell)*0.4d0
                crtr0(1)=9999.d0 ! 0.40d0
                crtr1(1)=dble(npart_per_cell)*0.4d0
            endif
            if(LvMax.eq.2) then
                crtr0(0)=999.d0 ! 0.40d0
                crtr1(0)=dble(npart_per_cell)*0.7
                crtr0(1)=999.d0 ! 0.40d0
                crtr1(1)=dble(npart_per_cell)*2.0
            !    crtr0(2)=999.d0 ! 0.40d0
            !    crtr1(2)=dble(npart_per_cell)*0.0
            endif
            if(LvMax.eq.3) then
                crtr0(0)=999.d0 ! 0.40d0
                crtr1(0)=dble(npart_per_cell)*2.0
                crtr0(1)=999.d0 ! 0.40d0
                crtr1(1)=dble(npart_per_cell)*4.0
!                crtr0(2)=999.d0 ! 0.40d0
!                crtr1(2)=dble(npart_per_cell)*8.0
            endif
          else if(AMRcondition==2)then
            if(LvMax.eq.1) then
                crtr0(0)=9999.d0 ! 0.40d0
                crtr1(0)=2.3d0
                crtr0(1)=9999.d0 ! 0.40d0
                crtr1(1)=9999.d0
            endif
            if(LvMax.eq.2) then
                crtr0(0)=999.d0 ! 0.40d0
                crtr1(0)=dble(npart_per_cell)*0.9
                crtr0(1)=999.d0 ! 0.40d0
                crtr1(1)=dble(npart_per_cell)*1.1
            !    crtr0(2)=999.d0 ! 0.40d0
            !    crtr1(2)=dble(npart_per_cell)*0.0
            endif
          else if(AMRcondition==3)then
            if(LvMax.eq.1) then
                crtr0(0)=dx(1)*NXR*NX*0.5d0-8.d0*dx(1) !Dpos(1)-36.d0*dx(1) ! 0.40d0
                crtr1(0)=dx(1)*NXR*NX*0.5d0+8.d0*dx(1) !Dpos(1)
                crtr0(1)=9999.0d0 ! 0.40d0
                crtr1(1)=9999.0d0
            endif
            if(LvMax.eq.2) then
                crtr0(0)=2.59d0 ! 0.40d0
                crtr1(0)=5.09d0
                crtr0(1)=3.2d0 ! 0.40d0
                crtr1(1)=4.48d0
!                crtr0(2)=3.2d0 ! 0.40d0
!                crtr1(2)=4.48d0
            endif
          end if
        end subroutine set_AMR_Criterion

        !************************************************************************
        subroutine set_evolution
          ! +---------------------------------------------------------------------+
          ! |     define some working array for the computation of time evolution |
          ! +---------------------------------------------------------------------+
          !************************************************************************
          use const
          use param
          use time_evolution
          ! --------------------------------------------------
          implicit none
          integer(kind=4)   :: i,itemp
          !
          jtemp(0)=0
          jtemp(1)=1
          do i=2,LvMax+1
             jtemp(i)=jtemp(i-1)+2**(i-1)
          end do
          do i=1,LvMax+1
             jtemp(i)=jtemp(i)*2
          end do

          itemp=jtemp(LvMax+1)
          call make_step(0,itemp)

          return
        end subroutine set_evolution
        !
        !************************************************************************
        recursive subroutine make_step(iLv,ii)
          ! +---------------------------------------------------------------------+
          ! |     define some working array for the computation of time evolution |
          ! +---------------------------------------------------------------------+
          !************************************************************************
          use const
          use param
          use time_evolution
          ! --------------------------------------------------
          implicit none
          integer(kind=4)  :: iLv,ii
          !
          jstep1(ii)=iLv
          jstep2(ii)=2
          if(iLv/=LvMax) call make_step(iLv+1,ii-1)

          ii=ii-jtemp(LvMax-iLv)-1
          jstep1(ii)=iLv
          jstep2(ii)=1
          if(iLv/=LvMax) call make_step(iLv+1,ii-1)

          return
        end subroutine make_step
        !************************************************************************

!====================================================================
!----------------------- Oct Setting Routines -----------------------
!====================================================================

        subroutine set_BMesh(icon)
          ! +-------------------------------------------------------------------+
          ! |                                                                   |
          ! | Get Morton index from each cell  and 
          ! | create Besh in that order
          ! |                                                                   |
          ! |  Morton ordering:                                                 |
          ! |   From bits of each cell index (i,j,k), ....                      |
          ! |                                                                   |
          ! |   (Example)                                                       |
          ! |     cell index = (3,4) in two-dimensional 2^3 x 2^3 meshes        |
          ! |     (3,4) => (011,100) => L = 100101 = 37.                        |
          ! |                                                                   |
          ! +-------------------------------------------------------------------+
          !***********************************************************************
          implicit none
          integer(kind=4),parameter :: N2Max=10
          !        real(kind=8) :: xx,yy,zz,xx2,yy2,zz2,rr
          integer(kind=4) :: ix,iy,iz,index
          integer(kind=4) :: i
          !        integer(kind=4) :: N2All
          integer(kind=4) :: n2x(1:N2Max),n2y(1:N2Max),n2z(1:N2Max)
          integer(kind=4) :: Nall_8
          integer(kind=4) :: sID,eID
          !---Variables necessary for BMesh
          integer(kind=4)    :: octN,octLv,octP
          integer(kind=4)    :: iFLG(3),Csort
          integer(kind=4)    :: iPOS(3)
          real(kind=8)       :: rPOS(3)
          integer(kind=4)    :: MrtN
          integer(kind=4)  :: prc_bndry, ptcl_loops
          integer(kind=4)    :: AllocateStatus
          integer(kind=4)    :: iC(iomp0)
          real(kind=8)       :: F(18),C(6*iomp0),Z(IonSorts),G(1:18),D(1:3),O(1:16)
          integer(kind=4),intent(in)::icon
          ! -- pointers --
          type(oct), pointer   :: newp
          ! -------------------------------------------
          !        print *, ' in Morton_number ..."'

          n2x=0
          n2y=0
          n2z=0
          Nall_8=NXB_2*NYB_2*NZB_2 
          !- Preparation for Making BMesh

          allocate(newp, stat=AllocateStatus)
          if(AllocateStatus /=0) then
             print *, '!!! No enough memory !!!'
             stop
          endif
          nullify(newp%octPrt)
          nullify(newp%octNb1) ; nullify(newp%octNb2)
          nullify(newp%octNb3) ; nullify(newp%octNb4)
          nullify(newp%octNb5) ; nullify(newp%octNb6)
          nullify(newp%octCh1) ; nullify(newp%octCh2)
          nullify(newp%octCh3) ; nullify(newp%octCh4)
          nullify(newp%octCh5) ; nullify(newp%octCh6)
          nullify(newp%octCh7) ; nullify(newp%octCh8)
          nullify(newp%Psort)
          nullify(newp%ptcl)   ; nullify(newp%ptclA) 

          octP  = 0
          octLv = -1
          iFLG  = 0
          iFLG(1) = 4
          Csort = 0
          F = 0.d0
          C = 0.d0
          Z = 0.d0
          G = 0.d0
          D=0.d0
          O=0.d0
          iC= 0
          ptcl_loops=0
          prc_bndry=0

          if(icon==0)then
             sID= rank   * Nall_8
             eID=(rank+1)* Nall_8 -1
          else
             sID=MinMn
             eID=MaxMn
          endif

!!$          if(debugMode>=3)then
!!$             print *,"setting BMesh sID=",sID,"eID=",eID,"MnDisp=",MnDisp,rank
!!$             print *,"MinID(1,-1)=",MinID(1,-1),"MaxID(1,-1)=",MaxID(1,-1),rank
!!$          endif

          do i=sID,eID

             MrtN=i
             !get Position
             call get_InverseIndex(MrtN,iPOS,1)
             !call Inverse_MortonN(MrtN,iPOS)

             ix=iPOS(1)
             iy=iPOS(2)
             iz=iPOS(3)

             !get hierarchy position
             iPOS(1) = 2*ix*LvMax2 - LvMax2
             iPOS(2) = 2*iy*LvMax2 - LvMax2
             iPOS(3) = 2*iz*LvMax2 - LvMax2

             !get real position
             rPOS(1) = real(ix)*dx(1)*2 - dx(1)
             rPOS(2) = real(iy)*dx(2)*2 - dx(2)
             rPOS(3) = real(iz)*dx(3)*2 - dx(3)

             index = MrtN-MnDisp
             octN  = index
             
             Mesh(index)=     &
                  oct(octN,octLv,octP,Csort,iFLG,     &
                  iPOS,rPOS,                          &
                  MrtN,iC,                     &
                  prc_bndry,                    &
                  -1,   &
                  F,C,Z,G,D,O,ptcl_loops,       & 
                  newp%octPrt,                        &
                  newp%octNb1, newp%octNb2,           &
                  newp%octNb3, newp%octNb4,           &
                  newp%octNb5, newp%octNb6,           &
                  newp%octCh1, newp%octCh2,           &
                  newp%octCh3, newp%octCh4,           &
                  newp%octCh5, newp%octCh6,           &
                  newp%octCh7, newp%octCh8,           &
                  newp%Psort ,                        &
                  newp%ptcl  , newp%ptclA)

          enddo

          deallocate(newp)

          return

        end subroutine set_BMesh
        !**********************************************************************

        subroutine set_Mesh
          ! +-------------------------------------------------------------------+
          ! |                                                                   |
          ! | Set & create initial Mesh for  Lv=0                               |
          ! |                                                                   |
          ! +-------------------------------------------------------------------+
          !***********************************************************************
          use particle_set
          ! ---------------------------------------------------------------
          implicit none
          integer(kind=4)    :: octN,octLv,octP,ich
          integer(kind=4)    :: iFLG(3),Csort
          integer(kind=4)    :: iPOS(3), Nint, Nint2, indexP
          real(kind=8)       :: ixx, iyy, izz, rPOS(3), R(9)
          integer(kind=4)    :: MrtN, ix, iy, iz, n1, n2, Isort
          integer(kind=4)  :: prc_bndry,octType, ptcl_loops
          integer(kind=4) :: index, index2, iii
          !        real(kind=8)       :: ixx,iyy,izz
          integer(kind=4)    :: AllocateStatus
          integer(kind=4)    :: iC(iomp0)
          real(kind=8)       :: F(18),C(6*iomp0),Z(IonSorts),G(1:18),D(1:3),O(1:16)
          type(oct), pointer :: newp, p0
          type(prtcl), pointer :: PrtList
          real(kind=8)    :: rx,ry,rz,rr,mr
          nullify(PrtList)
          ! ---------------------------------------------------------------

          Nall = Nall_ini 

          ! ------Set initial variables of the oct Lv =0--------------------
          !        LvMax2 = 2**LvMax ! for Lv 0

          octLv = 0
          octP  = 0
          iFLG  = 0
          Csort = 0

          Isort = -1
          F = 0.d0
          C = 0.d0
          Z = 0.d0
          G = 0.d0
          D=0.d0
          O=0.d0
          iC= 0
          prc_bndry=0
          octType = 0   !0 for Mesh
          ptcl_loops = 0

          MaxIP(:) = 0
          indexP = 0
          !     iPOS and rPOS for Mesh are determined from BMesh in initial_connection

          allocate(newp,stat=AllocateStatus)
          Nint2 = 2**LvMax
          index2 = MinID(1,0)-1
          do index=MinID(1,-1),MaxID(1,-1)
             p0 => Mesh(index)
             iFLG = 0  
             octLv   = p0%octLv  +1
             Nint  = 2**(LvMax-octLv)
             ! -- each directed child octs
             do ich=1,8
                nullify(newp%octPrt)
                nullify(newp%octNb1) ; nullify(newp%octNb2)
                nullify(newp%octNb3) ; nullify(newp%octNb4)
                nullify(newp%octNb5) ; nullify(newp%octNb6)
                nullify(newp%octCh1) ; nullify(newp%octCh2)
                nullify(newp%octCh3) ; nullify(newp%octCh4)
                nullify(newp%octCh5) ; nullify(newp%octCh6)
                nullify(newp%octCh7) ; nullify(newp%octCh8)
                nullify(newp%Psort)
                nullify(newp%ptcl)   ; nullify(newp%ptclA)

                index2 = index2 + 1
                MrtN   = p0%MrtN

                ! -- Extract position index for child-octs --
                n1 = ich
                iz = int((n1-1)/4) + 1
                n2 = n1 - 4*(iz-1)
                iy = int((n2-1)/2) + 1
                ix = n2 - 2*(iy-1)

                ! -- Position index --

                ix = p0%iPOS(1) + (2*ix-3)*Nint
                iy = p0%iPOS(2) + (2*iy-3)*Nint
                iz = p0%iPOS(3) + (2*iz-3)*Nint
                iPOS(1) = ix
                iPOS(2) = iy
                iPOS(3) = iz

                ! -- Position --
                ixx = real(ix + Nint2)/real(2*Nint2)
                iyy = real(iy + Nint2)/real(2*Nint2)
                izz = real(iz + Nint2)/real(2*Nint2)
                rPOS(1) = ixx*dx(1) - HALF*dx(1)
                rPOS(2) = iyy*dx(2) - HALF*dx(2)
                rPOS(3) = izz*dx(3) - HALF*dx(3)
                R(1:3) = rPOS(:)
                R(4:9) = 0d0


                Csort = ich
                octN = index2
                if(dipoleFlag==1)then
                    rx=rPOS(1)-Dpos(1)
                    ry=rPOS(2)-Dpos(2)
                    rz=rPOS(3)-Dpos(3)
                    rr=sqrt(rx**2+ry**2+rz**2)
                    if(rr.gt.dx(1)*2.0d0) then
                        mr=mx*rx+my*ry+mz*rz
                        D(1)=(-mx/(rr**3)+3*rx*mr/(rr**5))/(4*PI)
                        D(2)=(-my/(rr**3)+3*ry*mr/(rr**5))/(4*PI)
                        D(3)=(-mz/(rr**3)+3*rz*mr/(rr**5))/(4*PI)
                    endif
                end if
!                print *,"rPOS=",rPOS(1:3)," D=",D(1:3)
                ! -- Connect parent oct & next oct --
                newp%octPrt => p0
                nullify(newp%Psort)
                !
                ! -- representative particle 
                indexP = indexP + 1

                Pesh(indexP,0) =                  &
                     prtcl(R ,Isort,octN,indexP,PrtList,0)
                newp%ptcl => Pesh(indexP,0)

                ! -- Input informaions into a target variable Mesh(:) --
                Mesh(index2) =                        &
                     oct(octN,octLv,octP,Csort,iFLG,     &
                     iPOS,rPOS,                       &
                     MrtN,iC,                         &
                     prc_bndry, octType,            &
                     F,C,Z,G,D,O,ptcl_loops,    & 
                     newp%octPrt,                     &
                     newp%octNb1, newp%octNb2,        &
                     newp%octNb3, newp%octNb4,        &
                     newp%octNb5, newp%octNb6,        &
                     newp%octCh1, newp%octCh2,        &
                     newp%octCh3, newp%octCh4,        &
                     newp%octCh5, newp%octCh6,        &
                     newp%octCh7, newp%octCh8,        &
                     newp%Psort ,                     &
                     newp%ptcl  , newp%ptclA)

                ! -- Connect eight child-octs of parent oct --
                if(ich==1) then
                   p0%octCh1 => Mesh(index2)
                elseif(ich==2) then
                   p0%octCh2 => Mesh(index2)
                elseif(ich==3) then
                   p0%octCh3 => Mesh(index2)
                elseif(ich==4) then
                   p0%octCh4 => Mesh(index2)
                elseif(ich==5) then
                   p0%octCh5 => Mesh(index2)
                elseif(ich==6) then
                   p0%octCh6 => Mesh(index2)
                elseif(ich==7) then
                   p0%octCh7 => Mesh(index2)
                elseif(ich==8) then 
                   p0%octCh8 => Mesh(index2)
                endif
                !
             enddo  ! for do ich=1,8
          end do
          !--------------------------
          deallocate(newp,stat=iii)
          MaxIP(0) = indexP


          return
          !               
        end subroutine set_Mesh

        !***********************************************************************
        subroutine initial_connection(icon)
          ! +-------------------------------------------------------------------+
          ! |                       
          ! ! 1)Level starts with -1. Level -1 mesh is used for computing current
          ! ! density where neighbouring cells are required for computing
          ! ! the momentum.
          ! !
          ! ! 2)There are Nall_ini/8 grids in Lv=-1. 
          ! ! index for Level -1 : 1,2,....Nall_ini/8
          ! ! index for Level 0 : Nall_ini/8+1, Nall_ini/8+2,...,Nall_ini/8+Nall_ini
          ! !
          ! ! 3)Index and Position index are different,
          ! !   but corresponds one-to-one
          ! !                                           
          ! | Create initial pointer connection                                 |
          ! |                                                                   |
          ! |  Order of neighbor index for octNb#1-6 :                          |
          ! |                                                                   |
          ! |   Y                            Z                                  |
          ! |   |                            |                                  |
          ! |   |          octNb4            |          octNb6                  |
          ! |   +---X        |               +---X        |                     |
          ! |                |                            |                     |
          ! |    octNb1 -- octN -- octNb2     octNb1 -- octN -- octNb2          |
          ! |                |                            |                     |
          ! |                |                            |                     |
          ! |              octNb3                       octNb5                  |
          ! |                                                                   |
          ! +-------------------------------------------------------------------+
          !***********************************************************************
          ! --------------------------------------
          implicit none
          type(oct), pointer :: p0
          integer(kind=4) :: index
          integer(kind=4) :: nbindx(6), Mnarr(27)
          integer(kind=4) :: iPOS(3), iPOS2(3)
          integer(kind=4),intent(in) :: icon
          ! --------------------------------------
          !
          !       Initializing LV=-1 Mesh. There are Nall_ini/8 grids in LV=-1.
          !
          !if(debugMode>=1)print *,"initial_connection start icon=",icon,rank

!!$          if(debugMode>=3)then
!!$             do i=0,nprocs-1
!!$                if(rank==i)then
!!$                   print *,"MnMin =",MinMn,rank
!!$                   print *,"MnMax =",MaxMn,rank
!!$                   print *,"MnDisp=",MnDisp,rank
!!$                endif
!!$                call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$             enddo
!!$          endif

          if(icon==0)then !for normal initialize
             do index=MinID(1,-1),MaxID(1,-1) 
                p0 => Mesh(index)
                !
                iPOS = p0% iPOS

                !Get Morton number of neighbours
                call get_neighbourIndex(iPOS,Mnarr)
                iPOS2=(iPOS+LvMax2)/(2*LvMax2)
                nbindx(1:6)=Mnarr(1:6) - MnDisp

                !Note:: See Morton_neighbours on why octNbi don't correspond to Mnarr(i)          
                !  -- for -x neighbor --
                if(mod(iPOS2(1)-1,NXB/2) /= 0) then
                   p0%octNb1 => Mesh(nbindx(1))        
                endif

                !  -- for +x neighbor --
                if(mod(iPOS2(1),NXB/2) /= 0) then
                   p0%octNb2 => Mesh(nbindx(4))
                endif

                !  -- for -y neighbor --
                if(mod(iPOS2(2)-1,NYB/2) /= 0) then
                   p0%octNb3 => Mesh(nbindx(2))
                endif

                !  -- for +y neighbor --
                if(mod(iPOS2(2),NYB/2) /= 0) then
                   p0%octNb4 => Mesh(nbindx(5))
                endif

                !  -- for -z neighbor --
                if(mod(iPOS2(3)-1,NZB/2) /= 0) then
                   p0%octNb5 => Mesh(nbindx(3))
                endif

                !  -- for +z neighbor --
                if(mod(iPOS2(3),NZB/2) /= 0) then
                   p0%octNb6 => Mesh(nbindx(6))
                endif

                p0%Psort => p0
             enddo
          elseif(icon==1)then !for DDD
             do index=MinID(1,-1),MaxID(1,-1)
                p0 => Mesh(index)
                !
                iPOS = p0% iPOS

                !Get Morton number of neighbours
                call get_neighbourIndex(iPOS,Mnarr)
                iPOS2=(iPOS+LvMax2)/(2*LvMax2)
                nbindx(1:6)=Mnarr(1:6) - MnDisp

                !Note:: See Morton_neighbours on why octNbi don't correspond to Mnarr(i)          
                !  -- for -x neighbor --
                if(nbindx(1)>0 .and. nbindx(1)<=MaxID(1,-1)) then
                   p0%octNb1 => Mesh(nbindx(1))        
                endif

                !  -- for +x neighbor --
                if(nbindx(4)>0 .and. nbindx(4)<=MaxID(1,-1)) then
                   p0%octNb2 => Mesh(nbindx(4))
                endif

                !  -- for -y neighbor --
                if(nbindx(2)>0 .and. nbindx(2)<=MaxID(1,-1)) then
                   p0%octNb3 => Mesh(nbindx(2))
                endif

                !  -- for +y neighbor --
                if(nbindx(5)>0 .and. nbindx(5)<=MaxID(1,-1)) then
                   p0%octNb4 => Mesh(nbindx(5))
                endif

                !  -- for -z neighbor --
                if(nbindx(3)>0 .and. nbindx(3)<=MaxID(1,-1)) then
                   p0%octNb5 => Mesh(nbindx(3))
                endif

                !  -- for +z neighbor --
                if(nbindx(6)>0 .and. nbindx(6)<=MaxID(1,-1)) then
                   p0%octNb6 => Mesh(nbindx(6))
                endif

                p0%Psort => p0
             enddo
          endif


          do index = MinID(1,-1), MaxID(1,-1)

             p0 => Mesh(index)
             ! - for octNb1 -
             if(associated(p0%octNb1))then
                p0%octCh1%octNb1 => p0%octNb1%octCh2     
                p0%octCh3%octNb1 => p0%octNb1%octCh4   
                p0%octCh5%octNb1 => p0%octNb1%octCh6 
                p0%octCh7%octNb1 => p0%octNb1%octCh8
             endif
             ! - for octNb2 -                
             if(associated(p0%octNb2))then
                p0%octCh2%octNb2 => p0%octNb2%octCh1
                p0%octCh4%octNb2 => p0%octNb2%octCh3
                p0%octCh6%octNb2 => p0%octNb2%octCh5
                p0%octCh8%octNb2 => p0%octNb2%octCh7
             endif
             ! - for octNb3 -
             if(associated(p0%octNb3))then
                p0%octCh1%octNb3 => p0%octNb3%octCh3
                p0%octCh2%octNb3 => p0%octNb3%octCh4
                p0%octCh5%octNb3 => p0%octNb3%octCh7
                p0%octCh6%octNb3 => p0%octNb3%octCh8
             endif
             ! - for octNb4 -
             if(associated(p0%octNb4))then
                p0%octCh3%octNb4 => p0%octNb4%octCh1
                p0%octCh4%octNb4 => p0%octNb4%octCh2
                p0%octCh7%octNb4 => p0%octNb4%octCh5
                p0%octCh8%octNb4 => p0%octNb4%octCh6
             endif
             ! - for octNb5 -
             if(associated(p0%octNb5))then
                p0%octCh1%octNb5 => p0%octNb5%octCh5
                p0%octCh2%octNb5 => p0%octNb5%octCh6
                p0%octCh3%octNb5 => p0%octNb5%octCh7
                p0%octCh4%octNb5 => p0%octNb5%octCh8
             endif
             ! - for octNb6 -
             if(associated(p0%octNb6))then
                p0%octCh5%octNb6 => p0%octNb6%octCh1
                p0%octCh6%octNb6 => p0%octNb6%octCh2
                p0%octCh7%octNb6 => p0%octNb6%octCh3
                p0%octCh8%octNb6 => p0%octNb6%octCh4
             endif
             ! - for octCh1 -
             p0%octCh1%octNb2 => p0       %octCh2
             p0%octCh1%octNb4 => p0       %octCh3
             p0%octCh1%octNb6 => p0       %octCh5
             ! - for octCh2 -
             p0%octCh2%octNb1 => p0       %octCh1
             p0%octCh2%octNb4 => p0       %octCh4
             p0%octCh2%octNb6 => p0       %octCh6
             ! - for octCh3 -
             p0%octCh3%octNb2 => p0       %octCh4
             p0%octCh3%octNb3 => p0       %octCh1
             p0%octCh3%octNb6 => p0       %octCh7
             ! - for octCh4 -
             p0%octCh4%octNb1 => p0       %octCh3
             p0%octCh4%octNb3 => p0       %octCh2
             p0%octCh4%octNb6 => p0       %octCh8
             ! - for octCh5 -
             p0%octCh5%octNb2 => p0       %octCh6
             p0%octCh5%octNb4 => p0       %octCh7
             p0%octCh5%octNb5 => p0       %octCh1
             ! - for octCh6 -
             p0%octCh6%octNb1 => p0       %octCh5           
             p0%octCh6%octNb4 => p0       %octCh8
             p0%octCh6%octNb5 => p0       %octCh2
             ! - for octCh7 -
             p0%octCh7%octNb2 => p0       %octCh8
             p0%octCh7%octNb3 => p0       %octCh5
             p0%octCh7%octNb5 => p0       %octCh3
             ! - for octCh8 -
             p0%octCh8%octNb1 => p0       %octCh7
             p0%octCh8%octNb3 => p0       %octCh6
             p0%octCh8%octNb5 => p0       %octCh4

             !p0%iFLG(2)=1
          enddo

          return
          !
        end subroutine initial_connection
        !
        !***********************************************************************

      end subroutine initial_setting
!
!***********************************************************************
!!$      subroutine restart_setting
!!$! +-------------------------------------------------------------------+
!!$! +    set up initial parameters and initial conditions               +
!!$! +-------------------------------------------------------------------+
!!$!***********************************************************************
!!$        use param
!!$        use init_condition
!!$        use const
!!$        implicit none
!!$        integer(kind=4) :: i
!!$
!!$! -- Preparation --
!!$        call preparation
!!$        call CFL_check
!!$!
!!$! -- Initial setting for octs --
!!$!        call set_oct
!!$!        call initial_connection
!!$!        call periodic_boundary
!!$!
!!$! -- initial setting for time evolution
!!$        call set_evolution
!!$!
!!$        do i=1,3
!!$          dtdx(i)  = dt/dx(i)
!!$          dt2dx(i) = HALF*dt/dx(i)
!!$          odx(i)   = ONE/dx(i)
!!$        end do
!!$!
!!$! -- criterion setting
!!$        if(LvMax.eq.0) crtr0(0)=-999.d0
!!$        if(LvMax.eq.1) crtr0(0)=-999.d0 ! 0.40d0
!!$        if(LvMax.eq.2) then 
!!$           crtr0(0)=-999.d0 ; crtr0(1)=-999.d0
!!$        end if
!!$!
!!$        return
!!$        end subroutine restart_setting
!
!***********************************************************************
      subroutine refine_oct
! +--------------------------------------------------------------------+
! +     flagging for refinement and addition of refined new octs       +
! +                                                                    +
! +                   lower level : <present level>    : upper level   +
! +                                                                    +
! +     iFLG : <= -4      O       :       X            :      X        +
! +     iFLG :  = -3      O       =>      I  (edge)    :      X        +
! +     iFLG :  = -2      O       :       T  (overlap) :      X        +
! +     iFLG :  = -1      O       :       T  (overlap) :      X        +
! +     iFLG :  =  0      A      <=       O            :      X        +
! +     iFLG :  =  1      A       :       O            =>     I        +
! +     iFLG :  =  2      A       :       O            :      T        +
! +     iFLG :  =  3      A       :       O            :      T        +
! +     iFLG : >=  4      A       :       A  (refined)<=      O        +
! +                                                                    +
! +    condition of the oct                                            +
! +                                                                    +
! +          O : time evolution is computed in the level               +
! +          T : temporarily computed in the level                     +
! +          X : does not exsist or will be delated                    +
! +          A : averaged field data from child octs                   +
! +          I : interpolated field data from parent octs              +
! +                                                                    +
! +    information exchange                                            +
! +                                                                    +
! +          => : interpolation                                        +
! +          <= : averaging                                            +
! +                                                                    +
! +--------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
! ------------------------------------------
        implicit none
        integer(kind=4) :: iLv
!
        

        !if(debugMode>=1)print*,rank,'entering refine_oct', istep

!!$        if(debugMode>=3)then
!!$           print *,"debug call before refine_oct",rank
!!$           call output_octs_and_family(MinID(1,-1),MaxID(1,-1),LvMax,0,1,0)
!!$           call output_octs_and_family(MinID(3,-1),MaxID(3,-1),LvMax,1,1,0)
!!$           print *,"debug call before refine_oct ok",rank
!!$        endif

        !call reset_Gstate

        do iLv=0,LvMax-1
           call set_flag_block(iLv)
           call add_oct(iLv)
           call add_Goct(iLv)
        end do

        call set_flag_block(LvMax)

        do iLv=0,LvMax-1
           call connect_oct(iLv)
           call connect_Goct(iLv)
        enddo
    
        call connect_oct(LvMax)
        call connect_Goct(LvMax)

        do iLv=0,LvMax-1 !set_buffer_aveG is necessary. because ave is not work correct
                         !at process boundary
           call set_buffer_ave(iLv,4,6) ! 01234444443210 -> 01234666643210
           call set_buffer_in( iLv,5,4) ! 01234666643210 -> 01234566543210
           call refresh_iFLG(iLv)
        end do

!!$        if(LvMax>=1)then
!!$           do iLv=LvMax,1,-1
!!$              call set_flag_prt(iLv)
!!$              call set_flag_prtG(iLv)
!!$           end do
!!$        endif

        !if(debugMode>=1)print*,rank,'exiting refine_oct'

!!$        if(debugMode>=3)then
!!$           print *,"debug call after refine_oct",rank
!!$           call output_octs_and_family(MinID(1,-1),MaxID(1,-1),LvMax,0,1,0)
!!$           call output_octs_and_family(MinID(3,-1),MaxID(3,-1),LvMax,1,1,0)
!!$           print *,"debug call after refine_oct ok",rank
!!$        endif

        return
      end subroutine refine_oct

!2012/04/18 yagi added
      subroutine refine_oct_DDD
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
! ------------------------------------------
        implicit none
        integer(kind=4) :: iLv,i,ii,index
        type(oct),pointer::p0

        !if(debugMode>=3)print *,"enter refine_oct_DDD",rank

        !remake hierarchy Goct
        do iLv=0,LvMax-1
           if(minID(3,iLv)<MaxID(3,iLv))then
              do index=MinID(3,iLv),MaxID(3,iLv)
                 p0=>Mesh(index)
                 p0%iFLG(1)=min(p0%iFLG(1),4)
                 p0%iFLG(2)=p0%iFLG(1)
                 p0%iFLG(1)=0
              enddo
           endif

           if(minID(4,iLv)<MaxID(4,iLv))then
              do index=MinID(4,iLv),MaxID(4,iLv)
                 p0=>Mesh(index)
                 p0%iFLG(1)=min(p0%iFLG(1),4)
                 p0%iFLG(2)=p0%iFLG(1)
                 p0%iFLG(1)=0
              enddo
           endif

           call refresh_iFLG(iLv)

           do i=3,4
              if(minID(i,iLv).lt.maxID(i,iLv)) then 
                 do index=minID(i,iLv),maxID(i,iLv)
                    p0 => Mesh(index)
                    ii=p0%iFLG(2)
                    p0%iFLG(2)=0
                    if(ii>=4) then 
                       if((p0%iFLG(1)<=3).and.(p0%iFLG(1)>=0)) p0%iFLG(2)=-1
                    else if(ii<=0) then
                       if(p0%iFLG(1)>=1)p0%iFLG(2)=2
                    else if(ii<=3) then 
                       if(p0%iFLG(1)>=1) p0%iFLG(2)=1
                    endif
                 end do
              end if
           end do

           call add_Goct(iLv)

           call connect_oct(iLv)
           call connect_Goct(iLv)

        enddo

!!$        if(debugMode>=3)print *,"exit refine_oct_DDD",rank

      end subroutine refine_oct_DDD

!2011/09/13 yagi added
!============================================
      recursive subroutine pack_MeshFLG(index,sendLv,count,iLv,icon)
        use param
        use oct_set
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4)::count,iLv
        integer(kind=4),intent(in)::index,sendLv,icon
        type(oct),pointer::p0

        p0=>Mesh(index)

        if(p0%iFLG(1)>0)then
           if(iLv==sendLv)then
              sbuf_BFI(count*3-2,icon)=p0%octCh1%iFLG(3)
              sbuf_BFI(count*3-1,icon)=p0%octCh1%iFLG(1)
              sbuf_BFI(count*3  ,icon)=p0%octCh1%MrtN
              count=count+1
              sbuf_BFI(count*3-2,icon)=p0%octCh2%iFLG(3)
              sbuf_BFI(count*3-1,icon)=p0%octCh2%iFLG(1)
              sbuf_BFI(count*3  ,icon)=p0%octCh2%MrtN
              count=count+1
              sbuf_BFI(count*3-2,icon)=p0%octCh3%iFLG(3)
              sbuf_BFI(count*3-1,icon)=p0%octCh3%iFLG(1)
              sbuf_BFI(count*3  ,icon)=p0%octCh3%MrtN
              count=count+1
              sbuf_BFI(count*3-2,icon)=p0%octCh4%iFLG(3)
              sbuf_BFI(count*3-1,icon)=p0%octCh4%iFLG(1)
              sbuf_BFI(count*3  ,icon)=p0%octCh4%MrtN
              count=count+1
              sbuf_BFI(count*3-2,icon)=p0%octCh5%iFLG(3)
              sbuf_BFI(count*3-1,icon)=p0%octCh5%iFLG(1)
              sbuf_BFI(count*3  ,icon)=p0%octCh5%MrtN
              count=count+1
              sbuf_BFI(count*3-2,icon)=p0%octCh6%iFLG(3)
              sbuf_BFI(count*3-1,icon)=p0%octCh6%iFLG(1)
              sbuf_BFI(count*3  ,icon)=p0%octCh6%MrtN
              count=count+1
              sbuf_BFI(count*3-2,icon)=p0%octCh7%iFLG(3)
              sbuf_BFI(count*3-1,icon)=p0%octCh7%iFLG(1)
              sbuf_BFI(count*3  ,icon)=p0%octCh7%MrtN
              count=count+1
              sbuf_BFI(count*3-2,icon)=p0%octCh8%iFLG(3)
              sbuf_BFI(count*3-1,icon)=p0%octCh8%iFLG(1)
              sbuf_BFI(count*3  ,icon)=p0%octCh8%MrtN
              count=count+1
           else
              call pack_MeshFLG(p0%octCh1%octN,sendLv,count,iLv+1,icon)
              call pack_MeshFLG(p0%octCh2%octN,sendLv,count,iLv+1,icon)
              call pack_MeshFLG(p0%octCh3%octN,sendLv,count,iLv+1,icon)
              call pack_MeshFLG(p0%octCh4%octN,sendLv,count,iLv+1,icon)
              call pack_MeshFLG(p0%octCh5%octN,sendLv,count,iLv+1,icon)
              call pack_MeshFLG(p0%octCh6%octN,sendLv,count,iLv+1,icon)
              call pack_MeshFLG(p0%octCh7%octN,sendLv,count,iLv+1,icon)
              call pack_MeshFLG(p0%octCh8%octN,sendLv,count,iLv+1,icon) 
           endif
        endif

      end subroutine pack_MeshFLG
!-----------------------------------------------
      subroutine pack_BMeshFLGT(sdest,sendLv,icon)
        use param
        use oct_set
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4),intent(in)::sdest,sendLv,icon
        integer(kind=4)::i,iLv,index,sID,eID,count
        
        if(sendLv==-1)return
        if(n_rcells_proc(sdest+1)==0)return
        if(MinID(1,sendLv)>=MaxID(1,sendLv))return
        
        !set loop index
        if(sdest==0)then
           sID=1
        else
           sID=sum(n_rcells_proc(1:sdest))+1
        endif
        eID=sum(n_rcells_proc(1:sdest+1))
        
        count=1
        iLv  =-1
        
        do i=sID,eID
           index=iBMesh_arr(i)
           call pack_MeshFLG(index,sendLv,count,iLv+1,icon)
        enddo

        send_onum(icon)=count-1
        !if(istep==30)print *,"send_onum=",send_onum(icon),"sdest=",sdest,"rank=",rank
      end subroutine pack_BMeshFLGT
!============================================
      recursive subroutine set_received_GMeshFLG(index,sendLv,count,iLv,icon)
        use param
        use oct_set
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4),intent(in)::index,sendLv,icon
        integer(kind=4)::count,iLv
        type(oct),pointer::p0

        p0=>Mesh(index)
        if(p0%iFLG(1)>0)then
           if(iLv==sendLv)then
              !do not touch iFLG(2) here
              p0%octCh1%iFLG(1)=rbuf_GFI(count*3-1,icon)
              p0%octCh1%iFLG(3)=rbuf_GFI(count*3-2,icon)
              count=count+1
              p0%octCh2%iFLG(1)=rbuf_GFI(count*3-1,icon)
              p0%octCh2%iFLG(3)=rbuf_GFI(count*3-2,icon)
              count=count+1
              p0%octCh3%iFLG(1)=rbuf_GFI(count*3-1,icon)
              p0%octCh3%iFLG(3)=rbuf_GFI(count*3-2,icon)
              count=count+1
              p0%octCh4%iFLG(1)=rbuf_GFI(count*3-1,icon)
              p0%octCh4%iFLG(3)=rbuf_GFI(count*3-2,icon)
              count=count+1
              p0%octCh5%iFLG(1)=rbuf_GFI(count*3-1,icon)
              p0%octCh5%iFLG(3)=rbuf_GFI(count*3-2,icon)
              count=count+1
              p0%octCh6%iFLG(1)=rbuf_GFI(count*3-1,icon)
              p0%octCh6%iFLG(3)=rbuf_GFI(count*3-2,icon)
              count=count+1
              p0%octCh7%iFLG(1)=rbuf_GFI(count*3-1,icon)
              p0%octCh7%iFLG(3)=rbuf_GFI(count*3-2,icon)
              count=count+1
              p0%octCh8%iFLG(1)=rbuf_GFI(count*3-1,icon)
              p0%octCh8%iFLG(3)=rbuf_GFI(count*3-2,icon)
              count=count+1

           else
              call set_received_GMeshFLG(p0%octCh1%octN,sendLv,count,iLv+1,icon)
              call set_received_GMeshFLG(p0%octCh2%octN,sendLv,count,iLv+1,icon)
              call set_received_GMeshFLG(p0%octCh3%octN,sendLv,count,iLv+1,icon)
              call set_received_GMeshFLG(p0%octCh4%octN,sendLv,count,iLv+1,icon)
              call set_received_GMeshFLG(p0%octCh5%octN,sendLv,count,iLv+1,icon)
              call set_received_GMeshFLG(p0%octCh6%octN,sendLv,count,iLv+1,icon)
              call set_received_GMeshFLG(p0%octCh7%octN,sendLv,count,iLv+1,icon)
              call set_received_GMeshFLG(p0%octCh8%octN,sendLv,count,iLv+1,icon)     
           endif
        endif

      end subroutine set_received_GMeshFLG
!--------------------------------------------
      subroutine set_received_GMeshFLGT(sendLv,recvnum,icon)
        use message_passing_interface
        use oct_set
        use param
        use init_mesh_size
        implicit none
        integer(kind=4),intent(in)::sendLv,recvnum,icon
        integer(kind=4)::count,index,iLv

        if(sendLv==-1)return
        !this statement includes BUG when refresh_iFLG called at refine_oct_DDD
        !if(MinID(3,sendLv)>=MaxID(3,sendLv))return 

        if(MinID(3,sendLv)>=MaxID(3,sendLv) .and. MinID(4,sendLv)>=MaxID(4,sendLv))return
        if(recvnum<=0)return

        count=1
        iLv  =-1

        do while(count<=recvnum)
           index=GMn2octN(rbuf_GFI(count*3,icon))
           call set_received_GMeshFLG(index,sendLv,count,iLv+1,icon)
        enddo

      end subroutine set_received_GMeshFLGT

!===============================================================
      subroutine refresh_iFLG(sendLv)
        use message_passing_interface
        use param
        use oct_set
        implicit none
        integer(kind=4),intent(in)::sendLv
        integer(kind=4)::i,j,sdest(2),rdest(2),loopnum
        integer(kind=4)::req(2,2),istat(MPI_STATUS_SIZE,2,2)
        integer(kind=4)::stream1,stream2,streamtemp
        !if(debugMode>=1)print *,"entering refresh_iFLG",rank
    
        loopnum=rglistnum-1
        stream1=1
        stream2=2
        j=0

        !-----------stream=1 op=1 exchange onum
        j=j+1
        sdest    (stream1)=send_rlist(j)
        rdest    (stream1)=recv_glist(j)
        send_onum(stream1)=0
        recv_onum(stream1)=0
        if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_irecv (recv_onum(stream1),1,MPI_INTEGER,rdest(stream1),1,&
        &MPI_COMM_WORLD,req(1,stream1),ierr)
        call pack_BMeshFLGT(sdest(stream1),sendLv,stream1)
        if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_isend (send_onum(stream1),1,MPI_INTEGER,sdest(stream1),1,&
        &MPI_COMM_WORLD,req(2,stream1),ierr)
        if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_wait      (req(1,stream1),          istat(1,1,stream1),ierr)
        if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_wait      (req(2,stream1),          istat(1,2,stream1),ierr)
        !-----------stream=1 op=2 prepare to exchange oct
        if(recv_onum(stream1)>0)call mpi_irecv  (rBuf_GFI(:,stream1),recv_onum(stream1)*3,MPI_INTEGER,rdest(stream1),3,&
        &MPI_COMM_WORLD,req(1,stream1),ierr) 
        if(send_onum(stream1)>0)call mpi_isend  (sbuf_BFI(:,stream1),send_onum(stream1)*3,MPI_INTEGER,sdest(stream1),3,&
        &MPI_COMM_WORLD,req(2,stream1),ierr)

        do i=1,loopnum
           streamtemp=stream1
           stream1=stream2
           stream2=streamtemp
           !-----------stream=1 op=1 exchange onum
           j=j+1
           sdest    (stream1)=send_rlist(j)
           rdest    (stream1)=recv_glist(j)
           send_onum(stream1)=0
           recv_onum(stream1)=0
           if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_irecv (recv_onum(stream1),1,MPI_INTEGER,rdest(stream1),1,&
           &MPI_COMM_WORLD,req(1,stream1),ierr)
           call pack_BMeshFLGT(sdest(stream1),sendLv,stream1)
           if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_isend (send_onum(stream1),1,MPI_INTEGER,sdest(stream1),1,&
           &MPI_COMM_WORLD,req(2,stream1),ierr)
           if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_wait      (req(1,stream1),          istat(1,1,stream1),ierr)
           if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_wait      (req(2,stream1),          istat(1,2,stream1),ierr)
           !-----------stream=1 op=2 prepare to exchange oct
           if(recv_onum(stream1)>0)call mpi_irecv  (rBuf_GFI(:,stream1),recv_onum(stream1)*3,MPI_INTEGER,rdest(stream1),3,&
           &MPI_COMM_WORLD,req(1,stream1),ierr) 
           if(send_onum(stream1)>0)call mpi_isend  (sbuf_BFI(:,stream1),send_onum(stream1)*3,MPI_INTEGER,sdest(stream1),3,&
           &MPI_COMM_WORLD,req(2,stream1),ierr)
           !===========stream=2 op=3 complete to exchange oct
           if(recv_onum(stream2)>0)call mpi_wait      (req(1,stream2),istat(1,1,stream2),ierr)
           if(send_onum(stream2)>0)call mpi_wait      (req(2,stream2),istat(1,2,stream2),ierr)
           call set_received_GMeshFLGT(sendLv,recv_onum(stream2),stream2)      
        enddo

        !-----------stream=1 op=3 complete to exchange oct
        if(recv_onum(stream1)>0)call mpi_wait      (req(1,stream1),istat(1,1,stream1),ierr)
        if(send_onum(stream1)>0)call mpi_wait      (req(2,stream1),istat(1,2,stream1),ierr)
        call set_received_GMeshFLGT(sendLv,recv_onum(stream1),stream1)   

      end subroutine refresh_iFLG
!============================================
      subroutine reset_Gstate
        use param
        use init_mesh_size
        use oct_set
        implicit none
        integer(kind=4)::index
        type(oct),pointer::p0

        do index=MinID(3,-1),MaxID(3,-1)
           p0=>Mesh(index)
           p0%iFLG=(/4,0,0/)
        end do
        do index=MinID(3,0),MaxID(3,0)
           p0=>Mesh(index)
           p0%iFLG=(/0,0,0/)
        enddo
        if(LvMax>0)then
           do index=MinID(3,1),MaxID(3,LvMax)
              p0=>Mesh(index)
              p0%iFLG=(/-4,-4,0/)
           enddo
        endif


      end subroutine reset_Gstate


!============================================
!additional part end

!**********************************************************************
      subroutine advance_fieldT(iLv,ii)
! +-------------------------------------------------------------------+
! |                                                                   |
! | Compute next step for fields E,B with leap-frog scheme            |
! | (FDTD, Finite Differences Time Domain)                            |
! | [K. S. Yee, Numerical Solution of Initial Boundary Value          |
! | Problems Involving Maxwells Equations in Isotropic Media,         |
! | IEEE Trans. Antennas Prop., 14, 302 (1966)]                       |
! |                                                                   |
! |                                                                   |
! |   THE FIRST STEP:    B[n+1/2] - B[n]   = -dt/2 * ROT E[n]         |
! |        or                                                         |
! |   THE THIRD STEP:    B[n+1] - B[n+1/2] = -dt/2 * ROT E[n+1]       |
! |                                                                   |
! |  -- HALF STEP for B (implicit averaging of B over time step) --   |
! |         B[n+1/2] - B[n]     = -dt/2 * ROT E[n]                    |
! |              or                                                   |
! |         B[n+1]   - B[n+1/2] = -dt/2 * ROT E[n+1]                  |
! |                                                                   |
! |  --------------------------------------------------------------   |
! |                                                                   |
! |      iLv ; refinement level                                       |
! |                                                                   |
! |      ii = 1 ; first half of the calcuration                       |
! |         1) B(n) & E(n)     => B(n+1/2)                            |
! |         2)    get averaged field data from the upper level        |
! |         3)    put interpolated field data to the upper level      |
! |         4) E(n) & B(n+1/2) => E(n+1/2)                            |
! |         5)    get averaged field data from the upper level        |
! |         6)    put interpolated field data to the upper level      |
! |         7) E(n+1/2) & B(n+1/2) => E(n+1)                          |
! |                                                                   |
! |      ii = 2 ; second half of the calcuration                      |
! |         1)    get averaged field data from the upper level        |
! |         2)    put interpolated field data to the upper level      |
! |         3) B(n+1/2) & E(n+1)     => B(n+1)                        |
! |         4)    get averaged field data from the upper level        |
! |         5)    put interpolated field data to the upper level      |    
! +-------------------------------------------------------------------+
!***********************************************************************
 !       use mpi
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
!        use work_variable
        implicit none
        integer(kind=4) :: iLv,ii

        !if(debugMode/=0) print*, 'entering Advance_fieldT rank, iLv, irec(iLv)', rank, iLv, irec(ilv)
        if(ii==1) then

!------------ 
           if(irec(iLv).eq.0) then
!Finer levels (iLv>0) go through this if section, then the next else if section as well

              call copy_particleU( iLv-1    , 4,1, 1 )
              call copy_particleU( iLv-1    , 5,1, 1 )
              call copy_particleD( iLv-1    , 3,1, 1 )
              call copy_particleD( iLv-1    , 2,1, 1 )
              
!              call fipp_start()
!====this timing is not influence because reconnect_particle is more important=====
			  if(iLv == 0) then
			  call move_Vparticle
			  call reconnect_particleP(0,-2,5)
			  end if
!===================================================
              call move_particle      (iLv,-2,5)	
!              call advance_field     ( iLv,1.0d0,-2,4,"J")
!              call fipp_stop()
              
              call reconnect_particleP(iLv,-2,5)  
              if(ParticleInjectFlag(1)==1.and.iLv==0) then
                if(BackgroundParticle==1 .and. BeamParticle==1 .and. arraySize.gt.0)then
                    call deleteFieldInSatellite
                    call removeParticleForInjection
                    call particle_injctT_xs                
                else if(BackgroundParticle==1)then
                    call removeParticleForInjection
                    call particle_injctT_xs
                else if(BeamParticle==1 .and. arraySize.gt.0)then
				   if(satellite == 1)then
                    call deleteFieldInSatellite
				   end if
                    call removeParticleForInjection
                    call particle_injctT_xs
                endif
              end if
              call refresh_fields     ("J",iLv)
              call refresh_particles(iLv)

              irec(iLv)=abs(irec(iLv)-1)
           else if(irec(iLv).gt.0) then 
!iLv=0 always goes through only this portion (iLv>=0 as well)
!              call fipp_start()
!====this timing is not influence because reconnect_particle is more important=====
			  if(iLv == 0) then
			  call move_Vparticle
			  call reconnect_particleP(0,0,3)
			  end if
!===================================================
              call move_particle      (iLv,-2,5)
!              call advance_field     ( iLv,1.0d0,-2,4,"J")
!              call fipp_stop()
              
              call reconnect_particleP(iLv, 0,3)
              if(ParticleInjectFlag(1)==1.and.iLv==0) then
                if(BackgroundParticle==1 .and. BeamParticle==1 .and. arraySize.gt.0)then
                    call removeParticleForInjection
                    call particle_injctT_xs                
                else if(BackgroundParticle==1)then
                    call removeParticleForInjection
                    call particle_injctT_xs
                else if(BeamParticle==1 .and. arraySize.gt.0)then
                    call removeParticleForInjection
                    call particle_injctT_xs
                endif
              end if
              call refresh_fields     ("J",iLv)
              call refresh_particles(iLv)
              irec(iLv)=abs(irec(iLv)-1)
           end if
!------------
     !      print *,"start section3 rank=",rank
           if (CIP .eq. 0) then
             call advance_field     ( iLv,1.0d0,-2,4,"B")
             call refresh_Fields    ("B",iLv)
             call     average_fieldB( iLv      , 5,64)
             call interpolate_fieldB( iLv      , 1,1)
             call interpolate_fieldB( iLv      , 2,1)
             call interpolate_fieldB( iLv      , 3,1)
!             call     average_field( iLv      , 5,64,"B")
!             call interpolate_field( iLv      , 1,1,"B")
!             call interpolate_field( iLv      , 2,1,"B")
!             call interpolate_field( iLv      , 3,1,"B")
             call advance_field     ( iLv,0.5d0,-2,4,"E")
             call refresh_Fields    ("E",iLv)
             call interpolate_fieldE( iLv      , 1,1)
             call interpolate_fieldE( iLv      , 2,1)
!             call interpolate_field( iLv      , 1,1,"E")
!             call interpolate_field( iLv      , 2,1,"E")
           else 
             call non_advection(iLv,-2,6,0.25d0)
             call removeBoundaryECIP
             call removeBoundaryBCIP
              !call deleteFieldInSatellite
             call refresh_fields("E",iLv)
             call refresh_fields("G",iLv)
             
             call advection_x(iLv,-2,6,0.5d0)
             call removeBoundaryECIP
             call removeBoundaryBCIP
             !call deleteFieldInSatellite
             call refresh_fields("E",iLv)
             call refresh_fields("B",iLv)
             call refresh_fields("G",iLv)
             !
             call advection_y(iLv,-2,6,0.5d0)
             call removeBoundaryECIP
             call removeBoundaryBCIP
            ! call deleteFieldInSatellite
             call refresh_fields("E",iLv)
             call refresh_fields("B",iLv)
             call refresh_fields("G",iLv)
             !
             call advection_z(iLv,-2,6,0.5d0)
             call removeBoundaryECIP
             call removeBoundaryBCIP
             !call deleteFieldInSatellite
             call refresh_fields("E",iLv)
             call refresh_fields("B",iLv)
             call refresh_fields("G",iLv)
             !
             call non_advection(iLv,-2,6,0.25d0)
             call removeBoundaryECIP
             call removeBoundaryBCIP
             !call deleteFieldInSatellite
             call refresh_fields("E",iLv)
             call refresh_fields("G",iLv)
             !
             call interpolate_field(iLv,1,1,"E")
             call interpolate_field(iLv,1,1,"B")
             call average_field(iLv,5,64,"E")
             call average_field(iLv,5,64,"B")
           endif
        else if(ii==2) then  ! second half of Time advancement
           if (CIP .eq. 0) then
             call advance_field     ( iLv,0.5d0,-2,4,"E")
             call refresh_Fields    ("E",iLv)
             call     average_fieldE( iLv      , 4,4)
             call interpolate_fieldE( iLv      , 1,1)
             call interpolate_fieldE( iLv      , 2,1)
             call interpolate_fieldE( iLv      , 3,1)
!             call     average_field( iLv      , 4,4,"E")
!             call interpolate_field( iLv      , 1,1,"E")
!             call interpolate_field( iLv      , 2,1,"E")
!             call interpolate_field( iLv      , 3,1,"E")
             call advance_field     ( iLv,1.0d0,-2,4,"B")
             call refresh_fields    ("B",iLv)
             call     average_fieldE( iLv      , 4,64)
             call     average_fieldB( iLv      , 4,64)
             call interpolate_fieldE( iLv      , 5,2)
             call interpolate_fieldB( iLv      , 5,2)
!             call     average_field( iLv      , 4,64,"F")
!             call interpolate_field( iLv      , 5,2,"F")
           else
             call non_advection(iLv,-2,4,0.25d0)
             call removeBoundaryECIP
             call removeBoundaryBCIP
             call refresh_fields("E",iLv)
             call refresh_fields("G",iLv)

             call advection_z(iLv,-2,4,0.5d0)
             call removeBoundaryECIP
             call removeBoundaryBCIP
             call refresh_fields("E",iLv)
             call refresh_fields("B",iLv)
             call refresh_fields("G",iLv)
             !
             call advection_y(iLv,-2,4,0.5d0)
             call removeBoundaryECIP
             call removeBoundaryBCIP
             call refresh_fields("E",iLv)
             call refresh_fields("B",iLv)
             call refresh_fields("G",iLv)
             !
             call advection_x(iLv,-2,4,0.5d0)
             call removeBoundaryECIP
             call removeBoundaryBCIP
             call refresh_fields("E",iLv)
             call refresh_fields("B",iLv)
             call refresh_fields("G",iLv)
             !
             call non_advection(iLv,-2,4,0.25d0)
             call removeBoundaryECIP
             call removeBoundaryBCIP
             call refresh_fields("E",iLv)
             !
             call interpolate_field(iLv,5,2,"E")
             call interpolate_field(iLv,5,2,"B")
             call average_field(iLv,4,64,"E")
             call average_field(iLv,4,64,"B")
           endif
        end if

        return
        end subroutine advance_fieldT
!***********************************************************************
      subroutine sort_oct
! +-------------------------------------------------------------------+
! |     sorting / elimination of refined octs in iLv-th level mesh    |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4) :: iLv,ii,ii0,i,ip,indexP,indexQ
        integer(kind=4) :: iiG,iiG0
        integer(kind=4) :: index,nprt,iFLG0
        integer(kind=4) :: minIF(0:Lvmax),maxIF(0:Lvmax),GminIF(0:Lvmax),GmaxIF(0:Lvmax)
        type(oct), pointer :: p0,p1
        type(prtcl), pointer :: pp
!-------------------------------------------
!
        !if(debugMode>=1)print*,rank, 'entering sort_oct'

        if(debugMode>=3)then
           print *,"debug call before sort_oct",rank
           call full_oct_check
           print *,"debug call before sort_oct ok",rank
        endif
        call full_oct_check

        call reset_sortBuffer
 

        do iLv=LvMax-1,0,-1
           call copy_particleU(   iLv  ,-1, 2, 1 )
           call delete_particleT( iLv+1,-4,-4)
!           print *,'delete_particleT_1'
        end do

        ii  =BMeshBound
!        iiG =MaxID(3,-1)
        iiG =minID(3,0)-1   !==

        do iLv=0,LvMax 
           ii0 = ii
           iiG0= iiG
           do i=1,2
              if(maxID(i,iLv).gt.minID(i,iLv)) then
                 do index=MinID(i,iLv),MaxID(i,iLv)
                    p0 => Mesh(index)
                    if(p0%iFLG(1)>=-3)then
                       ii=ii+1 !count the number of remaining cells
                       p0%Psort => Mesh(ii)
                       Mesh2(ii)=Mesh(index) !Remaining cells are kept in Mesh2
                    end if
                 end do
              end if

              if(maxID(i+2,iLv).gt.minID(i+2,iLv)) then
                 do index=MinID(i+2,iLv),MaxID(i+2,iLv)
                    p0 => Mesh(index)
                    if(p0%iFLG(1)>=-3)then
                       iiG=iiG+1 
                       p0%Psort => Mesh(iiG)
                       Mesh2(iiG)=Mesh(index)
                    endif
                 end do
              end if

           end do
        
           if(ii.eq.ii0) then ! This case corresponds to the fact that there is no surviving cells. All of them are eliminated.
              minIF(iLv)=ii0
              maxIF(iLv)=ii0
           else
!!$omp parallel do private(index,p0,p1,nprt,indexP,indexQ,pp,ip) shared(ii0,ii,iLv,Mesh,Mesh2)
              do index=ii0+1,ii ! This case corresponds to the situation where surviving cells exist.

                 p0 => Mesh2(index)
                 p0%octN = index

                 p0%octNb1 => p0%octNb1%Psort
                 p0%octNb2 => p0%octNb2%Psort
                 p0%octNb3 => p0%octNb3%Psort
                 p0%octNb4 => p0%octNb4%Psort
                 p0%octNb5 => p0%octNb5%Psort
                 p0%octNb6 => p0%octNb6%Psort

                 if(iLv.gt.-1) then
                    p1 => p0%octPrt
                    nprt = p1%octN
                  
                    p1 => Mesh2(nprt)
               
                    if(p0%Csort<=4) then 
                       if(p0%Csort==1) then 
                          p1%octCh1 => p0%Psort
                       else if(p0%Csort==2) then 
                          p1%octCh2 => p0%Psort
                       else if(p0%Csort==3) then 
                          p1%octCh3 => p0%Psort
                       else 
                          p1%octCh4 => p0%Psort
                       end if
                    else
                       if(p0%Csort==5) then 
                          p1%octCh5 => p0%Psort
                       else if(p0%Csort==6) then 
                          p1%octCh6 => p0%Psort
                       else if(p0%Csort==7) then 
                          p1%octCh7 => p0%Psort
                       else 
                          p1%octCh8 => p0%Psort
                       end if
                    end if
                 end if

                 if(iLv.lt.LvMax) then ! if the level of the present oct is coarser than the finest level,
                    if(p0%iFLG(1).gt.0 ) then 
                       p0%octCh1%octPrt => p0%Psort
                       p0%octCh2%octPrt => p0%Psort
                       p0%octCh3%octPrt => p0%Psort
                       p0%octCh4%octPrt => p0%Psort
                       p0%octCh5%octPrt => p0%Psort
                       p0%octCh6%octPrt => p0%Psort
                       p0%octCh7%octPrt => p0%Psort
                       p0%octCh8%octPrt => p0%Psort
                    end if
                 end if

!print *,"set OctCh%octPrt form Psort",rank
!
                 indexP=p0%octP
                 indexQ=p0%Psort%octN
                 pp => p0%ptcl
                 do ip=1,indexP+1
                    pp%Ioct=indexQ
                    pp => pp%prtnxt
                 end do

              end do

!!$omp end parallel do 
              minIF(iLv)=ii0+1
              maxIF(iLv)=ii
           end if

           !sorting for GMesh
           if(iiG.eq.iiG0) then
              GminIF(iLv)=iiG0
              GmaxIF(iLv)=iiG0
           else
!!$omp parallel do private(index,p0,p1,nprt,indexP,indexQ,pp,ip) shared(ii0,ii,iLv,Mesh,Mesh2)
              do index=iiG0+1,iiG

                 p0 => Mesh2(index)
                 p0%octN = index

                 if(associated(p0%octNb1))p0%octNb1 => p0%octNb1%Psort
                 if(associated(p0%octNb2))p0%octNb2 => p0%octNb2%Psort
                 if(associated(p0%octNb3))p0%octNb3 => p0%octNb3%Psort
                 if(associated(p0%octNb4))p0%octNb4 => p0%octNb4%Psort
                 if(associated(p0%octNb5))p0%octNb5 => p0%octNb5%Psort
                 if(associated(p0%octNb6))p0%octNb6 => p0%octNb6%Psort
               
                 if(iLv.gt.-1) then
                    p1 => p0%octPrt
                    nprt = p1%octN
                    p1 => Mesh2(nprt)

                    if(p0%Csort<=4) then 
                       if(p0%Csort==1) then 
                          p1%octCh1 => p0%Psort
                       else if(p0%Csort==2) then 
                          p1%octCh2 => p0%Psort
                       else if(p0%Csort==3) then 
                          p1%octCh3 => p0%Psort
                       else 
                          p1%octCh4 => p0%Psort
                       end if
                    else
                       if(p0%Csort==5) then 
                          p1%octCh5 => p0%Psort
                       else if(p0%Csort==6) then 
                          p1%octCh6 => p0%Psort
                       else if(p0%Csort==7) then 
                          p1%octCh7 => p0%Psort
                       else 
                          p1%octCh8 => p0%Psort
                       end if
                    end if
                 end if
                 
                 if(iLv.lt.LvMax) then
                    if(p0%iFLG(1).gt.0 ) then 
                       p0%octCh1%octPrt => p0%Psort
                       p0%octCh2%octPrt => p0%Psort
                       p0%octCh3%octPrt => p0%Psort
                       p0%octCh4%octPrt => p0%Psort
                       p0%octCh5%octPrt => p0%Psort
                       p0%octCh6%octPrt => p0%Psort
                       p0%octCh7%octPrt => p0%Psort
                       p0%octCh8%octPrt => p0%Psort
                    end if
                 end if

                 indexP=p0%octP
                 indexQ=p0%Psort%octN
                
                 pp => p0%ptcl
                 do ip=1,indexP+1
                    pp%Ioct=indexQ
                    pp => pp%prtnxt
                 end do

              end do

!!$omp end parallel do 
              GminIF(iLv)=iiG0+1
              GmaxIF(iLv)=iiG
           end if

        end do


!!$omp parallel do private(index) shared(ii,Mesh,Mesh2)
        do index=BMeshBound+1,ii
           Mesh(index)=Mesh2(index)
        end do
!!$omp end parallel do
!!$omp parallel do private(index) shared(iiG,Mesh,Mesh2)
        do index=minID(3,0),iiG  !==
           Mesh(index)=Mesh2(index)
        enddo
!!$omp end parallel do 

 
!!$        if(debugMode>=2)then
!!$           print *,"=======================",rank
!!$           print *,"Before MinMaxID MeshRange=",1,"-",MeshBound,rank        
!!$           do iLv=0,LvMax
!!$              print "(A,I2,A,I7,$)","MinID(1,",iLv,")=",MinID(1,iLv)
!!$              print "(A,I2,A,I7,$)","|MaxID(1,",iLv,")=",MaxID(1,iLv)
!!$              print "(A,I2,A,I7,$)","|MinID(2,",iLv,")=",MinID(2,iLv)
!!$              print "(A,I2,A,I7,I4)","|MaxID(2,",iLv,")=",MaxID(2,iLv),rank
!!$           enddo
!!$        endif
!!$
!!$        if(debugMode>=2)then
!!$           print *,"=======================",rank
!!$           print *,"Before GMinMaxID GMeshRange=",MeshBound+1,"-",size(Mesh),rank
!!$           do iLv=-1,LvMax
!!$              print "(A,I2,A,I7,$)","MinID(3,",iLv,")=",MinID(3,iLv)
!!$              print "(A,I2,A,I7,$)","|MaxID(3,",iLv,")=",MaxID(3,iLv)
!!$              print "(A,I2,A,I7,$)","|MinID(4,",iLv,")=",MinID(4,iLv)
!!$              print "(A,I2,A,I7,I4)","|MaxID(4,",iLv,")=",MaxID(4,iLv),rank
!!$           enddo
!!$        endif

        do iLv=0,LvMax
           minID(1,iLv)=minIF(iLv)
           maxID(1,iLv)=maxIF(iLv)
        end do

        do iLv=0,LvMax
           minID(3,iLv)=GminIF(iLv)
           maxID(3,iLv)=GmaxIF(iLv)
        end do

        do iLv=0,LvMax-1
           call copy_particleD(   iLv, 1, 2, 1 )
!           if(CIP.eq.0)then
!            call interpolate_fieldE(iLv, 1, 2)
!            call interpolate_fieldB(iLv, 1, 2)
!            call interpolate_fieldJ(iLv, 1, 2,2)  !???
!           else
            call interpolate_field(iLv, 1, 2,"F")
            call interpolate_fieldC(iLv, 1, 2,2)  !???
!           end if
           minID(2,iLv)=maxIF(LvMax)
           maxID(2,iLv)=maxIF(LvMax)
        end do

        do iLv=0,LvMax-1
           minID(4,iLv)=GmaxIF(LvMax)
           maxID(4,iLv)=GmaxIF(LvMax)
        end do

        do iLv=0,LvMax
           call delete_particleT( iLv, 4, 64)
!           print *,'delete_particleT_2'
           call delete_particleT( iLv,-3,-1)
!           print *,'delete_particleT_3'
        end do

        minID(2,LvMax)=maxIF(LvMax)
        maxID(2,LvMax)=maxIF(LvMax)
        minID(4,LvMax)=GmaxIF(LvMax)
        maxID(4,LvMax)=GmaxIF(LvMax)

!!$        if(debugMode>=2)then
!!$           print *,"||||||||||||||||",rank
!!$           print *,"VVVVVVVVVVVVVVVV",rank
!!$           print *,"After MinMaxID MeshRange=",1,"-",MeshBound,rank
!!$           do iLv=0,LvMax
!!$              print "(A,I2,A,I7,$)","MinID(1,",iLv,")=",MinID(1,iLv)
!!$              print "(A,I2,A,I7,$)","|MaxID(1,",iLv,")=",MaxID(1,iLv)
!!$              print "(A,I2,A,I7,$)","|MinID(2,",iLv,")=",MinID(2,iLv)
!!$              print "(A,I2,A,I7,I4)","|MaxID(2,",iLv,")=",MaxID(2,iLv),rank
!!$           enddo
!!$           print *,"=======================",rank
!!$        endif
!!$
!!$        if(debugMode>=2)then
!!$           print *,"After GMinMaxID GMeshRange=",MeshBound+1,"-",size(Mesh),rank
!!$           do iLv=-1,LvMax
!!$              print "(A,I2,A,I7,$)","MinID(3,",iLv,")=",MinID(3,iLv)
!!$              print "(A,I2,A,I7,$)","|MaxID(3,",iLv,")=",MaxID(3,iLv)
!!$              print "(A,I2,A,I7,$)","|MinID(4,",iLv,")=",MinID(4,iLv)
!!$              print "(A,I2,A,I7,I4)","|MaxID(4,",iLv,")=",MaxID(4,iLv),rank
!!$           enddo
!!$           print *,"=======================",rank
!!$        endif

!!$omp parallel do private(index,iFLG0,p0) shared(ii,Mesh,BMeshBound)
        do index=BMeshBound+1,ii
           p0 => Mesh(index)
           iFLG0 = p0%iFLG(1)
           if((iFLG0.ge.1).and.(iFLG0.le.3)) then 
              p0%iFLG(2)=5
           else
              p0%iFLG(2)=0
           endif
        end do
!!$omp end parallel do
!!$omp parallel do private(index,iFLG0,p0) shared(iiG,Mesh,MeshBound)
        do index=MeshBound+1,iiG
           p0 => Mesh(index)
           iFLG0 = p0%iFLG(1)
           if((iFLG0.ge.1).and.(iFLG0.le.3)) then 
              p0%iFLG(2)=5
           else
              p0%iFLG(2)=0
           endif
        end do
!!$omp end parallel do
!
        
        if(debugMode>=3)then
           print *,"debug call after sort_oct",rank
           call full_oct_check
           print *,"debug call after sort_oct ok",rank
        endif


        !if(debugMode>=1)print*, rank, 'exiting sort_oct'
        return
        end subroutine sort_oct

        subroutine reset_sortBuffer
          use oct_set
          use init_mesh_size
          use param
          implicit none
          integer(kind=4)::index
          integer(kind=4)::size_mesh
          type(oct),pointer::p0
          
          size_mesh=size(Mesh2)
!!$omp parallel do private(index,p0) shared(size_mesh,Mesh)
          do index=1,size_mesh
             p0=>Mesh2(index)

             p0%iFLG=(/-4,0,0/)
             p0%octP=0
             nullify(p0%ptcl)
             nullify(p0%octNb1)
             nullify(p0%octNb2)
             nullify(p0%octNb3)
             nullify(p0%octNb4)
             nullify(p0%octNb5)
             nullify(p0%octNb6)
        
             nullify(p0%octCh1)
             nullify(p0%octCh2)
             nullify(p0%octCh3)
             nullify(p0%octCh4)
             nullify(p0%octCh5)
             nullify(p0%octCh6)
             nullify(p0%octCh7)
             nullify(p0%octCh8)
                       
          enddo
!!$omp end parallel do

        end subroutine reset_sortBuffer

! ===================================================================-
!
!2011/09/21 yagi added
!=========================================
  subroutine prepare_Gsort
    use oct_set
    use init_mesh_size
    use param
    implicit none
    integer(kind=4)::index,iLv,i
    type(oct),pointer::p0

    do iLv=-1,LvMax
       do i=3,4
          do index=MinID(i,iLv),MaxID(i,iLv)
             p0=>Mesh(index)
             p0%Psort=>p0
          enddo
       enddo
    enddo

  end subroutine prepare_Gsort
!***********************************************************************
      subroutine set_oct
! +-------------------------------------------------------------------+
! |                                                                   |
! | Set & create initial self oct with initial values                 |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        use param
        use message_passing_interface
! ---------------------------------------------------------------
        implicit none
        integer(kind=4)    :: octN,octLv,octP
        integer(kind=4)    :: iFLG(3),Csort
        integer(kind=4)    :: iPOS(3)
        real(kind=8)       :: rPOS(3)
        integer(kind=4)    :: MrtN
        integer(kind=4)  :: prc_bndry, octType, ptcl_loops
        integer(kind=4)    :: j,n1,n2,ix,iy,iz,index,iii
        real(kind=8)       :: ixx,iyy,izz
        integer(kind=4)    :: Nint
        integer(kind=4)    :: AllocateStatus
        integer(kind=4)    :: iC(iomp0)
        real(kind=8)       :: F(18),C(6*iomp0),Z(IonSorts),G(1:18),D(1:3),O(1:16)
        type(oct), pointer :: newp
! ---------------------------------------------------------------
!
        Nall = Nall_ini !Tatsuki Added 0401/2011
!
!        Nint = 2**LvMax
!        print *, 'Nint=',Nint
!
! -- allocate & nullify to a new pointer 'newp' --
        allocate(newp, stat=AllocateStatus)
        if(AllocateStatus /=0) then
          print *, '!!! Not enough memory !!!'
          stop
        endif
        nullify(newp%octPrt)
        nullify(newp%octNb1) ; nullify(newp%octNb2)
        nullify(newp%octNb3) ; nullify(newp%octNb4)
        nullify(newp%octNb5) ; nullify(newp%octNb6)
        nullify(newp%octCh1) ; nullify(newp%octCh2)
        nullify(newp%octCh3) ; nullify(newp%octCh4)
        nullify(newp%octCh5) ; nullify(newp%octCh6)
        nullify(newp%octCh7) ; nullify(newp%octCh8)
        nullify(newp%Psort)
        nullify(newp%ptcl)   ; nullify(newp%ptclA) 
!
! -- Set initial variables of the oct ---
!
        Nint = 2**(LvMax+1)
        octP  = 0
        octLv = -1
        iFLG  = 0
        Csort = 0

        F = 0.d0
        C = 0.d0
        Z = 0.d0
        G = 0.d0
        D=0.d0
        O=0.d0
        iC= 0

!------------------------------------------------
!       Level -1 Operation
        do j = 1, Nall_ini/8
          octN  = j
          index = j
          MrtN  = j*8
!
! -- Extract position index --
          n1 = j
          iz = int((n1-1)/(NXB*NYB/4)) + 1
          n2 = n1 - NXB*NYB/4*(iz-1)
          iy = int((n2-1)/(NXB/2)) + 1
          ix = n2 - NXB/2*(iy-1)
          prc_bndry = 0
          octType = 0
          ptcl_loops=0
! -- Position Index --
          ix = 2*ix*Nint - Nint
          iy = 2*iy*Nint - Nint
          iz = 2*iz*Nint - Nint
!!$          iPOS(1) = ix
!!$          iPOS(2) = iy
!!$          iPOS(3) = iz
!
! -- Position --
          ixx = real(ix + Nint)/real(2*Nint)
          iyy = real(iy + Nint)/real(2*Nint)
          izz = real(iz + Nint)/real(2*Nint)


          
! How about this!
          call i2pos(j, rPOS, iPOS, -1, 0)
if( ipos(1) /= ix .or. ipos(2) /= iy .or. ipos(3) /= iz) then 
   print*, 'there is a difference b/w ix, iy, iz and ipos, in Lv -1', ix, iy, iz, ipos, index
stop
endif
!
! -- Mesh Array --
          Mesh(index)=                               &
                 oct(octN,octLv,octP,Csort,iFLG,     &
                 iPOS,rPOS,                          &
                 MrtN,iC,                     &
                 prc_bndry,                    &
                 octType,   &
                 F,C,Z,G,D,O,ptcl_loops,       & 
                 newp%octPrt,                        &
                 newp%octNb1, newp%octNb2,           &
                 newp%octNb3, newp%octNb4,           &
                 newp%octNb5, newp%octNb6,           &
                 newp%octCh1, newp%octCh2,           &
                 newp%octCh3, newp%octCh4,           &
                 newp%octCh5, newp%octCh6,           &
                 newp%octCh7, newp%octCh8,           &
                 newp%Psort ,                        &
                 newp%ptcl  , newp%ptclA)
        enddo
! ------------------------------------------------------------------------
        Nint = 2**LvMax ! for Lv 0
        octLv = 0
!
! -- Loop for ganeration of Level 0 ---
!
!!$open(88, file="position0.dat"  , form=  'formatted')
!!$open(89, file="position1.dat"  , form=  'formatted')
!!$open(87, file="position2.dat"  , form=  'formatted')
!!$open(86, file="position3.dat"  , form=  'formatted')
        do j = 1, Nall_ini
!        do j = 1, NXB
          octN  = j+Nall_ini/8
          index = j+Nall_ini/8
          MrtN  = j+Nall_ini/8
          prc_bndry=0
!
! -- Extract position index --
          n1 = j
          iz = int((n1-1)/(NXB*NYB)) + 1
          n2 = n1 - NXB*NYB*(iz-1)
          iy = int((n2-1)/NXB) + 1
          ix = n2 - NXB*(iy-1)
!
! -- Position Index --
          ix = 2*ix*Nint - Nint
          iy = 2*iy*Nint - Nint
          iz = 2*iz*Nint - Nint
!!$          iPOS(1) = ix
!!$          iPOS(2) = iy
!!$          iPOS(3) = iz
!
! -- Position --
          ixx = real(ix + Nint)/real(2*Nint)
          iyy = real(iy + Nint)/real(2*Nint)
          izz = real(iz + Nint)/real(2*Nint)

! Obrain the position from index using i2pos 
call i2pos(index, rPOS, iPOS, 0, 0)

if( ipos(1) /= ix .or. ipos(2) /= iy .or. ipos(3) /= iz) then 
   print*, 'there is a difference b/w ix, iy, iz and ipos Lv 0 '
stop
endif

!if(rank==1)         print*,  ix, iy, iz, iPOS
!print*, rpos(1)
!!$if (rank==1) write(89,*),   index,  rpos(1), ixx*dx(1) - HALF*dx(1) + real(rank * NXB-(buf_width+1)*(rank+1)) * dx(1) 
!!$if (rank==0) write(88,*),   index,  rpos(1), ixx*dx(1) - HALF*dx(1) + real(rank * NXB-(buf_width+1)*(rank+1)) * dx(1) 
!!$if (rank==2) write(87,*),   index,  rpos(1), ixx*dx(1) - HALF*dx(1) + real(rank * NXB-(buf_width+1)*(rank+1)) * dx(1) 
!!$if (rank==3) write(86,*),   index,  rpos(1), ixx*dx(1) - HALF*dx(1) + real(rank * NXB-(buf_width+1)*(rank+1)) * dx(1) 

!
! -- Mesh Array --
          Mesh(index) =                              &
                 oct(octN,octLv,octP,Csort,iFLG,     &
                 iPOS,rPOS,                          &
                 MrtN,iC,                      &
                prc_bndry,                      &
                octType,                        &
                 F,C,Z,G,D,O, ptcl_loops,      & 
                 newp%octPrt,                        &
                 newp%octNb1, newp%octNb2,           &
                 newp%octNb3, newp%octNb4,           &
                 newp%octNb5, newp%octNb6,           &
                 newp%octCh1, newp%octCh2,           &
                 newp%octCh3, newp%octCh4,           &
                 newp%octCh5, newp%octCh6,           &
                 newp%octCh7, newp%octCh8,           &
                 newp%Psort ,                        &
                 newp%ptcl  , newp%ptclA)
        enddo

!!$close(88)
!!$close(89)
!!$close(87)
!!$close(86)
!!$stop
!--------------------------------------------------------------------
!----------------------------
! Level -1
        MinID(1,-1) = 1
        MaxID(1,-1) = Nall_ini/8
        MinID(2,-1) = Nall_ini/8
        MaxID(2,-1) = Nall_ini/8
!----------------------------
! Level 0

        MinID(1,0)  = 1+Nall_ini/8
        MaxID(1,0)  = Nall_ini+Nall_ini/8
        MinID(2,0)  = Nall_ini+Nall_ini/8
        MaxID(2,0)  = Nall_ini+Nall_ini/8

!----------------------------
! Level > 1   We do not create Mesh for Lv > 1 in this routine

        if(LvMax>0) then 
           do j=1,LvMax
              MinID(1,j) = Nall_ini+Nall_ini/8
              MaxID(1,j) = Nall_ini+Nall_ini/8
              MinID(2,j) = Nall_ini+Nall_ini/8
              MaxID(2,j) = Nall_ini+Nall_ini/8
           enddo
        end if
!
!        print *, ' MinID(1,0), MaxID(1,0) = ',MinID(1,0),MaxID(1,0)
!       
        deallocate(newp,stat=iii)

        return
!               
      end subroutine set_oct
!
!***********************************************************************
      subroutine periodic_boundary
! +-------------------------------------------------------------------+
! |                                                                   |
! | Create initial pointer connection for periodic boundary condition |
! |                                                                   |
! |  Order of neighbor index for octNb#1-#6                           |
! |                                                                   |
! |   Y                            Z                                  |
! |   |                            |                                  |
! |   |          octNb4            |          octNb6                  |
! |   +---X        |               +---X        |                     |
! |                |                            |                     |
! |    octNb1 -- octN -- octNb2     octNb1 -- octN -- octNb2          |
! |                |                            |                     |
! |                |                            |                     |
! |              octNb3                       octNb5                  |
! |                                                                   |
! |                                                                   |
! |  (ex.) for x-direction                                            |
! |                                                                   |
! |        ---+------+------+---...---+------+------+------+---       |
! |     ix:   1      2      3   ... NXB-1   NXB   NXB+1  NXB+2        |
! |           |      |      |         |      |      |      |          |
! |         NXB+1  NXB+2  NXB+3 ...  -1      0      1      2          |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use init_mesh_size
        use const
! ----------------------------------------------------
        implicit none
        type(oct), pointer :: p0, p1
!
        integer(kind=4) :: n1,n2,ii
        integer(kind=4) :: ix,iy,iz,index,index0
        integer(kind=4) :: jx,jy,jz,jxx,jyy,jzz,jnm
! ----------------------------------------------------
!
        do index=1,Nall_ini/8
          p0 => Mesh(index)
!
! -- extract position index --
          n1 = index
          iz = int((n1-1)/(NXB*NYB/4)) + 1
          n2 = n1 - NXB*NYB/4*(iz-1)
          iy = int((n2-1)/(NXB/2)) + 1
          ix = n2 - NXB/2*(iy-1)
!
          if( ix==1 .or. ix==NXB/2 .or.       &
              iy==1 .or. iy==NYB/2 .or.       &
              iz==1 .or. iz==NZB/2 ) then
!
! -- periodic boundary condition--
!  --- for -x neighbor ---
            if(ix==1) then
              p0%octNb1 => Mesh(index+(NXB/2-1))
            endif
!  --- for +x neighbor ---
            if(ix==NXB/2) then
              p0%octNb2 => Mesh(index-(NXB/2-1))
            endif
!  --- for -y neighbor ---
            if(iy==1) then
              p0%octNb3 => Mesh(index+(NXB/2)*(NYB/2-1))
            endif
!  --- for +y neighbor ---
            if(iy==NYB/2) then
              p0%octNb4 => Mesh(index-(NXB/2)*(NYB/2-1))
            endif
!  --- for -z neighbor ---
            if(iz==1) then
              p0%octNb5 => Mesh(index+(NXB*NYB/4)*(NZB/2-1))
            endif
!  --- for +z neighbor ---
            if(iz==NZB/2) then
              p0%octNb6 => Mesh(index-(NXB*NYB/4)*(NZB/2-1))
            endif
!
          endif
        enddo
!
!================================================================
!
        do index=1,Nall_ini
           index0 = index+Nall_ini/8
          p0 => Mesh(index0)
!
! -- extract position index --
          n1 = index
          iz = int((n1-1)/(NXB*NYB)) + 1
          n2 = n1 - NXB*NYB*(iz-1)
          iy = int((n2-1)/NXB) + 1
          ix = n2 - NXB*(iy-1) 
!
          if( ix==1 .or. ix==NXB .or.       &
              iy==1 .or. iy==NYB .or.       &
              iz==1 .or. iz==NZB ) then
!
! -- periodic boundary condition--
!  --- for -x neighbor ---
            if(ix==1) then
              p0%octNb1 => Mesh(index0+(NXB-1))
            endif
!  --- for +x neighbor ---
            if(ix==NXB) then
              p0%octNb2 => Mesh(index0-(NXB-1))
            endif
!  --- for -y neighbor ---
            if(iy==1) then
              p0%octNb3 => Mesh(index0+NXB*(NYB-1))
            endif
!  --- for +y neighbor ---
            if(iy==NYB) then
              p0%octNb4 => Mesh(index0-NXB*(NYB-1))
            endif
!  --- for -z neighbor ---
            if(iz==1) then
              p0%octNb5 => Mesh(index0+NXB*NYB*(NZB-1))
            endif
!  --- for +z neighbor ---
            if(iz==NZB) then
              p0%octNb6 => Mesh(index0-NXB*NYB*(NZB-1))
            endif
!
          endif
        enddo
!
!----------------------------------------------------------------
!
       ii=2**(Lvmax+1)

       do index=minID(1,-1),maxID(1,-1)
          p0 => Mesh(index)
          jx = p0%iPOS(1)
          jy = p0%iPOS(2)
          jz = p0%iPOS(3)
!-- ch1 --
          jxx= jx-1 ; jyy= jy-1 ; jzz= jz-1
          jnm= jxx/ii+NXB*(jyy/ii)+(NXB*NYB)*(jzz/ii)+1
          p1 => Mesh(jnm+Nall_ini/8)
          p0%octCh1 => p1
          p1%octPrt => p0
          p1%Csort = 1
!-- ch2 --
          jxx= jx+1 ; jyy= jy-1 ; jzz= jz-1
          jnm= jxx/ii+NXB*(jyy/ii)+(NXB*NYB)*(jzz/ii)+1
          p1 => Mesh(jnm+Nall_ini/8)
          p0%octCh2 => p1
          p1%octPrt => p0
          p1%Csort = 2
!-- ch3 --
          jxx= jx-1 ; jyy= jy+1 ; jzz= jz-1
          jnm= jxx/ii+NXB*(jyy/ii)+(NXB*NYB)*(jzz/ii)+1
          p1 => Mesh(jnm+Nall_ini/8)
          p0%octCh3 => p1
          p1%octPrt => p0
          p1%Csort = 3
!-- ch4 --
          jxx= jx+1 ; jyy= jy+1 ; jzz= jz-1
          jnm= jxx/ii+NXB*(jyy/ii)+(NXB*NYB)*(jzz/ii)+1
          p1 => Mesh(jnm+Nall_ini/8)
          p0%octCh4 => p1
          p1%octPrt => p0
          p1%Csort = 4
!-- ch5 --
          jxx= jx-1 ; jyy= jy-1 ; jzz= jz+1
          jnm= jxx/ii+NXB*(jyy/ii)+(NXB*NYB)*(jzz/ii)+1
          p1 => Mesh(jnm+Nall_ini/8)
          p0%octCh5 => p1
          p1%octPrt => p0
          p1%Csort = 5
!-- ch6 --
          jxx= jx+1 ; jyy= jy-1 ; jzz= jz+1
          jnm= jxx/ii+NXB*(jyy/ii)+(NXB*NYB)*(jzz/ii)+1
          p1 => Mesh(jnm+Nall_ini/8)
          p0%octCh6 => p1
          p1%octPrt => p0
          p1%Csort = 6
!-- ch7 --
          jxx= jx-1 ; jyy= jy+1 ; jzz= jz+1
          jnm= jxx/ii+NXB*(jyy/ii)+(NXB*NYB)*(jzz/ii)+1
          p1 => Mesh(jnm+Nall_ini/8)
          p0%octCh7 => p1
          p1%octPrt => p0
          p1%Csort = 7
!-- ch8 --
          jxx= jx+1 ; jyy= jy+1 ; jzz= jz+1
          jnm= jxx/ii+NXB*(jyy/ii)+(NXB*NYB)*(jzz/ii)+1
          p1 => Mesh(jnm+Nall_ini/8)
          p0%octCh8 => p1
          p1%octPrt => p0
          p1%Csort = 8
       end do

!  555   format(I5,X,6(F15.7,X))
!
        return
!
      end subroutine periodic_boundary
!***********************************************************************
      subroutine set_flag_block(iLv)
! +-------------------------------------------------------------------+
! |                  estimation of flagging condition                 |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
		use Vparam
! ------------------------------------------
        implicit none
        integer(kind=4) :: iLv
        integer(kind=4) :: index,maxindex,minindex,Gmaxindex,Gminindex,i,ii
        integer(kind=4) :: x,y,z,cstep,value0,iyz
        real(kind=8)    :: Fp1,Fp2,Fp3!,Fp0,Fp4,Fp5,Fp6,Fp7,Fp8
        real(kind=8)    :: Fq1,Fq2,Fq3!,Fq0,Fq4,Fq5,Fq6,Fq7,Fq8
    !    real(kind=8)    :: Fr0,Fr1,Fr2,Fr3,Fr4,Fr5,Fr6,Fr7,Fr8
        real(kind=8)    :: crtr(LvMax)
     !   real(kind=8)    :: xc,yc,zc
!-------------------------------------------
        type(oct), pointer :: p0,p1,p2,p3,p4,p5
!-------------------------------------------
        !if(debugMode>=1)print*,'Entering set_flag_block iLv=',iLv,rank

        nullify(p0)
        nullify(p1)
        nullify(p2)
        nullify(p3)

        maxindex = MaxID(1,iLv)
        minindex = MinID(1,iLv)
        Gmaxindex= MaxID(3,iLv)
        Gminindex= MinID(3,iLv)

        if(iLv==0) then
           !set iFLG(3)
           do index=minindex,maxindex
              p0=>Mesh(index)
              if(p0%iFLG(3)==1)then
                 if(p0%octCh1%iFLG(1)>=1)cycle
                 if(p0%octCh2%iFLG(1)>=1)cycle
                 if(p0%octCh3%iFLG(1)>=1)cycle
                 if(p0%octCh4%iFLG(1)>=1)cycle
                 if(p0%octCh5%iFLG(1)>=1)cycle
                 if(p0%octCh6%iFLG(1)>=1)cycle
                 if(p0%octCh7%iFLG(1)>=1)cycle
                 if(p0%octCh8%iFLG(1)>=1)cycle
              endif

              if(AMRcondition==1)then
                    Fp1 = dble(p0%Z(1))
                    if(Fp1>crtr1(iLv))then
                        p0%iFLG(3)=1
                    else
                        p0%iFLG(3)=0
                    endif
              else if(AMRcondition==2)then
                    Fp1 = dsqrt(p0%F(4)**2+p0%F(5)**2+p0%F(6)**2)
                    if(Fp1>crtr1(iLv))then
                        p0%iFLG(3)=1
                    else
                        p0%iFLG(3)=0
                    endif
              else if(AMRcondition==3)then
                    Fp1 = p0%rPos(1)
                    Fp2 = p0%rPos(2)
                    Fp3 = p0%rPos(3)
!                    if((Fp1.le.crtr0(iLv)).or.(Fp1.ge.crtr1(iLv)))then
!                        p0%iFLG(3)=0
!                    else 
!                        p0%iFLG(3)=1
!                    endif
                    p0%iFLG(3)=0
                    cstep=istep+int(0.25d0*DataIntvl2)
                    value0=mod(cstep,DataIntvl2)
!                    if((value0>int(0.5d0*DataIntvl2)).and.(value0<=DataIntvl2))then
!                        crtr(iLv)=crtr0(iLv)-14.d0*dx(1)
!                    else
                        crtr(iLv)=crtr0(iLv)
!                    end if
                    if((Fp1>=crtr(iLv)).and.(Fp1<=crtr1(iLv)))then
!                    if((Fp2>=crtr0(iLv)).and.(Fp2<=crtr1(iLv)))then
!                    if((Fp3>=crtr0(iLv)).and.(Fp3<=crtr1(iLv)))then
                        p0%iFLG(3)=1
!                    end if
!                    end if
                    end if
                    if(Model==7)then
                    if(Fp1>=Dpos(1)+3*dx(1))then
                        p0%iFLG(3)=0
                    end if
                    end if
              endif
           enddo

		   
!==================Nakano2016/02/13===================================
if(iLv == 0 .and. beamParticle == 1)then
	do iyz=1,arraySize
        p4 => VMesh(iyz)
		do index=minindex,maxindex
			p5=>Mesh(index)
			if(p4%rPos(1) == p5%rPos(1) .and. p4%rPos(2) == p5%rPos(2) .and. p4%rPos(3) == p5%rPos(3) ) then
				p5%octNb2%iFLG(3) = 0
				!call set_inject_flag(p5%octNb2%octN,3)
				call set_inject_flag2(p5%octNb2%octN)
				exit
			end if 
		end do
	end do
end if
!========================================================



           !set iFLG(1)
           do index=minindex,maxindex
              p0=>Mesh(index)
              p0%iFLG(1)=min(p0%iFLG(1),4)
              p0%iFLG(2)=p0%iFLG(1)
              p0%iFLG(1)=0
           enddo

           do index=Gminindex,Gmaxindex
              p0=>Mesh(index)
              p0%iFLG(1)=min(p0%iFLG(1),4)
              p0%iFLG(2)=p0%iFLG(1)
              p0%iFLG(1)=0
           enddo

           !check 7*7*7 oct
           do index=minindex,maxindex
              p0=>Mesh(index)
              if(p0%iFLG(3)/=0)then
                 p1=>p0%octNb1%octNb1%octNb1%octNb3%octNb3%octNb3%octNb5%octNb5%octNb5
                 p2=>p1
                 p3=>p2
                 do z=1,7
                    p2=>p1
                    do y=1,7
                       p3=>p2
                       do x=1,7
                          p3%iFLG(1)=4
                          p3=>p3%octNb2
                       enddo
                       p2=>p2%octNb4
                    enddo
                    p1=>p1%octNb6
                 enddo
              endif
           enddo

!!$omp parallel do private(index,Fp0,Fp1,Fp2,Fp3,Fp4,Fp5,Fp6,Fp7,Fp8) 

!--- If you want to enable AMR, turn this section on.
!--- crtr1(0) = 5 is advised. 
!--- Refinement criterion is the particle number per cell.

!--- AMR section End

        else if(iLv==LvMax) then 
           if(minindex.lt.maxindex) then 
              do index=minindex,maxindex
                 p0 => Mesh(index)
                 p0%iFLG(1)=p0%octPrt%iFLG(1)-4
              end do
           end if
           
           !for Goct
           if(MinID(3,iLv)<MaxID(3,iLv))then
              do index=MinID(3,iLv),MaxID(3,iLv)
                 p0=> Mesh(index)
                 p0%iFLG(1)=p0%octPrt%iFLG(1)-4
              enddo
           endif
           !if(debugMode>=1)print*,'Exiting set_flag_block iLv=',iLv,rank
           return         
        else
           if(minindex .lt. maxindex ) then
              do index=minindex,maxindex
                 p0=>Mesh(index)
                 if(p0%octPrt%iFLG(3)==1)then
                    if(p0%iFLG(3)==1)then
                       if(p0%octCh1%iFLG(1)>=1)cycle
                       if(p0%octCh2%iFLG(1)>=1)cycle
                       if(p0%octCh3%iFLG(1)>=1)cycle
                       if(p0%octCh4%iFLG(1)>=1)cycle
                       if(p0%octCh5%iFLG(1)>=1)cycle
                       if(p0%octCh6%iFLG(1)>=1)cycle
                       if(p0%octCh7%iFLG(1)>=1)cycle
                       if(p0%octCh8%iFLG(1)>=1)cycle
                    endif
                    
                    if(AMRcondition==1)then
                        Fq1 = dble(p0%Z(1))
                        if(Fq1>crtr1(iLv))then
                            p0%iFLG(3)=1
                        else
                            p0%iFLG(3)=0
                        endif
                    else if(AMRcondition==2)then
                        Fq1 = dsqrt(p0%F(4)**2+p0%F(5)**2+p0%F(6)**2)
                        if(Fq1>crtr1(iLv))then
                            p0%iFLG(3)=1
                        else
                            p0%iFLG(3)=0
                        endif
                    else if(AMRcondition==3)then
                        Fq1 = p0%rPos(1)
                        Fq2 = p0%rPos(2)
                        Fq3 = p0%rPos(3)
!                        if((Fq1.le.crtr0(iLv)).or.(Fq1.ge.crtr1(iLv)))then
!                            p0%iFLG(3)=0
!                        else 
!                            p0%iFLG(3)=1
!                        end if
                        p0%iFLG(3)=0
                        cstep=istep+int(0.25d0*DataIntvl2)
                        value0=mod(cstep,DataIntvl2)
!                        if((value0>int(0.5d0*DataIntvl2)).and.(value0<=DataIntvl2))then
!                            crtr(iLv)=crtr0(iLv)-14.d0*dx(1)
!                        else
                            crtr(iLv)=crtr0(iLv)
!                        end if
                        if((Fq1>=crtr(iLv)).and.(Fq1<=crtr1(iLv)))then
!                        if((Fq2>=crtr0(iLv)).and.(Fq2<=crtr1(iLv)))then
!                        if((Fq3>=crtr0(iLv)).and.(Fq3<=crtr1(iLv)))then
                            p0%iFLG(3)=1
!                        end if
!                        end if
                        end if
                        if(Model==7)then
                        if(Fq1>=Dpos(1)+3*dx(1))then
                            p0%iFLG(3)=0
                        end if
                        end if
                    endif
                 endif
              enddo

              do index=minindex,maxindex
                 p0=>Mesh(index)
                 p0%iFLG(1)=min(p0%iFLG(1),4)
                 p0%iFLG(2)=p0%iFLG(1)
                 p0%iFLG(1)=p0%octPrt%iFLG(1)-4
              enddo
              
              do index=Gminindex,Gmaxindex
                 p0=>Mesh(index)
                 p0%iFLG(1)=min(p0%iFLG(1),4)
                 p0%iFLG(2)=p0%iFLG(1)
                 p0%iFLG(1)=p0%octPrt%iFLG(1)-4
              enddo

              !check 7*7*7 oct
              do index=minindex,maxindex
                 p0=>Mesh(index)
                 if(p0%iFLG(3)/=0)then
                    p1=>p0%octNb1%octNb1%octNb1%octNb3%octNb3%octNb3%octNb5%octNb5%octNb5
                    p2=>p1
                    p3=>p2
                    do z=1,7
                       p2=>p1
                       do y=1,7
                          p3=>p2
                          do x=1,7
                             p3%iFLG(1)=4
                             p3=>p3%octNb2
                          enddo
                          p2=>p2%octNb4
                       enddo
                       p1=>p1%octNb6
                    enddo
                 endif
              enddo

           endif
        end if
!!$omp end parallel do
        
        call refresh_iFLG(iLv)
        do index=MinID(3,iLv),MaxID(3,iLv) !set iFLG(1)=4 from Goct
           p0=>Mesh(index)
           if(p0%iFLG(3)/=0)then
              !1: -x -y -z
              p1=>p0
              do z=1,4
                 p2=>p1
                 do y=1,4
                    p3=>p2
                    do x=1,4
                       p3%iFLG(1)=4
                       p3=>p3%octNb1
                    enddo
                    p2=>p2%octNb3
                 enddo
                 p1=>p1%octNb5
              enddo
              !2: +x -y -z
              p1=>p0
              do z=1,4
                 p2=>p1
                 do y=1,4
                    p3=>p2
                    do x=1,4
                       p3%iFLG(1)=4
                       p3=>p3%octNb2
                    enddo
                    p2=>p2%octNb3
                 enddo
                 p1=>p1%octNb5
              enddo
              !3: -x +y -z
              p1=>p0
              do z=1,4
                 p2=>p1
                 do y=1,4
                    p3=>p2
                    do x=1,4
                       p3%iFLG(1)=4
                       p3=>p3%octNb1
                    enddo
                    p2=>p2%octNb4
                 enddo
                 p1=>p1%octNb5
              enddo
              !4: +x +y -z
              p1=>p0
              do z=1,4
                 p2=>p1
                 do y=1,4
                    p3=>p2
                    do x=1,4
                       p3%iFLG(1)=4
                       p3=>p3%octNb2
                    enddo
                    p2=>p2%octNb4
                 enddo
                 p1=>p1%octNb5
              enddo
              !5: -x -y +z
              p1=>p0
              do z=1,4
                 p2=>p1
                 do y=1,4
                    p3=>p2
                    do x=1,4
                       p3%iFLG(1)=4
                       p3=>p3%octNb1
                    enddo
                    p2=>p2%octNb3
                 enddo
                 p1=>p1%octNb6
              enddo
              !6: +x -y +z
              p1=>p0
              do z=1,4
                 p2=>p1
                 do y=1,4
                    p3=>p2
                    do x=1,4
                       p3%iFLG(1)=4
                       p3=>p3%octNb2
                    enddo
                    p2=>p2%octNb3
                 enddo
                 p1=>p1%octNb6
              enddo
              !7: -x +y +z
              p1=>p0
              do z=1,4
                 p2=>p1
                 do y=1,4
                    p3=>p2
                    do x=1,4
                       p3%iFLG(1)=4
                       p3=>p3%octNb1
                    enddo
                    p2=>p2%octNb4
                 enddo
                 p1=>p1%octNb6
              enddo
              !8: +x +y +z
              p1=>p0
              do z=1,4
                 p2=>p1
                 do y=1,4
                    p3=>p2
                    do x=1,4
                       p3%iFLG(1)=4
                       p3=>p3%octNb2
                    enddo
                    p2=>p2%octNb4
                 enddo
                 p1=>p1%octNb6
              enddo
           endif
        enddo
        call refresh_iFLG(iLv)
        call set_buffer_out(iLv,4,3)
        call set_buffer_out(iLv,3,2)
        call set_buffer_out(iLv,2,1)
        call refresh_iFLG(iLv)

         do i=1,2
            if(minID(i,iLv).lt.maxID(i,iLv)) then 
!!$omp parallel do private(index,ii)  
               do index=minID(i,iLv),maxID(i,iLv)
                  p0 => Mesh(index)
                  ii=p0%iFLG(2)
                  p0%iFLG(2)=0
                  if(ii>=4) then 
                     if((p0%iFLG(1)<=3).and.(p0%iFLG(1)>=0)) p0%iFLG(2)=-1
                  else if(ii<=0) then
                     if(p0%iFLG(1)>=1) p0%iFLG(2)=2
                  else if(ii<=3) then 
                     if(p0%iFLG(1)>=1) p0%iFLG(2)=1
                  endif
               end do
!!$omp end parallel do
            end if

            if(MinID(i+2,iLv)<MaxID(i+2,iLv))then
               do index=minID(i+2,iLv),maxID(i+2,iLv)
                  p0 => Mesh(index)
                  ii=p0%iFLG(2)
                  p0%iFLG(2)=0
                  if(ii>=4) then 
                     if((p0%iFLG(1)<=3).and.(p0%iFLG(1)>=0)) p0%iFLG(2)=-1
                  else if(ii<=0) then
                     if(p0%iFLG(1)>=1) p0%iFLG(2)=2
                  else if(ii<=3) then 
                     if(p0%iFLG(1)>=1) p0%iFLG(2)=1
                  endif
               end do
            endif
         end do


        !if(debugMode>=1)print*,'Exiting set_flag_block iLv=',iLv,rank
    return
end subroutine set_flag_block
!
recursive subroutine set_inject_flag(octN,index)
		use oct_set
        use param
        use const
        use init_mesh_size
        implicit none
		integer(kind=4),intent(in) :: index,octN
		type(oct), pointer :: p0
		if(index == 0) return 
		p0 => Mesh(octN)
		p0%octNb1%iFLG(3) = 0
		if(index > 0) call set_inject_flag(p0%octNb1%octN,index-1)
		p0%octNb2%iFLG(3) = 0
		if(index > 0) call set_inject_flag(p0%octNb2%octN,index-1)
		p0%octNb3%iFLG(3) = 0
		if(index > 0) call set_inject_flag(p0%octNb3%octN,index-1)
		p0%octNb4%iFLG(3) = 0
		if(index > 0) call set_inject_flag(p0%octNb4%octN,index-1)
		p0%octNb5%iFLG(3) = 0
		if(index > 0) call set_inject_flag(p0%octNb5%octN,index-1)
		p0%octNb6%iFLG(3) = 0
		if(index > 0) call set_inject_flag(p0%octNb6%octN,index-1)

end subroutine set_inject_flag


subroutine set_inject_flag2(octN)
		use oct_set
        use param
        use const
        use init_mesh_size
        implicit none
		integer(kind=4),intent(in) :: octN
		integer(kind=4) :: x,y,z
		type(oct), pointer :: p0,p1,p2,p3
		p0 => Mesh(octN)
		p1=>p0%octNb1%octNb1%octNb1%octNb3%octNb3%octNb3%octNb5%octNb5%octNb5
                    p2=>p1
                    p3=>p2
                    do z=1,7
                       p2=>p1
                       do y=1,7
                          p3=>p2
                          do x=1,7
                             p3%iFLG(3)=0
                             p3=>p3%octNb2
                          enddo
                          p2=>p2%octNb4
                       enddo
                       p1=>p1%octNb6
                    enddo
end subroutine set_inject_flag2





!***********************************************************************
      subroutine set_flag_prt(iLv)
! +-------------------------------------------------------------------+
! |                  estimation of flagging condition                 |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
! ------------------------------------------
        implicit none
        integer(kind=4) :: iLv
        integer(kind=4) :: index,maxindex,minindex,iFLG0,iFLG1
!-------------------------------------------
        type(oct), pointer :: p0
!-------------------------------------------

        maxindex = MaxID(1,iLv)
        minindex = MinID(1,iLv)
        if(maxindex.eq.minindex) return 

        if(iLv.eq.0) then 
!!$omp parallel do private(index,iFLG0,iFLG1,p0) 
           do index=minindex,maxindex
              p0 => Mesh(index)
              if(p0%octPrt%iFLG(1).lt.4) p0%octPrt%iFLG(1)=4
              if(p0%iFLG(1).ge.0) then 
                 !
                 iFLG0 = 4+p0%iFLG(1)
                 iFLG1 = p0%octPrt%iFLG(1)

                 p0%octPrt%iFLG(1)=iFLG0
                 if(iFLG0.lt.iFLG1)   p0%octPrt%iFLG(1)=iFLG0
              endif
           end do
!!$omp end parallel do
!
           return
        end if

!!$omp parallel do private(index,iFLG0,iFLG1,p0) 
           do index=minindex,maxindex
              p0 => Mesh(index)
!
              if(p0%iFLG(1).gt.0) then 
                 iFLG0 = 4+p0%iFLG(1)
                 iFLG1 = p0%octPrt%iFLG(1)
                 if(iFLG1.le.6) then 
                    p0%octPrt%iFLG(1)=iFLG0
                 else if(iFLG0.lt.iFLG1) then 
                    p0%octPrt%iFLG(1)=iFLG0
                 end if

              endif
           end do
!!$omp end parallel do
!
        return
        end subroutine set_flag_prt
!
!***********************************************************************
      subroutine set_flag_prtG(iLv)
! +-------------------------------------------------------------------+
! |                  estimation of flagging condition                 |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
! ------------------------------------------
        implicit none
        integer(kind=4) :: iLv
        integer(kind=4) :: index,maxindex,minindex,iFLG0,iFLG1
!-------------------------------------------
        type(oct), pointer :: p0
!-------------------------------------------

        maxindex = MaxID(3,iLv)
        minindex = MinID(3,iLv)
        if(maxindex.eq.minindex) return 

        if(iLv.eq.0) then 
!!$omp parallel do private(index,iFLG0,iFLG1,p0) 
           do index=minindex,maxindex
              p0 => Mesh(index)

              if(p0%octPrt%iFLG(1).lt.4) p0%octPrt%iFLG(1)=4
              if(p0%iFLG(1).ge.0) then 
                 !
                 iFLG0 = 4+p0%iFLG(1)
                 iFLG1 = p0%octPrt%iFLG(1)

                 p0%octPrt%iFLG(1)=iFLG0
                 if(iFLG0.lt.iFLG1)   p0%octPrt%iFLG(1)=iFLG0
              end if
           end do
!!$omp end parallel do
!
           return
        end if

!!$omp parallel do private(index,iFLG0,iFLG1,p0) 
           do index=minindex,maxindex
              p0 => Mesh(index)
!
      
              if(p0%iFLG(1).gt.0) then 
                 iFLG0 = 4+p0%iFLG(1)
                 iFLG1 = p0%octPrt%iFLG(1)
                 if(iFLG1.le.6) then 
                    p0%octPrt%iFLG(1)=iFLG0
                 else if(iFLG0.lt.iFLG1) then 
                    p0%octPrt%iFLG(1)=iFLG0
                 end if

              endif
            
           end do
!!$omp end parallel do
!
        return
      end subroutine set_flag_prtG
!

!***********************************************************************
      subroutine add_oct(iLv)
! +-------------------------------------------------------------------+
! | Make eight child octs (level = iLv+1)                              |
! |                  at the flagged oct (level = iLv)                  |
! |                                                                    |
! |  Order of child index for octCh#1-8                                |
! |                                                                    |
! |         +-------+-------+                                          |
! |        /   7   /   8   /|                                          |
! |       +-------+-------+ |                                          |
! |      /   5   /   6   /|8|                                          |
! |     +-------+-------+ | +      Z                                   |
! |     |       |       |6|/|      |   Y                               |
! |     |   5   |   6   | + |      |  /                                |
! |     |       |       |/|4|      | /                                 |
! |     +-------+-------+ | +      |/                                  |
! |     |       |       |2|/       o--------X                          |
! |     |   1   |   2   | +                                            |
! |     |       |       |/                                             |
! |     +-------+-------+                                              |
! |                                                                    |
! |                                                        
! |          MinID(1,iLv)  = MinID(2, iLv) is for the internal task    !
! |  Mesh index  - one dimensional index -                             |
! |                                                                    |
! |              +-----+-----+-----+-----+-----+                       |
! |    Level=0   |  1  |  2  | ... | N-1 |  N  |   MinID(1,0) = 1      |
! |              +-----+-----+-----+-----+-----+   MaxID(1,0) = N      |
! |                                                                    |
! |              +-----+-----+-----+-----+-----+                       |
! |    Level=1   | N+1 | N+2 | ... | M-1 |  M  |   MinID(1,1) = N+1    |
! |              +-----+-----+-----+-----+-----+   MaxID(1,1) = M      |
! |                                                                    |
! |              +-----+-----+-----+-----+-----+                       |
! |    Level=2   | M+1 | M+2 | ... | L-1 |  L  |   MinID(1,2) = M+1    |
! |              +-----+-----+-----+-----+-----+   MaxID(1,2) = L      |
! |                                                                    |
! +-------------------------------------------------------------------+
!!***********************************************************************
        use oct_set
        use particle_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
! ------------------------------------------
        implicit none
        integer(kind=4) :: index,index2,ich,k,iLv,i
        integer(kind=4) :: octN0,octLv,octP,iFLG(3),Nint,Nint2,Csort
        integer(kind=4) :: Isort
        integer(kind=4) :: MrtN
        integer(kind=4) :: n1,n2,ix,iy,iz,iii
        integer(kind=4) :: iPOS(3)
        real(kind=8)    :: rPOS(3)
        real(kind=8)    :: ixx,iyy,izz
        real(kind=8)    :: F(18),C(6*iomp0),Z(IonSorts),R(9),G(1:18),D(1:3),O(1:16)
        integer(kind=4) :: AllocateStatus
        integer(kind=4) :: iC(iomp0)
        integer(kind=4) :: indexP, octType, prc_bndry, ptcl_loops
        real(kind=8)    :: rx,ry,rz,rr,mr
!-------------------------------------------
        type(oct), pointer :: p0,newp
        type(prtcl), pointer :: PrtList
!-------------------------------------------
        !if(debugMode>=1)print *,"start add_oct iLv=",iLv,rank
!        print*, 'add_oct iLv=',iLv
         allocate(newp,stat=AllocateStatus)

         Nint2 = 2**LvMax
         MinID(2,iLv+1)=max(MaxID(1,LvMax),MaxID(2,iLv))
         indexP = MaxIP(iLv+1)
         index2 = MinID(2,iLv+1)
         p0 => Mesh(index2)
         octN0=p0%octN

         Isort = -1
         R = 0.d0
         C = 0.d0
         Z = 0.d0
         G = 0.d0
         D=0.d0
         O=0.d0
         octP= 0
         iC  = 0
         octType = 0
         ptcl_loops=0
         iFLG=0

!         print*,'entering add oct:rank',rank,iLv
!         print*,'MinID, MaxID',MinID(1,iLv),MaxID(1,iLv)

         do i=1,2
            if(minID(i,iLv).lt.maxID(i,iLv)) then 
               do index=MinID(i,iLv),MaxID(i,iLv)
                  p0 => Mesh(index)
                  if(p0%iFLG(2) == 2)then
         
                     iFLG(2) = 0  
                     iFLG(1) = p0%iFLG(1)-4  
                     iFLG(3) = 0
                     octLv   = p0%octLv  +1
                     Nint  = 2**(LvMax-octLv)
!
! -- each directed child octs

                     do ich=1,8
                        nullify(newp%octPrt)
                        nullify(newp%octNb1) ; nullify(newp%octNb2)
                        nullify(newp%octNb3) ; nullify(newp%octNb4)
                        nullify(newp%octNb5) ; nullify(newp%octNb6)
                        nullify(newp%octCh1) ; nullify(newp%octCh2)
                        nullify(newp%octCh3) ; nullify(newp%octCh4)
                        nullify(newp%octCh5) ; nullify(newp%octCh6)
                        nullify(newp%octCh7) ; nullify(newp%octCh8)
                        nullify(newp%Psort)
                        nullify(newp%ptcl)   ; nullify(newp%ptclA)
!
                        nullify(PrtList)
!
                        index2 = index2 + 1
                        if(index2>=MinID(3,-1))then
                           print *,"==============================="
                           print *,"Mesh overflow index2=",index2,"size(Mesh)=",MinID(3,-1)-1,rank 
                           print *,"==============================="
                           stop
                        endif
                        MrtN   = p0%MrtN
! -- Extract position index for child-octs --
                        n1 = ich
                        iz = int((n1-1)/4) + 1
                        n2 = n1 - 4*(iz-1)
                        iy = int((n2-1)/2) + 1
                        ix = n2 - 2*(iy-1)
! -- Position index --
                        ix = p0%iPOS(1) + (2*ix-3)*Nint
                        iy = p0%iPOS(2) + (2*iy-3)*Nint
                        iz = p0%iPOS(3) + (2*iz-3)*Nint
                        iPOS(1) = ix
                        iPOS(2) = iy
                        iPOS(3) = iz

! -- Position --
                        ixx = real(ix + Nint2)/real(2*Nint2)
                        iyy = real(iy + Nint2)/real(2*Nint2)
                        izz = real(iz + Nint2)/real(2*Nint2)
                        rPOS(1) = R_lim(0,1)+ixx*dx(1) - HALF*dx(1)
                        rPOS(2) = R_lim(0,2)+iyy*dx(2) - HALF*dx(2)
                        rPOS(3) = R_lim(0,3)+izz*dx(3) - HALF*dx(3)
! -- Field values input --
                        do k=1,12
                           F(k) = p0%F(k)
                        enddo
 
                        do k=1,IonSorts
                           Z(k) = p0%Z(k)
                        enddo
                        
                        prc_bndry = p0% prc_bndry

                        Csort = ich
                        octN0 = index2
                        if(dipoleFlag==1)then
                            rx=rPOS(1)-Dpos(1)
                            ry=rPOS(2)-Dpos(2)
                            rz=rPOS(3)-Dpos(3)
                            rr=sqrt(rx**2+ry**2+rz**2)
                            if(rr.gt.dx(1)*2.0d0) then
                                mr=mx*rx+my*ry+mz*rz
                                D(1)=(-mx/(rr**3)+3*rx*mr/(rr**5))/(4*PI)
                                D(2)=(-my/(rr**3)+3*ry*mr/(rr**5))/(4*PI)
                                D(3)=(-mz/(rr**3)+3*rz*mr/(rr**5))/(4*PI)
                            endif
                        end if
!
! -- Connect parent oct & next oct --
                        newp%octPrt => p0
                        newp%Psort  => p0
! -- representative particle 
                        indexP = indexP + 1
                        R(1)=rPOS(1)
                        R(2)=rPOS(2)
                        R(3)=rPOS(3)
                        Pesh(indexP,iLv+1) =                  &
                         prtcl(R,Isort,octN0,indexP,PrtList,0)
                        newp%ptcl => Pesh(indexP,iLv+1)
! -- Input informaions into a target variable Mesh(:) --
                        Mesh(index2) =                        &
                         oct(octN0,octLv,octP,Csort,iFLG,     &
                             iPOS,rPOS,                       &
                             MrtN,iC,                    &
                             prc_bndry,                  &
                             octType,                      &
                             F,C,Z,G,D,O,ptcl_loops,    & 
                             newp%octPrt,                     &
                             newp%octNb1, newp%octNb2,        &
                             newp%octNb3, newp%octNb4,        &
                             newp%octNb5, newp%octNb6,        &
                             newp%octCh1, newp%octCh2,        &
                             newp%octCh3, newp%octCh4,        &
                             newp%octCh5, newp%octCh6,        &
                             newp%octCh7, newp%octCh8,        &
                             newp%Psort ,                     &
                             newp%ptcl  , newp%ptclA)

! -- Connect eight child-octs of parent oct --
                        if(ich==1) then
                           p0%octCh1 => Mesh(index2)
                        elseif(ich==2) then
                           p0%octCh2 => Mesh(index2)
                        elseif(ich==3) then
                           p0%octCh3 => Mesh(index2)
                        elseif(ich==4) then
                           p0%octCh4 => Mesh(index2)
                        elseif(ich==5) then
                           p0%octCh5 => Mesh(index2)
                        elseif(ich==6) then
                           p0%octCh6 => Mesh(index2)
                        elseif(ich==7) then
                           p0%octCh7 => Mesh(index2)
                        elseif(ich==8) then 
                           p0%octCh8 => Mesh(index2)
                        endif
                     enddo  ! for do ich=1,8
                     p0%iFLG(2) = 1
                  endif
               end do
            end if
!
         end do
!
         if(index2.eq.MinID(2,iLv+1)) then 
            MinID(2,iLv+1) = MaxID(2,iLv)
            MaxID(2,iLv+1) = MaxID(2,iLv)
         else
            MinID(2,iLv+1) = MaxID(2,iLv) + 1
            MaxID(2,iLv+1) = index2 
         end if
         maxIP(iLv+1)=indexP
!
         deallocate(newp,stat=iii)
! 
         !if(debugMode>=1)print *,"exit add_oct iLv=",iLv,rank
         return
         end subroutine add_oct
!
!***********************************************************************

!***********************************************************************
      subroutine add_Goct(iLv)
! +-------------------------------------------------------------------+
! | Make eight child octs (level = iLv+1)                              |
! |                  at the flagged oct (level = iLv)                  |
! |                                                                    |
! |  Order of child index for octCh#1-8                                |
! |                                                                    |
! |         +-------+-------+                                          |
! |        /   7   /   8   /|                                          |
! |       +-------+-------+ |                                          |
! |      /   5   /   6   /|8|                                          |
! |     +-------+-------+ | +      Z                                   |
! |     |       |       |6|/|      |   Y                               |
! |     |   5   |   6   | + |      |  /                                |
! |     |       |       |/|4|      | /                                 |
! |     +-------+-------+ | +      |/                                  |
! |     |       |       |2|/       o--------X                          |
! |     |   1   |   2   | +                                            |
! |     |       |       |/                                             |
! |     +-------+-------+                                              |
! |                                                                    |
! |                                                        
! |          MinID(1,iLv)  = MinID(2, iLv) is for the internal task    !
! |  Mesh index  - one dimensional index -                             |
! |                                                                    |
! |              +-----+-----+-----+-----+-----+                       |
! |    Level=0   |  1  |  2  | ... | N-1 |  N  |   MinID(1,0) = 1      |
! |              +-----+-----+-----+-----+-----+   MaxID(1,0) = N      |
! |                                                                    |
! |              +-----+-----+-----+-----+-----+                       |
! |    Level=1   | N+1 | N+2 | ... | M-1 |  M  |   MinID(1,1) = N+1    |
! |              +-----+-----+-----+-----+-----+   MaxID(1,1) = M      |
! |                                                                    |
! |              +-----+-----+-----+-----+-----+                       |
! |    Level=2   | M+1 | M+2 | ... | L-1 |  L  |   MinID(1,2) = M+1    |
! |              +-----+-----+-----+-----+-----+   MaxID(1,2) = L      |
! |                                                                    |
! +-------------------------------------------------------------------+
!!***********************************************************************
        use oct_set
        use particle_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
! ------------------------------------------
        implicit none
        integer(kind=4) :: index,index2,ich,k,iLv,i
        integer(kind=4) :: octN0,octLv,octP,iFLG(3),Nint,Nint2,Csort
        integer(kind=4) :: Isort
        integer(kind=4) :: MrtN
        integer(kind=4) :: n1,n2,ix,iy,iz,iii
        integer(kind=4) :: iPOS(3)
        real(kind=8)    :: rPOS(3)
        real(kind=8)    :: ixx,iyy,izz
        real(kind=8)    :: F(18),C(6*iomp0),Z(IonSorts),R(9),G(1:18),D(1:3),O(1:16)
        integer(kind=4) :: AllocateStatus
        integer(kind=4) :: iC(iomp0)
        integer(kind=4) :: indexP, octType, prc_bndry, ptcl_loops
        real(kind=8)    :: rx,ry,rz,rr,mr
!-------------------------------------------
        type(oct), pointer :: p0,newp
        type(prtcl), pointer :: PrtList
!-------------------------------------------
        !if(debugMode>=1)print *,"start add_Goct iLv=",iLv,rank
!        print*, 'add_oct iLv=',iLv
         allocate(newp,stat=AllocateStatus)

         Nint2 = 2**LvMax
         MinID(4,iLv+1)=max(MaxID(3,LvMax),MaxID(4,iLv))
         indexP = MaxIP(iLv+1)
         index2=MinID(4,iLv+1)
         p0 => Mesh(index2)
         octN0=p0%octN

         Isort = -1
         R = 0.d0
         C = 0.d0
         Z = 0.d0
         G = 0.d0
         D=0.d0
         O=0.d0
         octP= 0
         iC  = 0
         octType = 1 !== Goct
         ptcl_loops=0
         iFLG=0

!         print*,'entering add oct:rank',rank,iLv
!         print*,'MinID, MaxID',MinID(1,iLv),MaxID(1,iLv)

         do i=3,4
            if(minID(i,iLv).lt.maxID(i,iLv)) then 
               do index=MinID(i,iLv),MaxID(i,iLv)
                  p0 => Mesh(index)
                  if(p0%iFLG(2) == 2) then 

                     iFLG(2) = 0
                     iFLG(1) = min(p0%iFLG(1),4)-4  
                     iFLG(3) = 0
                     octLv   = p0%octLv  +1
                     Nint  = 2**(LvMax-octLv)
!
! -- each directed child octs

                     do ich=1,8
                        nullify(newp%octPrt)
                        nullify(newp%octNb1) ; nullify(newp%octNb2)
                        nullify(newp%octNb3) ; nullify(newp%octNb4)
                        nullify(newp%octNb5) ; nullify(newp%octNb6)
                        nullify(newp%octCh1) ; nullify(newp%octCh2)
                        nullify(newp%octCh3) ; nullify(newp%octCh4)
                        nullify(newp%octCh5) ; nullify(newp%octCh6)
                        nullify(newp%octCh7) ; nullify(newp%octCh8)
                        nullify(newp%Psort)
                        nullify(newp%ptcl)   ; nullify(newp%ptclA)
!
                        nullify(PrtList)
!
                        index2 = index2 + 1

                        if(index2>size(Mesh))then
                           print *,"==============================="
                           print *,"GMesh overflow index2=",index2,"size(GMesh)=",size(Mesh),rank 
                           print *,"==============================="
                           stop
                        endif

                        MrtN   = p0%MrtN
! -- Extract position index for child-octs --
                        n1 = ich
                        iz = int((n1-1)/4) + 1
                        n2 = n1 - 4*(iz-1)
                        iy = int((n2-1)/2) + 1
                        ix = n2 - 2*(iy-1)
! -- Position index --
                        ix = p0%iPOS(1) + (2*ix-3)*Nint
                        iy = p0%iPOS(2) + (2*iy-3)*Nint
                        iz = p0%iPOS(3) + (2*iz-3)*Nint
                        iPOS(1) = ix
                        iPOS(2) = iy
                        iPOS(3) = iz

! -- Position --
                        ixx = real(ix + Nint2)/real(2*Nint2)
                        iyy = real(iy + Nint2)/real(2*Nint2)
                        izz = real(iz + Nint2)/real(2*Nint2)
                        rPOS(1) = R_lim(0,1)+ixx*dx(1) - HALF*dx(1)
                        rPOS(2) = R_lim(0,2)+iyy*dx(2) - HALF*dx(2)
                        rPOS(3) = R_lim(0,3)+izz*dx(3) - HALF*dx(3)
! -- Field values input --
                        do k=1,12
                           F(k) = p0%F(k)
                        enddo
 
                        do k=1,IonSorts
                           Z(k) = p0%Z(k)
                        enddo
                        
                        prc_bndry = p0% prc_bndry

                        Csort = ich
                        octN0 = index2
                        
                        if(dipoleFlag==1)then
                            rx=rPOS(1)-Dpos(1)
                            ry=rPOS(2)-Dpos(2)
                            rz=rPOS(3)-Dpos(3)
                            rr=sqrt(rx**2+ry**2+rz**2)
                            if(rr.gt.dx(1)*2.0d0) then
                                mr=mx*rx+my*ry+mz*rz
                                D(1)=(-mx/(rr**3)+3*rx*mr/(rr**5))/(4*PI)
                                D(2)=(-my/(rr**3)+3*ry*mr/(rr**5))/(4*PI)
                                D(3)=(-mz/(rr**3)+3*rz*mr/(rr**5))/(4*PI)
                            endif
                        end if
!
! -- Connect parent oct & next oct --
                        newp%octPrt => p0
                        newp%Psort  => p0
!
! -- representative particle 
                        indexP = indexP + 1
                        R(1)=rPOS(1)
                        R(2)=rPOS(2)
                        R(3)=rPOS(3)
                        Pesh(indexP,iLv+1) =                  &
                         prtcl(R,Isort,octN0,indexP,PrtList,0)
                        newp%ptcl => Pesh(indexP,iLv+1)

! -- Input informaions into a target variable Mesh(:) --
                        Mesh(index2) =                        &
                         oct(octN0,octLv,octP,Csort,iFLG,     &
                             iPOS,rPOS,                       &
                             MrtN,iC,                    &
                             prc_bndry,                  &
                             octType,                      &
                             F,C,Z,G,D,O,ptcl_loops,    & 
                             newp%octPrt,                     &
                             newp%octNb1, newp%octNb2,        &
                             newp%octNb3, newp%octNb4,        &
                             newp%octNb5, newp%octNb6,        &
                             newp%octCh1, newp%octCh2,        &
                             newp%octCh3, newp%octCh4,        &
                             newp%octCh5, newp%octCh6,        &
                             newp%octCh7, newp%octCh8,        &
                             newp%Psort ,                     &
                             newp%ptcl  , newp%ptclA)
!
! -- Connect eight child-octs of parent oct --
                        if(ich==1) then
                           p0%octCh1 => Mesh(index2)
                        elseif(ich==2) then
                           p0%octCh2 => Mesh(index2)
                        elseif(ich==3) then
                           p0%octCh3 => Mesh(index2)
                        elseif(ich==4) then
                           p0%octCh4 => Mesh(index2)
                        elseif(ich==5) then
                           p0%octCh5 => Mesh(index2)
                        elseif(ich==6) then
                           p0%octCh6 => Mesh(index2)
                        elseif(ich==7) then
                           p0%octCh7 => Mesh(index2)
                        elseif(ich==8) then 
                           p0%octCh8 => Mesh(index2)
                        endif
!
                     enddo  ! for do ich=1,8
!
                     p0%iFLG(2) = 1
                  endif
               end do
            end if
!
         end do
!
         if(index2.eq.MinID(4,iLv+1)) then 
            MinID(4,iLv+1) = MaxID(4,iLv)
            MaxID(4,iLv+1) = MaxID(4,iLv)
         else
            MinID(4,iLv+1) = MaxID(4,iLv) + 1
            MaxID(4,iLv+1) = index2 
         end if
         maxIP(iLv+1)=indexP
!
         deallocate(newp,stat=iii)
! 
        ! if(debugMode>=1)print *,"exit add_Goct iLv=",iLv,rank
         return
       end subroutine add_Goct
!
!***********************************************************************
      subroutine set_buffer_out(iLv,iF1,iF2) !need to modify
! +-------------------------------------------------------------------+
! |     if iFLG >= iF1 ; iFLG of the neighbor octs should be > iF2    | 
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
! ------------------------------------------
        implicit none
        integer(kind=4) :: iF1,iF2
        integer(kind=4) :: i,iLv,index,maxindex,minindex  
!-------------------------------------------
        type(oct), pointer :: p0,p1
!-------------------------------------------
        !if(debugMode>=1)print *,"entering set buffer out iLv=",iLv,rank

        do i=1,4
        maxindex = MaxID(i,iLv)
        minindex = MinID(i,iLv)
        if(maxindex.gt.minindex) then 
!!$omp parallel do private(index,p0,p1)  
           do index=minindex,maxindex
              p0 => Mesh(index)

              if(p0%iFLG(1) >= iF1 ) then
                 if(p0%octNb1%iFLG(1) < iF2) then
                    p0%octNb1%iFLG(1) = iF2
                 endif
                 if(p0%octNb2%iFLG(1) < iF2) then
                    p0%octNb2%iFLG(1) = iF2
                 endif
                 if(p0%octNb3%iFLG(1) < iF2) then
                    p0%octNb3%iFLG(1) = iF2
                 endif
                 if(p0%octNb4%iFLG(1) < iF2) then
                    p0%octNb4%iFLG(1) = iF2
                 endif
                 if(p0%octNb1%octNb3%iFLG(1) < iF2) then
                    p0%octNb1%octNb3%iFLG(1) = iF2
                 endif
                 if(p0%octNb2%octNb4%iFLG(1) < iF2) then
                    p0%octNb2%octNb4%iFLG(1) = iF2
                 endif
                 if(p0%octNb3%octNb2%iFLG(1) < iF2) then
                    p0%octNb3%octNb2%iFLG(1) = iF2
                 endif
                 if(p0%octNb4%octNb1%iFLG(1) < iF2) then
                    p0%octNb4%octNb1%iFLG(1) = iF2
                 endif
!--------------------------------------
                 p1 => p0%octNb5
!--------------------------------------
                 if(p1%iFLG(1) < iF2) then
                    p1%iFLG(1) = iF2
                 endif
                 if(p1%octNb1%iFLG(1) < iF2) then
                    p1%octNb1%iFLG(1) = iF2
                 endif
                 if(p1%octNb2%iFLG(1) < iF2) then
                    p1%octNb2%iFLG(1) = iF2
                 endif
                 if(p1%octNb3%iFLG(1) < iF2) then
                    p1%octNb3%iFLG(1) = iF2
                 endif
                 if(p1%octNb4%iFLG(1) < iF2) then
                    p1%octNb4%iFLG(1) = iF2
                 endif
                 if(p1%octNb1%octNb3%iFLG(1) < iF2) then
                    p1%octNb1%octNb3%iFLG(1) = iF2
                 endif
                 if(p1%octNb2%octNb4%iFLG(1) < iF2) then
                    p1%octNb2%octNb4%iFLG(1) = iF2
                 endif
                 if(p1%octNb3%octNb2%iFLG(1) < iF2) then
                    p1%octNb3%octNb2%iFLG(1) = iF2
                 endif
                 if(p1%octNb4%octNb1%iFLG(1) < iF2) then
                    p1%octNb4%octNb1%iFLG(1) = iF2
                 endif
!--------------------------------------
                 p1 => p0%octNb6
!--------------------------------------
                 if(p1%iFLG(1) < iF2) then
                    p1%iFLG(1) = iF2
                 endif
                 if(p1%octNb1%iFLG(1) < iF2) then
                    p1%octNb1%iFLG(1) = iF2
                 endif
                 if(p1%octNb2%iFLG(1) < iF2) then
                    p1%octNb2%iFLG(1) = iF2
                 endif
                 if(p1%octNb3%iFLG(1) < iF2) then
                    p1%octNb3%iFLG(1) = iF2
                 endif
                 if(p1%octNb4%iFLG(1) < iF2) then
                    p1%octNb4%iFLG(1) = iF2
                 endif
                 if(p1%octNb1%octNb3%iFLG(1) < iF2) then
                    p1%octNb1%octNb3%iFLG(1) = iF2
                 endif
                 if(p1%octNb2%octNb4%iFLG(1) < iF2) then
                    p1%octNb2%octNb4%iFLG(1) = iF2
                 endif
                 if(p1%octNb3%octNb2%iFLG(1) < iF2) then
                    p1%octNb3%octNb2%iFLG(1) = iF2
                 endif
                 if(p1%octNb4%octNb1%iFLG(1) < iF2) then
                    p1%octNb4%octNb1%iFLG(1) = iF2
                 endif

              endif
           end do
!!$omp end parallel do
        end if
        end do
!

        !if(debugMode>=1)print *,"exiting set buffer out iLv=",iLv,rank
        return
        end subroutine set_buffer_out

!***********************************************************************
      subroutine set_buffer_prt(iLv,iF1,iF2)
! +-------------------------------------------------------------------+
! |     if iFLG >= iF1 ; iFLG of the neighbor octs should be > iF2    | 
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
! ------------------------------------------
        implicit none
        integer(kind=4) :: iF1,iF2
        integer(kind=4) :: i,iLv,index,maxindex,minindex  
!-------------------------------------------
        type(oct), pointer :: p0,p1
!-------------------------------------------

        do i=1,2
        maxindex = MaxID(i,iLv-1)
        minindex = MinID(i,iLv-1)
        if(maxindex.gt.minindex) then 
!!$omp parallel do private(index,p0,p1)  
           do index=minindex,maxindex
              p0 => Mesh(index)
              if(p0%iFLG(3) >= iF1 ) then
                 if(p0%octNb1%iFLG(3) < iF2) then
                    p0%octNb1%iFLG(3) = iF2
                 endif
                 if(p0%octNb2%iFLG(3) < iF2) then
                    p0%octNb2%iFLG(3) = iF2
                 endif
                 if(p0%octNb3%iFLG(3) < iF2) then
                    p0%octNb3%iFLG(3) = iF2
                 endif
                 if(p0%octNb4%iFLG(3) < iF2) then
                    p0%octNb4%iFLG(3) = iF2
                 endif
                 if(p0%octNb1%octNb3%iFLG(3) < iF2) then
                    p0%octNb1%octNb3%iFLG(3) = iF2
                 endif
                 if(p0%octNb2%octNb4%iFLG(3) < iF2) then
                    p0%octNb2%octNb4%iFLG(3) = iF2
                 endif
                 if(p0%octNb3%octNb2%iFLG(3) < iF2) then
                    p0%octNb3%octNb2%iFLG(3) = iF2
                 endif
                 if(p0%octNb4%octNb1%iFLG(3) < iF2) then
                    p0%octNb4%octNb1%iFLG(3) = iF2
                 endif
!--------------------------------------
                 p1 => p0%octNb5
!--------------------------------------
                 if(p1%iFLG(3) < iF2) then
                    p1%iFLG(3) = iF2
                 endif
                 if(p1%octNb1%iFLG(3) < iF2) then
                    p1%octNb1%iFLG(3) = iF2
                 endif
                 if(p1%octNb2%iFLG(3) < iF2) then
                    p1%octNb2%iFLG(3) = iF2
                 endif
                 if(p1%octNb3%iFLG(3) < iF2) then
                    p1%octNb3%iFLG(3) = iF2
                 endif
                 if(p1%octNb4%iFLG(3) < iF2) then
                    p1%octNb4%iFLG(3) = iF2
                 endif
                 if(p1%octNb1%octNb3%iFLG(3) < iF2) then
                    p1%octNb1%octNb3%iFLG(3) = iF2
                 endif
                 if(p1%octNb2%octNb4%iFLG(3) < iF2) then
                    p1%octNb2%octNb4%iFLG(3) = iF2
                 endif
                 if(p1%octNb3%octNb2%iFLG(3) < iF2) then
                    p1%octNb3%octNb2%iFLG(3) = iF2
                 endif
                 if(p1%octNb4%octNb1%iFLG(3) < iF2) then
                    p1%octNb4%octNb1%iFLG(3) = iF2
                 endif
!--------------------------------------
                 p1 => p0%octNb6
!--------------------------------------
                 if(p1%iFLG(3) < iF2) then
                    p1%iFLG(3) = iF2
                 endif
                 if(p1%octNb1%iFLG(3) < iF2) then
                    p1%octNb1%iFLG(3) = iF2
                 endif
                 if(p1%octNb2%iFLG(3) < iF2) then
                    p1%octNb2%iFLG(3) = iF2
                 endif
                 if(p1%octNb3%iFLG(3) < iF2) then
                    p1%octNb3%iFLG(3) = iF2
                 endif
                 if(p1%octNb4%iFLG(3) < iF2) then
                    p1%octNb4%iFLG(3) = iF2
                 endif
                 if(p1%octNb1%octNb3%iFLG(3) < iF2) then
                    p1%octNb1%octNb3%iFLG(3) = iF2
                 endif
                 if(p1%octNb2%octNb4%iFLG(3) < iF2) then
                    p1%octNb2%octNb4%iFLG(3) = iF2
                 endif
                 if(p1%octNb3%octNb2%iFLG(3) < iF2) then
                    p1%octNb3%octNb2%iFLG(3) = iF2
                 endif
                 if(p1%octNb4%octNb1%iFLG(3) < iF2) then
                    p1%octNb4%octNb1%iFLG(3) = iF2
                 endif

              endif
           end do
!!$omp end parallel do
        end if
        end do
!
        return
        end subroutine set_buffer_prt
!
!************************************************************************
      subroutine set_buffer_in(iLv,iF1,iF2)
! +---------------------------------------------------------------------+
! |  if iFLG of neighbor octs <= iF2 : iFLG of the oct should be <= iF1 | 
! +---------------------------------------------------------------------+
!************************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
! ------------------------------------------
        implicit none
        integer(kind=4) :: iF1,iF2,iFLG0,iFLG1,iFLG2,iFLG3
        integer(kind=4) :: i,iLv,index,maxindex,minindex
!-------------------------------------------
        type(oct), pointer :: p0,p1
!-------------------------------------------
!
! This routine sets iFLG(1) = iF1 if there is a neighbouring cell with iFLG(1) <= iF2

        do i=1,4
        maxindex = MaxID(i,iLv)
        minindex = MinID(i,iLv)
        if(maxindex.gt.minindex) then 
!!$omp parallel do private(index,p0,p1,iFLG0,iFLG1,iFLG2,iFLG3)  
        do index=minindex,maxindex
          p0 => Mesh(index)
          if(p0%iFLG(1) > iF1 ) then

             iFLG2 = min(p0%octNb1%octNb3%iFLG(1),p0%octNb2%octNb4%iFLG(1))
             iFLG3 = min(p0%octNb3%octNb2%iFLG(1),p0%octNb4%octNb1%iFLG(1))
             iFLG2 = min(iFLG2,iFLG3)
             if(iFLG2 <= iF2) then 
                p0%iFLG(1) = iF1
             else
                iFLG0 = min(p0%octNb1%iFLG(1),p0%octNb2%iFLG(1))
                iFLG1 = min(p0%octNb3%iFLG(1),p0%octNb4%iFLG(1))
                iFLG0 = min(iFLG0,iFLG1)
                if(iFLG0 <= iF2)  p0%iFLG(1) = iF1
             end if

             if(p0%iFLG(1) > iF1) then 
                p1 => p0%octNb5
                iFLG2 = min(p1%octNb1%octNb3%iFLG(1),p1%octNb2%octNb4%iFLG(1))
                iFLG3 = min(p1%octNb3%octNb2%iFLG(1),p1%octNb4%octNb1%iFLG(1))
                iFLG2 = min(iFLG2,iFLG3)
                if(iFLG2 <= iF2) then 
                   p0%iFLG(1) = iF1
                else
                   iFLG0 = min(p1%octNb1%iFLG(1),p1%octNb2%iFLG(1))
                   iFLG1 = min(p1%octNb3%iFLG(1),p1%octNb4%iFLG(1))
                   iFLG0 = min(iFLG0,iFLG1)
                   iFLG0 = min(iFLG0,p1%iFLG(1))
                   if(iFLG0 <= iF2)  p0%iFLG(1) = iF1
                end if
             end if

             if(p0%iFLG(1) > iF1) then 
                p1 => p0%octNb6
                iFLG2 = min(p1%octNb1%octNb3%iFLG(1),p1%octNb2%octNb4%iFLG(1))
                iFLG3 = min(p1%octNb3%octNb2%iFLG(1),p1%octNb4%octNb1%iFLG(1))
                iFLG2 = min(iFLG2,iFLG3)
                if(iFLG2 <= iF2) then 
                   p0%iFLG(1) = iF1
                else
                   iFLG0 = min(p1%octNb1%iFLG(1),p1%octNb2%iFLG(1))
                   iFLG1 = min(p1%octNb3%iFLG(1),p1%octNb4%iFLG(1))
                   iFLG0 = min(iFLG0,iFLG1)
                   iFLG0 = min(iFLG0,p1%iFLG(1))
                   if(iFLG0 <= iF2)  p0%iFLG(1) = iF1
                end if
             end if

          endif
        end do
!!$omp end parallel do
     end if
        end do
!
        return
        end subroutine set_buffer_in

!************************************************************************
      subroutine set_buffer_ave(iLv,iF1,iF2)
! +---------------------------------------------------------------------+
! |  if iFLG of all neighbor octs >= iF2 : iFLG of should be = iF1      | 
! +---------------------------------------------------------------------+
!************************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
! ------------------------------------------
        implicit none
        integer(kind=4) :: iF1,iF2
        integer(kind=4) :: i,iLv,index,maxindex,minindex
!-------------------------------------------
        type(oct), pointer :: p0
!-------------------------------------------
! This routine sets iFLG(1) = iF2 for those of iFLG(1) = iF1  
! but keep as it is if there is a neighbouring cell with iFLG<iF1
! For example if p0 % iFLG(1)= 4 then change it to p0 % iFLG(1) = 6 only 
! if there is a neighbour with iFLG(1) > 4

        do i=1,4
        maxindex = MaxID(i,iLv)
        minindex = MinID(i,iLv)
        if(maxindex.gt.minindex) then 
!!$omp parallel do private(index,p0)  
        do index=minindex,maxindex
          p0 => Mesh(index)
           if(p0%iFLG(1) == iF1 ) then
              p0%iFLG(1) = iF2
              if(p0%octNb1%iFLG(1) < iF1) then
                 p0%iFLG(1) = iF1
              endif
              if(p0%octNb2%iFLG(1) < iF1) then
                 p0%iFLG(1) = iF1
              endif
              if(p0%octNb3%iFLG(1) < iF1) then
                 p0%iFLG(1) = iF1
              endif
              if(p0%octNb4%iFLG(1) < iF1) then
                 p0%iFLG(1) = iF1
              endif
              if(p0%octNb5%iFLG(1) < iF1) then
                 p0%iFLG(1) = iF1
              endif
              if(p0%octNb6%iFLG(1) < iF1) then
                 p0%iFLG(1) = iF1
              endif
           endif
        end do
!!$omp end parallel do
     end if
  end do
  !
  return
end subroutine set_buffer_ave

!***********************************************************************
      subroutine connect_oct(iLv)
! +-------------------------------------------------------------------+
! |     create connections for additinal octs                         |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4)    :: iLv,i,index
        type(oct), pointer :: p0

        !if(debugMode>=2)print *,"start connect_oct iLv=",iLv,rank
        
        if(iLv.ge.LvMax .or. iLv<0) return 
        do i=1,2
           if(maxID(i,iLv).gt.minID(i,iLv)) then 
!!$omp parallel do private(index,p0)
              do index=minID(i,iLv),maxID(i,iLv)
                 p0 => Mesh(index)
                 !if(p0%iFLG(2)>=1) then 
                 if(p0%iFLG(1)>0)then
                    !initialize connection
                    !for octCh1
                    nullify(p0%octCh1%octNb1)
                    nullify(p0%octCh1%octNb2)
                    nullify(p0%octCh1%octNb3)
                    nullify(p0%octCh1%octNb4)
                    nullify(p0%octCh1%octNb5)
                    nullify(p0%octCh1%octNb6)
                    !for octCh2
                    nullify(p0%octCh2%octNb1)
                    nullify(p0%octCh2%octNb2)
                    nullify(p0%octCh2%octNb3)
                    nullify(p0%octCh2%octNb4)
                    nullify(p0%octCh2%octNb5)
                    nullify(p0%octCh2%octNb6)
                    !for octCh3
                    nullify(p0%octCh3%octNb1)
                    nullify(p0%octCh3%octNb2)
                    nullify(p0%octCh3%octNb3)
                    nullify(p0%octCh3%octNb4)
                    nullify(p0%octCh3%octNb5)
                    nullify(p0%octCh3%octNb6)
                    !for octCh4
                    nullify(p0%octCh4%octNb1)
                    nullify(p0%octCh4%octNb2)
                    nullify(p0%octCh4%octNb3)
                    nullify(p0%octCh4%octNb4)
                    nullify(p0%octCh4%octNb5)
                    nullify(p0%octCh4%octNb6)
                    !for octCh5
                    nullify(p0%octCh5%octNb1)
                    nullify(p0%octCh5%octNb2)
                    nullify(p0%octCh5%octNb3)
                    nullify(p0%octCh5%octNb4)
                    nullify(p0%octCh5%octNb5)
                    nullify(p0%octCh5%octNb6)
                    !for octCh6
                    nullify(p0%octCh6%octNb1)
                    nullify(p0%octCh6%octNb2)
                    nullify(p0%octCh6%octNb3)
                    nullify(p0%octCh6%octNb4)
                    nullify(p0%octCh6%octNb5)
                    nullify(p0%octCh6%octNb6)
                    !for octCh7
                    nullify(p0%octCh7%octNb1)
                    nullify(p0%octCh7%octNb2)
                    nullify(p0%octCh7%octNb3)
                    nullify(p0%octCh7%octNb4)
                    nullify(p0%octCh7%octNb5)
                    nullify(p0%octCh7%octNb6)
                    !for octCh8
                    nullify(p0%octCh8%octNb1)
                    nullify(p0%octCh8%octNb2)
                    nullify(p0%octCh8%octNb3)
                    nullify(p0%octCh8%octNb4)
                    nullify(p0%octCh8%octNb5)
                    nullify(p0%octCh8%octNb6)

                    if(p0%iFLG(1)>=2) then 
! - for octCh1 -
                       p0%octCh1%octNb1 => p0%octNb1%octCh2 
                       p0%octCh1%octNb2 => p0       %octCh2
                       p0%octCh1%octNb3 => p0%octNb3%octCh3
                       p0%octCh1%octNb4 => p0       %octCh3
                       p0%octCh1%octNb5 => p0%octNb5%octCh5
                       p0%octCh1%octNb6 => p0       %octCh5
! - for octCh2 -
                       p0%octCh2%octNb1 => p0       %octCh1
                       p0%octCh2%octNb2 => p0%octNb2%octCh1
                       p0%octCh2%octNb3 => p0%octNb3%octCh4
                       p0%octCh2%octNb4 => p0       %octCh4
                       p0%octCh2%octNb5 => p0%octNb5%octCh6
                       p0%octCh2%octNb6 => p0       %octCh6
! - for octCh3 -
                       p0%octCh3%octNb1 => p0%octNb1%octCh4
                       p0%octCh3%octNb2 => p0       %octCh4
                       p0%octCh3%octNb3 => p0       %octCh1
                       p0%octCh3%octNb4 => p0%octNb4%octCh1
                       p0%octCh3%octNb5 => p0%octNb5%octCh7
                       p0%octCh3%octNb6 => p0       %octCh7
! - for octCh4 -
                       p0%octCh4%octNb1 => p0       %octCh3
                       p0%octCh4%octNb2 => p0%octNb2%octCh3
                       p0%octCh4%octNb3 => p0       %octCh2
                       p0%octCh4%octNb4 => p0%octNb4%octCh2
                       p0%octCh4%octNb5 => p0%octNb5%octCh8
                       p0%octCh4%octNb6 => p0       %octCh8
! - for octCh5 -
                       p0%octCh5%octNb1 => p0%octNb1%octCh6 
                       p0%octCh5%octNb2 => p0       %octCh6
                       p0%octCh5%octNb3 => p0%octNb3%octCh7
                       p0%octCh5%octNb4 => p0       %octCh7
                       p0%octCh5%octNb5 => p0       %octCh1
                       p0%octCh5%octNb6 => p0%octNb6%octCh1
! - for octCh6 -
                       p0%octCh6%octNb1 => p0       %octCh5
                       p0%octCh6%octNb2 => p0%octNb2%octCh5
                       p0%octCh6%octNb3 => p0%octNb3%octCh8
                       p0%octCh6%octNb4 => p0       %octCh8
                       p0%octCh6%octNb5 => p0       %octCh2
                       p0%octCh6%octNb6 => p0%octNb6%octCh2
! - for octCh7 -
                       p0%octCh7%octNb1 => p0%octNb1%octCh8
                       p0%octCh7%octNb2 => p0       %octCh8
                       p0%octCh7%octNb3 => p0       %octCh5
                       p0%octCh7%octNb4 => p0%octNb4%octCh5
                       p0%octCh7%octNb5 => p0       %octCh3
                       p0%octCh7%octNb6 => p0%octNb6%octCh3
! - for octCh8 -
                       p0%octCh8%octNb1 => p0       %octCh7
                       p0%octCh8%octNb2 => p0%octNb2%octCh7
                       p0%octCh8%octNb3 => p0       %octCh6
                       p0%octCh8%octNb4 => p0%octNb4%octCh6
                       p0%octCh8%octNb5 => p0       %octCh4
                       p0%octCh8%octNb6 => p0%octNb6%octCh4
                    endif
                    if(p0%iFLG(1)==1) then
! - for octCh1 -
                       p0%octCh1%octNb1 => p0%octCh1
                       p0%octCh1%octNb2 => p0%octCh2 !
                       p0%octCh1%octNb3 => p0%octCh1 
                       p0%octCh1%octNb4 => p0%octCh3 !
                       p0%octCh1%octNb5 => p0%octCh1 
                       p0%octCh1%octNb6 => p0%octCh5 !
! - for octCh2 -
                       p0%octCh2%octNb1 => p0%octCh1 !
                       p0%octCh2%octNb2 => p0%octCh2
                       p0%octCh2%octNb3 => p0%octCh2
                       p0%octCh2%octNb4 => p0%octCh4 !
                       p0%octCh2%octNb5 => p0%octCh2
                       p0%octCh2%octNb6 => p0%octCh6 !
! - for octCh3 -
                       p0%octCh3%octNb1 => p0%octCh3 
                       p0%octCh3%octNb2 => p0%octCh4 !
                       p0%octCh3%octNb3 => p0%octCh1 !
                       p0%octCh3%octNb4 => p0%octCh3
                       p0%octCh3%octNb5 => p0%octCh3
                       p0%octCh3%octNb6 => p0%octCh7 !
! - for octCh4 -
                       p0%octCh4%octNb1 => p0%octCh3 !
                       p0%octCh4%octNb2 => p0%octCh4
                       p0%octCh4%octNb3 => p0%octCh2 !
                       p0%octCh4%octNb4 => p0%octCh4
                       p0%octCh4%octNb5 => p0%octCh4
                       p0%octCh4%octNb6 => p0%octCh8 !
! - for octCh5 -
                       p0%octCh5%octNb1 => p0%octCh5   
                       p0%octCh5%octNb2 => p0%octCh6 !
                       p0%octCh5%octNb3 => p0%octCh5
                       p0%octCh5%octNb4 => p0%octCh7 !
                       p0%octCh5%octNb5 => p0%octCh1 !
                       p0%octCh5%octNb6 => p0%octCh5
! - for octCh6 -
                       p0%octCh6%octNb1 => p0%octCh5 !
                       p0%octCh6%octNb2 => p0%octCh6
                       p0%octCh6%octNb3 => p0%octCh6
                       p0%octCh6%octNb4 => p0%octCh8 !
                       p0%octCh6%octNb5 => p0%octCh2 !
                       p0%octCh6%octNb6 => p0%octCh6
! - for octCh7 -
                       p0%octCh7%octNb1 => p0%octCh7 
                       p0%octCh7%octNb2 => p0%octCh8 !
                       p0%octCh7%octNb3 => p0%octCh5 !
                       p0%octCh7%octNb4 => p0%octCh7
                       p0%octCh7%octNb5 => p0%octCh3 !
                       p0%octCh7%octNb6 => p0%octCh7
! - for octCh8 -
                       p0%octCh8%octNb1 => p0%octCh7 !
                       p0%octCh8%octNb2 => p0%octCh8
                       p0%octCh8%octNb3 => p0%octCh6 !
                       p0%octCh8%octNb4 => p0%octCh8
                       p0%octCh8%octNb5 => p0%octCh4 !
                       p0%octCh8%octNb6 => p0%octCh8
!
                       if(p0%octNb1%iFLG(1) >= 1) then
                          p0%octCh1%octNb1 => p0%octNb1%octCh2
                          p0%octCh3%octNb1 => p0%octNb1%octCh4
                          p0%octCh5%octNb1 => p0%octNb1%octCh6
                          p0%octCh7%octNb1 => p0%octNb1%octCh8
                       endif
                       if(p0%octNb2%iFLG(1) >= 1) then
                          p0%octCh2%octNb2 => p0%octNb2%octCh1
                          p0%octCh4%octNb2 => p0%octNb2%octCh3
                          p0%octCh6%octNb2 => p0%octNb2%octCh5
                          p0%octCh8%octNb2 => p0%octNb2%octCh7
                       endif
                       if(p0%octNb3%iFLG(1) >= 1) then
                          p0%octCh1%octNb3 => p0%octNb3%octCh3
                          p0%octCh2%octNb3 => p0%octNb3%octCh4
                          p0%octCh5%octNb3 => p0%octNb3%octCh7
                          p0%octCh6%octNb3 => p0%octNb3%octCh8
                       endif
                       if(p0%octNb4%iFLG(1) >= 1) then
                          p0%octCh3%octNb4 => p0%octNb4%octCh1
                          p0%octCh4%octNb4 => p0%octNb4%octCh2
                          p0%octCh7%octNb4 => p0%octNb4%octCh5
                          p0%octCh8%octNb4 => p0%octNb4%octCh6
                       endif
                       if(p0%octNb5%iFLG(1) >= 1) then
                          p0%octCh1%octNb5 => p0%octNb5%octCh5
                          p0%octCh2%octNb5 => p0%octNb5%octCh6
                          p0%octCh3%octNb5 => p0%octNb5%octCh7
                          p0%octCh4%octNb5 => p0%octNb5%octCh8
                       endif
                       if(p0%octNb6%iFLG(1) >= 1) then
                          p0%octCh5%octNb6 => p0%octNb6%octCh1
                          p0%octCh6%octNb6 => p0%octNb6%octCh2
                          p0%octCh7%octNb6 => p0%octNb6%octCh3
                          p0%octCh8%octNb6 => p0%octNb6%octCh4
                       endif
                    endif
                    !p0%iFLG(2)=1
                 endif
              end do
!!$omp end parallel do
           endif
        end do
     
        !if(debugMode>=2)print *,"end connect_oct iLv=",iLv,rank

        return
        end subroutine connect_oct
!
!***********************************************************************
      subroutine connect_Goct(iLv)
! +-------------------------------------------------------------------+
! |     create connections for additinal octs                         |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4)    :: iLv,i,index,myoctN
        type(oct), pointer :: p0

        !if(debugMode>=2)print *,"start connect_Goct iLv=",iLv,rank

        if(iLv.ge.LvMax) return 
        if(iLv<=-1)return

        do i=3,4
           if(maxID(i,iLv).gt.minID(i,iLv)) then 
!!$omp parallel do private(index,p0,myoctN)
              do index=minID(i,iLv),maxID(i,iLv)
                 p0 => Mesh(index)
                 myoctN=p0%octN
                 !if(p0%iFLG(2)>=1) then
                 if(p0%iFLG(1)>0)then
                    if(p0%iFLG(1)>=2) then ! This means that the cell is not at the edge.
                       !! - for octNb1 -                   
                       if(p0%octNb1%octN==myoctN)then
                          p0%octCh1%octNb1 => p0%octCh1 
                          p0%octCh3%octNb1 => p0%octCh3 
                          p0%octCh5%octNb1 => p0%octCh5   
                          p0%octCh7%octNb1 => p0%octCh7 
                       else
                          p0%octCh1%octNb1 => p0%octNb1%octCh2     
                          p0%octCh3%octNb1 => p0%octNb1%octCh4   
                          p0%octCh5%octNb1 => p0%octNb1%octCh6 
                          p0%octCh7%octNb1 => p0%octNb1%octCh8
                       endif
                       ! - for octNb2 -                
                       if(p0%octNb2%octN==myoctN)then
                          p0%octCh2%octNb2 => p0%octCh2
                          p0%octCh4%octNb2 => p0%octCh4
                          p0%octCh6%octNb2 => p0%octCh6
                          p0%octCh8%octNb2 => p0%octCh8
                       else
                          p0%octCh2%octNb2 => p0%octNb2%octCh1
                          p0%octCh4%octNb2 => p0%octNb2%octCh3
                          p0%octCh6%octNb2 => p0%octNb2%octCh5
                          p0%octCh8%octNb2 => p0%octNb2%octCh7
                       endif
                       ! - for octNb3 -
                       if(p0%octNb3%octN==myoctN)then
                          p0%octCh1%octNb3 => p0%octCh1 
                          p0%octCh2%octNb3 => p0%octCh2
                          p0%octCh5%octNb3 => p0%octCh5
                          p0%octCh6%octNb3 => p0%octCh6
                       else
                          p0%octCh1%octNb3 => p0%octNb3%octCh3
                          p0%octCh2%octNb3 => p0%octNb3%octCh4
                          p0%octCh5%octNb3 => p0%octNb3%octCh7
                          p0%octCh6%octNb3 => p0%octNb3%octCh8
                       endif
                       ! - for octNb4 -
                       if(p0%octNb4%octN==myoctN)then
                          p0%octCh3%octNb4 => p0%octCh3
                          p0%octCh4%octNb4 => p0%octCh4
                          p0%octCh7%octNb4 => p0%octCh7
                          p0%octCh8%octNb4 => p0%octCh8
                       else
                          p0%octCh3%octNb4 => p0%octNb4%octCh1
                          p0%octCh4%octNb4 => p0%octNb4%octCh2
                          p0%octCh7%octNb4 => p0%octNb4%octCh5
                          p0%octCh8%octNb4 => p0%octNb4%octCh6
                       endif
                       ! - for octNb5 -
                       if(p0%octNb5%octN==myoctN)then
                          p0%octCh1%octNb5 => p0%octCh1 
                          p0%octCh2%octNb5 => p0%octCh2
                          p0%octCh3%octNb5 => p0%octCh3
                          p0%octCh4%octNb5 => p0%octCh4
                       else
                          p0%octCh1%octNb5 => p0%octNb5%octCh5
                          p0%octCh2%octNb5 => p0%octNb5%octCh6
                          p0%octCh3%octNb5 => p0%octNb5%octCh7
                          p0%octCh4%octNb5 => p0%octNb5%octCh8
                       endif
                       ! - for octNb6 -
                       if(p0%octNb6%octN==myoctN)then
                          p0%octCh5%octNb6 => p0%octCh5
                          p0%octCh6%octNb6 => p0%octCh6
                          p0%octCh7%octNb6 => p0%octCh7
                          p0%octCh8%octNb6 => p0%octCh8
                       else
                          p0%octCh5%octNb6 => p0%octNb6%octCh1
                          p0%octCh6%octNb6 => p0%octNb6%octCh2
                          p0%octCh7%octNb6 => p0%octNb6%octCh3
                          p0%octCh8%octNb6 => p0%octNb6%octCh4
                       endif
                       ! - for octCh1 -
                       p0%octCh1%octNb2 => p0       %octCh2
                       p0%octCh1%octNb4 => p0       %octCh3
                       p0%octCh1%octNb6 => p0       %octCh5
                       ! - for octCh2 -
                       p0%octCh2%octNb1 => p0       %octCh1
                       p0%octCh2%octNb4 => p0       %octCh4
                       p0%octCh2%octNb6 => p0       %octCh6
                       ! - for octCh3 -
                       p0%octCh3%octNb2 => p0       %octCh4
                       p0%octCh3%octNb3 => p0       %octCh1
                       p0%octCh3%octNb6 => p0       %octCh7
                       ! - for octCh4 -
                       p0%octCh4%octNb1 => p0       %octCh3
                       p0%octCh4%octNb3 => p0       %octCh2
                       p0%octCh4%octNb6 => p0       %octCh8
                       ! - for octCh5 -
                       p0%octCh5%octNb2 => p0       %octCh6
                       p0%octCh5%octNb4 => p0       %octCh7
                       p0%octCh5%octNb5 => p0       %octCh1
                       ! - for octCh6 -
                       p0%octCh6%octNb1 => p0       %octCh5           
                       p0%octCh6%octNb4 => p0       %octCh8
                       p0%octCh6%octNb5 => p0       %octCh2
                       ! - for octCh7 -
                       p0%octCh7%octNb2 => p0       %octCh8
                       p0%octCh7%octNb3 => p0       %octCh5
                       p0%octCh7%octNb5 => p0       %octCh3
                       ! - for octCh8 -
                       p0%octCh8%octNb1 => p0       %octCh7
                       p0%octCh8%octNb3 => p0       %octCh6
                       p0%octCh8%octNb5 => p0       %octCh4

                    endif
                    if(p0%iFLG(1)==1) then !This is at the edge of overlap region
                       ! - for octCh1 -
                       p0%octCh1%octNb1 => p0%octCh1 
                       p0%octCh1%octNb2 => p0%octCh2 !
                       p0%octCh1%octNb3 => p0%octCh1 
                       p0%octCh1%octNb4 => p0%octCh3 !
                       p0%octCh1%octNb5 => p0%octCh1 
                       p0%octCh1%octNb6 => p0%octCh5 !
                       ! - for octCh2 -
                       p0%octCh2%octNb1 => p0%octCh1 !
                       p0%octCh2%octNb2 => p0%octCh2
                       p0%octCh2%octNb3 => p0%octCh2
                       p0%octCh2%octNb4 => p0%octCh4 !
                       p0%octCh2%octNb5 => p0%octCh2
                       p0%octCh2%octNb6 => p0%octCh6 !
                       ! - for octCh3 -
                       p0%octCh3%octNb1 => p0%octCh3 
                       p0%octCh3%octNb2 => p0%octCh4 !
                       p0%octCh3%octNb3 => p0%octCh1 !
                       p0%octCh3%octNb4 => p0%octCh3
                       p0%octCh3%octNb5 => p0%octCh3
                       p0%octCh3%octNb6 => p0%octCh7 !
                       ! - for octCh4 -
                       p0%octCh4%octNb1 => p0%octCh3 !
                       p0%octCh4%octNb2 => p0%octCh4
                       p0%octCh4%octNb3 => p0%octCh2 !
                       p0%octCh4%octNb4 => p0%octCh4
                       p0%octCh4%octNb5 => p0%octCh4
                       p0%octCh4%octNb6 => p0%octCh8 !
                       ! - for octCh5 -
                       p0%octCh5%octNb1 => p0%octCh5   
                       p0%octCh5%octNb2 => p0%octCh6 !
                       p0%octCh5%octNb3 => p0%octCh5
                       p0%octCh5%octNb4 => p0%octCh7 !
                       p0%octCh5%octNb5 => p0%octCh1 !
                       p0%octCh5%octNb6 => p0%octCh5
                       ! - for octCh6 -
                       p0%octCh6%octNb1 => p0%octCh5 !
                       p0%octCh6%octNb2 => p0%octCh6
                       p0%octCh6%octNb3 => p0%octCh6
                       p0%octCh6%octNb4 => p0%octCh8 !
                       p0%octCh6%octNb5 => p0%octCh2 !
                       p0%octCh6%octNb6 => p0%octCh6
                       ! - for octCh7 -
                       p0%octCh7%octNb1 => p0%octCh7 
                       p0%octCh7%octNb2 => p0%octCh8 !
                       p0%octCh7%octNb3 => p0%octCh5 !
                       p0%octCh7%octNb4 => p0%octCh7
                       p0%octCh7%octNb5 => p0%octCh3 !
                       p0%octCh7%octNb6 => p0%octCh7
                       ! - for octCh8 -
                       p0%octCh8%octNb1 => p0%octCh7 !
                       p0%octCh8%octNb2 => p0%octCh8
                       p0%octCh8%octNb3 => p0%octCh6 !
                       p0%octCh8%octNb4 => p0%octCh8
                       p0%octCh8%octNb5 => p0%octCh4 !
                       p0%octCh8%octNb6 => p0%octCh8
                       !! In case there are children in octNb1
                       if(p0%octNb1%iFLG(1) >= 1 .and. p0%octNb1%octN/=myoctN) then
                          p0%octCh1%octNb1 => p0%octNb1%octCh2
                          p0%octCh3%octNb1 => p0%octNb1%octCh4
                          p0%octCh5%octNb1 => p0%octNb1%octCh6
                          p0%octCh7%octNb1 => p0%octNb1%octCh8
                       endif
                       if(p0%octNb2%iFLG(1) >= 1 .and. p0%octNb2%octN/=myoctN) then !
                          p0%octCh2%octNb2 => p0%octNb2%octCh1
                          p0%octCh4%octNb2 => p0%octNb2%octCh3
                          p0%octCh6%octNb2 => p0%octNb2%octCh5
                          p0%octCh8%octNb2 => p0%octNb2%octCh7
                       endif
                       if(p0%octNb3%iFLG(1) >= 1 .and. p0%octNb3%octN/=myoctN) then
                          p0%octCh1%octNb3 => p0%octNb3%octCh3
                          p0%octCh2%octNb3 => p0%octNb3%octCh4
                          p0%octCh5%octNb3 => p0%octNb3%octCh7
                          p0%octCh6%octNb3 => p0%octNb3%octCh8
                       endif
                       if(p0%octNb4%iFLG(1) >= 1 .and. p0%octNb4%octN/=myoctN) then
                          p0%octCh3%octNb4 => p0%octNb4%octCh1
                          p0%octCh4%octNb4 => p0%octNb4%octCh2
                          p0%octCh7%octNb4 => p0%octNb4%octCh5
                          p0%octCh8%octNb4 => p0%octNb4%octCh6
                       endif
                       if(p0%octNb5%iFLG(1) >= 1 .and. p0%octNb5%octN/=myoctN) then
                          p0%octCh1%octNb5 => p0%octNb5%octCh5
                          p0%octCh2%octNb5 => p0%octNb5%octCh6
                          p0%octCh3%octNb5 => p0%octNb5%octCh7
                          p0%octCh4%octNb5 => p0%octNb5%octCh8
                       endif
                       if(p0%octNb6%iFLG(1) >= 1 .and. p0%octNb6%octN/=myoctN) then
                          p0%octCh5%octNb6 => p0%octNb6%octCh1
                          p0%octCh6%octNb6 => p0%octNb6%octCh2
                          p0%octCh7%octNb6 => p0%octNb6%octCh3
                          p0%octCh8%octNb6 => p0%octNb6%octCh4
                       endif
                    endif
                    !p0%iFLG(2)=1
                 endif
              end do
!!$omp end parallel do
           endif
        end do
     
        !if(debugMode>=2)print *,"end connect_Goct iLv=",iLv,rank

        return
      end subroutine connect_Goct
!***********************************************************************
      subroutine connect_oct0(iLv)
! +-------------------------------------------------------------------+
! |     create connections for additinal octs                         |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
        implicit none
        integer(kind=4)    :: iLv,i,index,index0
        type(oct), pointer :: p0,p1

        if(iLv.ge.LvMax) return 
        do i=1,2
           if(maxID(i,iLv).gt.minID(i,iLv)) then 
!!$omp parallel do private(index,index0,p0,p1)
              do index=minID(i,iLv),maxID(i,iLv)
                 p0 => Mesh(index)
                 index0 = mod(index,Nall_ini/8)
                 p1 => Mesh(index0)
                 if(p0%iFLG(2)>=1) then 
                    if(p0%iFLG(1)>=2) then 
! - for octCh1 -
                       p0%octCh1%octNb1 => p0%octNb1%octCh2 
                       p0%octCh1%octNb2 => p0       %octCh2
                       p0%octCh1%octNb3 => p0%octNb3%octCh3
                       p0%octCh1%octNb4 => p0       %octCh3
                       p0%octCh1%octNb5 => p0%octNb5%octCh5
                       p0%octCh1%octNb6 => p0       %octCh5
! - for octCh2 -
                       p0%octCh2%octNb1 => p0       %octCh1
                       p0%octCh2%octNb2 => p0%octNb2%octCh1
                       p0%octCh2%octNb3 => p0%octNb3%octCh4
                       p0%octCh2%octNb4 => p0       %octCh4
                       p0%octCh2%octNb5 => p0%octNb5%octCh6
                       p0%octCh2%octNb6 => p0       %octCh6
! - for octCh3 -
                       p0%octCh3%octNb1 => p0%octNb1%octCh4
                       p0%octCh3%octNb2 => p0       %octCh4
                       p0%octCh3%octNb3 => p0       %octCh1
                       p0%octCh3%octNb4 => p0%octNb4%octCh1
                       p0%octCh3%octNb5 => p0%octNb5%octCh7
                       p0%octCh3%octNb6 => p0       %octCh7
! - for octCh4 -
                       p0%octCh4%octNb1 => p0       %octCh3
                       p0%octCh4%octNb2 => p0%octNb2%octCh3
                       p0%octCh4%octNb3 => p0       %octCh2
                       p0%octCh4%octNb4 => p0%octNb4%octCh2
                       p0%octCh4%octNb5 => p0%octNb5%octCh8
                       p0%octCh4%octNb6 => p0       %octCh8
! - for octCh5 -
                       p0%octCh5%octNb1 => p0%octNb1%octCh6 
                       p0%octCh5%octNb2 => p0       %octCh6
                       p0%octCh5%octNb3 => p0%octNb3%octCh7
                       p0%octCh5%octNb4 => p0       %octCh7
                       p0%octCh5%octNb5 => p0       %octCh1
                       p0%octCh5%octNb6 => p0%octNb6%octCh1
! - for octCh6 -
                       p0%octCh6%octNb1 => p0       %octCh5
                       p0%octCh6%octNb2 => p0%octNb2%octCh5
                       p0%octCh6%octNb3 => p0%octNb3%octCh8
                       p0%octCh6%octNb4 => p0       %octCh8
                       p0%octCh6%octNb5 => p0       %octCh2
                       p0%octCh6%octNb6 => p0%octNb6%octCh2
! - for octCh7 -
                       p0%octCh7%octNb1 => p0%octNb1%octCh8
                       p0%octCh7%octNb2 => p0       %octCh8
                       p0%octCh7%octNb3 => p0       %octCh5
                       p0%octCh7%octNb4 => p0%octNb4%octCh5
                       p0%octCh7%octNb5 => p0       %octCh3
                       p0%octCh7%octNb6 => p0%octNb6%octCh3
! - for octCh8 -
                       p0%octCh8%octNb1 => p0       %octCh7
                       p0%octCh8%octNb2 => p0%octNb2%octCh7
                       p0%octCh8%octNb3 => p0       %octCh6
                       p0%octCh8%octNb4 => p0%octNb4%octCh6
                       p0%octCh8%octNb5 => p0       %octCh4
                       p0%octCh8%octNb6 => p0%octNb6%octCh4
                    endif
                    if(p0%iFLG(1)==1) then
! - for octCh1 -
                       p0%octCh1%octNb1 => p1%octNb1 
                       p0%octCh1%octNb2 => p0%octCh2 !
                       p0%octCh1%octNb3 => p1%octNb1 
                       p0%octCh1%octNb4 => p0%octCh3 !
                       p0%octCh1%octNb5 => p1%octNb1 
                       p0%octCh1%octNb6 => p0%octCh5 !
! - for octCh2 -
                       p0%octCh2%octNb1 => p0%octCh1 !
                       p0%octCh2%octNb2 => p1%octNb2 
                       p0%octCh2%octNb3 => p1%octNb2 
                       p0%octCh2%octNb4 => p0%octCh4 !
                       p0%octCh2%octNb5 => p1%octNb2 
                       p0%octCh2%octNb6 => p0%octCh6 !
! - for octCh3 -
                       p0%octCh3%octNb1 => p1%octNb3  
                       p0%octCh3%octNb2 => p0%octCh4 !
                       p0%octCh3%octNb3 => p0%octCh1 !
                       p0%octCh3%octNb4 => p1%octNb3 
                       p0%octCh3%octNb5 => p1%octNb3 
                       p0%octCh3%octNb6 => p0%octCh7 !
! - for octCh4 -
                       p0%octCh4%octNb1 => p0%octCh3 !
                       p0%octCh4%octNb2 => p1%octNb4 
                       p0%octCh4%octNb3 => p0%octCh2 !
                       p0%octCh4%octNb4 => p1%octNb4 
                       p0%octCh4%octNb5 => p1%octNb4 
                       p0%octCh4%octNb6 => p0%octCh8 !
! - for octCh5 -
                       p0%octCh5%octNb1 => p1%octNb5    
                       p0%octCh5%octNb2 => p0%octCh6 !
                       p0%octCh5%octNb3 => p1%octNb5 
                       p0%octCh5%octNb4 => p0%octCh7 !
                       p0%octCh5%octNb5 => p0%octCh1 !
                       p0%octCh5%octNb6 => p1%octNb5 
! - for octCh6 -
                       p0%octCh6%octNb1 => p0%octCh5 !
                       p0%octCh6%octNb2 => p1%octNb6 
                       p0%octCh6%octNb3 => p1%octNb6 
                       p0%octCh6%octNb4 => p0%octCh8 !
                       p0%octCh6%octNb5 => p0%octCh2 !
                       p0%octCh6%octNb6 => p1%octNb6 
! - for octCh7 -
                       p0%octCh7%octNb1 => p1%octNb1 
                       p0%octCh7%octNb2 => p0%octCh8 !
                       p0%octCh7%octNb3 => p0%octCh5 !
                       p0%octCh7%octNb4 => p1%octNb1 
                       p0%octCh7%octNb5 => p0%octCh3 !
                       p0%octCh7%octNb6 => p1%octNb1 
! - for octCh8 -
                       p0%octCh8%octNb1 => p0%octCh7 !
                       p0%octCh8%octNb2 => p1%octNb3 
                       p0%octCh8%octNb3 => p0%octCh6 !
                       p0%octCh8%octNb4 => p1%octNb3 
                       p0%octCh8%octNb5 => p0%octCh4 !
                       p0%octCh8%octNb6 => p1%octNb3 
!
                       if(p0%octNb1%iFLG(1) >= 1) then
                          p0%octCh1%octNb1 => p0%octNb1%octCh2
                          p0%octCh3%octNb1 => p0%octNb1%octCh4
                          p0%octCh5%octNb1 => p0%octNb1%octCh6
                          p0%octCh7%octNb1 => p0%octNb1%octCh8
                       endif
                       if(p0%octNb2%iFLG(1) >= 1) then
                          p0%octCh2%octNb2 => p0%octNb2%octCh1
                          p0%octCh4%octNb2 => p0%octNb2%octCh3
                          p0%octCh6%octNb2 => p0%octNb2%octCh5
                          p0%octCh8%octNb2 => p0%octNb2%octCh7
                       endif
                       if(p0%octNb3%iFLG(1) >= 1) then
                          p0%octCh1%octNb3 => p0%octNb3%octCh3
                          p0%octCh2%octNb3 => p0%octNb3%octCh4
                          p0%octCh5%octNb3 => p0%octNb3%octCh7
                          p0%octCh6%octNb3 => p0%octNb3%octCh8
                       endif
                       if(p0%octNb4%iFLG(1) >= 1) then
                          p0%octCh3%octNb4 => p0%octNb4%octCh1
                          p0%octCh4%octNb4 => p0%octNb4%octCh2
                          p0%octCh7%octNb4 => p0%octNb4%octCh5
                          p0%octCh8%octNb4 => p0%octNb4%octCh6
                       endif
                       if(p0%octNb5%iFLG(1) >= 1) then
                          p0%octCh1%octNb5 => p0%octNb5%octCh5
                          p0%octCh2%octNb5 => p0%octNb5%octCh6
                          p0%octCh3%octNb5 => p0%octNb5%octCh7
                          p0%octCh4%octNb5 => p0%octNb5%octCh8
                       endif
                       if(p0%octNb6%iFLG(1) >= 1) then
                          p0%octCh5%octNb6 => p0%octNb6%octCh1
                          p0%octCh6%octNb6 => p0%octNb6%octCh2
                          p0%octCh7%octNb6 => p0%octNb6%octCh3
                          p0%octCh8%octNb6 => p0%octNb6%octCh4
                       endif
                    endif
                    p0%iFLG(2)=1
                 endif
              end do
!!$omp end parallel do
           endif
        end do
     
        return
        end subroutine connect_oct0
!
!***********************************************************************
      subroutine interpolate_field(iLv,iF,iK,ch)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ; put interpolated field data to upper level  |
! +-------------------------------------------------------------------+
!***********************************************************************
        use message_passing_interface,only:rank, ierr
        use oct_set
        use param
        use const
        use init_mesh_size
        implicit none
        integer(kind=4)    :: iLv,iF,iK,is,ie
        integer(kind=4)    :: i,index
        real(kind=8)       :: b1,b2,b3,b4,c1,c2
        character*1 ch
        type(oct), pointer :: p0
        type(oct), pointer :: new1,new2,new3,new4,new5,new6,new7,new8
!        
        if(iLv==LvMax) return
        !if(maxID(1,iLv).le.minID(1,iLv)) return
        if((maxID(1,iLv)<=minID(1,iLv)).and.(maxID(3,iLv)<=minID(3,iLv))) return 
        !if(debugMode>=1)print *,"entering interpolate_field iLv=",iLv,"ch=",ch,rank
!
        call set_index(ch,is,ie)
!

        
        if(maxID(1,iLv)>minID(1,iLv))then
!!$omp parallel do private(index,p0,i,b1,b2,b3,b4,c1,c2,new1,new2,new3,new4,new5,new6,new7,new8)
           do index=minID(1,iLv),maxID(1,iLv)
              
              p0 => Mesh(index)

              if(p0%iFLG(iK)==iF) then 
              
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb1
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb1
                 new7 => p0%octNb5%octNb3
                 new8 => p0%octNb5%octNb3%octNb1

                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh1%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb2
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb2
                 new7 => p0%octNb5%octNb3
                 new8 => p0%octNb5%octNb3%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh2%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb1
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb1
                 new7 => p0%octNb5%octNb4
                 new8 => p0%octNb5%octNb4%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh3%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb2
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb2
                 new7 => p0%octNb5%octNb4
                 new8 => p0%octNb5%octNb4%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh4%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb1
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb1
                 new7 => p0%octNb6%octNb3
                 new8 => p0%octNb6%octNb3%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh5%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb2
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb2
                 new7 => p0%octNb6%octNb3
                 new8 => p0%octNb6%octNb3%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh6%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb1
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb1
                 new7 => p0%octNb6%octNb4
                 new8 => p0%octNb6%octNb4%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh7%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb2
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb2
                 new7 => p0%octNb6%octNb4
                 new8 => p0%octNb6%octNb4%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh8%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo

                 nullify(new1) ; nullify(new2)
                 nullify(new3) ; nullify(new4)
                 nullify(new5) ; nullify(new6)
                 nullify(new7) ; nullify(new8)

                 !if(index==2625 .and. rank==0)print *,"check inter",ch," end index=",index,p0%F(2)!&
                      !,p0%octCh1%F(2),p0%octCh2%F(2),p0%octCh3%F(2),p0%octCh4%F(2),&
                      !p0%octCh5%F(2),p0%octCh6%F(2),p0%octCh7%F(2),p0%octCh8%F(2),rank
              endif
           end do
!!$omp end parallel do
        endif

        

!
 
        if(maxID(3,iLv)>minID(3,iLv))then
!!$omp parallel do private(index,p0,i,b1,b2,b3,b4,c1,c2,new1,new2,new3,new4,new5,new6,new7,new8)
           do index=minID(3,iLv),maxID(3,iLv)

              p0 => Mesh(index)

              if(p0%iFLG(iK)==iF) then 
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb1
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb1
                 new7 => p0%octNb5%octNb3
                 new8 => p0%octNb5%octNb3%octNb1

                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh1%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb2
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb2
                 new7 => p0%octNb5%octNb3
                 new8 => p0%octNb5%octNb3%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh2%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb1
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb1
                 new7 => p0%octNb5%octNb4
                 new8 => p0%octNb5%octNb4%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh3%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb2
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb2
                 new7 => p0%octNb5%octNb4
                 new8 => p0%octNb5%octNb4%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh4%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb1
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb1
                 new7 => p0%octNb6%octNb3
                 new8 => p0%octNb6%octNb3%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh5%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb2
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb2
                 new7 => p0%octNb6%octNb3
                 new8 => p0%octNb6%octNb3%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh6%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb1
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb1
                 new7 => p0%octNb6%octNb4
                 new8 => p0%octNb6%octNb4%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh7%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb2
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb2
                 new7 => p0%octNb6%octNb4
                 new8 => p0%octNb6%octNb4%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh8%F(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo

                 nullify(new1) ; nullify(new2)
                 nullify(new3) ; nullify(new4)
                 nullify(new5) ; nullify(new6)
                 nullify(new7) ; nullify(new8)
              endif
           end do
!!$omp end parallel do
        endif

        !if(debugMode>=1)print *,"exiting interpolate_field iLv=",iLv,"ch=",ch,rank
        return
        end subroutine interpolate_field
        
        subroutine interpolate_fieldE(iLv,iF,iK)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ; put interpolated field data to upper level  |
! +-------------------------------------------------------------------+
!***********************************************************************
        use message_passing_interface,only:rank, ierr
        use oct_set
        use param
        use const
        use init_mesh_size
        implicit none
        integer(kind=4)    :: iLv,iF,iK
        integer(kind=4)    :: index
        real(kind=8)       :: a,b,a1,a2,b1,b2
        type(oct), pointer :: p0
!        
        if(iLv==LvMax) return
        if((maxID(1,iLv)<=minID(1,iLv)).and.(maxID(3,iLv)<=minID(3,iLv))) return 

        if(maxID(1,iLv)>minID(1,iLv))then
           do index=minID(1,iLv),maxID(1,iLv)
              p0 => Mesh(index)
              if(p0%iFLG(iK)==iF) then 
                !===  Ex  ===
                !---p0%octCh1---
                a = THR_FOURTH* p0                     %F(1) +QUARTER* p0       %octNb3       %F(1)
                b = THR_FOURTH* p0              %octNb5%F(1) +QUARTER* p0       %octNb3%octNb5%F(1)
                p0%octCh1%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a1=       HALF* p0                     %F(1) +   HALF* p0%octNb2              %F(1)
                a2=       HALF* p0       %octNb3       %F(1) +   HALF* p0%octNb2%octNb3       %F(1)
                a = THR_FOURTH* a1 + QUARTER* a2
                b1=       HALF* p0              %octNb5%F(1) +   HALF* p0%octNb2       %octNb5%F(1)
                b2=       HALF* p0       %octNb3%octNb5%F(1) +   HALF* p0%octNb2%octNb3%octNb5%F(1)
                b = THR_FOURTH* b1 + QUARTER* b2
                p0%octCh2%F(1)= THR_FOURTH* a + QUARTER* b
                !---p0%octCh3---
                a = THR_FOURTH* p0                     %F(1) +QUARTER* p0       %octNb4       %F(1)
                b = THR_FOURTH* p0              %octNb5%F(1) +QUARTER* p0       %octNb4%octNb5%F(1)
                p0%octCh3%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a1=       HALF* p0                     %F(1) +   HALF* p0%octNb2              %F(1)
                a2=       HALF* p0       %octNb4       %F(1) +   HALF* p0%octNb2%octNb4       %F(1)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1=       HALF* p0              %octNb5%F(1) +   HALF* p0%octNb2       %octNb5%F(1)
                b2=       HALF* p0       %octNb4%octNb5%F(1) +   HALF* p0%octNb2%octNb4%octNb5%F(1)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh4%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a = THR_FOURTH* p0                     %F(1) +QUARTER* p0       %octNb3       %F(1)
                b = THR_FOURTH* p0              %octNb6%F(1) +QUARTER* p0       %octNb3%octNb6%F(1)
                p0%octCh5%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a1=       HALF* p0                     %F(1) +   HALF* p0%octNb2              %F(1)
                a2=       HALF* p0       %octNb3       %F(1) +   HALF* p0%octNb2%octNb3       %F(1)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1=       HALF* p0              %octNb6%F(1) +   HALF* p0%octNb2       %octNb6%F(1)
                b2=       HALF* p0       %octNb3%octNb6%F(1) +   HALF* p0%octNb2%octNb3%octNb6%F(1)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh6%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a = THR_FOURTH* p0                     %F(1) +QUARTER* p0       %octNb4       %F(1)
                b = THR_FOURTH* p0              %octNb6%F(1) +QUARTER* p0       %octNb4%octNb6%F(1)
                p0%octCh7%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1=       HALF* p0                     %F(1) +   HALF* p0%octNb2              %F(1)
                a2=       HALF* p0       %octNb4       %F(1) +   HALF* p0%octNb2%octNb4       %F(1)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1=       HALF* p0              %octNb6%F(1) +   HALF* p0%octNb2       %octNb6%F(1)
                b2=       HALF* p0       %octNb4%octNb6%F(1) +   HALF* p0%octNb2%octNb4%octNb6%F(1)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh8%F(1)= THR_FOURTH*a + QUARTER*b
                
                !===  Ey  ===
                !---p0%octCh1---
                a = THR_FOURTH* p0                     %F(2) +QUARTER* p0              %octNb5%F(2)
                b = THR_FOURTH* p0%octNb1              %F(2) +QUARTER* p0%octNb1       %octNb5%F(2)
                p0%octCh1%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a = THR_FOURTH* p0                     %F(2) +QUARTER* p0              %octNb5%F(2)
                b = THR_FOURTH* p0%octNb2              %F(2) +QUARTER* p0%octNb2       %octNb5%F(2)
                p0%octCh2%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh3---
                a1 =      HALF* p0                     %F(2) +   HALF* p0       %octNb4       %F(2)
                a2 =      HALF* p0%octNb1              %F(2) +   HALF* p0%octNb1%octNb4       %F(2)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0              %octNb5%F(2) +   HALF* p0       %octNb4%octNb5%F(2)
                b2 =      HALF* p0%octNb1       %octNb5%F(2) +   HALF* p0%octNb1%octNb4%octNb5%F(2)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh3%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a1 =      HALF* p0                     %F(2) +   HALF* p0       %octNb4       %F(2)
                a2 =      HALF* p0%octNb2              %F(2) +   HALF* p0%octNb2%octNb4       %F(2)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0              %octNb5%F(2) +   HALF* p0       %octNb4%octNb5%F(2)
                b2 =      HALF* p0%octNb2       %octNb5%F(2) +   HALF* p0%octNb2%octNb4%octNb5%F(2)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh4%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a = THR_FOURTH* p0                     %F(2) +QUARTER* p0              %octNb6%F(2)
                b = THR_FOURTH* p0%octNb1              %F(2) +QUARTER* p0%octNb1       %octNb6%F(2)
                p0%octCh5%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a = THR_FOURTH* p0                     %F(2) +QUARTER* p0              %octNb6%F(2)
                b = THR_FOURTH* p0%octNb2              %F(2) +QUARTER* p0%octNb2       %octNb6%F(2)
                p0%octCh6%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a1 =      HALF* p0                     %F(2) +   HALF* p0       %octNb4       %F(2)
                a2 =      HALF* p0%octNb1              %F(2) +   HALF* p0%octNb1%octNb4       %F(2)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0              %octNb6%F(2) +   HALF* p0       %octNb4%octNb6%F(2)
                b2 =      HALF* p0%octNb1       %octNb6%F(2) +   HALF* p0%octNb1%octNb4%octNb6%F(2)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh7%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1 =      HALF* p0                     %F(2) +   HALF* p0       %octNb4       %F(2)
                a2 =      HALF* p0%octNb2              %F(2) +   HALF* p0%octNb2%octNb4       %F(2)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0              %octNb6%F(2) +   HALF* p0       %octNb4%octNb6%F(2)
                b2 =      HALF* p0%octNb2       %octNb6%F(2) +   HALF* p0%octNb2%octNb4%octNb6%F(2)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh8%F(2)= THR_FOURTH*a + QUARTER*b
                
                !===  Ez  ===
                !---p0%octCh1---
                a = THR_FOURTH* p0                     %F(3) +QUARTER* p0%octNb1              %F(3)
                b = THR_FOURTH* p0       %octNb3       %F(3) +QUARTER* p0%octNb1%octNb3       %F(3)
                p0%octCh1%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a = THR_FOURTH* p0                     %F(3) +QUARTER* p0%octNb2              %F(3)
                b = THR_FOURTH* p0       %octNb3       %F(3) +QUARTER* p0%octNb2%octNb3       %F(3)
                p0%octCh2%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh3---
                a = THR_FOURTH* p0                     %F(3) +QUARTER* p0%octNb1              %F(3)
                b = THR_FOURTH* p0       %octNb4       %F(3) +QUARTER* p0%octNb1%octNb4       %F(3)
                p0%octCh3%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a = THR_FOURTH* p0                     %F(3) +QUARTER* p0%octNb2              %F(3)
                b = THR_FOURTH* p0       %octNb4       %F(3) +QUARTER* p0%octNb2%octNb4       %F(3)
                p0%octCh4%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a1 =      HALF* p0                     %F(3) +   HALF* p0              %octNb6%F(3)
                a2 =      HALF* p0%octNb1              %F(3) +   HALF* p0%octNb1       %octNb6%F(3)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0       %octNb3       %F(3) +   HALF* p0       %octNb3%octNb6%F(3)
                b2 =      HALF* p0%octNb1%octNb3       %F(3) +   HALF* p0%octNb1%octNb3%octNb6%F(3)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh5%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a1 =      HALF* p0                     %F(3) +   HALF* p0              %octNb6%F(3)
                a2 =      HALF* p0%octNb2              %F(3) +   HALF* p0%octNb2       %octNb6%F(3)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0       %octNb3       %F(3) +   HALF* p0       %octNb3%octNb6%F(3)
                b2 =      HALF* p0%octNb2%octNb3       %F(3) +   HALF* p0%octNb2%octNb3%octNb6%F(3)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh6%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a1 =      HALF* p0                     %F(3) +   HALF* p0              %octNb6%F(3)
                a2 =      HALF* p0%octNb1              %F(3) +   HALF* p0%octNb1       %octNb6%F(3)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0       %octNb4       %F(3) +   HALF* p0       %octNb4%octNb6%F(3)
                b2 =      HALF* p0%octNb1%octNb4       %F(3) +   HALF* p0%octNb1%octNb4%octNb6%F(3)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh7%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1 =      HALF* p0                     %F(3) +   HALF* p0              %octNb6%F(3)
                a2 =      HALF* p0%octNb2              %F(3) +   HALF* p0%octNb2       %octNb6%F(3)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0       %octNb4       %F(3) +   HALF* p0       %octNb4%octNb6%F(3)
                b2 =      HALF* p0%octNb2%octNb4       %F(3) +   HALF* p0%octNb2%octNb4%octNb6%F(3)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh8%F(3)= THR_FOURTH*a + QUARTER*b
              endif
           end do
        endif
        
        if(maxID(3,iLv)>minID(3,iLv))then
           do index=minID(3,iLv),maxID(3,iLv)
              p0 => Mesh(index)
              if(p0%iFLG(iK)==iF) then 
                !===  Ex  ===
                !---p0%octCh1---
                a = THR_FOURTH* p0                     %F(1) +QUARTER* p0       %octNb3       %F(1)
                b = THR_FOURTH* p0              %octNb5%F(1) +QUARTER* p0       %octNb3%octNb5%F(1)
                p0%octCh1%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a1=       HALF* p0                     %F(1) +   HALF* p0%octNb2              %F(1)
                a2=       HALF* p0       %octNb3       %F(1) +   HALF* p0%octNb2%octNb3       %F(1)
                a = THR_FOURTH* a1 + QUARTER* a2
                b1=       HALF* p0              %octNb5%F(1) +   HALF* p0%octNb2       %octNb5%F(1)
                b2=       HALF* p0       %octNb3%octNb5%F(1) +   HALF* p0%octNb2%octNb3%octNb5%F(1)
                b = THR_FOURTH* b1 + QUARTER* b2
                p0%octCh2%F(1)= THR_FOURTH* a + QUARTER* b
                !---p0%octCh3---
                a = THR_FOURTH* p0                     %F(1) +QUARTER* p0       %octNb4       %F(1)
                b = THR_FOURTH* p0              %octNb5%F(1) +QUARTER* p0       %octNb4%octNb5%F(1)
                p0%octCh3%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a1=       HALF* p0                     %F(1) +   HALF* p0%octNb2              %F(1)
                a2=       HALF* p0       %octNb4       %F(1) +   HALF* p0%octNb2%octNb4       %F(1)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1=       HALF* p0              %octNb5%F(1) +   HALF* p0%octNb2       %octNb5%F(1)
                b2=       HALF* p0       %octNb4%octNb5%F(1) +   HALF* p0%octNb2%octNb4%octNb5%F(1)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh4%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a = THR_FOURTH* p0                     %F(1) +QUARTER* p0       %octNb3       %F(1)
                b = THR_FOURTH* p0              %octNb6%F(1) +QUARTER* p0       %octNb3%octNb6%F(1)
                p0%octCh5%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a1=       HALF* p0                     %F(1) +   HALF* p0%octNb2              %F(1)
                a2=       HALF* p0       %octNb3       %F(1) +   HALF* p0%octNb2%octNb3       %F(1)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1=       HALF* p0              %octNb6%F(1) +   HALF* p0%octNb2       %octNb6%F(1)
                b2=       HALF* p0       %octNb3%octNb6%F(1) +   HALF* p0%octNb2%octNb3%octNb6%F(1)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh6%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a = THR_FOURTH* p0                     %F(1) +QUARTER* p0       %octNb4       %F(1)
                b = THR_FOURTH* p0              %octNb6%F(1) +QUARTER* p0       %octNb4%octNb6%F(1)
                p0%octCh7%F(1)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1=       HALF* p0                     %F(1) +   HALF* p0%octNb2              %F(1)
                a2=       HALF* p0       %octNb4       %F(1) +   HALF* p0%octNb2%octNb4       %F(1)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1=       HALF* p0              %octNb6%F(1) +   HALF* p0%octNb2       %octNb6%F(1)
                b2=       HALF* p0       %octNb4%octNb6%F(1) +   HALF* p0%octNb2%octNb4%octNb6%F(1)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh8%F(1)= THR_FOURTH*a + QUARTER*b
                
                !===  Ey  ===
                !---p0%octCh1---
                a = THR_FOURTH* p0                     %F(2) +QUARTER* p0              %octNb5%F(2)
                b = THR_FOURTH* p0%octNb1              %F(2) +QUARTER* p0%octNb1       %octNb5%F(2)
                p0%octCh1%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a = THR_FOURTH* p0                     %F(2) +QUARTER* p0              %octNb5%F(2)
                b = THR_FOURTH* p0%octNb2              %F(2) +QUARTER* p0%octNb2       %octNb5%F(2)
                p0%octCh2%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh3---
                a1 =      HALF* p0                     %F(2) +   HALF* p0       %octNb4       %F(2)
                a2 =      HALF* p0%octNb1              %F(2) +   HALF* p0%octNb1%octNb4       %F(2)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0              %octNb5%F(2) +   HALF* p0       %octNb4%octNb5%F(2)
                b2 =      HALF* p0%octNb1       %octNb5%F(2) +   HALF* p0%octNb1%octNb4%octNb5%F(2)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh3%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a1 =      HALF* p0                     %F(2) +   HALF* p0       %octNb4       %F(2)
                a2 =      HALF* p0%octNb2              %F(2) +   HALF* p0%octNb2%octNb4       %F(2)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0              %octNb5%F(2) +   HALF* p0       %octNb4%octNb5%F(2)
                b2 =      HALF* p0%octNb2       %octNb5%F(2) +   HALF* p0%octNb2%octNb4%octNb5%F(2)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh4%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a = THR_FOURTH* p0                     %F(2) +QUARTER* p0              %octNb6%F(2)
                b = THR_FOURTH* p0%octNb1              %F(2) +QUARTER* p0%octNb1       %octNb6%F(2)
                p0%octCh5%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a = THR_FOURTH* p0                     %F(2) +QUARTER* p0              %octNb6%F(2)
                b = THR_FOURTH* p0%octNb2              %F(2) +QUARTER* p0%octNb2       %octNb6%F(2)
                p0%octCh6%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a1 =      HALF* p0                     %F(2) +   HALF* p0       %octNb4       %F(2)
                a2 =      HALF* p0%octNb1              %F(2) +   HALF* p0%octNb1%octNb4       %F(2)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0              %octNb6%F(2) +   HALF* p0       %octNb4%octNb6%F(2)
                b2 =      HALF* p0%octNb1       %octNb6%F(2) +   HALF* p0%octNb1%octNb4%octNb6%F(2)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh7%F(2)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1 =      HALF* p0                     %F(2) +   HALF* p0       %octNb4       %F(2)
                a2 =      HALF* p0%octNb2              %F(2) +   HALF* p0%octNb2%octNb4       %F(2)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0              %octNb6%F(2) +   HALF* p0       %octNb4%octNb6%F(2)
                b2 =      HALF* p0%octNb2       %octNb6%F(2) +   HALF* p0%octNb2%octNb4%octNb6%F(2)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh8%F(2)= THR_FOURTH*a + QUARTER*b
                
                !===  Ez  ===
                !---p0%octCh1---
                a = THR_FOURTH* p0                     %F(3) +QUARTER* p0%octNb1              %F(3)
                b = THR_FOURTH* p0       %octNb3       %F(3) +QUARTER* p0%octNb1%octNb3       %F(3)
                p0%octCh1%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a = THR_FOURTH* p0                     %F(3) +QUARTER* p0%octNb2              %F(3)
                b = THR_FOURTH* p0       %octNb3       %F(3) +QUARTER* p0%octNb2%octNb3       %F(3)
                p0%octCh2%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh3---
                a = THR_FOURTH* p0                     %F(3) +QUARTER* p0%octNb1              %F(3)
                b = THR_FOURTH* p0       %octNb4       %F(3) +QUARTER* p0%octNb1%octNb4       %F(3)
                p0%octCh3%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a = THR_FOURTH* p0                     %F(3) +QUARTER* p0%octNb2              %F(3)
                b = THR_FOURTH* p0       %octNb4       %F(3) +QUARTER* p0%octNb2%octNb4       %F(3)
                p0%octCh4%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a1 =      HALF* p0                     %F(3) +   HALF* p0              %octNb6%F(3)
                a2 =      HALF* p0%octNb1              %F(3) +   HALF* p0%octNb1       %octNb6%F(3)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0       %octNb3       %F(3) +   HALF* p0       %octNb3%octNb6%F(3)
                b2 =      HALF* p0%octNb1%octNb3       %F(3) +   HALF* p0%octNb1%octNb3%octNb6%F(3)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh5%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a1 =      HALF* p0                     %F(3) +   HALF* p0              %octNb6%F(3)
                a2 =      HALF* p0%octNb2              %F(3) +   HALF* p0%octNb2       %octNb6%F(3)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0       %octNb3       %F(3) +   HALF* p0       %octNb3%octNb6%F(3)
                b2 =      HALF* p0%octNb2%octNb3       %F(3) +   HALF* p0%octNb2%octNb3%octNb6%F(3)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh6%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a1 =      HALF* p0                     %F(3) +   HALF* p0              %octNb6%F(3)
                a2 =      HALF* p0%octNb1              %F(3) +   HALF* p0%octNb1       %octNb6%F(3)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0       %octNb4       %F(3) +   HALF* p0       %octNb4%octNb6%F(3)
                b2 =      HALF* p0%octNb1%octNb4       %F(3) +   HALF* p0%octNb1%octNb4%octNb6%F(3)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh7%F(3)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1 =      HALF* p0                     %F(3) +   HALF* p0              %octNb6%F(3)
                a2 =      HALF* p0%octNb2              %F(3) +   HALF* p0%octNb2       %octNb6%F(3)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 =      HALF* p0       %octNb4       %F(3) +   HALF* p0       %octNb4%octNb6%F(3)
                b2 =      HALF* p0%octNb2%octNb4       %F(3) +   HALF* p0%octNb2%octNb4%octNb6%F(3)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh8%F(3)= THR_FOURTH*a + QUARTER*b
              endif
           end do
        endif

        return
        end subroutine interpolate_fieldE
        
        subroutine interpolate_fieldB(iLv,iF,iK)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ; put interpolated field data to upper level  |
! +-------------------------------------------------------------------+
!***********************************************************************
        use message_passing_interface,only:rank, ierr
        use oct_set
        use param
        use const
        use init_mesh_size
        implicit none
        integer(kind=4)    :: iLv,iF,iK
        integer(kind=4)    :: index
        real(kind=8)       :: a,b,a1,a2,b1,b2
        type(oct), pointer :: p0
!        
        if(iLv==LvMax) return
        if((maxID(1,iLv)<=minID(1,iLv)).and.(maxID(3,iLv)<=minID(3,iLv))) return 

        if(maxID(1,iLv)>minID(1,iLv))then
           do index=minID(1,iLv),maxID(1,iLv)
              p0 => Mesh(index)
              if(p0%iFLG(iK)==iF) then 
                !===  Bx  ===
                !---p0%octCh1---
                p0%octCh1%F(4)=THR_FOURTH*p0%F(4)+QUARTER*p0%octNb1%F(4)
                !---p0%octCh2---
                p0%octCh2%F(4)=THR_FOURTH*p0%F(4)+QUARTER*p0%octNb2%F(4)
                !---p0%octCh3---
                a  = HALF*p0       %F(4)+HALF*p0       %octNb4%F(4)
                b  = HALF*p0%octNb1%F(4)+HALF*p0%octNb1%octNb4%F(4)
                p0%octCh3%F(4)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh4---
                a  = HALF*p0       %F(4)+HALF*p0       %octNb4%F(4)
                b  = HALF*p0%octNb2%F(4)+HALF*p0%octNb2%octNb4%F(4)
                p0%octCh4%F(4)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh5---
                a  = HALF*p0       %F(4)+HALF*p0       %octNb6%F(4)
                b  = HALF*p0%octNb1%F(4)+HALF*p0%octNb1%octNb6%F(4)
                p0%octCh5%F(4)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh6---
                a  = HALF*p0       %F(4)+HALF*p0       %octNb6%F(4)
                b  = HALF*p0%octNb2%F(4)+HALF*p0%octNb2%octNb6%F(4)
                p0%octCh6%F(4)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh7---
                a1 = HALF*p0       %F(4)+HALF*p0       %octNb6%F(4)
                a2 = HALF*p0%octNb4%F(4)+HALF*p0%octNb4%octNb6%F(4)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0%octNb1       %F(4)+HALF*p0%octNb1       %octNb6%F(4)
                b2 = HALF*p0%octNb1%octNb4%F(4)+HALF*p0%octNb1%octNb4%octNb6%F(4)
                b  = HALF*b1+HALF*b2
                p0%octCh7%F(4)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh8---
                a1 = HALF*p0       %F(4)+HALF*p0       %octNb6%F(4)
                a2 = HALF*p0%octNb4%F(4)+HALF*p0%octNb4%octNb6%F(4)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0%octNb2       %F(4)+HALF*p0%octNb2       %octNb6%F(4)
                b2 = HALF*p0%octNb2%octNb4%F(4)+HALF*p0%octNb2%octNb4%octNb6%F(4)
                b  = HALF*b1+HALF*b2
                p0%octCh8%F(4)=THR_FOURTH*a+QUARTER*b
                
                !===  By  ===
                !---p0%octCh1---
                p0%octCh1%F(5)=THR_FOURTH*p0%F(5)+QUARTER*p0%octNb3%F(5)
                !---p0%octCh2---
                a  = HALF*p0       %F(5)+HALF*p0       %octNb2%F(5)
                b  = HALF*p0%octNb3%F(5)+HALF*p0%octNb3%octNb2%F(5)
                p0%octCh2%F(5)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh3---
                p0%octCh3%F(5)=THR_FOURTH*p0%F(5)+QUARTER*p0%octNb4%F(5)
                !---p0%octCh4---
                a  = HALF*p0       %F(5)+HALF*p0       %octNb2%F(5)
                b  = HALF*p0%octNb4%F(5)+HALF*p0%octNb4%octNb2%F(5)
                p0%octCh4%F(5)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh5---
                a  = HALF*p0       %F(5)+HALF*p0       %octNb6%F(5)
                b  = HALF*p0%octNb3%F(5)+HALF*p0%octNb3%octNb6%F(5)
                p0%octCh5%F(5)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh6---
                a1 = HALF*p0       %F(5)+HALF*p0       %octNb2%F(5)
                a2 = HALF*p0%octNb6%F(5)+HALF*p0%octNb6%octNb2%F(5)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0       %octNb3%F(5)+HALF*p0       %octNb3%octNb2%F(5)
                b2 = HALF*p0%octNb6%octNb3%F(5)+HALF*p0%octNb6%octNb3%octNb2%F(5)
                b  = HALF*b1+HALF*b2
                p0%octCh6%F(5)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh7---
                a  = HALF*p0       %F(5)+HALF*p0       %octNb6%F(5)
                b  = HALF*p0%octNb4%F(5)+HALF*p0%octNb4%octNb6%F(5)
                p0%octCh7%F(5)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh8---
                a1 = HALF*p0       %F(5)+HALF*p0       %octNb2%F(5)
                a2 = HALF*p0%octNb6%F(5)+HALF*p0%octNb6%octNb2%F(5)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0       %octNb4%F(5)+HALF*p0       %octNb4%octNb2%F(5)
                b2 = HALF*p0%octNb6%octNb4%F(5)+HALF*p0%octNb6%octNb4%octNb2%F(5)
                b  = HALF*b1+HALF*b2
                p0%octCh8%F(5)=THR_FOURTH*a+QUARTER*b
                
                !===  Bz  ===
                !---p0%octCh1---
                p0%octCh1%F(6)= THR_FOURTH*p0%F(6)+QUARTER*p0%octNb5%F(6)
                !---p0%octCh2---
                a = HALF*p0       %F(6)+HALF*p0       %octNb2%F(6)
                b = HALF*p0%octNb5%F(6)+HALF*p0%octNb5%octNb2%F(6)
                p0%octCh2%F(6)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh3---
                a = HALF*p0       %F(6)+HALF*p0       %octNb4%F(6)
                b = HALF*p0%octNb5%F(6)+HALF*p0%octNb5%octNb4%F(6)
                p0%octCh3%F(6)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh4---
                a1 = HALF*p0              %F(6)+HALF*p0              %octNb2%F(6)
                a2 = HALF*p0       %octNb4%F(6)+HALF*p0       %octNb4%octNb2%F(6)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0%octNb5       %F(6)+HALF*p0%octNb5       %octNb2%F(6)
                b2 = HALF*p0%octNb5%octNb4%F(6)+HALF*p0%octNb5%octNb4%octNb2%F(6)
                b  = HALF*b1+HALF*b2
                p0%octCh4%F(6)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh5---
                p0%octCh5%F(6)= THR_FOURTH*p0%F(6)+QUARTER*p0%octNb6%F(6)
                !---p0%octCh6---
                a = HALF*p0       %F(6)+HALF*p0       %octNb2%F(6)
                b = HALF*p0%octNb6%F(6)+HALF*p0%octNb6%octNb2%F(6)
                p0%octCh6%F(6)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh7---
                a = HALF*p0       %F(6)+HALF*p0       %octNb4%F(6)
                b = HALF*p0%octNb6%F(6)+HALF*p0%octNb6%octNb4%F(6)
                p0%octCh7%F(6)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh8---
                a1 = HALF*p0              %F(6)+HALF*p0              %octNb2%F(6)
                a2 = HALF*p0       %octNb4%F(6)+HALF*p0       %octNb4%octNb2%F(6)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0%octNb6       %F(6)+HALF*p0%octNb6       %octNb2%F(6)
                b2 = HALF*p0%octNb6%octNb4%F(6)+HALF*p0%octNb6%octNb4%octNb2%F(6)
                b  = HALF*b1+HALF*b2
                p0%octCh8%F(6)=THR_FOURTH*a+QUARTER*b
                
              endif
           end do
        endif
        
        if(maxID(3,iLv)>minID(3,iLv))then
           do index=minID(3,iLv),maxID(3,iLv)
              p0 => Mesh(index)
              if(p0%iFLG(iK)==iF) then 
                !===  Bx  ===
                !---p0%octCh1---
                p0%octCh1%F(4)=THR_FOURTH*p0%F(4)+QUARTER*p0%octNb1%F(4)
                !---p0%octCh2---
                p0%octCh2%F(4)=THR_FOURTH*p0%F(4)+QUARTER*p0%octNb2%F(4)
                !---p0%octCh3---
                a  = HALF*p0       %F(4)+HALF*p0       %octNb4%F(4)
                b  = HALF*p0%octNb1%F(4)+HALF*p0%octNb1%octNb4%F(4)
                p0%octCh3%F(4)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh4---
                a  = HALF*p0       %F(4)+HALF*p0       %octNb4%F(4)
                b  = HALF*p0%octNb2%F(4)+HALF*p0%octNb2%octNb4%F(4)
                p0%octCh4%F(4)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh5---
                a  = HALF*p0       %F(4)+HALF*p0       %octNb6%F(4)
                b  = HALF*p0%octNb1%F(4)+HALF*p0%octNb1%octNb6%F(4)
                p0%octCh5%F(4)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh6---
                a  = HALF*p0       %F(4)+HALF*p0       %octNb6%F(4)
                b  = HALF*p0%octNb2%F(4)+HALF*p0%octNb2%octNb6%F(4)
                p0%octCh6%F(4)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh7---
                a1 = HALF*p0       %F(4)+HALF*p0       %octNb6%F(4)
                a2 = HALF*p0%octNb4%F(4)+HALF*p0%octNb4%octNb6%F(4)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0%octNb1       %F(4)+HALF*p0%octNb1       %octNb6%F(4)
                b2 = HALF*p0%octNb1%octNb4%F(4)+HALF*p0%octNb1%octNb4%octNb6%F(4)
                b  = HALF*b1+HALF*b2
                p0%octCh7%F(4)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh8---
                a1 = HALF*p0       %F(4)+HALF*p0       %octNb6%F(4)
                a2 = HALF*p0%octNb4%F(4)+HALF*p0%octNb4%octNb6%F(4)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0%octNb2       %F(4)+HALF*p0%octNb2       %octNb6%F(4)
                b2 = HALF*p0%octNb2%octNb4%F(4)+HALF*p0%octNb2%octNb4%octNb6%F(4)
                b  = HALF*b1+HALF*b2
                p0%octCh8%F(4)=THR_FOURTH*a+QUARTER*b
                
                !===  By  ===
                !---p0%octCh1---
                p0%octCh1%F(5)=THR_FOURTH*p0%F(5)+QUARTER*p0%octNb3%F(5)
                !---p0%octCh2---
                a  = HALF*p0       %F(5)+HALF*p0       %octNb2%F(5)
                b  = HALF*p0%octNb3%F(5)+HALF*p0%octNb3%octNb2%F(5)
                p0%octCh2%F(5)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh3---
                p0%octCh3%F(5)=THR_FOURTH*p0%F(5)+QUARTER*p0%octNb4%F(5)
                !---p0%octCh4---
                a  = HALF*p0       %F(5)+HALF*p0       %octNb2%F(5)
                b  = HALF*p0%octNb4%F(5)+HALF*p0%octNb4%octNb2%F(5)
                p0%octCh4%F(5)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh5---
                a  = HALF*p0       %F(5)+HALF*p0       %octNb6%F(5)
                b  = HALF*p0%octNb3%F(5)+HALF*p0%octNb3%octNb6%F(5)
                p0%octCh5%F(5)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh6---
                a1 = HALF*p0       %F(5)+HALF*p0       %octNb2%F(5)
                a2 = HALF*p0%octNb6%F(5)+HALF*p0%octNb6%octNb2%F(5)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0       %octNb3%F(5)+HALF*p0       %octNb3%octNb2%F(5)
                b2 = HALF*p0%octNb6%octNb3%F(5)+HALF*p0%octNb6%octNb3%octNb2%F(5)
                b  = HALF*b1+HALF*b2
                p0%octCh6%F(5)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh7---
                a  = HALF*p0       %F(5)+HALF*p0       %octNb6%F(5)
                b  = HALF*p0%octNb4%F(5)+HALF*p0%octNb4%octNb6%F(5)
                p0%octCh7%F(5)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh8---
                a1 = HALF*p0       %F(5)+HALF*p0       %octNb2%F(5)
                a2 = HALF*p0%octNb6%F(5)+HALF*p0%octNb6%octNb2%F(5)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0       %octNb4%F(5)+HALF*p0       %octNb4%octNb2%F(5)
                b2 = HALF*p0%octNb6%octNb4%F(5)+HALF*p0%octNb6%octNb4%octNb2%F(5)
                b  = HALF*b1+HALF*b2
                p0%octCh8%F(5)=THR_FOURTH*a+QUARTER*b
                
                !===  Bz  ===
                !---p0%octCh1---
                p0%octCh1%F(6)= THR_FOURTH*p0%F(6)+QUARTER*p0%octNb5%F(6)
                !---p0%octCh2---
                a = HALF*p0       %F(6)+HALF*p0       %octNb2%F(6)
                b = HALF*p0%octNb5%F(6)+HALF*p0%octNb5%octNb2%F(6)
                p0%octCh2%F(6)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh3---
                a = HALF*p0       %F(6)+HALF*p0       %octNb4%F(6)
                b = HALF*p0%octNb5%F(6)+HALF*p0%octNb5%octNb4%F(6)
                p0%octCh3%F(6)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh4---
                a1 = HALF*p0              %F(6)+HALF*p0              %octNb2%F(6)
                a2 = HALF*p0       %octNb4%F(6)+HALF*p0       %octNb4%octNb2%F(6)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0%octNb5       %F(6)+HALF*p0%octNb5       %octNb2%F(6)
                b2 = HALF*p0%octNb5%octNb4%F(6)+HALF*p0%octNb5%octNb4%octNb2%F(6)
                b  = HALF*b1+HALF*b2
                p0%octCh4%F(6)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh5---
                p0%octCh5%F(6)= THR_FOURTH*p0%F(6)+QUARTER*p0%octNb6%F(6)
                !---p0%octCh6---
                a = HALF*p0       %F(6)+HALF*p0       %octNb2%F(6)
                b = HALF*p0%octNb6%F(6)+HALF*p0%octNb6%octNb2%F(6)
                p0%octCh6%F(6)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh7---
                a = HALF*p0       %F(6)+HALF*p0       %octNb4%F(6)
                b = HALF*p0%octNb6%F(6)+HALF*p0%octNb6%octNb4%F(6)
                p0%octCh7%F(6)=THR_FOURTH*a+QUARTER*b
                !---p0%octCh8---
                a1 = HALF*p0              %F(6)+HALF*p0              %octNb2%F(6)
                a2 = HALF*p0       %octNb4%F(6)+HALF*p0       %octNb4%octNb2%F(6)
                a  = HALF*a1+HALF*a2
                b1 = HALF*p0%octNb6       %F(6)+HALF*p0%octNb6       %octNb2%F(6)
                b2 = HALF*p0%octNb6%octNb4%F(6)+HALF*p0%octNb6%octNb4%octNb2%F(6)
                b  = HALF*b1+HALF*b2
                p0%octCh8%F(6)=THR_FOURTH*a+QUARTER*b
              endif
           end do
        endif

        return
        end subroutine interpolate_fieldB
        
        subroutine interpolate_fieldJ(iLv,iF,iK,icon)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ; put interpolated field data to upper level  |
! +-------------------------------------------------------------------+
!***********************************************************************
        use message_passing_interface,only:rank, ierr
        use oct_set
        use param
        use const
        use init_mesh_size
        implicit none
        integer(kind=4)    :: iLv,iF,iK,is,ie
        integer(kind=4)    :: i,index,icon
        real(kind=8)       :: a,b,a1,a2,b1,b2
        real(kind=8)       :: t1,t2,t3,t4
        type(oct), pointer :: p0
!        
        if(iLv.ge.LvMax) return
        if(iLv.lt.0) return 
        if((maxID(1,iLv)<=minID(1,iLv)).and.(maxID(3,iLv)<=minID(3,iLv))) return 
        is=7 ; ie=12
        if(icon.eq.0) then 
           t1 = -0.5d0 ; t2 = 1.5d0 ; t3 = -1.d0 ; t4 = 2.d0
        elseif(icon.eq.1) then 
           t1 =  1.5d0 ; t2 =-0.5d0 ; t3 =  1.d0 ; t4 = 0.d0
        elseif(icon.eq.2) then 
           t1 =  0.5d0 ; t2 = 0.5d0 ; t3 =  0.d0 ; t4 = 1.d0
        else
           t1 =  0.0d0 ; t2 = 0.0d0 ; t3 =  0.d0 ; t4 = 0.d0
        end if
        if(maxID(1,iLv)>minID(1,iLv))then
           do index=minID(1,iLv),maxID(1,iLv)
              p0 => Mesh(index)
              if(p0%iFLG(iK)==iF) then 
                do i=1,2
                !===  Jx  ===
                !---p0%octCh1---
                a = THR_FOURTH* p0       %F(3*i+4) +QUARTER* p0       %octNb3%F(3*i+4)
                b = THR_FOURTH* p0%octNb5%F(3*i+4) +QUARTER* p0%octNb5%octNb3%F(3*i+4)
                p0%octCh1%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a1= HALF* p0       %F(3*i+4) +HALF* p0       %octNb2%F(3*i+4)
                a2= HALF* p0%octNb3%F(3*i+4) +HALF* p0%octNb3%octNb2%F(3*i+4)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1= HALF* p0%octNb5       %F(3*i+4) +HALF* p0%octNb5       %octNb2%F(3*i+4)
                b2= HALF* p0%octNb5%octNb3%F(3*i+4) +HALF* p0%octNb5%octNb3%octNb2%F(3*i+4)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh2%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh3---
                a = THR_FOURTH* p0       %F(3*i+4) +QUARTER* p0       %octNb4%F(3*i+4)
                b = THR_FOURTH* p0%octNb5%F(3*i+4) +QUARTER* p0%octNb5%octNb4%F(3*i+4)
                p0%octCh3%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a1= HALF* p0       %F(3*i+4) +HALF* p0       %octNb2%F(3*i+4)
                a2= HALF* p0%octNb4%F(3*i+4) +HALF* p0%octNb4%octNb2%F(3*i+4)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1= HALF* p0%octNb5       %F(3*i+4) +HALF* p0%octNb5       %octNb2%F(3*i+4)
                b2= HALF* p0%octNb5%octNb4%F(3*i+4) +HALF* p0%octNb5%octNb4%octNb2%F(3*i+4)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh4%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a = THR_FOURTH* p0       %F(3*i+4) +QUARTER* p0       %octNb3%F(3*i+4)
                b = THR_FOURTH* p0%octNb6%F(3*i+4) +QUARTER* p0%octNb6%octNb3%F(3*i+4)
                p0%octCh5%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a1= HALF* p0       %F(3*i+4) +HALF* p0       %octNb2%F(3*i+4)
                a2= HALF* p0%octNb3%F(3*i+4) +HALF* p0%octNb3%octNb2%F(3*i+4)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1= HALF* p0%octNb6       %F(3*i+4) +HALF* p0%octNb6       %octNb2%F(3*i+4)
                b2= HALF* p0%octNb6%octNb3%F(3*i+4) +HALF* p0%octNb6%octNb3%octNb2%F(3*i+4)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh6%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a = THR_FOURTH* p0       %F(3*i+4) +QUARTER* p0       %octNb4%F(3*i+4)
                b = THR_FOURTH* p0%octNb6%F(3*i+4) +QUARTER* p0%octNb6%octNb4%F(3*i+4)
                p0%octCh7%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1= HALF* p0       %F(3*i+4) +HALF* p0       %octNb2%F(3*i+4)
                a2= HALF* p0%octNb4%F(3*i+4) +HALF* p0%octNb4%octNb2%F(3*i+4)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1= HALF* p0%octNb6       %F(3*i+4) +HALF* p0%octNb6       %octNb2%F(3*i+4)
                b2= HALF* p0%octNb6%octNb4%F(3*i+4) +HALF* p0%octNb6%octNb4%octNb2%F(3*i+4)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh8%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                
                !===  Jy  ===
                !---p0%octCh1---
                a = THR_FOURTH* p0       %F(3*i+5) +QUARTER* p0       %octNb5%F(3*i+5)
                b = THR_FOURTH* p0%octNb1%F(3*i+5) +QUARTER* p0%octNb1%octNb5%F(3*i+5)
                p0%octCh1%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a = THR_FOURTH* p0       %F(3*i+5) +QUARTER* p0       %octNb5%F(3*i+5)
                b = THR_FOURTH* p0%octNb2%F(3*i+5) +QUARTER* p0%octNb2%octNb5%F(3*i+5)
                p0%octCh2%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh3---
                a1 = HALF*p0       %F(3*i+5)+HALF*p0       %octNb4%F(3*i+5)
                a2 = HALF*p0%octNb1%F(3*i+5)+HALF*p0%octNb1%octNb4%F(3*i+5)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb5       %F(3*i+5)+HALF*p0%octNb5       %octNb4%F(3*i+5)
                b2 = HALF*p0%octNb5%octNb1%F(3*i+5)+HALF*p0%octNb5%octNb1%octNb4%F(3*i+5)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh3%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a1 = HALF*p0       %F(3*i+5)+HALF*p0       %octNb4%F(3*i+5)
                a2 = HALF*p0%octNb2%F(3*i+5)+HALF*p0%octNb2%octNb4%F(3*i+5)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb5       %F(3*i+5)+HALF*p0%octNb5       %octNb4%F(3*i+5)
                b2 = HALF*p0%octNb5%octNb2%F(3*i+5)+HALF*p0%octNb5%octNb2%octNb4%F(3*i+5)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh4%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a = THR_FOURTH* p0       %F(3*i+5) +QUARTER* p0       %octNb6%F(3*i+5)
                b = THR_FOURTH* p0%octNb1%F(3*i+5) +QUARTER* p0%octNb1%octNb6%F(3*i+5)
                p0%octCh5%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a = THR_FOURTH* p0       %F(3*i+5) +QUARTER* p0       %octNb6%F(3*i+5)
                b = THR_FOURTH* p0%octNb2%F(3*i+5) +QUARTER* p0%octNb2%octNb6%F(3*i+5)
                p0%octCh6%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a1 = HALF*p0       %F(3*i+5)+HALF*p0       %octNb4%F(3*i+5)
                a2 = HALF*p0%octNb1%F(3*i+5)+HALF*p0%octNb1%octNb4%F(3*i+5)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb6       %F(3*i+5)+HALF*p0%octNb6       %octNb4%F(3*i+5)
                b2 = HALF*p0%octNb6%octNb1%F(3*i+5)+HALF*p0%octNb6%octNb1%octNb4%F(3*i+5)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh7%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1 = HALF*p0       %F(3*i+5)+HALF*p0       %octNb4%F(3*i+5)
                a2 = HALF*p0%octNb2%F(3*i+5)+HALF*p0%octNb2%octNb4%F(3*i+5)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb6       %F(3*i+5)+HALF*p0%octNb6       %octNb4%F(3*i+5)
                b2 = HALF*p0%octNb6%octNb2%F(3*i+5)+HALF*p0%octNb6%octNb2%octNb4%F(3*i+5)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh8%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                
                !===  Jz  ===
                !---p0%octCh1---
                a = THR_FOURTH*p0       %F(3*i+6)+QUARTER*p0       %octNb1%F(3*i+6)
                b = THR_FOURTH*p0%octNb3%F(3*i+6)+QUARTER*p0%octNb3%octNb1%F(3*i+6)
                p0%octCh1%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a = THR_FOURTH*p0       %F(3*i+6)+QUARTER*p0       %octNb2%F(3*i+6)
                b = THR_FOURTH*p0%octNb3%F(3*i+6)+QUARTER*p0%octNb3%octNb2%F(3*i+6)
                p0%octCh2%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh3---
                a = THR_FOURTH*p0       %F(3*i+6)+QUARTER*p0       %octNb1%F(3*i+6)
                b = THR_FOURTH*p0%octNb4%F(3*i+6)+QUARTER*p0%octNb4%octNb1%F(3*i+6)
                p0%octCh3%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a = THR_FOURTH*p0       %F(3*i+6)+QUARTER*p0       %octNb2%F(3*i+6)
                b = THR_FOURTH*p0%octNb4%F(3*i+6)+QUARTER*p0%octNb4%octNb2%F(3*i+6)
                p0%octCh4%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a1 = HALF*p0              %F(3*i+6)+HALF*p0              %octNb6%F(3*i+6)
                a2 = HALF*p0       %octNb1%F(3*i+6)+HALF*p0       %octNb1%octNb6%F(3*i+6)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb3       %F(3*i+6)+HALF*p0%octNb3       %octNb6%F(3*i+6)
                b2 = HALF*p0%octNb3%octNb1%F(3*i+6)+HALF*p0%octNb3%octNb1%octNb6%F(3*i+6)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh5%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a1 = HALF*p0              %F(3*i+6)+HALF*p0              %octNb6%F(3*i+6)
                a2 = HALF*p0       %octNb2%F(3*i+6)+HALF*p0       %octNb2%octNb6%F(3*i+6)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb3       %F(3*i+6)+HALF*p0%octNb3       %octNb6%F(3*i+6)
                b2 = HALF*p0%octNb3%octNb2%F(3*i+6)+HALF*p0%octNb3%octNb2%octNb6%F(3*i+6)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh6%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a1 = HALF*p0              %F(3*i+6)+HALF*p0              %octNb6%F(3*i+6)
                a2 = HALF*p0       %octNb1%F(3*i+6)+HALF*p0       %octNb1%octNb6%F(3*i+6)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb4       %F(3*i+6)+HALF*p0%octNb4       %octNb6%F(3*i+6)
                b2 = HALF*p0%octNb4%octNb1%F(3*i+6)+HALF*p0%octNb4%octNb1%octNb6%F(3*i+6)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh7%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1 = HALF*p0              %F(3*i+6)+HALF*p0              %octNb6%F(3*i+6)
                a2 = HALF*p0       %octNb2%F(3*i+6)+HALF*p0       %octNb2%octNb6%F(3*i+6)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb4       %F(3*i+6)+HALF*p0%octNb4       %octNb6%F(3*i+6)
                b2 = HALF*p0%octNb4%octNb2%F(3*i+6)+HALF*p0%octNb4%octNb2%octNb6%F(3*i+6)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh8%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                end do
                do i=7,9
                    p0%octCh1%F(  i) = (t1*p0%octCh1%C(i)+t2*p0%octCh1%C(i+3))
                    p0%octCh2%F(  i) = (t1*p0%octCh2%C(i)+t2*p0%octCh2%C(i+3))
                    p0%octCh3%F(  i) = (t1*p0%octCh3%C(i)+t2*p0%octCh3%C(i+3))
                    p0%octCh4%F(  i) = (t1*p0%octCh4%C(i)+t2*p0%octCh4%C(i+3))
                    p0%octCh5%F(  i) = (t1*p0%octCh5%C(i)+t2*p0%octCh5%C(i+3))
                    p0%octCh6%F(  i) = (t1*p0%octCh6%C(i)+t2*p0%octCh6%C(i+3))
                    p0%octCh7%F(  i) = (t1*p0%octCh7%C(i)+t2*p0%octCh7%C(i+3))
                    p0%octCh8%F(  i) = (t1*p0%octCh8%C(i)+t2*p0%octCh8%C(i+3))

                    p0%octCh1%F(i+3) = (t3*p0%octCh1%C(i)+t4*p0%octCh1%C(i+3))
                    p0%octCh2%F(i+3) = (t3*p0%octCh2%C(i)+t4*p0%octCh2%C(i+3))
                    p0%octCh3%F(i+3) = (t3*p0%octCh3%C(i)+t4*p0%octCh3%C(i+3))
                    p0%octCh4%F(i+3) = (t3*p0%octCh4%C(i)+t4*p0%octCh4%C(i+3))
                    p0%octCh5%F(i+3) = (t3*p0%octCh5%C(i)+t4*p0%octCh5%C(i+3))
                    p0%octCh6%F(i+3) = (t3*p0%octCh6%C(i)+t4*p0%octCh6%C(i+3))
                    p0%octCh7%F(i+3) = (t3*p0%octCh7%C(i)+t4*p0%octCh7%C(i+3))
                    p0%octCh8%F(i+3) = (t3*p0%octCh8%C(i)+t4*p0%octCh8%C(i+3))
                end do
              endif
           end do
        endif
        
        if(maxID(3,iLv)>minID(3,iLv))then
           do index=minID(3,iLv),maxID(3,iLv)
              p0 => Mesh(index)
              if(p0%iFLG(iK)==iF) then 
                do i=1,2
                !===  Jx  ===
                !---p0%octCh1---
                a = THR_FOURTH* p0       %F(3*i+4) +QUARTER* p0       %octNb3%F(3*i+4)
                b = THR_FOURTH* p0%octNb5%F(3*i+4) +QUARTER* p0%octNb5%octNb3%F(3*i+4)
                p0%octCh1%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a1= HALF* p0       %F(3*i+4) +HALF* p0       %octNb2%F(3*i+4)
                a2= HALF* p0%octNb3%F(3*i+4) +HALF* p0%octNb3%octNb2%F(3*i+4)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1= HALF* p0%octNb5       %F(3*i+4) +HALF* p0%octNb5       %octNb2%F(3*i+4)
                b2= HALF* p0%octNb5%octNb3%F(3*i+4) +HALF* p0%octNb5%octNb3%octNb2%F(3*i+4)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh2%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh3---
                a = THR_FOURTH* p0       %F(3*i+4) +QUARTER* p0       %octNb4%F(3*i+4)
                b = THR_FOURTH* p0%octNb5%F(3*i+4) +QUARTER* p0%octNb5%octNb4%F(3*i+4)
                p0%octCh3%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a1= HALF* p0       %F(3*i+4) +HALF* p0       %octNb2%F(3*i+4)
                a2= HALF* p0%octNb4%F(3*i+4) +HALF* p0%octNb4%octNb2%F(3*i+4)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1= HALF* p0%octNb5       %F(3*i+4) +HALF* p0%octNb5       %octNb2%F(3*i+4)
                b2= HALF* p0%octNb5%octNb4%F(3*i+4) +HALF* p0%octNb5%octNb4%octNb2%F(3*i+4)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh4%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a = THR_FOURTH* p0       %F(3*i+4) +QUARTER* p0       %octNb3%F(3*i+4)
                b = THR_FOURTH* p0%octNb6%F(3*i+4) +QUARTER* p0%octNb6%octNb3%F(3*i+4)
                p0%octCh5%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a1= HALF* p0       %F(3*i+4) +HALF* p0       %octNb2%F(3*i+4)
                a2= HALF* p0%octNb3%F(3*i+4) +HALF* p0%octNb3%octNb2%F(3*i+4)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1= HALF* p0%octNb6       %F(3*i+4) +HALF* p0%octNb6       %octNb2%F(3*i+4)
                b2= HALF* p0%octNb6%octNb3%F(3*i+4) +HALF* p0%octNb6%octNb3%octNb2%F(3*i+4)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh6%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a = THR_FOURTH* p0       %F(3*i+4) +QUARTER* p0       %octNb4%F(3*i+4)
                b = THR_FOURTH* p0%octNb6%F(3*i+4) +QUARTER* p0%octNb6%octNb4%F(3*i+4)
                p0%octCh7%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1= HALF* p0       %F(3*i+4) +HALF* p0       %octNb2%F(3*i+4)
                a2= HALF* p0%octNb4%F(3*i+4) +HALF* p0%octNb4%octNb2%F(3*i+4)
                a = THR_FOURTH*a1 + QUARTER*a2
                b1= HALF* p0%octNb6       %F(3*i+4) +HALF* p0%octNb6       %octNb2%F(3*i+4)
                b2= HALF* p0%octNb6%octNb4%F(3*i+4) +HALF* p0%octNb6%octNb4%octNb2%F(3*i+4)
                b = THR_FOURTH*b1 + QUARTER*b2
                p0%octCh8%C(3*i+4)= THR_FOURTH*a + QUARTER*b
                
                !===  Jy  ===
                !---p0%octCh1---
                a = THR_FOURTH* p0       %F(3*i+5) +QUARTER* p0       %octNb5%F(3*i+5)
                b = THR_FOURTH* p0%octNb1%F(3*i+5) +QUARTER* p0%octNb1%octNb5%F(3*i+5)
                p0%octCh1%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a = THR_FOURTH* p0       %F(3*i+5) +QUARTER* p0       %octNb5%F(3*i+5)
                b = THR_FOURTH* p0%octNb2%F(3*i+5) +QUARTER* p0%octNb2%octNb5%F(3*i+5)
                p0%octCh2%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh3---
                a1 = HALF*p0       %F(3*i+5)+HALF*p0       %octNb4%F(3*i+5)
                a2 = HALF*p0%octNb1%F(3*i+5)+HALF*p0%octNb1%octNb4%F(3*i+5)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb5       %F(3*i+5)+HALF*p0%octNb5       %octNb4%F(3*i+5)
                b2 = HALF*p0%octNb5%octNb1%F(3*i+5)+HALF*p0%octNb5%octNb1%octNb4%F(3*i+5)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh3%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a1 = HALF*p0       %F(3*i+5)+HALF*p0       %octNb4%F(3*i+5)
                a2 = HALF*p0%octNb2%F(3*i+5)+HALF*p0%octNb2%octNb4%F(3*i+5)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb5       %F(3*i+5)+HALF*p0%octNb5       %octNb4%F(3*i+5)
                b2 = HALF*p0%octNb5%octNb2%F(3*i+5)+HALF*p0%octNb5%octNb2%octNb4%F(3*i+5)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh4%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a = THR_FOURTH* p0       %F(3*i+5) +QUARTER* p0       %octNb6%F(3*i+5)
                b = THR_FOURTH* p0%octNb1%F(3*i+5) +QUARTER* p0%octNb1%octNb6%F(3*i+5)
                p0%octCh5%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a = THR_FOURTH* p0       %F(3*i+5) +QUARTER* p0       %octNb6%F(3*i+5)
                b = THR_FOURTH* p0%octNb2%F(3*i+5) +QUARTER* p0%octNb2%octNb6%F(3*i+5)
                p0%octCh6%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a1 = HALF*p0       %F(3*i+5)+HALF*p0       %octNb4%F(3*i+5)
                a2 = HALF*p0%octNb1%F(3*i+5)+HALF*p0%octNb1%octNb4%F(3*i+5)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb6       %F(3*i+5)+HALF*p0%octNb6       %octNb4%F(3*i+5)
                b2 = HALF*p0%octNb6%octNb1%F(3*i+5)+HALF*p0%octNb6%octNb1%octNb4%F(3*i+5)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh7%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1 = HALF*p0       %F(3*i+5)+HALF*p0       %octNb4%F(3*i+5)
                a2 = HALF*p0%octNb2%F(3*i+5)+HALF*p0%octNb2%octNb4%F(3*i+5)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb6       %F(3*i+5)+HALF*p0%octNb6       %octNb4%F(3*i+5)
                b2 = HALF*p0%octNb6%octNb2%F(3*i+5)+HALF*p0%octNb6%octNb2%octNb4%F(3*i+5)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh8%C(3*i+5)= THR_FOURTH*a + QUARTER*b
                
                !===  Jz  ===
                !---p0%octCh1---
                a = THR_FOURTH*p0       %F(3*i+6)+QUARTER*p0       %octNb1%F(3*i+6)
                b = THR_FOURTH*p0%octNb3%F(3*i+6)+QUARTER*p0%octNb3%octNb1%F(3*i+6)
                p0%octCh1%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh2---
                a = THR_FOURTH*p0       %F(3*i+6)+QUARTER*p0       %octNb2%F(3*i+6)
                b = THR_FOURTH*p0%octNb3%F(3*i+6)+QUARTER*p0%octNb3%octNb2%F(3*i+6)
                p0%octCh2%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh3---
                a = THR_FOURTH*p0       %F(3*i+6)+QUARTER*p0       %octNb1%F(3*i+6)
                b = THR_FOURTH*p0%octNb4%F(3*i+6)+QUARTER*p0%octNb4%octNb1%F(3*i+6)
                p0%octCh3%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh4---
                a = THR_FOURTH*p0       %F(3*i+6)+QUARTER*p0       %octNb2%F(3*i+6)
                b = THR_FOURTH*p0%octNb4%F(3*i+6)+QUARTER*p0%octNb4%octNb2%F(3*i+6)
                p0%octCh4%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh5---
                a1 = HALF*p0              %F(3*i+6)+HALF*p0              %octNb6%F(3*i+6)
                a2 = HALF*p0       %octNb1%F(3*i+6)+HALF*p0       %octNb1%octNb6%F(3*i+6)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb3       %F(3*i+6)+HALF*p0%octNb3       %octNb6%F(3*i+6)
                b2 = HALF*p0%octNb3%octNb1%F(3*i+6)+HALF*p0%octNb3%octNb1%octNb6%F(3*i+6)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh5%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh6---
                a1 = HALF*p0              %F(3*i+6)+HALF*p0              %octNb6%F(3*i+6)
                a2 = HALF*p0       %octNb2%F(3*i+6)+HALF*p0       %octNb2%octNb6%F(3*i+6)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb3       %F(3*i+6)+HALF*p0%octNb3       %octNb6%F(3*i+6)
                b2 = HALF*p0%octNb3%octNb2%F(3*i+6)+HALF*p0%octNb3%octNb2%octNb6%F(3*i+6)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh6%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh7---
                a1 = HALF*p0              %F(3*i+6)+HALF*p0              %octNb6%F(3*i+6)
                a2 = HALF*p0       %octNb1%F(3*i+6)+HALF*p0       %octNb1%octNb6%F(3*i+6)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb4       %F(3*i+6)+HALF*p0%octNb4       %octNb6%F(3*i+6)
                b2 = HALF*p0%octNb4%octNb1%F(3*i+6)+HALF*p0%octNb4%octNb1%octNb6%F(3*i+6)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh7%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                !---p0%octCh8---
                a1 = HALF*p0              %F(3*i+6)+HALF*p0              %octNb6%F(3*i+6)
                a2 = HALF*p0       %octNb2%F(3*i+6)+HALF*p0       %octNb2%octNb6%F(3*i+6)
                a  = THR_FOURTH*a1+QUARTER*a2
                b1 = HALF*p0%octNb4       %F(3*i+6)+HALF*p0%octNb4       %octNb6%F(3*i+6)
                b2 = HALF*p0%octNb4%octNb2%F(3*i+6)+HALF*p0%octNb4%octNb2%octNb6%F(3*i+6)
                b  = THR_FOURTH*b1+QUARTER*b2
                p0%octCh8%C(3*i+6)= THR_FOURTH*a + QUARTER*b
                end do
                do i=7,9
                    p0%octCh1%F(  i) = (t1*p0%octCh1%C(i)+t2*p0%octCh1%C(i+3))
                    p0%octCh2%F(  i) = (t1*p0%octCh2%C(i)+t2*p0%octCh2%C(i+3))
                    p0%octCh3%F(  i) = (t1*p0%octCh3%C(i)+t2*p0%octCh3%C(i+3))
                    p0%octCh4%F(  i) = (t1*p0%octCh4%C(i)+t2*p0%octCh4%C(i+3))
                    p0%octCh5%F(  i) = (t1*p0%octCh5%C(i)+t2*p0%octCh5%C(i+3))
                    p0%octCh6%F(  i) = (t1*p0%octCh6%C(i)+t2*p0%octCh6%C(i+3))
                    p0%octCh7%F(  i) = (t1*p0%octCh7%C(i)+t2*p0%octCh7%C(i+3))
                    p0%octCh8%F(  i) = (t1*p0%octCh8%C(i)+t2*p0%octCh8%C(i+3))

                    p0%octCh1%F(i+3) = (t3*p0%octCh1%C(i)+t4*p0%octCh1%C(i+3))
                    p0%octCh2%F(i+3) = (t3*p0%octCh2%C(i)+t4*p0%octCh2%C(i+3))
                    p0%octCh3%F(i+3) = (t3*p0%octCh3%C(i)+t4*p0%octCh3%C(i+3))
                    p0%octCh4%F(i+3) = (t3*p0%octCh4%C(i)+t4*p0%octCh4%C(i+3))
                    p0%octCh5%F(i+3) = (t3*p0%octCh5%C(i)+t4*p0%octCh5%C(i+3))
                    p0%octCh6%F(i+3) = (t3*p0%octCh6%C(i)+t4*p0%octCh6%C(i+3))
                    p0%octCh7%F(i+3) = (t3*p0%octCh7%C(i)+t4*p0%octCh7%C(i+3))
                    p0%octCh8%F(i+3) = (t3*p0%octCh8%C(i)+t4*p0%octCh8%C(i+3))
                 end do
              endif
           end do
        endif

        return
        end subroutine interpolate_fieldJ
!
!***********************************************************************
      subroutine interpolate_fieldC(iLv,iF,iK,icon)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ; put interpolated field data to upper level  |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
!        use message_passing_interface,only:rank, ierr
!        use mpi
        implicit none
        integer(kind=4)    :: iLv,iF,iK,is,ie
        integer(kind=4)    :: i,index,icon
        real(kind=8)       :: b1,b2,b3,b4,c1,c2
        real(kind=8)       :: t1,t2,t3,t4
        type(oct), pointer :: p0
        type(oct), pointer :: new1,new2,new3,new4,new5,new6,new7,new8
!        

!        print*, 'entering interpolate_fieldC'

        if(iLv.ge.LvMax) return
        if(iLv.lt.0) return 
        if((maxID(1,iLv).le.minID(1,iLv)) .and. (maxID(3,iLv).le.minID(3,iLv))) return 
!
        is=7 ; ie=12
        if(icon.eq.0) then 
           t1 = -0.5d0 ; t2 = 1.5d0 ; t3 = -1.d0 ; t4 = 2.d0
        elseif(icon.eq.1) then 
           t1 =  1.5d0 ; t2 =-0.5d0 ; t3 =  1.d0 ; t4 = 0.d0
        elseif(icon.eq.2) then 
           t1 =  0.5d0 ; t2 = 0.5d0 ; t3 =  0.d0 ; t4 = 1.d0
        else
           t1 =  0.0d0 ; t2 = 0.0d0 ; t3 =  0.d0 ; t4 = 0.d0
        end if
!

        if(maxID(1,iLv)>minID(1,iLv))then
!!$omp parallel do private(index,p0,i,b1,b2,b3,b4,c1,c2,new1,new2,new3,new4,new5,new6,new7,new8)
           do index=minID(1,iLv),maxID(1,iLv)
              p0 => Mesh(index)
              if(p0%iFLG(iK)==iF) then 

                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb1
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb1
                 new7 => p0%octNb5%octNb3
                 new8 => p0%octNb5%octNb3%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh1%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb2
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb2
                 new7 => p0%octNb5%octNb3
                 new8 => p0%octNb5%octNb3%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh2%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb1
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb1
                 new7 => p0%octNb5%octNb4
                 new8 => p0%octNb5%octNb4%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh3%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb2
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb2
                 new7 => p0%octNb5%octNb4
                 new8 => p0%octNb5%octNb4%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh4%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb1
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb1
                 new7 => p0%octNb6%octNb3
                 new8 => p0%octNb6%octNb3%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh5%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb2
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb2
                 new7 => p0%octNb6%octNb3
                 new8 => p0%octNb6%octNb3%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh6%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb1
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb1
                 new7 => p0%octNb6%octNb4
                 new8 => p0%octNb6%octNb4%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh7%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb2
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb2
                 new7 => p0%octNb6%octNb4
                 new8 => p0%octNb6%octNb4%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh8%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo

                 do i=7,9
                    p0%octCh1%F(  i) = (t1*p0%octCh1%C(i)+t2*p0%octCh1%C(i+3))
                    p0%octCh2%F(  i) = (t1*p0%octCh2%C(i)+t2*p0%octCh2%C(i+3))
                    p0%octCh3%F(  i) = (t1*p0%octCh3%C(i)+t2*p0%octCh3%C(i+3))
                    p0%octCh4%F(  i) = (t1*p0%octCh4%C(i)+t2*p0%octCh4%C(i+3))
                    p0%octCh5%F(  i) = (t1*p0%octCh5%C(i)+t2*p0%octCh5%C(i+3))
                    p0%octCh6%F(  i) = (t1*p0%octCh6%C(i)+t2*p0%octCh6%C(i+3))
                    p0%octCh7%F(  i) = (t1*p0%octCh7%C(i)+t2*p0%octCh7%C(i+3))
                    p0%octCh8%F(  i) = (t1*p0%octCh8%C(i)+t2*p0%octCh8%C(i+3))

                    p0%octCh1%F(i+3) = (t3*p0%octCh1%C(i)+t4*p0%octCh1%C(i+3))
                    p0%octCh2%F(i+3) = (t3*p0%octCh2%C(i)+t4*p0%octCh2%C(i+3))
                    p0%octCh3%F(i+3) = (t3*p0%octCh3%C(i)+t4*p0%octCh3%C(i+3))
                    p0%octCh4%F(i+3) = (t3*p0%octCh4%C(i)+t4*p0%octCh4%C(i+3))
                    p0%octCh5%F(i+3) = (t3*p0%octCh5%C(i)+t4*p0%octCh5%C(i+3))
                    p0%octCh6%F(i+3) = (t3*p0%octCh6%C(i)+t4*p0%octCh6%C(i+3))
                    p0%octCh7%F(i+3) = (t3*p0%octCh7%C(i)+t4*p0%octCh7%C(i+3))
                    p0%octCh8%F(i+3) = (t3*p0%octCh8%C(i)+t4*p0%octCh8%C(i+3))
                 end do
                 nullify(new1) ; nullify(new2)
                 nullify(new3) ; nullify(new4)
                 nullify(new5) ; nullify(new6)
                 nullify(new7) ; nullify(new8)

              endif
           end do
!!$omp end parallel do
        endif
 
        if(maxID(3,iLv)>minID(3,iLv))then
!!$omp parallel do private(index,p0,i,b1,b2,b3,b4,c1,c2,new1,new2,new3,new4,new5,new6,new7,new8)
           do index=minID(3,iLv),maxID(3,iLv)
              p0 => Mesh(index)
              if(p0%iFLG(iK)==iF) then 

                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb1
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb1
                 new7 => p0%octNb5%octNb3
                 new8 => p0%octNb5%octNb3%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh1%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb2
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb2
                 new7 => p0%octNb5%octNb3
                 new8 => p0%octNb5%octNb3%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh2%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb1
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb1
                 new7 => p0%octNb5%octNb4
                 new8 => p0%octNb5%octNb4%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh3%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb2
                 new5 => p0%octNb5
                 new6 => p0%octNb5       %octNb2
                 new7 => p0%octNb5%octNb4
                 new8 => p0%octNb5%octNb4%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh4%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb1
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb1
                 new7 => p0%octNb6%octNb3
                 new8 => p0%octNb6%octNb3%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh5%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb3
                 new4 => p0       %octNb3%octNb2
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb2
                 new7 => p0%octNb6%octNb3
                 new8 => p0%octNb6%octNb3%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh6%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb1
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb1
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb1
                 new7 => p0%octNb6%octNb4
                 new8 => p0%octNb6%octNb4%octNb1
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh7%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo
                 new1 => p0
                 new2 => p0              %octNb2
                 new3 => p0       %octNb4
                 new4 => p0       %octNb4%octNb2
                 new5 => p0%octNb6
                 new6 => p0%octNb6       %octNb2
                 new7 => p0%octNb6%octNb4
                 new8 => p0%octNb6%octNb4%octNb2
                 do i=is,ie
                    b1 = THR_FOURTH*new1%F(i) + QUARTER*new2%F(i)
                    b2 = THR_FOURTH*new3%F(i) + QUARTER*new4%F(i)
                    b3 = THR_FOURTH*new5%F(i) + QUARTER*new6%F(i)
                    b4 = THR_FOURTH*new7%F(i) + QUARTER*new8%F(i)
                    c1 = THR_FOURTH*b1        + QUARTER*b2
                    c2 = THR_FOURTH*b3        + QUARTER*b4
                    p0%octCh8%C(i) = THR_FOURTH*c1 + QUARTER*c2
                 enddo

                 do i=7,9
                    p0%octCh1%F(  i) = (t1*p0%octCh1%C(i)+t2*p0%octCh1%C(i+3))
                    p0%octCh2%F(  i) = (t1*p0%octCh2%C(i)+t2*p0%octCh2%C(i+3))
                    p0%octCh3%F(  i) = (t1*p0%octCh3%C(i)+t2*p0%octCh3%C(i+3))
                    p0%octCh4%F(  i) = (t1*p0%octCh4%C(i)+t2*p0%octCh4%C(i+3))
                    p0%octCh5%F(  i) = (t1*p0%octCh5%C(i)+t2*p0%octCh5%C(i+3))
                    p0%octCh6%F(  i) = (t1*p0%octCh6%C(i)+t2*p0%octCh6%C(i+3))
                    p0%octCh7%F(  i) = (t1*p0%octCh7%C(i)+t2*p0%octCh7%C(i+3))
                    p0%octCh8%F(  i) = (t1*p0%octCh8%C(i)+t2*p0%octCh8%C(i+3))

                    p0%octCh1%F(i+3) = (t3*p0%octCh1%C(i)+t4*p0%octCh1%C(i+3))
                    p0%octCh2%F(i+3) = (t3*p0%octCh2%C(i)+t4*p0%octCh2%C(i+3))
                    p0%octCh3%F(i+3) = (t3*p0%octCh3%C(i)+t4*p0%octCh3%C(i+3))
                    p0%octCh4%F(i+3) = (t3*p0%octCh4%C(i)+t4*p0%octCh4%C(i+3))
                    p0%octCh5%F(i+3) = (t3*p0%octCh5%C(i)+t4*p0%octCh5%C(i+3))
                    p0%octCh6%F(i+3) = (t3*p0%octCh6%C(i)+t4*p0%octCh6%C(i+3))
                    p0%octCh7%F(i+3) = (t3*p0%octCh7%C(i)+t4*p0%octCh7%C(i+3))
                    p0%octCh8%F(i+3) = (t3*p0%octCh8%C(i)+t4*p0%octCh8%C(i+3))
                 end do
                 nullify(new1) ; nullify(new2)
                 nullify(new3) ; nullify(new4)
                 nullify(new5) ; nullify(new6)
                 nullify(new7) ; nullify(new8)

              endif
           end do
!!$omp end parallel do
        endif

        return
        end subroutine interpolate_fieldC
!***********************************************************************
subroutine advance_field(iLv,dstep,iFs,iFe,ch)
  use mpi
  use message_passing_interface
  use oct_set
  use param
  use const
  use init_mesh_size
  implicit none
  integer(kind=4) :: iLv,iFs,iFe,k,i
  integer(kind=4) :: index,minindex,maxindex
  real(kind=8)    :: dstep
  real(kind=8)    :: dvlevel,charge_ch
  real(kind=8)    :: cdxx(3)
  real(kind=8)    :: lamda,ome,tt,xm,ym,zm,xp,yp,zp
  character*1     :: ch
  type(oct), pointer :: p0
  
  !if(debugMode>=1)print *,"entering advance_field ",ch,rank
  
  maxindex = MaxID(1,iLv)
  minindex = MinID(1,iLV)  
  if(maxindex.le.minindex) return 

  if(ch=="E") then 
  ! -- Set spacing parameters --
!     call removeBoundaryJ(iLv)
     dvlevel   = 0.5d0**(iLv)
     PICnumber(1) = npart_per_cell
     charge = 1/omega
     do i=1,3
        cdx(i)=charge*dx(i)
     end do
     charge_ch = charge
     do k=1,3
        cdxx(k)   = charge_ch*dvlevel*dt
     enddo

     if(iLv>=1) then 
!!$omp parallel do private(index,p0) shared(minindex,maxindex,Mesh,dstep,dtdx,cdxx,iFs,iFe)
        do index=minindex,maxindex
           p0 => Mesh(index)
           if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
             p0%F(1) = p0%F(1)                                         &
                     + dstep*dtdx(2) * ( p0%octNb4%F(6) - p0%F(6) )    &
                     - dstep*dtdx(3) * ( p0%octNb6%F(5) - p0%F(5) )    &
                     - dstep*cdxx(1) *   p0%F(7)
             p0%F(2) = p0%F(2)                                         &
                     + dstep*dtdx(3) * ( p0%octNb6%F(4) - p0%F(4) )    &
                     - dstep*dtdx(1) * ( p0%octNb2%F(6) - p0%F(6) )    &
                     - dstep*cdxx(2) *   p0%F(8)
             p0%F(3) = p0%F(3)                                         &
                     + dstep*dtdx(1) * ( p0%octNb2%F(5) - p0%F(5) )    &
                     - dstep*dtdx(2) * ( p0%octNb4%F(4) - p0%F(4) )    &
                     - dstep*cdxx(3) *   p0%F(9)
          endif
       enddo
!!$omp end parallel do 
    else
!!$omp parallel do private(index,p0) shared(minindex,maxindex,Mesh,dstep,dtdx,cdxx,iFe)
       do index=minindex,maxindex
          p0 => Mesh(index)
          if(p0%iFLG(1)<=iFe) then
             p0%F(1) = p0%F(1)                                         &
                     + dstep*dtdx(2) * ( p0%octNb4%F(6) - p0%F(6) )    &
                     - dstep*dtdx(3) * ( p0%octNb6%F(5) - p0%F(5) )    &
                     - dstep*cdxx(1) *   p0%F(7)
             p0%F(2) = p0%F(2)                                         &
                     + dstep*dtdx(3) * ( p0%octNb6%F(4) - p0%F(4) )    &
                     - dstep*dtdx(1) * ( p0%octNb2%F(6) - p0%F(6) )    & 
                     - dstep*cdxx(2) *   p0%F(8)
             p0%F(3) = p0%F(3)                                         &
                     + dstep*dtdx(1) * ( p0%octNb2%F(5) - p0%F(5) )    &
                     - dstep*dtdx(2) * ( p0%octNb4%F(4) - p0%F(4) )    &    
                     - dstep*cdxx(3) *   p0%F(9)
          endif
       enddo
!!$omp end parallel do
    end if
    call removeBoundaryE(iLv)

 else if(ch=="B") then 
    if(iLv>=1) then 
!!$omp parallel do private(index,p0) shared(minindex,maxindex,Mesh,dstep,dtdx,cdxx,iFs,iFe)
       do index=minindex,maxindex
          p0 => Mesh(index)
          if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
             p0%F(4) = p0%F(4)                                         &
                     - dstep*dt2dx(2) * ( p0%F(3) - p0%octNb3%F(3) )   &
                     + dstep*dt2dx(3) * ( p0%F(2) - p0%octNb5%F(2) )  
             p0%F(5) = p0%F(5)                                         &
                     - dstep*dt2dx(3) * ( p0%F(1) - p0%octNb5%F(1) )   &
                     + dstep*dt2dx(1) * ( p0%F(3) - p0%octNb1%F(3) )
             p0%F(6) = p0%F(6)                                         &
                     - dstep*dt2dx(1) * ( p0%F(2) - p0%octNb1%F(2) )   &
                     + dstep*dt2dx(2) * ( p0%F(1) - p0%octNb3%F(1) )
          endif
       enddo
!!$omp end parallel do
    else
!!$omp parallel do private(index,p0) shared(minindex,maxindex,Mesh,dstep,dtdx,cdxx,iFe)   
       do index=minindex, maxindex
          p0 => Mesh(index)
          if(p0%iFLG(1)<=iFe) then
             p0%F(4) = p0%F(4)                                         &
                     - dstep*dt2dx(2) * ( p0%F(3) - p0%octNb3%F(3) )   &
                     + dstep*dt2dx(3) * ( p0%F(2) - p0%octNb5%F(2) )  
             p0%F(5) = p0%F(5)                                         &
                    - dstep*dt2dx(3) * ( p0%F(1) - p0%octNb5%F(1) )    &
                    + dstep*dt2dx(1) * ( p0%F(3) - p0%octNb1%F(3) )
             p0%F(6) = p0%F(6)                                         &
                     - dstep*dt2dx(1) * ( p0%F(2) - p0%octNb1%F(2) )   &
                     + dstep*dt2dx(2) * ( p0%F(1) - p0%octNb3%F(1) )
          endif
       enddo
!!$omp end parallel do
    end if
    call removeBoundaryB(iLv)

 else if(ch=="J") then 
    lamda=10.d0*dx(1)
    ome=2*PI/lamda
    tt=dble(istep)*dt
!    if(rank==0)then
!        print *,"ome=",ome," tt=",tt
!    end if
    xm=(0.5d0*dble(NXB*NXR)+0.5d0)*dx(1)-0.3d0*dx(1)
    xp=(0.5d0*dble(NXB*NXR)+0.5d0)*dx(1)+0.3d0*dx(1)
    ym=(0.5d0*dble(NYB*NYR)+0.5d0)*dx(2)-0.3d0*dx(2)
    yp=(0.5d0*dble(NYB*NYR)+0.5d0)*dx(2)+0.3d0*dx(2)
    zm=(0.5d0*dble(NZB*NZR)+0.5d0)*dx(3)-0.3d0*dx(3)
    zp=(0.5d0*dble(NZB*NZR)+0.5d0)*dx(3)+0.3d0*dx(3)
    do index=minindex, maxindex
        p0 => Mesh(index)
        if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
            if((p0%rPos(1)>=xm).and.(p0%rPos(1)<xp))then
            if((p0%rPos(2)>=ym).and.(p0%rPos(2)<yp))then
!            if((p0%rPos(3)>=zm).and.(p0%rPos(3)<zp))then
                p0%F(7)=0.d0!1.d0*dsin(ome*tt)
                p0%F(8)=0.d0!1.d0*dsin(ome*tt)
                p0%F(9)=1.d0*dsin(ome*tt)
!            end if
            end if
            end if
        endif
    end do
 end if
 return
end subroutine advance_field

!***********************************************************************
      subroutine average_field(iLv,iFs,iFe,ch)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ; get averaged field data from finer levels    |
! +-------------------------------------------------------------------+
!***********************************************************************
!This routine is for averaging (1-6)
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4)    :: iLv,iFs,iFe,is,ie
        integer(kind=4)    :: k,index
        character*1 ch
        type(oct), pointer :: p0

        if(iLv==LvMax) return
        !if(maxID(1,iLv).le.minID(1,iLv)) return 
        if(maxID(1,iLv)<=minID(1,iLv).and.maxID(3,iLv)<=minID(3,iLv)) return 
!
        !if(debugMode>=1)print *,"entering average_field",ch,rank

        call set_index(ch,is,ie)
! For E, k varies from 1-3 while for B, k varies from 4-6
! Averaging happens when iFLG > 0

        if(iFs.ne.iFe) then 
           if(maxID(1,iLv)>minID(1,iLv))then
!!$omp parallel do private(index,k,p0)
                do index=minID(1,iLv),maxID(1,iLv)
                 p0 => Mesh(index)
                 
                 if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
                    do k=is,ie 
                       p0%F(k) = ( p0%octCh1%F(k) + p0%octCh2%F(k)       &
                                 + p0%octCh3%F(k) + p0%octCh4%F(k)       &
                                 + p0%octCh5%F(k) + p0%octCh6%F(k)       &
                                 + p0%octCh7%F(k) + p0%octCh8%F(k) )/8.d0

                    enddo
                 endif
              end do
!!$omp end parallel do
              endif

           if(maxID(3,iLv)>minID(3,iLv))then
!!$omp parallel do private(index,k,p0)
              do index=minID(3,iLv),maxID(3,iLv)
                 p0 => Mesh(index)
                 if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
                    do k=is,ie 
                       p0%F(k) = ( p0%octCh1%F(k) + p0%octCh2%F(k)       &
                                 + p0%octCh3%F(k) + p0%octCh4%F(k)       &
                                 + p0%octCh5%F(k) + p0%octCh6%F(k)       &
                                 + p0%octCh7%F(k) + p0%octCh8%F(k) )/8.d0

                    enddo
                 endif
              end do
!!$omp end parallel do
           endif

        else
         
           if(maxID(1,iLv)>minID(1,iLv))then
!!$omp parallel do private(index,k,p0)
               do index=minID(1,iLv),maxID(1,iLv)
                 p0 => Mesh(index)
                 if(p0%iFLG(1)==iFs) then 
                    do k=is,ie 
                       p0%F(k) = ( p0%octCh1%F(k) + p0%octCh2%F(k)       &
                                 + p0%octCh3%F(k) + p0%octCh4%F(k)       &
                                 + p0%octCh5%F(k) + p0%octCh6%F(k)       &
                                 + p0%octCh7%F(k) + p0%octCh8%F(k) )/8.d0

                    enddo
                 endif
              end do
!!$omp end parallel do
           endif

           if(maxID(3,iLv)>minID(3,iLv))then
!!$omp parallel do private(index,k,p0)
                do index=minID(3,iLv),maxID(3,iLv)
                 p0 => Mesh(index)
                 if(p0%iFLG(1)==iFs) then 
                    do k=is,ie 
                       p0%F(k) = ( p0%octCh1%F(k) + p0%octCh2%F(k)       &
                                 + p0%octCh3%F(k) + p0%octCh4%F(k)       &
                                 + p0%octCh5%F(k) + p0%octCh6%F(k)       &
                                 + p0%octCh7%F(k) + p0%octCh8%F(k) )/8.d0
                    enddo
                 endif
              end do
!!$omp end parallel do
           endif

        end if

        !if(debugMode>=1)print *,"exiting average_field",ch,rank
        return
        end subroutine average_field
        
        subroutine average_fieldE(iLv,iFs,iFe)
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4)    :: iLv,iFs,iFe
        integer(kind=4)    :: index
        type(oct), pointer :: p0

        if(iLv==LvMax) return
        if(maxID(1,iLv)<=minID(1,iLv).and.maxID(3,iLv)<=minID(3,iLv)) return 

        if(iFs.ne.iFe) then 
           if(maxID(1,iLv)>minID(1,iLv))then
                do index=minID(1,iLv),maxID(1,iLv)
                 p0 => Mesh(index)
                 if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
                    !===Ex===
                    p0%F(1)=QUARTER*(p0%octCh1%F(1)+p0%octCh3%F(1)+p0%octCh5%F(1)+p0%octCh7%F(1))
                    !===Ey===
                    p0%F(2)=QUARTER*(p0%octCh1%F(2)+p0%octCh2%F(2)+p0%octCh5%F(2)+p0%octCh6%F(2))
                    !===Ez===
                    p0%F(3)=QUARTER*(p0%octCh1%F(3)+p0%octCh2%F(3)+p0%octCh3%F(3)+p0%octCh4%F(3))
                 endif
              end do
           endif
           if(maxID(3,iLv)>minID(3,iLv))then
              do index=minID(3,iLv),maxID(3,iLv)
                 p0 => Mesh(index)
                 if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
                    !===Ex===
                    p0%F(1)=QUARTER*(p0%octCh1%F(1)+p0%octCh3%F(1)+p0%octCh5%F(1)+p0%octCh7%F(1))
                    !===Ey===
                    p0%F(2)=QUARTER*(p0%octCh1%F(2)+p0%octCh2%F(2)+p0%octCh5%F(2)+p0%octCh6%F(2))
                    !===Ez===
                    p0%F(3)=QUARTER*(p0%octCh1%F(3)+p0%octCh2%F(3)+p0%octCh3%F(3)+p0%octCh4%F(3))
                 endif
              end do
           endif
        else
           if(maxID(1,iLv)>minID(1,iLv))then
               do index=minID(1,iLv),maxID(1,iLv)
                 p0 => Mesh(index)
                 if(p0%iFLG(1)==iFs) then 
                    !===Ex===
                    p0%F(1)=QUARTER*(p0%octCh1%F(1)+p0%octCh3%F(1)+p0%octCh5%F(1)+p0%octCh7%F(1))
                    !===Ey===
                    p0%F(2)=QUARTER*(p0%octCh1%F(2)+p0%octCh2%F(2)+p0%octCh5%F(2)+p0%octCh6%F(2))
                    !===Ez===
                    p0%F(3)=QUARTER*(p0%octCh1%F(3)+p0%octCh2%F(3)+p0%octCh3%F(3)+p0%octCh4%F(3))
                 endif
              end do
           endif
           if(maxID(3,iLv)>minID(3,iLv))then
                do index=minID(3,iLv),maxID(3,iLv)
                 p0 => Mesh(index)
                 if(p0%iFLG(1)==iFs) then 
                    !===Ex===
                    p0%F(1)=QUARTER*(p0%octCh1%F(1)+p0%octCh3%F(1)+p0%octCh5%F(1)+p0%octCh7%F(1))
                    !===Ey===
                    p0%F(2)=QUARTER*(p0%octCh1%F(2)+p0%octCh2%F(2)+p0%octCh5%F(2)+p0%octCh6%F(2))
                    !===Ez===
                    p0%F(3)=QUARTER*(p0%octCh1%F(3)+p0%octCh2%F(3)+p0%octCh3%F(3)+p0%octCh4%F(3))
                 endif
              end do
           endif
        end if

        return
        end subroutine average_fieldE
        
        subroutine average_fieldB(iLv,iFs,iFe)
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4)    :: iLv,iFs,iFe
        integer(kind=4)    :: index
        type(oct), pointer :: p0

        if(iLv==LvMax) return
        if(maxID(1,iLv)<=minID(1,iLv).and.maxID(3,iLv)<=minID(3,iLv)) return 

        if(iFs.ne.iFe) then 
           if(maxID(1,iLv)>minID(1,iLv))then
                do index=minID(1,iLv),maxID(1,iLv)
                 p0 => Mesh(index)
                 if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
                    !===Bx===
                    p0%F(4)=HALF*(p0%octCh1%F(4)+p0%octCh2%F(4))
                    !===By===
                    p0%F(5)=HALF*(p0%octCh1%F(5)+p0%octCh3%F(5))
                    !===Bz===
                    p0%F(6)=HALF*(p0%octCh1%F(6)+p0%octCh5%F(6))
                 endif
              end do
           endif
           if(maxID(3,iLv)>minID(3,iLv))then
              do index=minID(3,iLv),maxID(3,iLv)
                 p0 => Mesh(index)
                 if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
                    !===Bx===
                    p0%F(4)=HALF*(p0%octCh1%F(4)+p0%octCh2%F(4))
                    !===By===
                    p0%F(5)=HALF*(p0%octCh1%F(5)+p0%octCh3%F(5))
                    !===Bz===
                    p0%F(6)=HALF*(p0%octCh1%F(6)+p0%octCh5%F(6))
                 endif
              end do
           endif
        else
           if(maxID(1,iLv)>minID(1,iLv))then
               do index=minID(1,iLv),maxID(1,iLv)
                 p0 => Mesh(index)
                 if(p0%iFLG(1)==iFs) then 
                    !===Bx===
                    p0%F(4)=HALF*(p0%octCh1%F(4)+p0%octCh2%F(4))
                    !===By===
                    p0%F(5)=HALF*(p0%octCh1%F(5)+p0%octCh3%F(5))
                    !===Bz===
                    p0%F(6)=HALF*(p0%octCh1%F(6)+p0%octCh5%F(6))
                 endif
              end do
           endif
           if(maxID(3,iLv)>minID(3,iLv))then
                do index=minID(3,iLv),maxID(3,iLv)
                 p0 => Mesh(index)
                 if(p0%iFLG(1)==iFs) then 
                    !===Bx===
                    p0%F(4)=HALF*(p0%octCh1%F(4)+p0%octCh2%F(4))
                    !===By===
                    p0%F(5)=HALF*(p0%octCh1%F(5)+p0%octCh3%F(5))
                    !===Bz===
                    p0%F(6)=HALF*(p0%octCh1%F(6)+p0%octCh5%F(6))
                 endif
              end do
           endif
        end if

        return
        end subroutine average_fieldB
        
        subroutine average_fieldJ(iLv,iFs,iFe,icon)
        use oct_set
        use param
        use const
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4)    :: iLv,iFs,iFe,is,ie
        integer(kind=4)    :: index,icon,ii
        type(oct), pointer :: p0
        
        is=0; ie=0; ii=0
        if(iLv.ge.LvMax) return
        if(iLv.lt.0) return 
        if(maxID(1,iLv)<=minID(1,iLv).and.maxID(3,iLv)<=minID(3,iLv)) return 
        
        if(icon.eq.0) then 
           is=7 ; ie=9 ; ii=3
        elseif(icon.eq.1) then 
           is=10; ie=12 ; ii=0
        end if
        
        if(iFs.ne.iFe) then 
           if(maxID(1,iLv)>minID(1,iLv))then
                do index=minID(1,iLv),maxID(1,iLv)
                 p0 => Mesh(index)
                 if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
                    if(icon==0)then
                        !===Jx===
                        p0%F(7)=QUARTER*(p0%octCh1%F(10)+p0%octCh3%F(10)+p0%octCh5%F(10)+p0%octCh7%F(10))
                        !===Jy===
                        p0%F(8)=QUARTER*(p0%octCh1%F(11)+p0%octCh2%F(11)+p0%octCh5%F(11)+p0%octCh6%F(11))
                        !===Jz===
                        p0%F(9)=QUARTER*(p0%octCh1%F(12)+p0%octCh2%F(12)+p0%octCh3%F(12)+p0%octCh4%F(12))
                    else if(icon==1)then
                        !===Jx===
                        p0%F(10)=QUARTER*(p0%octCh1%F(10)+p0%octCh3%F(10)+p0%octCh5%F(10)+p0%octCh7%F(10))
                        !===Jy===
                        p0%F(11)=QUARTER*(p0%octCh1%F(11)+p0%octCh2%F(11)+p0%octCh5%F(11)+p0%octCh6%F(11))
                        !===Jz===
                        p0%F(12)=QUARTER*(p0%octCh1%F(12)+p0%octCh2%F(12)+p0%octCh3%F(12)+p0%octCh4%F(12))
                    end if
                 endif
              end do
           endif
           if(maxID(3,iLv)>minID(3,iLv))then
              do index=minID(3,iLv),maxID(3,iLv)
                 p0 => Mesh(index)
                 if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
                   if(icon==0)then
                        !===Jx===
                        p0%F(7)=QUARTER*(p0%octCh1%F(10)+p0%octCh3%F(10)+p0%octCh5%F(10)+p0%octCh7%F(10))
                        !===Jy===
                        p0%F(8)=QUARTER*(p0%octCh1%F(11)+p0%octCh2%F(11)+p0%octCh5%F(11)+p0%octCh6%F(11))
                        !===Jz===
                        p0%F(9)=QUARTER*(p0%octCh1%F(12)+p0%octCh2%F(12)+p0%octCh3%F(12)+p0%octCh4%F(12))
                    else if(icon==1)then
                        !===Jx===
                        p0%F(10)=QUARTER*(p0%octCh1%F(10)+p0%octCh3%F(10)+p0%octCh5%F(10)+p0%octCh7%F(10))
                        !===Jy===
                        p0%F(11)=QUARTER*(p0%octCh1%F(11)+p0%octCh2%F(11)+p0%octCh5%F(11)+p0%octCh6%F(11))
                        !===Jz===
                        p0%F(12)=QUARTER*(p0%octCh1%F(12)+p0%octCh2%F(12)+p0%octCh3%F(12)+p0%octCh4%F(12))
                    end if
                 endif
              end do
           endif
        else
           if(maxID(1,iLv)>minID(1,iLv))then
               do index=minID(1,iLv),maxID(1,iLv)
                 p0 => Mesh(index)
                 if(p0%iFLG(1)==iFs) then 
                    if(icon==0)then
                        !===Jx===
                        p0%F(7)=QUARTER*(p0%octCh1%F(10)+p0%octCh3%F(10)+p0%octCh5%F(10)+p0%octCh7%F(10))
                        !===Jy===
                        p0%F(8)=QUARTER*(p0%octCh1%F(11)+p0%octCh2%F(11)+p0%octCh5%F(11)+p0%octCh6%F(11))
                        !===Jz===
                        p0%F(9)=QUARTER*(p0%octCh1%F(12)+p0%octCh2%F(12)+p0%octCh3%F(12)+p0%octCh4%F(12))
                    else if(icon==1)then
                        !===Jx===
                        p0%F(10)=QUARTER*(p0%octCh1%F(10)+p0%octCh3%F(10)+p0%octCh5%F(10)+p0%octCh7%F(10))
                        !===Jy===
                        p0%F(11)=QUARTER*(p0%octCh1%F(11)+p0%octCh2%F(11)+p0%octCh5%F(11)+p0%octCh6%F(11))
                        !===Jz===
                        p0%F(12)=QUARTER*(p0%octCh1%F(12)+p0%octCh2%F(12)+p0%octCh3%F(12)+p0%octCh4%F(12))
                    end if
                 endif
              end do
           endif
           if(maxID(3,iLv)>minID(3,iLv))then
                do index=minID(3,iLv),maxID(3,iLv)
                 p0 => Mesh(index)
                 if(p0%iFLG(1)==iFs) then 
                   if(icon==0)then
                        !===Jx===
                        p0%F(7)=QUARTER*(p0%octCh1%F(10)+p0%octCh3%F(10)+p0%octCh5%F(10)+p0%octCh7%F(10))
                        !===Jy===
                        p0%F(8)=QUARTER*(p0%octCh1%F(11)+p0%octCh2%F(11)+p0%octCh5%F(11)+p0%octCh6%F(11))
                        !===Jz===
                        p0%F(9)=QUARTER*(p0%octCh1%F(12)+p0%octCh2%F(12)+p0%octCh3%F(12)+p0%octCh4%F(12))
                    else if(icon==1)then
                        !===Jx===
                        p0%F(10)=QUARTER*(p0%octCh1%F(10)+p0%octCh3%F(10)+p0%octCh5%F(10)+p0%octCh7%F(10))
                        !===Jy===
                        p0%F(11)=QUARTER*(p0%octCh1%F(11)+p0%octCh2%F(11)+p0%octCh5%F(11)+p0%octCh6%F(11))
                        !===Jz===
                        p0%F(12)=QUARTER*(p0%octCh1%F(12)+p0%octCh2%F(12)+p0%octCh3%F(12)+p0%octCh4%F(12))
                    end if
                 endif
              end do
           endif
        end if

        return
        end subroutine average_fieldJ
!
!***********************************************************************
      subroutine average_fieldC(iLv,iFs,iFe,icon)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ; get averaged field data from upper level    |
! +-------------------------------------------------------------------+
!***********************************************************************
!This routine is for averaging F(7-12)
        use oct_set
        use param
        use const
        use init_mesh_size
        implicit none
        integer(kind=4)    :: iLv,iFs,iFe,is,ie,icon
        integer(kind=4)    :: k,index,i,ii
        type(oct), pointer :: p0

        is=0; ie=0; ii=0

        if(iLv.ge.LvMax) return
        if(iLv.lt.0) return 
        !if(maxID(1,iLv).le.minID(1,iLv)) return 
        if((maxID(1,iLv).le.minID(1,iLv)) .and. (maxID(3,iLv).le.minID(3,iLv))) return 
!
        if(icon.eq.0) then 
           is=7 ; ie=9 ; ii=3
        elseif(icon.eq.1) then 
           is=10; ie=12 ; ii=0
        end if
!
        if(iFs.ne.iFe) then 
          
           if(maxID(1,iLv)>minID(1,iLv))then
!!$omp parallel do private(index,k,p0,i)
              do index=minID(1,iLv),maxID(1,iLv)
                 p0 => Mesh(index)
                 if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
                    do k=is,ie 
                       i=k+ii
                       p0%F(k) = ( p0%octCh1%F(i) + p0%octCh2%F(i)       &
                                 + p0%octCh3%F(i) + p0%octCh4%F(i)       &
                                 + p0%octCh5%F(i) + p0%octCh6%F(i)       &
                                 + p0%octCh7%F(i) + p0%octCh8%F(i) )/8.d0
                    enddo
                 endif
              end do
!!$omp end parallel do
           endif

           if(maxID(3,iLv)>minID(3,iLv))then
!!$omp parallel do private(index,k,p0,i)
              do index=minID(3,iLv),maxID(3,iLv)
                 p0 => Mesh(index)
                 if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
                    do k=is,ie 
                       i=k+ii
                       p0%F(k) = ( p0%octCh1%F(i) + p0%octCh2%F(i)       &
                                 + p0%octCh3%F(i) + p0%octCh4%F(i)       &
                                 + p0%octCh5%F(i) + p0%octCh6%F(i)       &
                                 + p0%octCh7%F(i) + p0%octCh8%F(i) )/8.d0
                    enddo
                 endif
              end do
!!$omp end parallel do
           endif

        else
      
           if(maxID(1,iLv)>minID(1,iLv))then
!!$omp parallel do private(index,k,p0,i)
              do index=minID(3,iLv),maxID(3,iLv)
                 p0 => Mesh(index)
                 if(p0%iFLG(1)==iFs) then
                    do k=is,ie 
                       i=k+ii
                       p0%F(k) = ( p0%octCh1%F(i) + p0%octCh2%F(i)       &
                            + p0%octCh3%F(i) + p0%octCh4%F(i)       &
                            + p0%octCh5%F(i) + p0%octCh6%F(i)       &
                            + p0%octCh7%F(i) + p0%octCh8%F(i) )/8.d0
                    enddo
                 endif
              end do
!!$omp end parallel do
           endif

           if(maxID(1,iLv)>minID(1,iLv))then
!!$omp parallel do private(index,k,p0,i)
              do index=minID(1,iLv),maxID(1,iLv)
                 p0 => Mesh(index)
                 if(p0%iFLG(1)==iFs) then
                    do k=is,ie 
                       i=k+ii
                       p0%F(k) = ( p0%octCh1%F(i) + p0%octCh2%F(i)       &
                            + p0%octCh3%F(i) + p0%octCh4%F(i)       &
                            + p0%octCh5%F(i) + p0%octCh6%F(i)       &
                            + p0%octCh7%F(i) + p0%octCh8%F(i) )/8.d0
                    enddo
                 endif
              end do
!!$omp end parallel do
           endif

        end if

        return
        end subroutine average_fieldC
!
!***********************************************************************
      subroutine average_fieldZ(iLv,iFs,iFe)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ; get averaged field data from upper level    |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
        implicit none
        integer(kind=4)    :: iLv,iFs,iFe
        integer(kind=4)    :: k,index
        type(oct), pointer :: p0

        if(iLv==LvMax) return
        if(maxID(1,iLv).le.minID(1,iLv)) return 
!
        if(iFs.ne.iFe) then 
!!$omp parallel do private(index,k,p0)
        do index=minID(1,iLv),maxID(1,iLv)
           p0 => Mesh(index)
           if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
              do k=1,Ionsorts 
                 p0%Z(k) = ( p0%octCh1%Z(k) + p0%octCh2%Z(k)       &
                           + p0%octCh3%Z(k) + p0%octCh4%Z(k)       &
                           + p0%octCh5%Z(k) + p0%octCh6%Z(k)       &
                           + p0%octCh7%Z(k) + p0%octCh8%Z(k) )/8.d0
              enddo
           endif
        end do
!!$omp end parallel do
        else
!!$omp parallel do private(index,k,p0)
        do index=minID(1,iLv),maxID(1,iLv)
           p0 => Mesh(index)
           if(p0%iFLG(1)==iFs) then
              do k=1,Ionsorts 
                 p0%Z(k) = ( p0%octCh1%Z(k) + p0%octCh2%Z(k)       &
                           + p0%octCh3%Z(k) + p0%octCh4%Z(k)       &
                           + p0%octCh5%Z(k) + p0%octCh6%Z(k)       &
                           + p0%octCh7%Z(k) + p0%octCh8%Z(k) )/8.d0
              enddo
           endif
        end do
!!$omp end parallel do
        end if

        return
        end subroutine average_fieldZ
!
!************************************************************************
      subroutine set_index(ch,is,ie)
!************************************************************************
        implicit none
        character*1 ch
        integer(kind=4) :: is,ie
!
        if(ch=="E") then 
           is=1 ; ie=3
        else if(ch=="B") then 
           is=4 ; ie=6
        else if(ch=="C") then 
           is=7 ; ie=9
        else if(ch=="D") then 
           is=10; ie=12
        else if(ch=="F") then 
           is=1 ; ie=6
        end if
!
        return 
        end subroutine set_index
!
!============================================================================
!      subroutine restart(time)
! +-------------------------------------------------------------------+
! |                                                                   |
! |  field correction at istep   (for debug)                          |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************
!!$        use oct_set
!!$        use param
!!$        use const
!!$        use init_mesh_size
!!$        use init_condition
!-----------------------------------------
!!$        implicit none
!!$        integer(kind=4)   :: index,i
!!$        real(kind=8)      :: time
!!$        real(kind=8)      :: pos(3)
!!$        real(kind=8)      :: phase
!!$!------------------------------------------
!!$        type(oct), pointer :: p0
!!$!------------------------------------------ 
!!$!
!!$! == BEGIN: ============================================================
!!$!
!!$        do index=MinID(1,0),MaxID(1,LvMax)
!!$           p0 => Mesh(index)
!!$!
!!$           pos(1) = p0%rPOS(1)
!!$           pos(2) = p0%rPOS(2)
!!$           pos(3) = p0%rPOS(3)
!!$!
!!$           phase = pi2*(kiniN(1)*pos(1)/R_lim(1,1)       &
!!$                       +kiniN(2)*pos(2)/R_lim(1,2)       &
!!$                       +kiniN(3)*pos(3)/R_lim(1,3))
!!$!
!!$           p0%F(1) = EiniN(1)*dcos(phase+phi)
!!$           p0%F(2) = EiniN(2)*dcos(phase+phi)
!!$           p0%F(3) = EiniN(3)*dcos(phase+phi)
!!$           p0%F(4) = BiniN(1)*dcos(phase+phi)
!!$           p0%F(5) = BiniN(2)*dcos(phase+phi)
!!$           p0%F(6) = BiniN(3)*dcos(phase+phi)
!!$        enddo
!!$!
!!$        do i=1,3
!!$          dtdx(i)  = dt/dx(i)
!!$          dt2dx(i) = HALF*dt/dx(i)
!!$          odx(i)   = ONE/dx(i)
!!$        end do
!!$!
!!$        time = 0.d0
!!$!
!!$        print *, ' '
!!$        print *, '************** START THE MAIN LOOP **************'
!!$        print *, ' '
!!$!
!!$        return
!!$      end subroutine restart
!
! ================================================================================
!     subroutines for particles
! ================================================================================
!
!***********************************************************************
      subroutine copy_particleD(iLv,iF,iK,icon)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ;  copy particles to the child octs           |
! |     if icon==0 ; particle is deleted after copy process (move)    |
! |     if icon==1 ; particle is not deleted after copy process (copy)|
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use particle_set
        use param
        use init_mesh_size
        use const
        use message_passing_interface
        implicit none
        integer(kind=4)    :: iLv,iF,iK,icon,k
        integer(kind=4)    :: index,indexP,indexR
        integer(kind=4)    :: ipx,ipy,ipz,ich,ILv2
        integer(kind=4)    :: octP1,octP2,octP3,octP4
        integer(kind=4)    :: octP5,octP6,octP7,octP8
        integer(kind=4)    :: Isort,Ioct
        real(kind=8)       :: R(9),pos(3)
! -- pointers --
        type(  oct), pointer :: p0
        type(prtcl), pointer :: pptc
        type(prtcl), pointer :: PrtList1,PrtList2,PrtList3,PrtList4
        type(prtcl), pointer :: PrtList5,PrtList6,PrtList7,PrtList8

        if(iLv.ge.LvMax) return
        if(iLv.lt.0) return 

        if(maxID(1,iLv+1).eq.minID(1,iLv+1)) return 
        if(maxID(1,iLv  ).eq.minID(1,iLv  )) return 

        iLv2 = iLv+1
        indexP = MaxIP(iLv2)
!!! $omp parallel do private(index,p0,octP1,octP2,octP3,octP4,octP5,octP6,octP7,octP8 &
!!! $omp&                   ,PrtList1,PrtList2,PrtList3,PrtList4,PrtList5,PrtList6,PrtList7,PrtList8 &
!!! $omp&                   ,pos,indexR,pptc,k,R,Isort,ipx,ipy,ipz,ich,Ioct) &
!!! $omp&             shared(minID,maxID,iK,iF,Pesh,Mesh,iLv,iLv2,icon) &
!!! $omp&          reduction(+: indexP)
!!$        if(debugMode>=1)then
!!$           print*,'entering copy_paricleD',rank,iLv 
!!$           print*,"range min-max",MinID(1,iLv),"--",MaxID(1,iLv),rank
!!$        endif
        do index=minID(1,iLv),maxID(1,iLv)
           p0 => Mesh(index)
           if(p0%iFLG(iK)==iF) then 
              !if(p0%prc_bndry==0) then !Not to bring particles down to finer cells at process boundary
!
              octP1 = 0 ; octP2 = 0
              octP3 = 0 ; octP4 = 0
              octP5 = 0 ; octP6 = 0
              octP7 = 0 ; octP8 = 0
              nullify(PrtList1) ; nullify(PrtList2) 
              nullify(PrtList3) ; nullify(PrtList4) 
              nullify(PrtList5) ; nullify(PrtList6) 
              nullify(PrtList7) ; nullify(PrtList8) 
!
              pos(1) = p0%rPOS(1)
              pos(2) = p0%rPOS(2)
              pos(3) = p0%rPOS(3)
              indexR = p0%octP
!
              pptc => p0%ptcl%prtnxt
!
             
              do k=1,indexR
!
                 if(pptc%isort.le.0) cycle
!
                 R(1) = pptc%R(1) ; R(2) = pptc%R(2) ; R(3) = pptc%R(3)
                 R(4) =(pptc%R(4)+pptc%R(7))*0.5d0
                 R(5) =(pptc%R(5)+pptc%R(8))*0.5d0
                 R(6) =(pptc%R(6)+pptc%R(9))*0.5d0
                 R(7) = pptc%R(7) ; R(8) = pptc%R(8) ; R(9) = pptc%R(9)
                 Isort= pptc%Isort
                 ipx=(idnint(sign(1.d0,R(1)-pos(1)))+1)/2
                 ipy=(idnint(sign(2.d0,R(2)-pos(2)))+2)/2
                 ipz=(idnint(sign(4.d0,R(3)-pos(3)))+4)/2
                 ich=ipx+ipy+ipz+1
!
                 indexP = indexP + 1
                 if(ich.le.4) then 
                    if(ich == 1) then
                       Ioct = p0%octCh1%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList1,0)
                       PrtList1 => Pesh(indexP,iLv2)
                       octP1 = octP1 + 1
                    elseif(ich == 2) then
                       Ioct = p0%octCh2%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList2,0)
                       PrtList2 => Pesh(indexP,iLv2)
                       octP2 = octP2 + 1
                    elseif(ich == 3) then
                       Ioct = p0%octCh3%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList3,0)
                       PrtList3 => Pesh(indexP,iLv2)
                       octP3 = octP3 + 1
                    else
                       Ioct = p0%octCh4%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList4,0)
                       PrtList4 => Pesh(indexP,iLv2)
                       octP4 = octP4 + 1
                    endif
                 else
                    if(ich == 5) then
                       Ioct = p0%octCh5%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList5,0)
                       PrtList5 => Pesh(indexP,iLv2)
                       octP5 = octP5 + 1
                    elseif(ich == 6) then
                       Ioct = p0%octCh6%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList6,0)
                       PrtList6 => Pesh(indexP,iLv2)
                       octP6 = octP6 + 1
                    elseif(ich == 7) then
                       Ioct = p0%octCh7%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList7,0)
                       PrtList7 => Pesh(indexP,iLv2)
                       octP7 = octP7 + 1
                    else
                       Ioct = p0%octCh8%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList8,0)
                       PrtList8 => Pesh(indexP,iLv2)
                       octP8 = octP8 + 1
                    endif
                 endif
               
                 pptc%isort = icon*(pptc%isort)
                 pptc => pptc%prtnxt
              enddo
              p0%octCh1%ptcl%PrtNxt => PrtList1 ; p0%octCh1%octP = octP1
              p0%octCh2%ptcl%PrtNxt => PrtList2 ; p0%octCh2%octP = octP2
              p0%octCh3%ptcl%PrtNxt => PrtList3 ; p0%octCh3%octP = octP3
              p0%octCh4%ptcl%PrtNxt => PrtList4 ; p0%octCh4%octP = octP4
              p0%octCh5%ptcl%PrtNxt => PrtList5 ; p0%octCh5%octP = octP5
              p0%octCh6%ptcl%PrtNxt => PrtList6 ; p0%octCh6%octP = octP6
              p0%octCh7%ptcl%PrtNxt => PrtList7 ; p0%octCh7%octP = octP7
              p0%octCh8%ptcl%PrtNxt => PrtList8 ; p0%octCh8%octP = octP8
!
              p0%octP = icon*(p0%octP)
              !if(rank==1)print *,"move particle complete",index,p0%octLv,rank
              !endif
           endif
        end do
!!! $omp end parallel do
!
        MaxIP(iLv2)=indexP
        !if(debugMode>=1)print*,'exiting copy_particleD',rank
!
        return
        end subroutine copy_particleD
!
!***********************************************************************
      subroutine copy_particleU(iLv,iF,iK,icon)
! +-------------------------------------------------------------------+
! |     Bring particles up to a Coaser level                          |
! |     if iFLG(iK)==iF ;  copy particles to the parent octs          |
! |     if icon==0 ; particle is deleted after copy process (move)    |
! |     if icon==1 ; particle is not deleted after copy process (copy)|
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use particle_set
        use param
        use init_mesh_size
        use const
        use message_passing_interface
        implicit none
        integer(kind=4)    :: iLv,iF,iK,icon,k,ii
        integer(kind=4)    :: index,indexP,indexR,octP1
        integer(kind=4)    :: Isort,Ioct
        real(kind=8)       :: R(9)
        integer(kind=4)    :: index8(8)
! -- pointers --
        type(  oct), pointer :: p0,p1
        type(prtcl), pointer :: pptc,PrtList1
!        
        if(iLv.ge.LvMax) return
        if(iLv.lt.0) return 
        !if(debugMode>0)print*,'entering copy_particleU',rank,iLv

        if(maxID(1,iLv+1).eq.minID(1,iLv+1)) return 
        if(maxID(1,iLv  ).eq.minID(1,iLv  )) return 

        indexP = MaxIP(iLv)
!! $omp parallel do private(index,p0,octP1,Ioct,index8,ii,indexR,p1,pptc,R,Isort,PrtList1,k) &
!! $omp&             shared(minID,maxID,Mesh,Pesh,icon,iLv,iK,iF) reduction(+: indexP)

        do index=minID(1,iLv),maxID(1,iLv)
           p0 => Mesh(index)
           !           if(p0%iFLG(iK)==iF) then 

           if(p0%iFLG(iK)==iF)  then 
              !if(p0%prc_bndry==0) then! Skip the process boundary buffer as well 
              ! because we don't want to erase particles in Lv=0
              ! by setting octP1=0
                 octP1 = 0
                 nullify(PrtList1)
                 Ioct  = p0%octN
            !
                 index8(1)=p0%octCh1%octN ; index8(2)=p0%octCh2%octN
                 index8(3)=p0%octCh3%octN ; index8(4)=p0%octCh4%octN
                 index8(5)=p0%octCh5%octN ; index8(6)=p0%octCh6%octN
                 index8(7)=p0%octCh7%octN ; index8(8)=p0%octCh8%octN
                 do ii=1,8
                    p1 => Mesh(index8(ii))
                    indexR = p1%octP
                    pptc => p1%ptcl%prtnxt
                    do k=1,indexR
                       if(pptc%isort.le.0) cycle
                       R(1) = pptc%R(1) ; R(2) = pptc%R(2) ; R(3) = pptc%R(3)
                       R(4) = 2.d0*pptc%R(4)-pptc%R(7)
                       R(5) = 2.d0*pptc%R(5)-pptc%R(8)
                       R(6) = 2.d0*pptc%R(6)-pptc%R(9)
                       R(7) = pptc%R(7) ; R(8) = pptc%R(8) ; R(9) = pptc%R(9)
                       Isort= pptc%Isort
                       indexP= indexP+ 1
                       Pesh(indexP,iLv) = prtcl(R,Isort,Ioct,indexP,PrtList1,0)
                       PrtList1 => Pesh(indexP,iLv)
                       octP1 = octP1 + 1
                       pptc%Isort = icon*(pptc%Isort)
                       pptc => pptc%prtnxt
                    end do
                    p1%octP = icon*(p1%octP)
                 end do
                 p0%ptcl%PrtNxt => PrtList1 ; p0%octP = octP1
              !endif
           end if
        end do
!! $omp end parallel do
!

        MaxIP(iLv)=indexP
!
        return
        end subroutine copy_particleU
!
!***********************************************************************
      subroutine delete_particleT(iLv,iFs,iFe)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ;  delete particles of the own oct            |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use particle_set
        use param
        use init_mesh_size

        use message_passing_interface
        implicit none
        integer(kind=4)    :: iLv,iFs,iFe
        integer(kind=4)    :: index,indexR,k
! -- pointers --
        type(oct), pointer :: p0
        type(prtcl), pointer :: pptc

        if(iLv.lt.0) return 
        if(minID(1,iLv)>=maxID(1,iLv))return
        !if(debugMode>=1)print *,"entering delete_particleT iLv=",iLv,"iFs=",iFs,"iFe=",iFe,rank
!!$omp parallel do private(index,p0,indexR,pptc,k) shared(minID,maxID,Mesh,iFs,iFe)
        do index=minID(1,iLv),maxID(1,iLv)
           p0 => Mesh(index)
           if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then             
              indexR = p0%octP
              pptc   => p0%ptcl%prtnxt
              do k=1,indexR
                 pptc%isort = 0
                 pptc => pptc%prtnxt
              enddo
              p0%octP = 0
           endif
        end do
!!$omp end parallel do
!
        !if(debugMode>=1)print *,"exiting delete_particleT iLv=",iLv,"iFs=",iFs,"iFe=",iFe,rank
        return
        end subroutine delete_particleT
!
!***********************************************************************
      subroutine sort_particle
!***********************************************************************
        use oct_set
        use particle_set
        use param
        use init_mesh_size
        use const
        use message_passing_interface
		use Vparam
        implicit none
        integer(kind=4)    :: iLv,i
        integer(kind=4)    :: index,indexP,indexR,index2,indexS
        integer(kind=4)    :: minindex,maxindex,Gminindex,Gmaxindex!,Lmax
! -- pointers --
        type(  oct), pointer :: p0
        type(prtcl), pointer :: pptc
!
        !if(debugMode>0)print*, 'sorting particles...', rank
      
        do iLv=0,LvMax

           indexP = 0
           minindex = minID(1,iLv)               
           maxindex = maxID(1,iLv)
           Gminindex= minID(3,iLv)
           Gmaxindex= maxID(3,iLv)
		  
		   
           !debug
           !print *,"L=",L,"MeshLen=",MeshLength,"GMeshLen=",GMeshLength,rank
!!$        if((maxindex.le.minindex)) return
           if((maxindex.le.minindex) .and. (Gmaxindex .le. Gminindex)) cycle  ! modify
 
           if(maxindex>minindex)then

              !The following loop organizes the particles in Pesh2 in order.
!=======Hybrid Error===========
!!$omp parallel do private(index,p0,indexR,pptc,indexS,i) & 
!!$omp&             shared(iLv,Mesh,Pesh,Pesh2,minindex,maxindex) &
!!$omp&          reduction(+: indexP)
              do index=minindex,maxindex 
                 p0 => Mesh(index)              ! p0 scans Mesh from minindex to maxindex
                 indexR = p0%octP                           
                 pptc => p0%ptcl                                      ! pptc points to a rep particle in the cell 

                 indexP = indexP+1                 ! indexP is the particle number counted from 1st cell in Mesh

                 indexS = pptc%inum                               ! inum is the index of PRE-SORTED particles in Pesh 
                 pptc%inum = indexP                      ! Line up particles in Pesh according to indexP, which is the particle number counted from the 1st cell in Mesh
                 Pesh2(indexP) = Pesh(indexS,iLv)         ! 
                 p0%ptcl => Pesh(indexP,iLv)
                 pptc => pptc%prtnxt
                 do i=1,indexR                                          ! For each cell, put the particle info into Pesh2 in the according order with indexP  
                    if(pptc%isort.gt.0) then                         ! For real particles
                       indexP = indexP + 1                          ! If a real particle is found, count the number
                       indexS = pptc%inum                         ! where is indexP'th particle located in the original Pesh BEFORE SORTING
                       pptc%inum = indexP                          ! new location of the particle in Pesh is indexP
                       Pesh2(indexP) = Pesh(indexS,iLv)    ! Extract particle info from the old Pesh and put it into Pesh 2 in the order of indexP 
                       Pesh2(indexP-1)%prtnxt=>Pesh(indexP,iLv)
                    end if
                    pptc => pptc%prtnxt
                 end do
              end do
!!$omp end parallel do
           endif

           !for GMesh domain
           if(Gmaxindex>Gminindex)then
!!$omp parallel do private(index,p0,indexR,pptc,indexP,indexS,i) & 
!!$omp&             shared(iLv,Mesh,Pesh,Pesh2,Gminindex,Gmaxindex)
              do index=Gminindex,GMaxindex
                 p0 => Mesh(index)
                 indexR = p0%octP
                 pptc => p0%ptcl
                 indexP = indexP+1
                 indexS = pptc%inum
                 pptc%inum = indexP
                 Pesh2(indexP) = Pesh(indexS,iLv)
                 p0%ptcl => Pesh(indexP,iLv)
                 pptc => pptc%prtnxt
                 do i=1,indexR
                    if(pptc%isort>0)then
                       indexP = indexP + 1
                       indexS = pptc%inum
                       pptc%inum = indexP
                       Pesh2(indexP) = Pesh(indexS,iLv)
                       Pesh2(indexP-1)%prtnxt=>Pesh(indexP,iLv)
                    endif
                    pptc => pptc%prtnxt
                 enddo
              enddo
!!$omp end parallel do
           endif

!====================================2016/1/22 Nakano=======================

	        if(iLV == 0)then

!!$omp parallel do private(index,p0,indexR,pptc,indexS,i) & 
!!$omp&             shared(iLv,Mesh,Pesh,Pesh2,minindex,maxindex) &
!!$omp&          reduction(+: indexP)
                do index=1,VMeshSize 
                 p0 => VMesh(index)              
                 indexR = p0%octP                           
                 pptc => p0%ptcl                                      
                 indexP = indexP+1                
                 indexS = pptc%inum                             
                 pptc%inum = indexP                      
                 Pesh2(indexP) = Pesh(indexS,iLv)        
                 p0%ptcl => Pesh(indexP,iLv)
                 pptc => pptc%prtnxt
                   do i=1,indexR                                        
                    if(pptc%isort.gt.0) then                         
                       indexP = indexP + 1                         
                       indexS = pptc%inum                      
                       pptc%inum = indexP                        
                       Pesh2(indexP) = Pesh(indexS,iLv)    
                       Pesh2(indexP-1)%prtnxt=>Pesh(indexP,iLv)
                    end if
                    pptc => pptc%prtnxt
                    end do
                end do
!!$omp end parallel do
           endif

!===========================================================================		   
		   		   
		   
		   
		   
		   
           index2 = indexP   ! index 2 is the total number of particles in the new 

           !!$omp parallel do private(i,pp,p0,inum,pptc) shared(iLv,index2,Mesh,Pesh)


           ! The following Do loop re-connect the particle list elements in Mesh in the sorted order.
!!$           do i=1,index2                                            ! loop in Pesh
!!$              pp => Pesh2(i)                                      ! sample particle in sorted order
!!$              if(pp%isort.eq.-1) then                          ! If the particle is of rep kind
!!$                 inum = pp%ioct                                  !Which cell does this rep particle belong to?
!!$                 p0 => Mesh(inum)                             
!!$                 inum = pp%inum                               ! Now where is the particle located in Pesh? 
!!$                 p0%ptcl => Pesh(inum,iLv)               ! Connect Mesh and Pesh
!!$              end if
!!$              if(associated(pp%prtnxt))then
!!$                 pptc => pp%prtnxt                               ! Next one!
!!$                 inum = pptc%inum                            ! where are you in Pesh?
!!$                 pp%prtnxt => Pesh(inum,iLv)            ! Connect Mesh and Pesh
!!$              end if
!!$           end do

           !!$omp end parallel do
           !
!!$omp parallel do private(i) shared(index2,iLv,Pesh,Pesh2)
           do i=1,index2
              Pesh(i,iLv)=Pesh2(i)                                 
           end do
!!$omp end parallel do
           ! Erase the rest tail of Pesh. (Added by Tatsuki)
!!$omp parallel do private(i) shared(index2,maxip,iLv,Pesh)
           do i=index2+1,maxip(iLv)
              !           if(Pesh(i,iLv) % isort > 0) Mesh(Pesh(i,iLv) % ioct) % octP  = Mesh(Pesh(i,iLv) % ioct) % octP  -1 
              Pesh(i,iLv) % isort =0
           enddo
!!$omp end parallel do
           maxIP(iLv)=index2
		   

		   
		   

        enddo

        !if(debugMode>0)print *,"sort_particle Done, iLv, maxIP",iLv,maxIP(iLv)

        return
        end subroutine sort_particle

!***********************************************************************
        subroutine reconnect_particleP(iLv,iFs,iFe)
! +-------------------------------------------------------------------+
! |                                                                   |
! | Reconnect list structure for particles                            |
! |                                                                   |
! |                                                                   |
! |  Order of octs which the particle escape to.                      |
! |   isgn = [1-27]                                                   |
! |                                                                   |
! |  Original oct : isgn=14                                           |
! |                                                                   |
! |           +-------+-------+-------+                               |
! |          /  25   /  26   /  27   /|                               |
! |         +-------+-------+-------+ |                               |
! |        /  22   /  23   /  24   /| |                               |
! |       +-------+-------+-------+ | +                               |
! |      /  19   /  20   /  21   /| |/|                               |
! |     +-------+-------+-------+ | + |                               |
! |     |       |       |       | |/| |                               |
! |     |  19   |  20   |  21   | + | +                               |
! |     |       |       |       |/| |/|        Z                      |
! |     +-------+-------+-------+ | + |        |   Y                  |
! |     |       |       |       | |/| |        |  /                   |
! |     |  10   |  11   |  12   | + | +        | /                    |
! |     |       |       |       |/| |/         |/                     |
! |     +-------+-------+-------+ | +          o--------X             |
! |     |       |       |       | |/                                  |
! |     |   1   |   2   |   3   | +                                   |
! |     |       |       |       |/                                    |
! |     +-------+-------+-------+                                     |
! |                                                                   |
! |                                                                   |
! |           +-------+-------+-------+                               |
! |          /  19   /  22   /  25   /|                               |
! |         +-------+-------+-------+ |                               |
! |        /  20   /  23   /  26   /| |                               |
! |       +-------+-------+-------+ | +                               |
! |      /  21   /  24   /  27   /| |/|                               |
! |     +-------+-------+-------+ | + |                               |
! |     |       |       |       | |/| |                               |
! |     |  21   |  24   |  27   | + | +                               |
! |     |       |       |       |/| |/|        Z                      |
! |     +-------+-------+-------+ | + |        |                      |
! |     |       |       |       | |/| |        |                      |
! |     |  12   |  15   |  18   | + | +        |                      |
! |     |       |       |       |/| |/         |                      |
! |     +-------+-------+-------+ | +          o--------Y             |
! |     |       |       |       | |/          /                       |
! |     |   3   |   6   |   9   | +          /                        |
! |     |       |       |       |/          /                         |
! |     +-------+-------+-------+          X                          |
! |                                                                   |
! |                                                                   |
! |                                                                   |
! |                                                                   |
! |   Argument parameters:                                            |
! |       i2step  : iLv+1 levels time step                            |
! |       istep   : iLv   levels time step                            |
! |       iLv     : target hierarchical level                         |
! |                                                                   |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        use param
        use message_passing_interface
		use Vparam
!---------------------------------------------------------------
        implicit none
        integer(kind=4) :: indx(3),iLv,iFs,iFe
        integer(kind=4) :: j,k,index,indexP
        integer(kind=4) :: octType, newoctN
        integer(kind=4) :: maxindex,minindex,minP,maxP
        integer(kind=4) :: Gmaxindex,Gminindex
        real(kind=8)    :: dxx(3)
        real(kind=8)    :: dvlevel,Xout,Xin,Yout,Yin,Zout,Zin
        integer(kind=4) :: i,iomp2,iC0

! -- pointers --
        type(oct),   pointer :: p0,ptmp,ptmp1,ptmp2,p1
        type(prtcl), pointer :: pp,pptc
!!$        include 'omp_lib.h'
!---------------------------------------------------------------
        dvlevel  = HALF**iLv
        do k=1,3
           dxx(k) = dvlevel * dx(k)  ! dxx= dx for iLv=0
        enddo
!  reset oct infomation
        maxindex = MaxID(1,iLv)
        minindex = MinID(1,iLv)
        Gmaxindex= MaxID(3,iLv)
        Gminindex= MinID(3,iLv)
        if(maxindex<=minindex .and. Gmaxindex<=Gminindex)return
! At this point all the particles are threaded to Mesh, including the ones who are 
! physically in the GMesh Zone.   
        if(minindex<maxindex)then
!!$omp parallel do private(index,p0) shared(minindex,maxindex,Mesh)
           do index=minindex,maxindex
              p0 => Mesh(index)   
              p0%octP = 0
              do i=1,iomp0
                 p0%iC(i) = 0 
              end do
           end do
!!$omp end parallel do
        endif

        if(Gminindex<Gmaxindex)then
!!$omp parallel do private(index,p0) shared(minindex,maxindex,Mesh)
           do index=Gminindex,Gmaxindex
              p0 => Mesh(index)   
              p0%octP = 0
              do i=1,iomp0
                 p0%iC(i) = 0 
              end do
           end do
!!$omp end parallel do
        endif
		
!=================2016/1/23Nakano======================
	if(iLv == 0) then
		do index=1 , VMeshSize
			p0 => VMesh(index)
			p0%octP = 0
			do i=1,iomp0
                 p0%iC(i) = 0 
            end do
		end do
	end if
!======================================================
        
        maxP = MaxIP(iLv)
        minP = 1
        
        !ADD 20140630
        Xin =R_lim(0,1)+dx(1)!dx(1)*0.5d0**intLv
        if(Model==7)then
            Xout=Dpos(1)-7.d0*dx(1)
        else
            Xout=R_lim(1,1)-dx(1)!*0.5d0**intLv !0.5d0*R_lim(1,1)-dx(1)*0.5d0**intLv
        end if
        Yin =R_lim(0,2)+dx(2)!*0.5d0**intLv
        Yout=R_lim(1,2)-dx(2)!*0.5d0**intLv
        Zin =R_lim(0,3)+dx(3)!*0.5d0**intLv
        Zout=R_lim(1,3)-dx(3)!*0.5d0**intLv

!!$omp parallel do private(indexP,pp,p0,k,indx,ptmp1,ptmp2,ptmp,p1,octType,newoctN) &
!!$omp&             shared(minP,maxP,Pesh,dxx,R_lim,iLv,iFs,iFe)   
        do indexP=minP,maxP
		
           do i=1,iomp0
              iCp(i,indexP) = 0 
              !iCp2(i,indexP) = 0
           end do

		   
		   pp => Pesh(indexP,iLv)
           if(pp%isort <= 0) cycle
!=================2016/1/23Nakano======================
		if(pp%Vflg == 1 .and. iLv == 0 ) then
			p0 => VMesh(pp%Ioct)
		else 
			p0 => Mesh(pp%Ioct) ! bring up the cell that pp belongs to
		end if
!======================================================

           do k=1,3
              indx(k) = int(pp%R(k)/dxx(k)-p0%rPOS(k)/dxx(k)+2.5d0)-1 ! compute where the particle 
           enddo
           if(indx(1)==0) then 
              ptmp1 => p0%octNb1
           else if(indx(1)==1) then 
              ptmp1 => p0
           else
              ptmp1 => p0%octNb2
			  pp%Vflg = 0
           end if
           if(indx(2)==0) then 
              ptmp2 => ptmp1%octNb3
           else if(indx(2)==1) then 
              ptmp2 => ptmp1
           else
              ptmp2 => ptmp1%octNb4
           end if
           if(indx(3)==0) then 
              ptmp  => ptmp2%octNb5
           else if(indx(3)==1) then 
              ptmp  => ptmp2
           else
              ptmp  => ptmp2%octNb6
           end if
		if(pp%Vflg == 1) then
			pp%isort = 0
			cycle
		end if
           octType = ptmp%octType
           newoctN = ptmp%octN
           p1 => Mesh(newoctN)
		
! How about GMesh 
           pp%Ioct = p1%octN
! Periodic Boundary Condition
           if((p1%iFLG(1).ge.iFs).and.(p1%iFLG(1).le.iFe)) then
            ! --- remove particles outside of boundary ---
            if(ParticleDeleteFlag(1)==1)then
            ! -- for x direction --
!                if(pp%R(1)<Xin) then
!                    pp%isort=0
!                endif
                if(pp%R(1)>Xout) then
                    pp%isort=0
                endif
            else if(ParticleDeleteFlag(2)==1)then
            ! -- for y direction --
!                if(pp%R(2)<Yin) then
!                   pp%isort=0
!                endif
                if(pp%R(2)>Yout) then
                   pp%isort=0
                endif
            else if(ParticleDeleteFlag(3)==1)then
            ! -- for z direction --
!                if(pp%R(3)<Zin) then
!                   pp%isort=0
!                endif
                if(pp%R(3)>Zout) then
                   pp%isort=0
                endif
            end if
            ! -- Periodic boundary for particles --
              do k=1,3
                 if(pp%R(k) < R_lim(0,k)) then
                    pp%R(k) = pp%R(k) + R_lim(1,k)
                 endif
                 if(pp%R(k) > R_lim(1,k)) then
                    pp%R(k) = pp%R(k) - R_lim(1,k)
                 endif
              enddo
           else
              pp%isort = 0 
           endif
!--Pedriodic B.C END
!--erace particles which are outside of overlap regions           
        end do
!!$omp end parallel do


!-------------------Parallel Thread Operation Begin------------------------------
!!$omp parallel do private(j,iomp2,pp,p0) shared(minP,maxP,Pesh,Mesh,iCp)  
        do j=minP,maxP
           pp => Pesh(j,iLv)
           iomp2 = 1
           if(pp%isort.le.0)cycle
!!$          iomp2 = 1+omp_get_thread_num()
           p0 => Mesh(pp%ioct)
           iCp(iomp2,j)=p0%iC(iomp2) ; p0%iC(iomp2)=j
        end do 
!!$omp end parallel do

        if(minindex<maxindex)then
!!$omp parallel do private(index,p0,pptc,i,iC0) shared(minindex,maxindex,Mesh,Pesh)
           do index=minindex,maxindex
              p0   => Mesh(index)
              pptc => p0%ptcl
              do i=1,iomp0
                 iC0 = p0%iC(i)
                 if(iC0.gt.0) then 
                    pptc%prtnxt => Pesh(iC0,iLv)
                    do while(iC0.gt.0)
                       p0%octP = p0%octP+1
                       pptc%prtnxt => Pesh(iC0,iLv)
                       pptc => pptc%prtnxt
                       iC0 = iCp(i,pptc%inum)
                    end do
                 end if
              end do
           end do
!!$omp end parallel do
        endif

        if(Gminindex<Gmaxindex)then
!!$omp parallel do private(index,p0,pptc,i,iC0) shared(Gminindex,Gmaxindex,Mesh,Pesh)
           do index=Gminindex,Gmaxindex
              p0   => Mesh(index)
              pptc => p0%ptcl
              do i=1,iomp0
                 iC0 = p0%iC(i)
                 if(iC0.gt.0) then 
                    pptc%prtnxt => Pesh(iC0,iLv)
                    do while(iC0.gt.0)
                       p0%octP = p0%octP+1
                       pptc%prtnxt => Pesh(iC0,iLv)
                       pptc => pptc%prtnxt
                       iC0 = iCp(i,pptc%inum)
                    end do
                 end if
              end do
           end do
!!$omp end parallel do
        endif
		
		if(iLv == 0)then
		 do index=1,VMeshSize
              p0   => VMesh(index)
              pptc => p0%ptcl
              do i=1,iomp0
                 iC0 = p0%iC(i)
                 if(iC0.gt.0) then 
                    pptc%prtnxt => Pesh(iC0,iLv)
                    do while(iC0.gt.0)
                       p0%octP = p0%octP+1
                       pptc%prtnxt => Pesh(iC0,iLv)
                       pptc => pptc%prtnxt
                       iC0 = iCp(i,pptc%inum)
                    end do
                 end if
              end do
           end do
		end if
		
		
!-------------------Parallel Thread Operation END--------------------------------
        return
end subroutine reconnect_particleP
!
!
!***********************************************************************
      subroutine move_particle(iLv,iFs,iFe)
! +-------------------------------------------------------------------+
! |                                                                   |
! |  Computes force interpolation, advances particles                 |
! |   [C. K. Birdsall and A. B. Langdon, Plasma Physics Via           |
! |    Computer Simulation (Adam-Hilger, 1991)].                      |
! |                                                                   |
! |  Computes current density by density decomposition                |
! |  [T.Zh.Esirkepov, Exact charge conservation scheme for            |
! |   particle-in-cell simulation with an arbitrary form-factor,      |
! |   Computer Physics Communications, v.135 (2001), pp.144-153]      |
! |                                                                   |
! |                                                                   |
! |  2nd order form-factor. Particle occupies 27*dx*dy*dz cells.      |
! |                                                                   |
! |  W(X_{j}-x0) = 3/4 - (|X_{j}-x0|/dx)^2, if |X_{j}-x0|/dx <= 1/2 ; |
! |   or 1/2*( 3/2-|X_{j}-x0|/dx )^2, if 1/2 < |X_{j}-x0|/dx <= 3/2 ; |
! |   or 0, if |X_{j}-x0|/dx > 3/2 .                                  |
! |                                                                   |
! |  dx=dy=dz is strongly recommended!                                |
! |                                                                   |
! |                                                                   |
! |                                                                   |
! |   Argument parameters:                                            |
! |       i2step  : iLv+1 levels time step                            |
! |       istep   : iLv   levels time step                            |
! |       iLv     : target hierarchical level                         |
! |                                                                   |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************
        use message_passing_interface
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        use param
!-----------------------------------------------------------------------
        implicit none
        integer(kind=4) :: iLv,iFs,iFe,iomp2
        real(kind=8)    :: dvlevel, deltaT,dvlevel3
!---------------
        integer(kind=4)      :: index,minindex,maxindex,minindex0,maxindex0
        integer(kind=4)      :: Nin,itotal
        real(kind=4)         :: Rinin
        real(kind=8)         :: R1(10)
        integer(kind=4)      :: i,j,k,ii,jj,kk,kp,indexP,particle,ipp,imf
        type(oct), pointer   :: p0,p1
        integer(kind=4)      :: ip(-1:1,-1:1,-1:1),jp(2,2,2)
        !integer(kind=4)      :: ip2(-1:1,-1:1,-1:1) ! ip2 plays a role in differentiating BMesh, Mesh and GMesh
!  It is not efficient to do in this way but we get to write something that works ASAP.
        real(kind=8)         :: Ffld(6,6,6,6)
        integer(kind=4)      :: j0(10,3),k1,m1,n1,k0,m0,n0
        type(prtcl), pointer :: pptc
        real(kind=8)         :: R(23),R0(-2:20,3),dipoleMagnetic(1:3),rx,ry,rz,rr,mr
        integer(kind=4)      :: j1(7)
        real(kind=8)         :: Cfld(3,8,8,8,IonSorts)
        type(oct), pointer   :: p2,p3,p4,p5,p6,p7,p8
!------------------------
!$        include 'omp_lib.h'
!
! == BEGIN: ============================================================

        !if(debugMode>=1)print*,'entering move particles:rank, iLv',rank, ilv
        if((maxID(1,iLv).eq.minID(1,iLv)) .and. (maxID(3,iLv).eq.minID(3,iLv))) return 
!
        if(maxID(1,iLv)>minID(1,iLv))then
!$omp parallel do private(index,p0,k) shared(minID,maxID,iFs,iFe,Mesh,iLv)
! Move old data into F(10,11,12)
           do index=MinID(1,iLv),MaxID(1,iLv)
              p0 => Mesh(index)
              if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then 
                 do k=7,9
                    p0%F(k+3) = p0%F(k)
                    p0%F(k  ) = ZERO
                 enddo
                 do k=13,18
                    p0%F(k)=ZERO
                 end do
                 do k=1,6*iomp0
                    p0%C(k  ) = ZERO
                 end do
              endif
           end do
!$omp end parallel do
        endif
!
        if(maxID(3,iLv)>minID(3,iLv))then
!$omp parallel do private(index,p0,k) shared(MinID,MaxID,iLv,Mesh,iFs,iFe)
           do index=MinID(3,iLv),MaxID(3,iLv)
              p0 => Mesh(index)
              if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then 
                 do k=7,9
                    p0%F(k+3) = p0%F(k)
                    p0%F(k  ) = ZERO
                 enddo
                 do k=13,18
                    p0%F(k)=ZERO
                 end do
                 do k=1,6*iomp0
                    p0%C(k  ) = ZERO
                 end do
              endif
           end do
!$omp end parallel do
        endif

!=======================================================
        R1=0.d0
! -- Set spacing parameters --
        dvlevel = 0.5d0**(iLv)
        dvlevel3=(dvlevel**3)
        deltaT  = dvlevel * dt
        Nin   = 2**(LvMax-iLv)
        Rinin = 1.d0/real(2*Nin)
        PICnumber(1) = npart_per_cell
        PICnumber(2) = npart_per_cell
!
        !R1( 1) = 0.5d0*QM*deltaT*Rcharge(1)/Rmass(1)
        R1( 1) = 0.5d0*QM*deltaT*Rcharge(1)/Rmass(1)*omega
        R1( 2) = Rcharge(1)/dvlevel3*(-1.d0)/dble(PICnumber(1))
        R1( 3) = ONE/(dvlevel*dx(1))
!                      if(rank==1)print*,'R1(3)',R1(3), dvlevel, dx(1)
        !R1( 4) = 0.5d0*QM*deltaT*Rcharge(2)/Rmass(2)
        R1( 4) = 0.5d0*QM*deltaT*Rcharge(2)/Rmass(2)*omega
        R1( 5) = Rcharge(2)/dvlevel3*(-1.d0)/dble(PICnumber(2))
        R1( 6) = ONE/(dvlevel*dx(2))
        R1( 7)=  Rinin
        R1( 8)=  deltaT
        R1( 9)=-dx(1)*dx(2)*BD
        R1(10) = ONE/(dvlevel*dx(3))
! -- parallel loop

        maxindex = MaxID(1,iLv)
        minindex = MinID(1,iLv)
        maxindex0= MaxID(1,iLv-1)
        minindex0= MinID(1,iLv-1)

!$omp parallel do private(iomp2,index,p0,ip,i,j,k,ii,jj,kk,Ffld,Cfld,jp    &
!$omp&        ,R0,indexP,pptc,kp,R,k0,m0,n0,j0,j1,p1,p2,p3,p4,p5,p6,p7,p8,itotal,k1,m1,n1 & 
!$omp&        ,rx,ry,rz,rr,mr,dipoleMagnetic,particle) &
!$omp&             shared(minindex0,maxindex0,Mesh,iFs,iFe,Pesh,R1)        
!!$omp&           schedule(dynamic,2)

        do index = minindex0, maxindex0

           p0 => Mesh(index)
            iomp2 = 0
!$          iomp2 = 3*(omp_get_thread_num())
! ENTERED BY TATSUKI

! TATSUKI END

          R0(-1,1)=R1( 1) ; R0( 0,1)=R1( 2) ; R0( 1,1)=R1( 3)
          R0(-1,2)=R1( 4) ; R0( 0,2)=R1( 5) ; R0( 1,2)=R1( 6)
          R0(-2,3)=R1( 7) ; R0(-1,3)=R1( 8) ; R0( 0,3)=R1( 9) ; R0( 1,3)=R1(10)

!          print*,'I am ',rank,'Current p0 is :octType', p0%octType, 'index',index 
!          print*,'[',(p0%iPOS+2**(LvMax-(iLv-1)))/(2**(LvMax-(iLv-1)+1)),']'
          
!print*, 'came here' 
          if((p0%iFLG(1)<iFs+4).or.(p0%iFLG(1)>iFe+4)) cycle
!            if(p0%prc_bndry/=0) cycle
        ! -- for vacuum background cases

!-- added by TATSUKI for multi-process run. We do not move particles if they are at the boundaries

          itotal = p0%octCh1%octP+p0%octCh2%octP+p0%octCh3%octP+p0%octCh4%octP &
                  +p0%octCh5%octP+p0%octCh6%octP+p0%octCh7%octP+p0%octCh8%octP

          if(itotal.le.0)then
             cycle
          endif

          ip(-1,-1,-1) = p0%octNb5%octNb3%octNb1%octN
          ip( 0,-1,-1) = p0%octNb5%octNb3       %octN
          ip( 1,-1,-1) = p0%octNb5%octNb3%octNb2%octN
          ip(-1, 0,-1) = p0%octNb5       %octNb1%octN
          ip( 0, 0,-1) = p0%octNb5              %octN
          ip( 1, 0,-1) = p0%octNb5       %octNb2%octN
          ip(-1, 1,-1) = p0%octNb5%octNb4%octNb1%octN
          ip( 0, 1,-1) = p0%octNb5%octNb4       %octN
          ip( 1, 1,-1) = p0%octNb5%octNb4%octNb2%octN

          ip(-1,-1, 0) = p0       %octNb3%octNb1%octN
          ip( 0,-1, 0) = p0       %octNb3       %octN
          ip( 1,-1, 0) = p0       %octNb3%octNb2%octN
          ip(-1, 0, 0) = p0              %octNb1%octN
          ip( 0, 0, 0) = p0                     %octN
          ip( 1, 0, 0) = p0              %octNb2%octN
          ip(-1, 1, 0) = p0       %octNb4%octNb1%octN
          ip( 0, 1, 0) = p0       %octNb4       %octN
          ip( 1, 1, 0) = p0       %octNb4%octNb2%octN

          ip(-1,-1, 1) = p0%octNb6%octNb3%octNb1%octN
          ip( 0,-1, 1) = p0%octNb6%octNb3       %octN
          ip( 1,-1, 1) = p0%octNb6%octNb3%octNb2%octN
          ip(-1, 0, 1) = p0%octNb6       %octNb1%octN
          ip( 0, 0, 1) = p0%octNb6              %octN
          ip( 1, 0, 1) = p0%octNb6       %octNb2%octN
          ip(-1, 1, 1) = p0%octNb6%octNb4%octNb1%octN
          ip( 0, 1, 1) = p0%octNb6%octNb4       %octN
          ip( 1, 1, 1) = p0%octNb6%octNb4%octNb2%octN

! Copying Fields onto Ffld array

          do k=-1,1
             do j=-1,1
                do i=-1,1
                   ! This is a dumb but the easiest way to select a pointer
                  
                   p1 => Mesh(ip(i,j,k))

                   ii = 3+i*2
                   jj = 3+j*2
                   kk = 3+k*2
                   
                   if (dipoleFlag == 1) then
                      rx=p1%octCh1%rPOS(1)-Dpos(1)
                      ry=p1%octCh1%rPOS(2)-Dpos(2)
                      rz=p1%octCh1%rPOS(3)-Dpos(3)
                      rr=sqrt(rx**2+ry**2+rz**2)
                      mr=mx*rx+my*ry+mz*rz
                      dipoleMagnetic(1) = (-mx/(rr**3)+3*mr*rx/(rr**5))/(4*PI)
                      dipoleMagnetic(2) = (-my/(rr**3)+3*mr*ry/(rr**5))/(4*PI)
                      dipoleMagnetic(3) = (-mz/(rr**3)+3*mr*rz/(rr**5))/(4*PI)
                   endif
                   
                   Ffld( 1,ii  ,jj,kk)=p1%octCh1%F(1)
                   Ffld( 2,ii  ,jj,kk)=p1%octCh1%F(2)
                   Ffld( 3,ii  ,jj,kk)=p1%octCh1%F(3)
                   Ffld( 4,ii  ,jj,kk)=p1%octCh1%F(4)+dipoleMagnetic(1)
                   Ffld( 5,ii  ,jj,kk)=p1%octCh1%F(5)+dipoleMagnetic(2)
                   Ffld( 6,ii  ,jj,kk)=p1%octCh1%F(6)+dipoleMagnetic(3)
                   
                   if (dipoleFlag == 1) then
                      rx=p1%octCh2%rPOS(1)-Dpos(1)
                      ry=p1%octCh2%rPOS(2)-Dpos(2)
                      rz=p1%octCh2%rPOS(3)-Dpos(3)
                      rr=sqrt(rx**2+ry**2+rz**2)
                      mr=mx*rx+my*ry+mz*rz
                      dipoleMagnetic(1) = (-mx/(rr**3)+3*mr*rx/(rr**5))/(4*PI)
                      dipoleMagnetic(2) = (-my/(rr**3)+3*mr*ry/(rr**5))/(4*PI)
                      dipoleMagnetic(3) = (-mz/(rr**3)+3*mr*rz/(rr**5))/(4*PI)
                   endif
                      
                   Ffld( 1,ii+1,jj,kk)=p1%octCh2%F(1)
                   Ffld( 2,ii+1,jj,kk)=p1%octCh2%F(2)
                   Ffld( 3,ii+1,jj,kk)=p1%octCh2%F(3)
                   Ffld( 4,ii+1,jj,kk)=p1%octCh2%F(4)+dipoleMagnetic(1)
                   Ffld( 5,ii+1,jj,kk)=p1%octCh2%F(5)+dipoleMagnetic(2)
                   Ffld( 6,ii+1,jj,kk)=p1%octCh2%F(6)+dipoleMagnetic(3)
                   
                   jj=jj+1
                   
                   if (dipoleFlag == 1) then
                      rx=p1%octCh3%rPOS(1)-Dpos(1)
                      ry=p1%octCh3%rPOS(2)-Dpos(2)
                      rz=p1%octCh3%rPOS(3)-Dpos(3)
                      rr=sqrt(rx**2+ry**2+rz**2)
                      mr=mx*rx+my*ry+mz*rz
                      dipoleMagnetic(1) = (-mx/(rr**3)+3*mr*rx/(rr**5))/(4*PI)
                      dipoleMagnetic(2) = (-my/(rr**3)+3*mr*ry/(rr**5))/(4*PI)
                      dipoleMagnetic(3) = (-mz/(rr**3)+3*mr*rz/(rr**5))/(4*PI)
                   endif
                   Ffld( 1,ii  ,jj,kk)=p1%octCh3%F(1)
                   Ffld( 2,ii  ,jj,kk)=p1%octCh3%F(2)
                   Ffld( 3,ii  ,jj,kk)=p1%octCh3%F(3)
                   Ffld( 4,ii  ,jj,kk)=p1%octCh3%F(4)+dipoleMagnetic(1)
                   Ffld( 5,ii  ,jj,kk)=p1%octCh3%F(5)+dipoleMagnetic(2)
                   Ffld( 6,ii  ,jj,kk)=p1%octCh3%F(6)+dipoleMagnetic(3)
                   
                   if (dipoleFlag == 1) then
                      rx=p1%octCh4%rPOS(1)-Dpos(1)
                      ry=p1%octCh4%rPOS(2)-Dpos(2)
                      rz=p1%octCh4%rPOS(3)-Dpos(3)
                      rr=sqrt(rx**2+ry**2+rz**2)
                      mr=mx*rx+my*ry+mz*rz
                      dipoleMagnetic(1) = (-mx/(rr**3)+3*mr*rx/(rr**5))/(4*PI)
                      dipoleMagnetic(2) = (-my/(rr**3)+3*mr*ry/(rr**5))/(4*PI)
                      dipoleMagnetic(3) = (-mz/(rr**3)+3*mr*rz/(rr**5))/(4*PI)
                   endif
                   
                   Ffld( 1,ii+1,jj,kk)=p1%octCh4%F(1)
                   Ffld( 2,ii+1,jj,kk)=p1%octCh4%F(2)
                   Ffld( 3,ii+1,jj,kk)=p1%octCh4%F(3)
                   Ffld( 4,ii+1,jj,kk)=p1%octCh4%F(4)+dipoleMagnetic(1)
                   Ffld( 5,ii+1,jj,kk)=p1%octCh4%F(5)+dipoleMagnetic(2)
                   Ffld( 6,ii+1,jj,kk)=p1%octCh4%F(6)+dipoleMagnetic(3)
                end do
             end do
          end do

!          print*, 'came here4' , rank,index
! Copying Fields onto Ffld array
          do k=-1,1
             do j=-1,1
                do i=-1,1

                   p1 => Mesh(ip(i,j,k))

                   ii = 3+i*2
                   jj = 3+j*2
                   kk = 4+k*2
                   
                   if (dipoleFlag == 1) then
                      rx=p1%octCh5%rPOS(1)-Dpos(1)
                      ry=p1%octCh5%rPOS(2)-Dpos(2)
                      rz=p1%octCh5%rPOS(3)-Dpos(3)
                      rr=sqrt(rx**2+ry**2+rz**2)
                      mr=mx*rx+my*ry+mz*rz
                      dipoleMagnetic(1) = (-mx/(rr**3)+3*mr*rx/(rr**5))/(4*PI)
                      dipoleMagnetic(2) = (-my/(rr**3)+3*mr*ry/(rr**5))/(4*PI)
                      dipoleMagnetic(3) = (-mz/(rr**3)+3*mr*rz/(rr**5))/(4*PI)
                   endif
                   
                   Ffld( 1,ii  ,jj,kk)=p1%octCh5%F(1)
                   Ffld( 2,ii  ,jj,kk)=p1%octCh5%F(2)
                   Ffld( 3,ii  ,jj,kk)=p1%octCh5%F(3)
                   Ffld( 4,ii  ,jj,kk)=p1%octCh5%F(4)+dipoleMagnetic(1)
                   Ffld( 5,ii  ,jj,kk)=p1%octCh5%F(5)+dipoleMagnetic(2)
                   Ffld( 6,ii  ,jj,kk)=p1%octCh5%F(6)+dipoleMagnetic(3)
                   
                   if (dipoleFlag == 1) then
                      rx=p1%octCh6%rPOS(1)-Dpos(1)
                      ry=p1%octCh6%rPOS(2)-Dpos(2)
                      rz=p1%octCh6%rPOS(3)-Dpos(3)
                      rr=sqrt(rx**2+ry**2+rz**2)
                      mr=mx*rx+my*ry+mz*rz
                      dipoleMagnetic(1) = (-mx/(rr**3)+3*mr*rx/(rr**5))/(4*PI)
                      dipoleMagnetic(2) = (-my/(rr**3)+3*mr*ry/(rr**5))/(4*PI)
                      dipoleMagnetic(3) = (-mz/(rr**3)+3*mr*rz/(rr**5))/(4*PI)
                   endif
                   
                   Ffld( 1,ii+1,jj,kk)=p1%octCh6%F(1)
                   Ffld( 2,ii+1,jj,kk)=p1%octCh6%F(2)
                   Ffld( 3,ii+1,jj,kk)=p1%octCh6%F(3)
                   Ffld( 4,ii+1,jj,kk)=p1%octCh6%F(4)+dipoleMagnetic(1)
                   Ffld( 5,ii+1,jj,kk)=p1%octCh6%F(5)+dipoleMagnetic(2)
                   Ffld( 6,ii+1,jj,kk)=p1%octCh6%F(6)+dipoleMagnetic(3)
                   
                   jj=jj+1
                   
                   if (dipoleFlag == 1) then
                      rx=p1%octCh7%rPOS(1)-Dpos(1)
                      ry=p1%octCh7%rPOS(2)-Dpos(2)
                      rz=p1%octCh7%rPOS(3)-Dpos(3)
                      rr=sqrt(rx**2+ry**2+rz**2)
                      mr=mx*rx+my*ry+mz*rz
                      dipoleMagnetic(1) = (-mx/(rr**3)+3*mr*rx/(rr**5))/(4*PI)
                      dipoleMagnetic(2) = (-my/(rr**3)+3*mr*ry/(rr**5))/(4*PI)
                      dipoleMagnetic(3) = (-mz/(rr**3)+3*mr*rz/(rr**5))/(4*PI)
                   endif
                   
                   Ffld( 1,ii  ,jj,kk)=p1%octCh7%F(1)
                   Ffld( 2,ii  ,jj,kk)=p1%octCh7%F(2)
                   Ffld( 3,ii  ,jj,kk)=p1%octCh7%F(3)
                   Ffld( 4,ii  ,jj,kk)=p1%octCh7%F(4)+dipoleMagnetic(1)
                   Ffld( 5,ii  ,jj,kk)=p1%octCh7%F(5)+dipoleMagnetic(2)
                   Ffld( 6,ii  ,jj,kk)=p1%octCh7%F(6)+dipoleMagnetic(3)
                   
                   if (dipoleFlag == 1) then
                      rx=p1%octCh8%rPOS(1)-Dpos(1)
                      ry=p1%octCh8%rPOS(2)-Dpos(2)
                      rz=p1%octCh8%rPOS(3)-Dpos(3)
                      rr=sqrt(rx**2+ry**2+rz**2)
                      mr=mx*rx+my*ry+mz*rz
                      dipoleMagnetic(1) = (-mx/(rr**3)+3*mr*rx/(rr**5))/(4*PI)
                      dipoleMagnetic(2) = (-my/(rr**3)+3*mr*ry/(rr**5))/(4*PI)
                      dipoleMagnetic(3) = (-mz/(rr**3)+3*mr*rz/(rr**5))/(4*PI)
                   endif
                   
                   Ffld( 1,ii+1,jj,kk)=p1%octCh8%F(1)
                   Ffld( 2,ii+1,jj,kk)=p1%octCh8%F(2)
                   Ffld( 3,ii+1,jj,kk)=p1%octCh8%F(3)
                   Ffld( 4,ii+1,jj,kk)=p1%octCh8%F(4)+dipoleMagnetic(1)
                   Ffld( 5,ii+1,jj,kk)=p1%octCh8%F(5)+dipoleMagnetic(2)
                   Ffld( 6,ii+1,jj,kk)=p1%octCh8%F(6)+dipoleMagnetic(3)
                end do
             end do
          end do
  
!*soption PREFETCH_ZERO(Cfld)
          Cfld(:,:,:,:,:) = 0d0
!
          jp(1,1,1)=p0%octCh1%octN !NOTE: p0 is moving along either BMesh or Mesh so jp are ALWAYS 
          jp(2,1,1)=p0%octCh2%octN ! pointing at Mesh (i.e., NOT BMesh,or GMesh) 
          jp(1,2,1)=p0%octCh3%octN
          jp(2,2,1)=p0%octCh4%octN
          jp(1,1,2)=p0%octCh5%octN
          jp(2,1,2)=p0%octCh6%octN
          jp(1,2,2)=p0%octCh7%octN
          jp(2,2,2)=p0%octCh8%octN
!

!call Display_Grid(0)

          do kk=1,2
             do jj=1,2
                do ii=1,2

                   
                   p1 => Mesh(jp(ii,jj,kk))
!
! -- Get cell index --
                   R0(2,1) = real(p1%iPOS(1))*R0(-2,3) + 0.5d0
                   R0(2,2) = real(p1%iPOS(2))*R0(-2,3) + 0.5d0
                   R0(2,3) = real(p1%iPOS(3))*R0(-2,3) + 0.5d0
!
! -- pointer loop for particle --
!
                   indexP =  p1%octP 
                   if(indexP==0) cycle

                   pptc   => p1%ptcl%prtnxt

!    Pushing particles INSIDE Mesh 
                   do kp=1,indexP
!            print*, 'I am ',rank, pptc % ioct
!
                      R( 1)= pptc%R(1)
                      R( 2)= pptc%R(2)
                      R( 3)= pptc%R(3)
                      R( 4)= pptc%R(4)
                      R( 5)= pptc%R(5)
                      R( 6)= pptc%R(6)

! --- Line Dipole --- kozawa
! --------------------------------------------------------------------------
!            R( 7)= R0(-2,3)+R(1)
!            R( 8)= R0(-1,3)+R(2)
!            R( 9)= R0( 0,3)/((R(7)**2+R(8)**2)**2)
!            R(10)=(   R(7)*R(8)   )*R(9)
!            R(11)=(R(8)**2-R(7)**2)*R(9)*0.5d0
!            R(12)= 0.d0
! --------------------------------------------------------------------------
                      R( 7)= 0.d0
                      R( 8)= 0.d0
                      R( 9)= 0.d0
                      R(10)= 0.d0
                      R(11)= 0.d0
                      R(12)= 0.d0
                      R(13)= R(4)
                      R(14)= R(5)
                      R(15)= R(6)
                      R(16)= R0(-1,pptc%Isort)
                      R(23)= R0( 0,pptc%Isort)
                      particle= pptc%Isort
!
! =====================================================
!    Form-factor calculation before particle movement
!    R_lim(0,k) should be 0
! =====================================================
!
! -- compute 1D form-factors --


                      do k=1,3
                         R0( 3,k) = R(k)*R0(1,k) + 0.5d0
! ---  base grid (+1) ---
                         R0( 4,k) = R0(3,k) + 0.5d0                  
                         j0( 1,k) = NINT( R0(4,k) )        ! j=1 corresponds to x=0
                         j0( 2,k) = j0(1,k) - R0(2,k)
                         R0( 4,k) = R0(4,k) - dble(j0(1,k))
                         R0( 5,k) = R0(4,k)*R0(4,k)
                         R0( 6,k) = 0.125d0 + 0.5d0*R0(5,k)
                         R0( 4,k) =           0.5d0*R0(4,k)
                         R0( 7,k) = R0(6,k)       - R0(4,k)      ! 0.5*(0.5-h1)**2 ; W1(-1,k)
                         R0( 8,k) =  0.75d0       - R0(5,k)   ! 0.75 - h1**2    ; W1( 0,k)
                         R0( 9,k) = R0(6,k)       + R0(4,k)      !                 ; W1( 1,k)
                         R0(10,k) = 0.d0                   !                 ; W1( 2,k)
! --- shifted grid (+1/2) ---
                         R0( 4,k) = R0(3,k)
                         j0( 3,k) = NINT( R0(4,k) )        ! j=1 corresponds to x=0
                         j0( 4,k) = j0(3,k) - R0(2,k)
                         R0( 4,k) = R0(4,k) - dble(j0(3,k))
                         R0( 5,k) = R0(4,k)*R0(4,k)
                         R0( 6,k) = 0.125d0 + 0.5d0*R0(5,k)
                         R0( 4,k) =           0.5d0*R0(4,k)
                         R0(11,k) = R0(6,k)       - R0(4,k)           ! 0.5*(0.5-h1)**2 ; W0(-1,k)
                         R0(12,k) =  0.75d0       - R0(5,k)   ! 0.75 - h1**2    ; W0( 0,k)
                         R0(13,k) = R0(6,k)       + R0(4,k)            ! 0.5*(0.5+h1)**2 ; W0( 1,k)
                         R0(14,k) = 0.d0              !                 ; W0( 2,k)
                         
                      enddo

! =====================================================
!    Force calculation
! =====================================================
                      
                      k1=j0(2,1)+ii+2
                      m1=j0(2,2)+jj+2
                      n1=j0(2,3)+kk+2
                      k0=j0(4,1)+ii+2
                      m0=j0(4,2)+jj+2
                      n0=j0(4,3)+kk+2


                      R( 7)=R( 7) +(( Ffld(1,k1-1,m0-1,n0-1)*R0( 7,1)                     &
                                    + Ffld(1,k1  ,m0-1,n0-1)*R0( 8,1)                     &
                                    + Ffld(1,k1+1,m0-1,n0-1)*R0( 9,1))*R0(11,2)           &
                                   +( Ffld(1,k1-1,m0  ,n0-1)*R0( 7,1)                     &
                                    + Ffld(1,k1  ,m0  ,n0-1)*R0( 8,1)                     &
                                    + Ffld(1,k1+1,m0  ,n0-1)*R0( 9,1))*R0(12,2)           &
                                   +( Ffld(1,k1-1,m0+1,n0-1)*R0( 7,1)                     &
                                    + Ffld(1,k1  ,m0+1,n0-1)*R0( 8,1)                     &
                                    + Ffld(1,k1+1,m0+1,n0-1)*R0( 9,1))*R0(13,2))*R0(11,3) &
                                  +(( Ffld(1,k1-1,m0-1,n0  )*R0( 7,1)                     &
                                    + Ffld(1,k1  ,m0-1,n0  )*R0( 8,1)                     &
                                    + Ffld(1,k1+1,m0-1,n0  )*R0( 9,1))*R0(11,2)           &
                                   +( Ffld(1,k1-1,m0  ,n0  )*R0( 7,1)                     &
                                    + Ffld(1,k1  ,m0  ,n0  )*R0( 8,1)                     &
                                    + Ffld(1,k1+1,m0  ,n0  )*R0( 9,1))*R0(12,2)           &
                                   +( Ffld(1,k1-1,m0+1,n0  )*R0( 7,1)                     &
                                    + Ffld(1,k1  ,m0+1,n0  )*R0( 8,1)                     &
                                    + Ffld(1,k1+1,m0+1,n0  )*R0( 9,1))*R0(13,2))*R0(12,3) &
                                  +(( Ffld(1,k1-1,m0-1,n0+1)*R0( 7,1)                     &
                                    + Ffld(1,k1  ,m0-1,n0+1)*R0( 8,1)                     &
                                    + Ffld(1,k1+1,m0-1,n0+1)*R0( 9,1))*R0(11,2)           &
                                   +( Ffld(1,k1-1,m0  ,n0+1)*R0( 7,1)                     &
                                    + Ffld(1,k1  ,m0  ,n0+1)*R0( 8,1)                     &
                                    + Ffld(1,k1+1,m0  ,n0+1)*R0( 9,1))*R0(12,2)           &
                                   +( Ffld(1,k1-1,m0+1,n0+1)*R0( 7,1)                     &
                                    + Ffld(1,k1  ,m0+1,n0+1)*R0( 8,1)                     &
                                    + Ffld(1,k1+1,m0+1,n0+1)*R0( 9,1))*R0(13,2))*R0(13,3) 

!-----------------------------------------------------------------------------------------
                      R( 8)=R( 8) +(( Ffld(2,k0-1,m1-1,n0-1)*R0(11,1)                     &
                                    + Ffld(2,k0  ,m1-1,n0-1)*R0(12,1)                     &
                                    + Ffld(2,k0+1,m1-1,n0-1)*R0(13,1))*R0( 7,2)           &
                                   +( Ffld(2,k0-1,m1  ,n0-1)*R0(11,1)                     &
                                    + Ffld(2,k0  ,m1  ,n0-1)*R0(12,1)                     &
                                    + Ffld(2,k0+1,m1  ,n0-1)*R0(13,1))*R0( 8,2)           &
                                   +( Ffld(2,k0-1,m1+1,n0-1)*R0(11,1)                     &
                                    + Ffld(2,k0  ,m1+1,n0-1)*R0(12,1)                     &
                                    + Ffld(2,k0+1,m1+1,n0-1)*R0(13,1))*R0( 9,2))*R0(11,3) &
                                  +(( Ffld(2,k0-1,m1-1,n0  )*R0(11,1)                     &
                                    + Ffld(2,k0  ,m1-1,n0  )*R0(12,1)                     &
                                    + Ffld(2,k0+1,m1-1,n0  )*R0(13,1))*R0( 7,2)           &
                                   +( Ffld(2,k0-1,m1  ,n0  )*R0(11,1)                     &
                                    + Ffld(2,k0  ,m1  ,n0  )*R0(12,1)                     &
                                    + Ffld(2,k0+1,m1  ,n0  )*R0(13,1))*R0( 8,2)           &
                                   +( Ffld(2,k0-1,m1+1,n0  )*R0(11,1)                     &
                                    + Ffld(2,k0  ,m1+1,n0  )*R0(12,1)                     &
                                    + Ffld(2,k0+1,m1+1,n0  )*R0(13,1))*R0( 9,2))*R0(12,3) &
                                  +(( Ffld(2,k0-1,m1-1,n0+1)*R0(11,1)                     &
                                    + Ffld(2,k0  ,m1-1,n0+1)*R0(12,1)                     &
                                    + Ffld(2,k0+1,m1-1,n0+1)*R0(13,1))*R0( 7,2)           &
                                   +( Ffld(2,k0-1,m1  ,n0+1)*R0(11,1)                     &
                                    + Ffld(2,k0  ,m1  ,n0+1)*R0(12,1)                     &
                                    + Ffld(2,k0+1,m1  ,n0+1)*R0(13,1))*R0( 8,2)           &
                                   +( Ffld(2,k0-1,m1+1,n0+1)*R0(11,1)                     &
                                    + Ffld(2,k0  ,m1+1,n0+1)*R0(12,1)                     &
                                    + Ffld(2,k0+1,m1+1,n0+1)*R0(13,1))*R0( 9,2))*R0(13,3) 
!-----------------------------------------------------------------------------------------
                      R( 9)=R( 9) +(( Ffld(3,k0-1,m0-1,n1-1)*R0(11,1)                     &
                                    + Ffld(3,k0  ,m0-1,n1-1)*R0(12,1)                     &
                                    + Ffld(3,k0+1,m0-1,n1-1)*R0(13,1))*R0(11,2)           &
                                   +( Ffld(3,k0-1,m0  ,n1-1)*R0(11,1)                     &
                                    + Ffld(3,k0  ,m0  ,n1-1)*R0(12,1)                     &
                                    + Ffld(3,k0+1,m0  ,n1-1)*R0(13,1))*R0(12,2)           &
                                   +( Ffld(3,k0-1,m0+1,n1-1)*R0(11,1)                     &
                                    + Ffld(3,k0  ,m0+1,n1-1)*R0(12,1)                     &
                                    + Ffld(3,k0+1,m0+1,n1-1)*R0(13,1))*R0(13,2))*R0( 7,3) &
                                  +(( Ffld(3,k0-1,m0-1,n1  )*R0(11,1)                     &
                                    + Ffld(3,k0  ,m0-1,n1  )*R0(12,1)                     &
                                    + Ffld(3,k0+1,m0-1,n1  )*R0(13,1))*R0(11,2)           &
                                   +( Ffld(3,k0-1,m0  ,n1  )*R0(11,1)                     &
                                    + Ffld(3,k0  ,m0  ,n1  )*R0(12,1)                     &
                                    + Ffld(3,k0+1,m0  ,n1  )*R0(13,1))*R0(12,2)           &
                                   +( Ffld(3,k0-1,m0+1,n1  )*R0(11,1)                     &
                                    + Ffld(3,k0  ,m0+1,n1  )*R0(12,1)                     &
                                    + Ffld(3,k0+1,m0+1,n1  )*R0(13,1))*R0(13,2))*R0( 8,3) &
                                  +(( Ffld(3,k0-1,m0-1,n1+1)*R0(11,1)                     &
                                    + Ffld(3,k0  ,m0-1,n1+1)*R0(12,1)                     &
                                    + Ffld(3,k0+1,m0-1,n1+1)*R0(13,1))*R0(11,2)           &
                                   +( Ffld(3,k0-1,m0  ,n1+1)*R0(11,1)                     &
                                    + Ffld(3,k0  ,m0  ,n1+1)*R0(12,1)                     &
                                    + Ffld(3,k0+1,m0  ,n1+1)*R0(13,1))*R0(12,2)           &
                                   +( Ffld(3,k0-1,m0+1,n1+1)*R0(11,1)                     &
                                    + Ffld(3,k0  ,m0+1,n1+1)*R0(12,1)                     &
                                    + Ffld(3,k0+1,m0+1,n1+1)*R0(13,1))*R0(13,2))*R0( 9,3) 
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
                      R(10)=R(10) +(( Ffld(4,k0-1,m1-1,n1-1)*R0(11,1)                     &
                                    + Ffld(4,k0  ,m1-1,n1-1)*R0(12,1)                     &
                                    + Ffld(4,k0+1,m1-1,n1-1)*R0(13,1))*R0( 7,2)           &
                                   +( Ffld(4,k0-1,m1  ,n1-1)*R0(11,1)                     &
                                    + Ffld(4,k0  ,m1  ,n1-1)*R0(12,1)                     &
                                    + Ffld(4,k0+1,m1  ,n1-1)*R0(13,1))*R0( 8,2)           &
                                   +( Ffld(4,k0-1,m1+1,n1-1)*R0(11,1)                     &
                                    + Ffld(4,k0  ,m1+1,n1-1)*R0(12,1)                     &
                                    + Ffld(4,k0+1,m1+1,n1-1)*R0(13,1))*R0( 9,2))*R0( 7,3) &
                                  +(( Ffld(4,k0-1,m1-1,n1  )*R0(11,1)                     &
                                    + Ffld(4,k0  ,m1-1,n1  )*R0(12,1)                     &
                                    + Ffld(4,k0+1,m1-1,n1  )*R0(13,1))*R0( 7,2)           &
                                   +( Ffld(4,k0-1,m1  ,n1  )*R0(11,1)                     &
                                    + Ffld(4,k0  ,m1  ,n1  )*R0(12,1)                     &
                                    + Ffld(4,k0+1,m1  ,n1  )*R0(13,1))*R0( 8,2)           &
                                   +( Ffld(4,k0-1,m1+1,n1  )*R0(11,1)                     &
                                    + Ffld(4,k0  ,m1+1,n1  )*R0(12,1)                     &
                                    + Ffld(4,k0+1,m1+1,n1  )*R0(13,1))*R0( 9,2))*R0( 8,3) &
                                  +(( Ffld(4,k0-1,m1-1,n1+1)*R0(11,1)                     &
                                    + Ffld(4,k0  ,m1-1,n1+1)*R0(12,1)                     &
                                    + Ffld(4,k0+1,m1-1,n1+1)*R0(13,1))*R0( 7,2)           &
                                   +( Ffld(4,k0-1,m1  ,n1+1)*R0(11,1)                     &
                                    + Ffld(4,k0  ,m1  ,n1+1)*R0(12,1)                     &
                                    + Ffld(4,k0+1,m1  ,n1+1)*R0(13,1))*R0( 8,2)           &
                                   +( Ffld(4,k0-1,m1+1,n1+1)*R0(11,1)                     &
                                    + Ffld(4,k0  ,m1+1,n1+1)*R0(12,1)                     &
                                    + Ffld(4,k0+1,m1+1,n1+1)*R0(13,1))*R0( 9,2))*R0( 9,3) 
!-----------------------------------------------------------------------------------------
                      R(11)=R(11) +(( Ffld(5,k1-1,m0-1,n1-1)*R0( 7,1)                     &
                                    + Ffld(5,k1  ,m0-1,n1-1)*R0( 8,1)                     &
                                    + Ffld(5,k1+1,m0-1,n1-1)*R0( 9,1))*R0(11,2)           &
                                   +( Ffld(5,k1-1,m0  ,n1-1)*R0( 7,1)                     &
                                    + Ffld(5,k1  ,m0  ,n1-1)*R0( 8,1)                     &
                                    + Ffld(5,k1+1,m0  ,n1-1)*R0( 9,1))*R0(12,2)           &
                                   +( Ffld(5,k1-1,m0+1,n1-1)*R0( 7,1)                     &
                                    + Ffld(5,k1  ,m0+1,n1-1)*R0( 8,1)                     &
                                    + Ffld(5,k1+1,m0+1,n1-1)*R0( 9,1))*R0(13,2))*R0( 7,3) &
                                  +(( Ffld(5,k1-1,m0-1,n1  )*R0( 7,1)                     &
                                    + Ffld(5,k1  ,m0-1,n1  )*R0( 8,1)                     &
                                    + Ffld(5,k1+1,m0-1,n1  )*R0( 9,1))*R0(11,2)           &
                                   +( Ffld(5,k1-1,m0  ,n1  )*R0( 7,1)                     &
                                    + Ffld(5,k1  ,m0  ,n1  )*R0( 8,1)                     &
                                    + Ffld(5,k1+1,m0  ,n1  )*R0( 9,1))*R0(12,2)           &
                                   +( Ffld(5,k1-1,m0+1,n1  )*R0( 7,1)                     &
                                    + Ffld(5,k1  ,m0+1,n1  )*R0( 8,1)                     &
                                    + Ffld(5,k1+1,m0+1,n1  )*R0( 9,1))*R0(13,2))*R0( 8,3) &
                                  +(( Ffld(5,k1-1,m0-1,n1+1)*R0( 7,1)                     &
                                    + Ffld(5,k1  ,m0-1,n1+1)*R0( 8,1)                     &
                                    + Ffld(5,k1+1,m0-1,n1+1)*R0( 9,1))*R0(11,2)           &
                                   +( Ffld(5,k1-1,m0  ,n1+1)*R0( 7,1)                     &
                                    + Ffld(5,k1  ,m0  ,n1+1)*R0( 8,1)                     &
                                    + Ffld(5,k1+1,m0  ,n1+1)*R0( 9,1))*R0(12,2)           &
                                   +( Ffld(5,k1-1,m0+1,n1+1)*R0( 7,1)                     &
                                    + Ffld(5,k1  ,m0+1,n1+1)*R0( 8,1)                     &
                                    + Ffld(5,k1+1,m0+1,n1+1)*R0( 9,1))*R0(13,2))*R0( 9,3) 
!-----------------------------------------------------------------------------------------
                      R(12)=R(12) +(( Ffld(6,k1-1,m1-1,n0-1)*R0( 7,1)                     &
                                    + Ffld(6,k1  ,m1-1,n0-1)*R0( 8,1)                     &
                                    + Ffld(6,k1+1,m1-1,n0-1)*R0( 9,1))*R0( 7,2)           &
                                   +( Ffld(6,k1-1,m1  ,n0-1)*R0( 7,1)                     &
                                    + Ffld(6,k1  ,m1  ,n0-1)*R0( 8,1)                     &
                                    + Ffld(6,k1+1,m1  ,n0-1)*R0( 9,1))*R0( 8,2)           &
                                   +( Ffld(6,k1-1,m1+1,n0-1)*R0( 7,1)                     &
                                    + Ffld(6,k1  ,m1+1,n0-1)*R0( 8,1)                     &
                                    + Ffld(6,k1+1,m1+1,n0-1)*R0( 9,1))*R0( 9,2))*R0(11,3) &
                                  +(( Ffld(6,k1-1,m1-1,n0  )*R0( 7,1)                     &
                                    + Ffld(6,k1  ,m1-1,n0  )*R0( 8,1)                     &
                                    + Ffld(6,k1+1,m1-1,n0  )*R0( 9,1))*R0( 7,2)           &
                                   +( Ffld(6,k1-1,m1  ,n0  )*R0( 7,1)                     &
                                     + Ffld(6,k1  ,m1  ,n0  )*R0( 8,1)                    &
                                    + Ffld(6,k1+1,m1  ,n0  )*R0( 9,1))*R0( 8,2)           &
                                   +( Ffld(6,k1-1,m1+1,n0  )*R0( 7,1)                     &
                                    + Ffld(6,k1  ,m1+1,n0  )*R0( 8,1)                     &
                                    + Ffld(6,k1+1,m1+1,n0  )*R0( 9,1))*R0( 9,2))*R0(12,3) &
                                  +(( Ffld(6,k1-1,m1-1,n0+1)*R0( 7,1)                     &
                                    + Ffld(6,k1  ,m1-1,n0+1)*R0( 8,1)                     &
                                    + Ffld(6,k1+1,m1-1,n0+1)*R0( 9,1))*R0( 7,2)           &
                                   +( Ffld(6,k1-1,m1  ,n0+1)*R0( 7,1)                     &
                                    + Ffld(6,k1  ,m1  ,n0+1)*R0( 8,1)                     &
                                    + Ffld(6,k1+1,m1  ,n0+1)*R0( 9,1))*R0( 8,2)           &
                                   +( Ffld(6,k1-1,m1+1,n0+1)*R0( 7,1)                     &
                                    + Ffld(6,k1  ,m1+1,n0+1)*R0( 8,1)                     &
                                    + Ffld(6,k1+1,m1+1,n0+1)*R0( 9,1))*R0( 9,2))*R0(13,3) 

! =====================================================
!    Calculation of particle movement
! =====================================================
!
                    !---IMF---
                      do imf=1,6
                        R(imf+6)=R(imf+6)+FIMF(imf)
                      end do
                      R( 7)=R( 7)*R(16)
                      R( 8)=R( 8)*R(16)
                      R( 9)=R( 9)*R(16)
                      R(10)=R(10)*R(16)
                      R(11)=R(11)*R(16)
                      R(12)=R(12)*R(16)
!
! ---  u^- = u^{n-1/2} + (qdt/2m)*E  ----
!TATSUKI altered this part. NEED TO CHANGE BACK TO ORIGINAL WHEN INCLIDING E
                      R(17)=R(13)+R(7)
                      R(18)=R(14)+R(8)
                      R(19)=R(15)+R(9)
!!$            R(17)=R(13)
!!$            R(18)=R(14)
!!$            R(19)=R(15)
!
! ---  gamma = sqrt{ 1 + |u^-|^2 }  ---
                      R(20)=1.d0/SQRT( ONE + R(17)*R(17) + R(18)*R(18) + R(19)*R(19) )
!
! ---  h = (qdt/2m)*(1/gamma)*B  ---
!TATSUKI altered this part. NEED TO CHANGE BACK TO ORIGINAL WHEN INCLIDING E
                      R(10)=R(10)*R(20)
                      R(11)=R(11)*R(20)
                      R(12)=R(12)*R(20)
!            R(10)=0d0
!            R(11)=0d0
!            R(12)=0d0
!
! ---  coefficients for u^{n+1/2}  ---
                      R(20) = R(10)*R(10) + R(11)*R(11) + R(12)*R(12)
                      R(21) = (ONE - R(20))
                      R(22) = ONE/(ONE + R(20))
                      R(20) = R(10)*R(17) + R(11)*R(18) + R(12)*R(19) 
!
! ---  update velocities
!      u^{n+1/2} = (1-h^2)/(1+h^2)*u^-
!                + 2/(1+h^2)*( (h,u^-) + [u^- x h] ) + (qdt/2m)*E
                      R(13) = R(7) + R(22)*(R(21)*R(17) + 2.d0*(R(20)*R(10) + R(18)*R(12)-R(19)*R(11)) ) 
                      R(14) = R(8) + R(22)*(R(21)*R(18) + 2.d0*(R(20)*R(11) + R(19)*R(10)-R(17)*R(12)) ) 
                      R(15) = R(9) + R(22)*(R(21)*R(19) + 2.d0*(R(20)*R(12) + R(17)*R(11)-R(18)*R(10)) )

!
! ---  update coordinates
!      x^{n+1} = x^{n} + u^{n+1/2}*dt/gamma^{n+1/2}
                      R(20)= R0(-1,3)/SQRT( ONE + R(13)*R(13) + R(14)*R(14) + R(15)*R(15) )
!
                      R(1) = R(1) + R(13)*R(20)
                      R(2) = R(2) + R(14)*R(20)
                      R(3) = R(3) + R(15)*R(20) 

!
! =====================================================
!    Form-factor calculation after particle movement
! =====================================================
! -- compute new 1D form-factors --
!
!print *,"compute new 1D form-factors start",rank
                      do k=1,3
                         R0( 4,k) = R(k)*R0(1,k) + HALF
                         j0( 2,k) = NINT( R0(4,k) )                ! j=1 corresponds to x=0
                         R0( 4,k) = R0(4,k) - j0(2,k)
                         j0( 1,k) = j0(2,k) - j0(3,k)              ! absolute shift of a particle = -1,0,1
                         j0( 4,k) = j0(3,k) + j0(1,k)*(1-j0(1,k))/2 - R0(2,k) ! *OR* = j0(k) + MIN(jd,0)
                         j0( 2,k) = j0(1,k)*(1+j0(1,k))/2 + 7                ! *OR* = MAX(jd,0) - 1 + 2
                         R0( 5,k) = R0(4,k)*R0(4,k)
                         R0( 6,k) = EIGHTH + HALF*R0(5,k)
                         R0( 4,k) = HALF*R0(4,k)
                         R0( 7,k) = ZERO
                         R0(10,k) = ZERO
!
                         j0( 1,k) = (sign(1,j0(1,k))-1)/2
                         R0(14,k) = R0(14+j0(1,k),k)
                         R0(13,k) = R0(13+j0(1,k),k)
                         R0(12,k) = R0(12+j0(1,k),k)
                         R0(11,k) = R0(11+j0(1,k),k)               ! W0(-2,k)
!
                         R0(j0(2,k)  ,k) = R0(6,k) - R0(4,k)         ! 0.5*(0.5-h1)**2
                         R0(j0(2,k)+1,k) = THR_FOURTH - R0(5,k)      ! 0.75 - h1**2
                         R0(j0(2,k)+2,k) = R0(6,k) + R0(4,k)         ! 0.5*(0.5+h1)**2
!
                         R0( 7,k) = R0( 7,k) - R0(11,k)
                         R0( 8,k) = R0( 8,k) - R0(12,k)
                         R0( 9,k) = R0( 9,k) - R0(13,k)
                         R0(10,k) = R0(10,k) - R0(14,k)

                         R0( 4,k)  = R0(7,k)
                         R0( 5,k)  = R0(7,k) + R0(8,k)
                         R0( 6,k)  = R0(7,k) + R0(8,k) + R0(9,k)
                      end do
!
! ==========================================================
!    Compute the decomposition of the total density change
! ==========================================================
! -- Contribute to current density of corresponding oct --
! -- Dont forget to multiply the current by cdx(*) in the Maxwell solver! --


                      j1(1) = j0(4,1)+ii+1
                      j1(2) = j0(4,2)+jj+1
                      j1(3) = j0(4,3)+kk+1

                      R0(3,3) = R(23)
!print *,"compute the decomposition of the total density change start",rank
                      do j=1,4
                         R0( 3,2) = (     R0(j+10,2) + HALF *R0(j+6,2))*R0(3,3)
                         R0(16,2) = (HALF*R0(j+10,2) + THIRD*R0(j+6,2))*R0(3,3)
                         R0(15,3) = (     R0(j+10,3) + HALF *R0(j+6,3))*R0(3,3)
                         R0(16,3) = (HALF*R0(j+10,3) + THIRD*R0(j+6,3))*R0(3,3)
                         j1(4) = j+j1(3)
                         j1(5) = j+j1(2)
                         do i=1,4
                            R0(15,2) = R0(15,3)*R0(i+10,2) + R0(16,3)*R0(i+6,2)
                            R0(15,1) = R0(15,3)*R0(i+10,1) + R0(16,3)*R0(i+6,1)
                            R0(16,1) = R0( 3,2)*R0(i+10,1) + R0(16,2)*R0(i+6,1)
                            j1(6) = i+j1(1)
                            j1(7) = i+j1(2)
                            Cfld(1, 2+j1(1),j1(7),j1(4),particle)= Cfld(1, 2+j1(1),j1(7),j1(4),particle)+ R0(15,2)*R0(4,1)
                            Cfld(1, 3+j1(1),j1(7),j1(4),particle)= Cfld(1, 3+j1(1),j1(7),j1(4),particle)+ R0(15,2)*R0(5,1)
                            Cfld(1, 4+j1(1),j1(7),j1(4),particle)= Cfld(1, 4+j1(1),j1(7),j1(4),particle)+ R0(15,2)*R0(6,1)
                            !
                            Cfld(2,j1(6), 2+j1(2),j1(4),particle)= Cfld(2,j1(6), 2+j1(2),j1(4),particle)+ R0(15,1)*R0(4,2)
                            Cfld(2,j1(6), 3+j1(2),j1(4),particle)= Cfld(2,j1(6), 3+j1(2),j1(4),particle)+ R0(15,1)*R0(5,2)
                            Cfld(2,j1(6), 4+j1(2),j1(4),particle)= Cfld(2,j1(6), 4+j1(2),j1(4),particle)+ R0(15,1)*R0(6,2)
                            !
                            Cfld(3,j1(6),j1(5), 2+j1(3),particle)= Cfld(3,j1(6),j1(5), 2+j1(3),particle)+ R0(16,1)*R0(4,3)
                            Cfld(3,j1(6),j1(5), 3+j1(3),particle)= Cfld(3,j1(6),j1(5), 3+j1(3),particle)+ R0(16,1)*R0(5,3)
                            Cfld(3,j1(6),j1(5), 4+j1(3),particle)= Cfld(3,j1(6),j1(5), 4+j1(3),particle)+ R0(16,1)*R0(6,3)
                         end do
                      end do
!
                      pptc%R(1) = R(1)
                      pptc%R(2) = R(2)
                      pptc%R(3) = R(3)


                      pptc%R(4) = R(13)
                      pptc%R(5) = R(14)
                      pptc%R(6) = R(15)
                      pptc%R(7) = (R(13)-R(4))*0.5d0+R(13)
                      pptc%R(8) = (R(14)-R(5))*0.5d0+R(14)
                      pptc%R(9) = (R(15)-R(6))*0.5d0+R(15)

                      pptc => pptc%prtNxt
            !
                   end do ! particle loop
                end do
             end do
          end do

          cfld(1,:,:,:,:)=cfld(1,:,:,:,:)*dx(1)/dt
          cfld(2,:,:,:,:)=cfld(2,:,:,:,:)*dx(2)/dt
          cfld(3,:,:,:,:)=cfld(3,:,:,:,:)*dx(3)/dt
          
          do ipp=1,Ionsorts
           do k=-1,1
             do j=-1,1
                do i=-1,1

                   p0 => Mesh(ip(i,j,k))

                   ii = 4+i*2
                   jj = 4+j*2
                   kk = 4+k*2

                   p0%octCh1%C(1+iomp2) = p0%octCh1%C(1+iomp2) + Cfld(1,ii  ,jj,kk,ipp)
                   p0%octCh1%C(2+iomp2) = p0%octCh1%C(2+iomp2) + Cfld(2,ii  ,jj,kk,ipp)
                   p0%octCh1%C(3+iomp2) = p0%octCh1%C(3+iomp2) + Cfld(3,ii  ,jj,kk,ipp)
                   p0%octCh2%C(1+iomp2) = p0%octCh2%C(1+iomp2) + Cfld(1,ii+1,jj,kk,ipp) ! ch2
                   p0%octCh2%C(2+iomp2) = p0%octCh2%C(2+iomp2) + Cfld(2,ii+1,jj,kk,ipp)
                   p0%octCh2%C(3+iomp2) = p0%octCh2%C(3+iomp2) + Cfld(3,ii+1,jj,kk,ipp)
                   
                   if(eleJoutput==1 .or. ionJoutput==1)then
                   if(ipp==1)then
                        p0%octCh1%C(1+iomp2+3*iomp0) = p0%octCh1%C(1+iomp2+3*iomp0) + Cfld(1,ii  ,jj,kk,ipp)
                        p0%octCh1%C(2+iomp2+3*iomp0) = p0%octCh1%C(2+iomp2+3*iomp0) + Cfld(2,ii  ,jj,kk,ipp)
                        p0%octCh1%C(3+iomp2+3*iomp0) = p0%octCh1%C(3+iomp2+3*iomp0) + Cfld(3,ii  ,jj,kk,ipp)
                        p0%octCh2%C(1+iomp2+3*iomp0) = p0%octCh2%C(1+iomp2+3*iomp0) + Cfld(1,ii+1,jj,kk,ipp) ! ch2
                        p0%octCh2%C(2+iomp2+3*iomp0) = p0%octCh2%C(2+iomp2+3*iomp0) + Cfld(2,ii+1,jj,kk,ipp)
                        p0%octCh2%C(3+iomp2+3*iomp0) = p0%octCh2%C(3+iomp2+3*iomp0) + Cfld(3,ii+1,jj,kk,ipp)
                   end if
                   end if

                   jj=jj+1

                   p0%octCh3%C(1+iomp2) = p0%octCh3%C(1+iomp2) + Cfld(1,ii  ,jj,kk,ipp)
                   p0%octCh3%C(2+iomp2) = p0%octCh3%C(2+iomp2) + Cfld(2,ii  ,jj,kk,ipp)
                   p0%octCh3%C(3+iomp2) = p0%octCh3%C(3+iomp2) + Cfld(3,ii  ,jj,kk,ipp)
                   p0%octCh4%C(1+iomp2) = p0%octCh4%C(1+iomp2) + Cfld(1,ii+1,jj,kk,ipp) ! ch4 
                   p0%octCh4%C(2+iomp2) = p0%octCh4%C(2+iomp2) + Cfld(2,ii+1,jj,kk,ipp)
                   p0%octCh4%C(3+iomp2) = p0%octCh4%C(3+iomp2) + Cfld(3,ii+1,jj,kk,ipp)
                   
                   if(eleJoutput==1 .or. ionJoutput==1)then
                   if(ipp==1)then
                        p0%octCh3%C(1+iomp2+3*iomp0) = p0%octCh3%C(1+iomp2+3*iomp0) + Cfld(1,ii  ,jj,kk,ipp)
                        p0%octCh3%C(2+iomp2+3*iomp0) = p0%octCh3%C(2+iomp2+3*iomp0) + Cfld(2,ii  ,jj,kk,ipp)
                        p0%octCh3%C(3+iomp2+3*iomp0) = p0%octCh3%C(3+iomp2+3*iomp0) + Cfld(3,ii  ,jj,kk,ipp)
                        p0%octCh4%C(1+iomp2+3*iomp0) = p0%octCh4%C(1+iomp2+3*iomp0) + Cfld(1,ii+1,jj,kk,ipp) ! ch4 
                        p0%octCh4%C(2+iomp2+3*iomp0) = p0%octCh4%C(2+iomp2+3*iomp0) + Cfld(2,ii+1,jj,kk,ipp)
                        p0%octCh4%C(3+iomp2+3*iomp0) = p0%octCh4%C(3+iomp2+3*iomp0) + Cfld(3,ii+1,jj,kk,ipp)
                   end if
                   end if
                   jj=jj-1
                   kk=kk+1
                   p0%octCh5%C(1+iomp2) = p0%octCh5%C(1+iomp2) + Cfld(1,ii  ,jj,kk,ipp)
                   p0%octCh5%C(2+iomp2) = p0%octCh5%C(2+iomp2) + Cfld(2,ii  ,jj,kk,ipp)
                   p0%octCh5%C(3+iomp2) = p0%octCh5%C(3+iomp2) + Cfld(3,ii  ,jj,kk,ipp)
                   p0%octCh6%C(1+iomp2) = p0%octCh6%C(1+iomp2) + Cfld(1,ii+1,jj,kk,ipp) ! ch6
                   p0%octCh6%C(2+iomp2) = p0%octCh6%C(2+iomp2) + Cfld(2,ii+1,jj,kk,ipp)
                   p0%octCh6%C(3+iomp2) = p0%octCh6%C(3+iomp2) + Cfld(3,ii+1,jj,kk,ipp)
                   
                   if(eleJoutput==1 .or. ionJoutput==1)then
                   if(ipp==1)then
                        p0%octCh5%C(1+iomp2+3*iomp0) = p0%octCh5%C(1+iomp2+3*iomp0) + Cfld(1,ii  ,jj,kk,ipp)
                        p0%octCh5%C(2+iomp2+3*iomp0) = p0%octCh5%C(2+iomp2+3*iomp0) + Cfld(2,ii  ,jj,kk,ipp)
                        p0%octCh5%C(3+iomp2+3*iomp0) = p0%octCh5%C(3+iomp2+3*iomp0) + Cfld(3,ii  ,jj,kk,ipp)
                        p0%octCh6%C(1+iomp2+3*iomp0) = p0%octCh6%C(1+iomp2+3*iomp0) + Cfld(1,ii+1,jj,kk,ipp) ! ch6
                        p0%octCh6%C(2+iomp2+3*iomp0) = p0%octCh6%C(2+iomp2+3*iomp0) + Cfld(2,ii+1,jj,kk,ipp)
                        p0%octCh6%C(3+iomp2+3*iomp0) = p0%octCh6%C(3+iomp2+3*iomp0) + Cfld(3,ii+1,jj,kk,ipp)
                   end if
                   end if
                   jj=jj+1
                   p0%octCh7%C(1+iomp2) = p0%octCh7%C(1+iomp2) + Cfld(1,ii  ,jj,kk,ipp)
                   p0%octCh7%C(2+iomp2) = p0%octCh7%C(2+iomp2) + Cfld(2,ii  ,jj,kk,ipp)
                   p0%octCh7%C(3+iomp2) = p0%octCh7%C(3+iomp2) + Cfld(3,ii  ,jj,kk,ipp)
                   p0%octCh8%C(1+iomp2) = p0%octCh8%C(1+iomp2) + Cfld(1,ii+1,jj,kk,ipp) ! ch8
                   p0%octCh8%C(2+iomp2) = p0%octCh8%C(2+iomp2) + Cfld(2,ii+1,jj,kk,ipp)
                   p0%octCh8%C(3+iomp2) = p0%octCh8%C(3+iomp2) + Cfld(3,ii+1,jj,kk,ipp)
                   
                   if(eleJoutput==1 .or. ionJoutput==1)then
                   if(ipp==1)then
                        p0%octCh7%C(1+iomp2+3*iomp0) = p0%octCh7%C(1+iomp2+3*iomp0) + Cfld(1,ii  ,jj,kk,ipp)
                        p0%octCh7%C(2+iomp2+3*iomp0) = p0%octCh7%C(2+iomp2+3*iomp0) + Cfld(2,ii  ,jj,kk,ipp)
                        p0%octCh7%C(3+iomp2+3*iomp0) = p0%octCh7%C(3+iomp2+3*iomp0) + Cfld(3,ii  ,jj,kk,ipp)
                        p0%octCh8%C(1+iomp2+3*iomp0) = p0%octCh8%C(1+iomp2+3*iomp0) + Cfld(1,ii+1,jj,kk,ipp) ! ch8
                        p0%octCh8%C(2+iomp2+3*iomp0) = p0%octCh8%C(2+iomp2+3*iomp0) + Cfld(2,ii+1,jj,kk,ipp)
                        p0%octCh8%C(3+iomp2+3*iomp0) = p0%octCh8%C(3+iomp2+3*iomp0) + Cfld(3,ii+1,jj,kk,ipp)
                    end if
                    end if
                end do
             end do
          end do

          do k=-1,1
             do j=-1,1

                p1 => Mesh(ip(-1,j,k))
                p2 => Mesh(ip(1,j,k))

                jj = 4+j*2
                kk = 4+k*2

                do i=1+iomp2,3+iomp2
                   ii=i-iomp2
                   p1%octCh1%octNb1%C(i) = p1%octCh1%octNb1%C(i) + Cfld(ii,1,jj  ,kk  ,ipp)
                   p2%octCh2%octNb2%C(i) = p2%octCh2%octNb2%C(i) + Cfld(ii,8,jj  ,kk  ,ipp)
                   p1%octCh3%octNb1%C(i) = p1%octCh3%octNb1%C(i) + Cfld(ii,1,jj+1,kk  ,ipp)
                   p2%octCh4%octNb2%C(i) = p2%octCh4%octNb2%C(i) + Cfld(ii,8,jj+1,kk  ,ipp)
                   p1%octCh5%octNb1%C(i) = p1%octCh5%octNb1%C(i) + Cfld(ii,1,jj  ,kk+1,ipp)
                   p2%octCh6%octNb2%C(i) = p2%octCh6%octNb2%C(i) + Cfld(ii,8,jj  ,kk+1,ipp)
                   p1%octCh7%octNb1%C(i) = p1%octCh7%octNb1%C(i) + Cfld(ii,1,jj+1,kk+1,ipp)
                   p2%octCh8%octNb2%C(i) = p2%octCh8%octNb2%C(i) + Cfld(ii,8,jj+1,kk+1,ipp)
                   
                   if(eleJoutput==1 .or. ionJoutput==1)then
                   if(ipp==1)then
                        p1%octCh1%octNb1%C(i+3*iomp0) = p1%octCh1%octNb1%C(i+3*iomp0) + Cfld(ii,1,jj  ,kk  ,ipp)
                        p2%octCh2%octNb2%C(i+3*iomp0) = p2%octCh2%octNb2%C(i+3*iomp0) + Cfld(ii,8,jj  ,kk  ,ipp)
                        p1%octCh3%octNb1%C(i+3*iomp0) = p1%octCh3%octNb1%C(i+3*iomp0) + Cfld(ii,1,jj+1,kk  ,ipp)
                        p2%octCh4%octNb2%C(i+3*iomp0) = p2%octCh4%octNb2%C(i+3*iomp0) + Cfld(ii,8,jj+1,kk  ,ipp)
                        p1%octCh5%octNb1%C(i+3*iomp0) = p1%octCh5%octNb1%C(i+3*iomp0) + Cfld(ii,1,jj  ,kk+1,ipp)
                        p2%octCh6%octNb2%C(i+3*iomp0) = p2%octCh6%octNb2%C(i+3*iomp0) + Cfld(ii,8,jj  ,kk+1,ipp)
                        p1%octCh7%octNb1%C(i+3*iomp0) = p1%octCh7%octNb1%C(i+3*iomp0) + Cfld(ii,1,jj+1,kk+1,ipp)
                        p2%octCh8%octNb2%C(i+3*iomp0) = p2%octCh8%octNb2%C(i+3*iomp0) + Cfld(ii,8,jj+1,kk+1,ipp)
                   end if
                   end if
                end do

                p1 => Mesh(ip(j,-1,k))         
                p2 => Mesh(ip(j,1,k))

                do i=1+iomp2,3+iomp2
                   ii=i-iomp2
                   p1%octCh1%octNb3%C(i) = p1%octCh1%octNb3%C(i) + Cfld(ii,jj  ,1,kk  ,ipp)
                   p1%octCh2%octNb3%C(i) = p1%octCh2%octNb3%C(i) + Cfld(ii,jj  ,1,kk  ,ipp)
                   p2%octCh3%octNb4%C(i) = p2%octCh3%octNb4%C(i) + Cfld(ii,jj+1,8,kk  ,ipp)
                   p2%octCh4%octNb4%C(i) = p2%octCh4%octNb4%C(i) + Cfld(ii,jj+1,8,kk  ,ipp)
                   p1%octCh5%octNb3%C(i) = p1%octCh5%octNb3%C(i) + Cfld(ii,jj  ,1,kk+1,ipp)
                   p1%octCh6%octNb3%C(i) = p1%octCh6%octNb3%C(i) + Cfld(ii,jj  ,1,kk+1,ipp)
                   p2%octCh7%octNb4%C(i) = p2%octCh7%octNb4%C(i) + Cfld(ii,jj+1,8,kk+1,ipp)
                   p2%octCh8%octNb4%C(i) = p2%octCh8%octNb4%C(i) + Cfld(ii,jj+1,8,kk+1,ipp)
                   
                   if(eleJoutput==1 .or. ionJoutput==1)then
                   if(ipp==1)then
                        p1%octCh1%octNb3%C(i+3*iomp0) = p1%octCh1%octNb3%C(i+3*iomp0) + Cfld(ii,jj  ,1,kk  ,ipp)
                        p1%octCh2%octNb3%C(i+3*iomp0) = p1%octCh2%octNb3%C(i+3*iomp0) + Cfld(ii,jj  ,1,kk  ,ipp)
                        p2%octCh3%octNb4%C(i+3*iomp0) = p2%octCh3%octNb4%C(i+3*iomp0) + Cfld(ii,jj+1,8,kk  ,ipp)
                        p2%octCh4%octNb4%C(i+3*iomp0) = p2%octCh4%octNb4%C(i+3*iomp0) + Cfld(ii,jj+1,8,kk  ,ipp)
                        p1%octCh5%octNb3%C(i+3*iomp0) = p1%octCh5%octNb3%C(i+3*iomp0) + Cfld(ii,jj  ,1,kk+1,ipp)
                        p1%octCh6%octNb3%C(i+3*iomp0) = p1%octCh6%octNb3%C(i+3*iomp0) + Cfld(ii,jj  ,1,kk+1,ipp)
                        p2%octCh7%octNb4%C(i+3*iomp0) = p2%octCh7%octNb4%C(i+3*iomp0) + Cfld(ii,jj+1,8,kk+1,ipp)
                        p2%octCh8%octNb4%C(i+3*iomp0) = p2%octCh8%octNb4%C(i+3*iomp0) + Cfld(ii,jj+1,8,kk+1,ipp)
                   end if
                   end if
                end do

                p1 => Mesh(ip(j,k,-1))
                p2 => Mesh(ip(j,k,1))

                do i=1+iomp2,3+iomp2
                   ii=i-iomp2
                   p1%octCh1%octNb5%C(i) = p1%octCh1%octNb5%C(i) + Cfld(ii,jj  ,kk  ,1,ipp)
                   p1%octCh2%octNb5%C(i) = p1%octCh2%octNb5%C(i) + Cfld(ii,jj  ,kk  ,1,ipp)
                   p1%octCh3%octNb5%C(i) = p1%octCh3%octNb5%C(i) + Cfld(ii,jj+1,kk  ,1,ipp)
                   p1%octCh4%octNb5%C(i) = p1%octCh4%octNb5%C(i) + Cfld(ii,jj+1,kk  ,1,ipp)
                   p2%octCh5%octNb6%C(i) = p2%octCh5%octNb6%C(i) + Cfld(ii,jj  ,kk+1,8,ipp)
                   p2%octCh6%octNb6%C(i) = p2%octCh6%octNb6%C(i) + Cfld(ii,jj  ,kk+1,8,ipp)
                   p2%octCh7%octNb6%C(i) = p2%octCh7%octNb6%C(i) + Cfld(ii,jj+1,kk+1,8,ipp)
                   p2%octCh8%octNb6%C(i) = p2%octCh8%octNb6%C(i) + Cfld(ii,jj+1,kk+1,8,ipp)
                   
                   if(eleJoutput==1 .or. ionJoutput==1)then
                   if(ipp==1)then
                        p1%octCh1%octNb5%C(i+3*iomp0) = p1%octCh1%octNb5%C(i+3*iomp0) + Cfld(ii,jj  ,kk  ,1,ipp)
                        p1%octCh2%octNb5%C(i+3*iomp0) = p1%octCh2%octNb5%C(i+3*iomp0) + Cfld(ii,jj  ,kk  ,1,ipp)
                        p1%octCh3%octNb5%C(i+3*iomp0) = p1%octCh3%octNb5%C(i+3*iomp0) + Cfld(ii,jj+1,kk  ,1,ipp)
                        p1%octCh4%octNb5%C(i+3*iomp0) = p1%octCh4%octNb5%C(i+3*iomp0) + Cfld(ii,jj+1,kk  ,1,ipp)
                        p2%octCh5%octNb6%C(i+3*iomp0) = p2%octCh5%octNb6%C(i+3*iomp0) + Cfld(ii,jj  ,kk+1,8,ipp)
                        p2%octCh6%octNb6%C(i+3*iomp0) = p2%octCh6%octNb6%C(i+3*iomp0) + Cfld(ii,jj  ,kk+1,8,ipp)
                        p2%octCh7%octNb6%C(i+3*iomp0) = p2%octCh7%octNb6%C(i+3*iomp0) + Cfld(ii,jj+1,kk+1,8,ipp)
                        p2%octCh8%octNb6%C(i+3*iomp0) = p2%octCh8%octNb6%C(i+3*iomp0) + Cfld(ii,jj+1,kk+1,8,ipp)
                   end if
                   end if
                end do
             end do
          end do

          do k=-1,1

             p1 => Mesh(ip(-1,-1,k))
             p2 => Mesh(ip(1,-1,k)) 
             p3 => Mesh(ip(-1,1,k))
             p4 => Mesh(ip(1,1,k))

             kk = 4+k*2

             do i=1+iomp2,3+iomp2
                ii=i-iomp2
                p1%octCh1%octNb1%octNb3%C(i) = p1%octCh1%octNb1%octNb3%C(i) + Cfld(ii,1,1,kk  ,ipp)
                p2%octCh2%octNb2%octNb3%C(i) = p2%octCh2%octNb2%octNb3%C(i) + Cfld(ii,8,1,kk  ,ipp)
                p3%octCh3%octNb1%octNb4%C(i) = p3%octCh3%octNb1%octNb4%C(i) + Cfld(ii,1,8,kk  ,ipp)
                p4%octCh4%octNb2%octNb4%C(i) = p4%octCh4%octNb2%octNb4%C(i) + Cfld(ii,8,8,kk  ,ipp)
                p1%octCh5%octNb1%octNb3%C(i) = p1%octCh5%octNb1%octNb3%C(i) + Cfld(ii,1,1,kk+1,ipp)
                p2%octCh6%octNb2%octNb3%C(i) = p2%octCh6%octNb2%octNb3%C(i) + Cfld(ii,8,1,kk+1,ipp)
                p3%octCh7%octNb1%octNb4%C(i) = p3%octCh7%octNb1%octNb4%C(i) + Cfld(ii,1,8,kk+1,ipp)
                p4%octCh8%octNb2%octNb4%C(i) = p4%octCh8%octNb2%octNb4%C(i) + Cfld(ii,8,8,kk+1,ipp)
                
                if(eleJoutput==1 .or. ionJoutput==1)then
                if(ipp==1)then
                    p1%octCh1%octNb1%octNb3%C(i+3*iomp0) = p1%octCh1%octNb1%octNb3%C(i+3*iomp0) + Cfld(ii,1,1,kk  ,ipp)
                    p2%octCh2%octNb2%octNb3%C(i+3*iomp0) = p2%octCh2%octNb2%octNb3%C(i+3*iomp0) + Cfld(ii,8,1,kk  ,ipp)
                    p3%octCh3%octNb1%octNb4%C(i+3*iomp0) = p3%octCh3%octNb1%octNb4%C(i+3*iomp0) + Cfld(ii,1,8,kk  ,ipp)
                    p4%octCh4%octNb2%octNb4%C(i+3*iomp0) = p4%octCh4%octNb2%octNb4%C(i+3*iomp0) + Cfld(ii,8,8,kk  ,ipp)
                    p1%octCh5%octNb1%octNb3%C(i+3*iomp0) = p1%octCh5%octNb1%octNb3%C(i+3*iomp0) + Cfld(ii,1,1,kk+1,ipp)
                    p2%octCh6%octNb2%octNb3%C(i+3*iomp0) = p2%octCh6%octNb2%octNb3%C(i+3*iomp0) + Cfld(ii,8,1,kk+1,ipp)
                    p3%octCh7%octNb1%octNb4%C(i+3*iomp0) = p3%octCh7%octNb1%octNb4%C(i+3*iomp0) + Cfld(ii,1,8,kk+1,ipp)
                    p4%octCh8%octNb2%octNb4%C(i+3*iomp0) = p4%octCh8%octNb2%octNb4%C(i+3*iomp0) + Cfld(ii,8,8,kk+1,ipp)
                end if
                end if
             end do


             p1 => Mesh(ip(-1,k,-1))
             p2 => Mesh(ip(1,k,-1))
             p3 => Mesh(ip(-1,k,1))
             p4 => Mesh(ip(1,k,1))

             do i=1+iomp2,3+iomp2
                ii=i-iomp2
                p1%octCh1%octNb1%octNb5%C(i) = p1%octCh1%octNb1%octNb5%C(i) + Cfld(ii,1,kk  ,1,ipp)
                p2%octCh2%octNb2%octNb5%C(i) = p2%octCh2%octNb2%octNb5%C(i) + Cfld(ii,8,kk  ,1,ipp)
                p1%octCh3%octNb1%octNb5%C(i) = p1%octCh3%octNb1%octNb5%C(i) + Cfld(ii,1,kk+1,1,ipp)
                p2%octCh4%octNb2%octNb5%C(i) = p2%octCh4%octNb2%octNb5%C(i) + Cfld(ii,8,kk+1,1,ipp)
                p3%octCh5%octNb1%octNb6%C(i) = p3%octCh5%octNb1%octNb6%C(i) + Cfld(ii,1,kk  ,8,ipp)
                p4%octCh6%octNb2%octNb6%C(i) = p4%octCh6%octNb2%octNb6%C(i) + Cfld(ii,8,kk  ,8,ipp)
                p3%octCh7%octNb1%octNb6%C(i) = p3%octCh7%octNb1%octNb6%C(i) + Cfld(ii,1,kk+1,8,ipp)
                p4%octCh8%octNb2%octNb6%C(i) = p4%octCh8%octNb2%octNb6%C(i) + Cfld(ii,8,kk+1,8,ipp)
                
                if(eleJoutput==1 .or. ionJoutput==1)then
                if(ipp==1)then
                    p1%octCh1%octNb1%octNb5%C(i+3*iomp0) = p1%octCh1%octNb1%octNb5%C(i+3*iomp0) + Cfld(ii,1,kk  ,1,ipp)
                    p2%octCh2%octNb2%octNb5%C(i+3*iomp0) = p2%octCh2%octNb2%octNb5%C(i+3*iomp0) + Cfld(ii,8,kk  ,1,ipp)
                    p1%octCh3%octNb1%octNb5%C(i+3*iomp0) = p1%octCh3%octNb1%octNb5%C(i+3*iomp0) + Cfld(ii,1,kk+1,1,ipp)
                    p2%octCh4%octNb2%octNb5%C(i+3*iomp0) = p2%octCh4%octNb2%octNb5%C(i+3*iomp0) + Cfld(ii,8,kk+1,1,ipp)
                    p3%octCh5%octNb1%octNb6%C(i+3*iomp0) = p3%octCh5%octNb1%octNb6%C(i+3*iomp0) + Cfld(ii,1,kk  ,8,ipp)
                    p4%octCh6%octNb2%octNb6%C(i+3*iomp0) = p4%octCh6%octNb2%octNb6%C(i+3*iomp0) + Cfld(ii,8,kk  ,8,ipp)
                    p3%octCh7%octNb1%octNb6%C(i+3*iomp0) = p3%octCh7%octNb1%octNb6%C(i+3*iomp0) + Cfld(ii,1,kk+1,8,ipp)
                    p4%octCh8%octNb2%octNb6%C(i+3*iomp0) = p4%octCh8%octNb2%octNb6%C(i+3*iomp0) + Cfld(ii,8,kk+1,8,ipp)
                end if
                end if
             end do

             p1 => Mesh(ip(k,-1,-1))       
             p2 => Mesh(ip(k,1,-1))         
             p3 => Mesh(ip(k,-1,1))   
             p4 => Mesh(ip(k,1,1))

             do i=1+iomp2,3+iomp2
                ii=i-iomp2
                p1%octCh1%octNb3%octNb5%C(i) = p1%octCh1%octNb3%octNb5%C(i) + Cfld(ii,kk  ,1,1,ipp)
                p1%octCh2%octNb3%octNb5%C(i) = p1%octCh2%octNb3%octNb5%C(i) + Cfld(ii,kk+1,1,1,ipp)
                p2%octCh3%octNb4%octNb5%C(i) = p2%octCh3%octNb4%octNb5%C(i) + Cfld(ii,kk  ,8,1,ipp)
                p2%octCh4%octNb4%octNb5%C(i) = p2%octCh4%octNb4%octNb5%C(i) + Cfld(ii,kk+1,8,1,ipp)
                p3%octCh5%octNb3%octNb6%C(i) = p3%octCh5%octNb3%octNb6%C(i) + Cfld(ii,kk  ,1,8,ipp)
                p3%octCh6%octNb3%octNb6%C(i) = p3%octCh6%octNb3%octNb6%C(i) + Cfld(ii,kk+1,1,8,ipp)
                p4%octCh7%octNb4%octNb6%C(i) = p4%octCh7%octNb4%octNb6%C(i) + Cfld(ii,kk  ,8,8,ipp)
                p4%octCh8%octNb4%octNb6%C(i) = p4%octCh8%octNb4%octNb6%C(i) + Cfld(ii,kk+1,8,8,ipp)
                
                if(eleJoutput==1 .or. ionJoutput==1)then
                if(ipp==1)then
                    p1%octCh1%octNb3%octNb5%C(i+3*iomp0) = p1%octCh1%octNb3%octNb5%C(i+3*iomp0) + Cfld(ii,kk  ,1,1,ipp)
                    p1%octCh2%octNb3%octNb5%C(i+3*iomp0) = p1%octCh2%octNb3%octNb5%C(i+3*iomp0) + Cfld(ii,kk+1,1,1,ipp)
                    p2%octCh3%octNb4%octNb5%C(i+3*iomp0) = p2%octCh3%octNb4%octNb5%C(i+3*iomp0) + Cfld(ii,kk  ,8,1,ipp)
                    p2%octCh4%octNb4%octNb5%C(i+3*iomp0) = p2%octCh4%octNb4%octNb5%C(i+3*iomp0) + Cfld(ii,kk+1,8,1,ipp)
                    p3%octCh5%octNb3%octNb6%C(i+3*iomp0) = p3%octCh5%octNb3%octNb6%C(i+3*iomp0) + Cfld(ii,kk  ,1,8,ipp)
                    p3%octCh6%octNb3%octNb6%C(i+3*iomp0) = p3%octCh6%octNb3%octNb6%C(i+3*iomp0) + Cfld(ii,kk+1,1,8,ipp)
                    p4%octCh7%octNb4%octNb6%C(i+3*iomp0) = p4%octCh7%octNb4%octNb6%C(i+3*iomp0) + Cfld(ii,kk  ,8,8,ipp)
                    p4%octCh8%octNb4%octNb6%C(i+3*iomp0) = p4%octCh8%octNb4%octNb6%C(i+3*iomp0) + Cfld(ii,kk+1,8,8,ipp)
                end if
                end if
             end do
          end do


      
          p1 => Mesh(ip(-1,-1,-1))      
          p2 => Mesh(ip(1,-1,-1))
          p3 => Mesh(ip(-1,1,-1))
          p4 => Mesh(ip(1,1,-1))     
          p5 => Mesh(ip(-1,-1,1))      
          p6 => Mesh(ip(1,-1,1))     
          p7 => Mesh(ip(-1,1,1))    
          p8 => Mesh(ip(1,1,1))

          do i=1+iomp2,3+iomp2
             ii=i-iomp2
             p1%octCh1%octNb1%octNb3%octNb5%C(i) = p1%octCh1%octNb1%octNb3%octNb5%C(i) + Cfld(ii,1,1,1,ipp)
             p2%octCh2%octNb2%octNb3%octNb5%C(i) = p2%octCh2%octNb2%octNb3%octNb5%C(i) + Cfld(ii,8,1,1,ipp)
             p3%octCh3%octNb1%octNb4%octNb5%C(i) = p3%octCh3%octNb1%octNb4%octNb5%C(i) + Cfld(ii,1,8,1,ipp)
             p4%octCh4%octNb2%octNb4%octNb5%C(i) = p4%octCh4%octNb2%octNb4%octNb5%C(i) + Cfld(ii,8,8,1,ipp)
             p5%octCh5%octNb1%octNb3%octNb6%C(i) = p5%octCh5%octNb1%octNb3%octNb6%C(i) + Cfld(ii,1,1,8,ipp)
             p6%octCh6%octNb2%octNb3%octNb6%C(i) = p6%octCh6%octNb2%octNb3%octNb6%C(i) + Cfld(ii,8,1,8,ipp)
             p7%octCh7%octNb1%octNb4%octNb6%C(i) = p7%octCh7%octNb1%octNb4%octNb6%C(i) + Cfld(ii,1,8,8,ipp)
             p8%octCh8%octNb2%octNb4%octNb6%C(i) = p8%octCh8%octNb2%octNb4%octNb6%C(i) + Cfld(ii,8,8,8,ipp)
             
             if(eleJoutput==1 .or. ionJoutput==1)then
             if(ipp==1)then
                p1%octCh1%octNb1%octNb3%octNb5%C(i+3*iomp0) = p1%octCh1%octNb1%octNb3%octNb5%C(i+3*iomp0) + Cfld(ii,1,1,1,ipp)
                p2%octCh2%octNb2%octNb3%octNb5%C(i+3*iomp0) = p2%octCh2%octNb2%octNb3%octNb5%C(i+3*iomp0) + Cfld(ii,8,1,1,ipp)
                p3%octCh3%octNb1%octNb4%octNb5%C(i+3*iomp0) = p3%octCh3%octNb1%octNb4%octNb5%C(i+3*iomp0) + Cfld(ii,1,8,1,ipp)
                p4%octCh4%octNb2%octNb4%octNb5%C(i+3*iomp0) = p4%octCh4%octNb2%octNb4%octNb5%C(i+3*iomp0) + Cfld(ii,8,8,1,ipp)
                p5%octCh5%octNb1%octNb3%octNb6%C(i+3*iomp0) = p5%octCh5%octNb1%octNb3%octNb6%C(i+3*iomp0) + Cfld(ii,1,1,8,ipp)
                p6%octCh6%octNb2%octNb3%octNb6%C(i+3*iomp0) = p6%octCh6%octNb2%octNb3%octNb6%C(i+3*iomp0) + Cfld(ii,8,1,8,ipp)
                p7%octCh7%octNb1%octNb4%octNb6%C(i+3*iomp0) = p7%octCh7%octNb1%octNb4%octNb6%C(i+3*iomp0) + Cfld(ii,1,8,8,ipp)
                p8%octCh8%octNb2%octNb4%octNb6%C(i+3*iomp0) = p8%octCh8%octNb2%octNb4%octNb6%C(i+3*iomp0) + Cfld(ii,8,8,8,ipp)
             end if
             end if
          end do
         end do
       end do




       if(MaxID(1,iLv)>MinID(1,iLv))then
!$omp parallel do private(index,p0,k) shared(minindex,maxindex,iFs,iFe,Mesh)
          do index=MinID(1,iLv),MaxID(1,iLv)
             !do index=minindex,maxindex ! This index runs for iLv <=0 and inside the real space so p0 always
             ! points to Mesh
             p0 => Mesh(index)
             if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then 
                do k=1,3
                    do i=0,iomp0-1
                        p0%F(k+6) = p0%F(k+6) + p0%C(k+3*i)
                        p0%F(k+12) = p0%F(k+12) + p0%C(k+3*i+3*iomp0)
                    end do
                    p0%F(k+9) = (p0%F(k+6)-p0%F(k+9))*0.5d0+p0%F(k+6)
                    p0%F(k+15) = p0%F(k+6) - p0%F(k+12)
                end do
             endif
          end do
!$omp end parallel do 
       endif




       if(maxID(3,iLv)>minID(3,iLv))then
!$omp parallel do private(index,p0,k) shared(minindex,maxindex,iFs,iFe,Mesh)
          do index=minID(3,iLv),maxID(3,iLv)
             p0 => Mesh(index)
             if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then 
                do k=1,3
                   do i=0,iomp0-1
                        p0%F(k+6) = p0%F(k+6) + p0%C(k+3*i)
                        p0%F(k+12) = p0%F(k+12) + p0%C(k+3*i+3*iomp0)
                   end do
                   p0%F(k+9) = (p0%F(k+6)-p0%F(k+9))*0.5d0+p0%F(k+6)
                   p0%F(k+15) = p0%F(k+6) - p0%F(k+12)
                end do
             endif
          end do
!$omp end parallel do 
       endif
!
       !if(debugMode>=1)print*, rank, 'After move_particle'
       return
       !
     end subroutine move_particle
!
!**********************************************************************
      subroutine issue(itype, id, message)
! +-------------------------------------------------------------------+
! |                                                                   |
! | if itype=0, print WARNING, if >0 print ERROR and exit.            |
! |                                                                   |
! +-------------------------------------------------------------------+
!**********************************************************************
        implicit none
        integer(kind=4) :: itype
        character*(*)   :: id, message
!
        goto (100,999) itype+1
        goto 999

100     continue        ! Just a warning: itype=0
        write(*,*) '<WARNING',itype,'>',id,'|',message
        return

999     continue        ! Fatal error: itype>0
        write(*,*) '<ERROR',itype,'>',id,'|',message
!
        stop
!
      end subroutine issue
!
!----------------------------------------------------------------------
!    Subroutines for debug
!----------------------------------------------------------------------
!
!**********************************************************************
      subroutine density(iLv)
! +-------------------------------------------------------------------+
! |                                                                   |
! |   get particle number in the oct                                  |
! |                                                                   |
! +-------------------------------------------------------------------+
!**********************************************************************
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        use param
!---------------------------------------------------------------
        implicit none
        integer(kind=4) :: index,indexP,k
        integer(kind=4) :: iLv,is
        integer(kind=4) :: maxindex, minindex
        real(kind=8)::dxx,dvlevel
! -- pointers -
        type(oct),   pointer :: p0
        type(prtcl), pointer :: pptc
!---------------------------------------------------------------
!
        if(iLv.gt.Lvmax) return 
        maxindex = MaxID(1,iLv)
        minindex = MinID(1,iLv)
        if(maxindex.eq.minindex) return 

        do index=minindex,maxindex
           p0 => Mesh(index)
           do k=1,IonSorts
              p0%Z(k) = 0.d0
           end do
        end do
        dvlevel=0.5d0**iLv
        dxx=dvlevel**3.0d0
        do index=minindex,maxindex
           p0 => Mesh(index)
           if((p0%iFLG(1).le.3).and.(p0%iFLG(1).ge.0)) then 
              indexP = p0%octP
              pptc => p0%ptcl%prtnxt
              do k=1,indexP
                 is = pptc%Isort
                 if(is.gt.0) then 
                    p0%Z(is) = p0%Z(is) + 1.d0
                 end if
                 pptc => pptc%prtNxt
              end do
              do k=1,IonSorts
                p0%Z(k)=p0%Z(k)/dxx
              end do
           else if(p0%iFLG(1).gt.3) then 
            do k=1,IonSorts
                p0%Z(k)=(p0%octCh1%Z(k)+p0%octCh2%Z(k)  &
                        +p0%octCh3%Z(k)+p0%octCh4%Z(k)  &
                        +p0%octCh5%Z(k)+p0%octCh6%Z(k)  &
                        +p0%octCh7%Z(k)+p0%octCh8%Z(k) )*0.125d0
              end do
           end if
        end do
!
        return 
        end subroutine density

!**********************************************************************
      subroutine collect_ptcl_loops(ptcl_lps_per_prc)
! +-------------------------------------------------------------------+
! |                                                                   |
! |   get particle number in the oct                                  |
! |                                                                   |
! +-------------------------------------------------------------------+
!**********************************************************************
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        use param
        use mpi
        use message_passing_interface
!---------------------------------------------------------------
        implicit none
        integer(kind=4) :: index,indexP,k
        integer(kind=4) :: iLv,is, temp
        integer(kind=4) :: maxindex, minindex
        integer(kind=4) ,intent(out):: ptcl_lps_per_prc
! -- pointers -
        type(oct),   pointer :: p0
        type(prtcl), pointer :: pptc
!---------------------------------------------------------------
!
!        ptcl_loops_per_process = 0

        !if(debugMode>=1)print*,'Entering collect_ptcl_loops:',rank
        do iLv = LvMax, 0, -1
           maxindex = MaxID(1,iLv)
           minindex = MinID(1,iLv)

           temp=2**iLv
           if(maxindex .gt. minindex) then
              do index=minindex,maxindex
                 p0 => Mesh(index)
                 p0%ptcl_loops=0
!if(rank==5 .and. iLv==1) print*,'came here now'
!           do k=1,IonSorts
!              p0%Z(k) = 0.d0
!           end do
                 if((p0%iFLG(1).le.3).and.(p0%iFLG(1).ge.0)) then !there is no child below
                    indexP = p0%octP
                    !if(rank==0)print *,"p0%octP",p0%octP,rank
                    pptc => p0%ptcl%prtnxt

                    do k=1,indexP
                       is = pptc%Isort
                       if(is.gt.0) then 
                          p0%ptcl_loops = p0%ptcl_loops + temp
                         ! if(rank==0)print *,"p0%ptcl_loops",p0%ptcl_loops,rank
                       end if
                       pptc => pptc%prtNxt
                    end do
                 else if(p0%iFLG(1).gt.3) then !if there were children, then collect them up
                    p0%ptcl_loops=p0%ptcl_loops+p0%octCh1%ptcl_loops+p0%octCh2%ptcl_loops  &
                         +p0%octCh3%ptcl_loops+p0%octCh4%ptcl_loops  &
                         +p0%octCh5%ptcl_loops+p0%octCh6%ptcl_loops  &
                         +p0%octCh7%ptcl_loops+p0%octCh8%ptcl_loops  
                 end if
              end do
           endif
        enddo

        ptcl_lps_per_prc = 0
        do index = MinID(1, -1), MaxID(1, -1)
           !p0 => BMesh(index)
           p0=>Mesh(index)
!!$           p0 % octP = p0%octCh1% octP + p0%octCh2% octP &
!!$                + p0%octCh3% octP + p0%octCh4% octP &
!!$                + p0%octCh5% octP + p0%octCh6% octP &
!!$                + p0%octCh7% octP + p0%octCh8% octP 

           p0%ptcl_loops=p0%octCh1%ptcl_loops+p0%octCh2%ptcl_loops  &
                +p0%octCh3%ptcl_loops+p0%octCh4%ptcl_loops  &
                +p0%octCh5%ptcl_loops+p0%octCh6%ptcl_loops  &
                +p0%octCh7%ptcl_loops+p0%octCh8%ptcl_loops  
           ptcl_lps_per_prc = ptcl_lps_per_prc + p0%ptcl_loops
        enddo

        return 
    end subroutine collect_ptcl_loops

!***********************************************************************
subroutine back_push2
! +-------------------------------------------------------------------+
! |                                                                   |
! |  Computes force interpolation, advances particles                 |
! |   [C. K. Birdsall and A. B. Langdon, Plasma Physics Via           |
! |    Computer Simulation (Adam-Hilger, 1991)].                      |
! |                                                                   |
! |  Computes current density by density decomposition                |
! |  [T.Zh.Esirkepov, Exact charge conservation scheme for            |
! |   particle-in-cell simulation with an arbitrary form-factor,      |
! |   Computer Physics Communications, v.135 (2001), pp.144-153]      |
! |                                                                   |
! |                                                                   |
! |  2nd order form-factor. Particle occupies 27*dx*dy*dz cells.      |
! |                                                                   |
! |  W(X_{j}-x0) = 3/4 - (|X_{j}-x0|/dx)^2, if |X_{j}-x0|/dx <= 1/2 ; |
! |   or 1/2*( 3/2-|X_{j}-x0|/dx )^2, if 1/2 < |X_{j}-x0|/dx <= 3/2 ; |
! |   or 0, if |X_{j}-x0|/dx > 3/2 .                                  |
! |                                                                   |
! |  dx=dy=dz is strongly recommended!                                |
! |                                                                   |
! |                                                                   |
! |                                                                   |
! |   Argument parameters:                                            |
! |       i2step  : iLv+1 levels time step                            |
! |       istep   : iLv   levels time step                            |
! |       iLv     : target hierarchical level                         |
! |                                                                   |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        use param
!-----------------------------------------------------------------------
        implicit none
        real(kind=8) :: dvlevel, deltaT,dvlevel3
        real(kind=8) :: time_factor, charge_factor
        real(kind=8) :: R(9), Rchrg, Rmss
        integer(kind=4) :: i1,i2,i3,k,j,jd,index,iLv!,iFe,iFs,kk
        integer(kind=4) :: ii,ij,ik
        integer(kind=4) :: j0(3), j1(3)
        integer(kind=4) :: j00(3),j11(3)
        real(kind=8) :: h1,h2,h3, u1,u2,u3, v1,v2,v3, gamma, c1,c2
        real(kind=8) :: W0(-1:3,3) ! (-1:2) is enough, but bank conflict
        real(kind=8) :: W1(-1:3,3) ! (-1:2) is enough, but bank conflict
        real(kind=8) :: FF(6)
        real(kind=8) :: Ffld(6,-2:2,-2:2,-2:2)
        real(kind=8) :: Cfld(3,-1:2,-1:2,-1:2)
!
        real(kind=8) :: odxx(3), icel(3)
    !    real(kind=8) :: hh(3),jj(3)
!
        integer(kind=4) :: ix,iy,iz
        integer(kind=4) :: Nin!,ipx,ipy,ipz,index2
        integer(kind=4) :: ImovT,IsortT
        real(kind=8)    :: RT(6)
    !    real(kind=8)    :: ixx,iyy,izz
        integer(kind=8) :: indexPT,indexRT,indexP,k0
        integer(kind=8) :: minindex,maxindex,kp!,index0,istart,iend
!
! -- pointers --
        type(oct), pointer   :: p0!,pp0
        type(oct), pointer   :: ptmp
        type(prtcl), pointer :: pptc!,tempptr
!-----------------------------------------------------------------------
!
! == BEGIN: ============================================================
!
        print *, 'In "move_particle"...'
!
! -- Set temporal representative particles
        ImovT = 0
        IsortT= 0
        RT    = 0.d0
        indexPT= 0
        indexRT= 0
!
        iLv=0
! -- initialize / current density
        do index=MinID(1,iLv),MaxID(1,iLv)
           p0 => Mesh(index)
           do k=7,9
              p0%F(k+3) = ZERO
              p0%F(k  ) = ZERO
           enddo
           do k=1,48
              p0%C(k  ) = ZERO
           end do
        end do
!
! -- Set spacing parameters --
        dvlevel = 0.5d0**(iLv)
        dvlevel3=(dvlevel**3)
        deltaT  = dvlevel * dt
        do k=1,3
          odxx(k) = ONE/(dvlevel*dx(k))
        end do
!
        print *, ' iLv=',iLv
        print *, ' dt =',deltaT
!
! -- parallel loop
        maxindex = MaxID(1,iLv)
        minindex = MinID(1,iLv)
        do index=minindex,maxindex
          p0 => Mesh(index)                      
!
! -- Get cell index --
          ix = p0%iPOS(1)
          iy = p0%iPOS(2)
          iz = p0%iPOS(3)
          Nin     = 2**(LvMax-iLv)
          icel(1) = real(ix + Nin)/real(2*Nin)
          icel(2) = real(iy + Nin)/real(2*Nin)
          icel(3) = real(iz + Nin)/real(2*Nin)
!
! -- Field preparation --
          do k=1,6
!            ---------------      %----  X  ----%----  Y  ----%----  Z  ----%
             Ffld(k,-2,-2,-2) = p0%octNb1%octNb1%octNb3%octNb3%octNb5%octNb5%F(k)
             Ffld(k,-1,-2,-2) = p0%octNb1       %octNb3%octNb3%octNb5%octNb5%F(k)
             Ffld(k, 0,-2,-2) = p0              %octNb3%octNb3%octNb5%octNb5%F(k)
             Ffld(k, 1,-2,-2) = p0%octNb2       %octNb3%octNb3%octNb5%octNb5%F(k)
             Ffld(k, 2,-2,-2) = p0%octNb2%octNb2%octNb3%octNb3%octNb5%octNb5%F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2,-1,-2) = p0%octNb1%octNb1%octNb3       %octNb5%octNb5%F(k)
             Ffld(k,-1,-1,-2) = p0%octNb1       %octNb3       %octNb5%octNb5%F(k)
             Ffld(k, 0,-1,-2) = p0              %octNb3       %octNb5%octNb5%F(k)
             Ffld(k, 1,-1,-2) = p0%octNb2       %octNb3       %octNb5%octNb5%F(k)
             Ffld(k, 2,-1,-2) = p0%octNb2%octNb2%octNb3       %octNb5%octNb5%F(k)
!             --------------------------------------------------------------------
             Ffld(k,-2, 0,-2) = p0%octNb1%octNb1              %octNb5%octNb5%F(k)
             Ffld(k,-1, 0,-2) = p0%octNb1                     %octNb5%octNb5%F(k)
             Ffld(k, 0, 0,-2) = p0                            %octNb5%octNb5%F(k)
             Ffld(k, 1, 0,-2) = p0%octNb2                     %octNb5%octNb5%F(k)
             Ffld(k, 2, 0,-2) = p0%octNb2%octNb2              %octNb5%octNb5%F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 1,-2) = p0%octNb1%octNb1%octNb4       %octNb5%octNb5%F(k)
             Ffld(k,-1, 1,-2) = p0%octNb1       %octNb4       %octNb5%octNb5%F(k)
             Ffld(k, 0, 1,-2) = p0              %octNb4       %octNb5%octNb5%F(k)
             Ffld(k, 1, 1,-2) = p0%octNb2       %octNb4       %octNb5%octNb5%F(k)
             Ffld(k, 2, 1,-2) = p0%octNb2%octNb2%octNb4       %octNb5%octNb5%F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 2,-2) = p0%octNb1%octNb1%octNb4%octNb4%octNb5%octNb5%F(k)
             Ffld(k,-1, 2,-2) = p0%octNb1       %octNb4%octNb4%octNb5%octNb5%F(k)
             Ffld(k, 0, 2,-2) = p0              %octNb4%octNb4%octNb5%octNb5%F(k)
             Ffld(k, 1, 2,-2) = p0%octNb2       %octNb4%octNb4%octNb5%octNb5%F(k)
             Ffld(k, 2, 2,-2) = p0%octNb2%octNb2%octNb4%octNb4%octNb5%octNb5%F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2,-2,-1) = p0%octNb1%octNb1%octNb3%octNb3%octNb5       %F(k)
             Ffld(k,-1,-2,-1) = p0%octNb1       %octNb3%octNb3%octNb5       %F(k)
             Ffld(k, 0,-2,-1) = p0              %octNb3%octNb3%octNb5       %F(k)
             Ffld(k, 1,-2,-1) = p0%octNb2       %octNb3%octNb3%octNb5       %F(k)
             Ffld(k, 2,-2,-1) = p0%octNb2%octNb2%octNb3%octNb3%octNb5       %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2,-1,-1) = p0%octNb1%octNb1%octNb3       %octNb5       %F(k)
             Ffld(k,-1,-1,-1) = p0%octNb1       %octNb3       %octNb5       %F(k)
             Ffld(k, 0,-1,-1) = p0              %octNb3       %octNb5       %F(k)
             Ffld(k, 1,-1,-1) = p0%octNb2       %octNb3       %octNb5       %F(k)
             Ffld(k, 2,-1,-1) = p0%octNb2%octNb2%octNb3       %octNb5       %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 0,-1) = p0%octNb1%octNb1              %octNb5       %F(k)
             Ffld(k,-1, 0,-1) = p0%octNb1                     %octNb5       %F(k)
             Ffld(k, 0, 0,-1) = p0                            %octNb5       %F(k)
             Ffld(k, 1, 0,-1) = p0%octNb2                     %octNb5       %F(k)
             Ffld(k, 2, 0,-1) = p0%octNb2%octNb2              %octNb5       %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 1,-1) = p0%octNb1%octNb1%octNb4       %octNb5       %F(k)
             Ffld(k,-1, 1,-1) = p0%octNb1       %octNb4       %octNb5       %F(k)
             Ffld(k, 0, 1,-1) = p0              %octNb4       %octNb5       %F(k)
             Ffld(k, 1, 1,-1) = p0%octNb2       %octNb4       %octNb5       %F(k)
             Ffld(k, 2, 1,-1) = p0%octNb2%octNb2%octNb4       %octNb5       %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 2,-1) = p0%octNb1%octNb1%octNb4%octNb4%octNb5       %F(k)
             Ffld(k,-1, 2,-1) = p0%octNb1       %octNb4%octNb4%octNb5       %F(k)
             Ffld(k, 0, 2,-1) = p0              %octNb4%octNb4%octNb5       %F(k)
             Ffld(k, 1, 2,-1) = p0%octNb2       %octNb4%octNb4%octNb5       %F(k)
             Ffld(k, 2, 2,-1) = p0%octNb2%octNb2%octNb4%octNb4%octNb5       %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2,-2, 0) = p0%octNb1%octNb1%octNb3%octNb3              %F(k)
             Ffld(k,-1,-2, 0) = p0%octNb1       %octNb3%octNb3              %F(k)
             Ffld(k, 0,-2, 0) = p0              %octNb3%octNb3              %F(k)
             Ffld(k, 1,-2, 0) = p0%octNb2       %octNb3%octNb3              %F(k)
             Ffld(k, 2,-2, 0) = p0%octNb2%octNb2%octNb3%octNb3              %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2,-1, 0) = p0%octNb1%octNb1%octNb3                     %F(k)
             Ffld(k,-1,-1, 0) = p0%octNb1       %octNb3                     %F(k)
             Ffld(k, 0,-1, 0) = p0              %octNb3                     %F(k)
             Ffld(k, 1,-1, 0) = p0%octNb2       %octNb3                     %F(k)
             Ffld(k, 2,-1, 0) = p0%octNb2%octNb2%octNb3                     %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 0, 0) = p0%octNb1%octNb1                            %F(k)
             Ffld(k,-1, 0, 0) = p0%octNb1                                   %F(k)
             Ffld(k, 0, 0, 0) = p0                                          %F(k)
             Ffld(k, 1, 0, 0) = p0%octNb2                                   %F(k)
             Ffld(k, 2, 0, 0) = p0%octNb2%octNb2                            %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 1, 0) = p0%octNb1%octNb1%octNb4                     %F(k)
             Ffld(k,-1, 1, 0) = p0%octNb1       %octNb4                     %F(k)
             Ffld(k, 0, 1, 0) = p0              %octNb4                     %F(k)
             Ffld(k, 1, 1, 0) = p0%octNb2       %octNb4                     %F(k)
             Ffld(k, 2, 1, 0) = p0%octNb2%octNb2%octNb4                     %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 2, 0) = p0%octNb1%octNb1%octNb4%octNb4              %F(k)
             Ffld(k,-1, 2, 0) = p0%octNb1       %octNb4%octNb4              %F(k)
             Ffld(k, 0, 2, 0) = p0              %octNb4%octNb4              %F(k)
             Ffld(k, 1, 2, 0) = p0%octNb2       %octNb4%octNb4              %F(k)
             Ffld(k, 2, 2, 0) = p0%octNb2%octNb2%octNb4%octNb4              %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2,-2, 1) = p0%octNb1%octNb1%octNb3%octNb3%octNb6       %F(k)
             Ffld(k,-1,-2, 1) = p0%octNb1       %octNb3%octNb3%octNb6       %F(k)
             Ffld(k, 0,-2, 1) = p0              %octNb3%octNb3%octNb6       %F(k)
             Ffld(k, 1,-2, 1) = p0%octNb2       %octNb3%octNb3%octNb6       %F(k)
             Ffld(k, 2,-2, 1) = p0%octNb2%octNb2%octNb3%octNb3%octNb6       %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2,-1, 1) = p0%octNb1%octNb1%octNb3       %octNb6       %F(k)
             Ffld(k,-1,-1, 1) = p0%octNb1       %octNb3       %octNb6       %F(k)
             Ffld(k, 0,-1, 1) = p0              %octNb3       %octNb6       %F(k)
             Ffld(k, 1,-1, 1) = p0%octNb2       %octNb3       %octNb6       %F(k)
             Ffld(k, 2,-1, 1) = p0%octNb2%octNb2%octNb3       %octNb6       %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 0, 1) = p0%octNb1%octNb1              %octNb6       %F(k)
             Ffld(k,-1, 0, 1) = p0%octNb1                     %octNb6       %F(k)
             Ffld(k, 0, 0, 1) = p0                            %octNb6       %F(k)
             Ffld(k, 1, 0, 1) = p0%octNb2                     %octNb6       %F(k)
             Ffld(k, 2, 0, 1) = p0%octNb2%octNb2              %octNb6       %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 1, 1) = p0%octNb1%octNb1%octNb4       %octNb6       %F(k)
             Ffld(k,-1, 1, 1) = p0%octNb1       %octNb4       %octNb6       %F(k)
             Ffld(k, 0, 1, 1) = p0              %octNb4       %octNb6       %F(k)
             Ffld(k, 1, 1, 1) = p0%octNb2       %octNb4       %octNb6       %F(k)
             Ffld(k, 2, 1, 1) = p0%octNb2%octNb2%octNb4       %octNb6       %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 2, 1) = p0%octNb1%octNb1%octNb4%octNb4%octNb6       %F(k)
             Ffld(k,-1, 2, 1) = p0%octNb1       %octNb4%octNb4%octNb6       %F(k)
             Ffld(k, 0, 2, 1) = p0              %octNb4%octNb4%octNb6       %F(k)
             Ffld(k, 1, 2, 1) = p0%octNb2       %octNb4%octNb4%octNb6       %F(k)
             Ffld(k, 2, 2, 1) = p0%octNb2%octNb2%octNb4%octNb4%octNb6       %F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2,-2, 2) = p0%octNb1%octNb1%octNb3%octNb3%octNb6%octNb6%F(k)
             Ffld(k,-1,-2, 2) = p0%octNb1       %octNb3%octNb3%octNb6%octNb6%F(k)
             Ffld(k, 0,-2, 2) = p0              %octNb3%octNb3%octNb6%octNb6%F(k)
             Ffld(k, 1,-2, 2) = p0%octNb2       %octNb3%octNb3%octNb6%octNb6%F(k)
             Ffld(k, 2,-2, 2) = p0%octNb2%octNb2%octNb3%octNb3%octNb6%octNb6%F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2,-1, 2) = p0%octNb1%octNb1%octNb3       %octNb6%octNb6%F(k)
             Ffld(k,-1,-1, 2) = p0%octNb1       %octNb3       %octNb6%octNb6%F(k)
             Ffld(k, 0,-1, 2) = p0              %octNb3       %octNb6%octNb6%F(k)
             Ffld(k, 1,-1, 2) = p0%octNb2       %octNb3       %octNb6%octNb6%F(k)
             Ffld(k, 2,-1, 2) = p0%octNb2%octNb2%octNb3       %octNb6%octNb6%F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 0, 2) = p0%octNb1%octNb1              %octNb6%octNb6%F(k)
             Ffld(k,-1, 0, 2) = p0%octNb1                     %octNb6%octNb6%F(k)
             Ffld(k, 0, 0, 2) = p0                            %octNb6%octNb6%F(k)
             Ffld(k, 1, 0, 2) = p0%octNb2                     %octNb6%octNb6%F(k)
             Ffld(k, 2, 0, 2) = p0%octNb2%octNb2              %octNb6%octNb6%F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 1, 2) = p0%octNb1%octNb1%octNb4       %octNb6%octNb6%F(k)
             Ffld(k,-1, 1, 2) = p0%octNb1       %octNb4       %octNb6%octNb6%F(k)
             Ffld(k, 0, 1, 2) = p0              %octNb4       %octNb6%octNb6%F(k)
             Ffld(k, 1, 1, 2) = p0%octNb2       %octNb4       %octNb6%octNb6%F(k)
             Ffld(k, 2, 1, 2) = p0%octNb2%octNb2%octNb4       %octNb6%octNb6%F(k)
!            --------------------------------------------------------------------
             Ffld(k,-2, 2, 2) = p0%octNb1%octNb1%octNb4%octNb4%octNb6%octNb6%F(k)
             Ffld(k,-1, 2, 2) = p0%octNb1       %octNb4%octNb4%octNb6%octNb6%F(k)
             Ffld(k, 0, 2, 2) = p0              %octNb4%octNb4%octNb6%octNb6%F(k)
             Ffld(k, 1, 2, 2) = p0%octNb2       %octNb4%octNb4%octNb6%octNb6%F(k)
             Ffld(k, 2, 2, 2) = p0%octNb2%octNb2%octNb4%octNb4%octNb6%octNb6%F(k)
          enddo
!
! -- pointer loop for particle --
!
          indexP =  p0%octP 
          pptc   => p0%ptcl%prtnxt
!
          do kp=1,indexP
!
            if(pptc%Isort.le.0) cycle
!
            Rchrg = Rcharge(pptc%Isort)
            Rmss  = Rmass(pptc%Isort)
!            time_factor    =  PI*QM*deltaT*Rchrg/Rmss
! toseo
            time_factor    = -0.25d0*QM*deltaT*Rchrg/Rmss
            charge_factor  = -dble(Rchrg)/dvlevel3  ! minus for current assignment
!
            do k=1,6
              R(k) = pptc%R(k)
            enddo

            pptc%R(7)=pptc%R(4)
            pptc%R(8)=pptc%R(5)
            pptc%R(9)=pptc%R(6)
!
! =====================================================
!    Form-factor calculation before particle movement
! =====================================================
!
! -- compute 1D form-factors --
            do k=1,3
!              gamma    = ( R(k) - R_lim(0,k) )*odxx(k) + HALF
              gamma    = R(k)*odxx(k) + HALF
!
! --- shifted grid (+1/2) ---
              h1       = gamma
              j        = NINT( h1 )        ! j=1 corresponds to x=0
              h1       = h1 - j
              j0(k)    = j
              h2       = h1*h1
              h3       = EIGHTH + HALF*h2
              h1       = HALF*h1
              W0(-1,k) = h3 - h1           ! 0.5*(0.5-h1)**2
              W0( 0,k) = THR_FOURTH - h2   ! 0.75 - h1**2
              W0( 1,k) = h3 + h1           ! 0.5*(0.5+h1)**2
              W0( 2,k) = ZERO
!
! ---  base grid (+1) ---
              h1       = gamma + HALF
              j        = NINT( h1 )        ! j=1 corresponds to x=0
              h1       = h1 - j
              j1(k)    = j
              h2       = h1*h1
              h3       = EIGHTH + HALF*h2
              h1       = HALF*h1
              W1(-1,k) = h3 - h1           ! 0.5*(0.5-h1)**2
              W1( 0,k) = THR_FOURTH - h2   ! 0.75 - h1**2
              W1( 1,k) = h3 + h1
            enddo
!
            do k=1,3
              j00(k) = j0(k)
              j11(k) = j1(k)
            enddo
!
            do k=1,3
              j0(k) = j0(k) - icel(k)
              j1(k) = j1(k) - icel(k)
            enddo
!
! =====================================================
!    Force calculation
! =====================================================
!
            do k=1,6
               FF(k) = (( Ffld(k,j1(1)-1,j0(2)-1,j0(3)-1) * W1(-1,1)             &
                        + Ffld(k,j1(1)  ,j0(2)-1,j0(3)-1) * W1( 0,1)             &
                        + Ffld(k,j1(1)+1,j0(2)-1,j0(3)-1) * W1( 1,1)) * W0(-1,2) &
                      + ( Ffld(k,j1(1)-1,j0(2)  ,j0(3)-1) * W1(-1,1)             &
                        + Ffld(k,j1(1)  ,j0(2)  ,j0(3)-1) * W1( 0,1)             &
                        + Ffld(k,j1(1)+1,j0(2)  ,j0(3)-1) * W1( 1,1)) * W0( 0,2) &
                      + ( Ffld(k,j1(1)-1,j0(2)+1,j0(3)-1) * W1(-1,1)             &
                        + Ffld(k,j1(1)  ,j0(2)+1,j0(3)-1) * W1( 0,1)             &
                        + Ffld(k,j1(1)+1,j0(2)+1,j0(3)-1) * W1( 1,1)) * W0( 1,2) &
                       ) * W0(-1,3)                                              &
!                -------------------------------------------------------------
                     + (( Ffld(k,j1(1)-1,j0(2)-1,j0(3)  ) * W1(-1,1)             &
                        + Ffld(k,j1(1)  ,j0(2)-1,j0(3)  ) * W1( 0,1)             &
                        + Ffld(k,j1(1)+1,j0(2)-1,j0(3)  ) * W1( 1,1)) * W0(-1,2) &
                      + ( Ffld(k,j1(1)-1,j0(2)  ,j0(3)  ) * W1(-1,1)             &
                        + Ffld(k,j1(1)  ,j0(2)  ,j0(3)  ) * W1( 0,1)             &
                        + Ffld(k,j1(1)+1,j0(2)  ,j0(3)  ) * W1( 1,1)) * W0( 0,2) &
                      + ( Ffld(k,j1(1)-1,j0(2)+1,j0(3)  ) * W1(-1,1)             &
                        + Ffld(k,j1(1)  ,j0(2)+1,j0(3)  ) * W1( 0,1)             &
                        + Ffld(k,j1(1)+1,j0(2)+1,j0(3)  ) * W1( 1,1)) * W0( 1,2) &
                       ) * W0( 0,3)                                              &
!                -------------------------------------------------------------
                     + (( Ffld(k,j1(1)-1,j0(2)-1,j0(3)+1) * W1(-1,1)             &
                        + Ffld(k,j1(1)  ,j0(2)-1,j0(3)+1) * W1( 0,1)             &
                        + Ffld(k,j1(1)+1,j0(2)-1,j0(3)+1) * W1( 1,1)) * W0(-1,2) &
                      + ( Ffld(k,j1(1)-1,j0(2)  ,j0(3)+1) * W1(-1,1)             &
                        + Ffld(k,j1(1)  ,j0(2)  ,j0(3)+1) * W1( 0,1)             &
                        + Ffld(k,j1(1)+1,j0(2)  ,j0(3)+1) * W1( 1,1)) * W0( 0,2) &
                      + ( Ffld(k,j1(1)-1,j0(2)+1,j0(3)+1) * W1(-1,1)             &
                        + Ffld(k,j1(1)  ,j0(2)+1,j0(3)+1) * W1( 0,1)             &
                        + Ffld(k,j1(1)+1,j0(2)+1,j0(3)+1) * W1( 1,1)) * W0( 1,2) &
                       ) * W0( 1,3)
            enddo
!
! =====================================================
!    Calculation of particle movement
! =====================================================
!
            do k=1,6
              FF(k) = FF(k) * time_factor
            enddo
!
! ---  u^- = u^{n-1/2} + (qdt/2m)*E  ----
            u1 = R(4) + FF(1)
            u2 = R(5) + FF(2)
            u3 = R(6) + FF(3)
!
! ---  gamma = sqrt{ 1 + |u^-|^2 }  ---
            gamma = SQRT( ONE + u1*u1 + u2*u2 + u3*u3 )
!
! ---  h = (qdt/2m)*(1/gamma)*B  ---
            h1 = FF(4)/gamma
            h2 = FF(5)/gamma
            h3 = FF(6)/gamma
!
! ---  coefficients for u^{n+1/2}  ---
            gamma = h1*h1 + h2*h2 + h3*h3
            c1    = (ONE - gamma)
            c2    = ONE/(ONE + gamma)
            gamma = h1*u1 + h2*u2 + h3*u3 
!
! ---  update velocities
!      u^{n+1/2} = (1-h^2)/(1+h^2)*u^-
!                + 2/(1+h^2)*( (h,u^-) + [u^- x h] ) + (qdt/2m)*E
            v1 = c2*(c1*u1 + 2*(gamma*h1 + u2*h3-u3*h2) ) + FF(1)
            v2 = c2*(c1*u2 + 2*(gamma*h2 + u3*h1-u1*h3) ) + FF(2)
            v3 = c2*(c1*u3 + 2*(gamma*h3 + u1*h2-u2*h1) ) + FF(3)
!
            pptc%R(4)= v1
            pptc%R(5)= v2
            pptc%R(6)= v3
!
! ---  update coordinates
!      x^{n+1} = x^{n} + u^{n+1/2}*dt/gamma^{n+1/2}
            gamma = deltaT/SQRT( ONE + v1*v1 + v2*v2 + v3*v3 )
!
            R(1) = R(1) + gamma*v1
            R(2) = R(2) + gamma*v2
            R(3) = R(3) + gamma*v3  
!
!            pptc%R(1) = R(1)
!            pptc%R(2) = R(2)
!            pptc%R(3) = R(3)
!
! =====================================================
!    Form-factor calculation after particle movement
! =====================================================
!
! -- Initialization for summation --
            do k=1,3
              do ik=-1,2
                do ij=-1,2
                  do ii=-1,2
                    Cfld(k,ii,ij,ik) = ZERO
                  enddo
                enddo
              enddo
            enddo
!
! -- compute new 1D form-factors --
!
            do k=1,3
              W1(-1,k) = ZERO
              W1( 2,k) = ZERO
            enddo
!
            do k=1,3
!              h1        = ( R(k) - R_lim(0,k) )*odxx(k) + HALF
              h1        =  R(k)*odxx(k) + HALF
!
              j         = NINT( h1 )            ! j=1 corresponds to x=0
              h1        = h1 - j
              jd        = j - j00(k)            ! absolute shift of a particle = -1,0,1
              j0(k)     = j00(k) + jd*(1-jd)/2  ! *OR* = j0(k) + MIN(jd,0)
              j         = jd*(1+jd)/2 - 1       ! *OR* = MAX(jd,0) - 1
              h2        = h1*h1
              h3        = EIGHTH + HALF*h2
              h1        = HALF*h1
!
              W1(j  ,k) = h3 - h1              ! 0.5*(0.5-h1)**2
              W1(j+1,k) = THR_FOURTH - h2      ! 0.75 - h1**2
              W1(j+2,k) = h3 + h1              ! 0.5*(0.5+h1)**2
!
!
              if (jd < 0) then
                W0( 2,k) = W0( 1,k)
                W0( 1,k) = W0( 0,k)
                W0( 0,k) = W0(-1,k)
                W0(-1,k) = ZERO                ! W0(-2,k)
              end if
!
            enddo
!
            do k=1,3
              j0(k) = j0(k) - icel(k)
            enddo
!
! ---  Compute changes of 1D form-factors  ---
            do k=1,3
              W1(-1,k) = W1(-1,k) - W0(-1,k)
              W1( 0,k) = W1( 0,k) - W0( 0,k)
              W1( 1,k) = W1( 1,k) - W0( 1,k)
              W1( 2,k) = W1( 2,k) - W0( 2,k)
            enddo
!
! ==========================================================
!    Compute the decomposition of the total density change
! ==========================================================
! -- Contribute to current density of corresponding oct --
! -- Dont forget to multiply the current by cdx(*) in the Maxwell solver! --
!
            do i3=-1,2
              h1 = charge_factor*(      W0(i3,3) + HALF *W1(i3,3) )
              h2 = charge_factor*( HALF*W0(i3,3) + THIRD*W1(i3,3) )
              do i2=-1,2
                h3 = h1*W0(i2,2) + h2*W1(i2,2)
                Cfld(1,0,i2,i3) = h3 * W1(-1,1)
                Cfld(1,1,i2,i3) = h3*( W1(-1,1) + W1(0,1) )
                Cfld(1,2,i2,i3) = h3*( W1(-1,1) + W1(0,1) + W1(1,1) )
              enddo
            enddo
!           ----------------------------------------------------------
            do i3=-1,2
              h1 = charge_factor*(      W0(i3,3) + HALF *W1(i3,3) )
              h2 = charge_factor*( HALF*W0(i3,3) + THIRD*W1(i3,3) )
              do i1=-1,2
                h3 = h1*W0(i1,1) + h2*W1(i1,1)
                Cfld(2,i1,0,i3) = h3 * W1(-1,2)
                Cfld(2,i1,1,i3) = h3*( W1(-1,2) + W1(0,2) )
                Cfld(2,i1,2,i3) = h3*( W1(-1,2) + W1(0,2) + W1(1,2) )
              enddo
            enddo
!           ----------------------------------------------------------
            do i2=-1,2
              h1 = charge_factor*(      W0(i2,2) + HALF *W1(i2,2) )
              h2 = charge_factor*( HALF*W0(i2,2) + THIRD*W1(i2,2) )
              do i1=-1,2
                h3 = h1*W0(i1,1) + h2*W1(i1,1)
                Cfld(3,i1,i2,0) = h3 * W1(-1,3)
                Cfld(3,i1,i2,1) = h3*( W1(-1,3) + W1(0,3) )
                Cfld(3,i1,i2,2) = h3*( W1(-1,3) + W1(0,3) + W1(1,3) )
              enddo
            enddo
!
! -- Pointer for target oct for summation of current --
            ptmp => p0
!
!  -- X-direction --
            if(j0(1) == -1) then
              ptmp => ptmp%octNb1
            elseif(j0(1) == 1) then
              ptmp => ptmp%octNb2
            else
              ptmp => ptmp
            endif
!  -- Y-direction --
            if(j0(2) == -1) then
              ptmp => ptmp%octNb3
            elseif(j0(2) == 1) then
              ptmp => ptmp%octNb4
            else
              ptmp => ptmp
            endif
!  -- Z-direction --
            if(j0(3) == -1) then
              ptmp => ptmp%octNb5
            elseif(j0(3) == 1) then
              ptmp => ptmp%octNb6
            else
              ptmp => ptmp
            endif
!
! ---  Sum-up current density  ---
!            k1=1+(kk-1)*3
!            k3=k1+2
            do k=1,3
!              k0 = k+(kk-1)*3
              k0=k
!!            ptmp%----  X  ----%----  Y  ----%----  Z  ----%C(k0) ------------------
              ptmp%octNb1       %octNb3       %octNb5       %C(k0) = Cfld(k,-1,-1,-1) &
            + ptmp%octNb1       %octNb3       %octNb5       %C(k0)
              ptmp              %octNb3       %octNb5       %C(k0) = Cfld(k, 0,-1,-1) &
            + ptmp              %octNb3       %octNb5       %C(k0)
              ptmp%octNb2       %octNb3       %octNb5       %C(k0) = Cfld(k, 1,-1,-1) &
            + ptmp%octNb2       %octNb3       %octNb5       %C(k0)
              ptmp%octNb2%octNb2%octNb3       %octNb5       %C(k0) = Cfld(k, 2,-1,-1) &
            + ptmp%octNb2%octNb2%octNb3       %octNb5       %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1                     %octNb5       %C(k0) = Cfld(k,-1, 0,-1) &
            + ptmp%octNb1                     %octNb5       %C(k0)
              ptmp                            %octNb5       %C(k0) = Cfld(k, 0, 0,-1) &
            + ptmp                            %octNb5       %C(k0)
              ptmp%octNb2                     %octNb5       %C(k0) = Cfld(k, 1, 0,-1) &
            + ptmp%octNb2                     %octNb5       %C(k0)
              ptmp%octNb2%octNb2              %octNb5       %C(k0) = Cfld(k, 2, 0,-1) &
            + ptmp%octNb2%octNb2              %octNb5       %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb4       %octNb5       %C(k0) = Cfld(k,-1, 1,-1) &
            + ptmp%octNb1       %octNb4       %octNb5       %C(k0)
              ptmp              %octNb4       %octNb5       %C(k0) = Cfld(k, 0, 1,-1) &
            + ptmp              %octNb4       %octNb5       %C(k0)
              ptmp%octNb2       %octNb4       %octNb5       %C(k0) = Cfld(k, 1, 1,-1) &
            + ptmp%octNb2       %octNb4       %octNb5       %C(k0)
              ptmp%octNb2%octNb2%octNb4       %octNb5       %C(k0) = Cfld(k, 2, 1,-1) &
            + ptmp%octNb2%octNb2%octNb4       %octNb5       %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb4%octNb4%octNb5       %C(k0) = Cfld(k,-1, 2,-1) &
            + ptmp%octNb1       %octNb4%octNb4%octNb5       %C(k0)
              ptmp              %octNb4%octNb4%octNb5       %C(k0) = Cfld(k, 0, 2,-1) &
            + ptmp              %octNb4%octNb4%octNb5       %C(k0)
              ptmp%octNb2       %octNb4%octNb4%octNb5       %C(k0) = Cfld(k, 1, 2,-1) &
            + ptmp%octNb2       %octNb4%octNb4%octNb5       %C(k0)
              ptmp%octNb2%octNb2%octNb4%octNb4%octNb5       %C(k0) = Cfld(k, 2, 2,-1) &
            + ptmp%octNb2%octNb2%octNb4%octNb4%octNb5       %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb3                     %C(k0) = Cfld(k,-1,-1, 0) &
            + ptmp%octNb1       %octNb3                     %C(k0)
              ptmp              %octNb3                     %C(k0) = Cfld(k, 0,-1, 0) &
            + ptmp              %octNb3                     %C(k0)
              ptmp%octNb2       %octNb3                     %C(k0) = Cfld(k, 1,-1, 0) &
            + ptmp%octNb2       %octNb3                     %C(k0)
              ptmp%octNb2%octNb2%octNb3                     %C(k0) = Cfld(k, 2,-1, 0) &
            + ptmp%octNb2%octNb2%octNb3                     %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1                                   %C(k0) = Cfld(k,-1, 0, 0) &
            + ptmp%octNb1                                   %C(k0)
              ptmp                                          %C(k0) = Cfld(k, 0, 0, 0) &
            + ptmp                                          %C(k0)
              ptmp%octNb2                                   %C(k0) = Cfld(k, 1, 0, 0) &
            + ptmp%octNb2                                   %C(k0)
              ptmp%octNb2%octNb2                            %C(k0) = Cfld(k, 2, 0, 0) &
            + ptmp%octNb2%octNb2                            %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb4                     %C(k0) = Cfld(k,-1, 1, 0) &
            + ptmp%octNb1       %octNb4                     %C(k0)
              ptmp              %octNb4                     %C(k0) = Cfld(k, 0, 1, 0) &
            + ptmp              %octNb4                     %C(k0)
              ptmp%octNb2       %octNb4                     %C(k0) = Cfld(k, 1, 1, 0) &
            + ptmp%octNb2       %octNb4                     %C(k0)
              ptmp%octNb2%octNb2%octNb4                     %C(k0) = Cfld(k, 2, 1, 0) &
            + ptmp%octNb2%octNb2%octNb4                     %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb4%octNb4              %C(k0) = Cfld(k,-1, 2, 0) &
            + ptmp%octNb1       %octNb4%octNb4              %C(k0)
              ptmp              %octNb4%octNb4              %C(k0) = Cfld(k, 0, 2, 0) &
            + ptmp              %octNb4%octNb4              %C(k0)
              ptmp%octNb2       %octNb4%octNb4              %C(k0) = Cfld(k, 1, 2, 0) &
            + ptmp%octNb2       %octNb4%octNb4              %C(k0)
              ptmp%octNb2%octNb2%octNb4%octNb4              %C(k0) = Cfld(k, 2, 2, 0) &
            + ptmp%octNb2%octNb2%octNb4%octNb4              %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb3       %octNb6       %C(k0) = Cfld(k,-1,-1, 1) &
            + ptmp%octNb1       %octNb3       %octNb6       %C(k0)
              ptmp              %octNb3       %octNb6       %C(k0) = Cfld(k, 0,-1, 1) &
            + ptmp              %octNb3       %octNb6       %C(k0)
              ptmp%octNb2       %octNb3       %octNb6       %C(k0) = Cfld(k, 1,-1, 1) &
            + ptmp%octNb2       %octNb3       %octNb6       %C(k0)
              ptmp%octNb2%octNb2%octNb3       %octNb6       %C(k0) = Cfld(k, 2,-1, 1) &
            + ptmp%octNb2%octNb2%octNb3       %octNb6       %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1                     %octNb6       %C(k0) = Cfld(k,-1, 0, 1) &
            + ptmp%octNb1                     %octNb6       %C(k0)
              ptmp                            %octNb6       %C(k0) = Cfld(k, 0, 0, 1) &
            + ptmp                            %octNb6       %C(k0)
              ptmp%octNb2                     %octNb6       %C(k0) = Cfld(k, 1, 0, 1) &
            + ptmp%octNb2                     %octNb6       %C(k0)
              ptmp%octNb2%octNb2              %octNb6       %C(k0) = Cfld(k, 2, 0, 1) &
            + ptmp%octNb2%octNb2              %octNb6       %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb4       %octNb6       %C(k0) = Cfld(k,-1, 1, 1) &
            + ptmp%octNb1       %octNb4       %octNb6       %C(k0)
              ptmp              %octNb4       %octNb6       %C(k0) = Cfld(k, 0, 1, 1) &
            + ptmp              %octNb4       %octNb6       %C(k0)
              ptmp%octNb2       %octNb4       %octNb6       %C(k0) = Cfld(k, 1, 1, 1) &
            + ptmp%octNb2       %octNb4       %octNb6       %C(k0)
              ptmp%octNb2%octNb2%octNb4       %octNb6       %C(k0) = Cfld(k, 2, 1, 1) &
            + ptmp%octNb2%octNb2%octNb4       %octNb6       %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb4%octNb4%octNb6       %C(k0) = Cfld(k,-1, 2, 1) &
            + ptmp%octNb1       %octNb4%octNb4%octNb6       %C(k0)
              ptmp              %octNb4%octNb4%octNb6       %C(k0) = Cfld(k, 0, 2, 1) &
            + ptmp              %octNb4%octNb4%octNb6       %C(k0)
              ptmp%octNb2       %octNb4%octNb4%octNb6       %C(k0) = Cfld(k, 1, 2, 1) &
            + ptmp%octNb2       %octNb4%octNb4%octNb6       %C(k0)
              ptmp%octNb2%octNb2%octNb4%octNb4%octNb6       %C(k0) = Cfld(k, 2, 2, 1) &
            + ptmp%octNb2%octNb2%octNb4%octNb4%octNb6       %C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb3       %octNb6%octNb6%C(k0) = Cfld(k,-1,-1, 2) &
            + ptmp%octNb1       %octNb3       %octNb6%octNb6%C(k0)
              ptmp              %octNb3       %octNb6%octNb6%C(k0) = Cfld(k, 0,-1, 2) &
            + ptmp              %octNb3       %octNb6%octNb6%C(k0)
              ptmp%octNb2       %octNb3       %octNb6%octNb6%C(k0) = Cfld(k, 1,-1, 2) &
            + ptmp%octNb2       %octNb3       %octNb6%octNb6%C(k0)
              ptmp%octNb2%octNb2%octNb3       %octNb6%octNb6%C(k0) = Cfld(k, 2,-1, 2) &
            + ptmp%octNb2%octNb2%octNb3       %octNb6%octNb6%C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1                     %octNb6%octNb6%C(k0) = Cfld(k,-1, 0, 2) &
            + ptmp%octNb1                     %octNb6%octNb6%C(k0)
              ptmp                            %octNb6%octNb6%C(k0) = Cfld(k, 0, 0, 2) &
            + ptmp                            %octNb6%octNb6%C(k0)
              ptmp%octNb2                     %octNb6%octNb6%C(k0) = Cfld(k, 1, 0, 2) &
            + ptmp%octNb2                     %octNb6%octNb6%C(k0)
              ptmp%octNb2%octNb2              %octNb6%octNb6%C(k0) = Cfld(k, 2, 0, 2) &
            + ptmp%octNb2%octNb2              %octNb6%octNb6%C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb4       %octNb6%octNb6%C(k0) = Cfld(k,-1, 1, 2) &
            + ptmp%octNb1       %octNb4       %octNb6%octNb6%C(k0)
              ptmp              %octNb4       %octNb6%octNb6%C(k0) = Cfld(k, 0, 1, 2) &
            + ptmp              %octNb4       %octNb6%octNb6%C(k0)
              ptmp%octNb2       %octNb4       %octNb6%octNb6%C(k0) = Cfld(k, 1, 1, 2) &
            + ptmp%octNb2       %octNb4       %octNb6%octNb6%C(k0)
              ptmp%octNb2%octNb2%octNb4       %octNb6%octNb6%C(k0) = Cfld(k, 2, 1, 2) &
            + ptmp%octNb2%octNb2%octNb4       %octNb6%octNb6%C(k0)
!             --------------------------------------------------------------------
              ptmp%octNb1       %octNb4%octNb4%octNb6%octNb6%C(k0) = Cfld(k,-1, 2, 2) &
            + ptmp%octNb1       %octNb4%octNb4%octNb6%octNb6%C(k0)
              ptmp              %octNb4%octNb4%octNb6%octNb6%C(k0) = Cfld(k, 0, 2, 2) &
            + ptmp              %octNb4%octNb4%octNb6%octNb6%C(k0)
              ptmp%octNb2       %octNb4%octNb4%octNb6%octNb6%C(k0) = Cfld(k, 1, 2, 2) &
            + ptmp%octNb2       %octNb4%octNb4%octNb6%octNb6%C(k0)
              ptmp%octNb2%octNb2%octNb4%octNb4%octNb6%octNb6%C(k0) = Cfld(k, 2, 2, 2) &
            + ptmp%octNb2%octNb2%octNb4%octNb4%octNb6%octNb6%C(k0)
            enddo
!
            pptc => pptc%prtNxt  ! to next particle
!
          end do ! particle loop
        end do   ! oct loop
!        end do   ! parallel loop
!
        do index=minindex,maxindex
           p0 => Mesh(index) 
              do k=1,3
                 p0%F(k+6) =  (p0%C(k))*(-1.d0)
              end do
        end do
!
        return
!
end subroutine back_push2
!
!***************************************************************
!yagi 2011/12/26
      subroutine setup_model
        use param
        use init_mesh_size
        implicit none
        integer(kind=4)::icon,injct,iLv
        real(kind=8)  :: XinR,YinR!,XinL,YinL
        real(kind=8)  :: VXinRe,VYinRe,VXinRi,VYinRi!,VXinLe,VYinLe,VXinLi,VYinLi
        
        icon=-999
        injct=200
        XinR= R_lim(1,1)*(0.5d0+1.d0/dble(NXB))
        YinR= R_lim(1,2)*0.5
        VXinRe=0.03d0  ;   VYinRe= 0.00d0
        VXinRi=0.03d0  ;   VYinRi= 0.00d0
        
        
        call particle_injct(XinR,YinR,VXinRe,VYinRe,VXinRi,VYinRi,injct,icon)

        !send data to Goct
        do iLv=0,LvMax
           call refresh_Fields("E",iLv)
           call refresh_Fields("B",iLv)
	if(CIP == 1) then
           call refresh_Fields("G",iLv)
	end if
           call refresh_Fields("J",iLv)
           call refresh_particles(iLv)
        enddo

      end subroutine setup_model

!***************************************************************
      subroutine particle_injct(Xc,Yc,VXc,VYc,VXd,VYd,injct,icon)
!***************************************************************
        use mpi
        use message_passing_interface
        use param
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        use init_condition
        use mtmod
        implicit none
        integer(kind=4)               :: i,j,index,IR,IR1,IR2,IR3
        integer(kind=4)               :: icon,injct
        integer(kind=4):: npart
        real(kind=8)                  :: Xc,Yc
        real(kind=8)                  :: VXc,VYc,VXd,VYd
        real(kind=8),dimension(Niall*IonSorts) :: PEwork
        real(kind=8),dimension(Niall) ::PEinitX,PEinitY,PEinitZ,PIinitX,PIinitY,PIinitZ
        real(kind=8)                  :: ar,wavek0,ratio0,vx0,Ex0
        integer(kind=4)               :: indexR,ecount,icount
        real(kind=8),dimension(9)     :: R
        real(kind=8),dimension(6)     :: rArea
        real(kind=8),dimension(3)     :: ptclVelc
        real(kind=8),dimension(3)     :: VT
        real(kind=8),dimension(32)    :: distParam
        integer(kind=4),dimension(3)  :: distType
        integer(kind=4)               :: ptclpercell
        integer(kind=4),dimension(IonSorts) :: ptclType
        integer(kind=4),dimension(IonSorts)::injctNum
        real(kind=8)                  :: xx,yy,zz,phase
        integer(kind=4)               ::factor
        real(kind=8)                  ::inject

!--- Debug Variables
    !    integer(kind=4) :: nptcl
! -- pointers --
        type(oct), pointer   :: p0
        type(prtcl), pointer :: PrtList!,pp

        !unused variables
        icon=0
        injct=0
        VXd=0
        VXc=0
        VYc=0
        VYd=0
        Xc=0
        Yc=0
        !initialize
        indexR = 0
        j=0
        distParam=0.0d0
        distType=0
        injctNum=0
        ptclType=0
        ptclpercell=0
        wavek0 = 0.d0
        ratio0 = 0.d0
        vx0   = 0.d0
        VT = 0.d0
! At initial, all particles are inserted in Lv 0
        do i=1,Lvmax
           MaxIP(i)=0
        end do

        npart=0
        
        ar=0.04d0
        IR=3207+rank !don't add rank
        call sgrnd(IR)
        IR1=3207
        IR2=4217
        IR3=5227

        if(Model==1 .or. Model==3 .or. Model==4 .or.Model==7)then!don't use R_lim_local!
           rArea(1)=R_lim(0,1)
           rArea(2)=R_lim(0,2)
           rArea(3)=R_lim(0,3)
           rArea(4)=R_lim(1,1)
           rArea(5)=R_lim(1,2)
           rArea(6)=R_lim(1,3)
           distParam=0.d0
           distType=0
           ptclType(1)=1
           ptclType(2)=2
           ptclVelc(1)=0.0d0
           ptclVelc(2)=0.0d0
           ptclVelc(3)=0.0d0
           call setParticleSquare(rArea,&
                distType,distParam,&
                ptclType,npart_per_cell,ptclVelc,&
                IR+rank,injctNum)
        else if(model==6)then
           rArea(1)=R_lim(0,1)
           rArea(2)=R_lim(1,2)*0.40d0
           rArea(3)=R_lim(1,3)*0.40d0
           rArea(4)=R_lim(1,1)
           rArea(5)=R_lim(1,2)*0.60d0
           rArea(6)=R_lim(1,3)*0.60d0
           distParam=0.d0
           distType(1:3)=2
           ptclType(1)=1
           ptclType(2)=2
           ptclVelc(1)=0.0d0
           ptclVelc(2)=0.0d0
           ptclVelc(3)=0.0d0
           call setParticleSquare(rArea,&
                distType,distParam,&
                ptclType,npart_per_cell,ptclVelc,&
                IR+rank,injctNum)
        else if(model==8)then
           factor=int(wall/10)+1
           inject=(factor*10.d0)/(NXR*NXB)
           rArea(1)=R_lim(0,1)
           rArea(2)=R_lim(1,2)
           rArea(3)=R_lim(1,3)
           rArea(4)=R_lim(1,1)*inject
           rArea(5)=R_lim(1,2)
           rArea(6)=R_lim(1,3)
           distParam=0.d0
           distType(1:3)=2
           ptclType(1)=1
           ptclType(2)=2
           ptclVelc(1)=0.0d0
           ptclVelc(2)=0.0d0
           ptclVelc(3)=0.0d0
           call setParticleSquare(rArea,&
                distType,distParam,&
                ptclType,npart_per_cell,ptclVelc,&
                IR+rank,injctNum)        
        else if(model==2)then
           rArea(1)=R_lim(0,1)
           rArea(2)=R_lim(0,2)
           rArea(3)=R_lim(0,3)
           rArea(4)=R_lim(1,1)
           rArea(5)=R_lim(1,2)
           rArea(6)=R_lim(1,3)
           distParam=0.d0
           distType=0
           ptclpercell=Int(0.5*npart_per_cell)
           ptclType(1)=1
           ptclType(2)=2
           ptclVelc(1)=0.0d0
           ptclVelc(2)=0.0d0
           ptclVelc(3)=0.0d0
           call setParticleSquare(rArea,&
                distType,distParam,&
                ptclType,ptclpercell,ptclVelc,&
                IR+rank,injctNum)
           !-y cloud
           rArea(1)=R_lim(1,1)/2
           rArea(2)=R_lim(0,2)+R_lim(1,2)/4
           rArea(3)=R_lim(1,3)/2
           rArea(4)=R_lim(1,1)/8     
           distParam=0.d0
           distType=0
           ptclpercell=Int(0.5*npart_per_cell)
           ptclType(1)=1
           ptclType(2)=2
           ptclVelc=0.d0
           ptclVelc(1)=0.d0
           ptclVelc(2)=1.d0
           call setParticleSphere(rArea(1:3),rArea(4),&
                distType,distParam,&
                ptclType,ptclpercell,ptclVelc,&
                IR+rank,injctNum)    

           !+y cloud
           rArea(1)=R_lim(1,1)/2
           rArea(2)=R_lim(1,2)-R_lim(1,2)/4
           rArea(3)=R_lim(1,3)/2
           rArea(4)=R_lim(1,1)/8     
           distParam=0.d0
           distType=0
           ptclpercell=Int(0.5*npart_per_cell)
           ptclType(1)=1
           ptclType(2)=2
           ptclVelc=0.d0
           ptclVelc(1)=0.d0
           ptclVelc(2)=-1.d0
           call setParticleSphere(rArea(1:3),rArea(4),&
                distType,distParam,&
                ptclType,ptclpercell,ptclVelc,&
                IR+rank,injctNum) 

           !+x cloud
           rArea(1)=R_lim(0,1)+R_lim(1,1)/4
           rArea(2)=R_lim(1,2)/2
           rArea(3)=R_lim(1,3)/2
           rArea(4)=R_lim(1,1)/8     
           distParam=0.d0
           distType=0
           ptclpercell=Int(0.5*npart_per_cell)
           ptclType(1)=1
           ptclType(2)=2
           ptclVelc=0.d0
           ptclVelc(1)=1.d0
           ptclVelc(2)=0.d0
           call setParticleSphere(rArea(1:3),rArea(4),&
                distType,distParam,&
                ptclType,ptclpercell,ptclVelc,&
                IR+rank,injctNum) 

           !-x cloud
           rArea(1)=R_lim(1,1)-R_lim(1,1)/4
           rArea(2)=R_lim(1,2)/2
           rArea(3)=R_lim(1,3)/2
           rArea(4)=R_lim(1,1)/8     
           distParam=0.d0
           distType=0
           ptclpercell=Int(0.5*npart_per_cell)
           ptclType(1)=1
           ptclType(2)=2
           ptclVelc=0.d0
           ptclVelc(1)=-1.d0
           ptclVelc(2)=0.d0
           call setParticleSphere(rArea(1:3),rArea(4),&
                distType,distParam,&
                ptclType,ptclpercell,ptclVelc,&
                IR+rank,injctNum) 
        endif


!!$!!!!!!!!!!!!!!!!!! Structure of Pesh Array !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!               !-----Rep particles-----!     
!!$!     Pesh  !                           |  e  |   i   |   e   |   i   |
!!$!    index  1                        indexR                 j  
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$        if(debugMode>0)print *,"particle inject complete",injctNum,rank
!!$        if(debugMode>0)call output_octs_and_family(MinID(1,-1),MaxID(1,-1),LvMax,0,1,0)
!!$        if(debugMode>0)print*,'MaxIP is',MaxIP(0),rank
  

! Create Maxwell D
! Further deploy more particles with Maxwell Distribution function

!-- Particle coordinates loading finished
     

!------------------------------------------------------------------
! Record the number of particles per cell and store in octP
        do index=MinID(1,0), MaxID(1,0)
           indexR = 0 
           p0 => Mesh(index)
           PrtList => p0%ptcl

           if (associated(PrtList % prtnxt)) then
              do while(associated(PrtList))
                 if (PrtList % isort > 0 ) then
                    indexR = indexR+1
                 endif
                 Prtlist => PrtList%prtnxt
              end do
           endif
           p0%octP = indexR         
        end do
!----------------------------------------------------------------- 

        !if(debugMode>1)print *,"set octP complete",rank

 ! Invoke Maxwellian velocity distribution
       
        if(Model==1)then
           !if(debugMode>=3)print *,"invoke maxwellian",injctNum(1),maxIP(0),rank

           wavek0 = 1.d0*PI/(R_lim_local(1,1)-R_lim_local(0,1))
           ratio0 = omega/wavek0
           vx0   = ar*ratio0
           !print *,"QSTART start",rank
           if(TVELE.gt.0.d0) then 
              if(injctNum(1)+injctNum(2)>0)then
                 CALL QSTART( IR+rank, PEinitX, PEwork, injctNum(1)+injctNum(2), TVELE, 0.d0)
                 CALL QSTART( IR+rank, PEinitY, PEwork, injctNum(1)+injctNum(2), TVELE, 0.d0)
                 CALL QSTART( IR+rank, PEinitZ, PEwork, injctNum(1)+injctNum(2), TVELE, 0.d0)
              endif
           end if
           !print *,"QSTART ok",rank
           j=1
           do i=1,MaxIP(0)
              PrtList => Pesh(i,0)
              if(Prtlist%isort<=0)cycle 
              !print *,"j chosen=",j,rank
              R(1) = PrtList %R(1)
              R(2) = PrtList %R(2)
              R(3) = PrtList %R(3)
              R(4) = PEinitX(j)
              R(5) = PEinitY(j)
              R(6) = PEinitZ(j)
              PrtList%R(4) = R(4)+0.0
              PrtList%R(5) = R(5)! +vx0*dsin((R(2)-R_lim_local(0,2))*2.d0*PI/(R_lim_local(1,2)-R_lim_local(0,2)))
              PrtList%R(6) = R(6)
              PrtList%R(7) = R(4)
              PrtList%R(8) = R(5)
              PrtList%R(9) = R(6)
              j=j+1
           end do
           !if(debugMode>=3)print *,"maxwellian comp",rank
        else if(Model==3)then
           VT(1)=TVELE
           VT(2)=TVELE
           VT(3)=5.0d0*TVELE
           if(TVELE.gt.0.d0) then 
              if(injctNum(1)+injctNum(2)>0)then
                 CALL QSTART( IR+rank, PEinitX, PEwork, injctNum(1)+injctNum(2), VT(1), 0.d0)
                 CALL QSTART( IR+rank, PEinitY, PEwork, injctNum(1)+injctNum(2), VT(2), 0.d0)
                 CALL QSTART( IR+rank, PEinitZ, PEwork, injctNum(1)+injctNum(2), VT(3), 0.d0)
              endif
           end if
           j=1
           do i=1,MaxIP(0)
              PrtList => Pesh(i,0)
              if(Prtlist%isort<=0)cycle 
              R(1) = PrtList %R(1)
              R(2) = PrtList %R(2)
              R(3) = PrtList %R(3)
              R(4) = PEinitX(j)
              R(5) = PEinitY(j)
              R(6) = PEinitZ(j)
              PrtList%R(4) = R(4)+0.0
              PrtList%R(5) = R(5)
              PrtList%R(6) = R(6)
              PrtList%R(7) = R(4)
              PrtList%R(8) = R(5)
              PrtList%R(9) = R(6)
              j=j+1
           end do
        else if(Model==6)then
           !if(debugMode>=3)print *,"invoke maxwellian",injctNum(1),maxIP(0),rank
           !print *,"QSTART start",rank
           VT(1)=TVELE
           VT(2)=TVELE
           VT(3)=TVELE
           if(TVELE.gt.0.d0) then 
              if(injctNum(1)+injctNum(2)>0)then
                 call set_thermal_velocity(VT(1),PEinitX,Niall,injctNum(1)+injctNum(2))
                 call set_thermal_velocity(VT(2),PEinitY,Niall,injctNum(1)+injctNum(2))
                 call set_thermal_velocity(VT(3),PEinitZ,Niall,injctNum(1)+injctNum(2))
              endif
           end if
           !print *,"QSTART ok",rank
           j=1
           do i=1,MaxIP(0)
              PrtList => Pesh(i,0)

              if(Prtlist%isort<=0)cycle 
              !print *,"j chosen=",j,rank
              R(1) = PrtList %R(1)
              R(2) = PrtList %R(2)
              R(3) = PrtList %R(3)
              if (Prtlist%isort==1) then
                R(4) = PEinitX(j)*5.0d0
                R(5) = PEinitY(j)*5.0d0
                R(6) = PEinitZ(j)*5.0d0
              else
                R(4) = PEinitX(j)
                R(5) = PEinitY(j)
                R(6) = PEinitZ(j)
              end if
              PrtList%R(4) = R(4)
              PrtList%R(5) = R(5)! +vx0*dsin((R(2)-R_lim_local(0,2))*2.d0*PI/(R_lim_local(1,2)-R_lim_local(0,2)))
              PrtList%R(6) = R(6)
              PrtList%R(7) = R(4)
              PrtList%R(8) = R(5)! +vx0*dsin((R(2)-R_lim_local(0,2))*2.d0*PI/(R_lim_local(1,2)-R_lim_local(0,2)))
              PrtList%R(9) = R(6)
              j=j+1
           end do
           !if(debugMode>=3)print *,"maxwellian comp",rank
        else if(Model==4)then
           VT(1)=TVELE
           VT(2)=TVION
           if(TVELE.gt.0.d0) then 
                call set_thermal_velocity(VT(1),PEinitX,Niall,injctNum(1))
                call set_thermal_velocity(VT(1),PEinitY,Niall,injctNum(1))
                call set_thermal_velocity(VT(1),PEinitZ,Niall,injctNum(1))
                call set_thermal_velocity(VT(2),PIinitX,Niall,injctNum(2))
                call set_thermal_velocity(VT(2),PIinitY,Niall,injctNum(2))
                call set_thermal_velocity(VT(2),PIinitZ,Niall,injctNum(2))
           end if
           ecount=1
           icount=1
           do i=1,MaxIP(0)
              PrtList => Pesh(i,0)
              R(1) = PrtList %R(1)
              R(2) = PrtList %R(2)
              R(3) = PrtList %R(3)
              if(Prtlist%isort==1)then
                R(4) = PEinitX(ecount)+flow(1)
                R(5) = PEinitY(ecount)+flow(2)
                R(6) = PEinitZ(ecount)+flow(3)
                ecount=ecount+1
              else if(Prtlist%isort==2)then
                R(4) = PIinitX(icount)+flow(1)
                R(5) = PIinitY(icount)+flow(2)
                R(6) = PIinitZ(icount)+flow(3)
                icount=icount+1
              end if
                PrtList%R(4) = R(4)
                PrtList%R(5) = R(5)
                PrtList%R(6) = R(6)
                PrtList%R(7) = R(4)
                PrtList%R(8) = R(5)
                PrtList%R(9) = R(6)
           end do
        else if(Model==7)then
           VT(1)=TVELE
           VT(2)=TVION
           if(TVELE.gt.0.d0) then 
                call set_thermal_velocity(VT(1),PEinitX,Niall,injctNum(1))
                call set_thermal_velocity(VT(1),PEinitY,Niall,injctNum(1))
                call set_thermal_velocity(VT(1),PEinitZ,Niall,injctNum(1))
                call set_thermal_velocity(VT(2),PIinitX,Niall,injctNum(2))
                call set_thermal_velocity(VT(2),PIinitY,Niall,injctNum(2))
                call set_thermal_velocity(VT(2),PIinitZ,Niall,injctNum(2))
           end if
           ecount=1
           icount=1
           do i=1,MaxIP(0)
              PrtList => Pesh(i,0)
              R(1) = PrtList %R(1)
              R(2) = PrtList %R(2)
              R(3) = PrtList %R(3)
              if(Prtlist%isort==1)then
                R(4) = PEinitX(ecount)+flow(1)
                R(5) = PEinitY(ecount)+flow(2)
                R(6) = PEinitZ(ecount)+flow(3)
                ecount=ecount+1
              else if(Prtlist%isort==2)then
                R(4) = PIinitX(icount)+flow(1)
                R(5) = PIinitY(icount)+flow(2)
                R(6) = PIinitZ(icount)+flow(3)
                icount=icount+1
              end if
                PrtList%R(4) = R(4)
                PrtList%R(5) = R(5)
                PrtList%R(6) = R(6)
                PrtList%R(7) = R(4)
                PrtList%R(8) = R(5)
                PrtList%R(9) = R(6)
           end do   
        end if

        wavek0 = 1.d0*PI/(R_lim_local(1,1)-R_lim_local(0,1))
        ratio0 = omega/wavek0
        vx0   = ar*ratio0

        if(Model == 1) then
           Ex0 = ar*ratio0*omega
           BiniN=(/Ex0, 0.d0, 0.d0/)
           EiniN=(/0.d0, 0.d0, 0.d0/)
           kiniN=(/1.d0, 0.d0, 0.d0/)
           phi=0.d0
           do index=minID(1,0),maxID(1,0)
              p0 => Mesh(index)
              xx = p0 % rPOS(1)
              yy = p0 % rPOS(2)
              zz = p0 % rPOS(3)

              phase = pi2*(kiniN(1)*xx/(R_lim_local(1,1)-R_lim_local(0,1)) &
                   +kiniN(2)*yy/(R_lim_local(1,2)-R_lim_local(0,2)) &
                   +kiniN(3)*zz/(R_lim_local(1,3)-R_lim_local(0,3)))

              p0 % F(1:3) = 0d0 ! EiniN*dcos(phase+phi)
              p0 % F(4:6) = BiniN ! *dcos(phase+phi)
              p0 % F(7:12)= 0d0
           enddo
        else if(Model==3) then
           do index=minID(1,0),maxID(1,0)
              p0 => Mesh(index)
              p0 % F(1:3) = 0.d0
              p0 % F(4:6) = 0.d0
              p0 % F(7:12)= 0.d0
           enddo
        else if(Model==0)then
           BiniN=(/0.0d0, 0.0d0, 1.0d0/)
           EiniN=(/0.0d0, 1.0d0, 0.0d0/)
           kiniN=(/1.0d0, 0.0d0, 0.0d0/)
           phi=0.d0
           do index=minID(1,0),maxID(1,0)
              p0 => Mesh(index)
              xx = p0 % rPOS(1)
              yy = p0 % rPOS(2)
              zz = p0 % rPOS(3)

              phase = pi2*(kiniN(1)*xx/(R_lim_local(1,1)-R_lim_local(0,1)) &
                   +kiniN(2)*yy/(R_lim_local(1,2)-R_lim_local(0,2)) &
                   +kiniN(3)*zz/(R_lim_local(1,3)-R_lim_local(0,3)))
              !debug
              phase = pi2*(kiniN(1)*xx/(R_lim(1,1)-R_lim(0,1)) &
                   +kiniN(2)*yy/(R_lim(1,2)-R_lim(0,2)) &
                   +kiniN(3)*zz/(R_lim(1,3)-R_lim(0,3)))

              p0 % F(1:3) = EiniN*dsin(phase+phi)
              p0 % F(4:6) = BiniN*dsin(phase+phi)             
              p0 % F(7:12)=0d0
           enddo
        else if(Model==6)then
           BiniN=(/0.0d0, 0.0d0, 1.0d0/)
           EiniN=(/0.0d0, 1.0d0, 0.0d0/)
           kiniN=(/1.0d0, 0.0d0, 0.0d0/)
           phi=0.d0
           do index=minID(1,0),maxID(1,0)
              p0 => Mesh(index)
              xx = p0 % rPOS(1)
              yy = p0 % rPOS(2)
              zz = p0 % rPOS(3)

              phase = pi2*(kiniN(1)*xx/(R_lim_local(1,1)-R_lim_local(0,1)) &
                   +kiniN(2)*yy/(R_lim_local(1,2)-R_lim_local(0,2)) &
                   +kiniN(3)*zz/(R_lim_local(1,3)-R_lim_local(0,3)))
              !debug
              phase = pi2*(kiniN(1)*xx/(R_lim(1,1)-R_lim(0,1)) &
                   +kiniN(2)*yy/(R_lim(1,2)-R_lim(0,2)) &
                   +kiniN(3)*zz/(R_lim(1,3)-R_lim(0,3)))

              p0 % F(1:3) = 0.d0
              p0 % F(4:6) = 0.d0             
              p0 % F(7:9) = kiniN*dsin(phase+phi)
              p0 % F(10:12)=kiniN*dsin(phase+phi)
           enddo  
        else
           do index=minID(1,0),maxID(1,0)
              p0 => Mesh(index)

              p0 % F(1:3) = 0.d0
              p0 % F(4:6) = 0.d0             
              p0 % F(7:12)=0d0
           enddo
        endif
        return 
        
      contains
           !yagi 2012/02/13
subroutine injctParticle(R,Isort)
          implicit none
          integer(kind=4),intent(in)::Isort
          real(kind=8),dimension(9),intent(in)::R
          integer(kind=4),dimension(3)::iPosTemp
          integer(kind=4)::Mn,IPinit,j
          integer(kind=4)::chdx,chdy,chdz,iC
          type(oct),pointer::p0
          type(prtcl),pointer::PrtList

          
          MaxIP(0)=MaxIP(0)+1
          j=MaxIP(0)
          nullify(PrtList)

          if(j>=PeshSize)then
             print '(A,I3,A,I12)',"[injctParticle] In rank=",rank," :Pesh is full MaxIP(0)=",MaxIP(0)
             stop
          endif

          iPosTemp(1)=int(R(1)/(dx(1)*2))+1
          iPosTemp(2)=int(R(2)/(dx(2)*2))+1
          iPosTemp(3)=int(R(3)/(dx(3)*2))+1      
          call get_indexNumberN(iPosTemp,Mn)

          p0=>Mesh(Mn-MnDisp)

          chdx=int((R(1)-p0%rPos(1))/dx(1)+1)
          chdy=int((R(2)-p0%rPos(2))/dx(2)+1)
          chdz=int((R(3)-p0%rPos(3))/dx(3)+1)
          iC=4*chdz+2*chdy+chdx+1

          if(iC<1 .or. iC>8)then
             print *,"[injctParticle] particle position error! ",rank,iC,iposTemp,R(1:3)
             stop
          endif

          if(iC<=4)then
             if(iC<=2)then
                if(iC==1)then 
                   p0 => p0%octCh1
                else
                   p0 => p0%octCh2
                endif
             else
                if(iC==3)then
                   p0 => p0%octCh3
                else
                   p0 => p0%octCh4
                endif
             endif
          else
             if(iC<=6)then
                if(iC==5)then
                   p0 => p0%octCh5
                else
                   p0 => p0%octCh6
                endif
             else
                if(iC==7)then
                   p0 => p0%octCh7
                else
                   p0 => p0%octCh8
                endif
             endif
          endif

          IPinit=p0%octN
          PrtList=>p0%ptcl%prtnxt
          Pesh(j,0) = prtcl(R,isort,IPinit,j,PrtList,0)
          p0%ptcl%prtnxt=> Pesh(j,0)

end subroutine injctParticle

subroutine setParticleSquare(rArea,ptclDistType,distParam,ptclType,ptclpercell,ptclVelc,IR,injctNum)
          use message_passing_interface
          use init_mesh_size
          implicit none
          !define the area to inject particles.
          real(kind=8)   ,intent(in)::rArea(6),distParam(32),ptclVelc(3)
          !if the value of rArea was 0, it would be neglected. 
          !if the value of rArea was negative, particles would set in whole region. 
          integer(kind=4),intent(in)::ptclDistType(3),ptclType(IonSorts),ptclpercell,IR
          !16 types of particles can be injected at once by same position. 
          !particle distribution type = 0 : white noise
          !                           = 1 : sine random
          integer(kind=4),intent(inout),dimension(IonSorts)::injctNum
          integer(kind=4)::ptclNum,octNum,num,start

          real(kind=8)::rPosS(3),rPosE(3),rPosLS(3),rPosLE(3)
          real(kind=8),dimension(ptclpercell)::PinitX,PinitY,PinitZ
          real(kind=8),dimension(ptclpercell/8)::temp
          real(kind=8),dimension(ptclpercell*NXR*NYR*NZR*NX*NY*NZ)::randomX,randomY,randomZ
          real(kind=8),dimension(ptclpercell)::PEwork
          real(kind=8),dimension(9)::R
          integer(kind=4)::i,j,k
          type(oct),pointer::p0

          !if(debugMode>0)print *,"Starting setParticleSquare",rank
          rPosS=rArea(1:3)
          rPosE=rArea(4:6)
          octNum=0

          !modify inputs
          if(rArea(1)<0)rPosS(1)=R_lim(0,1)
          if(rArea(2)<0)rPosS(2)=R_lim(0,2)
          if(rArea(3)<0)rPosS(3)=R_lim(0,3)
          if(rArea(4)<0)rPosE(1)=R_lim(1,1)
          if(rArea(5)<0)rPosE(2)=R_lim(1,2)
          if(rArea(6)<0)rPosE(3)=R_lim(1,3)
          
          !confine injecting area in region of myrank
          if(rPosS(1)<R_lim_local(0,1))rPosS(1)=R_lim_local(0,1)
          if(rPosS(2)<R_lim_local(0,2))rPosS(2)=R_lim_local(0,2)
          if(rPosS(3)<R_lim_local(0,3))rPosS(3)=R_lim_local(0,3)
          if(rPosS(1)>R_lim_local(1,1))return
          if(rPosS(2)>R_lim_local(1,2))return
          if(rPosS(3)>R_lim_local(1,3))return

          if(rPosE(1)>R_lim_local(1,1))rPosE(1)=R_lim_local(1,1)
          if(rPosE(2)>R_lim_local(1,2))rPosE(2)=R_lim_local(1,2)
          if(rPosE(3)>R_lim_local(1,3))rPosE(3)=R_lim_local(1,3)
          if(rPosE(1)<R_lim_local(0,1))return
          if(rPosE(2)<R_lim_local(0,2))return
          if(rPosE(3)<R_lim_local(0,3))return

          if(ptclDistType(1)==2)then
             CALL setRandom(12345,randomX, ptclpercell*NXR*NYR*NZR*NX*NY*NZ, 0.0d0, 1.0d0)
             CALL setRandom(37567,randomY, ptclpercell*NXR*NYR*NZR*NX*NY*NZ, 0.0d0, 1.0d0)
             CALL setRandom(56789,randomZ, ptclpercell*NXR*NYR*NZR*NX*NY*NZ, 0.0d0, 1.0d0)
          endif

          do index=MinID(1,0),MaxID(1,0)
             p0=>Mesh(index)
             if(p0%rPos(1)<rPosS(1) .or. p0%rPos(2)<rPosS(2) .or. p0%rPos(3)<rPosS(3) .or.&
                p0%rPos(1)>rPosE(1) .or. p0%rPos(2)>rPosE(2) .or. p0%rPos(3)>rPosE(3))cycle
             octNum=octNum+1
             ptclNum=ptclpercell

             do i=1,IonSorts
                if(ptclType(i)<=0)cycle
                injctNum(ptclType(i))=injctNum(ptclType(i))+ptclNum
                if(injctNum(ptclType(i))>Niall)then
                   print *,"[Error]setParticleSphere, sum of the ptcl num isort=",ptclType,rank
                   print *,"[Error]violates Niall. ptclNum=",injctNum(ptclType(i)),"Niall=",Niall,rank
                   stop
                endif
             enddo
             rPosLS=p0%rPos-dx*HALF
             rPosLE=p0%rPOs+dx*HALF
           
             rPosLS=rPosLS+0.00000001d0
             rPosLE=rPosLE-0.00000001d0
             if(ptclDistType(1)==0)then
                CALL RPBMP1( IR+p0%MrtN, PinitX, PEwork, ptclNum, rPosLS(1), rPosLE(1), 5)
             elseif(ptclDistType(1)==1)then
                CALL SINRAN( IR+p0%MrtN, PinitX, PEwork, ptclNum, rPosLE(1), rPosLS(1),distParam(1))
             elseif(ptclDistType(1)==2)then
                num = ptclNum/8
                start = index-BMeshBound+Nall_ini*rank
CALL setAxis(randomX, temp, num, 1+ptclNum*(start-1), num+ptclNum*(start-1), rPosLE(1)-dx(1)*HALF, rPosLS(1))
do k=1,num
  PinitX(k) = temp(k)
enddo
CALL setAxis(randomX, temp, num, 1+ptclNum*(start-1)+num, num*2+ptclNum*(start-1), rPosLE(1), rPosLS(1)+dx(1)*HALF)
do k=1,num
  PinitX(k+num) = temp(k)
enddo
CALL setAxis(randomX, temp, num, 1+ptclNum*(start-1)+num*2, num*3+ptclNum*(start-1), rPosLE(1)-dx(1)*HALF, rPosLS(1))
do k=1,num
  PinitX(k+num*2) = temp(k)
enddo
CALL setAxis(randomX, temp, num, 1+ptclNum*(start-1)+num*3, num*4+ptclNum*(start-1), rPosLE(1), rPosLS(1)+dx(1)*HALF)
do k=1,num
  PinitX(k+num*3) = temp(k)
enddo
CALL setAxis(randomX, temp, num, 1+ptclNum*(start-1)+num*4, num*5+ptclNum*(start-1), rPosLE(1)-dx(1)*HALF, rPosLS(1))
do k=1,num
  PinitX(k+num*4) = temp(k)
enddo
CALL setAxis(randomX, temp, num, 1+ptclNum*(start-1)+num*5, num*6+ptclNum*(start-1), rPosLE(1), rPosLS(1)+dx(1)*HALF)
do k=1,num
  PinitX(k+num*5) = temp(k)
enddo
CALL setAxis(randomX, temp, num, 1+ptclNum*(start-1)+num*6, num*7+ptclNum*(start-1), rPosLE(1)-dx(1)*HALF, rPosLS(1))
do k=1,num
  PinitX(k+num*6) = temp(k)
enddo
CALL setAxis(randomX, temp, num, 1+ptclNum*(start-1)+num*7, ptclNum+ptclNum*(start-1), rPosLE(1), rPosLS(1)+dx(1)*HALF)
do k=1,num
  PinitX(k+num*7) = temp(k)
enddo

             endif
             if(ptclDistType(2)==0)then
                CALL RPBMP1( IR+p0%MrtN+65789, PinitY, PEwork, ptclNum, rPosLS(2), rPosLE(2), 5)
             elseif(ptclDistType(2)==1)then
                CALL SINRAN( IR+p0%MrtN, PinitY, PEwork, ptclNum, rPosLE(2), rPosLS(2),distParam(1))
             elseif(ptclDistType(2)==2)then
                num = ptclNum/8
                start = index-BMeshBound+Nall_ini*rank
CALL setAxis(randomY, temp, num, 1+ptclNum*(start-1), num+ptclNum*(start-1), rPosLE(2)-dx(2)*HALF, rPosLS(2))
do k=1,num
  PinitY(k) = temp(k)
enddo
CALL setAxis(randomY, temp, num, 1+ptclNum*(start-1)+num, num*2+ptclNum*(start-1), rPosLE(2)-dx(2)*HALF, rPosLS(2))
do k=1,num
  PinitY(k+num) = temp(k)
enddo
CALL setAxis(randomY, temp, num, 1+ptclNum*(start-1)+num*2, num*3+ptclNum*(start-1), rPosLE(2), rPosLS(2)+dx(2)*HALF)
do k=1,num
  PinitY(k+num*2) = temp(k)
enddo
CALL setAxis(randomY, temp, num, 1+ptclNum*(start-1)+num*3, num*4+ptclNum*(start-1), rPosLE(2), rPosLS(2)+dx(2)*HALF)
do k=1,num
  PinitY(k+num*3) = temp(k)
enddo
CALL setAxis(randomY, temp, num, 1+ptclNum*(start-1)+num*4, num*5+ptclNum*(start-1), rPosLE(2)-dx(2)*HALF, rPosLS(2))
do k=1,num
  PinitY(k+num*4) = temp(k)
enddo
CALL setAxis(randomY, temp, num, 1+ptclNum*(start-1)+num*5, num*6+ptclNum*(start-1), rPosLE(2)-dx(2)*HALF, rPosLS(2))
do k=1,num
  PinitY(k+num*5) = temp(k)
enddo
CALL setAxis(randomY, temp, num, 1+ptclNum*(start-1)+num*6, num*7+ptclNum*(start-1), rPosLE(2), rPosLS(2)+dx(2)*HALF)
do k=1,num
  PinitY(k+num*6) = temp(k)
enddo
CALL setAxis(randomY, temp, num, 1+ptclNum*(start-1)+num*7, ptclNum+ptclNum*(start-1), rPosLE(2), rPosLS(2)+dx(2)*HALF)
do k=1,num
  PinitY(k+num*7) = temp(k)
enddo
             endif
             if(ptclDistType(3)==0)then
                CALL RPBMP1( IR+p0%MrtN+37567, PinitZ, PEwork, ptclNum, rPosLS(3), rPosLE(3), 5)
             elseif(ptclDistType(3)==1)then
                CALL SINRAN( IR+p0%MrtN, PinitZ, PEwork, ptclNum, rPosLE(3), rPosLS(3),distParam(1))
             elseif(ptclDistType(3)==2)then
                num = ptclNum/8
                start = index-BMeshBound+Nall_ini*rank
CALL setAxis(randomZ, temp, num, 1+ptclNum*(start-1), num+ptclNum*(start-1), rPosLE(3)-dx(3)*HALF, rPosLS(3))
do k=1,num
  PinitZ(k) = temp(k)
enddo
CALL setAxis(randomZ, temp, num, 1+ptclNum*(start-1)+num, num*2+ptclNum*(start-1), rPosLE(3)-dx(3)*HALF, rPosLS(3))
do k=1,num
  PinitZ(k+num) = temp(k)
enddo
CALL setAxis(randomZ, temp, num, 1+ptclNum*(start-1)+num*2, num*3+ptclNum*(start-1), rPosLE(3)-dx(3)*HALF, rPosLS(3))
do k=1,num
  PinitZ(k+num*2) = temp(k)
enddo
CALL setAxis(randomZ, temp, num, 1+ptclNum*(start-1)+num*3, num*4+ptclNum*(start-1), rPosLE(3)-dx(3)*HALF, rPosLS(3))
do k=1,num
  PinitZ(k+num*3) = temp(k)
enddo
CALL setAxis(randomZ, temp, num, 1+ptclNum*(start-1)+num*4, num*5+ptclNum*(start-1), rPosLE(3), rPosLS(3)+dx(3)*HALF)
do k=1,num
  PinitZ(k+num*4) = temp(k)
enddo
CALL setAxis(randomZ, temp, num, 1+ptclNum*(start-1)+num*5, num*6+ptclNum*(start-1), rPosLE(3), rPosLS(3)+dx(3)*HALF)
do k=1,num
  PinitZ(k+num*5) = temp(k)
enddo
CALL setAxis(randomZ, temp, num, 1+ptclNum*(start-1)+num*6, num*7+ptclNum*(start-1), rPosLE(3), rPosLS(3)+dx(3)*HALF)
do k=1,num
  PinitZ(k+num*6) = temp(k)
enddo
CALL setAxis(randomZ, temp, num, 1+ptclNum*(start-1)+num*7, ptclNum+ptclNum*(start-1), rPosLE(3), rPosLS(3)+dx(3)*HALF)
do k=1,num
  PinitZ(k+num*7) = temp(k)
enddo
             endif

           
             do j=1,IonSorts
                if(ptclType(j)<=0)cycle
                do i=1,ptclNum
                   R=0.d0
                   R(1)=PinitX(i)
                   R(2)=PinitY(i)
                   R(3)=PinitZ(i)
                   R(4)=ptclVelc(1)
                   R(5)=ptclVelc(2)
                   R(6)=ptclVelc(3)
               
                   call injctParticle(R,ptclType(j))
      
                enddo
             enddo
          enddo
      
          !result output
!!$          do i=1,IonSorts
!!$             if(ptclType(i)<=0)cycle
!!$             print *,"[setParticleSquare]:injecting isort=",ptclType(i),"num=",injctNum(ptclType(i)),"injected octNum=",octNum,"rank=",rank
!!$          enddo
   

        end subroutine setParticleSquare

        subroutine setParticleSphere(centerPos,radius,ptclDistType,distParam,ptclType,ptclpercell,ptclVelc,IR,injctNum)
          implicit none
          !define the area to inject particles.
          real(kind=8)   ,intent(in)::centerPos(3),radius,distParam(32),ptclVelc(3)
          !if the value of rArea was 0, it would be neglected. 
          !if the value of rArea was negative, particles would set in whole region. 
          integer(kind=4),intent(in)::ptclDistType(3),ptclType(IonSorts),ptclpercell,IR
          !particle distribution type = 0 : white noise
          !                           = 1 : sine random
          integer(kind=4),intent(inout),dimension(IonSorts)::injctNum
          integer(kind=4)::ptclNum,octNum
          real(kind=8),dimension(ptclpercell)::PinitX,PinitY,PinitZ
          real(kind=8),dimension(ptclpercell)::PEwork
          real(kind=8),dimension(9)::R
          real(kind=8)::distToCenter,rPosS(3),rPosE(3)
          integer(kind=4)::i,index
          type(oct),pointer::p0

         ! if(debugMode>0)print *,"Starting setParticleSphere",rank
          ptclNum=0
          octNum=0

          if(radius>dx(1)*dble(NX*NXR/2) .or. radius>dx(2)*dble(NY*NYR/2) .or. &
               radius>dx(3)*dble(NZ*NZR/2))then
             print *,"Error in setParticleSphere radius is too large : radius=",radius,"space size=",dx(1)*dble(NX*NXR/2)
             stop
          endif

          do index=MinID(1,0),MaxID(1,0)
             p0=>Mesh(index)
             distToCenter=sqrt((p0%rPos(1)-centerPos(1))**2+&
                               (p0%rPos(2)-centerPos(2))**2+&
                               (p0%rPos(3)-centerPos(3))**2)
             if(distToCenter>radius)cycle
             octNum=octNum+1
             ptclNum=ptclpercell

             !add summation of particle number
             do i=1,IonSorts
                if(ptclType(i)<=0)cycle
                injctNum(ptclType(i))=injctNum(ptclType(i))+ptclNum
                if(injctNum(ptclType(i))>Niall)then
                   print *,"[Error]setParticleSphere, sum of the ptcl num isort=",ptclType,rank
                   print *,"[Error]violates Niall. ptclNum=",injctNum(ptclType(i)),"Niall=",Niall,rank
                   stop
                endif
             enddo

             rPosS=p0%rPos-dx*HALF
             rPosE=p0%rPos+dx*HALF
   
       
             rPosS=rPosS+0.00000001d0
             rPosE=rPosE-0.00000001d0
             if(ptclDistType(1)==0)then
                CALL RPBMP1( IR+p0%MrtN, PinitX, PEwork, ptclNum, rPosS(1), rPosE(1), 5)
             elseif(ptclDistType(1)==1)then
                CALL SINRAN( IR+p0%MrtN, PinitX, PEwork, ptclNum, rPosE(1), rPosS(1),distParam(1))
             endif
             if(ptclDistType(2)==0)then
                CALL RPBMP1( IR+p0%MrtN, PinitY, PEwork, ptclNum, rPosS(2), rPosE(2), 5)
             elseif(ptclDistType(2)==1)then
                CALL SINRAN( IR+p0%MrtN, PinitY, PEwork, ptclNum, rPosE(2), rPosS(2),distParam(1))
             endif
             if(ptclDistType(3)==0)then
                CALL RPBMP1( IR+p0%MrtN, PinitZ, PEwork, ptclNum, rPosS(3), rPosE(3), 5)
             elseif(ptclDistType(3)==1)then
                CALL SINRAN( IR+p0%MrtN, PinitZ, PEwork, ptclNum, rPosE(3), rPosS(3),distParam(1))
             endif
       
             do j=1,IonSorts
                if(ptclType(j)<=0)cycle
                do i=1,ptclNum
                   R=0.d0
                   R(1)=PinitX(i)
                   R(2)=PinitY(i)
                   R(3)=PinitZ(i)
                   R(4)=ptclVelc(1)
                   R(5)=ptclVelc(2)
                   R(6)=ptclVelc(3)
              
                   call injctParticle(R,ptclType(j))
          
                enddo
             enddo
          end do

          !result output
!!$          do i=1,IonSorts
!!$             if(ptclType(i)<=0)cycle
!!$             print *,"[setParticleSphere]:injecting isort=",ptclType(i),"num=",injctNum(ptclType(i)),"injected octNum=",octNum,"rank=",rank
!!$          enddo
 

        end subroutine setParticleSphere

      end subroutine particle_injct
!
      subroutine get_MaxwellianV(RandNum,MxwllV,Vt,RSeed1,RSeed2)
        use mtmod
        use const
        use message_passing_interface
        implicit none
        integer(kind=4),intent(in)::RandNum,RSeed1,RSeed2
        real(kind=8),intent(in)::Vt
        real(kind=8),dimension(RandNum),intent(inout)::MxwllV
        real(kind=8)::al
        real(kind=8),dimension(:),allocatable::drand1,drand2
        integer(kind=4)::errD1,errD2,i

        allocate(drand1(1:RandNum),stat=errD1)
        allocate(drand2(1:RandNum),stat=errD2)   
   
        call sgrnd(RSeed1)
        do i=1,RandNum
           drand1(i)=grnd()
        end do
        call sgrnd(RSeed2)
        do i=1,RandNum
           drand2(i)=grnd()
        end do

        do i=1,RandNum
           al=Vt*dsqrt(-2.d0*dlog(drand1(i)))
           MxwllV(i)=al*dcos(drand2(i)*PI2)
        end do

        deallocate(drand2)
        deallocate(drand1)
  
      end subroutine get_MaxwellianV

!***********************************************************************
      subroutine copy_particle0(iLv,iF,iK,icon)
! +-------------------------------------------------------------------+
! |     if iFLG(iK)==iF ;  copy particles to the child octs           |
! |     if icon==0 ; particle is deleted after copy process (move)    |
! |     if icon==1 ; particle is not deleted after copy process (copy)|
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use particle_set
        use param
        use init_mesh_size
        use const
        implicit none
        integer(kind=4)    :: iLv,iF,iK,icon,k
        integer(kind=4)    :: index,indexP,indexR,indexQ
        integer(kind=4)    :: ipx,ipy,ipz,ich,ILv2
        integer(kind=4)    :: octP1,octP2,octP3,octP4
        integer(kind=4)    :: octP5,octP6,octP7,octP8
        integer(kind=4)    :: octP01,octP02,octP03,octP04
        integer(kind=4)    :: octP05,octP06,octP07,octP08
        integer(kind=4)    :: Isort,Ioct
        real(kind=8)       :: R(9),pos(3)
! -- pointers --
        type(  oct), pointer :: p0
        type(prtcl), pointer :: pptc
        type(prtcl), pointer :: PrtList1,PrtList2,PrtList3,PrtList4
        type(prtcl), pointer :: PrtList5,PrtList6,PrtList7,PrtList8

        icon=0 !unused
        
        if(iLv.eq.LvMax) return
        if(maxID(1,iLv+1).eq.minID(1,iLv+1)) return 
        if(maxID(1,iLv  ).eq.minID(1,iLv  )) return 

        iLv2 = iLv+1
        indexP = maxIP(iLv2)
        indexQ = 0
        do index=minID(1,iLv),maxID(1,iLv)
           p0 => Mesh(index)
           if(p0%iFLG(iK)==iF) then 
!
              PrtList1 => p0%octCh1%ptcl%prtnxt ; octP1 = p0%octCh1%octP
              PrtList2 => p0%octCh2%ptcl%prtnxt ; octP2 = p0%octCh2%octP
              PrtList3 => p0%octCh3%ptcl%prtnxt ; octP3 = p0%octCh3%octP
              PrtList4 => p0%octCh4%ptcl%prtnxt ; octP4 = p0%octCh4%octP
              PrtList5 => p0%octCh5%ptcl%prtnxt ; octP5 = p0%octCh5%octP
              PrtList6 => p0%octCh6%ptcl%prtnxt ; octP6 = p0%octCh6%octP
              PrtList7 => p0%octCh7%ptcl%prtnxt ; octP7 = p0%octCh7%octP
              PrtList8 => p0%octCh8%ptcl%prtnxt ; octP8 = p0%octCh8%octP
!
              octP01 = octP1 ; octP02 = octP2 ; octP03 = octP3 ; octP04 = octP4
              octP05 = octP5 ; octP06 = octP6 ; octP07 = octP7 ; octP08 = octP8
!
              pos(1) = p0%rPOS(1)
              pos(2) = p0%rPOS(2)
              pos(3) = p0%rPOS(3)
              indexR = p0%octP
!
              pptc => p0%ptcl%prtnxt
!
              do k=1,indexR
!
                 if(pptc%isort.le.0) cycle
!
                 R(1) = pptc%R(1) ; R(2) = pptc%R(2) ; R(3) = pptc%R(3)
                 R(4) =(pptc%R(4)+pptc%R(7))*0.5d0
                 R(5) =(pptc%R(5)+pptc%R(8))*0.5d0
                 R(6) =(pptc%R(6)+pptc%R(9))*0.5d0
                 R(7) = pptc%R(7) ; R(8) = pptc%R(8) ; R(9) = pptc%R(9)
                 Isort= pptc%Isort
                 ipx=(idnint(sign(1.d0,R(1)-pos(1)))+1)/2
                 ipy=(idnint(sign(2.d0,R(2)-pos(2)))+2)/2
                 ipz=(idnint(sign(4.d0,R(3)-pos(3)))+4)/2
                 ich=ipx+ipy+ipz+1
!
                 indexP = indexP + 1
                 indexQ = indexQ + 1
                 if(ich.le.4) then 
                    if(ich == 1) then
                       Ioct = p0%octCh1%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList1,0)
                       PrtList1 => Pesh(indexP,iLv2)
                       octP1 = octP1 + 1
                    elseif(ich == 2) then
                       Ioct = p0%octCh2%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList2,0)
                       PrtList2 => Pesh(indexP,iLv2)
                       octP2 = octP2 + 1
                    elseif(ich == 3) then
                       Ioct = p0%octCh3%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList3,0)
                       PrtList3 => Pesh(indexP,iLv2)
                       octP3 = octP3 + 1
                    else
                       Ioct = p0%octCh4%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList4,0)
                       PrtList4 => Pesh(indexP,iLv2)
                       octP4 = octP4 + 1
                    endif
                 else
                    if(ich == 5) then
                       Ioct = p0%octCh5%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList5,0)
                       PrtList5 => Pesh(indexP,iLv2)
                       octP5 = octP5 + 1
                    elseif(ich == 6) then
                       Ioct = p0%octCh6%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList6,0)
                       PrtList6 => Pesh(indexP,iLv2)
                       octP6 = octP6 + 1
                    elseif(ich == 7) then
                       Ioct = p0%octCh7%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList7,0)
                       PrtList7 => Pesh(indexP,iLv2)
                       octP7 = octP7 + 1
                    else
                       Ioct = p0%octCh8%octN
                       Pesh(indexP,iLv2) = prtcl(R,Isort,Ioct,indexP,PrtList8,0)
                       PrtList8 => Pesh(indexP,iLv2)
                       octP8 = octP8 + 1
                    endif
                 endif
                 pptc%isort = 0
                 pptc => pptc%prtnxt
              enddo
              p0%octCh1%ptcl%PrtNxt => PrtList1 ; p0%octCh1%octP = octP1
              p0%octCh2%ptcl%PrtNxt => PrtList2 ; p0%octCh2%octP = octP2
              p0%octCh3%ptcl%PrtNxt => PrtList3 ; p0%octCh3%octP = octP3
              p0%octCh4%ptcl%PrtNxt => PrtList4 ; p0%octCh4%octP = octP4
              p0%octCh5%ptcl%PrtNxt => PrtList5 ; p0%octCh5%octP = octP5
              p0%octCh6%ptcl%PrtNxt => PrtList6 ; p0%octCh6%octP = octP6
              p0%octCh7%ptcl%PrtNxt => PrtList7 ; p0%octCh7%octP = octP7
              p0%octCh8%ptcl%PrtNxt => PrtList8 ; p0%octCh8%octP = octP8
!
              if((octP1.gt.octP01).and.(p0%octCh1%iFLG(1).ge.4)) p0%octCh1%iFLG(3)=8 
              if((octP2.gt.octP02).and.(p0%octCh2%iFLG(1).ge.4)) p0%octCh2%iFLG(3)=8 
              if((octP3.gt.octP03).and.(p0%octCh3%iFLG(1).ge.4)) p0%octCh3%iFLG(3)=8 
              if((octP4.gt.octP04).and.(p0%octCh4%iFLG(1).ge.4)) p0%octCh4%iFLG(3)=8 
              if((octP5.gt.octP05).and.(p0%octCh5%iFLG(1).ge.4)) p0%octCh5%iFLG(3)=8 
              if((octP6.gt.octP06).and.(p0%octCh6%iFLG(1).ge.4)) p0%octCh6%iFLG(3)=8 
              if((octP7.gt.octP07).and.(p0%octCh7%iFLG(1).ge.4)) p0%octCh7%iFLG(3)=8 
              if((octP8.gt.octP08).and.(p0%octCh8%iFLG(1).ge.4)) p0%octCh8%iFLG(3)=8 

              p0%octP = 0
           endif
        end do
!
        MaxIP(iLv2)=MaxIP(iLv2)+indexQ
        MaxIP( iLv)=MaxIP( iLv)-indexQ
!
        return
        end subroutine copy_particle0
!
!======================================================================
      subroutine RPBMP1(IR,XP,WK,NPT,XMN,XMX,ICHG)
!======================================================================
      implicit real*8(A-H,O-Z)
      dimension  XP(NPT),WK(NPT)
!
        DG = (XMX-XMN)/ DBLE(NPT)
        WK(1) = 0.
      do 250  I = 1 , NPT
        XP(I) = (DBLE(I)) * DG + XMN
  250 continue
!
      if(ICHG.eq.0) return
        CALL EXCHNG( IR , XP , WK , NPT , ICHG )
!
      return
      end subroutine RPBMP1
!
!======================================================================
      SUBROUTINE SINRAN(IR,A,WK,NPT,RMX,RMN,ar)
!======================================================================
      use const 
      IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION WK(NPT), A(NPT)
!
      if(ar.lt.0.01d0) ar=1.d0
      br = 1.d0/ar

        DG = (RMX-RMN)/DBLE(NPT-1)
        WK(1) = 0.
      do I = 2 , NPT
        GG = (DBLE(I)-0.5d0) * DG + RMN
        WK(I) = WK(I-1) + DSIN( 2.d0*PI*GG/(RMX-RMN) ) + br
      end do
!
!        DG = (RMX-RMN)/2.d0/PI
!
        DF = WK(NPT) / NPT
        II = 1
        J  = 1
      DO 251  I = 1 , NPT
          FV = (DBLE(I) - 5.D-1) * DF
!
  252   CONTINUE
        IF(FV.LT.WK(J+1)) GO TO 253
!
          J = J + 1
        IF(J.GT.NPT-1) GO TO 254
          GO TO 252
!
  253   CONTINUE
          A(II) = DG * ( DBLE(J-1)                       &
                       + (FV-WK(J))/(WK(J+1)-WK(J)) )
          II = II + 1
  251 CONTINUE
!
  254 CONTINUE
!
      DO 256  I = 1 , NPT
        A(I) = A(I) + RMN
  256 CONTINUE
!
      CALL EXCHNG( IR , A , WK , NPT , 5 )
!
      return
      end subroutine SINRAN
!
!======================================================================
      SUBROUTINE EXCHNG(IR,A,W,NPT,IMX)
!======================================================================
      IMPLICIT  REAL*8(A-H,O-Z)
      INTEGER J1S,J2S,IAD1U(19994)
      REAL QQQ(19999)
      EQUIVALENCE (QQQ(1),IAD1U)
      DIMENSION  A(NPT), W(NPT)
!
      INTEGER J1
      ISS = 0

! --------------
!!$  120 CONTINUE
! --------------
      J1 = MAX(IMX - 1,0)
      DO ISS = 0, J1
         CALL RANDN2 (IR, W, NPT)
!
         IF (NPT .GE. 5) THEN
            DO J1S = 0, NPT - 1, 19993
               J2S = MIN0(NPT - J1S,19993)
               DO I = 1, J2S
                  IAD1U(I) = IDNINT(W(J1S+I)*DBLE(NPT-1)) + 1
               END DO
               DO I = 1, J2S
                  AAA = A(J1S+I)
                  A(J1S+I) = A(IAD1U(I))
                  A(IAD1U(I)) = AAA
               END DO
            END DO
         ELSE
            DO I = 1, NPT
               IAD = IDNINT(W(I)*DBLE(NPT-1)) + 1
               AAA = A(I)
               A(I) = A(IAD)
               A(IAD) = AAA
            END DO
         ENDIF
      END DO
!
      return
      end subroutine exchng
!
!======================================================================
      SUBROUTINE RANDN(IR , W , NPT )
!======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(NPT)
!
      DO 100  I = 1 , 500
        CALL RANDUN(IR,RN1)
  100 CONTINUE
!
      DO 200  I = 1 , NPT
        CALL RANDUN( IR , RN1 )
        W(I) = RN1
  200 CONTINUE
!
      RETURN
      end subroutine RANDN
!
!======================================================================
      SUBROUTINE RANDUN(IX,X)
!======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
!
      RMOD = 2147483648.
      RLAMDA = 65539.
!
      IF(IX) 2, 4, 6
    2 IX = -IX
      GO TO 6
    4 IX = 3907
    6 W = IX
      W = RLAMDA * W
      IF(W-RMOD) 20, 10, 10
   10 IOVER = W / RMOD
      W = W - DBLE(IOVER)*RMOD
   20 IX = W
      X = W / RMOD
!
      return 
      end subroutine RANDUN
!
!======================================================================
      SUBROUTINE RANDN2(IR , W , NPT )
!======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      integer(kind=4) :: i
      integer(kind=4) :: IR,NPT,iover
      real(kind=8)    :: W(NPT)
      real(kind=8),parameter :: RMOD = 2147483648.d0
      real(kind=8),parameter :: RLAMDA = 65539.d0
!
      do I = 1 , 500
         if(IR.lt.0) IR=-IR
         if(IR.eq.0) IR=3907
         A=dble(IR)*RLAMDA
         if(A.ge.RMOD) then 
            iover=A/RMOD
            A=A-dble(iover)*RMOD
         end if
         IR=A
      end do
!
      do I = 1 , NPT
         if(IR.lt.0) IR=-IR
         if(IR.eq.0) IR=3907
         A=dble(IR)*RLAMDA
         if(A.ge.RMOD) then 
            iover=A/RMOD
            A=A-dble(iover)*RMOD
         end if
         IR=A
         W(i)=A/RMOD
      end do
!
      return 
      end subroutine RANDN2
!
!======================================================================
      SUBROUTINE QSTART(IR,A,WK,NPT,SD,AV)
!======================================================================
      IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION WK(NPT), A(NPT)
!
        GMAX = 5.D0 * SD
        DG = 2.D0 * GMAX / DBLE(NPT-1)
        WK(1) = 0.D0
      DO 250  I = 2 , NPT
        GG = ((DBLE(I) - 1.5D0) * DG - GMAX) / SD
        WK(I) = WK(I-1) + DEXP( -GG**2 / 2.D0)
  250 CONTINUE
!
        DF = WK(NPT) / NPT
        II = 1
        J  = 1
      DO 251  I = 1 , NPT
          FV = (DBLE(I) - 0.5D0) * DF
!
  252   CONTINUE
        IF(FV.LT.WK(J+1)) GO TO 253
!
          J = J + 1
        IF(J.GT.NPT-1) GO TO 254
          GO TO 252
!
  253   CONTINUE
!
          A(II) = DG * ( DBLE(J-1)            &
                  + (FV-WK(J))/(WK(J+1)-WK(J)) ) - GMAX
          II = II + 1
  251 CONTINUE
!
!---------------------*
  254 CONTINUE
!---------------------*
!
        SGM = 0.D+0
      DO 255  I = 1 , NPT
        SGM = SGM + A(I)
  255 CONTINUE
!
        SGM = SGM / DBLE(NPT)
      DO 256  I = 1 , NPT
        A(I) = A(I) - SGM + AV
  256 CONTINUE
!
        CALL EXCHNG( IR , A , WK , NPT , 5 )
!
      RETURN
      END subroutine QSTART
!

!***************************************************************
! FFT OF SINGLE REAL FUNCTION
!   INPUT DATA HAVE 2N ELEMENTS
!     ISIGN = 1 FOR FOURIER TRANSFORM
!          TRANSFORMD DATA SHOULD BE MULTIPLIED BY 1/N
!          WITH SEQUANCE OF
!          C(0),C(N),C(1),S(1),C(2),S(2),......C(N-1),S(N-1)
!       WHILE X(J) : J=1,2,....,2N IS EXPRESSED AS
!          X(J)= 0.5*C(0) + SUM(C(K)COS(..)+S(K)SIN(..)) + 0.5*C(N)
!
!     ISIGN = -1 FOR INVERSE FOURIER TRANSFORM
!     
!     REFERENCE: NUMERICAL RECIPES BY W.H. PRESS ET AL., CAMBRIDGE 1986
!        MODIFIED BY Y. OMURA, SEPTEMBER, 1989
!****************************************************************
      SUBROUTINE REALFT(DATA,N2,ISIGN)
        implicit none
!        implicit double precision(a-h,o-z)
        integer(kind=4)::ISIGN, N, N2, I1, N2P3
        integer(kind=4)::H1R, H1I, H2R, H2I
        REAL(kind=8) ::  WR,WI,WPR,WPI,WTEMP,THETA
        real(kind=8) ::  DATA(N2), C1,C2, WRS, WIS
        N=N2/2
        THETA=3.141592653589793D0/DBLE(N)
        C1=0.5
        if(ISIGN.eq.1) then
           C2=-0.5
           call FOUR1(DATA,N2,1)
        else
           C2=0.5
           THETA=-THETA
        end if
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0+WPR
        WI=WPI
        N2P3=2*N+3
        do I1=3,N-1,2
           WRS=SNGL(WR)
           WIS=SNGL(WI)
           H1R=C1*(DATA(I1)+DATA(N2P3-I1-1))
           H1I=C1*(DATA(I1+1)-DATA(N2P3-I1))
           H2R=-C2*(DATA(I1+1)+DATA(N2P3-I1))
           H2I=C2*(DATA(I1)-DATA(N2P3-I1-1))
           DATA(I1)=H1R+WRS*H2R-WIS*H2I
           DATA(I1+1)=H1I+WRS*H2I+WIS*H2R
           DATA(N2P3-I1-1)=H1R-WRS*H2R+WIS*H2I
           DATA(N2P3-I1)=-H1I+WRS*H2I+WIS*H2R
           WTEMP=WR
           WR=WR*WPR-WI*WPI+WR
           WI=WI*WPR+WTEMP*WPI+WI
        end do
      IF(ISIGN.EQ.1) THEN
        H1R=DATA(1)
        DATA(1)=H1R+DATA(2)
        DATA(2)=H1R-DATA(2)
      ELSE
        H1R=DATA(1)
        DATA(1)=C1*(H1R+DATA(2))
        DATA(2)=C1*(H1R-DATA(2))
        CALL FOUR1(DATA,N2,-1)
      ENDIF
      RETURN
      end subroutine REALFT

!***************************************************************
      SUBROUTINE FOUR1(DATA, N, ISIGN)
!***************************************************************
        implicit none
!      implicit double precision(a-h,o-z)
      integer(kind=4) :: I,J,N,ISTEP, ISIGN
      integer(kind=4) :: M, MMAX
      REAL(kind=8) :: WR,WI,WPR,WPI,WTEMP,THETA
      real(kind=8) :: DATA(N)
      real(kind=8) ::TEMPR, TEMPI
      J=1
      DO I=1,N,2
        IF(J.GT.I) THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
        do while((M.GE.2).AND.(J.GT.M))
          J=J-M
          M=M/2
        end do
        J=J+M
      end do
      MMAX=2
      do while(N.GT.MMAX)
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO M=1,MMAX,2
          DO I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
          end do
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
        end do
        MMAX=ISTEP
      end do
!
      RETURN
      end subroutine FOUR1
!
!****************************************************************
      subroutine wrtout(octmax,ipmax1,ipmax2,ipmax3)
!****************************************************************
        use param
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        implicit none 
        integer(kind=4) :: octmax,ipmax1,ipmax2,ipmax3
        integer(kind=4) :: i,iLv
        integer(kind=4) :: octLvw(  octmax),octPw( octmax),Csortw( octmax)
        integer(kind=4) :: iPOSw( 3,octmax),iFLGw( octmax)
        real(kind=8)    :: rPOSw( 3,octmax)
        integer(kind=4) :: octNbw(6,octmax),octChw(8,octmax),octPrtw(octmax)
        real(kind=8)    :: Fw(12,octmax)
        integer(kind=4) :: octptcl(octmax)
        integer(kind=4) :: Ioctw1( ipmax1),Ioctw2( ipmax2),Ioctw3( ipmax3)
        integer(kind=4) :: nxtw1(  ipmax1),nxtw2(  ipmax2),nxtw3(  ipmax3)
        integer(kind=4) :: Isortw1(ipmax1),Isortw2(ipmax2),Isortw3(ipmax3)
        real(kind=8)    ::   Rw1(9,ipmax1),  Rw2(9,ipmax2),  Rw3(9,ipmax3)
        character(LEN=32) :: fo_name
!-- pointers
        type(  oct), pointer :: p0
        type(prtcl), pointer :: pp

!  ---  output data for oct
        do i=1,octmax
           p0 => Mesh(i)
           octLvw(  i)=p0%octLv
            octPw(  i)=p0%octP
           Csortw(  i)=p0%Csort
            iPOSw(1,i)=p0%iPOS(1)   
            iPOSw(2,i)=p0%iPOS(2) 
            iPOSw(3,i)=p0%iPOS(3)
            rPOSw(1,i)=p0%rPOS(1)   
            rPOSw(2,i)=p0%rPOS(2) 
            rPOSw(3,i)=p0%rPOS(3)
            iFLGw(  i)=p0%iFLG(1)
!
           octNbw(1,i)=p0%octNb1%octN  ;octNbw(2,i)=p0%octNb2%octN  
           octNbw(3,i)=p0%octNb3%octN  ;octNbw(4,i)=p0%octNb4%octN  
           octNbw(5,i)=p0%octNb5%octN  ;octNbw(6,i)=p0%octNb6%octN
!
           if((p0%octLv.eq.-1).or.(p0%iFLG(1).gt.0)) then 
           octChw(1,i)=p0%octCh1%octN  ;octChw(2,i)=p0%octCh2%octN
           octChw(3,i)=p0%octCh3%octN  ;octChw(4,i)=p0%octCh4%octN
           octChw(5,i)=p0%octCh5%octN  ;octChw(6,i)=p0%octCh6%octN
           octChw(7,i)=p0%octCh7%octN  ;octChw(8,i)=p0%octCh8%octN
           else 
           octChw(1,i)=0               ;octChw(2,i)=0
           octChw(3,i)=0               ;octChw(4,i)=0
           octChw(5,i)=0               ;octChw(6,i)=0
           octChw(7,i)=0               ;octChw(8,i)=0
           end if
!
           if(p0%octLv.ne.-1) then 
           octPrtw( i)=p0%octPrt%octN
           else
           octPrtw( i)=0
           end if
!
           Fw(1,i)=p0%F(1) ; Fw( 2,i)=p0%F( 2) ; Fw( 3,i)=p0%F( 3) ; Fw( 4,i)=p0%F( 4) 
           Fw(5,i)=p0%F(5) ; Fw( 6,i)=p0%F( 6) ; Fw( 7,i)=p0%F( 7) ; Fw( 8,i)=p0%F( 8) 
           Fw(9,i)=p0%F(9) ; Fw(10,i)=p0%F(10) ; Fw(11,i)=p0%F(11) ; Fw(12,i)=p0%F(12) 
!
           octptcl(i)=p0%ptcl%Inum
!
        end do

!  ---  output data for particle (Lv=0)
        if(ipmax1.gt.0) then 
           do i=1,ipmax1
              pp => Pesh(i,0)
              Isortw1(i)=pp%Isort
               Ioctw1(i)=pp%Ioct
              Rw1(1,i)=pp%R(1) ; Rw1(2,i)=pp%R(2) ; Rw1(3,i)=pp%R(3)
              Rw1(4,i)=pp%R(4) ; Rw1(5,i)=pp%R(5) ; Rw1(6,i)=pp%R(6)
              Rw1(7,i)=pp%R(7) ; Rw1(8,i)=pp%R(8) ; Rw1(9,i)=pp%R(9)
              if(associated(pp%prtnxt)) then 
                 nxtw1(i)=pp%prtnxt%Inum
              else
                 nxtw1(i)=0
              end if
           end do
        end if
!  ---  output data for particle (Lv=1)
        if(ipmax2.gt.0) then 
           do i=1,ipmax2
              pp => Pesh(i,1)
              Isortw2(i)=pp%Isort
               Ioctw2(i)=pp%Ioct
              Rw2(1,i)=pp%R(1) ; Rw2(2,i)=pp%R(2) ; Rw2(3,i)=pp%R(3)
              Rw2(4,i)=pp%R(4) ; Rw2(5,i)=pp%R(5) ; Rw2(6,i)=pp%R(6)
              Rw2(7,i)=pp%R(7) ; Rw2(8,i)=pp%R(8) ; Rw2(9,i)=pp%R(9)
              if(associated(pp%prtnxt)) then 
                 nxtw2(i)=pp%prtnxt%Inum
              else
                 nxtw2(i)=0
              end if

           end do
        end if
!  ---  output data for particle (Lv=2)
        if(ipmax3.gt.0) then 
           do i=1,ipmax3
              pp => Pesh(i,2)
              Isortw3(i)=pp%Isort
               Ioctw3(i)=pp%Ioct
              Rw3(1,i)=pp%R(1) ; Rw3(2,i)=pp%R(2) ; Rw3(3,i)=pp%R(3)
              Rw3(4,i)=pp%R(4) ; Rw3(5,i)=pp%R(5) ; Rw3(6,i)=pp%R(6)
              Rw3(7,i)=pp%R(7) ; Rw3(8,i)=pp%R(8) ; Rw3(9,i)=pp%R(9)
              if(associated(pp%prtnxt)) then 
                 nxtw3(i)=pp%prtnxt%Inum
              else
                 nxtw3(i)=0
              end if
           end do
        end if
!  ---

        iLv=0
!        write(fo_name,555) 'wrtmain',iLv,istep
!        open (94, file=fo_name   , form=  'unformatted')
!          write(94) istep,octmax,ipmax1,ipmax2,ipmax3
!        close(94)

        write(fo_name,555) 'wrtoct_',iLv,istep
        open (95, file=fo_name   , form=  'unformatted')
          write(95) octLvw,octPw,Csortw,iPOSw,rPOSw,iFLGw,octNbw,octChw,octPrtw,Fw,octptcl
        close(95)

        iLv=0
        if(ipmax1.gt.0) then 
           write(fo_name,555) 'wrtprt_',iLv,istep
           open (97, file=fo_name   , form=  'unformatted')
             write(97) Isortw1,Ioctw1,nxtw1,Rw1
           close(97)
        end if
        iLv=1
        if(ipmax2.gt.0) then 
           write(fo_name,555) 'wrtprt_',iLv,istep
           open (98, file=fo_name   , form=  'unformatted')
             write(98) Isortw2,Ioctw2,nxtw2,Rw2
           close(98)
        end if
        iLv=2
        if(ipmax3.gt.0) then 
           write(fo_name,555) 'wrtprt_',iLv,istep
           open (99, file=fo_name   , form=  'unformatted')
             write(99) Isortw3,Ioctw3,nxtw3,Rw3
           close(99)
        end if
!
555     format(A,I2.2,'_',I6.6,'.datPL2') 
!
        return 
        end subroutine wrtout
!
!****************************************************************
      subroutine restart(time)
!****************************************************************
        use param
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        implicit none 
        real(kind=8)       :: time
        time=0 !temp
!!$        integer(kind=4)    :: octmax,ipmax1,ipmax2,ipmax3,iLv,isort,i
!!$        integer(kind=4)    :: octmax1,octmax2,octmax3
!!$        character(LEN=32) :: fo_name
!!$
!!$        call restart_setting
!!$        iLv=0
!!$        istep=Nstep0
!!$        write(fo_name,555) 'wrtmain',iLv,istep
!!$        open (94, file=fo_name   , form=  'unformatted')
!!$         read(94) istep,octmax,octmax1,octmax2,octmax3    &
!!$                               ,ipmax1, ipmax2, ipmax3
!!$        close(94)
!!$        print *, istep,octmax,octmax1,octmax2,octmax3    &
!!$                             , ipmax1, ipmax2, ipmax3
!!$!
!!$        minID(1,-1) = 1
!!$        maxID(1,-1) = Nall_ini/8
!!$        maxID(1, 0) = octmax1
!!$        maxID(1, 1) = octmax2
!!$        maxID(1, 2) = octmax3
!!$        if(maxID(1,0).eq.maxID(1,-1)) then 
!!$           minID(1,0)=maxID(1,0)
!!$        else
!!$           minID(1,0)=maxID(1,-1)+1
!!$        end if
!!$        if(maxID(1,1).eq.maxID(1,0)) then 
!!$           minID(1,1)=maxID(1,1)
!!$        else
!!$           minID(1,1)=maxID(1,0)+1
!!$        end if
!!$        if(maxID(1,2).eq.maxID(1,1)) then 
!!$           minID(1,2)=maxID(1,2)
!!$        else
!!$           minID(1,2)=maxID(1,1)+1
!!$        end if
!!$        
!!$        minID(2,-1)=maxID(1,-1) ; maxID(2,-1)=maxID(1,-1)
!!$        minID(2, 0)=maxID(1, 0) ; maxID(2, 0)=maxID(1, 0)
!!$        minID(2, 1)=maxID(1, 1) ; maxID(2, 1)=maxID(1, 1)
!!$        minID(2, 2)=maxID(1, 2) ; maxID(2, 2)=maxID(1, 2)
!!$        
!!$        maxIP(  0)=ipmax1
!!$        maxIP(  1)=ipmax2
!!$        maxIP(  2)=ipmax3
!!$        call restart_data(octmax,ipmax1,ipmax2,ipmax3)
!!$!
!!$        Isort = 0 ; Rcharge(Isort) =  0.d0 ; Rmass(Isort)   =  0.d0
!!$        isort = 1 ; Rcharge(isort) = -1.d0 ; Rmass(isort)   =  1.d0
!!$        isort = 2 ; Rcharge(isort) =  1.d0 ; Rmass(isort)   = 25.d0
!!$        
!!$        PICnumber(1) = dble(Niall)/dble(Nall_ini)
!!$        PICnumber(2) = dble(Niall)/dble(Nall_ini)
!!$        charge = (omega**2)/dble(PICnumber(1))/(4.d0*PI)
!!$        do i=1,3
!!$           cdx(i)=charge*dx(i)
!!$        end do
!!$
!!$        time = dble(Nstep0)*dt
!!$
!!$555     format(A,I2.2,'_',I6.6,'.datPL2') 

        return 
        end subroutine restart
!
!****************************************************************
      subroutine restart_data(octmax,ipmax1,ipmax2,ipmax3)
!****************************************************************
        use param
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        implicit none 
        integer(kind=4) :: octN,octLv,octP
        integer(kind=4) :: iFLG(3),Csort
        integer(kind=4) :: iPOS(3)
        real(kind=8)    :: rPOS(3)
        integer(kind=4) :: MrtN
        integer(kind=4) :: prc_bndry, octType, ptcl_loops
        integer(kind=4) :: index
        integer(kind=4) :: AllocateStatus
        integer(kind=4) :: iC(iomp0)
        real(kind=8)    :: F(18),C(6*iomp0),Z(IonSorts),G(1:18),D(1:3),O(1:16)
        integer(kind=4) :: octmax,ipmax1,ipmax2,ipmax3
        integer(kind=4) :: i,iLv
        integer(kind=4) :: octLvw(  octmax),octPw( octmax),Csortw( octmax)
        integer(kind=4) :: iPOSw( 3,octmax),iFLGw( octmax)
        real(kind=8)    :: rPOSw( 3,octmax)
        integer(kind=4) :: octNbw(6,octmax),octChw(8,octmax),octPrtw(octmax)
        real(kind=8)    :: Fw(12,octmax)
        integer(kind=4) :: octptcl(octmax)
        integer(kind=4) :: Isort,Ioct
        integer(kind=4) :: Ioctw1( ipmax1),Ioctw2( ipmax2),Ioctw3( ipmax3)
        integer(kind=4) :: nxtw1(  ipmax1),nxtw2(  ipmax2),nxtw3(  ipmax3)
        integer(kind=4) :: Isortw1(ipmax1),Isortw2(ipmax2),Isortw3(ipmax3)
        real(kind=8)    ::   Rw1(9,ipmax1),  Rw2(9,ipmax2),  Rw3(9,ipmax3)
        real(kind=8)    :: R(9)
        character(LEN=32) :: fo_name
!-- pointers
        type(  oct), pointer :: newp
        type(prtcl), pointer :: prtlist
!
!-----------------------------------------------------
        allocate(newp, stat=AllocateStatus)
        nullify(newp%octPrt)
        nullify(newp%octNb1) ; nullify(newp%octNb2)
        nullify(newp%octNb3) ; nullify(newp%octNb4)
        nullify(newp%octNb5) ; nullify(newp%octNb6)
        nullify(newp%octCh1) ; nullify(newp%octCh2)
        nullify(newp%octCh3) ; nullify(newp%octCh4)
        nullify(newp%octCh5) ; nullify(newp%octCh6)
        nullify(newp%octCh7) ; nullify(newp%octCh8)
        nullify(newp%Psort)
        nullify(newp%ptcl)   ; nullify(newp%ptclA) 
!-----------------------------------------------------

        octP  = 0
        Csort = 0

        C = 0.d0
        Z = 0.d0
        G = 0.d0
        D=0.d0
        O=0.d0
        iC= 0
        prc_bndry=0
        octType = 0
        ptcl_loops=0
        iLv=0
        write(fo_name,555) 'wrtoct_',iLv,istep
        open (95, file=fo_name   , form=  'unformatted')
          read(95) octLvw,octPw,Csortw,iPOSw,rPOSw,iFLGw,octNbw,octChw,octPrtw,Fw,octptcl
        close(95)

        do i = 1, octmax
          octN   = i
          octLv  = octLvw(i)
          octP   = octPw(i)
          Csort  = Csortw(i)
          index  = i
          MrtN   = i*8
          iFLG(1)= iFLGw(i)
          iFLG(2)= 0
          iFLG(3)= 0
!
          iPOS(1) = iPOSw(1,i)
          iPOS(2) = iPOSw(2,i)
          iPOS(3) = iPOSw(3,i)
          rPOS(1) = rPOSw(1,i)
          rPOS(2) = rPOSw(2,i)
          rPOS(3) = rPOSw(3,i)
!
          F(1)=Fw(1,i) ; F( 2)=Fw( 2,i) ; F( 3)=Fw( 3,i) ; F( 4)=Fw( 4,i)
          F(5)=Fw(5,i) ; F( 6)=Fw( 6,i) ; F( 7)=Fw( 7,i) ; F( 8)=Fw( 8,i)
          F(9)=Fw(9,i) ; F(10)=Fw(10,i) ; F(11)=Fw(11,i) ; F(12)=Fw(12,i)
!
          if(octLv.ne.-1) then 
             newp%octPrt => Mesh(octPrtw(i))
          else
             nullify(newp%octPrt)
          end if
          newp%octNb1 => Mesh(octNbw(1,i))
          newp%octNb2 => Mesh(octNbw(2,i))
          newp%octNb3 => Mesh(octNbw(3,i))
          newp%octNb4 => Mesh(octNbw(4,i))
          newp%octNb5 => Mesh(octNbw(5,i))
          newp%octNb6 => Mesh(octNbw(6,i))
!
          if((iFLG(1).gt.0).or.(octLv.eq.-1)) then  
             newp%octCh1 => Mesh(octChw(1,i))
             newp%octCh2 => Mesh(octChw(2,i))
             newp%octCh3 => Mesh(octChw(3,i))
             newp%octCh4 => Mesh(octChw(4,i))
             newp%octCh5 => Mesh(octChw(5,i))
             newp%octCh6 => Mesh(octChw(6,i))
             newp%octCh7 => Mesh(octChw(7,i))
             newp%octCh8 => Mesh(octChw(8,i))
          else
             nullify(newp%octCh1)
             nullify(newp%octCh2)
             nullify(newp%octCh3)
             nullify(newp%octCh4)
             nullify(newp%octCh5)
             nullify(newp%octCh6)
             nullify(newp%octCh7)
             nullify(newp%octCh8)
          end if

          if(octLv.ne.-1) then 
             newp%ptcl => Pesh(octPtcl(i),octLv)
          else
             nullify(newp%ptcl)
          end if
          
! -- Mesh Array --
          Mesh(index)=                               &
                 oct(octN,octLv,octP,Csort,iFLG,     &
                 iPOS,rPOS,                          &
                 MrtN, iC,                       &
                prc_bndry,                 & 
                octType,                       &
                F,C,Z,G,D,O,ptcl_loops,       & 
                newp%octPrt,                        &
                 newp%octNb1, newp%octNb2,           &
                 newp%octNb3, newp%octNb4,           &
                 newp%octNb5, newp%octNb6,           &
                 newp%octCh1, newp%octCh2,           &
                 newp%octCh3, newp%octCh4,           &
                 newp%octCh5, newp%octCh6,           &
                 newp%octCh7, newp%octCh8,           &
                 newp%Psort ,                        &
                 newp%ptcl  , newp%ptclA)
        enddo
!
        if(ipmax1.gt.0) then 

           iLv=0
           write(fo_name,555) 'wrtprt_',iLv,istep
           open (97, file=fo_name   , form=  'unformatted')
            read(97) Isortw1,Ioctw1,nxtw1,Rw1
           close(97)
           
           do i=1,ipmax1
              R(1)=Rw1(1,i) ; R(2)=Rw1(2,i) ; R(3)=Rw1(3,i)
              R(4)=Rw1(4,i) ; R(5)=Rw1(5,i) ; R(6)=Rw1(6,i)
              R(7)=Rw1(7,i) ; R(8)=Rw1(8,i) ; R(9)=Rw1(9,i)
              Isort=Isortw1(i)
              Ioct = Ioctw1(i)
              if(nxtw1(i).ne.0) then  
                 PrtList => Pesh(nxtw1(i),0)
              else
                 nullify(PrtList)
              end if
              Pesh(i,0) = prtcl(R,Isort,Ioct,i  ,PrtList,0)
           end do
              
        end if
!
        if(ipmax2.gt.0) then 

           iLv=1
           write(fo_name,555) 'wrtprt_',iLv,istep
           open (97, file=fo_name   , form=  'unformatted')
            read(97) Isortw2,Ioctw2,nxtw2,Rw2
           close(97)
           
           do i=1,ipmax2
              R(1)=Rw2(1,i) ; R(2)=Rw2(2,i) ; R(3)=Rw2(3,i)
              R(4)=Rw2(4,i) ; R(5)=Rw2(5,i) ; R(6)=Rw2(6,i)
              R(7)=Rw2(7,i) ; R(8)=Rw2(8,i) ; R(9)=Rw2(9,i)
              Isort=Isortw2(i)
              Ioct = Ioctw2(i)
              if(nxtw2(i).ne.0) then  
                 PrtList => Pesh(nxtw2(i),1)
              else
                 nullify(PrtList)
              end if
              Pesh(i,1) = prtcl(R,Isort,Ioct,i  ,PrtList,0)
           end do
              
        end if
!
        if(ipmax3.gt.0) then 

           iLv=2
           write(fo_name,555) 'wrtprt_',iLv,istep
           open (97, file=fo_name   , form=  'unformatted')
            read(97) Isortw3,Ioctw3,nxtw3,Rw3
           close(97)
           
           do i=1,ipmax3
              R(1)=Rw3(1,i) ; R(2)=Rw3(2,i) ; R(3)=Rw3(3,i)
              R(4)=Rw3(4,i) ; R(5)=Rw3(5,i) ; R(6)=Rw3(6,i)
              R(7)=Rw3(7,i) ; R(8)=Rw3(8,i) ; R(9)=Rw3(9,i)
              Isort=Isortw3(i)
              Ioct = Ioctw3(i)
              if(nxtw3(i).ne.0) then  
                 PrtList => Pesh(nxtw3(i),2)
              else
                 nullify(PrtList)
              end if
              Pesh(i,2) = prtcl(R,Isort,Ioct,i  ,PrtList,0)
           end do
              
        end if
!
555     format(A,I2.2,'_',I6.6,'.datPL2') 
!
        return 
        end subroutine restart_data

subroutine i2pos(index, x, ipos, iLv, global)
use param
use const
use init_mesh_size
use message_passing_interface

implicit none
integer (kind=4):: index, index2, iLv
integer (kind=4):: ix, iy, iz, Nint
integer (kind = 4) ::n1, n2, n3
real (kind = 8) :: factor
integer (kind = 4) :: iRx, iRy, iRz
real (kind = 8):: xstart, ystart, zstart, x(3), dR(3)
integer(kind=4) :: iPOS(3), global !global flag returns global ipos
!-------------------------- This routine converts from local index to x_i-------------------
! Note that x_i is the position of i'th cell measured at the center of the cell cube
! Rank is another input, together with index
!
!iRz = int(rank /(NXR*NYR))
!Needs to add global ipos function

index2=0

!index must be larger than Nall_ini/8
if (iLv ==0)  index2 = index-Nall_ini/8
if( iLv == -1) index2 = index

iRz = int(rank / (NXR* NYR))
n3 =  rank - iRz *NXR * NYR 
iRy = int(n3 / NXR) 
iRx = n3 - iRy *  NXR

dR(1) = (NXB  -2 *buf_width_J) * dx(1)
dR(2) = (NYB  -2 *buf_width_J) * dx(2)
dR(3) = (NZB  -2 *buf_width_J) * dx(3)

! Having the expansion to multiple levels in mind
factor = 2d0**(-iLv)
Nint = 2**LvMax*factor
n1 = NXB / factor
n2 = NYB / factor

iz = int(( index2 - 1 )/(n1*n2)) + 1
n3 = index2 - n1 * n2 *(iz-1)
iy = int(( n3 - 1) / n1)+1
ix = n3 - n1 * (iy -1) 

xstart = dR(1) * iRx - (buf_width_J + HALF) * dx(1) 
ystart = dR(2) * iRy - (buf_width_J + HALF) * dx(2)
zstart = dR(3) * iRz - (buf_width_J + HALF) * dx(3)

x(1) = ix * dx(1) * factor + xstart
x(2) = iy * dx(2) * factor + ystart
x(3) = iz * dx(3) * factor + zstart
!print*,'in i2pos', ix, iy, iz

if(global ==1) then

ipos(1) = 2 * ix * factor - factor
ipos(2) = 2 * iy  * factor - factor
ipos(3) = 2 * iz * factor - factor

else 

ipos(1) = 2 * ix * Nint - Nint
ipos(2) = 2 * iy  *  Nint - Nint
ipos(3) = 2 * iz * Nint - Nint

endif

return
end subroutine i2pos

subroutine pos2i(x, index, ipos)
!returns local ipos and CPU number from position

use param
!use const
use init_mesh_size
use message_passing_interface

implicit none
integer (kind=4):: index
integer (kind=4):: ix, iy, iz
integer(kind = 4)::ipos(3)
real (kind = 8):: x(3)
!---------------------------------------------------------------------
! this is usable for Lv 0 only at this point  of time. 09/10/10
ix = int((x(1) - R_lim(0,1))/dx(1)) +1
iy = int((x(2) - R_lim(0,2))/dx(2)) +1
iz = int((x(3) - R_lim(0,3))/dx(3)) +1

!ix = int((x(1) - Min_R_lim(1))/dx(1)) +1   !GLOBAL iPOS(1)
!iy = int((x(2) - Min_R_lim(2))/dx(2)) +1   !GLOBAL iPOS(2)
!iz = int((x(3) - Min_R_lim(3))/dx(3)) +1   !GLOBAL iPOS(3)

!Mtn_temp=MnPOS(ix, iy, iz)

!!$i=1
!!$do while(Nall_arr <= Mtn_temp)
!!$   
!!$   i=i+1      
!!$   
!!$enddo

index = (iz-1)*(NXB*NYB)+(iy-1)*(NXB)+ix+Nall_ini/8

ipos(1) =  ix 
ipos(2) =  iy 
ipos(3) =  iz 

return
end subroutine pos2i

!yagi 2012/02/13 plz don't use for tuning
subroutine cnv_HtoN(iPOSH,iLv,iPOSN)
  use param
  implicit none
  integer(kind=4),dimension(3),intent(in)::iPOSH
  integer(kind=4),dimension(3),intent(out)::iPOSN
  integer(kind=4),intent(in)::iLv

  iPOSN=(iPOSH-sConst(iLv))/intNxt(iLv)

end subroutine cnv_HtoN

subroutine cnv_NtoH(iPOSN,iLv,iPOSH)
  use param
  implicit none
  integer(kind=4),dimension(3),intent(out)::iPOSH
  integer(kind=4),dimension(3),intent(in)::iPOSN
  integer(kind=4),intent(in)::iLv

  iPOSH=iPOSN*intNxt(iLv)+sConst(iLv)

end subroutine cnv_NtoH

!---------------------------------------------------------------------

!-----------------------Morton related routines------------------------
!***********************************************************************
!2012/01/10 yagi
      subroutine get_indexNumber(iPOS,index)
        use param
        implicit none
        integer(kind=4),intent(in),dimension(3)::iPOS
        integer(kind=4)                        ::index 

        if(MrtNisUsable/=0)then
           call Morton_number(iPOS,index)
        else
           call Raster_number(iPOS,index)
        endif

      end subroutine get_indexNumber

!2012/02/08 yagi
      subroutine get_indexNumberN(iPOS,index)
        use param
        implicit none
        integer(kind=4),intent(in),dimension(3)::iPOS
        integer(kind=4)                        ::index

        if(MrtNisUsable/=0)then
           call Morton_numberN(iPOS,index)
        else
           call Raster_numberN(iPOS,index)
        endif

      end subroutine get_indexNumberN
!***********************************************************************

!***********************************************************************
subroutine Morton_number(iPOS,Mn)
! +-------------------------------------------------------------------+
! |                                                                   |
! | Get Morton index from each cell  and 
! | create Besh in that order
! |                                                                   |
! |  Morton ordering:                                                 |
! |   From bits of each cell index (i,j,k), ....                      |
! |                                                                   |
! |   (Example)                                                       |
! |     cell index = (3,4) in two-dimensional 2^3 x 2^3 meshes        |
! |     (3,4) => (011,100) => L = 100101 = 37.                        |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************

        use message_passing_interface
        use init_mesh_size
        use const
        use param
! -------------------------------------------
        implicit none
        integer(kind=4),parameter :: N2Max=10
        integer(kind=4), intent(in) :: iPOS(3)
        integer(kind=4):: iPOS2(3)
        integer(kind=4) :: i,ixx,iyy,izz
        integer(kind=4) :: n2x(1:N2Max),n2y(1:N2Max),n2z(1:N2Max)

!---Variables necessary for BMesh
        integer(kind=4)    :: Mn

        Mn=0
        n2x=0
        n2y=0
        n2z=0

!    Convert from Hieralchial to Normal ordering
        iPOS2=(iPOS+LvMax2)/(2*LvMax2)

        ixx=iPOS2(1)-1
        iyy=iPOS2(2)-1
        izz=iPOS2(3)-1

!----Periodic B. C.-----
        if( ixx < 0) then
           ixx = ixx + NX_2t
        else if (ixx > NX_2t)then
           ixx = ixx - Nx_2t
        endif

        if( iyy < 0) then
           iyy = iyy + NY_2t
        else if (iyy > NY_2t )then
           iyy = iyy - NY_2t
        endif

        if( izz < 0) then
!if(rank==0) print*,'B:',ixx, iyy,izz
           izz = izz + NZ_2t
!if(rank==0) print*,'A:',ixx, iyy,izz
        else if (izz > NZ_2t)then
           izz = izz - NZ_2t
        endif
!--------------------------------------

!!$         if( ixx < 0) ixx = NXR*NXB-1
!!$          if( iyy < 0) iyy = NYR*NYB-1
!!$          if( izz < 0) izz = NZR*NZB-1
!!$          
!!$          if( ixx > NXR*NXB-1) ixx = 0
!!$          if( iyy >  NYR*NYB-1) iyy= 0
!!$          if( izz > NZR*NZB-1 ) izz = 0

        do i=1,N2Max
           n2x(i) = mod(ixx ,2)
           n2y(i) = mod(iyy ,2)
           n2z(i) = mod(izz ,2)
           ixx = ixx/2
           iyy = iyy/2
           izz = izz/2
        enddo
        do i=1,N2Max
           Mn = Mn             &
                + n2x(i) * (2**(3*i-3))    &
                + n2y(i) * (2**(3*i-2))    &
                + n2z(i) * (2**(3*i-1))    
        enddo
        return
      end subroutine Morton_number

!***********************************************************************
      subroutine Morton_numberN(iPOS,Mn) !Morton_number for Normal iPos
        use message_passing_interface
        use init_mesh_size
        use const
        use param
! -------------------------------------------
        implicit none
        integer(kind=4),parameter :: N2Max=10
        integer(kind=4), intent(in) :: iPOS(3)
        integer(kind=4):: iPOS2(3)
        integer(kind=4) :: i,ixx,iyy,izz
        integer(kind=4) :: n2x(1:N2Max),n2y(1:N2Max),n2z(1:N2Max)

!---Variables necessary for BMesh
        integer(kind=4)    :: Mn

        iPOS2=iPOS

        Mn=0
        n2x=0
        n2y=0
        n2z=0

        ixx=iPOS2(1)-1
        iyy=iPOS2(2)-1
        izz=iPOS2(3)-1

!----Periodic B. C.-----
        if( ixx < 0) then
           ixx = ixx + NX_2t
        else if (ixx > NX_2t)then
           ixx = ixx - Nx_2t
        endif

        if( iyy < 0) then
           iyy = iyy + NY_2t
        else if (iyy > NY_2t )then
           iyy = iyy - NY_2t
        endif

        if( izz < 0) then
!if(rank==0) print*,'B:',ixx, iyy,izz
           izz = izz + NZ_2t
!if(rank==0) print*,'A:',ixx, iyy,izz
        else if (izz > NZ_2t)then
           izz = izz - NZ_2t
        endif
!--------------------------------------
        do i=1,N2Max
           n2x(i) = mod(ixx ,2)
           n2y(i) = mod(iyy ,2)
           n2z(i) = mod(izz ,2)
           ixx = ixx/2
           iyy = iyy/2
           izz = izz/2
        enddo
        do i=1,N2Max
           Mn = Mn             &
                + n2x(i) * (2**(3*i-3))    &
                + n2y(i) * (2**(3*i-2))    &
                + n2z(i) * (2**(3*i-1))    
        enddo
        return
      end subroutine Morton_numberN

!***********************************************************************
!2012/01/10 yagi
      subroutine Raster_number(iPOS,index)
        use init_mesh_size
        use param
        use message_passing_interface
        implicit none
        integer(kind=4),intent(in)::iPOS(3)
        integer(kind=4)           ::iPOS2(3),index
        integer(kind=4)           ::rankiPos(3),iPosP(3)
        
        iPOS2=(iPOS+LvMax2)/(2*LvMax2)

        rankiPos(1)=iPos2(1)/NXB_2
        rankiPos(2)=iPos2(1)/NXB_2
        rankiPos(3)=iPos2(1)/NXB_2
        iPosP(1) = iPos2(1) - rankiPos(1)*NXB_2
        iPosP(2) = iPos2(2) - rankiPos(2)*NYB_2
        iPosP(3) = iPos2(3) - rankiPos(3)*NZB_2
        index = (iPosP(1)-1) + (iPosP(2)-1)*NXB_2 + (iPosP(3)-1)*NXB_2*NYB_2 + rank*NXB_2*NYB_2*NZB_2
        print*,"Raster_number rank=",rank,"iPos2=",iPos2,"index=",index

      end subroutine Raster_number
!***********************************************************************
!***********************************************************************
      subroutine Raster_numberN(iPOS,index)
        use init_mesh_size
        use param
        use message_passing_interface
        implicit none
        integer(kind=4),intent(in)::iPOS(3)
        integer(kind=4)           ::iPOS2(3),index,myrank
        integer(kind=4)           ::rankiPos(3),iPosP(3)

        iPOS2=iPOS

        rankiPos(:)=0
        iPosP(:)=0
             
        rankiPos(1)=(iPos2(1)-1)/NXB_2
        rankiPos(2)=(iPos2(2)-1)/NYB_2
        rankiPos(3)=(iPos2(3)-1)/NZB_2
        iPosP(1) = iPos2(1) - rankiPos(1)*NXB_2
        iPosP(2) = iPos2(2) - rankiPos(2)*NYB_2
        iPosP(3) = iPos2(3) - rankiPos(3)*NZB_2
        
        myrank = rankiPos(1) + rankiPos(2)*NXR + rankiPos(3)*NXR*NYR
        
        index = (iPosP(1)-1) + (iPosP(2)-1)*NXB_2 + (iPosP(3)-1)*NXB_2*NYB_2 + myrank*NXB_2*NYB_2*NZB_2

!        print*,"Raster_numberN rank=",rank," myrank=",myrank,"| iPos2=",iPos2,"| index=",index
      end subroutine Raster_numberN

!***********************************************************************
!2012/01/10 yagi
      subroutine get_inverseIndex(index, iPOS, iPOSNFlg)
        use param
        implicit none
        integer(kind=4),intent(in)  ::index,iPOSNFlg
        integer(kind=4),dimension(3)::iPOS
        
        if(MrtNisUsable/=0)then
           if(iPOSNFlg/=0)then
              call Inverse_MortonN(index,iPOS)
           else
              call Inverse_Morton(index,iPOS)
           endif
        else
           if(iPOSNFlg/=0)then
              call Inverse_RasterN(index,iPOS)
           else
              call Inverse_Raster(index,iPOS)
           endif
        endif

      end subroutine get_InverseIndex
!***********************************************************************

!***********************************************************************
    subroutine Inverse_Morton(Mn, iPOS)
      use param
        implicit none
        integer(kind=4),parameter :: N2Max=10
        integer(kind=4) :: a, b, c
        integer(kind=4), intent(in)::Mn
        integer(kind=4):: Mn2
        integer(kind=4), intent(out):: iPOS(3)
        integer(kind=4) :: i,ixx,iyy,izz
        integer(kind=4) :: n2x(1:N2Max),n2y(1:N2Max),n2z(1:N2Max)

        n2x=0
        n2y=0
        n2z=0
        ixx=0
        iyy=0
        izz=0        
        
        Mn2=Mn
        
        do i=N2Max, 1, -1

           a = 2**(3*i-1)
           b = 2**(3*i-2)
           c = 2**(3*i-3)

           if(a==0 .or. b==0 .or. c==0) then
              print*,'divided by zero. Decrease N2Max'
              stop
           endif
           n2z(i) = Mn2 / a

           Mn2 = Mn2 - n2z(i) * a
           n2y(i) = Mn2 / b
           Mn2 = Mn2 - n2y(i) * b
           n2x(i) = Mn2 / c
           Mn2 = Mn2 - n2x(i) * c

        enddo

        do i=N2Max, 1, -1
           a = 2**(i-1)
           ixx = ixx + n2x(i) * a
           iyy = iyy + n2y(i) * a
           izz = izz + n2z(i) * a
        enddo
        iPOS(1)=ixx+1
        iPOS(2)=iyy+1
        iPOS(3)=izz+1
 !Convert iPOS from Normal to Hieralchial ordering
        iPOS = 2*LvMax2*iPOS-LvMax2
        return
      end subroutine Inverse_Morton 

!***********************************************************************

!***********************************************************************
    subroutine Inverse_MortonN(Mn, iPOS)
      use param
        implicit none
        integer(kind=4),parameter :: N2Max=10
        integer(kind=4) :: a, b, c
        integer(kind=4), intent(in)::Mn
        integer(kind=4):: Mn2
        integer(kind=4), intent(out):: iPOS(3)
        integer(kind=4) :: i,ixx,iyy,izz
        integer(kind=4) :: n2x(1:N2Max),n2y(1:N2Max),n2z(1:N2Max)

        n2x=0
        n2y=0
        n2z=0
        ixx=0
        iyy=0
        izz=0        
        
        Mn2=Mn
        
        do i=N2Max, 1, -1

           a = 2**(3*i-1)
           b = 2**(3*i-2)
           c = 2**(3*i-3)

           if(a==0 .or. b==0 .or. c==0) then
              print*,'divided by zero. Decrease N2Max'
              stop
           endif
           n2z(i) = Mn2 / a

           Mn2 = Mn2 - n2z(i) * a
           n2y(i) = Mn2 / b
           Mn2 = Mn2 - n2y(i) * b
           n2x(i) = Mn2 / c
           Mn2 = Mn2 - n2x(i) * c

        enddo

        do i=N2Max, 1, -1
           a = 2**(i-1)
           ixx = ixx + n2x(i) * a
           iyy = iyy + n2y(i) * a
           izz = izz + n2z(i) * a
        enddo
        iPOS(1)=ixx+1
        iPOS(2)=iyy+1
        iPOS(3)=izz+1

        return
      end subroutine Inverse_MortonN 
!***********************************************************************
!***********************************************************************
      subroutine Inverse_Raster(index,iPOS)
        use param
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4),intent(in)::index
        integer(kind=4),dimension(3)::iPOS,iPosP,rankiPos
        integer(kind=4)::ix,iy,iz,myrank

        myrank = int(index/(NXB_2*NYB_2*NZB_2))
        
        ix=myrank
        iy=myrank/NXR
        iz=myrank/(NXR*NYR)
        rankiPos(1)=mod(ix,NXR)
        rankiPos(2)=mod(iy,NYR)
        rankiPos(3)=mod(iz,NZR)
        
        iPosP(1) = mod(index,NXB_2)+1
        iPosP(2) = mod(Int(index/NXB_2),NYB_2)+1
        iPosP(3) = mod(Int(index/(NXB_2*NYB_2)),NZB_2)+1
        
        iPos(1) = iPosP(1) + rankiPos(1)*NXB_2
        iPos(2) = iPosP(2) + rankiPos(2)*NYB_2
        iPos(3) = iPosP(3) + rankiPos(3)*NZB_2

        iPOS = 2*LvMax2*iPOS-LvMax2
        return
      end subroutine Inverse_Raster
!***********************************************************************
!***********************************************************************
      subroutine Inverse_RasterN(index,iPOS)
        use param
        use init_mesh_size
        use message_passing_interface
        implicit none
        integer(kind=4),intent(in)::index
        integer(kind=4),dimension(3)::iPOS,iPosP,rankiPos
        integer(kind=4)::ix,iy,iz,myrank

        rankiPos(:)=0
        iPosP(:)=0
        
        myrank = int(index/(NXB_2*NYB_2*NZB_2))
        ix=myrank
        iy=myrank/NXR
        iz=myrank/(NXR*NYR)
        rankiPos(1)=mod(ix,NXR)
        rankiPos(2)=mod(iy,NYR)
        rankiPos(3)=mod(iz,NZR)
        
        iPosP(1) = mod(index,NXB_2)+1
        iPosP(2) = mod(Int(index/NXB_2),NYB_2)+1
        iPosP(3) = mod(Int(index/(NXB_2*NYB_2)),NZB_2)+1
        
        iPos(1) = iPosP(1) + rankiPos(1)*NXB_2
        iPos(2) = iPosP(2) + rankiPos(2)*NYB_2
        iPos(3) = iPosP(3) + rankiPos(3)*NZB_2

!        print '(3i3,A,i4,A,2i4)',iPOS,"|",index,"|",rank,MnDisp
!        print *,"rankiPos",rankiPos
        return
      end subroutine Inverse_RasterN
!***********************************************************************
!***********************************************************************
!2012/01/10 yagi
      subroutine get_neighbourIndex(iPOS,IDarr)
        use param
        implicit none
        integer(kind=4),intent(in),dimension(3)::iPOS
        integer(kind=4),dimension(27)          ::IDarr
        
        if(MrtNisUsable/=0)then
           call Morton_neighbours(iPOS,IDarr,1)
        else
           call Raster_neighbours(iPOS,IDarr)
        endif

      end subroutine get_neighbourIndex
!***********************************************************************

!***********************************************************************
      subroutine Morton_neighbours(iPOS, Mnarr, mnnb_flag)
! +-------------------------------------------------------------------+
! |                                                                   |
! | Get all the adjacent neighbours' Morton indecies from iPOS
! |                                                                   |
! |  Morton ordering:                                                 |
! |   From bits of each cell index (i,j,k), ....                      |
! |                                                                   |
! |   (Example)                                                       |
! |     cell index = (3,4) in two-dimensional 2^3 x 2^3 meshes        |
! |     (3,4) => (011,100) => L = 100101 = 37.                        |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************
        use message_passing_interface
        use init_mesh_size
        use const
        use param
! -------------------------------------------
        implicit none
        integer(kind=4),parameter :: N2Max=20
        integer(kind=4) :: i,ix,iy,iz, ix2, iy2, iz2
        integer(kind=4) :: j,k,m
        integer(kind=4),intent(in) :: iPOS(3)
        integer(kind=4) :: ixx(3),iyy(3),izz(3)
        integer(kind=4) :: n2x(1:N2Max),n2y(1:N2Max),n2z(1:N2Max)

!---Variables necessary for BMesh
        integer(kind=4)    :: Mnarr(27) !There are 26 neighbours (3^3-1(self)), actually
! Mnarr(14) is self (0, 0, 0)
! If Mnnb_flag==1, this program computes only 6 faces
        integer(kind=4) :: Mnnb_flag

!        print *, ' in Morton_number ..."'
        Mnarr(:)=0
        n2x=0
        n2y=0
        n2z=0
        m=1

!Convert iPOS system from hieralchial to Normal ordering

        ix= (iPOS(1)+LvMax2)/(2*LvMax2)
        iy= (iPOS(2)+LvMax2)/(2*LvMax2)
        iz= (iPOS(3)+LvMax2)/(2*LvMax2)

        ix=ix-1
        iy=iy-1
        iz=iz-1

        if(Mnnb_flag ==0) then        
! Remember,
! octNb1 <=> Mnarr(13)
! octNb2 <=> Mnarr(15)
! octNb3 <=> Mnarr(11)
! octNb4 <=> Mnarr(17)
! octNb5 <=> Mnarr(5)
! octNb6 <=> Mnarr(23)
!!$
!!$           ixx = (/ix -1, ix, ix+1/) !-1, 0, +1
!!$           iyy = (/iy -1, iy, iy+1/)
!!$           izz = (/iz -1, iz, iz+1/)
!!$        
!!$! For the case that index points outside of the simulation box
!!$! Points to the other end  ===> Periodic Boundary condition
!!$           if( ixx(1) < 0) ixx(1) = NXR*NXB
!!$           if( iyy(1) < 0) iyy(1) = NYR*NYB
!!$           if( izz(1) < 0) izz(1) = NZR*NZB
!!$
!!$           if( ixx(3) > NXR*NXB-1) ixx(3) = 0
!!$           if( iyy(3) <  NYR*NYB-1) iyy(3)= 0
!!$           if( izz(3) < NZR*NZB-1 ) izz(3) = 0
!!$
!!$        do j=1, 3
!!$           do k = 1, 3
!!$              do l = 1, 3
!!$                     do i=1,N2Max
!!$                        n2x(i) = mod(ixx(l) ,2)
!!$                        n2y(i) = mod(iyy(k) ,2)
!!$                        n2z(i) = mod(izz(j) ,2)
!!$                        ixx = ixx/2
!!$                        iyy = iyy/2
!!$                        izz = izz/2
!!$                     enddo
!!$                     do i=1,N2Max
!!$                        Mnarr(m) = Mnarr(m)             &
!!$                             + n2x(i) * (2**(3*i-3))    &
!!$                             + n2y(i) * (2**(3*i-2))    &
!!$                             + n2z(i) * (2**(3*i-1))    
!!$                     enddo
!!$                     m=m+1
!!$                  enddo
!!$               enddo
!!$            enddo
!!$            
 else if(Mnnb_flag == 1) then
! octNb1 <=> Mnarr(1)
! octNb2 <=> Mnarr(4)
! octNb3 <=> Mnarr(2)
! octNb4 <=> Mnarr(5)
! octNb5 <=> Mnarr(3)
! octNb6 <=> Mnarr(6)
    m=1
    Mnarr=0
!    Mnarr(1)=Nb1,2=Nb3, 3,  
    do j=-1,1, 2
       
       do k=1,3 
          ixx  = (/1, 0, 0/)
          iyy  = (/0, 1, 0/)
          izz  = (/0, 0, 1/)
          
          ix2 = ix + ixx(k)*j
          iy2 = iy + iyy(k)*j
          iz2 = iz + izz(k)*j

!!$          if( ix2 < 0) ix2 = NXR*NXB-1
!!$          if( iy2 < 0) iy2 = NYR*NYB-1
!!$          if( iz2 < 0) iz2 = NZR*NZB-1
!!$          
!!$          if( ix2 > NXR*NXB-1) ix2 = 0
!!$          if( iy2 >  NYR*NYB-1) iy2 = 0
!!$          if( iz2 > NZR*NZB-1 ) iz2 = 0
          
          if( ix2 < 0) ix2 = NXR*NXB/2-1
          if( iy2 < 0) iy2 = NYR*NYB/2-1
          if( iz2 < 0) iz2 = NZR*NZB/2-1
          
          if( ix2 > NXR*NXB/2-1) ix2 = 0
          if( iy2 > NYR*NYB/2-1) iy2 = 0
          if( iz2 > NZR*NZB/2-1) iz2 = 0

          do i=1,N2Max
             n2x(i) = mod(ix2 ,2)
             n2y(i) = mod(iy2 ,2)
             n2z(i) = mod(iz2 ,2)
             ix2 = ix2/2
             iy2 = iy2/2
             iz2 = iz2/2
          enddo
          do i=1,N2Max
             Mnarr(m) = Mnarr(m)             &
                  + n2x(i) * (2**(3*i-3))    &
                  + n2y(i) * (2**(3*i-2))    &
                  + n2z(i) * (2**(3*i-1))    
          enddo
          m=m+1
       enddo
    enddo
 endif
 
! print*,Mnarr(1:6)
 return
end subroutine Morton_neighbours
!***********************************************************************
subroutine Raster_neighbours(iPOS,IDarr)
  use init_mesh_size
  use const
  use param
  !
  use message_passing_interface
  !
  implicit none
  integer(kind=4),intent(in)::iPOS(3)
  integer(kind=4),dimension(3)::iPOStemp,iPOS2!,iPOS3
  integer(kind=4),dimension(27)::IDarr!,temparr
  integer(kind=4)           ::i,ixx(6),iyy(6),izz(6)!,temp
  ! octNb1 <=> Mnarr(1)
  ! octNb2 <=> Mnarr(4)
  ! octNb3 <=> Mnarr(2)
  ! octNb4 <=> Mnarr(5)
  ! octNb5 <=> Mnarr(3)
  ! octNb6 <=> Mnarr(6)  
  ixx=(/-1, 0, 0,1,0,0/)
  iyy=(/ 0,-1, 0,0,1,0/)
  izz=(/ 0, 0,-1,0,0,1/)
  iPOS2=(iPOS+LvMax2)/(2*LvMax2)
  
  iPostemp(:)=0
  
  do i=1,6
     iPOStemp(1)=iPOS2(1)+ixx(i)
     iPOStemp(2)=iPOS2(2)+iyy(i)
     iPOStemp(3)=iPOS2(3)+izz(i)

     if(iPOStemp(1)<1)iPOStemp(1)=NXR*NXB_2
     if(iPOStemp(2)<1)iPOStemp(2)=NYR*NYB_2
     if(iPOStemp(3)<1)iPOStemp(3)=NZR*NZB_2

     if(iPOStemp(1)>NXR*NXB_2)iPOStemp(1)=1
     if(iPOStemp(2)>NYR*NYB_2)iPOStemp(2)=1
     if(iPOStemp(3)>NZR*NZB_2)iPOStemp(3)=1
     !
!     temparr(1+(i-1)*3:i*3)=iPOStemp
     !
     call Raster_numberN(iPOStemp,IDarr(i))

  enddo

!  call Raster_numberN(iPOS2,temp)
  !call Inverse_RasterN(temp,iPOS3)
  !print '(3i3,A,18i4,A,2i4)',iPOS2,"|",temparr(1:18),"|",rank,MnDisp
  !print '(3i3,A,18i4,A,2i4)',iPOS3,"|",temparr(1:18),"|",rank,MnDisp
!  print '(i9,A,6i12,A,2i4)',temp,"|",IDarr(1:6),"|",rank,MnDisp
end subroutine Raster_neighbours

!***********************************************************************
      subroutine sort_heap(n,arr)
! +-------------------------------------------------------------------+
! |                                                                   |
! | Sorting routine using heap-sort method                            |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************
        integer(kind=4) :: n
        integer(kind=4) :: arr(n)
        integer(kind=4) :: i
! -----------------------------------------------
!
        print *, ' in sort_heap ..."'
!
        do i=n/2,1,-1
          call shift_down(i,n)
        enddo
!
        do i=n,2,-1
          call swap(arr(1),arr(i))
          call shift_down(1,i-1)
        enddo
!
!  ***********************************************
        contains
          subroutine shift_down(l,r)
            integer(kind=4),intent(in) :: l,r
            integer(kind=4) :: j,jold
            integer(kind=4) :: a
            a = arr(l)
            jold = l
            j    = l+1
            do
              if(j > r) exit
              if(j < r) then
                if(arr(j) < arr(j+1)) j=j+1
              endif
              if(a >= arr(j)) exit
              arr(jold) = arr(j)
              jold = j
              j    = j+1
            enddo
            arr(jold) = a
          end subroutine shift_down
!
          subroutine swap(a,b)
            integer(kind=4),intent(inout) :: a,b
            integer(kind=4) :: dum
            dum = a
            a = b
            b = dum
          end subroutine swap
      end subroutine sort_heap


recursive subroutine quicksort(a, first, last)
  implicit none
  integer(kind=4)  a(*)
  integer(kind=4) first, last !first is the first element of a array( namely 1), last is the length of the array
  integer(kind=4) i, j
  real(kind=8)  x, t
  
  if(first>=last)return

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i - 1) call quicksort(a, first, i - 1)
  if (j + 1 < last)  call quicksort(a, j + 1, last)
end subroutine quicksort

!------------------------------Morton Related Routines END----------------------------

subroutine DDD_check
 use message_passing_interface
 use init_mesh_size
 use oct_set
 implicit none 
 integer(kind=4)::load_sum,flagsum,temp,iLv
 real(kind=8)::average_ratio,load_ratio
 
 !temptime=mpi_wtime()

 DDDflag=0
 if(DDD<0)DDDflag=1

 flagsum=0

 nloop=0

!simple check of particle loops
 do iLv=0,LvMax
    temp=2**iLv
    nloop=nloop+temp*(MaxIP(iLv)-(MaxID(1,iLv)-MinID(1,iLv)+MaxID(3,iLv)-MinID(3,iLv)))
 enddo

 !print *,"elapse time of collect_ptcl_loops=",mpi_wtime()-temptime
 
 call MPI_Allreduce(nloop,load_sum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
 
 average_ratio=1.d0/dble(nprocs)
 load_ratio=dble(nloop)/dble(load_sum)

 if(average_ratio*1.3 < load_ratio)DDDflag=1
 !if(nprocs<10)print *,"average_ratio:",average_ratio,"load_ratio:",load_ratio,rank
 !if(DDDflag==1)print *,"rank=",rank," needs DDD!! load_ratio=",load_ratio,"average_ratio=",average_ratio

 call MPI_Allreduce(DDDflag,flagsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
 
 DDDflag=flagsum
 if(DDDflag/=0)then
    DDDcount=DDDcount+1
    !if(rank==0)print *,"DDD should execute at this step DDDflag=",DDDflag
endif

!print *,"elapse time of DDDcheck=",mpi_wtime()-temptime


end subroutine DDD_check


subroutine load_balance
!This subroutine collects particle loops from all the other processes and compute the 
!average particle loops. It then re-compute a new Mn2CPU variable based on the average
!particle loops. The end product of this subroutine is Mn2CPU(nprocs) variable
!The value in the array indicates the largest Morton number in each process. 

use message_passing_interface
use init_mesh_size
use oct_set
implicit none

integer(kind=4):: equal_loop, equal_loop2, whole_ncells, whole_ptcl_loops
 !whole_ncells=Nall_ini*nprocs 
integer(kind=4):: i, j, k
integer(kind=4), allocatable, dimension(:):: Nall_arr !, rcounts2
integer(kind=4), allocatable, dimension(:):: loop, proc_total_loops, rdispls 
integer(kind=4), allocatable:: accum_loops(:)

integer(kind=4), allocatable, dimension(:):: loop_final ! for debug

type(oct), pointer::p0


!backup Mn2CPU
Mn2CPU_old=Mn2CPU

!nloop is calculated at DDD_check
call collect_ptcl_loops(nloop) ! Collect the particle loops of each level and bring them up to Lv =-1

!print *,"[load_balance]:nloop=",nloop,"istep=",istep,"rank=",rank

!if(rank==0)print*,'Mn2CPU Before',Mn2CPU

call MPI_Allreduce(nloop, whole_ptcl_loops, 1 , MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

whole_ncells=(Nall_ini/8)*nprocs ! this number is constant
Nall = MaxID(1,-1)-MinID(1,-1) + 1

allocate(Nall_arr(nprocs)) 
allocate(rdispls(nprocs))
allocate(accum_loops(whole_ncells))
allocate(loop(Nall)) !loop is the number of ptcl loops for each cell in the local process
allocate(proc_total_loops(Nall)) !total_number of accum_loops in the local process

k = 1
do j = MinID(1,-1), MaxID(1,-1)
   p0=>Mesh(j)
   loop(k)= p0% ptcl_loops
   k=k+1
end do

call MPI_allgather(Nall, 1, MPI_INTEGER, Nall_arr, 1, MPI_INTEGER, MPI_comm_world, ierr)

!Nall_arr records the number of Lv=-1 cells in each process
rdispls(1)=0
do i=2, nprocs
   rdispls(i)=rdispls(i-1)+Nall_arr(i-1)
enddo

proc_total_loops=loop
do i=2, nall
   proc_total_loops(i)=proc_total_loops(i)+proc_total_loops(i-1)
enddo

call MPI_allgatherV(proc_total_loops, nall, MPI_INTEGER, accum_loops, Nall_arr, rdispls, MPI_INTEGER, MPI_comm_world, ierr)

equal_loop=whole_ptcl_loops/nprocs
!if(rank==0)print *,"[load_balance]:equal_loop=",equal_loop

Mn2CPU=Nall_arr ! initially, Mn2CPU is the same as Nall_arr

do i=2, nprocs
   Mn2CPU(i)=Mn2CPU(i)+Mn2CPU(i-1)
   accum_loops(Mn2CPU(i-1)+1:Mn2CPU(i))= accum_loops(Mn2CPU(i-1)+1:Mn2CPU(i))+accum_loops(Mn2CPU(i-1))
enddo

!----Now Slice accum_loops by equal_loop 
Mn2CPU(:)=0
equal_loop2=equal_loop

!rcounts3(:)=0 !rcounts3 records the cell number after load balance in each process 

allocate(loop_final(nprocs))
j=1
do i=1, whole_ncells

   MN2CPU(j)=Mn2CPU(j)+1  ! add up cell number in each index which corresponds to process number
   if (accum_loops(i)>= equal_loop2) then
      j=j+1
      equal_loop2=equal_loop2+equal_loop
   endif
   if(j>nprocs)exit
enddo

if(istep/=0)then
   do i=2, nprocs-1
      if(Mn2CPU(i) < 4) Mn2CPU(i)= 4
      Mn2CPU(i) = Mn2CPU(i) + Mn2CPU(i-1) !Remember, Mn2CPU is cumulative.
      !if(Mn2CPU(i) <= Mn2CPU_old(i-1))Mn2CPU(i)=Mn2CPU_old(i-1)+4
      !if(Mn2CPU(i) >= Mn2CPU_old(i+1))Mn2CPU(i)=Mn2CPU_old(i+1)-4
   enddo
else
   do i=2, nprocs-1
      if(Mn2CPU(i) < 4) Mn2CPU(i)= 4
      Mn2CPU(i) = Mn2CPU(i) + Mn2CPU(i-1) !Remember, Mn2CPU is cumulative.
   enddo
endif

Mn2CPU(nprocs) = whole_ncells ! The end guy takes care of the rest of the cells unconditionally.

Mn2CPU = Mn2CPU - 1 ! remember, Morton number starts from zero, so it has to be subtracted by 1.

!if(rank==0)print*, 'Mn2CPU after:',Mn2CPU

!--- Specify the Min and Max Morton on each process.
if(rank==0)then
   MinMn=0
else
   MinMn=Mn2CPU(rank)+1
endif
MaxMn=Mn2CPU(rank+1)

deallocate(accum_loops)
deallocate(rdispls)
deallocate(loop)
deallocate(proc_total_loops)
deallocate(loop_final)

nullify(p0)
return
end subroutine load_balance

!--------------------------

!yagi 2011/12/26 remake and added
subroutine set_BufferLength
  use param
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::i

  !=========set MeshLength===========
  LevelSize(-1)=BMeshSize/20 !temporary
  GLevelSize(-1)=((NXB_2+4)*(NYB_2+4)*(NZB_2+4)-Nall_ini/8)
  do i=0,10
     LevelSize(i)=int(dble(LevelSize(i-1))*refineRatio(i-1))*8
     GLevelSize(i)=int(dble(GLevelSize(i-1))*refineRatio(i-1))*8
  end do

  LevelSize(-1)=LevelSize(-1)*20 !temporary

  MeshSizeIni =sum(LevelSize)
  GMeshSizeIni=sum(GLevelSize)
  MeshSize    =MeshSizeIni+GMeshSizeIni
  BMeshBound  =LevelSize(-1)!BMeshSize
  MeshBound   =MeshSizeIni
  !PeshSize    =Nall_ini*npart_per_cell*IonSorts*2*2
  PeshSize    =Nall_ini*npart_per_cell*IonSorts*2*4

  if(rank==0)then
     do i=-1,10
        print '(A,i3,A,f6.3,A,i10,A,i10)',"Lv=",i," refineRatio: ",refineRatio(i),"->LevelSize Mesh: ",LevelSize(i)," GMesh: ",GLevelSize(i)
     enddo

     print *,"BMeshSizeIni=",BMeshSize
     print *,"MeshSizeIni =", MeshSizeIni
     print *,"GMeshSizeIni=",GMeshSizeIni
     print *,"MeshSize    =", MeshSize
     print *,"BMeshBound  =",BMeshBound
     print *,"MeshBound   =",MeshBound
     print *,"PeshSize    =",PeshSize     
  endif
  !===================================
  
  !======= set MPIbufferLength =======
  MPIbufsize=GLevelSize(LvMax)
  DDDbufsize=MeshSizeIni
  if(rank==0)then
     print *,"MPIbufMaxLength=",MPIbufsize
     print *,"DDDbufMaxLength=",DDDbufsize
  endif
  !===================================
end subroutine set_BufferLength

subroutine allocate_Memory
  use param
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::sizeofINT,sizeofDBL,sizeofOct,sizeofPtcl
  integer(kind=4)::INTnum,DBLnum
  integer(kind=4)::errM,errM2,errP,errP2,erriCp
  integer(kind=4)::errD(16),errC(10)

  integer(kind=8)::sizeofMesh
  integer(kind=8)::sizeofPesh
  integer(kind=8)::sizeofiCp
  integer(kind=8)::sizeofMPIbuf,sizeofDDDbuf

  integer(kind=4)::i
!  print *,"Allocating Memories in rank:",rank
  
  call MPI_type_size(MPI_INTEGER,sizeofINT,ierr)
  call MPI_type_size(MPI_DOUBLE_PRECISION,sizeofDBL,ierr)

  !for OCT
  INTnum=32+iomp0
  DBLnum=63+IonSorts
  sizeofOct=INTnum*sizeofINT+DBLnum*sizeofDBL

  !for PTCL
  INTnum=4
  DBLnum=9
  sizeofPtcl=INTnum*sizeofINT+DBLnum*sizeofDBL
  

  sizeofMesh=MeshSize/(1024**2)*sizeofOct*2
  sizeofPesh=PeshSize/(1024**2)*sizeofPtcl*(LvMax+2)
  sizeofiCp =PeshSize/(1024**2)*iomp0*sizeofPtcl
  sizeofMPIbuf=MPIbufsize*12*sizeofInt/(1024**2)&
              +MPIbufsize*12*sizeofDBL/(1024**2)&
              +Nall_ini*2 *npart_per_cell*IonSorts*sizeofInt/(1024**2)&
              +Nall_ini*18*npart_per_cell*IonSorts*sizeofDBL/(1024**2)
  sizeofDDDbuf=DDDbufsize*12*sizeofInt/(1024**2)&
              +DDDbufsize*48*sizeofDBL/(1024**2)&
              +Nall_ini*4 *npart_per_cell*IonSorts*sizeofInt/(1024**2)&
              +Nall_ini*36*npart_per_cell*IonSorts*sizeofDBL/(1024**2)

  if(rank==0)then
     print *,"***************************"
     print '(A,i7,A)',"Oct Size  :",sizeofOct,"Byte"
     print '(A,i7,A)',"Prtcl Size:",sizeofPtcl,"Byte"
     print *,"---------------------------"
     print '(A,i7,A)',"Mesh Size:",sizeofMesh,"MByte"
     print '(A,i7,A)',"Pesh Size:",sizeofPesh,"MByte"
     print '(A,i7,A)',"iCp  Size:",sizeofiCp ,"MByte"
     print *,"---------------------------"
     print '(A,i7,A)',"MPI Buffer Size:",sizeofMPIbuf *2,"MByte"
     print '(A,i7,A)',"DDD Buffer Size:",sizeofDDDbuf   ,"MByte"
     print *,"***************************"
  endif

!------- allocate MPI array -------
  if(debugMode>=3 .and. rank==0)print *,"allocating",nprocs*sizeofINT," B for Mn2CPU"
  allocate(Mn2CPU(nprocs+mempat))
  if(debugMode>=3 .and. rank==0)print *,"allocating",nprocs*sizeofINT," B for Mn2CPU_old"
  allocate(Mn2CPU_old(nprocs+mempat))
  if(debugMode>=3 .and. rank==0)print *,"allocating",nprocs*sizeofINT," B for n_rcells_proc"
  allocate(n_rcells_proc(nprocs+mempat))
  if(debugMode>=3 .and. rank==0)print *,"allocating",nprocs*sizeofINT," B for n_gcells_proc"
  allocate(n_gcells_proc(nprocs+mempat))
  if(debugMode>=3 .and. rank==0)print *,"allocating",nprocs*sizeofINT," B for n_rcells_sum"
  allocate(n_rcells_sum(nprocs+mempat))
  if(debugMode>=3 .and. rank==0)print *,"allocating",nprocs*sizeofINT," B for n_gcells_sum"
  allocate(n_gcells_sum(nprocs+mempat))
!--------- allocate Mesh ----------
  if(debugMode>=3 .and. rank==0)print *,"allocating",MeshSize*sizeofOct/(1024**2),"MB for Mesh"
  allocate(Mesh(MeshSize+mempat),stat=errM)
  if(debugMode>=3 .and. rank==0)print *,"allocating",MeshSize*sizeofOct/(1024**2),"MB for Mesh2"
  allocate(Mesh2(MeshSize+mempat),stat=errM2)
!--------- allocate Pesh ----------
  if(debugMode>=3 .and. rank==0)print *,"allocating",PeshSize*sizeofPtcl*(LvMax+1)/(1024**2),"MB for Pesh"
  allocate(Pesh(PeshSize+mempat,0:LvMax),stat=errP)
  if(debugMode>=3 .and. rank==0)print *,"allocating",PeshSize*sizeofPtcl/(1024**2),"MB for Pesh2"
  allocate(Pesh2(PeshSize+mempat),stat=errP2)
  if(debugMode>=3 .and. rank==0)print *,"allocating",sizeofiCp,"MB for iCp"
  allocate(iCp(iomp0,PeshSize+mempat),stat=erriCp)
!------ allocate MPI buffer -------

  if(debugMode>=3 .and. rank==0)print *,"allocating",MPIbufsize*6*sizeofDBL/(1024**2),"MB for sbuf_BFR"
  allocate(sbuf_BFR(MPIbufsize*6+mempat,2),stat=errC(1))
  if(debugMode>=3 .and. rank==0)print *,"allocating",MPIbufsize*6*sizeofDBL/(1024**2),"MB for sbuf_GFR"
  allocate(rbuf_GFR(MPIbufsize*18+mempat,2),stat=errC(2))
  if(debugMode>=3 .and. rank==0)print *,"allocating",MPIbufsize*18*sizeofDBL/(1024**2),"MB for rbuf_G"
  allocate(rbuf_G(MPIbufsize*18+mempat,2),stat=errC(2))
  if(debugMode>=3 .and. rank==0)print *,"allocating",MPIbufsize*18*sizeofDBL/(1024**2),"MB for sbuf_G"
  allocate(sbuf_G(MPIbufsize*18+mempat,2),stat=errC(2))
  if(debugMode>=3 .and. rank==0)print *,"allocating",MPIbufsize*6*sizeofINT/(1024**2),"MB for sbuf_BFI"
  allocate(sbuf_BFI(MPIbufsize*3+mempat,2),stat=errC(3))
  if(debugMode>=3 .and. rank==0)print *,"allocating",MPIbufsize*6*sizeofINT/(1024**2),"MB for sbuf_GFI"
  allocate(rbuf_GFI(MPIbufsize*3+mempat,2),stat=errC(4))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*9*IonSorts*npart_per_cell*sizeofDBL/(1024**2),"MB for sbuf_GPR"
  allocate(sbuf_GPR(Nall_ini*9*IonSorts*npart_per_cell+mempat,2),stat=errC(5))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*9*IonSorts*npart_per_cell*sizeofDBL/(1024**2),"MB for sbuf_BPR"
  allocate(rbuf_BPR(Nall_ini*9*IonSorts*npart_per_cell+mempat,2),stat=errC(6))
  if(debugMode>=3 .and. rank==0)print *,"allocating",MPIbufsize*3*sizeofINT/(1024**2),"MB for sbuf_GPO"
  allocate(sbuf_GPO(MPIbufsize*3+mempat,2),stat=errC(7))
  if(debugMode>=3 .and. rank==0)print *,"allocating",MPIbufsize*3*sizeofINT/(1024**2),"MB for sbuf_BPO"
  allocate(rbuf_BPO(MPIbufsize*3+mempat,2),stat=errC(8))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*1*IonSorts*npart_per_cell*sizeofINT/(1024**2),"MB for sbuf_GPI"
  allocate(sbuf_GPI(Nall_ini*1*IonSorts*npart_per_cell+mempat,2),stat=errC(9))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*1*IonSorts*npart_per_cell*sizeofINT/(1024**2),"MB for sbuf_BPI"
  allocate(rbuf_BPI(Nall_ini*1*IonSorts*npart_per_cell+mempat,2),stat=errC(10))
!------ allocate DDD buffer ------
  if(debugMode>=3 .and. rank==0)print *,"allocating",DDDbufsize*12*sizeofDBL/(1024**2),"MB for sbufL_FR"
  allocate(sBufL_FR(DDDbufsize*12+mempat),stat=errD( 1))
  if(debugMode>=3 .and. rank==0)print *,"allocating",DDDbufsize*12*sizeofDBL/(1024**2),"MB for sbufR_FR"
  allocate(sBufR_FR(DDDbufsize*12+mempat),stat=errD( 2))
  if(debugMode>=3 .and. rank==0)print *,"allocating",DDDbufsize*12*sizeofDBL/(1024**2),"MB for rbufL_FR"
  allocate(rBufL_FR(DDDbufsize*12+mempat),stat=errD( 3))
  if(debugMode>=3 .and. rank==0)print *,"allocating",DDDbufsize*12*sizeofDBL/(1024**2),"MB for rbufR_FR"
  allocate(rBufR_FR(DDDbufsize*12+mempat),stat=errD( 4))
  if(debugMode>=3 .and. rank==0)print *,"allocating",DDDbufsize*3*sizeofINT/(1024**2),"MB for sbufL_FI"
  allocate(sBufL_FI(DDDbufsize*3 +mempat),stat=errD( 5))
  if(debugMode>=3 .and. rank==0)print *,"allocating",DDDbufsize*3*sizeofINT/(1024**2),"MB for sbufR_FI"
  allocate(sBufR_FI(DDDbufsize*3 +mempat),stat=errD( 6))
  if(debugMode>=3 .and. rank==0)print *,"allocating",DDDbufsize*3*sizeofINT/(1024**2),"MB for rbufL_FI"
  allocate(rBufL_FI(DDDbufsize*3 +mempat),stat=errD( 7))
  if(debugMode>=3 .and. rank==0)print *,"allocating",DDDbufsize*3*sizeofINT/(1024**2),"MB for rbufR_FI"
  allocate(rBufR_FI(DDDbufsize*3 +mempat),stat=errD( 8))

  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*9*IonSorts*npart_per_cell*sizeofDBL/(1024**2),"MB for sbufL_PR"
  allocate(sBufL_PR(Nall_ini*9 *IonSorts*npart_per_cell +mempat),stat=errD( 9))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*9*IonSorts*npart_per_cell*sizeofDBL/(1024**2),"MB for sbufR_PR"
  allocate(sBufR_PR(Nall_ini*9 *IonSorts*npart_per_cell +mempat),stat=errD(10))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*9*IonSorts*npart_per_cell*sizeofDBL/(1024**2),"MB for rbufL_PR"
  allocate(rBufL_PR(Nall_ini*9 *IonSorts*npart_per_cell +mempat),stat=errD(11))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*9*IonSorts*npart_per_cell*sizeofDBL/(1024**2),"MB for rbufR_PR"
  allocate(rBufR_PR(Nall_ini*9 *IonSorts*npart_per_cell +mempat),stat=errD(12))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*IonSorts*npart_per_cell*sizeofINT/(1024**2),"MB for sbufL_PI"
  allocate(sBufL_PI(Nall_ini   *IonSorts*npart_per_cell +mempat),stat=errD(13))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*IonSorts*npart_per_cell*sizeofINT/(1024**2),"MB for sbufR_PI"
  allocate(sBufR_PI(Nall_ini   *IonSorts*npart_per_cell +mempat),stat=errD(14))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*IonSorts*npart_per_cell*sizeofINT/(1024**2),"MB for rbufL_PI"
  allocate(rBufL_PI(Nall_ini   *IonSorts*npart_per_cell +mempat),stat=errD(15))
  if(debugMode>=3 .and. rank==0)print *,"allocating",Nall_ini*IonSorts*npart_per_cell*sizeofINT/(1024**2),"MB for rbufR_PI"
  allocate(rBufR_PI(Nall_ini   *IonSorts*npart_per_cell +mempat),stat=errD(16))

!error check
  if(errM/=0)then
     print *,"Mesh allocation error",rank ; stop
  endif
  if(errM2/=0)then
     print *,"Mesh2 allocation error",rank ; stop
  endif
  if(errP/=0)then
     print *,"Pesh allocation error",rank ; stop
  endif
  if(errP2/=0)then
     print *,"Pesh2 allocation error",rank ; stop
  endif
  if(erriCp/=0)then
     print *,"iCp allocation error",rank ; stop
  endif

  do i=1,10
     if(errC(i)/=0)then
        print *,"Buffer allocation error i=",i,rank ; stop
     endif
  enddo
  do i=1,16
     if(errD(i)/=0)then
        print *,"Buffer allocation error i=",i,rank ; stop
     endif
  enddo

  !if(debugMode>=3 .and. rank==0)print *,"AllocateMemory Completed"

end subroutine allocate_Memory

subroutine deallocate_Memory
  use param
  use init_mesh_size
  use message_passing_interface
  implicit none

  deallocate(rBufR_PI)
  deallocate(rBufL_PI)
  deallocate(sBufR_PI)
  deallocate(sBufL_PI)
  deallocate(rBufR_PR)
  deallocate(rBufL_PR)
  deallocate(sBufR_PR)
  deallocate(sBufL_PR)
  deallocate(rBufR_FI)
  deallocate(rBufL_FI)
  deallocate(sBufR_FI)
  deallocate(sBufL_FI)
  deallocate(rBufR_FR)
  deallocate(rBufL_FR)
  deallocate(sBufR_FR)
  deallocate(sBufL_FR)

  deallocate(rbuf_BPI)
  deallocate(sbuf_GPI)
  deallocate(rbuf_BPO)
  deallocate(sbuf_GPO)
  deallocate(rbuf_BPR)
  deallocate(sbuf_GPR)
  deallocate(rbuf_GFI)
  deallocate(sbuf_BFI)
  deallocate(rbuf_GFR)
  deallocate(sbuf_BFR)
  deallocate(sbuf_G)
  deallocate(rbuf_G)

  deallocate(iCp)
  deallocate(Pesh2)
  deallocate(Pesh)
  deallocate(Mesh2)
  deallocate(Mesh)
  deallocate(n_gcells_proc)
  deallocate(n_rcells_proc)
  deallocate(n_gcells_sum)
  deallocate(n_rcells_sum)
  deallocate(Mn2CPU)

end subroutine deallocate_Memory



!--------------------------------GMesh Handler Routines----------------------------------
!------            GMesh is an interprocess communication buffer
!------            which envelopes each process domain in 3D
!------            and overlaps onto the adjacent process domain.
!------            (This part is written without considering the performance, so it has 
!-----             a lot of room for optimization)

subroutine make_GMesh
  use mpi
  use param
  use message_passing_interface
  use init_mesh_size
  use const
  use oct_set
  use MOD_GMap
  implicit none
  integer(kind=4) :: i, index, index2, ix, iy, iz, ixx, iyy, izz, indexG
  integer(kind=4) :: GMn2Min,GMn2Max,GMrtN,com_range,GMnCount
  integer(kind=4) :: minCPU,maxCPU,minindex,maxindex

  integer(kind=4) :: BMesh_len
! --------------Ingredients for GMesh Structure--------------
        integer(kind=4) :: ich,iLv
        integer(kind=4) :: indexP
        integer(kind=4) :: octN0,octLv,octP,iFLG(3),Csort
        integer(kind=4) :: MrtN
        integer(kind=4) :: n1,n2, Nint, k
        real(kind=8)    :: rPOS(3), R(9)
        real(kind=8)    :: F(18),C(6*iomp0),Z(IonSorts),G(1:18),D(1:3),O(1:16)
        integer(kind=4) :: iC(iomp0)
        integer(kind=4) :: octType,prc_bndry, ptcl_loops
!-------------------------------------------

  integer(kind=4) :: Mn, ngst,debug
  integer(kind=4) :: iPOS(3), iPOSH(3), iPOSN(3), iPOSN_center(3), spt(3), ept(3)
  integer(kind=4) :: errG, CPU
  integer(kind=4) :: sdestl,rdestl
  real(kind=8)    :: rx,ry,rz,rr,mr
!  integer(kind=4) :: Gst_Cube((2*buf_width+1)**3)
  integer(kind=4)::Mnarr(27),nbindx(6)
  integer(kind=4), allocatable::n_rflag(:),n_rindex(:),iGMn2OctN(:)
  type(GMap_type), allocatable:: GMap(:)
  type(oct):: newp
  type(oct), pointer::p0,p1
  type(prtcl), pointer :: PrtList

  
!--------                 Obtain minimum and maximum iPos in the process  and create GMap ---------                                      ------
!--------           GMap is a structure that consists of iPOS and Mn number   ------
!-------- GMap acts like a "workbench" on which ghost cells and real cells are classified and recorded ----- 
!-------- GMap is lined up like conventional "raster scan sequence"  with incremental ix -> iy -> iz ------

if(debugMode>=1)print*,'**************Entering make_GMesh',rank
minindex=MinID(1,-1)
maxindex=MaxID(1,-1)

!allocate CPUarr
BMesh_len=maxindex-minindex+1
com_range=(buf_width*2+1)**3

errG=0
allocate(n_rflag(nprocs),stat=errG)
if(errG/=0)then
     print *,"Not enough memory! error no:",errG
     stop
endif
allocate(n_rindex(nprocs),stat=errG)
if(errG/=0)then
     print *,"Not enough memory! error no:",errG
     stop
endif
allocate(GMap(BMesh_len*com_range),stat=errG)
if(errG/=0)then
     print *,"Not enough memory! error no:",errG
     stop
endif

n_rcells_proc=0
n_gcells_proc=0
n_rcells_sum=0
n_gcells_sum=0

do i=1,BMesh_len*com_range
   GMap(i)%MrtN=-1
   GMap(i)%iPOSN=-1
   GMap(i)%proc=-1
   GMap(i)%octN=-1
enddo

GMn2Min=NXB*NYB*NZB !maximum number of morton is NXB*NYB*NZB/8
GMn2Max=-1          !minimum number of morton is 0

ngst=0
!ngst2=0
GMeshSize=0

!get Morton numeber of new ghost oct
do index=minindex,maxindex
   index2=(index-1)*com_range
   p0=>Mesh(index)
   ngst=0
   n_rflag=0
   
   iPOSH=p0%iPOS
   iPOSN_center=(iPOSH + LvMax2)/(2*LvMax2)

   spt=iPOSN_center-buf_width
   ept=iPOSN_center+buf_width

   do iz=spt(3),ept(3)
      do iy=spt(2),ept(2)
         do ix=spt(1),ept(1)
            ixx=ix
            iyy=iy
            izz=iz

            if(ixx<1)then
               ixx=ixx+NX_2t
            elseif(ixx>NX_2t)then
               ixx=ixx-NX_2t
            endif
            if(iyy<1)then
               iyy=iyy+NY_2t
            elseif(iyy>NY_2t)then
               iyy=iyy-NY_2t
            endif
            if(izz<1)then
               izz=izz+NZ_2t
            elseif(izz>NZ_2t)then
               izz=izz-NZ_2t
            endif

            iPOSN(1)=ixx
            iPOSN(2)=iyy
            iPOSN(3)=izz

            call get_indexNumberN(iPOSN,GMrtN)
            if(MinMn<=GMrtN .and. GMrtN<=MaxMn)then
               !target is oct. then do nothing.
            else
               GMn2Min=min(GMn2Min,GMrtN)
               GMn2Max=max(GMn2Max,GMrtN)
               ngst=ngst+1
               !ngst2=ngst2+1
               call findCPU(GMrtN,CPU)
               GMap(index2+ngst)%MrtN=GMrtN
               GMap(index2+ngst)%iPOSN=iPOSN
               GMap(index2+ngst)%proc=CPU+1
               n_rflag(CPU+1)=1
            endif
         end do
      end do
   end do
   
   if(ngst>0)then
      p0%prc_bndry=1
      if(associated(p0%octNb1).and.associated(p0%octNb2).and.associated(p0%octNb3).and.&
         associated(p0%octNb4).and.associated(p0%octNb5).and.associated(p0%octNb6)     )then
!		 print *,"Error"
      else
         p0%prc_bndry=2
      endif
   endif

   n_rcells_proc(1:nprocs)=n_rcells_proc(1:nprocs)+n_rflag(1:nprocs)
enddo

allocate(iBMesh_arr(sum(n_rcells_proc)),stat=errG)
if(errG/=0)then
     print *,"Not enough memory! error no:",errG
     stop
endif
allocate(GMn2OctN(GMn2Min:GMn2Max),stat=errG)
if(errG/=0)then
     print *,"Not enough memory! error no:",errG
     stop
endif
allocate(iGMn2OctN(BMesh_len*com_range),stat=errG)
if(errG/=0)then
     print *,"Not enough memory! error no:",errG
     stop
endif

!Create n_rcells_sum
do CPU=1,nprocs
   if(CPU==1)then
      n_rcells_sum(CPU)=0
   else
      n_rcells_sum(CPU)=sum(n_rcells_proc(1:CPU-1))
   endif
enddo
!!$print *,"rank=",rank,"n_rcells_proc=",n_rcells_proc
!!$print *,"rank=",rank,"n_rcells_sum =",n_rcells_sum

!print *,"size iBM=",sum(n_rcells_proc),rank
GMn2OctN=-1
iGMn2OctN=-1
iBMesh_arr=-1
n_rflag=0

!prepare parameter for new Goct
nullify(newp%octPrt)
nullify(newp%octNb1) ; nullify(newp%octNb2)
nullify(newp%octNb3) ; nullify(newp%octNb4)
nullify(newp%octNb5) ; nullify(newp%octNb6)
nullify(newp%octCh1) ; nullify(newp%octCh2)
nullify(newp%octCh3) ; nullify(newp%octCh4)
nullify(newp%octCh5) ; nullify(newp%octCh6)
nullify(newp%octCh7) ; nullify(newp%octCh8)
nullify(newp%Psort)
nullify(newp%ptcl)   ; nullify(newp%ptclA)
octLv = -1
octP = 0
iFLG = (/4,0,0/) !iFLG(1)=4 sets Lv=0 below, and iFLG(2) calls for creation of child cells
Csort = 0
C = 0d0
F=0.d0
Z=0.d0
G=0.d0
D=0.d0
O=0.d0
iC=0
octN0 = 0
ptcl_loops=0
prc_bndry=0
octType=1 !for GMesh

GMnCount=1
n_rindex=0

!set GOct and GMn2OctN
indexG=MeshBound+1
!print *,"indexG=",indexG,rank
debug=0
do index=minindex,maxindex
   p0=>Mesh(index)
   minCPU=(index-1)*com_range+1
   maxCPU=index*com_range

   n_rflag(1:nprocs)=0

   if(p0%prc_bndry>=1)then
      do index2=minCPU,maxCPU

         CPU=Gmap(index2)%proc
         if(CPU<=0) cycle
         !create iBMesh_arr
         if(n_rflag(CPU)==0)then
            iBMesh_arr(n_rcells_sum(CPU)+n_rindex(CPU)+1)=index
            n_rindex(CPU)=n_rindex(CPU)+1
            n_rflag(CPU)=1
         endif

         if(p0%prc_bndry>=2)then
            GMrtN=GMap(index2)%MrtN
            if(GMrtN>=0)then
               if(GMn2OctN(GMrtN)<0)then
                  GMn2OctN(GMrtN)=indexG
                  iGMn2OctN(GMnCount)=GMrtN
                  GMnCount=GMnCount+1

                  Gmap(index2)%octN=indexG

                  !add number to n_gcells_proc
                  n_gcells_proc(CPU)=n_gcells_proc(CPU)+1

                  !add Lv-1 Goct
                  iPOSN=Gmap(index2)%iPOSN
                  rPOS(:) = real(iPOSN(:))*dx(:)*2d0 - dx(:)
                  iPOSH = 2*LvMax2*iPOSN - LvMax2
                  octN0=indexG
                  prc_bndry=CPU !!!
                  MrtN=GMrtN
                  if(dipoleFlag==1)then
                    rx=rPOS(1)-Dpos(1)
                    ry=rPOS(2)-Dpos(2)
                    rz=rPOS(3)-Dpos(3)
                    rr=sqrt(rx**2+ry**2+rz**2)
                    if(rr.gt.dx(1)*2.0d0) then
                        mr=mx*rx+my*ry+mz*rz
                        D(1)=(-mx/(rr**3)+3*rx*mr/(rr**5))/(4*PI)
                        D(2)=(-my/(rr**3)+3*ry*mr/(rr**5))/(4*PI)
                        D(3)=(-mz/(rr**3)+3*rz*mr/(rr**5))/(4*PI)
                    endif
                  end if
                  Mesh(indexG) =                            &
                        oct(octN0,octLv,octP,Csort,iFLG,    &
                        iPOSH,rPOS,                         &
                        MrtN,iC,                            &
                        prc_bndry,                          &
                        octType,                            &
                        F,C,Z,G,D,O,ptcl_loops,       & 
                        newp%octPrt,                        &
                        newp%octNb1, newp%octNb2,           &
                        newp%octNb3, newp%octNb4,           &
                        newp%octNb5, newp%octNb6,           &
                        newp%octCh1, newp%octCh2,           &
                        newp%octCh3, newp%octCh4,           &
                        newp%octCh5, newp%octCh6,           &
                        newp%octCh7, newp%octCh8,           &
                        newp%Psort ,                        &
                        newp%ptcl  , newp%ptclA)  

                  indexG=indexG+1
                  GMeshSize=GMeshSize+1
               endif
            endif
         endif
      enddo   
   endif
enddo

!print *,"debug",debug,rank
!print *,"n_rindex iBM=",n_rindex,rank
!!$print *,"GMeshSize=",GMeshSize,rank
!!$print *,"indexGafter =",indexG,rank
!set Min Max ID
MinID(3,-1)=MeshBound+1
MaxID(3,-1)=MeshBound+GMeshSize
MinID(3,0) =MaxID(3,-1)+1
MaxID(3,0) =MaxID(3,-1)+8*GMeshSize
if(LvMax>0)then
   do iLv=1,LvMax
      MinID(3,iLv) = MaxID(3,iLv-1)
      MaxID(3,iLv) = MinID(3,iLv)
   enddo
endif
do iLv=-1, LvMax
   MinID(4,iLv) = MaxID(3,LvMax)
   MaxID(4,iLv) = MaxID(3,LvMax)
enddo

!create iGMesh_arr
allocate(iGMesh_arr(sum(n_gcells_proc)),stat=errG)
if(errG/=0)then
     print *,"Not enough memory! error no:",errG
     stop
endif
iGMesh_arr=-1

!create n_gcells_sum
do CPU=1,nprocs
   if(CPU==1)then
      n_gcells_sum(CPU)=0
   else
      n_gcells_sum(CPU)=sum(n_gcells_proc(1:CPU-1))
   endif
enddo
!!$print *,"rank=",rank,"n_gcells_proc=",n_gcells_proc
!!$print *,"rank=",rank,"n_gcells_sum=",n_gcells_sum

n_rindex=0
GMnCount=GMnCount-1
do index2=1,GMnCount
   GMrtN=iGMn2OctN(index2)
   index=GMn2OctN(GMrtN)

   CPU=Mesh(index)%prc_bndry

   iGMesh_arr(1+n_gcells_sum(CPU)+n_rindex(CPU))=index
   n_rindex(CPU)=n_rindex(CPU)+1   
enddo

!!$do GMrtN=GMn2Min,GMn2Max
!!$   index=GMn2OctN(GMrtN)
!!$   if(index>0)then
!!$      CPU=Mesh(index)%prc_bndry
!!$
!!$      iGMesh_arr(1+n_gcells_sum(CPU)+n_rindex(CPU))=index
!!$      n_rindex(CPU)=n_rindex(CPU)+1
!!$   endif
!!$enddo

!debug
!if(istep==9)then
!!$if(istep==30 .and. rank==6)then
!!$   print *,"rank=",rank,"n_rcells_proc=",n_rcells_proc(1:nprocs)
!!$   print *,"rank=",rank,"n_gcells_proc=",n_gcells_proc(1:nprocs)
!!$   print *,"rank=",rank,"n_rindex=",n_rindex(1:nprocs)
!!$   print *,"rank=",rank,"iBMesh_arr",iBMesh_arr(1:sum(n_rcells_proc))
!!$   print *,"rank=",rank,"size iBM=",sum(n_rcells_proc)
!!$   print *,"rank=",rank,"iGMesh_arr",iGMesh_arr(1:sum(n_gcells_proc))
!!$   print *,"rank=",rank,"size iGM=",sum(n_gcells_proc),rank
!!$   print *,"rank=",rank,"MinMn=",MinMn,"MaxMn",MaxMn,"GMn2Min",GMn2Min,"GMn2Max",GMn2Max
!!$   print *,"rank=",rank,"MnDisp=",MnDisp,"MaxMn",MaxMn,"GMn2Min",GMn2Min,"GMn2Max",GMn2Max
!!$endif
!endif

!connect Lv-1 Goct and Goct or Oct
maxindex=MaxID(3,-1)
minindex=MinID(3,-1)
do index=minindex,maxindex
   p0 => Mesh(index)
   iPOSH = p0%iPOS

   call get_neighbourIndex(iPOSH,Mnarr)
   iPOSN = (iPOSN+LvMax2)/(2*LvMax2)
   nbindx=-1
   if(Mnarr(1)>=GMn2Min .and. Mnarr(1)<=GMn2Max)nbindx(1)=GMn2OctN(Mnarr(1))
   if(Mnarr(4)>=GMn2Min .and. Mnarr(4)<=GMn2Max)nbindx(2)=GMn2OctN(Mnarr(4))
   if(Mnarr(2)>=GMn2Min .and. Mnarr(2)<=GMn2Max)nbindx(3)=GMn2OctN(Mnarr(2))
   if(Mnarr(5)>=GMn2Min .and. Mnarr(5)<=GMn2Max)nbindx(4)=GMn2OctN(Mnarr(5))
   if(Mnarr(3)>=GMn2Min .and. Mnarr(3)<=GMn2Max)nbindx(5)=GMn2OctN(Mnarr(3))
   if(Mnarr(6)>=GMn2Min .and. Mnarr(6)<=GMn2Max)nbindx(6)=GMn2OctN(Mnarr(6))

   !if(p0%MrtN==49)print *,"Mnarr",Mnarr(1:6),"nbindx",nbindx(1:6),rank

   if(nbindx(1)>=minindex .and. nbindx(1)<=maxindex)then !nb is ghost
      p0%octNb1=>Mesh(nbindx(1))
   else if(Mnarr(1)>=MinMn .and. Mnarr(1)<=MaxMn)then !nb is oct
      p1=>Mesh(Mnarr(1)-MnDisp)
      p0%octNb1=>p1
      p1%octNb2=>p0 
   else
      nullify(p0%octNb1)
   endif
   if(nbindx(2)>=minindex .and. nbindx(2)<=maxindex)then
      p0%octNb2=>Mesh(nbindx(2))
   else if(Mnarr(4)>=MinMn .and. Mnarr(4)<=MaxMn)then
      p1=>Mesh(Mnarr(4)-MnDisp)
      p0%octNb2=>p1
      p1%octNb1=>p0
   else
      nullify(p0%octNb2)
   endif
   if(nbindx(3)>=minindex .and. nbindx(3)<=maxindex)then
      p0%octNb3=>Mesh(nbindx(3))
   else if(Mnarr(2)>=MinMn .and. Mnarr(2)<=MaxMn)then
      p1=>Mesh(Mnarr(2)-MnDisp)
      p0%octNb3=>p1
      p1%octNb4=>p0
   else
      nullify(p0%octNb3)
   endif
   if(nbindx(4)>=minindex .and. nbindx(4)<=maxindex)then
      p0%octNb4=>Mesh(nbindx(4))
   else if(Mnarr(5)>=MinMn .and. Mnarr(5)<=MaxMn)then
      p1=>Mesh(Mnarr(5)-MnDisp)
      p0%octNb4=>p1
      p1%octNb3=>p0
   else
      nullify(p0%octNb4)
   endif
   if(nbindx(5)>=minindex .and. nbindx(5)<=maxindex)then
      p0%octNb5=>Mesh(nbindx(5))
   else if(Mnarr(3)>=MinMn .and. Mnarr(3)<=MaxMn)then
      p1=>Mesh(Mnarr(3)-MnDisp)
      p0%octNb5=>p1
      p1%octNb6=>p0
   else
      nullify(p0%octNb5)
   endif
   if(nbindx(6)>=minindex .and. nbindx(6)<=maxindex)then
      p0%octNb6=>Mesh(nbindx(6))
   else if(Mnarr(6)>=MinMn .and. Mnarr(6)<=MaxMn)then
      p1=>Mesh(Mnarr(6)-MnDisp)
      p0%octNb6=>p1
      p1%octNb5=>p0
   else
      nullify(p0%octNb6)
   endif
enddo

!set Lv0 Goct
indexP = MaxIP(0)
iFLG(:)= (/0,0,0/)  
octLv   = 0
Nint  = 2**(LvMax - octLv)
nullify(PrtList)
!print *,"MinID(3,:)=",MinID(3,:),"MaxID(3,:)=",MaxID(3,:),rank
do index=MinID(3,-1),MaxID(3,-1)
   p0 => Mesh(index)

   index2 = MinID(3,octLv) + (index-MinID(3,octLv-1))* 8 
   !if(rank==0)print *,"index2=",index2,rank
   ! -- each directed child octs
   do ich=1,8
      nullify(newp%octPrt)
      nullify(newp%octNb1) ; nullify(newp%octNb2)
      nullify(newp%octNb3) ; nullify(newp%octNb4)
      nullify(newp%octNb5) ; nullify(newp%octNb6)
      nullify(newp%octCh1) ; nullify(newp%octCh2)
      nullify(newp%octCh3) ; nullify(newp%octCh4)
      nullify(newp%octCh5) ; nullify(newp%octCh6)
      nullify(newp%octCh7) ; nullify(newp%octCh8)
      nullify(newp%Psort)
      nullify(newp%ptcl)   ; nullify(newp%ptclA)

      Mn = p0 % MrtN
      !if(rank==1)print *,"Mn=",Mn,"index=",index
      ! -- Extract position index for child-octs --
      n1 = ich
      iz = int((n1-1)/4) + 1
      n2 = n1 - 4*(iz-1)
      iy = int((n2-1)/2) + 1
      ix = n2 - 2*(iy-1)
      ! -- Position index --
      ix = p0%iPOS(1) + (2*ix-3)*Nint
      iy = p0%iPOS(2) + (2*iy-3)*Nint
      iz = p0%iPOS(3) + (2*iz-3)*Nint
      iPOS(1) = ix
      iPOS(2) = iy
      iPOS(3) = iz
      ! -- Position --
      ixx = real(ix + Nint)/real(2*Nint) !Converting from H to N
      iyy = real(iy + Nint)/real(2*Nint)
      izz = real(iz + Nint)/real(2*Nint)
      rPOS(1) = dx(1)*(ixx - HALF)* 2d0**(-octLv)
      rPOS(2) = dx(2)*(iyy - HALF)* 2d0**(-octLv)
      rPOS(3) = dx(3)*(izz - HALF)* 2d0**(-octLv)

      R(1:3)  = rPOS(:)
      R(4:9) = 0d0

      ! -- Field values input --
      do k=1,12
         F(k) = p0%F(k)
      enddo

      do k=1,IonSorts
         Z(k) = p0%Z(k)
      enddo
      if(dipoleFlag==1)then
        rx=rPOS(1)-Dpos(1)
        ry=rPOS(2)-Dpos(2)
        rz=rPOS(3)-Dpos(3)
        rr=sqrt(rx**2+ry**2+rz**2)
        if(rr.gt.dx(1)*2.0d0) then
            mr=mx*rx+my*ry+mz*rz
            D(1)=(-mx/(rr**3)+3*rx*mr/(rr**5))/(4*PI)
            D(2)=(-my/(rr**3)+3*ry*mr/(rr**5))/(4*PI)
            D(3)=(-mz/(rr**3)+3*rz*mr/(rr**5))/(4*PI)
        endif
      end if
      C(:) = 0d0

      prc_bndry = p0% prc_bndry

      octP = 0
      Csort = ich
      octN0 = index2
      octType = 1 !This is GMesh
      ! -- Connect parent oct & next oct --
      !                        newp%octPrt => p0
      newp%Psort  => p0

      ! -- representative particle 
      indexP = indexP + 1
      Pesh(indexP,octLv) =                  &
           prtcl( R ,-1,octN0,indexP,PrtList,0)
      newp%ptcl => Pesh(indexP,octLv)


      ! -- Input informaions into a target variable Mesh(:) --
      Mesh(index2) =                       &
           oct(octN0,octLv,octP,Csort,iFLG, &
           iPOS,rPOS,                       &
           Mn,iC,                           &
           prc_bndry,                       &
           octType,                         &
           F,C,Z,G,D,O,ptcl_loops,    & 
           newp%octPrt,                     &
           newp%octNb1, newp%octNb2,        &
           newp%octNb3, newp%octNb4,        &
           newp%octNb5, newp%octNb6,        &
           newp%octCh1, newp%octCh2,        &
           newp%octCh3, newp%octCh4,        &
           newp%octCh5, newp%octCh6,        &
           newp%octCh7, newp%octCh8,        &
           newp%Psort ,                     &
           newp%ptcl  , newp%ptclA)

      ! -- Connect eight child-octs of parent oct --
      if(ich==1) then
         p0%octCh1 => Mesh(index2)
      elseif(ich==2) then
         p0%octCh2 => Mesh(index2)
      elseif(ich==3) then
         p0%octCh3 => Mesh(index2)
      elseif(ich==4) then
         p0%octCh4 => Mesh(index2)
      elseif(ich==5) then
         p0%octCh5 => Mesh(index2)
      elseif(ich==6) then
         p0%octCh6 => Mesh(index2)
      elseif(ich==7) then
         p0%octCh7 => Mesh(index2)
      elseif(ich==8) then 
         p0%octCh8 => Mesh(index2)
      endif
      Mesh(index2) % octPrt => p0
      index2 = index2+1
   enddo  ! for do ich=1,8
   p0%iFLG(2) = 1 !Creation of 8 children are completed
enddo

MaxIP(0)=indexP

do i=MinID(1,-1),MaxID(1,-1)
   p0 => Mesh(i)

   if(.not. associated(p0 % octNb1) .or. .not. associated(p0 % octNb2) .or. &
        .not. associated(p0 % octNb3) .or. .not. associated(p0 % octNb4) .or. &
        .not. associated(p0 % octNb5) .or. .not. associated(p0 % octNb6)) then
      print*, "This BMesh is not connected to the neighbour",rank, p0%octN,&
           (p0% iPOS +LvMax2)/(2*LvMax2),&
           associated(p0 % octNb1),associated(p0 % octNb2), &
           associated(p0 % octNb3),associated(p0 % octNb4), &
           associated(p0 % octNb5),associated(p0 % octNb6)
   endif

enddo
! connect Lv0 octs
call connect_G_n_M

deallocate(n_rindex)
deallocate(n_rflag)
deallocate(GMap)



!create send/recv_list for refresh routines
send_glist=0
send_rlist=0
recv_glist=0
recv_rlist=0
sglistnum=0
srlistnum=0
rglistnum=0
rrlistnum=0

do i=1,nprocs-1
   sdestl=mod(rank+i,nprocs)
   rdestl=mod(rank-i+nprocs,nprocs)
 
   sglistnum=sglistnum+1
   send_glist(sglistnum)=sdestl

   srlistnum=srlistnum+1
   send_rlist(srlistnum)=sdestl

   rglistnum=rglistnum+1
   recv_glist(rglistnum)=rdestl

   rrlistnum=rrlistnum+1
   recv_rlist(rrlistnum)=rdestl

enddo

if(sglistnum/=rrlistnum .or. srlistnum/=rglistnum)then
   print *,"send_ghost_num and recv_oct_num is not equal",sglistnum,rrlistnum,rank
   print *,"send_oct_num and recv_ghost_num is not equal",srlistnum,rglistnum,rank
endif
 
return 
end subroutine make_GMesh

subroutine reset_Gst
  use message_passing_interface
  use param
  use oct_set
  implicit none
  integer(kind=4)::iLv,index,i
  type(oct),pointer::p0
  type(prtcl),pointer::pptc

 ! if(debugMode>=1)print *,"starting reset_Gst",rank

  do iLv=-1,LvMax
     do i=3,4
        if(MinID(i,iLv)>=MaxID(i,iLv))cycle
        do index=MinID(i,iLv),MaxID(i,iLv)

           p0=>Mesh(index)
           !Delete particles
           if(p0%octP/=0)then
              print *,"Error in reset_Gst : Particles are remaining in Goct"
              stop
           endif
           p0%octP=0
           if(iLv >=0) then  !Remember there is no rep particle in iLv=-1
              pptc => p0%ptcl
              pptc%isort = 0
              nullify(p0%ptcl)
           endif

          
           if(associated(p0%octNb1))nullify(p0%octNb1%octNb2)
           if(associated(p0%octNb2))nullify(p0%octNb2%octNb1)
           if(associated(p0%octNb3))nullify(p0%octNb3%octNb4)
           if(associated(p0%octNb4))nullify(p0%octNb4%octNb3)
           if(associated(p0%octNb5))nullify(p0%octNb5%octNb6)
           if(associated(p0%octNb6))nullify(p0%octNb6%octNb5)

           nullify(p0%octCh1)
           nullify(p0%octCh2)
           nullify(p0%octCh3)
           nullify(p0%octCh4)
           nullify(p0%octCh5)
           nullify(p0%octCh6)
           nullify(p0%octCh7)
           nullify(p0%octCh8)
           nullify(p0%octPrt)
           nullify(p0%ptcl)
        enddo
     enddo
  enddo

  !set MinMaxID
  MinID(3,-1)=MeshBound+1
  MaxID(3,-1)=MeshBound
  MinID(3,0) =MinID(3,-1)+1
  MaxID(3,0) =MinID(3,-1)
  if(LvMax>0)then
     do iLv=1,LvMax
        MinID(3,iLv) = MaxID(3,iLv-1)+1
        MaxID(3,iLv) = MaxID(3,iLv-1)
     end do
  endif
  do iLv=-1, LvMax
     MinID(4,iLv) = MaxID(3,LvMax)+1
     MaxID(4,iLv) = MaxID(3,LvMax)
  enddo

!  MaxID(0) = MaxID(0)-(GMaxID(0)-GMinID(0)) This sentence has to be implemented after sort_particle works

  deallocate(iBMesh_arr)
  deallocate(iGMesh_arr)
  deallocate(GMn2octN)

 ! if(debugMode>=1)print *,"exiting reset_Gst",rank

end subroutine reset_Gst

!------------ Interprocess comm of Fields for time evolution ------------

!new communication routine
recursive subroutine pack_GMeshJ(index,sendLv,count,iLv,icon)
  use param
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4),intent(in)::index,sendLv,iLv,icon
  integer(kind=4),intent(inout)::count
  type(oct),pointer::p0

  p0=>Mesh(index)
  if(p0%iFLG(1)>0)then
     if(iLv==sendLv)then
        rbuf_GFR(1+(count-1)*12:count*12,icon)=p0%octCh1%F(7:18)
        rbuf_GFI(   count               ,icon)=p0%octCh1%MrtN
                                               p0%octCh1%F(7:18)=0.0d0
        count=count+1
        rbuf_GFR(1+(count-1)*12:count*12,icon)=p0%octCh2%F(7:18)
        rbuf_GFI(   count               ,icon)=p0%octCh2%MrtN
                                               p0%octCh2%F(7:18)=0.0d0
        count=count+1
        rbuf_GFR(1+(count-1)*12:count*12,icon)=p0%octCh3%F(7:18)
        rbuf_GFI(   count               ,icon)=p0%octCh3%MrtN
                                               p0%octCh3%F(7:18)=0.0d0
        count=count+1
        rbuf_GFR(1+(count-1)*12:count*12,icon)=p0%octCh4%F(7:18)
        rbuf_GFI(   count               ,icon)=p0%octCh4%MrtN
                                               p0%octCh4%F(7:18)=0.0d0
        count=count+1
        rbuf_GFR(1+(count-1)*12:count*12,icon)=p0%octCh5%F(7:18)
        rbuf_GFI(   count               ,icon)=p0%octCh5%MrtN
                                               p0%octCh5%F(7:18)=0.0d0
        count=count+1
        rbuf_GFR(1+(count-1)*12:count*12,icon)=p0%octCh6%F(7:18)
        rbuf_GFI(   count               ,icon)=p0%octCh6%MrtN
                                               p0%octCh6%F(7:18)=0.0d0
        count=count+1
        rbuf_GFR(1+(count-1)*12:count*12,icon)=p0%octCh7%F(7:18)
        rbuf_GFI(   count               ,icon)=p0%octCh7%MrtN
                                               p0%octCh7%F(7:18)=0.0d0
        count=count+1
        rbuf_GFR(1+(count-1)*12:count*12,icon)=p0%octCh8%F(7:18)
        rbuf_GFI(   count               ,icon)=p0%octCh8%MrtN
                                               p0%octCh8%F(7:18)=0.0d0
        count=count+1

     else
        call pack_GMeshJ(p0%octCh1%octN,sendLv,count,iLv+1,icon)
        call pack_GMeshJ(p0%octCh2%octN,sendLv,count,iLv+1,icon)
        call pack_GMeshJ(p0%octCh3%octN,sendLv,count,iLv+1,icon)
        call pack_GMeshJ(p0%octCh4%octN,sendLv,count,iLv+1,icon)
        call pack_GMeshJ(p0%octCh5%octN,sendLv,count,iLv+1,icon)
        call pack_GMeshJ(p0%octCh6%octN,sendLv,count,iLv+1,icon)
        call pack_GMeshJ(p0%octCh7%octN,sendLv,count,iLv+1,icon)
        call pack_GMeshJ(p0%octCh8%octN,sendLv,count,iLv+1,icon)
     endif
  endif
end subroutine pack_GMeshJ

subroutine pack_GMeshJT(sdest,sendLv,icon)
  use param
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4),intent(in)::sdest,sendLv,icon
  integer(kind=4)::count
  integer(kind=4)::sID,eID,index,i,iLv

  if(MinID(3,sendLv)>=MaxID(3,sendLv))return
  if(n_gcells_proc(sdest+1)==0)return
  if(sendLv==-1)return

  if(sdest==0)then
     sID=1
  else
     sID=sum(n_gcells_proc(1:sdest))+1
  endif
  eID=sum(n_gcells_proc(1:sdest+1))

  count=1
  iLv  =-1

  do i=sID,eID
     index=iGMesh_arr(i)
     call pack_GMeshJ(index,sendLv,count,iLv+1,icon)
  enddo

  send_onum(icon)=count-1

end subroutine pack_GMeshJT

recursive subroutine pack_BMeshEB(index,sendLv,count,icon,iLv,bicon)
  use param
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::count   !count and count_p must be 1 when this subroutine called
  integer(kind=4),intent(in)::index,icon,iLv,sendLv,bicon
  type(oct),pointer::p0

  p0=>Mesh(index)
  if(p0%iFLG(1)>0)then
     if(iLv==sendLv)then
        sbuf_BFR(1+(count-1)*3:count*3,bicon)=p0%octCh1%F(1+(icon-1)*3:icon*3)
        sbuf_BFI(   count             ,bicon)=p0%octCh1%MrtN
        count=count+1
        sbuf_BFR(1+(count-1)*3:count*3,bicon)=p0%octCh2%F(1+(icon-1)*3:icon*3)
        sbuf_BFI(   count             ,bicon)=p0%octCh2%MrtN
        count=count+1
        sbuf_BFR(1+(count-1)*3:count*3,bicon)=p0%octCh3%F(1+(icon-1)*3:icon*3)
        sbuf_BFI(   count             ,bicon)=p0%octCh3%MrtN
        count=count+1
        sbuf_BFR(1+(count-1)*3:count*3,bicon)=p0%octCh4%F(1+(icon-1)*3:icon*3)
        sbuf_BFI(   count             ,bicon)=p0%octCh4%MrtN
        count=count+1
        sbuf_BFR(1+(count-1)*3:count*3,bicon)=p0%octCh5%F(1+(icon-1)*3:icon*3)
        sbuf_BFI(   count             ,bicon)=p0%octCh5%MrtN
        count=count+1
        sbuf_BFR(1+(count-1)*3:count*3,bicon)=p0%octCh6%F(1+(icon-1)*3:icon*3)
        sbuf_BFI(   count             ,bicon)=p0%octCh6%MrtN
        count=count+1
        sbuf_BFR(1+(count-1)*3:count*3,bicon)=p0%octCh7%F(1+(icon-1)*3:icon*3)
        sbuf_BFI(   count             ,bicon)=p0%octCh7%MrtN
        count=count+1
        sbuf_BFR(1+(count-1)*3:count*3,bicon)=p0%octCh8%F(1+(icon-1)*3:icon*3)
        sbuf_BFI(   count             ,bicon)=p0%octCh8%MrtN
        count=count+1
     else
        call pack_BMeshEB(p0%octCh1%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshEB(p0%octCh2%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshEB(p0%octCh3%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshEB(p0%octCh4%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshEB(p0%octCh5%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshEB(p0%octCh6%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshEB(p0%octCh7%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshEB(p0%octCh8%octN,sendLv,count,icon,iLv+1,bicon)
     endif
  endif
end subroutine pack_BMeshEB

subroutine pack_BMeshEBT(sdest,sendLv,icon,bicon)
  use param
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4),intent(in)::sdest,icon,sendLv,bicon
  integer(kind=4)::i,index,count,sID,eID,iLv

  if(MinID(1,sendLv)>=MaxID(1,sendLv))return
  if(n_rcells_proc(sdest+1)==0)return
  if(sendLv==-1)return

  !set loop index
  if(sdest==0)then
     sID=1
  else
     sID=sum(n_rcells_proc(1:sdest))+1
  endif
  eID=sum(n_rcells_proc(1:sdest+1))

  count=1
  iLv  =-1
  do i=sID,eID
     index=iBMesh_arr(i)
     call pack_BMeshEB(index,sendLv,count,icon,iLv+1,bicon)
  enddo

  send_onum(bicon)=count-1

end subroutine pack_BMeshEBT


recursive subroutine pack_BMeshG(index,sendLv,count,icon,iLv,bicon)
  use param
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::count   !count and count_p must be 1 when this subroutine called
  integer(kind=4),intent(in)::index,icon,iLv,sendLv,bicon
  type(oct),pointer::p0

  p0=>Mesh(index)
  if(p0%iFLG(1)>0)then
     if(iLv==sendLv)then
        sbuf_G(1+(count-1)*18:count*18,bicon)=p0%octCh1%G(1:18)
        count=count+1
        sbuf_G(1+(count-1)*18:count*18,bicon)=p0%octCh2%G(1:18)
        count=count+1
        sbuf_G(1+(count-1)*18:count*18,bicon)=p0%octCh3%G(1:18)
        count=count+1
        sbuf_G(1+(count-1)*18:count*18,bicon)=p0%octCh4%G(1:18)
        count=count+1
        sbuf_G(1+(count-1)*18:count*18,bicon)=p0%octCh5%G(1:18)
        count=count+1
        sbuf_G(1+(count-1)*18:count*18,bicon)=p0%octCh6%G(1:18)
        count=count+1
        sbuf_G(1+(count-1)*18:count*18,bicon)=p0%octCh7%G(1:18)
        count=count+1
        sbuf_G(1+(count-1)*18:count*18,bicon)=p0%octCh8%G(1:18)
        count=count+1
     else
        call pack_BMeshG(p0%octCh1%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshG(p0%octCh2%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshG(p0%octCh3%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshG(p0%octCh4%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshG(p0%octCh5%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshG(p0%octCh6%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshG(p0%octCh7%octN,sendLv,count,icon,iLv+1,bicon)
        call pack_BMeshG(p0%octCh8%octN,sendLv,count,icon,iLv+1,bicon)
     endif
  endif
end subroutine pack_BMeshG

subroutine pack_BMeshGT(sdest,sendLv,icon,bicon)
  use param
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4),intent(in)::sdest,icon,sendLv,bicon
  integer(kind=4)::i,index,count,sID,eID,iLv

  if(MinID(1,sendLv)>=MaxID(1,sendLv))return
  if(n_rcells_proc(sdest+1)==0)return
  if(sendLv==-1)return

  !set loop index
  if(sdest==0)then
     sID=1
  else
     sID=sum(n_rcells_proc(1:sdest))+1
  endif
  eID=sum(n_rcells_proc(1:sdest+1))

  count=1
  iLv  =-1
  do i=sID,eID
     index=iBMesh_arr(i)
     call pack_BMeshG(index,sendLv,count,icon,iLv+1,bicon)
  enddo

  send_onum(bicon)=count-1

end subroutine pack_BMeshGT


!new communication routine for J
recursive subroutine set_received_BMeshJ(index,recvLv,count,iLv,icon)
  use message_passing_interface
  use oct_set
  use param
  use init_mesh_size
  implicit none
  integer(kind=4),intent(in)::index,recvLv,iLv,icon
  integer(kind=4)::count
  type(oct),pointer::p0

  p0=>Mesh(index)
  if(p0%iFLG(1)>0)then
     if(iLv==recvLv)then
        p0%octCh1%F(7:18)=p0%octCh1%F(7:18)+sbuf_BFR(1+(count-1)*12:count*12,icon)
        count=count+1
        p0%octCh2%F(7:18)=p0%octCh2%F(7:18)+sbuf_BFR(1+(count-1)*12:count*12,icon)
        count=count+1
        p0%octCh3%F(7:18)=p0%octCh3%F(7:18)+sbuf_BFR(1+(count-1)*12:count*12,icon)
        count=count+1
        p0%octCh4%F(7:18)=p0%octCh4%F(7:18)+sbuf_BFR(1+(count-1)*12:count*12,icon)
        count=count+1
        p0%octCh5%F(7:18)=p0%octCh5%F(7:18)+sbuf_BFR(1+(count-1)*12:count*12,icon)
        count=count+1
        p0%octCh6%F(7:18)=p0%octCh6%F(7:18)+sbuf_BFR(1+(count-1)*12:count*12,icon)
        count=count+1
        p0%octCh7%F(7:18)=p0%octCh7%F(7:18)+sbuf_BFR(1+(count-1)*12:count*12,icon)
        count=count+1
        p0%octCh8%F(7:18)=p0%octCh8%F(7:18)+sbuf_BFR(1+(count-1)*12:count*12,icon)
        count=count+1
     else
        call set_received_BMeshJ(p0%octCh1%octN,recvLv,count,iLv+1,icon)
        call set_received_BMeshJ(p0%octCh2%octN,recvLv,count,iLv+1,icon)
        call set_received_BMeshJ(p0%octCh3%octN,recvLv,count,iLv+1,icon)
        call set_received_BMeshJ(p0%octCh4%octN,recvLv,count,iLv+1,icon)
        call set_received_BMeshJ(p0%octCh5%octN,recvLv,count,iLv+1,icon)
        call set_received_BMeshJ(p0%octCh6%octN,recvLv,count,iLv+1,icon)
        call set_received_BMeshJ(p0%octCh7%octN,recvLv,count,iLv+1,icon)
        call set_received_BMeshJ(p0%octCh8%octN,recvLv,count,iLv+1,icon)
     endif
  endif

end subroutine set_received_BMeshJ

subroutine set_received_BMeshJT(recvLv,recvnum,icon)
  use message_passing_interface
  use oct_set
  use param
  use init_mesh_size
  implicit none
  integer(kind=4),intent(in)::recvLv,recvnum,icon
  integer(kind=4)::count,index,iLv

  if(MinID(1,recvLv)>=MaxID(1,recvLv))return
  if(recvnum<=0)return  

  count=1
  iLv=-1

  do while(count<=recvnum)
     index=sbuf_BFI(count,icon)-MnDisp 
     call set_received_BMeshJ(index,recvLv,count,iLv+1,icon)
  end do

end subroutine set_received_BMeshJT

recursive subroutine set_received_GMeshEB(index,recvLv,count,icon,iLv,bicon)
  use param
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::count
  integer(kind=4),intent(in)::index,recvLv,iLv,bicon,icon
  type(oct),pointer::p0

  p0=>Mesh(index)
  if(p0%iFLG(1)>0)then
     if(iLv==recvLv)then
        p0%octCh1%F(1+(icon-1)*3:icon*3)=rbuf_GFR(1+(count-1)*3:count*3,bicon)
        count=count+1
        p0%octCh2%F(1+(icon-1)*3:icon*3)=rbuf_GFR(1+(count-1)*3:count*3,bicon)
        count=count+1
        p0%octCh3%F(1+(icon-1)*3:icon*3)=rbuf_GFR(1+(count-1)*3:count*3,bicon)
        count=count+1
        p0%octCh4%F(1+(icon-1)*3:icon*3)=rbuf_GFR(1+(count-1)*3:count*3,bicon)
        count=count+1
        p0%octCh5%F(1+(icon-1)*3:icon*3)=rbuf_GFR(1+(count-1)*3:count*3,bicon)
        count=count+1
        p0%octCh6%F(1+(icon-1)*3:icon*3)=rbuf_GFR(1+(count-1)*3:count*3,bicon)
        count=count+1
        p0%octCh7%F(1+(icon-1)*3:icon*3)=rbuf_GFR(1+(count-1)*3:count*3,bicon)
        count=count+1
        p0%octCh8%F(1+(icon-1)*3:icon*3)=rbuf_GFR(1+(count-1)*3:count*3,bicon)
        count=count+1
     else
        call set_received_GMeshEB(p0%octCh1%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshEB(p0%octCh2%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshEB(p0%octCh3%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshEB(p0%octCh4%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshEB(p0%octCh5%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshEB(p0%octCh6%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshEB(p0%octCh7%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshEB(p0%octCh8%octN,recvLv,count,icon,iLv+1,bicon)
     endif
  endif
end subroutine set_received_GMeshEB

subroutine set_received_GMeshEBT(recvLv,recvnum,icon,bicon)
  use message_passing_interface
  use oct_set
  use param
  use init_mesh_size
  implicit none
  integer(kind=4),intent(in)::recvLv,recvnum,icon,bicon !process num of sender
  integer(kind=4)::count,index,iLv

  if(MinID(3,recvLv)>=MaxID(3,recvLv))return
  if(recvnum<=0)return
  count=1
  iLv=-1

  do while(count<=recvnum)
     index=GMn2octN(rbuf_GFI(count,bicon))
    
     call set_received_GMeshEB(index,recvLv,count,icon,iLv+1,bicon)
  enddo
end subroutine set_received_GMeshEBT

recursive subroutine set_received_GMeshG(index,recvLv,count,icon,iLv,bicon)
  use param
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::count
  integer(kind=4),intent(in)::index,recvLv,iLv,bicon,icon
  type(oct),pointer::p0

  p0=>Mesh(index)
  if(p0%iFLG(1)>0)then
     if(iLv==recvLv)then
        p0%octCh1%G(1:18)=rbuf_G(1+(count-1)*18:count*18,bicon)
        count=count+1
        p0%octCh2%G(1:18)=rbuf_G(1+(count-1)*18:count*18,bicon)
        count=count+1
        p0%octCh3%G(1:18)=rbuf_G(1+(count-1)*18:count*18,bicon)
        count=count+1
        p0%octCh4%G(1:18)=rbuf_G(1+(count-1)*18:count*18,bicon)
        count=count+1
        p0%octCh5%G(1:18)=rbuf_G(1+(count-1)*18:count*18,bicon)
        count=count+1
        p0%octCh6%G(1:18)=rbuf_G(1+(count-1)*18:count*18,bicon)
        count=count+1
        p0%octCh7%G(1:18)=rbuf_G(1+(count-1)*18:count*18,bicon)
        count=count+1
        p0%octCh8%G(1:18)=rbuf_G(1+(count-1)*18:count*18,bicon)
        count=count+1
     else
        call set_received_GMeshG(p0%octCh1%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshG(p0%octCh2%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshG(p0%octCh3%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshG(p0%octCh4%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshG(p0%octCh5%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshG(p0%octCh6%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshG(p0%octCh7%octN,recvLv,count,icon,iLv+1,bicon)
        call set_received_GMeshG(p0%octCh8%octN,recvLv,count,icon,iLv+1,bicon)
     endif
  endif
end subroutine set_received_GMeshG

subroutine set_received_GMeshGT(recvLv,recvnum,icon,bicon)
  use message_passing_interface
  use oct_set
  use param
  use init_mesh_size
  implicit none
  integer(kind=4),intent(in)::recvLv,recvnum,icon,bicon !process num of sender
  integer(kind=4)::count,index,iLv

  if(MinID(3,recvLv)>=MaxID(3,recvLv))return
  if(recvnum<=0)return
  count=1
  iLv=-1

  do while(count<=recvnum)
     index=GMn2octN(rbuf_GFI(count,bicon))
    
     call set_received_GMeshG(index,recvLv,count,icon,iLv+1,bicon)
  enddo
end subroutine set_received_GMeshGT

subroutine refresh_Fields(type,sendLv)
  use message_passing_interface
  use param
  use oct_set
  !
  use init_mesh_size
  !
  implicit none
  integer(kind=4)::i,j,rdest(2),sdest(2),icon
  integer(kind=4)::req(4,2),istat(MPI_STATUS_SIZE,4,2)
  integer(kind=4),intent(in)::sendLv
  character(LEN=1),intent(in)::type
  integer(kind=4)::loopnum,stream1,stream2,streamtemp

  if(sendLv==-1)return

  if(type=="E")then
     icon=1
  else if(type=="B")then
     icon=2
  else if(type=="J")then
     icon=3
  else if(type=="G")then
     icon=4
  else
     print *,"type error type=",type,rank
     stop
  endif

  !first operation
  loopnum=rrlistnum-1
  stream1=1
  stream2=2
  j=0



  if(icon==3)then !J ==========================================================

     !first operation
     loopnum=rrlistnum-1
     stream1=1
     stream2=2
     j=0
     !-----stream1 op=1
     j=j+1
     sdest(stream1)=send_glist(j)
     rdest(stream1)=recv_rlist(j)
     send_onum(stream1)=0
     recv_onum(stream1)=0
     if(n_rcells_proc(rdest(stream1)+1)>0)call mpi_irecv(recv_onum(stream1),1,MPI_INTEGER,rdest(stream1),1,&
           &MPI_COMM_WORLD,req(1,stream1),ierr)
     call pack_GMeshJT (sdest(stream1),sendLv,stream1)
     if(n_gcells_proc(sdest(stream1)+1)>0)call mpi_isend(send_onum(stream1),1,MPI_INTEGER,sdest(stream1),1,&
           &MPI_COMM_WORLD,req(2,stream1),ierr)
     if(n_rcells_proc(rdest(stream1)+1)>0)call mpi_wait     (req(1,stream1),          istat(1,1,stream1),ierr)
     if(n_gcells_proc(sdest(stream1)+1)>0)call mpi_wait     (req(2,stream1),          istat(1,2,stream1),ierr)
     !-----stream1 op=2
     if(recv_onum(stream1)>0)then
        call mpi_irecv(sBuf_BFR(:,stream1),recv_onum(stream1)*12,MPI_DOUBLE_PRECISION,rdest(stream1),2,&
           &MPI_COMM_WORLD,req(1,stream1),ierr)
        call mpi_irecv(sBuf_BFI(:,stream1),recv_onum(stream1)   ,MPI_INTEGER         ,rdest(stream1),3,&
           &MPI_COMM_WORLD,req(2,stream1),ierr)
     endif
     if(send_onum(stream1)>0)then
        call mpi_isend(rbuf_GFR(:,stream1),send_onum(stream1)*12,MPI_DOUBLE_PRECISION,sdest(stream1),2,&
           &MPI_COMM_WORLD,req(3,stream1),ierr)
        call mpi_isend(rbuf_GFI(:,stream1),send_onum(stream1)   ,MPI_INTEGER         ,sdest(stream1),3,&
           &MPI_COMM_WORLD,req(4,stream1),ierr)
     endif
     
     do i=1,loopnum
        streamtemp=stream1
        stream1=stream2
        stream2=streamtemp
        !-------stream1 op=1
        j=j+1
        sdest(stream1)=send_glist(j)
        rdest(stream1)=recv_rlist(j)
        send_onum(stream1)=0
        recv_onum(stream1)=0
        if(n_rcells_proc(rdest(stream1)+1)>0)call mpi_irecv(recv_onum(stream1),1,MPI_INTEGER,rdest(stream1),1,&
           &MPI_COMM_WORLD,req(1,stream1),ierr)!prepare to receive recv_onum
        call pack_GMeshJT (sdest(stream1),sendLv,stream1)
        if(n_gcells_proc(sdest(stream1)+1)>0)call mpi_isend(send_onum(stream1),1,MPI_INTEGER,sdest(stream1),1,&
           &MPI_COMM_WORLD,req(2,stream1),ierr)!trade send_onum
        if(n_rcells_proc(rdest(stream1)+1)>0)call mpi_wait     (req(1,stream1),          istat(1,1,stream1),ierr)
        if(n_gcells_proc(sdest(stream1)+1)>0)call mpi_wait     (req(2,stream1),          istat(1,2,stream1),ierr)
        !-------stream1 op=2
        if(recv_onum(stream1)>0)then
           call mpi_irecv(sBuf_BFR(:,stream1),recv_onum(stream1)*18,MPI_DOUBLE_PRECISION,rdest(stream1),2,&
           &MPI_COMM_WORLD,req(1,stream1),ierr)
           call mpi_irecv(sBuf_BFI(:,stream1),recv_onum(stream1)   ,MPI_INTEGER         ,rdest(stream1),3,&
           &MPI_COMM_WORLD,req(2,stream1),ierr)
        endif
        if(send_onum(stream1)>0)then
           call mpi_isend(rbuf_GFR(:,stream1),send_onum(stream1)*18,MPI_DOUBLE_PRECISION,sdest(stream1),2,&
           &MPI_COMM_WORLD,req(3,stream1),ierr)
           call mpi_isend(rbuf_GFI(:,stream1),send_onum(stream1)   ,MPI_INTEGER         ,sdest(stream1),3,&
           &MPI_COMM_WORLD,req(4,stream1),ierr)
        endif
        !-------stream2 op=3
        if(recv_onum(stream2)>0)call mpi_waitall(2,req(1,stream2),istat(1,1,stream2),ierr)
        if(send_onum(stream2)>0)call mpi_waitall(2,req(3,stream2),istat(1,2,stream2),ierr)
        call set_received_BMeshJT (sendLv,recv_onum(stream2),stream2)
     enddo
     !---------stream1 op=3
     if(recv_onum(stream1)>0)call mpi_waitall(2,req(1,stream1),istat(1,1,stream1),ierr)
     if(send_onum(stream1)>0)call mpi_waitall(2,req(3,stream1),istat(1,2,stream1),ierr)
     call set_received_BMeshJT (sendLv,recv_onum(stream1),stream1)

  else if (icon==4) then
     !first operation
     loopnum=rrlistnum-1
     stream1=1
     stream2=2
     j=0
     !-----stream1 op=1
     j=j+1
     sdest(stream1)=send_rlist(j)
     rdest(stream1)=recv_glist(j)
     send_onum(stream1)=0
     recv_onum(stream1)=0
     if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_irecv(recv_onum(stream1),1,MPI_INTEGER,rdest(stream1),1,&
     &MPI_COMM_WORLD,req(1,stream1),ierr)
     call pack_BMeshGT(sdest(stream1),sendLv,icon,stream1)
     if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_isend(send_onum(stream1),1,MPI_INTEGER,sdest(stream1),1,&
     &MPI_COMM_WORLD,req(2,stream1),ierr)
     if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_wait(     req(1,stream1),          istat(1,1,stream1),ierr)
     if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_wait(     req(2,stream1),          istat(1,2,stream1),ierr)
     !-----stream1 op=2
     if(recv_onum(stream1)>0)then
        call mpi_irecv(rBuf_G(:,stream1),recv_onum(stream1)*18,MPI_DOUBLE_PRECISION,rdest(stream1),2,&
        &MPI_COMM_WORLD,req(1,stream1),ierr)
     endif
     if(send_onum(stream1)>0)then
        call mpi_isend(sbuf_G(:,stream1),send_onum(stream1)*18,MPI_DOUBLE_PRECISION,sdest(stream1),2,&
        &MPI_COMM_WORLD,req(3,stream1),ierr)
     endif

     do i=1,loopnum
        streamtemp=stream1
        stream1=stream2
        stream2=streamtemp
        !-----stream1 op=1
        j=j+1
        sdest(stream1)=send_rlist(j)
        rdest(stream1)=recv_glist(j)
        send_onum(stream1)=0
        recv_onum(stream1)=0
        if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_irecv(recv_onum(stream1),1,MPI_INTEGER,rdest(stream1),1,&
        &MPI_COMM_WORLD,req(1,stream1),ierr)
        call pack_BMeshGT(sdest(stream1),sendLv,icon,stream1)
        if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_isend(send_onum(stream1),1,MPI_INTEGER,sdest(stream1),1,&
        &MPI_COMM_WORLD,req(2,stream1),ierr)
        if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_wait(     req(1,stream1),          istat(1,1,stream1),ierr)
        if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_wait(     req(2,stream1),          istat(1,2,stream1),ierr)
        !-------stream1 op=2
        if(recv_onum(stream1)>0)then
           call mpi_irecv(rBuf_G(:,stream1),recv_onum(stream1)*18,MPI_DOUBLE_PRECISION,rdest(stream1),2,&
           &MPI_COMM_WORLD,req(1,stream1),ierr)
        endif
        if(send_onum(stream1)>0)then
           call mpi_isend(sbuf_G(:,stream1),send_onum(stream1)*18,MPI_DOUBLE_PRECISION,sdest(stream1),2,&
           &MPI_COMM_WORLD,req(3,stream1),ierr)
        endif
        !-------stream2 op=3
        if(recv_onum(stream2)>0)call mpi_waitall(2,req(1,stream2),istat(1,1,stream2),ierr)
        if(send_onum(stream2)>0)call mpi_waitall(2,req(3,stream2),istat(1,1,stream2),ierr)
        call set_received_GMeshGT(sendLv,recv_onum(stream2),icon,stream2)
     end do
     !-----stream1 op=3
     if(recv_onum(stream1)>0)call mpi_waitall(2,req(1,stream1),istat(1,1,stream1),ierr)
     if(send_onum(stream1)>0)call mpi_waitall(2,req(3,stream1),istat(1,1,stream1),ierr)
     call set_received_GMeshGT(sendLv,recv_onum(stream1),icon,stream1)

  else !E,B =================================================================

     !if(rank==0)print *,"refresh E start"
     !first operation
     loopnum=rrlistnum-1
     stream1=1
     stream2=2
     j=0
     !-----stream1 op=1
     j=j+1
     sdest(stream1)=send_rlist(j)
     rdest(stream1)=recv_glist(j)
     send_onum(stream1)=0
     recv_onum(stream1)=0
     if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_irecv(recv_onum(stream1),1,MPI_INTEGER,rdest(stream1),1,&
     &MPI_COMM_WORLD,req(1,stream1),ierr)
     call pack_BMeshEBT(sdest(stream1),sendLv,icon,stream1)
     if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_isend(send_onum(stream1),1,MPI_INTEGER,sdest(stream1),1,&
     &MPI_COMM_WORLD,req(2,stream1),ierr)
     if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_wait(     req(1,stream1),          istat(1,1,stream1),ierr)
     if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_wait(     req(2,stream1),          istat(1,2,stream1),ierr)
     !-----stream1 op=2
     if(recv_onum(stream1)>0)then
        call mpi_irecv(rBuf_GFR(:,stream1),recv_onum(stream1)*6,MPI_DOUBLE_PRECISION,rdest(stream1),2,&
        &MPI_COMM_WORLD,req(1,stream1),ierr)
        call mpi_irecv(rBuf_GFI(:,stream1),recv_onum(stream1)  ,MPI_INTEGER         ,rdest(stream1),3,&
        &MPI_COMM_WORLD,req(2,stream1),ierr)
     endif
     if(send_onum(stream1)>0)then
        call mpi_isend(sbuf_BFR(:,stream1),send_onum(stream1)*6,MPI_DOUBLE_PRECISION,sdest(stream1),2,&
        &MPI_COMM_WORLD,req(3,stream1),ierr)
        call mpi_isend(sbuf_BFI(:,stream1),send_onum(stream1)  ,MPI_INTEGER         ,sdest(stream1),3,&
        &MPI_COMM_WORLD,req(4,stream1),ierr)
     endif

     do i=1,loopnum
        streamtemp=stream1
        stream1=stream2
        stream2=streamtemp
        !-----stream1 op=1
        j=j+1
        sdest(stream1)=send_rlist(j)
        rdest(stream1)=recv_glist(j)
        send_onum(stream1)=0
        recv_onum(stream1)=0
        if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_irecv(recv_onum(stream1),1,MPI_INTEGER,rdest(stream1),1,&
        &MPI_COMM_WORLD,req(1,stream1),ierr)
        call pack_BMeshEBT(sdest(stream1),sendLv,icon,stream1)
        if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_isend(send_onum(stream1),1,MPI_INTEGER,sdest(stream1),1,&
        &MPI_COMM_WORLD,req(2,stream1),ierr)
        if(n_gcells_proc(rdest(stream1)+1)>0)call mpi_wait(     req(1,stream1),          istat(1,1,stream1),ierr)
        if(n_rcells_proc(sdest(stream1)+1)>0)call mpi_wait(     req(2,stream1),          istat(1,2,stream1),ierr)
        !-------stream1 op=2
        if(recv_onum(stream1)>0)then
           call mpi_irecv(rBuf_GFR(:,stream1),recv_onum(stream1)*6,MPI_DOUBLE_PRECISION,rdest(stream1),2,&
           &MPI_COMM_WORLD,req(1,stream1),ierr)
           call mpi_irecv(rBuf_GFI(:,stream1),recv_onum(stream1)  ,MPI_INTEGER         ,rdest(stream1),3,&
           &MPI_COMM_WORLD,req(2,stream1),ierr)
        endif
        if(send_onum(stream1)>0)then
           call mpi_isend(sbuf_BFR(:,stream1),send_onum(stream1)*6,MPI_DOUBLE_PRECISION,sdest(stream1),2,&
           &MPI_COMM_WORLD,req(3,stream1),ierr)
           call mpi_isend(sbuf_BFI(:,stream1),send_onum(stream1)  ,MPI_INTEGER         ,sdest(stream1),3,&
           &MPI_COMM_WORLD,req(4,stream1),ierr)
        endif
        !-------stream2 op=3
        if(recv_onum(stream2)>0)call mpi_waitall(2,req(1,stream2),istat(1,1,stream2),ierr)
        if(send_onum(stream2)>0)call mpi_waitall(2,req(3,stream2),istat(1,1,stream2),ierr)
        call set_received_GMeshEBT(sendLv,recv_onum(stream2),icon,stream2)
     end do
     !-----stream1 op=3
     if(recv_onum(stream1)>0)call mpi_waitall(2,req(1,stream1),istat(1,1,stream1),ierr)
     if(send_onum(stream1)>0)call mpi_waitall(2,req(3,stream1),istat(1,1,stream1),ierr)
     call set_received_GMeshEBT(sendLv,recv_onum(stream1),icon,stream1)

  endif

end subroutine refresh_Fields

!--------------------Interprocess comm of particles for time evolution------------- 
!
!----------------------------------------------------------------------------------

recursive subroutine set_received_BMeshP(index,countOct,countPtcl,recvLv,iLv,icon)
  use param
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::countOct,countPtcl,i,indexP
  integer(kind=4)::octP1,octP2,octP3,octP4,octP5,octP6,octP7,octP8
  integer(kind=4)::octN1,octN2,octN3,octN4,octN5,octN6,octN7,octN8
  integer(kind=4),intent(in)::index,recvLv,iLv,icon
  type(oct),pointer::p0
  type(prtcl), pointer:: pp,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8

  p0=>Mesh(index)

  if(p0%iFLG(1)>0)then
     if(iLv==recvLv)then
        octN1=p0%octCh1%octN
        octN2=p0%octCh2%octN
        octN3=p0%octCh3%octN
        octN4=p0%octCh4%octN
        octN5=p0%octCh5%octN
        octN6=p0%octCh6%octN
        octN7=p0%octCh7%octN
        octN8=p0%octCh8%octN

        octP1=rbuf_BPO(countOct*2-1,icon)
        countOct=countOct+1
        octP2=rbuf_BPO(countOct*2-1,icon)
        countOct=countOct+1
        octP3=rbuf_BPO(countOct*2-1,icon)
        countOct=countOct+1
        octP4=rbuf_BPO(countOct*2-1,icon)
        countOct=countOct+1
        octP5=rbuf_BPO(countOct*2-1,icon)
        countOct=countOct+1
        octP6=rbuf_BPO(countOct*2-1,icon)
        countOct=countOct+1
        octP7=rbuf_BPO(countOct*2-1,icon)
        countOct=countOct+1
        octP8=rbuf_BPO(countOct*2-1,icon)
        countOct=countOct+1

        p0%octCh1%octP=p0%octCh1%octP+octP1
        p0%octCh2%octP=p0%octCh2%octP+octP2
        p0%octCh3%octP=p0%octCh3%octP+octP3
        p0%octCh4%octP=p0%octCh4%octP+octP4
        p0%octCh5%octP=p0%octCh5%octP+octP5
        p0%octCh6%octP=p0%octCh6%octP+octP6
        p0%octCh7%octP=p0%octCh7%octP+octP7
        p0%octCh8%octP=p0%octCh8%octP+octP8      

        pp1=>p0%octCh1%ptcl
        pp2=>p0%octCh2%ptcl
        pp3=>p0%octCh3%ptcl
        pp4=>p0%octCh4%ptcl
        pp5=>p0%octCh5%ptcl
        pp6=>p0%octCh6%ptcl
        pp7=>p0%octCh7%ptcl
        pp8=>p0%octCh8%ptcl

        indexP=MaxIP(iLv)

        do i=1,octP1
           indexP=indexP+1
           pp=>pp1%prtnxt
           Pesh(indexP,iLv)=prtcl(rBuf_BPR(1+(countPtcl-1)*9:countPtcl*9,icon),rBuf_BPI(countPtcl,icon),octN1,indexP,pp,0)
           pp1%prtnxt => Pesh(indexP,iLv)
           countPtcl=countPtcl+1
        enddo
        do i=1,octP2
           indexP=indexP+1
           pp=>pp2%prtnxt
           Pesh(indexP,iLv)=prtcl(rBuf_BPR(1+(countPtcl-1)*9:countPtcl*9,icon),rBuf_BPI(countPtcl,icon),octN2,indexP,pp,0)
           pp2%prtnxt => Pesh(indexP,iLv)
           countPtcl=countPtcl+1
        enddo
        do i=1,octP3
           indexP=indexP+1
           pp=>pp3%prtnxt
           Pesh(indexP,iLv)=prtcl(rBuf_BPR(1+(countPtcl-1)*9:countPtcl*9,icon),rBuf_BPI(countPtcl,icon),octN3,indexP,pp,0)
           pp3%prtnxt => Pesh(indexP,iLv)
           countPtcl=countPtcl+1
        enddo
        do i=1,octP4
           indexP=indexP+1
           pp=>pp4%prtnxt
           Pesh(indexP,iLv)=prtcl(rBuf_BPR(1+(countPtcl-1)*9:countPtcl*9,icon),rBuf_BPI(countPtcl,icon),octN4,indexP,pp,0)
           pp4%prtnxt => Pesh(indexP,iLv)
           countPtcl=countPtcl+1
        enddo
        do i=1,octP5
           indexP=indexP+1
           pp=>pp5%prtnxt
           Pesh(indexP,iLv)=prtcl(rBuf_BPR(1+(countPtcl-1)*9:countPtcl*9,icon),rBuf_BPI(countPtcl,icon),octN5,indexP,pp,0)
           pp5%prtnxt => Pesh(indexP,iLv)
           countPtcl=countPtcl+1
        enddo
        do i=1,octP6
           indexP=indexP+1
           pp=>pp6%prtnxt
           Pesh(indexP,iLv)=prtcl(rBuf_BPR(1+(countPtcl-1)*9:countPtcl*9,icon),rBuf_BPI(countPtcl,icon),octN6,indexP,pp,0)
           pp6%prtnxt => Pesh(indexP,iLv)
           countPtcl=countPtcl+1
        enddo
        do i=1,octP7
           indexP=indexP+1
           pp=>pp7%prtnxt
           Pesh(indexP,iLv)=prtcl(rBuf_BPR(1+(countPtcl-1)*9:countPtcl*9,icon),rBuf_BPI(countPtcl,icon),octN7,indexP,pp,0)
           pp7%prtnxt => Pesh(indexP,iLv)
           countPtcl=countPtcl+1
        enddo
        do i=1,octP8
           indexP=indexP+1
           pp=>pp8%prtnxt
           Pesh(indexP,iLv)=prtcl(rBuf_BPR(1+(countPtcl-1)*9:countPtcl*9,icon),rBuf_BPI(countPtcl,icon),octN8,indexP,pp,0)
           pp8%prtnxt => Pesh(indexP,iLv)
           countPtcl=countPtcl+1
        enddo

        MaxIP(iLv)=indexP


     else
        call set_received_BMeshP(p0%octCh1%octN,countOct,countPtcl,recvLv,iLv+1,icon)
        call set_received_BMeshP(p0%octCh2%octN,countOct,countPtcl,recvLv,iLv+1,icon)
        call set_received_BMeshP(p0%octCh3%octN,countOct,countPtcl,recvLv,iLv+1,icon)
        call set_received_BMeshP(p0%octCh4%octN,countOct,countPtcl,recvLv,iLv+1,icon)
        call set_received_BMeshP(p0%octCh5%octN,countOct,countPtcl,recvLv,iLv+1,icon)
        call set_received_BMeshP(p0%octCh6%octN,countOct,countPtcl,recvLv,iLv+1,icon)
        call set_received_BMeshP(p0%octCh7%octN,countOct,countPtcl,recvLv,iLv+1,icon)
        call set_received_BMeshP(p0%octCh8%octN,countOct,countPtcl,recvLv,iLv+1,icon)
     endif
  endif

end subroutine set_received_BMeshP

subroutine set_received_BMeshPT(recvLv,recvonum,recvpnum,icon)
  use message_passing_interface
  use oct_set
  use param
  use init_mesh_size
  implicit none
  integer(kind=4),intent(in)::recvLv,recvonum,recvpnum,icon
  integer(kind=4)::countOct,countPtcl,index,iLv

  if(MinID(1,recvLv)>=MaxID(1,recvLv))return
  if(recvpnum<=0)return
  if(recvLv==-1)return

  countOct=1
  countPtcl=1
  iLv=-1

  do while(countOct<=recvonum)
     index=rbuf_BPO(countOct*2,icon)-MnDisp
     call set_received_BMeshP(index,countOct,countPtcl,recvLv,iLv+1,icon)
  enddo
end subroutine set_received_BMeshPT

subroutine refresh_Particles(sendLv)
  use message_passing_interface
  use param
  use oct_set
  implicit none
  integer(kind=4),intent(in)::sendLv
  integer(kind=4)::i,j,rdest(2),sdest(2)
  integer(kind=4)::req(6,2),istat(MPI_STATUS_SIZE,6,2)
  integer(kind=4)::stream1,stream2,streamtemp,loopnum

     loopnum=rrlistnum-1
     stream1=1
     stream2=2
     j=0
     !-----stream=1 op=1 exchange onum and pnum
     j=j+1
     sdest(stream1)=send_glist(j)
     rdest(stream1)=recv_rlist(j)
     send_onum(stream1)=0
     recv_onum(stream1)=0
     send_pnum(stream1)=0
     recv_pnum(stream1)=0

     if(n_rcells_proc(rdest(stream1)+1)>0)then
        call mpi_irecv(recv_pnum(stream1),1,MPI_INTEGER,rdest(stream1),1,MPI_COMM_WORLD,req(1,stream1),ierr)
        call mpi_irecv(recv_onum(stream1),1,MPI_INTEGER,rdest(stream1),2,MPI_COMM_WORLD,req(2,stream1),ierr)
     endif
     call pack_GMesh_and_PnumT(sdest(stream1),sendLv,stream1)
     if(n_gcells_proc(sdest(stream1)+1)>0)then
        call mpi_isend(send_pnum(stream1),1,MPI_INTEGER,sdest(stream1),1,MPI_COMM_WORLD,req(3,stream1),ierr)
        call mpi_isend(send_onum(stream1),1,MPI_INTEGER,sdest(stream1),2,MPI_COMM_WORLD,req(4,stream1),ierr)
     endif

     if(n_rcells_proc(rdest(stream1)+1)>0)call mpi_waitall(2,req(1,stream1),istat(1,1,stream1),ierr)
     if(n_gcells_proc(sdest(stream1)+1)>0)call mpi_waitall(2,req(3,stream1),istat(1,1,stream1),ierr)
     !-----stream=1 op=2 prepare to exchange
     if(recv_pnum(stream1)>0)then
        call mpi_irecv(rBuf_BPO(:,stream1),recv_onum(stream1)*2,MPI_INTEGER         ,rdest(stream1),3,MPI_COMM_WORLD,req(1,stream1),ierr)
        call mpi_irecv(rBuf_BPI(:,stream1),recv_pnum(stream1)*1,MPI_INTEGER         ,rdest(stream1),4,MPI_COMM_WORLD,req(2,stream1),ierr)
        call mpi_irecv(rBuf_BPR(:,stream1),recv_pnum(stream1)*9,MPI_DOUBLE_PRECISION,rdest(stream1),5,MPI_COMM_WORLD,req(3,stream1),ierr)
     endif
     if(send_pnum(stream1)>0)then
        call mpi_isend(sbuf_GPO(:,stream1),send_onum(stream1)*2,MPI_INTEGER         ,sdest(stream1),3,MPI_COMM_WORLD,req(4,stream1),ierr)
        call mpi_isend(sbuf_GPI(:,stream1),send_pnum(stream1)*1,MPI_INTEGER         ,sdest(stream1),4,MPI_COMM_WORLD,req(5,stream1),ierr)
        call mpi_isend(sbuf_GPR(:,stream1),send_pnum(stream1)*9,MPI_DOUBLE_PRECISION,sdest(stream1),5,MPI_COMM_WORLD,req(6,stream1),ierr)
     endif

     do i=1,loopnum
        streamtemp=stream1
        stream1=stream2
        stream2=streamtemp
        !-------stream1 op=1
        j=j+1
        sdest(stream1)=send_glist(j)
        rdest(stream1)=recv_rlist(j)
        send_onum(stream1)=0
        recv_onum(stream1)=0
        send_pnum(stream1)=0
        recv_pnum(stream1)=0
        if(n_rcells_proc(rdest(stream1)+1)>0)then
           call mpi_irecv(recv_pnum(stream1),1,MPI_INTEGER,rdest(stream1),1,MPI_COMM_WORLD,req(1,stream1),ierr)
           call mpi_irecv(recv_onum(stream1),1,MPI_INTEGER,rdest(stream1),2,MPI_COMM_WORLD,req(2,stream1),ierr)
        endif
        call pack_GMesh_and_PnumT(sdest(stream1),sendLv,stream1)
        if(n_gcells_proc(sdest(stream1)+1)>0)then
           call mpi_isend(send_pnum(stream1),1,MPI_INTEGER,sdest(stream1),1,MPI_COMM_WORLD,req(3,stream1),ierr)
           call mpi_isend(send_onum(stream1),1,MPI_INTEGER,sdest(stream1),2,MPI_COMM_WORLD,req(4,stream1),ierr)
        endif
        if(n_rcells_proc(rdest(stream1)+1)>0)call mpi_waitall(2,req(1,stream1),istat(1,1,stream1),ierr)
        if(n_gcells_proc(sdest(stream1)+1)>0)call mpi_waitall(2,req(3,stream1),istat(1,1,stream1),ierr)
        !-------stream1 op=2
        if(recv_pnum(stream1)>0)then
           call mpi_irecv(rBuf_BPO(:,stream1),recv_onum(stream1)*2,MPI_INTEGER         ,rdest(stream1),3,MPI_COMM_WORLD,req(1,stream1),ierr)
           call mpi_irecv(rBuf_BPI(:,stream1),recv_pnum(stream1)*1,MPI_INTEGER         ,rdest(stream1),4,MPI_COMM_WORLD,req(2,stream1),ierr)
           call mpi_irecv(rBuf_BPR(:,stream1),recv_pnum(stream1)*9,MPI_DOUBLE_PRECISION,rdest(stream1),5,MPI_COMM_WORLD,req(3,stream1),ierr)
        endif
        if(send_pnum(stream1)>0)then
           call mpi_isend(sbuf_GPO(:,stream1),send_onum(stream1)*2,MPI_INTEGER         ,sdest(stream1),3,MPI_COMM_WORLD,req(4,stream1),ierr)
           call mpi_isend(sbuf_GPI(:,stream1),send_pnum(stream1)*1,MPI_INTEGER         ,sdest(stream1),4,MPI_COMM_WORLD,req(5,stream1),ierr)
           call mpi_isend(sbuf_GPR(:,stream1),send_pnum(stream1)*9,MPI_DOUBLE_PRECISION,sdest(stream1),5,MPI_COMM_WORLD,req(6,stream1),ierr)
        endif
        !------stream2 op=3
        if(recv_pnum(stream2)>0)call mpi_waitall(3,req(1,stream2),istat(1,1,stream2),ierr)
        if(send_pnum(stream2)>0)call mpi_waitall(3,req(4,stream2),istat(1,1,stream2),ierr)
        call set_received_BMeshPT(sendLv,recv_onum(stream2),recv_pnum(stream2),stream2)
     end do

     !-----stream1 op=3
     if(recv_pnum(stream1)>0)call mpi_waitall(3,req(1,stream1),istat(1,1,stream1),ierr)
     if(send_pnum(stream1)>0)call mpi_waitall(3,req(4,stream1),istat(1,1,stream1),ierr)
     call set_received_BMeshPT(sendLv,recv_onum(stream1),recv_pnum(stream1),stream1)

  !if(debugMode>=1)print *,"refresh_Particle completed",rank
end subroutine refresh_Particles

!collect octData+pnum+onum
subroutine pack_GMesh_and_PnumT(sdest,sendLv,icon)
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::i,index,sID,eID,iLv
  integer(kind=4)::countOct,countPtcl
  integer(kind=4),intent(in)::sdest,sendLv,icon

  if(MinID(3,sendLv)>=MaxID(3,sendLv))return
  if(n_gcells_proc(sdest+1)==0)return

  if(sdest==0)then
     sID=1
  else
     sID=sum(n_gcells_proc(1:sdest))+1
  endif
  eID=sum(n_gcells_proc(1:sdest+1))
 
  countOct=1
  countPtcl=1
  iLv=-1

  do i=sID,eID
     index=iGMesh_arr(i)  
     call pack_GMesh_and_Pnum(index,countOct,countPtcl,sendLv,iLv+1,icon)  
  enddo

  send_onum(icon)=countOct-1
  send_pnum(icon)=countPtcl-1

  return
end subroutine pack_GMesh_and_PnumT

recursive subroutine pack_GMesh_and_Pnum(index,countOct,countPtcl,sendLv,iLv,icon)
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::i,index,iLv
  integer(kind=4)::countOct,countPtcl
  integer(kind=4),intent(in)::sendLv,icon
  integer(kind=4)::octP1,octP2,octP3,octP4,octP5,octP6,octP7,octP8
  type(oct),pointer::p0 
  type(prtcl),pointer::pp

  p0=>Mesh(index)
  if(p0%iFLG(1)>0)then
     if(iLv==sendLv)then
        octP1 = p0%octCh1%octP
        octP2 = p0%octCh2%octP
        octP3 = p0%octCh3%octP
        octP4 = p0%octCh4%octP
        octP5 = p0%octCh5%octP
        octP6 = p0%octCh6%octP
        octP7 = p0%octCh7%octP
        octP8 = p0%octCh8%octP

        !pack MrtN and octP
        sbuf_GPO(countOct*2  ,icon)=p0%octCh1%MrtN
        sbuf_GPO(countOct*2-1,icon)=octP1
        countOct=countOct+1
        sbuf_GPO(countOct*2  ,icon)=p0%octCh2%MrtN
        sbuf_GPO(countOct*2-1,icon)=octP2
        countOct=countOct+1
        sbuf_GPO(countOct*2  ,icon)=p0%octCh3%MrtN
        sbuf_GPO(countOct*2-1,icon)=octP3
        countOct=countOct+1
        sbuf_GPO(countOct*2  ,icon)=p0%octCh4%MrtN
        sbuf_GPO(countOct*2-1,icon)=octP4
        countOct=countOct+1
        sbuf_GPO(countOct*2  ,icon)=p0%octCh5%MrtN
        sbuf_GPO(countOct*2-1,icon)=octP5
        countOct=countOct+1
        sbuf_GPO(countOct*2  ,icon)=p0%octCh6%MrtN
        sbuf_GPO(countOct*2-1,icon)=octP6
        countOct=countOct+1
        sbuf_GPO(countOct*2  ,icon)=p0%octCh7%MrtN
        sbuf_GPO(countOct*2-1,icon)=octP7
        countOct=countOct+1
        sbuf_GPO(countOct*2  ,icon)=p0%octCh8%MrtN
        sbuf_GPO(countOct*2-1,icon)=octP8
        countOct=countOct+1

        !pack particle
        if(octP1>0)then
           pp=>p0%octCh1%ptcl
           do i=1,octP1
              pp=>pp%prtnxt
              
              sbuf_GPR(1+(countPtcl-1)*9:countPtcl*9,icon)=pp%R(1:9)
              sbuf_GPI(countPtcl                    ,icon)=pp%isort
              pp%isort=0
              countPtcl=countPtcl+1
           end do
        endif
        if(octP2>0)then
           pp=>p0%octCh2%ptcl
           do i=1,octP2
              pp=>pp%prtnxt
              
              sbuf_GPR(1+(countPtcl-1)*9:countPtcl*9,icon)=pp%R(1:9)
              sbuf_GPI(countPtcl                    ,icon)=pp%isort
              pp%isort=0
              countPtcl=countPtcl+1
           end do
        endif
        if(octP3>0)then
           pp=>p0%octCh3%ptcl
           do i=1,octP3
              pp=>pp%prtnxt
              
              sbuf_GPR(1+(countPtcl-1)*9:countPtcl*9,icon)=pp%R(1:9)
              sbuf_GPI(countPtcl                    ,icon)=pp%isort
              pp%isort=0
              countPtcl=countPtcl+1
           end do
        endif
        if(octP4>0)then
           pp=>p0%octCh4%ptcl
           do i=1,octP4
              pp=>pp%prtnxt
              
              sbuf_GPR(1+(countPtcl-1)*9:countPtcl*9,icon)=pp%R(1:9)
              sbuf_GPI(countPtcl                    ,icon)=pp%isort
              pp%isort=0
              countPtcl=countPtcl+1
           end do
        endif
        if(octP5>0)then
           pp=>p0%octCh5%ptcl
           do i=1,octP5
              pp=>pp%prtnxt
              
              sbuf_GPR(1+(countPtcl-1)*9:countPtcl*9,icon)=pp%R(1:9)
              sbuf_GPI(countPtcl                    ,icon)=pp%isort
              pp%isort=0
              countPtcl=countPtcl+1
           end do
        endif
        if(octP6>0)then
           pp=>p0%octCh6%ptcl
           do i=1,octP6
              pp=>pp%prtnxt
              
              sbuf_GPR(1+(countPtcl-1)*9:countPtcl*9,icon)=pp%R(1:9)
              sbuf_GPI(countPtcl                    ,icon)=pp%isort
              pp%isort=0
              countPtcl=countPtcl+1
           end do
        endif
        if(octP7>0)then
           pp=>p0%octCh7%ptcl
           do i=1,octP7
              pp=>pp%prtnxt
              
              sbuf_GPR(1+(countPtcl-1)*9:countPtcl*9,icon)=pp%R(1:9)
              sbuf_GPI(countPtcl                    ,icon)=pp%isort
              pp%isort=0
              countPtcl=countPtcl+1
           end do
        endif
        if(octP8>0)then
           pp=>p0%octCh8%ptcl
           do i=1,octP8
              pp=>pp%prtnxt
              
              sbuf_GPR(1+(countPtcl-1)*9:countPtcl*9,icon)=pp%R(1:9)
              sbuf_GPI(countPtcl                    ,icon)=pp%isort
              pp%isort=0
              countPtcl=countPtcl+1
           end do
        endif

        p0%octCh1%octP=0
        p0%octCh2%octP=0
        p0%octCh3%octP=0
        p0%octCh4%octP=0
        p0%octCh5%octP=0
        p0%octCh6%octP=0
        p0%octCh7%octP=0
        p0%octCh8%octP=0
     else
        call pack_GMesh_and_Pnum(p0%octCh1%octN,countOct,countPtcl,sendLv,iLv+1,icon)
        call pack_GMesh_and_Pnum(p0%octCh2%octN,countOct,countPtcl,sendLv,iLv+1,icon)
        call pack_GMesh_and_Pnum(p0%octCh3%octN,countOct,countPtcl,sendLv,iLv+1,icon)
        call pack_GMesh_and_Pnum(p0%octCh4%octN,countOct,countPtcl,sendLv,iLv+1,icon)
        call pack_GMesh_and_Pnum(p0%octCh5%octN,countOct,countPtcl,sendLv,iLv+1,icon)
        call pack_GMesh_and_Pnum(p0%octCh6%octN,countOct,countPtcl,sendLv,iLv+1,icon)
        call pack_GMesh_and_Pnum(p0%octCh7%octN,countOct,countPtcl,sendLv,iLv+1,icon)
        call pack_GMesh_and_Pnum(p0%octCh8%octN,countOct,countPtcl,sendLv,iLv+1,icon)
     endif
  endif

end subroutine pack_GMesh_and_Pnum

!---------------------------Routines for Refresh Particle END--------------------------

!-----------------------------Routines for DDD--------------------------------
!yagi added 2011/04/07 modify 2012/03/21
subroutine get_moved_octs_size
  use message_passing_interface
  use param
  implicit none
  integer(kind=4)::req(20),istat(MPI_STATUS_SIZE,20)

  !if(debugMode>=1)print *,"===== get_moved_octs_size start =====",rank

!init variables
  next=mod(rank+1,nprocs)
  prev=mod(rank-1+nprocs,nprocs)
  sendSize=0
  recvSize=0
!set parameters
  if(rank/=0)then
     leftDisp=Mn2CPU(prev+1)-(MnDisp)
  else
     leftDisp=0
  endif
  rightDisp=Mn2CPU(rank+1)-(MnDisp+MaxID(1,-1))

  LeftSendMin =1
  LeftSendMax =0
  if(leftDisp>0)LeftSendMax=LeftSendMax+leftDisp
  SavedMin    =MinID(1,-1)
  SavedMax    =MaxID(1,-1)
  RightSendMin=MaxID(1,-1)+1
  RightSendMax=MaxID(1,-1)
  if(leftDisp>0)then
     SavedMin    =SavedMin+leftDisp
  endif
  if(rightDisp<0)then
     SavedMax    =SavedMax+rightDisp
     RightSendMin=RightSendMin+rightDisp
 ! else if(rightDisp>0)then
 !    RightSendMax=RightSendMax+rightDisp
  endif

!!$  if(debugMode>=2)then
!!$     do i=0,nprocs-1
!!$        if(i==rank)then
!!$           print *,"*******************",rank
!!$           print *,"LeftDisp=",leftDisp,rank
!!$           print *,"RightDisp=",rightDisp,rank
!!$           print *,"-------------------",rank
!!$
!!$           print *,"LeftSendMin=",LeftSendMin,rank
!!$           print *,"LeftSendMax=",LeftSendMax,rank
!!$           print *,"SavedMin=",SavedMin,rank
!!$           print *,"SavedMax=",SavedMax,rank
!!$           print *,"RightSendMin=",RightSendMin,rank
!!$           print *,"RightSendMax=",RightSendMax,rank
!!$           print *,"*******************",rank
!!$        endif
!!$        call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$     enddo
!!$  endif

  call mpi_irecv(recvSize(1),2,MPI_INTEGER,prev,1,MPI_COMM_WORLD,req(1),ierr)
  call mpi_irecv(recvSize(3),2,MPI_INTEGER,next,2,MPI_COMM_WORLD,req(2),ierr)

!check num of octs for sending
  if(leftDisp>0)then
     call count_octsT(LeftSendMin ,LeftSendMax ,sendSize(1),sendSize(2))
  endif
  if(rightDisp<0)then
     call count_octsT(RightSendMin,RightSendMax,sendSize(3),sendSize(4))
  endif

!get receive size
  call mpi_isend(sendSize(1),2,MPI_INTEGER,prev,2,MPI_COMM_WORLD,req(3),ierr)
  call mpi_isend(sendSize(3),2,MPI_INTEGER,next,1,MPI_COMM_WORLD,req(4),ierr)

  call mpi_waitall(4,req(1),istat(1,1),ierr)
  call mpi_barrier(MPI_COMM_WORLD,ierr)

!check output
!!$  if(debugMode>=2)then
!!$     do i=0,nprocs-1
!!$        if(i==rank)then
!!$           print "(I4,A,I3,A,I4)",sendSize(1),"<-oct[",rank,"]->",sendSize(3)
!!$           print "(I4,A,I3,A,I4)",recvSize(1),"->oct[",rank,"]<-",recvSize(3)
!!$           print "(I4,A,I3,A,I4)",sendSize(2),"<-ptcl[",rank,"]->",sendSize(4)
!!$           print "(I4,A,I3,A,I4)",recvSize(2),"->ptcl[",rank,"]<-",recvSize(4)
!!$        endif
!!$        call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$     enddo
!!$  endif

 ! if(debugMode>=1)print *,"===== get_moved_octs_size end =====",rank

end subroutine get_moved_octs_size

subroutine exchange_moved_octs
  use message_passing_interface
  use param
  use oct_set
  implicit none
  integer(kind=4)::req(20),istat(MPI_STATUS_SIZE,20)

 ! if(debugMode>=1)print *,"===== exchange_moved_octs start =====",rank


!-----start sending oct(field + particle)-----
!setup irecv
  !field
  call mpi_irecv(rBufL_FR,recvSize(1)*12,MPI_DOUBLE_PRECISION,prev,1,MPI_COMM_WORLD,req(5),ierr)
  call mpi_irecv(rBufR_FR,recvSize(3)*12,MPI_DOUBLE_PRECISION,next,2,MPI_COMM_WORLD,req(6),ierr)
  call mpi_irecv(rBufL_FI,recvSize(1)*3 ,MPI_INTEGER   ,prev,3,MPI_COMM_WORLD,req(7),ierr)
  call mpi_irecv(rBufR_FI,recvSize(3)*3 ,MPI_INTEGER   ,next,4,MPI_COMM_WORLD,req(8),ierr)

  !prtcl
  call mpi_irecv(rbufL_PR,recvSize(2)*9,MPI_DOUBLE_PRECISION,prev,5,MPI_COMM_WORLD,req(9),ierr)
  call mpi_irecv(rbufR_PR,recvSize(4)*9,MPI_DOUBLE_PRECISION,next,6,MPI_COMM_WORLD,req(10),ierr)
  call mpi_irecv(rbufL_PI,recvSize(2)  ,MPI_INTEGER   ,prev,7,MPI_COMM_WORLD,req(11),ierr)
  call mpi_irecv(rbufR_PI,recvSize(4)  ,MPI_INTEGER   ,next,8,MPI_COMM_WORLD,req(12),ierr)

!pack data
  !prev <- me
  if(sendsize(1)>0)then
     
     call pack_octsTL
  endif

  call mpi_isend(sbufL_FR,sendSize(1)*12,MPI_DOUBLE_PRECISION,prev,2,MPI_COMM_WORLD,req(13),ierr)
  call mpi_isend(sbufL_FI,sendSize(1)*3 ,MPI_INTEGER   ,prev,4,MPI_COMM_WORLD,req(14),ierr)
  call mpi_isend(sbufL_PR,sendSize(2)*9 ,MPI_DOUBLE_PRECISION,prev,6,MPI_COMM_WORLD,req(15),ierr)
  call mpi_isend(sbufL_PI,sendSize(2)   ,MPI_INTEGER   ,prev,8,MPI_COMM_WORLD,req(16),ierr)

  !me -> next
  if(sendsize(3)>0)then
   
     call pack_octsTR
  endif
  
  call mpi_isend(sbufR_FR,sendSize(3)*12,MPI_DOUBLE_PRECISION,next,1,MPI_COMM_WORLD,req(17),ierr)
  call mpi_isend(sbufR_FI,sendSize(3)*3 ,MPI_INTEGER   ,next,3,MPI_COMM_WORLD,req(18),ierr)
  call mpi_isend(sbufR_PR,sendSize(4)*9 ,MPI_DOUBLE_PRECISION,next,5,MPI_COMM_WORLD,req(19),ierr)
  call mpi_isend(sbufR_PI,sendSize(4)   ,MPI_INTEGER   ,next,7,MPI_COMM_WORLD,req(20),ierr)
  call mpi_waitall(16,req(5),istat(1,5),ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  !
  !if(debugMode>=1)print *,"===== exchange_moved_octs end =====",rank

end subroutine exchange_moved_octs

!------------- Count Sending Data Number --------------
recursive subroutine count_octs(index,num,num_p)
  use param  
  use oct_set
  use init_mesh_size
  implicit none
  integer(kind=4)::num,num_p
  integer(kind=4),intent(in)::index
  type(oct),pointer::p0

  !print *,"checking index=",index
  p0=>Mesh(index)

  num=num+1
  num_p=num_p+p0%octP

  if(p0%iFLG(1)>0)then
     call count_octs(p0%octCh1%octN,num,num_p)
     call count_octs(p0%octCh2%octN,num,num_p)
     call count_octs(p0%octCh3%octN,num,num_p)
     call count_octs(p0%octCh4%octN,num,num_p)
     call count_octs(p0%octCh5%octN,num,num_p)
     call count_octs(p0%octCh6%octN,num,num_p)
     call count_octs(p0%octCh7%octN,num,num_p)
     call count_octs(p0%octCh8%octN,num,num_p)
  endif

end subroutine count_octs

subroutine count_octsT(sID,eID,num,num_p)
  use param
  use init_mesh_size
  use oct_set
  use message_passing_interface
  implicit none
  integer(kind=4)::num,num_p,index
  integer(kind=4),intent(in)::sID,eID !this variable must be Lv-1 index
  type(oct),pointer::p0

  if(sID<MinID(1,-1) .or. eID>MaxID(1,-1))then
     print *,"illegal sID & eID",sID,eID, rank
     print *,"MinID=",MinID(1,-1),"MaxID=",MaxID(1,-1)
     stop
  endif

  do index=sID,eID
     !p0=>BMesh(index)
     p0=>Mesh(index)
     num=num+1
     num_p=num_p+p0%octP
        call count_octs(p0%octCh1%octN,num,num_p)
        call count_octs(p0%octCh2%octN,num,num_p)
        call count_octs(p0%octCh3%octN,num,num_p)
        call count_octs(p0%octCh4%octN,num,num_p)
        call count_octs(p0%octCh5%octN,num,num_p)
        call count_octs(p0%octCh6%octN,num,num_p)
        call count_octs(p0%octCh7%octN,num,num_p)
        call count_octs(p0%octCh8%octN,num,num_p)
  enddo

end subroutine count_octsT

!------------ Pack Octs -----------

!i will add pack FR pack FI pack PR pack PI

recursive subroutine pack_octsR(index,count,count_p)
  use param
  use oct_set
  use particle_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::count,count_p,i  !count and count_p must be 1 when this subroutine called
  integer(kind=4),intent(in)::index
  type(oct),pointer::p0
  type(prtcl),pointer::pp


  p0=>Mesh(index)
             
  sbufR_FR(1+(count-1)*12:count*12)=p0%F(1:12)
  sbufR_FI(count*3-2)=p0%iFLG(3)
  sbufR_FI(count*3-1)=p0%iFLG(1)
  sbufR_FI(count*3  )=0

  if(p0%iFLG(1)>-4)then

     sbufR_FI(count*3  )=p0%octP    
     if(p0%octP>0)then
        pp=>p0%ptcl
        do i=1,p0%octP
           !if(pp%isort==0)cycle
           pp=>pp%prtnxt
           sbufR_PR(1+(count_p-1)*9:count_p*9)=pp%R(1:9)
           sbufR_PI(count_p)=pp%isort
           count_p=count_p+1

           pp%isort=0
        enddo
     endif
     p0%octP=0
  endif
  count=count+1

  if(p0%iFLG(1)>0)then
     call pack_octsR(p0%octCh1%octN,count,count_p)
     call pack_octsR(p0%octCh2%octN,count,count_p)
     call pack_octsR(p0%octCh3%octN,count,count_p)
     call pack_octsR(p0%octCh4%octN,count,count_p)
     call pack_octsR(p0%octCh5%octN,count,count_p)
     call pack_octsR(p0%octCh6%octN,count,count_p)
     call pack_octsR(p0%octCh7%octN,count,count_p)
     call pack_octsR(p0%octCh8%octN,count,count_p)
  endif

end subroutine pack_octsR

recursive subroutine pack_octsL(index,count,count_p)
  use param
  use oct_set
  use particle_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::count,count_p,i   !count and count_p must be 1 when this subroutine called
  integer(kind=4),intent(in)::index
  type(oct),pointer::p0
  type(prtcl),pointer::pp

  p0=>Mesh(index)
  sbufL_FR(1+(count-1)*12:count*12)=p0%F(1:12)
  sbufL_FI(count*3-2)=p0%iFLG(3)
  sbufL_FI(count*3-1)=p0%iFLG(1)
  sbufL_FI(count*3  )=0

  if(p0%iFLG(1)>-4)then
     sbufL_FI(count*3  )=p0%octP
     if(p0%octP>0)then
        pp=>p0%ptcl
        do i=1,p0%octP
           !if(pp%isort==0)cycle
           pp=>pp%prtnxt
           sbufL_PR(1+(count_p-1)*9:count_p*9)=pp%R(1:9)
           sbufL_PI(count_p)=pp%isort
           count_p=count_p+1

           pp%isort=0
        enddo
     endif
     p0%octP=0
  endif
  count=count+1

  if(p0%iFLG(1)>0)then
     call pack_octsL(p0%octCh1%octN,count,count_p)
     call pack_octsL(p0%octCh2%octN,count,count_p)
     call pack_octsL(p0%octCh3%octN,count,count_p)
     call pack_octsL(p0%octCh4%octN,count,count_p)
     call pack_octsL(p0%octCh5%octN,count,count_p)
     call pack_octsL(p0%octCh6%octN,count,count_p)
     call pack_octsL(p0%octCh7%octN,count,count_p)
     call pack_octsL(p0%octCh8%octN,count,count_p)
  endif

end subroutine pack_octsL

subroutine pack_octsTR
  use param
  use init_mesh_size
  use message_passing_interface
  use oct_set
  implicit none
  integer(kind=4)::index,count,count_p
  integer(kind=4)::sID,eID
  type(oct),pointer::p0

  sID=RightSendMin
  eID=RightSendMax

!!$  if(debugMode>=3)then
!!$     print *,"Start pack_octsT for R buf Datanum=",size(sbufR_FR),rank
!!$     print *,"sID=",sID,"eID=",eID,rank
!!$  endif

  if(sID<MinID(1,-1) .or. eID>MaxID(1,-1))then
     print *,"illegal sID & eID",sID,eID
     print *,"MinID=",MinID(1,-1),"MaxID=",MaxID(1,-1)
     stop
  endif

  count=1
  count_p=1

  do index=sID,eID
     !p0=>BMesh(index)
     p0=>Mesh(index)

     sbufR_FR(1+(count-1)*12:count*12)=p0%F(1:12)
     sbufR_FI(count*3-2)=p0%iFLG(3)
     sbufR_FI(count*3-1)=p0%iFLG(1)
     sbufR_FI(count*3  )=0
     count=count+1

     call pack_octsR(p0%octCh1%octN,count,count_p)
     call pack_octsR(p0%octCh2%octN,count,count_p)
     call pack_octsR(p0%octCh3%octN,count,count_p)
     call pack_octsR(p0%octCh4%octN,count,count_p)
     call pack_octsR(p0%octCh5%octN,count,count_p)
     call pack_octsR(p0%octCh6%octN,count,count_p)
     call pack_octsR(p0%octCh7%octN,count,count_p)
     call pack_octsR(p0%octCh8%octN,count,count_p)
  enddo

!!$  if(debugMode>=3)print *,"[Pack_octsTR] oct count=",count,rank
!!$  if(debugMode>=3)print *,"[Pack_octsTR] ptcl count=",count_p,rank

end subroutine pack_octsTR

subroutine pack_octsTL
  use param
  use init_mesh_size
  use message_passing_interface
  use oct_set
  implicit none
  integer(kind=4)::index,count,count_p
  integer(kind=4)::sID,eID
  type(oct),pointer::p0

  sID=LeftSendMin
  eID=LeftSendMax

!!$  if(debugMode>=3)then
!!$     print *,"Start pack_octsT for L buf Datanum=",size(sbufL_FR),rank
!!$     print *,"sID=",sID,"eID=",eID,rank
!!$  endif

  if(sID<MinID(1,-1) .or. eID>MaxID(1,-1))then
     print *,"illegal sID & eID",sID,eID
     print *,"MinID=",MinID(1,-1),"MaxID=",MaxID(1,-1)
     stop
  endif

  count=1
  count_p=1

  do index=sID,eID
     !p0=>BMesh(index)
     p0=>Mesh(index)

     sbufL_FR(1+(count-1)*12:count*12)=p0%F(1:12)
     sbufL_FI(count*3-2)=p0%iFLG(3)
     sbufL_FI(count*3-1)=p0%iFLG(1)
     sbufL_FI(count*3  )=0
     count=count+1

     call pack_octsL(p0%octCh1%octN,count,count_p)
     call pack_octsL(p0%octCh2%octN,count,count_p)
     call pack_octsL(p0%octCh3%octN,count,count_p)
     call pack_octsL(p0%octCh4%octN,count,count_p)
     call pack_octsL(p0%octCh5%octN,count,count_p)
     call pack_octsL(p0%octCh6%octN,count,count_p)
     call pack_octsL(p0%octCh7%octN,count,count_p)
     call pack_octsL(p0%octCh8%octN,count,count_p)
  enddo

!!$  if(debugMode>=3)print *,"[Pack_octsTL] oct count=",count,rank
!!$  if(debugMode>=3)print *,"[Pack_octsTL] ptcl count=",count_p,rank

end subroutine pack_octsTL
!------------------------------

!------ Set Received Data -----
!last mod. 2011/06/30
!now modifying 2012/03/21
subroutine set_moved_octs
  use param
  use message_passing_interface
  use oct_set
  use particle_set
  use init_mesh_size
  implicit none
  integer(kind=4)::index,i,count
  type(oct),pointer::p0,p1

  integer(kind=4)::index2,temp,CopyMin,CopyMax,intv
  integer(kind=4)::slideSize

  integer(kind=4),dimension(0:LvMax)::setNum

  !if(debugMode>=1)print *,"===== Start set_moved_octs =====",rank

  intv=1
  LeftNum=0
  RightNum=0

!count detail receive data number
  if(recvSize(1)>0)then
     call count_recv_dataLT
  endif
  if(recvSize(3)>0)then
     call count_recv_dataRT
  endif
  
!!$  if(debugMode>=2)then
!!$     call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$     do i=0,nprocs-1
!!$        if(i==rank)then
!!$           print *,"*******************",rank
!!$           do j=-1,LvMax
!!$              Print '(A,I2,A,I4,I3)',"Left Num(",j,")=",LeftNum(j),rank
!!$              print '(A,I2,A,I4,I3)',"RightNum(",j,")=",RightNum(j),rank
!!$           enddo
!!$           print *,"*******************",rank
!!$        endif
!!$        call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$     enddo
!!$  endif

!set slide amount of BMesh
  slideSize=leftDisp*(-1)
!set new Lv:-1 index of octs
  LeftRecvMin =1
  LeftRecvMax =LeftNum(-1)
  ExistMin    =LeftNum(-1)+MinID(1,-1)
  ExistMax    =LeftNum(-1)+MaxID(1,-1)
  if(leftDisp>0)ExistMax=ExistMax-leftDisp
  if(rightDisp<0)ExistMax=ExistMax+rightDisp
  RightRecvMin=ExistMax+1
  RightRecvMax=ExistMax+RightNum(-1)

!set slide area
  CopyMin=SavedMin
  CopyMax=SavedMax

!!$  if(debugMode>=2)then
!!$     call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$     do i=0,nprocs-1
!!$        if(i==rank)then
!!$           print *,"*******************",rank
!!$           print *,"slideSize=",slideSize,rank
!!$           print *,"-------------------",rank
!!$           print *,"LeftRecvMin=",LeftRecvMin,rank
!!$           print *,"LeftRecvMax=",LeftRecvMax,rank
!!$           print *,"ExistMin=",ExistMin,rank
!!$           print *,"ExistMax=",ExistMax,rank
!!$           print *,"RightRecvMin=",RightRecvMin,rank
!!$           print *,"RightRecvMax=",RightRecvMax,rank
!!$           print *,"-------------------",rank
!!$           print *,"CopyMin=",CopyMin,rank
!!$           print *,"CopyMax=",CopyMax,rank
!!$           print *,"*******************",rank
!!$        endif
!!$        call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$     enddo
!!$  endif


!delete non using domain
  call set_deleted_oct_flagT(CopyMin,CopyMax)
  !if(debugMode>=3)print *,"set_deleted_oct_flag completed",rank

!change size of BMesh  
!!$if(debugMode>=2)then
!!$   print *,"leftDisp=",leftDisp,"MaxID(1,-1)=",MaxID(1,-1),"rightDisp=",rightDisp
!!$   print *,"BMeshSize=",BMeshSize
!!$endif

if((MaxID(1,-1)+leftDisp+rightDisp)>BMeshSize)then
   print *,"BMeshSize overflow BMeshSize:",BMeshSize," < newBMeshSize",MaxID(1,-1)+leftDisp+rightDisp
   stop
endif


!check copying direction
  if(slideSize>0)then
     temp=CopyMin
     CopyMin=CopyMax
     CopyMax=temp
     intv=-1
  endif

  count=0

!start slide existing octs 
  if(slideSize/=0)then
     do index=CopyMin,CopyMax,intv
        index2=index+slideSize
        !p0=>BMesh(index)
        !p1=>BMesh(index2)
        p0=>Mesh(index)
        p1=>Mesh(index2)

      
!copy data     
        p1=p0
!refine flags
        p1%octN=index2
        p0%iFLG=-4

!delete non used data
        p0%iPOS=-1
        p0%rPOS=-1

!reset connection
        nullify(p0%octNb1)
        nullify(p0%octNb2)
        nullify(p0%octNb3)
        nullify(p0%octNb4)
        nullify(p0%octNb5)
        nullify(p0%octNb6)

        nullify(p0%octCh1)
        nullify(p0%octCh2)
        nullify(p0%octCh3)
        nullify(p0%octCh4)
        nullify(p0%octCh5)
        nullify(p0%octCh6)
        nullify(p0%octCh7)
        nullify(p0%octCh8)

        nullify(p0%ptcl)

!set new connection of child
        p1%octCh1%octPrt=>p1
        p1%octCh2%octPrt=>p1
        p1%octCh3%octPrt=>p1
        p1%octCh4%octPrt=>p1
        p1%octCh5%octPrt=>p1
        p1%octCh6%octPrt=>p1
        p1%octCh7%octPrt=>p1
        p1%octCh8%octPrt=>p1
     enddo
  else
     do index=copyMin,copyMax,intv
        !p0=>BMesh(index)
        p0=>Mesh(index)
        !reset connection
        p0%octCh1%octPrt=>p0
        p0%octCh2%octPrt=>p0
        p0%octCh3%octPrt=>p0
        p0%octCh4%octPrt=>p0
        p0%octCh5%octPrt=>p0
        p0%octCh6%octPrt=>p0
        p0%octCh7%octPrt=>p0
        p0%octCh8%octPrt=>p0
     enddo
  endif
!change new data domain
  ! I must modify this part
  MnDisp=MnDisp+leftDisp
  MinID(1,-1)=1
  MaxID(1,-1)=ExistMax+RightNum(-1)
  MinID(2,-1)=MaxID(1,-1)+1
  MaxID(2,-1)=MaxID(1,-1)


  MinID(2,0)=MaxID(1,LvMax)+1
  MaxID(2,0)=MinID(2,0)+LeftNum(0)+RightNum(0)-1

  if(LvMax>0)then
     do i=1,LvMax
        MinID(2,i)=MaxID(2,i-1)+1
        MaxID(2,i)=MinID(2,i)+LeftNum(i)+RightNum(i)-1  
     enddo
  endif



!!$  if(debugMode>=2)then
!!$     do i=0,nprocs-1
!!$        if(rank==i)then
!!$           print *,"MnDisp=",MnDisp,rank
!!$           print *,"-- now MinID,MaxID is --",rank
!!$           do j=-1,LvMax
!!$              print '(A,I2,A,I8,I3)',"minID(1,",j,")=",minID(1,j),rank
!!$              print '(A,I2,A,I8,I3)',"maxID(1,",j,")=",maxID(1,j),rank
!!$              print '(A,I2,A,I8,I3)',"minID(2,",j,")=",minID(2,j),rank
!!$              print '(A,I2,A,I8,I3)',"maxID(2,",j,")=",maxID(2,j),rank
!!$           enddo
!!$           print *,"------------------------",rank
!!$        endif
!!$        call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$     enddo
!!$  endif

  setNum=0

  if(leftDisp<0)call set_packed_octTL(setNum)
!!$  if(debugMode>=3)print *,"set_packed_octTLcomp",rank
!!$  if(debugMode>=2)then
!!$     print *,"*******************",rank
!!$     do k=0,LvMax
!!$        Print '(A,I2,A,I4,I3)',"Left Added Num(",k,")=",setNum(k),rank
!!$     enddo
!!$     print *,"*******************",rank
!!$  endif

  if(rightDisp>0)call set_packed_octTR(setNum)
!!$  if(debugMode>=3)print *,"set_packed_octTRcomp",rank
!!$  if(debugMode>=2)then
!!$     print *,"*******************",rank
!!$     do k=0,LvMax
!!$        Print '(A,I2,A,I4,I3)',"Right Added Num(",k,")=",setNum(k),rank
!!$     enddo
!!$     print *,"*******************",rank
!!$  endif

 ! if(debugMode>=1)print *,"===== set_moved_octs end =====",rank

end subroutine set_moved_octs


!yagi added 2011/05/06

recursive subroutine set_deleted_oct_flag(index)
  use message_passing_interface
  use init_mesh_size
  use oct_set
  implicit none
  integer(kind=4),intent(in)::index
  integer(kind=4)::i,pnum
  type(oct),pointer::p0
  type(prtcl),pointer::pp,ppnxt

  !print *,"checking index=",index,rank
  p0=>Mesh(index)
  
  if(p0%iFLG(1)>0)then
     call set_deleted_oct_flag(p0%octCh1%octN)
     call set_deleted_oct_flag(p0%octCh2%octN)
     call set_deleted_oct_flag(p0%octCh3%octN)
     call set_deleted_oct_flag(p0%octCh4%octN)
     call set_deleted_oct_flag(p0%octCh5%octN)
     call set_deleted_oct_flag(p0%octCh6%octN)
     call set_deleted_oct_flag(p0%octCh7%octN)
     call set_deleted_oct_flag(p0%octCh8%octN)
     nullify(p0%octCh1)
     nullify(p0%octCh2)
     nullify(p0%octCh3)
     nullify(p0%octCh4)
     nullify(p0%octCh5)
     nullify(p0%octCh6)
     nullify(p0%octCh7)
     nullify(p0%octCh8)
  endif

  nullify(p0%octPrt)

  p0%iFLG=-4
  p0%iPOS=-1
  p0%rPOS=-1
  
  !delete particle
  pnum=p0%octP
  if(associated(p0%ptcl))then
     if(associated(p0%ptcl%prtnxt))then
        pp=>p0%ptcl
        do i=1,pnum+1
           pp%isort=0
           ppnxt=>pp%prtnxt
           nullify(pp%prtnxt)
           pp=>ppnxt
        enddo
     endif
    
     nullify(p0%ptcl)
  endif
  p0%octP=0

  if(associated(p0%octNb1))then
     nullify(p0%octNb1%octNb2)
     nullify(p0%octNb1)
  endif
  if(associated(p0%octNb2))then
     nullify(p0%octNb2%octNb1)
     nullify(p0%octNb2)
  endif
  if(associated(p0%octNb3))then
     nullify(p0%octNb3%octNb4)
     nullify(p0%octNb3)
  endif
  if(associated(p0%octNb4))then
     nullify(p0%octNb4%octNb3)
     nullify(p0%octNb4)
  endif
  if(associated(p0%octNb5))then
     nullify(p0%octNb5%octNb6)
     nullify(p0%octNb5)
  endif
  if(associated(p0%octNb6))then
     nullify(p0%octNb6%octNb5)
     nullify(p0%octNb6)
  endif

end subroutine set_deleted_oct_flag

subroutine set_deleted_oct_flagT(SaveSID,SaveEID)
  use message_passing_interface
  use init_mesh_size
  use oct_set
  implicit none
  integer(kind=4),intent(in)::SaveSID,SaveEID
  integer(kind=4)::index
  type(oct),pointer::p0

!!$  if(debugMode>=3)then
!!$     print *,"clean BMesh",1,"-",BMeshSize,rank
!!$     print *,"except",SaveSID,"-",SaveEID,rank
!!$  endif

 ! do index=minID(1,-1),maxID(1,-1)
  do index=1,BMeshSize
     if(index>SaveeID .or. index<SavesID)then
        !print *,"cheking indexLv-1=",index,rank
        !p0=>BMesh(index)
        p0=>Mesh(index)

        if(p0%iFLG(1)>0)then
           call set_deleted_oct_flag(p0%octCh1%octN)
           call set_deleted_oct_flag(p0%octCh2%octN)
           call set_deleted_oct_flag(p0%octCh3%octN)
           call set_deleted_oct_flag(p0%octCh4%octN)
           call set_deleted_oct_flag(p0%octCh5%octN)
           call set_deleted_oct_flag(p0%octCh6%octN)
           call set_deleted_oct_flag(p0%octCh7%octN)
           call set_deleted_oct_flag(p0%octCh8%octN)
        endif



        if(associated(p0%octNb1))then
           nullify(p0%octNb1%octNb2)
           nullify(p0%octNb1)
        endif
        if(associated(p0%octNb2))then
           nullify(p0%octNb2%octNb1)
           nullify(p0%octNb2)
        endif
        if(associated(p0%octNb3))then
           nullify(p0%octNb3%octNb4)
           nullify(p0%octNb3)
        endif
        if(associated(p0%octNb4))then
           nullify(p0%octNb4%octNb3)
           nullify(p0%octNb4)
        endif
        if(associated(p0%octNb5))then
           nullify(p0%octNb5%octNb6)
           nullify(p0%octNb5)
        endif
        if(associated(p0%octNb6))then
           nullify(p0%octNb6%octNb5)
           nullify(p0%octNb6)
        endif

        nullify(p0%octCh1)
        nullify(p0%octCh2)
        nullify(p0%octCh3)
        nullify(p0%octCh4)
        nullify(p0%octCh5)
        nullify(p0%octCh6)
        nullify(p0%octCh7)
        nullify(p0%octCh8)


        !print *,"here!"
        p0%iFLG=-4
        p0%iPOS=-1
        p0%rPOS=-1
        !print *,"here!"
     endif
  enddo

end subroutine set_deleted_oct_flagT

!yagi added 2011/04/11
recursive subroutine count_recv_dataR(count,iLv)
  use const
  use param
  use message_passing_interface
  implicit none
  integer(kind=4),intent(in)::iLv
  integer(kind=4)::count

  RightNum(iLv)=RightNum(iLv)+1
  count=count+1
  if(rbufR_FI(count*3-1)>0)then
     call count_recv_dataR(count,iLv+1)
     call count_recv_dataR(count,iLv+1)
     call count_recv_dataR(count,iLv+1)
     call count_recv_dataR(count,iLv+1)
     call count_recv_dataR(count,iLv+1)
     call count_recv_dataR(count,iLv+1)
     call count_recv_dataR(count,iLv+1)
     call count_recv_dataR(count,iLv+1)
  endif
end subroutine count_recv_dataR

recursive subroutine count_recv_dataL(count,iLv)
  use const
  use param
  use message_passing_interface
  implicit none
  integer(kind=4),intent(in)::iLv
  integer(kind=4)::count

  LeftNum(iLv)=LeftNum(iLv)+1
  count=count+1
  if(rbufL_FI(count*3-1)>0)then
     call count_recv_dataL(count,iLv+1)
     call count_recv_dataL(count,iLv+1)
     call count_recv_dataL(count,iLv+1)
     call count_recv_dataL(count,iLv+1)
     call count_recv_dataL(count,iLv+1)
     call count_recv_dataL(count,iLv+1)
     call count_recv_dataL(count,iLv+1)
     call count_recv_dataL(count,iLv+1)
  endif
end subroutine count_recv_dataL
!!$
subroutine count_recv_dataRT
  use param
  use message_passing_interface
  implicit none
  integer(kind=4)::count
  count=0
  do while(count<recvSize(3))
     call count_recv_dataR(count,-1)
  enddo
end subroutine count_recv_dataRT

subroutine count_recv_dataLT
  use param
  use message_passing_interface
  implicit none
  integer(kind=4)::count
  count=0
  do while(count<recvSize(1))
     call count_recv_dataL(count,-1)
  enddo
end subroutine count_recv_dataLT

!refine 2011/5/10 modify 2012/04/18
recursive subroutine set_packed_octL(count,count_p,num,sID,iLv,Cs)
  use param
  use message_passing_interface
  use oct_set
  use particle_set
  use const
  implicit none
  integer(kind=4),intent(in)::iLv,Cs,sID
  integer(kind=4),dimension(0:LvMax)::num
  integer(kind=4)::count,count_p,index,i,indexP,Isort
  real(kind=8)::ixx,iyy,izz,R(9)
  type(oct),pointer::p0
  type(prtcl),pointer::pp,PrtList

  nullify(PrtList)


  index=sID+num(iLv)
  num(iLv)=num(iLv)+1

  p0=>Mesh(index)
  p0%octN   =index
  p0%Csort  =Cs
  p0%octLv  =iLv
  p0%MrtN   =p0%octPrt%MrtN
  !--set octP--
  p0%octP   =rBufL_FI(count*3  )
  !--set iFLG--
  p0%iFLG(1)=rBufL_FI(count*3-1)
  p0%iFLG(3)=rBufL_FI(count*3-2)
  if((p0%iFLG(1).ge.1).and.(p0%iFLG(1).le.3)) then 
     p0%iFLG(2)=5
  else
     p0%iFLG(2)=0
  endif
  !p0%iFLG(3)=0
  p0%octType = 0
  
  !--set iPos--
  p0%iPos=p0%octPrt%iPos

  do i=1,3
     p0%iPos(i)=p0%iPos(i)+intvLv(Cs,i,iLv)
  enddo

  !-----translate iPos to rPos------
  ixx = real(p0%iPos(1)+Nintv)/real(2*Nintv)
  iyy = real(p0%iPos(2)+Nintv)/real(2*Nintv)
  izz = real(p0%iPos(3)+Nintv)/real(2*Nintv) 
  p0%rPOS(1) = ixx*dx(1) - HALF*dx(1)
  p0%rPOS(2) = iyy*dx(2) - HALF*dx(2)
  p0%rPOS(3) = izz*dx(3) - HALF*dx(3)
  !---------------------------------


  !---set field parameter---
  p0%F(1:12)=rBufL_FR(1+(count-1)*12:count*12)

  count=count+1 !care to change

  !=================================
  indexP=MaxIP(iLv)
  !---set representative particle---
  if(associated(p0%ptcl))then
!!$     if(p0%ptcl%isort/=-1)then
!!$        print *,"oct initializing error in set moved oct",rank
!!$        stop
!!$     endif
     !temporary
     nullify(p0%ptcl)
  endif
  if(.not. associated(p0%ptcl))then
     indexP=indexP+1
     R(1:3)=p0%rPos(1:3)
     R(4:9)=0.0d0
     Isort=-1
     Pesh(indexP,iLv)=prtcl(R,Isort,p0%octN,indexP,PrtList,0)
     p0%ptcl=>Pesh(indexP,iLv)
  endif
  !--------------------------------

  !-----injecting particle----
  if(p0%iFLG(1)>=0 .and. p0%iFLG(1)<4)then !particle will exist only cells which has no child
     pp=>p0%ptcl%prtnxt
     do i=1,p0%octP
       
        indexP=indexP+1
        R(1:9)=rBufL_PR(1+(count_p-1)*9:count_p*9)
        Isort=rBufL_PI(count_p)
        Pesh(indexP,iLv)=prtcl(R,Isort,p0%octN,indexP,pp,0)
        p0%ptcl%prtnxt => Pesh(indexP,iLv)
        pp=>p0%ptcl%prtnxt
        count_p=count_p+1
     enddo
  endif
  MaxIP(iLv)=indexP
  !=================================


  !connect pointer
  if(p0%iFLG(1)>0)then
     p0%octCh1=>Mesh(MinID(2,iLv+1)+num(iLv+1)+0)
     p0%octCh2=>Mesh(MinID(2,iLv+1)+num(iLv+1)+1)
     p0%octCh3=>Mesh(MinID(2,iLv+1)+num(iLv+1)+2)
     p0%octCh4=>Mesh(MinID(2,iLv+1)+num(iLv+1)+3)
     p0%octCh5=>Mesh(MinID(2,iLv+1)+num(iLv+1)+4)
     p0%octCh6=>Mesh(MinID(2,iLv+1)+num(iLv+1)+5)
     p0%octCh7=>Mesh(MinID(2,iLv+1)+num(iLv+1)+6)
     p0%octCh8=>Mesh(MinID(2,iLv+1)+num(iLv+1)+7)
     p0%octCh1%octPrt=>p0
     p0%octCh2%octPrt=>p0
     p0%octCh3%octPrt=>p0
     p0%octCh4%octPrt=>p0
     p0%octCh5%octPrt=>p0
     p0%octCh6%octPrt=>p0
     p0%octCh7%octPrt=>p0
     p0%octCh8%octPrt=>p0
     call set_packed_octL(count,count_p,num,minID(2,iLv+1),iLv+1,1)
     call set_packed_octL(count,count_p,num,minID(2,iLv+1),iLv+1,2)
     call set_packed_octL(count,count_p,num,minID(2,iLv+1),iLv+1,3)
     call set_packed_octL(count,count_p,num,minID(2,iLv+1),iLv+1,4)
     call set_packed_octL(count,count_p,num,minID(2,iLv+1),iLv+1,5)
     call set_packed_octL(count,count_p,num,minID(2,iLv+1),iLv+1,6)
     call set_packed_octL(count,count_p,num,minID(2,iLv+1),iLv+1,7)
     call set_packed_octL(count,count_p,num,minID(2,iLv+1),iLv+1,8)
  endif

end subroutine set_packed_octL

recursive subroutine set_packed_octR(count,count_p,num,sID,iLv,Cs)
  use param
  use message_passing_interface
  use oct_set
  use particle_set
  use const
  implicit none
  integer(kind=4),intent(in)::iLv,Cs,sID
  integer(kind=4),dimension(0:LvMax)::num
  integer(kind=4)::count,count_p,index,i,indexP,Isort
  real(kind=8)::ixx,iyy,izz,R(9)
  type(oct),pointer::p0
  type(prtcl),pointer::pp,PrtList
  
  nullify(PrtList)

  index=sID+num(iLv)
  num(iLv)=num(iLv)+1

  p0=>Mesh(index)
  p0%octN   =index
  p0%Csort  =Cs
  p0%octLv  =iLv
  p0%MrtN   =p0%octPrt%MrtN 
  p0%octP   =rBufR_FI(count*3  )
  p0%iFLG(1)=rBufR_FI(count*3-1)
  p0%iFLG(3)=rBufR_FI(count*3-2)
  if((p0%iFLG(1).ge.1).and.(p0%iFLG(1).le.3)) then 
     p0%iFLG(2)=5
  else
     p0%iFLG(2)=0
  endif
  !p0%iFLG(3)=0
  p0%octType = 0
 
  !----set iPos----
  p0%iPos=p0%octPrt%iPos
  do i=1,3
     p0%iPos(i)=p0%iPos(i)+intvLv(Cs,i,iLv)
  enddo

  !-----translate iPos to rPos------
  ixx = real(p0%iPos(1)+Nintv)/real(2*Nintv)
  iyy = real(p0%iPos(2)+Nintv)/real(2*Nintv)
  izz = real(p0%iPos(3)+Nintv)/real(2*Nintv) 
  p0%rPOS(1) = ixx*dx(1) - HALF*dx(1)
  p0%rPOS(2) = iyy*dx(2) - HALF*dx(2)
  p0%rPOS(3) = izz*dx(3) - HALF*dx(3)
 
  p0%F(1:12)=rBufR_FR(1+(count-1)*12:count*12)
  
  count=count+1
  
  !==========set particles============
  indexP=MaxIP(iLv)
  !----set representative particle----
  if(associated(p0%ptcl))then
!!$     if(p0%ptcl%isort/=-1)then
!!$        print *,"oct initializing error in set moved oct",rank
!!$        stop
!!$     endif
     !temporary
     nullify(p0%ptcl)
  endif
  if(.not. associated(p0%ptcl))then
     indexP=indexP+1
     R(1:3)=p0%rPos(1:3)
     R(4:9)=0.0d0
     Isort=-1
     Pesh(indexP,iLv)=prtcl(R,Isort,p0%octN,indexP,PrtList,0)
     p0%ptcl=>Pesh(indexP,iLv)
  endif

  !----inject particle----
  if(p0%iFLG(1)>=0 .and. p0%iFLG(1)<4)then
     pp=>p0%ptcl%prtnxt
    
     do i=1,p0%octP
      
        indexP=indexP+1
        R(1:9)=rBufR_PR(1+(count_p-1)*9:count_p*9)
        Isort=rBufR_PI(count_p)
        Pesh(indexP,iLv)=prtcl(R,Isort,p0%octN,indexP,pp,0)
        p0%ptcl%prtnxt => Pesh(indexP,iLv)
        pp=>p0%ptcl%prtnxt
        count_p=count_p+1
     enddo
  endif
  
  MaxIP(iLv)=indexP
  !===================================

  !connect pointer
  if(p0%iFLG(1)>0)then
     p0%octCh1=>Mesh(MinID(2,iLv+1)+num(iLv+1)+0)
     p0%octCh2=>Mesh(MinID(2,iLv+1)+num(iLv+1)+1)
     p0%octCh3=>Mesh(MinID(2,iLv+1)+num(iLv+1)+2)
     p0%octCh4=>Mesh(MinID(2,iLv+1)+num(iLv+1)+3)
     p0%octCh5=>Mesh(MinID(2,iLv+1)+num(iLv+1)+4)
     p0%octCh6=>Mesh(MinID(2,iLv+1)+num(iLv+1)+5)
     p0%octCh7=>Mesh(MinID(2,iLv+1)+num(iLv+1)+6)
     p0%octCh8=>Mesh(MinID(2,iLv+1)+num(iLv+1)+7)
     p0%octCh1%octPrt=>p0
     p0%octCh2%octPrt=>p0
     p0%octCh3%octPrt=>p0
     p0%octCh4%octPrt=>p0
     p0%octCh5%octPrt=>p0
     p0%octCh6%octPrt=>p0
     p0%octCh7%octPrt=>p0
     p0%octCh8%octPrt=>p0
     call set_packed_octR(count,count_p,num,minID(2,iLv+1),iLv+1,1)
     call set_packed_octR(count,count_p,num,minID(2,iLv+1),iLv+1,2)
     call set_packed_octR(count,count_p,num,minID(2,iLv+1),iLv+1,3)
     call set_packed_octR(count,count_p,num,minID(2,iLv+1),iLv+1,4)
     call set_packed_octR(count,count_p,num,minID(2,iLv+1),iLv+1,5)
     call set_packed_octR(count,count_p,num,minID(2,iLv+1),iLv+1,6)
     call set_packed_octR(count,count_p,num,minID(2,iLv+1),iLv+1,7)
     call set_packed_octR(count,count_p,num,minID(2,iLv+1),iLv+1,8)
  endif

end subroutine set_packed_octR

subroutine set_packed_octTL(num)
  use const
  use param
  use message_passing_interface
  use init_mesh_size
  use oct_set
  implicit none
  integer(kind=4),dimension(0:LvMax)::num
  integer(kind=4)::count,count_p,i,sID,eID,index
  integer(kind=4),dimension(3)::iPos
  real(kind=8)::ixx,iyy,izz
  type(oct),pointer::p0

  count=1
  count_p=1

  sID=LeftRecvMin
  eID=LeftRecvMax

  do i=sID,eID
 
     index=i
     
     !p0=>BMesh(index)
     p0=>Mesh(index)
     p0%octN   =index
     p0%octLv  =-1
     p0%MrtN   =index+MnDisp
     !if(istep==230)print *,"Mn is set",p0%MrtN
     p0%octP   =rBufL_FI(count*3  )
     p0%iFLG(1)=rBufL_FI(count*3-1)
     p0%iFLG(3)=rBufL_FI(count*3-2)
     p0%iFLG(2)=0
     p0%octType = -1
     p0%F(1:12)=rBufL_FR(1+(count-1)*12:count*12)
     count=count+1

     !call Inverse_Morton(p0%MrtN,iPos)
     call get_InverseIndex(p0%MrtN,iPos,0)
     p0%iPos=iPos
          
     ixx = real(p0%iPos(1)+Nintv)/real(2*Nintv)
     iyy = real(p0%iPos(2)+Nintv)/real(2*Nintv)
     izz = real(p0%iPos(3)+Nintv)/real(2*Nintv) 
     p0%rPOS(1) = ixx*dx(1) - HALF*dx(1)
     p0%rPOS(2) = iyy*dx(2) - HALF*dx(2)
     p0%rPOS(3) = izz*dx(3) - HALF*dx(3)

     nullify(p0%ptcl)

     if(p0%iFLG(1)>0)then
        p0%octCh1=>Mesh(MinID(2,0)+num(0)+0)
        p0%octCh2=>Mesh(MinID(2,0)+num(0)+1)
        p0%octCh3=>Mesh(MinID(2,0)+num(0)+2)
        p0%octCh4=>Mesh(MinID(2,0)+num(0)+3)
        p0%octCh5=>Mesh(MinID(2,0)+num(0)+4)
        p0%octCh6=>Mesh(MinID(2,0)+num(0)+5)
        p0%octCh7=>Mesh(MinID(2,0)+num(0)+6)
        p0%octCh8=>Mesh(MinID(2,0)+num(0)+7)
        p0%octCh1%octPrt=>p0
        p0%octCh2%octPrt=>p0
        p0%octCh3%octPrt=>p0
        p0%octCh4%octPrt=>p0
        p0%octCh5%octPrt=>p0
        p0%octCh6%octPrt=>p0
        p0%octCh7%octPrt=>p0
        p0%octCh8%octPrt=>p0
        call set_packed_octL(count,count_p,num,minID(2,0),0,1)
        call set_packed_octL(count,count_p,num,minID(2,0),0,2)
        call set_packed_octL(count,count_p,num,minID(2,0),0,3)
        call set_packed_octL(count,count_p,num,minID(2,0),0,4)
        call set_packed_octL(count,count_p,num,minID(2,0),0,5)
        call set_packed_octL(count,count_p,num,minID(2,0),0,6)
        call set_packed_octL(count,count_p,num,minID(2,0),0,7)
        call set_packed_octL(count,count_p,num,minID(2,0),0,8)
     endif

  enddo

  !if(debugMode>=1)print *,"set_packed_octTL completed",rank
end subroutine set_packed_octTL

subroutine set_packed_octTR(num)
  use const
  use param
  use message_passing_interface
  use init_mesh_size
  use oct_set
  implicit none
  integer(kind=4),dimension(0:LvMax)::num
  integer(kind=4)::count,count_p,i,sID,eID,index
  integer(kind=4),dimension(3)::iPos
  real(kind=8)::ixx,iyy,izz
  type(oct),pointer::p0

  count=1
  count_p=1

  sID=RightRecvMin
  eID=RightRecvMax

  do i=sID,eID
     index=i
     
     !p0=>BMesh(index)
     p0=>Mesh(index)
     p0%octN   =index
     p0%octLv  =-1
     p0%MrtN   =index+MnDisp
     p0%octP   =rBufR_FI(count*3  )
     p0%iFLG(1)=rBufR_FI(count*3-1)
     p0%iFLG(3)=rBufR_FI(count*3-2)
     p0%iFLG(2)=0

     p0%octType = -1
     p0%F(1:12)=rBufR_FR(1+(count-1)*12:count*12)
     count=count+1

     !call Inverse_Morton(p0%MrtN,iPos)
     call get_InverseIndex(p0%MrtN,iPos,0)
     p0%iPos=iPos
     
     ixx = real(p0%iPos(1)+Nintv)/real(2*Nintv)
     iyy = real(p0%iPos(2)+Nintv)/real(2*Nintv)
     izz = real(p0%iPos(3)+Nintv)/real(2*Nintv) 
     p0%rPOS(1) = ixx*dx(1) - HALF*dx(1)
     p0%rPOS(2) = iyy*dx(2) - HALF*dx(2)
     p0%rPOS(3) = izz*dx(3) - HALF*dx(3)

     nullify(p0%ptcl)

     if(p0%iFLG(1)>0)then
        !if(rank==1)print *,"Right Parent octN=",p0%octN,"child octN=",MinID(2,0)+num(0)+0
        p0%octCh1=>Mesh(MinID(2,0)+num(0)+0)
        p0%octCh2=>Mesh(MinID(2,0)+num(0)+1)
        p0%octCh3=>Mesh(MinID(2,0)+num(0)+2)
        p0%octCh4=>Mesh(MinID(2,0)+num(0)+3)
        p0%octCh5=>Mesh(MinID(2,0)+num(0)+4)
        p0%octCh6=>Mesh(MinID(2,0)+num(0)+5)
        p0%octCh7=>Mesh(MinID(2,0)+num(0)+6)
        p0%octCh8=>Mesh(MinID(2,0)+num(0)+7)
        p0%octCh1%octPrt=>p0
        !if(rank==1)print *,"Right check Prt octN=",p0%octCh1%octPrt%octN
        p0%octCh2%octPrt=>p0
        p0%octCh3%octPrt=>p0
        p0%octCh4%octPrt=>p0
        p0%octCh5%octPrt=>p0
        p0%octCh6%octPrt=>p0
        p0%octCh7%octPrt=>p0
        p0%octCh8%octPrt=>p0
        call set_packed_octR(count,count_p,num,minID(2,0),0,1)
        call set_packed_octR(count,count_p,num,minID(2,0),0,2)
        call set_packed_octR(count,count_p,num,minID(2,0),0,3)
        call set_packed_octR(count,count_p,num,minID(2,0),0,4)
        call set_packed_octR(count,count_p,num,minID(2,0),0,5)
        call set_packed_octR(count,count_p,num,minID(2,0),0,6)
        call set_packed_octR(count,count_p,num,minID(2,0),0,7)
        call set_packed_octR(count,count_p,num,minID(2,0),0,8)
     endif
  enddo

  !if(debugMode>=1)print *,"set_packed_octTR completed",rank
end subroutine set_packed_octTR

!2011/05/16
subroutine connect_new_octs
  use oct_set
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::index,iLv
  integer(kind=4),dimension(27)::MtnNB
  integer(kind=4),dimension(6 )::tNBIndex
  type(oct),pointer::p0

!!$  if(debugMode>=3)then
!!$     do i=0,nprocs-1
!!$        if(rank==i)then
!!$           print *,"MnMin =",MinMn,rank
!!$           print *,"MnMax =",MaxMn,rank
!!$           print *,"MnDisp=",MnDisp,rank
!!$        endif
!!$        call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$     enddo
!!$  endif
  
  !connect new left octs

  do index=MinID(1,-1),MaxID(1,-1)
     !p0=>BMesh(index)
     p0=>Mesh(index)
     !call Morton_neighbours(p0%iPOS,MtnNB,1)
     call get_neighbourIndex(p0%iPOS,MtnNB)
     tNbIndex(1:6)=MtnNB(1:6) - MnDisp

     

     if(MtnNB(1)>=MinMn .and. MtnNB(1)<=MaxMn)then
        p0%octNb1=>Mesh(tNbIndex(1))
     else
        nullify(p0%octNb1)
     endif
     if(MtnNB(4)>=MinMn .and. MtnNB(4)<=MaxMn)then
        p0%octNb2=>Mesh(tNbIndex(4))
     else
        nullify(p0%octNb2)
     endif
     if(MtnNB(2)>=MinMn .and. MtnNB(2)<=MaxMn)then
        p0%octNb3=>Mesh(tNbIndex(2))
     else
        nullify(p0%octNb3)
     endif
     if(MtnNB(5)>=MinMn .and. MtnNB(5)<=MaxMn)then
        p0%octNb4=>Mesh(tNbIndex(5))
     else
        nullify(p0%octNb4)
     endif
     if(MtnNB(3)>=MinMn .and. MtnNB(3)<=MaxMn)then
        p0%octNb5=>Mesh(tNbIndex(3))
     else
        nullify(p0%octNb5)
     endif
     if(MtnNB(6)>=MinMn .and. MtnNB(6)<=MaxMn)then
        p0%octNb6=>Mesh(tNbIndex(6))
     else
        nullify(p0%octNb6)
     endif
  enddo

  !connect children
  !if(debugMode>=3)print *,"connecting children",rank
  do iLv=-1,LvMax
     call connect_oct_with_BMesh(iLv)
  enddo

!!$  if(debugMode>=1)then
!!$     print *,"connect_new_octs completed",rank
!!$  endif

end subroutine connect_new_octs
!2011/05/20

!***********************************************************************
subroutine connect_oct_with_BMesh(iLv)
! +-------------------------------------------------------------------+
! |     create connections for additinal octs                         |
! +-------------------------------------------------------------------+
!***********************************************************************
  use oct_set
  use param
  use const
  use init_mesh_size
  implicit none
  integer(kind=4)    :: iLv,i,index
  type(oct), pointer :: p0
  
  if(iLv.ge.LvMax) return 
  do i=1,2
     !if(i==2.and.iLv==-1)cycle
     if(maxID(i,iLv).gt.minID(i,iLv)) then 
        
!$omp parallel do private(index)
        do index=minID(i,iLv),maxID(i,iLv)                    
           p0 => Mesh(index)
         
           !if(p0%iFLG(2)>=1) then 
           if(p0%iFLG(1)>0)then
              if(p0%iFLG(1)>=2) then ! This means that the cell is not at the edge.
! - for octNb1 -
                 if(associated(p0%octNb1))then
                    p0%octCh1%octNb1 => p0%octNb1%octCh2     
                    p0%octCh3%octNb1 => p0%octNb1%octCh4   
                    p0%octCh5%octNb1 => p0%octNb1%octCh6 
                    p0%octCh7%octNb1 => p0%octNb1%octCh8
                 endif
! - for octNb2 -                
                 if(associated(p0%octNb2))then
                    p0%octCh2%octNb2 => p0%octNb2%octCh1
                    p0%octCh4%octNb2 => p0%octNb2%octCh3
                    p0%octCh6%octNb2 => p0%octNb2%octCh5
                    p0%octCh8%octNb2 => p0%octNb2%octCh7
                 endif
! - for octNb3 -
                 if(associated(p0%octNb3))then
                    p0%octCh1%octNb3 => p0%octNb3%octCh3
                    p0%octCh2%octNb3 => p0%octNb3%octCh4
                    p0%octCh5%octNb3 => p0%octNb3%octCh7
                    p0%octCh6%octNb3 => p0%octNb3%octCh8
                 endif
! - for octNb4 -
                 if(associated(p0%octNb4))then
                    p0%octCh3%octNb4 => p0%octNb4%octCh1
                    p0%octCh4%octNb4 => p0%octNb4%octCh2
                    p0%octCh7%octNb4 => p0%octNb4%octCh5
                    p0%octCh8%octNb4 => p0%octNb4%octCh6
                 endif
! - for octNb5 -
                 if(associated(p0%octNb5))then
                    p0%octCh1%octNb5 => p0%octNb5%octCh5
                    p0%octCh2%octNb5 => p0%octNb5%octCh6
                    p0%octCh3%octNb5 => p0%octNb5%octCh7
                    p0%octCh4%octNb5 => p0%octNb5%octCh8
                 endif
! - for octNb6 -
                 if(associated(p0%octNb6))then
                    p0%octCh5%octNb6 => p0%octNb6%octCh1
                    p0%octCh6%octNb6 => p0%octNb6%octCh2
                    p0%octCh7%octNb6 => p0%octNb6%octCh3
                    p0%octCh8%octNb6 => p0%octNb6%octCh4
                 endif
! - for octCh1 -
                 p0%octCh1%octNb2 => p0       %octCh2
                 p0%octCh1%octNb4 => p0       %octCh3
                 p0%octCh1%octNb6 => p0       %octCh5
! - for octCh2 -
                 p0%octCh2%octNb1 => p0       %octCh1
                 p0%octCh2%octNb4 => p0       %octCh4
                 p0%octCh2%octNb6 => p0       %octCh6
! - for octCh3 -
                 p0%octCh3%octNb2 => p0       %octCh4
                 p0%octCh3%octNb3 => p0       %octCh1
                 p0%octCh3%octNb6 => p0       %octCh7
! - for octCh4 -
                 p0%octCh4%octNb1 => p0       %octCh3
                 p0%octCh4%octNb3 => p0       %octCh2
                 p0%octCh4%octNb6 => p0       %octCh8
! - for octCh5 -
                 p0%octCh5%octNb2 => p0       %octCh6
                 p0%octCh5%octNb4 => p0       %octCh7
                 p0%octCh5%octNb5 => p0       %octCh1
! - for octCh6 -
                 p0%octCh6%octNb1 => p0       %octCh5           
                 p0%octCh6%octNb4 => p0       %octCh8
                 p0%octCh6%octNb5 => p0       %octCh2
! - for octCh7 -
                 p0%octCh7%octNb2 => p0       %octCh8
                 p0%octCh7%octNb3 => p0       %octCh5
                 p0%octCh7%octNb5 => p0       %octCh3
! - for octCh8 -
                 p0%octCh8%octNb1 => p0       %octCh7
                 p0%octCh8%octNb3 => p0       %octCh6
                 p0%octCh8%octNb5 => p0       %octCh4
              endif
              if(p0%iFLG(1)==1) then !This is at the edge of overlap region
! - for octCh1 -
                 p0%octCh1%octNb1 => p0%octCh1 
                 p0%octCh1%octNb2 => p0%octCh2 !
                 p0%octCh1%octNb3 => p0%octCh1 
                 p0%octCh1%octNb4 => p0%octCh3 !
                 p0%octCh1%octNb5 => p0%octCh1 
                 p0%octCh1%octNb6 => p0%octCh5 !
! - for octCh2 -
                 p0%octCh2%octNb1 => p0%octCh1 !
                 p0%octCh2%octNb2 => p0%octCh2
                 p0%octCh2%octNb3 => p0%octCh2
                 p0%octCh2%octNb4 => p0%octCh4 !
                 p0%octCh2%octNb5 => p0%octCh2
                 p0%octCh2%octNb6 => p0%octCh6 !
! - for octCh3 -
                 p0%octCh3%octNb1 => p0%octCh3 
                 p0%octCh3%octNb2 => p0%octCh4 !
                 p0%octCh3%octNb3 => p0%octCh1 !
                 p0%octCh3%octNb4 => p0%octCh3
                 p0%octCh3%octNb5 => p0%octCh3
                 p0%octCh3%octNb6 => p0%octCh7 !
! - for octCh4 -
                 p0%octCh4%octNb1 => p0%octCh3 !
                 p0%octCh4%octNb2 => p0%octCh4
                 p0%octCh4%octNb3 => p0%octCh2 !
                 p0%octCh4%octNb4 => p0%octCh4
                 p0%octCh4%octNb5 => p0%octCh4
                 p0%octCh4%octNb6 => p0%octCh8 !
! - for octCh5 -
                 p0%octCh5%octNb1 => p0%octCh5   
                 p0%octCh5%octNb2 => p0%octCh6 !
                 p0%octCh5%octNb3 => p0%octCh5
                 p0%octCh5%octNb4 => p0%octCh7 !
                 p0%octCh5%octNb5 => p0%octCh1 !
                 p0%octCh5%octNb6 => p0%octCh5
! - for octCh6 -
                 p0%octCh6%octNb1 => p0%octCh5 !
                 p0%octCh6%octNb2 => p0%octCh6
                 p0%octCh6%octNb3 => p0%octCh6
                 p0%octCh6%octNb4 => p0%octCh8 !
                 p0%octCh6%octNb5 => p0%octCh2 !
                 p0%octCh6%octNb6 => p0%octCh6
! - for octCh7 -
                 p0%octCh7%octNb1 => p0%octCh7 
                 p0%octCh7%octNb2 => p0%octCh8 !
                 p0%octCh7%octNb3 => p0%octCh5 !
                 p0%octCh7%octNb4 => p0%octCh7
                 p0%octCh7%octNb5 => p0%octCh3 !
                 p0%octCh7%octNb6 => p0%octCh7
! - for octCh8 -
                 p0%octCh8%octNb1 => p0%octCh7 !
                 p0%octCh8%octNb2 => p0%octCh8
                 p0%octCh8%octNb3 => p0%octCh6 !
                 p0%octCh8%octNb4 => p0%octCh8
                 p0%octCh8%octNb5 => p0%octCh4 !
                 p0%octCh8%octNb6 => p0%octCh8
!
                 if(associated(p0%octNb1))then
                    if(p0%octNb1%iFLG(1) >= 1) then ! In case there are children in octNb1
                       p0%octCh1%octNb1 => p0%octNb1%octCh2
                       p0%octCh3%octNb1 => p0%octNb1%octCh4
                       p0%octCh5%octNb1 => p0%octNb1%octCh6
                       p0%octCh7%octNb1 => p0%octNb1%octCh8
                    endif
                 endif
                 if(associated(p0%octNb2))then
                    if(p0%octNb2%iFLG(1) >= 1) then !
                       p0%octCh2%octNb2 => p0%octNb2%octCh1
                       p0%octCh4%octNb2 => p0%octNb2%octCh3
                       p0%octCh6%octNb2 => p0%octNb2%octCh5
                       p0%octCh8%octNb2 => p0%octNb2%octCh7
                    endif
                 endif
                 if(associated(p0%octNb3))then
                    if(p0%octNb3%iFLG(1) >= 1) then
                       p0%octCh1%octNb3 => p0%octNb3%octCh3
                       p0%octCh2%octNb3 => p0%octNb3%octCh4
                       p0%octCh5%octNb3 => p0%octNb3%octCh7
                       p0%octCh6%octNb3 => p0%octNb3%octCh8
                    endif
                 endif
                 if(associated(p0%octNb4))then
                    if(p0%octNb4%iFLG(1) >= 1) then
                       p0%octCh3%octNb4 => p0%octNb4%octCh1
                       p0%octCh4%octNb4 => p0%octNb4%octCh2
                       p0%octCh7%octNb4 => p0%octNb4%octCh5
                       p0%octCh8%octNb4 => p0%octNb4%octCh6
                    endif
                 endif
                 if(associated(p0%octNb5))then
                    if(p0%octNb5%iFLG(1) >= 1) then
                       p0%octCh1%octNb5 => p0%octNb5%octCh5
                       p0%octCh2%octNb5 => p0%octNb5%octCh6
                       p0%octCh3%octNb5 => p0%octNb5%octCh7
                       p0%octCh4%octNb5 => p0%octNb5%octCh8
                    endif
                 endif
                 if(associated(p0%octNb6))then
                    if(p0%octNb6%iFLG(1) >= 1) then
                       p0%octCh5%octNb6 => p0%octNb6%octCh1
                       p0%octCh6%octNb6 => p0%octNb6%octCh2
                       p0%octCh7%octNb6 => p0%octNb6%octCh3
                       p0%octCh8%octNb6 => p0%octNb6%octCh4
                    endif
                 endif
              endif
           endif
        end do
!$omp end parallel do
     endif
  end do
     
  return
end subroutine connect_oct_with_BMesh
!2011/07/05
!initialize Mesh for sort_set_oct
subroutine reset_Mesh
  use oct_set
  use param
  use init_mesh_size
  implicit none
  integer(kind=4)::index
  type(oct),pointer::p0

  do index=1,maxID(2,LvMax) !for all used octs
     p0=>Mesh(index)
     
     p0%iFLG=-4
     p0%iPOS=-1
     p0%rPOS=-1
     p0%octP=-1
     p0%octN=-1
     p0%MrtN=-1
     

     nullify(p0%octCh1)
     nullify(p0%octCh2)
     nullify(p0%octCh3)
     nullify(p0%octCh4)
     nullify(p0%octCh5)
     nullify(p0%octCh6)
     nullify(p0%octCh7)
     nullify(p0%octCh8)

     nullify(p0%octPrt)

     nullify(p0%octNb1)
     nullify(p0%octNb2)
     nullify(p0%octNb3)
     nullify(p0%octNb4)
     nullify(p0%octNb5)
     nullify(p0%octNb6)

     !don't delete instance of particle
     nullify(p0%ptcl)
     
  enddo

end subroutine reset_Mesh

!2011/06/30
!this subroutine sorts new received octs MinID(2,iLv)--MaxID(2,iLv)
subroutine sort_set_octs
  use oct_set
  use param
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4) :: iLv,ii,ii0,i,ip,indexP,indexQ
  integer(kind=4) :: index,nprt
  integer(kind=4) :: minIF(0:Lvmax),maxIF(0:Lvmax)
  type(oct), pointer :: p0,p1
  type(prtcl), pointer :: pp
!-------------------------------------------

  !if(debugMode>=1)print*,rank, 'entering sort_set_oct'
  
  call reset_sortBuffer

! ----------- special routine for Lv -1 ------------
  do index=MinID(1,-1),MaxID(1,-1)
     p0 => Mesh(index)
     p0%Psort => Mesh(index)
     Mesh2(index)=Mesh(index)
  enddo
! --------------------------------------------------
  !ii=Nall_ini/8
  !ii=0
  ii=BMeshBound
!
  do iLv=0,LvMax 
     ii0 = ii
     do i=1,2
        if(maxID(i,iLv).gt.minID(i,iLv)) then
           do index=MinID(i,iLv),MaxID(i,iLv)
              p0 => Mesh(index)

              if(p0%iFLG(1)>=-3) then
                 ii=ii+1 !count the number of remaining cells
                 p0%Psort => Mesh(ii)
                 Mesh2(ii)=Mesh(index) !Remaining cells are kept in Mesh2
!!$
              end if
           end do
        end if
     end do
!
     if(ii.eq.ii0) then
        minIF(iLv)=ii0
        maxIF(iLv)=ii0
     else

!$omp parallel do private(index,p0,p1,nprt,indexP,indexQ,pp,ip) shared(ii0,ii,iLv,Mesh,Mesh2)

        do index=ii0+1,ii ! This case corresponds to the situation where surviving cells exist.

           p0 => Mesh2(index)
           p0%octN = index


           if(associated(p0%octNb1))p0%octNb1 => p0%octNb1%Psort
           if(associated(p0%octNb2))p0%octNb2 => p0%octNb2%Psort
           if(associated(p0%octNb3))p0%octNb3 => p0%octNb3%Psort
           if(associated(p0%octNb4))p0%octNb4 => p0%octNb4%Psort
           if(associated(p0%octNb5))p0%octNb5 => p0%octNb5%Psort
           if(associated(p0%octNb6))p0%octNb6 => p0%octNb6%Psort

           if(iLv.gt.-1) then
              p1 => p0%octPrt
              nprt = p1%octN

              p1 => Mesh2(nprt)

              if(p0%Csort<=4) then 
                 if(p0%Csort==1) then 
                    p1%octCh1 => p0%Psort
                 else if(p0%Csort==2) then 
                    p1%octCh2 => p0%Psort
                 else if(p0%Csort==3) then 
                    p1%octCh3 => p0%Psort
                 else 
                    p1%octCh4 => p0%Psort
                 end if
              else
                 if(p0%Csort==5) then 
                    p1%octCh5 => p0%Psort
                 else if(p0%Csort==6) then 
                    p1%octCh6 => p0%Psort
                 else if(p0%Csort==7) then 
                    p1%octCh7 => p0%Psort
                 else 
                    p1%octCh8 => p0%Psort
                 end if
              end if
           end if

           !
           if(iLv.lt.LvMax) then ! if the level of the present oct is coarser than the finest level,
              if(p0%iFLG(1).gt.0 ) then 
                 p0%octCh1%octPrt => p0%Psort
                 p0%octCh2%octPrt => p0%Psort
                 p0%octCh3%octPrt => p0%Psort
                 p0%octCh4%octPrt => p0%Psort
                 p0%octCh5%octPrt => p0%Psort
                 p0%octCh6%octPrt => p0%Psort
                 p0%octCh7%octPrt => p0%Psort
                 p0%octCh8%octPrt => p0%Psort
              end if
           end if
           !
           indexP=p0%octP
           indexQ=p0%Psort%octN

           pp => p0%ptcl
           do ip=1,indexP+1
              pp%Ioct=indexQ
              pp => pp%prtnxt
           end do

        end do

!$omp end parallel do 
        minIF(iLv)=ii0+1
        maxIF(iLv)=ii
     end if
  end do

!reset octs
  call reset_Mesh


!$omp parallel do private(index) shared(ii,Mesh,Mesh2)
  do index=1,ii
     Mesh(index)=Mesh2(index)
  end do
!$omp end parallel do 


!!$  if(debugMode>=2)then
!!$     print *,"=======================",rank
!!$     print *,"Before MinMaxID",rank
!!$     do iLv=-1,LvMax
!!$        print "(A,I2,A,I5,I2)","MinID(1,",iLv,")=",MinID(1,iLv),rank
!!$        print "(A,I2,A,I5,I2)","MaxID(1,",iLv,")=",MaxID(1,iLv),rank
!!$        print "(A,I2,A,I5,I2)","MinID(2,",iLv,")=",MinID(2,iLv),rank
!!$        print "(A,I2,A,I5,I2)","MaxID(2,",iLv,")=",MaxID(2,iLv),rank
!!$     enddo
!!$  endif

  do iLv=0,LvMax
     minID(1,iLv)=minIF(iLv)
     maxID(1,iLv)=maxIF(iLv)
  end do

  do iLv=0,LvMax-1
     minID(2,iLv)=maxIF(LvMax)
     maxID(2,iLv)=maxIF(LvMax)
  end do

  minID(2,LvMax)=maxIF(LvMax)
  maxID(2,LvMax)=maxIF(LvMax)

!!$  if(debugMode>=2)then
!!$     print *,"||||||||||||||||",rank
!!$     print *,"VVVVVVVVVVVVVVVV",rank
!!$     print *,"After MinMaxID",rank
!!$     do iLv=-1,LvMax
!!$        print "(A,I2,A,I5,I2)","MinID(1,",iLv,")=",MinID(1,iLv),rank
!!$        print "(A,I2,A,I5,I2)","MaxID(1,",iLv,")=",MaxID(1,iLv),rank
!!$        print "(A,I2,A,I5,I2)","MinID(2,",iLv,")=",MinID(2,iLv),rank
!!$        print "(A,I2,A,I5,I2)","MaxID(2,",iLv,")=",MaxID(2,iLv),rank
!!$     enddo
!!$     print *,"=======================",rank
!!$  endif

 ! if(debugMode>=1)print*, rank, 'exiting sort_set_oct'
  return

end subroutine sort_set_octs

subroutine DDD_sendrecv
  use param
  use message_passing_interface
  implicit none

  if(debugMode>=1)then
     print*,"===========================================================",rank
     print*,"==================== DDD_sendrecv Start ===================",rank
     print*,"===========================================================",rank
  endif

  call get_moved_octs_size 
  call exchange_moved_octs

  call set_moved_octs
  call connect_new_octs

  !sort new octs  
  call sort_set_octs

  if(debugMode>=1)then
     print*,"============================================================",rank
     print*,'==================== DDD_sendrecv  End  ====================',rank
     print*,"============================================================",rank
  endif
end subroutine DDD_sendrecv

subroutine initial_AMR
  use param
  use message_passing_interface
  implicit none
  integer(kind=4)::i,j

 ! if(debugMode>=1)print *,"entering initial_AMR",rank
!create hierarchy grid to get particle loops
  if(LvMax>0)then
     do j=1,LvMax
        do i=j-1,0,-1
           call density(i)
        enddo
        call refine_oct
        call sort_oct
     enddo
  endif
 ! if(debugMode>=1)print *,"exiting inital_AMR",rank
end subroutine initial_AMR

subroutine initial_DDD
  use param
  use message_passing_interface
  implicit none
  integer(kind=4)::i,j

  !if(debugMode>=1)print *,"entering initial_DDD",rank

  call load_balance !get load balance
  call reset_all_mesh !reset Mesh and Pesh

  if(rank==0)print *,"OOOOOOOOOO reset all mesh completed OOOOOOOOOO"

  call initial_setting(1) !remake Mesh and particles
  call setup_model

!remake hierarchy grid
  if(LvMax>0)then
     do j=1,LvMax
        do i=j-1,0,-1
           call density(i)
        enddo
        call refine_oct
        call sort_oct
     enddo
  endif

  !if(debugMode>=1)print *,"exiting inital_DDD",rank
end subroutine initial_DDD

!additional part end
!---------------------------------DDD related Subroutines END ------------------------------

!-------------------------------------------------------------------------------------------
!-------------------------------- Non essential Routines ----------------------------------- 

subroutine nprocs_check
  use mpi
  use message_passing_interface
  implicit none
  if(rank==0)print *,"*******************Checking Num of Process******************"
  if(NXR*NYR*NZR/=nprocs)then
     print *,"Error(nprocs_check):NXR*NYR*NZR isnt equal to nprocs_param",rank,nprocs,NXR*NYR*NZR
     stop
  endif
end subroutine nprocs_check


subroutine findCPU(Mn, CPU)
use message_passing_interface
use init_mesh_size
implicit none

integer(kind=4)::i,sID,eID,midian,loopnum
integer(kind=4),intent(in)::Mn
integer(kind=4),intent(out):: CPU

CPU=-1
if (Mn < 0 ) then 
   print*,'Mn should be >=0', rank, Mn
   stop
endif

if( Mn <=Mn2CPU(1)) then
   CPU = 0
else
   sID=1
   eID=nprocs
   midian=sID+(eID-sID)/2
   loopnum=int(log(dble(nprocs))/log(2.d0))+1
 
   do i=1,loopnum
      if(Mn>=Mn2CPU(midian)+1 .and. Mn<=Mn2CPU(midian+1))then
         CPU=midian
      elseif(Mn>Mn2CPU(midian+1))then
         sID=midian
         eID=eID
         midian=sID+(eID-sID)/2
      elseif(Mn<Mn2CPU(midian)+1)then
         sID=sID
         eID=midian
         midian=sID+(eID-sID)/2
      endif
   enddo
endif

!--Debug BGN
if(CPU > nprocs .or. CPU<0) then
   print*,'something is wrong in findCPU',rank,Mn, CPU
   print *,"Mn2CPU=",Mn2CPU
   stop
endif
!--Debug END

return
end subroutine findCPU

!***********************************************************************************************
!---------------------- subtroutines for Debug
!***********************************************************************************************

subroutine Debug_Display_BMesh(iPOSN)
  Use init_mesh_size
  use param
  use message_passing_interface
  use const
  use oct_set
  implicit none
  integer(kind=4) :: iPOSN(3) ! iPOSN is the Normal ordering, while iPOSH is the hieralchal ordering
  integer(kind=4) :: Mn, index

!!$  iPOSH = 2*LvMax2*iPOSN-LvMax2
!!$  call Morton_number(iPOSH,Mn)
  call get_indexNumberN(iPOSN,Mn)

  if(rank==0) then
     if(Mn >Mn2CPU(1) .or. Mn < 0 ) then
        print*,"This BMesh does not exist in my process: rank,Mn", rank, Mn
        goto 115
     endif
  else 
     if(Mn < 0 .or. Mn >Mn2CPU(rank+1) .or. Mn <= Mn2CPU(rank) ) then
        print*,"This BMesh does not exist in my process: rank, Mn", rank, Mn
        goto 115
     endif
  endif
  
  index=Mn-Mndisp
  
  print '( A,1X, I4, 1X, I4, 1X, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh: index, Mn, iPOS, rPOS', &
index, Mn, Mesh(index) % iPOS, Mesh(index) % rPOS
  if(associated(Mesh(index) % octNb1))then
     print'(A,I4, I4, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Nb1', &
index, Mn, Mesh(index) % octNb1 % iPOS, Mesh(index) % octNb1 % rPOS
  else
     print*,'No Nb1'
  end if
  if(associated(Mesh(index) % octNb2))then
     print'(A,I4, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Nb2',&
 index, Mn, Mesh(index) % octNb2 % iPOS, Mesh(index) % octNb2 % rPOS
  else
     print*,'No Nb2'
  end if
  if(associated(Mesh(index) % octNb3))then
     print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Nb3', &
index, Mn, Mesh(index) % octNb3 % iPOS, Mesh(index) % octNb3 % rPOS
  else
     print*,'No Nb3'
  end if
  if(associated(Mesh(index) % octNb4))then
     print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Nb4', &
index, Mn, Mesh(index) % octNb4 % iPOS, Mesh(index) % octNb4 % rPOS
  else
     print*,'No Nb4'
  end if
  if(associated(Mesh(index) % octNb5))then
     print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Nb5', &
index, Mn, Mesh(index) % octNb5 % iPOS, Mesh(index) % octNb5 % rPOS
  else
     print*,'No Nb5'
  end if
  if(associated(Mesh(index) % octNb6))then
     print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Nb6', &
index, Mn, Mesh(index) % octNb6 % iPOS, Mesh(index) % octNb6 % rPOS
  else
     print*,'No Nb6'
  end if
  print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Ch1', index, Mn, Mesh(index) % octCh1 % iPOS, Mesh(index) % octCh1 % rPOS
  print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Ch2', index, Mn, Mesh(index) % octCh2 % iPOS, Mesh(index) % octCh2 % rPOS
  print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Ch3', index, Mn, Mesh(index) % octCh3 % iPOS, Mesh(index) % octCh3 % rPOS
  print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Ch4', index, Mn, Mesh(index) % octCh4 % iPOS, Mesh(index) % octCh4 % rPOS
  print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Ch5', index, Mn, Mesh(index) % octCh5 % iPOS, Mesh(index) % octCh5 % rPOS
  print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Ch6', index, Mn, Mesh(index) % octCh6 % iPOS, Mesh(index) % octCh6 % rPOS
  print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Ch7', index, Mn, Mesh(index) % octCh7 % iPOS, Mesh(index) % octCh7 % rPOS
  print'(A,I3, I3, 3(1X, I4),1X, 3(1X, F7.4))', 'BMesh Ch8', index, Mn, Mesh(index) % octCh8 % iPOS, Mesh(index) % octCh8 % rPOS
        
!!$print*, 'Nb1', p0 % octNb1 % rPOS, 'Nb2', p0 % octNb2 % rPOS, 'NB3', p0%octNb3% rPOS
!!$ print*, 'rank', rank, 'NB4:' ,p0 % octNb4 % rPOS,'Nb5', p0 % octNb5 % rPOS, 'NB6', p0%octNb6% rPOS
115 continue 
  return
end subroutine Debug_Display_BMesh

subroutine Debug_one_oct_connection(index,nowLv,MeshState,errFlag)
  use param
  use message_passing_interface
  use init_mesh_size
  implicit none
  integer(kind=4),intent(in)::index,nowLv,MeshState
  integer(kind=4)::errFlag,errConn
  integer(kind=4)::iLv, iPOSH(3), iPOSN(3)
  integer(kind=4)::iPOSH2(3), iPOSN2(3)

  type(oct),pointer::p0

!print*,'*************Mesh connection check start*************', rank

  p0 => Mesh(index)

  iLv=nowLv
  iPOSH = p0%iPOS
  iPOSN = (iPOSH + 2**(LvMax-iLv))/(2**(LvMax+1-iLv))
!--octNb connection Check (Oct only)
  if(p0%iFLG(1)>-3 .and. .not.(p0%octType==1 .and. iLv==-1) )then
     if(.not.  associated(p0%octNb1))then
        if(MeshState==1)then
           print*,'No octNb1 id=', index,"iPosN=",iPOSN,"Lv=",p0%octLv,"Type=",p0%octType,"FLG=",p0%iFLG,rank
           if(iLv>=0)print*,"PFLG=",p0%octPrt%iFLG,rank
           errFlag=1
        endif
     endif
     if(.not.  associated(p0%octNb2))then
        if(MeshState==1)then
           print*,'No octNb2 id=', index,"iPosN=",iPOSN,"Lv=",p0%octLv,"Type=",p0%octType,"FLG=",p0%iFLG,rank
           if(iLv>=0)print*,"PFLG=",p0%octPrt%iFLG,rank
           errFlag=1
        endif
     end if
     if(.not.  associated(p0%octNb3))then
        if(MeshState==1)then
           print*,'No octNb3 id=', index,"iPosN=",iPOSN,"Lv=",p0%octLv,"Type=",p0%octType,"FLG=",p0%iFLG,rank
           if(iLv>=0)print*,"PFLG=",p0%octPrt%iFLG,rank
           errFlag=1
        endif
     endif
     if(.not.  associated(p0%octNb4))then
        if(MeshState==1)then
           print*,'No octNb4 id=', index,"iPosN=",iPOSN,"Lv=",p0%octLv,"Type=",p0%octType,"FLG=",p0%iFLG,rank
           if(iLv>=0)print*,"PFLG=",p0%octPrt%iFLG,rank
           errFlag=1
        endif
     endif
     if(.not.  associated(p0%octNb5))then
        if(MeshState==1)then
           print*,'No octNb5 id=', index,"iPosN=",iPOSN,"Lv=",p0%octLv,"Type=",p0%octType,"FLG=",p0%iFLG,rank
           if(iLv>=0)print*,"PFLG=",p0%octPrt%iFLG,rank
           errFlag=1
        endif
     endif
     if(.not.  associated(p0%octNb6))then
        if(MeshState==1)then
           print*,'No octNb6 id=', index,"iPosN=",iPOSN,"Lv=",p0%octLv,"Type=",p0%octType,"FLG=",p0%iFLG,rank
           if(iLv>=0)print*,"PFLG=",p0%octPrt%iFLG,rank
           errFlag=1
        endif
     endif

  endif

!-----octLv Check   
  if(associated(p0%octNb1))then
     if(p0%octLv /= p0%octNb1% octLv)then
        print*,'octNb1%octLv Descrepancy ',"index=",index,"octN=",p0%octN,"iPos=",p0%iPOS,rank
        errFlag=1
        !stop
     endif
  endif
  if(associated(p0%octNb2))then
     if(p0%octLv /= p0%octNb2% octLv)then
        print*,'octNb2%octLv Descrepancy ',"index=",index,"octN=",p0%octN,"iPos=",p0%iPOS,rank
        errFlag=1
        !stop
     endif
  endif
  if(associated(p0%octNb3))then
     if(p0%octLv /= p0%octNb3% octLv)then
        print*,'octNb3%octLv Descrepancy ',"index=",index,"octN=",p0%octN,"iPos=",p0%iPOS,rank
        errFlag=1
        !stop
     endif
  endif
  if(associated(p0%octNb4))then
     if(p0%octLv /= p0%octNb4% octLv)then
        print*,'octNb4%octLv Descrepancy ',"index=",index,"octN=",p0%octN,"iPos=",p0%iPOS,rank
        errFlag=1
        !stop
     endif
  endif
  if(associated(p0%octNb5))then
     if(p0%octLv /= p0%octNb5% octLv)then
        print*,'octNb5%octLv Descrepancy ',"index=",index,"octN=",p0%octN,"iPos=",p0%iPOS,rank
        errFlag=1
        !stop
     endif
  endif
  if(associated(p0%octNb6))then
     if(p0%octLv /= p0%octNb6% octLv)then
        print*,'octNb6%octLv Descrepancy ',"index=",index,"octN=",p0%octN,"iPos=",p0%iPOS,rank
        errFlag=1
        !stop
     endif
  endif

!--Location Check 
!  if(p0%octType==0)then
     if(associated(p0%octNb1))then
        iPOSH2 = p0%octNb1%iPOS
        iPOSN2 = (iPOSH2 + 2**(LvMax-iLv))/(2**(LvMax+1-iLv))
        if(iLv==-1)then
           if(iPOSN2(1) == NX_2t)iPOSN2(1)=0
        else
           if(iPOSN2(1) == NXR*NXB*2**(iLv))iPOSN2(1) =0
        endif
        errConn=0
        if(iLv>0 .and. p0%octType/= 1        .and. iPOSN2(1)+1/=iPOSN(1) .and. p0%iFLG(1)/=-3       )errConn=1
        if(iLv==0.and. p0%octType/= 1        .and. iPOSN2(1)+1/=iPOSN(1)                            )errConn=1
        if(iLv>=0.and. p0%octType== 1        .and. iPOSN2(1)+1/=iPOSN(1) .and. iPOSN2(1)  /=iPOSN(1))errConn=1
        if(iLv<0 .and.                             iPOSN2(1)+1/=iPOSN(1)                            )errConn=1

        if(errConn/=0)then
           
           print*,'octNb1 iPOS error index=', index," iLv=",iLv,rank
           if(iLv/=-1)print*,"iFLG=",p0%iFLG,"prtiFLG=",p0%octPrt%iFLG,rank
           print*,"p0%iPos=",iPOSN,"NbiPos=", iPOSN2,rank
           print*,"octNb1 octN=",p0%octNb1%octN,"octType=",p0%octType,"Nb1 octType=",p0%octNb1%octType,rank
           if(iLv/=-1)print*,"octPrt=",p0%octPrt%octN,"Nb1 octPrt=",p0%octNb1%octPrt%octN,"Csort=",p0%Csort,"Nb1 Csort=",p0%octNb1%Csort
           print*,"---------------------------------------------------------------------",rank
           errFlag=1   
        endif
     endif
     if(associated(p0%octNb2))then
        iPOSH2 = p0%octNb2%iPOS
        iPOSN2 = (iPOSH2 + 2**(LvMax-iLv))/(2**(LvMax+1-iLv))
        if(iLv==-1)then
           if(iPOSN2(1) == 1)iPOSN2(1)=NX_2t+1
        else
           if(iPOSN2(1) == 1)iPOSN2(1) =NXR*NXB*2**(iLv)+1
        endif
        errConn=0
        if(iLv>0 .and. p0%octType/= 1        .and. iPOSN2(1)-1/=iPOSN(1) .and. p0%iFLG(1)/=-3       )errConn=1
        if(iLv==0.and. p0%octType/= 1        .and. iPOSN2(1)-1/=iPOSN(1)                            )errConn=1
        if(iLv>=0.and. p0%octType== 1        .and. iPOSN2(1)-1/=iPOSN(1) .and. iPOSN2(1)  /=iPOSN(1))errConn=1
        if(iLv<0 .and.                             iPOSN2(1)-1/=iPOSN(1)                            )errConn=1

        if(errConn/=0)then
           print*,'octNb2 iPOS error index=', index," iLv=",iLv,rank
           if(iLv/=-1)print*,"iFLG=",p0%iFLG,"prtiFLG=",p0%octPrt%iFLG,rank
           print*,"p0%iPos=",iPOSN,"NbiPos=", iPOSN2,rank
           print*,"octNb2 octN=",p0%octNb2%octN,"octType=",p0%octType,"Nb2 octType=",p0%octNb2%octType,rank
           if(iLv/=-1)print*,"octPrt=",p0%octPrt%octN,"Nb2 octPrt=",p0%octNb2%octPrt%octN,"Csort=",p0%Csort,"Nb2 Csort=",p0%octNb2%Csort
           print*,"---------------------------------------------------------------------",rank
           errFlag=1   
        endif
     endif
     if(associated(p0%octNb3))then
        iPOSH2 = p0%octNb3%iPOS
        iPOSN2 = (iPOSH2 + 2**(LvMax-iLv))/(2**(LvMax+1-iLv))
        if(iLv==-1)then
           if(iPOSN2(2) == NY_2t)iPOSN2(2)=0
        else
           if(iPOSN2(2) == NYR*NYB*2**(iLv))iPOSN2(2) =0
        endif
        errConn=0

        if(iLv>0 .and. p0%octType/= 1        .and. iPOSN2(2)+1/=iPOSN(2) .and. p0%iFLG(1)/=-3       )errConn=1
        if(iLv==0.and. p0%octType/= 1        .and. iPOSN2(2)+1/=iPOSN(2)                            )errConn=1
        if(iLv>=0.and. p0%octType== 1        .and. iPOSN2(2)+1/=iPOSN(2) .and. iPOSN2(2)  /=iPOSN(2))errConn=1
        if(iLv<0 .and.                             iPOSN2(2)+1/=iPOSN(2)                            )errConn=1

        if(errConn/=0)then
           print*,'octNb3 iPOS error index=', index," iLv=",iLv,rank
           if(iLv/=-1)print*,"iFLG=",p0%iFLG,"prtiFLG=",p0%octPrt%iFLG,rank
           print*,"p0%iPos=",iPOSN,"NbiPos=", iPOSN2,rank
           print*,"octNb3 octN=",p0%octNb3%octN,"octType=",p0%octType,"Nb3 octType=",p0%octNb3%octType,rank
           if(iLv/=-1)print*,"octPrt=",p0%octPrt%octN,"Nb3 octPrt=",p0%octNb3%octPrt%octN,"Csort=",p0%Csort,"Nb3 Csort=",p0%octNb3%Csort
           print*,"---------------------------------------------------------------------",rank
           errFlag=1 
        endif
     endif
     if(associated(p0%octNb4))then
        iPOSH2 = p0%octNb4%iPOS
        iPOSN2 = (iPOSH2 + 2**(LvMax-iLv))/(2**(LvMax+1-iLv))

        if(iLv==-1)then
           if(iPOSN2(2) == 1)iPOSN2(2)=NY_2t+1
        else
           if(iPOSN2(2) == 1)iPOSN2(2) =NYR*NYB*2**(iLv)+1
        endif
        errConn=0

        if(iLv>0 .and. p0%octType/= 1        .and. iPOSN2(2)-1/=iPOSN(2) .and. p0%iFLG(1)/=-3       )errConn=1
        if(iLv==0.and. p0%octType/= 1        .and. iPOSN2(2)-1/=iPOSN(2)                            )errConn=1
        if(iLv>=0.and. p0%octType== 1        .and. iPOSN2(2)-1/=iPOSN(2) .and. iPOSN2(2)  /=iPOSN(2))errConn=1
        if(iLv<0 .and.                             iPOSN2(2)-1/=iPOSN(2)                            )errConn=1

        if(errConn/=0)then
           print*,'octNb4 iPOS error index=', index," iLv=",iLv,rank
           if(iLv/=-1)print*,"iFLG=",p0%iFLG,"prtiFLG=",p0%octPrt%iFLG,rank
           print*,"p0%iPos=",iPOSN,"NbiPos=", iPOSN2,rank
           print*,"octNb4 octN=",p0%octNb4%octN,"octType=",p0%octType,"Nb4 octType=",p0%octNb4%octType,rank
           if(iLv/=-1)print*,"octPrt=",p0%octPrt%octN,"Nb4 octPrt=",p0%octNb4%octPrt%octN,"Csort=",p0%Csort,"Nb4 Csort=",p0%octNb4%Csort
           print*,"---------------------------------------------------------------------",rank
           errFlag=1
        endif
     endif
     if(associated(p0%octNb5))then
        iPOSH2 = p0%octNb5%iPOS
        iPOSN2 = (iPOSH2 + 2**(LvMax-iLv))/(2**(LvMax+1-iLv))

        if(iLv==-1)then
           if(iPOSN2(3) == NZ_2t)iPOSN2(3)=0
        else
           if(iPOSN2(3) == NZR*NZB*2**(iLv))iPOSN2(3) =0
        endif
        errConn=0

        if(iLv>0 .and. p0%octType/= 1        .and. iPOSN2(3)+1/=iPOSN(3) .and. p0%iFLG(1)/=-3       )errConn=1
        if(iLv==0.and. p0%octType/= 1        .and. iPOSN2(3)+1/=iPOSN(3)                            )errConn=1
        if(iLv>=0.and. p0%octType== 1        .and. iPOSN2(3)+1/=iPOSN(3) .and. iPOSN2(3)  /=iPOSN(3))errConn=1
        if(iLv<0 .and.                             iPOSN2(3)+1/=iPOSN(3)                            )errConn=1


        if(errConn/=0)then
           print*,'octNb5 iPOS error index=', index," iLv=",iLv,rank
           if(iLv/=-1)print*,"iFLG=",p0%iFLG,"prtiFLG=",p0%octPrt%iFLG,rank
           print*,"p0%iPos=",iPOSN,"NbiPos=", iPOSN2,rank
           print*,"octNb5 octN=",p0%octNb5%octN,"octType=",p0%octType,"Nb5 octType=",p0%octNb5%octType,rank
           if(iLv/=-1)print*,"octPrt=",p0%octPrt%octN,"Nb5 octPrt=",p0%octNb5%octPrt%octN,"Csort=",p0%Csort,"Nb5 Csort=",p0%octNb5%Csort
           print*,"---------------------------------------------------------------------",rank
           errFlag=1
        endif
     endif
     if(associated(p0%octNb6))then
        iPOSH2 = p0%octNb6%iPOS
        iPOSN2 = (iPOSH2 + 2**(LvMax-iLv))/(2**(LvMax+1-iLv))

        if(iLv==-1)then
           if(iPOSN2(3) == 1)iPOSN2(3)=NZ_2t+1
        else
           if(iPOSN2(3) == 1)iPOSN2(3) =NZR*NZB*2**(iLv)+1
        endif
        errConn=0
        if(iLv>0 .and. p0%octType/= 1        .and. iPOSN2(3)-1/=iPOSN(3) .and. p0%iFLG(1)/=-3       )errConn=1
        if(iLv==0.and. p0%octType/= 1        .and. iPOSN2(3)-1/=iPOSN(3)                            )errConn=1
        if(iLv>=0.and. p0%octType== 1        .and. iPOSN2(3)-1/=iPOSN(3) .and. iPOSN2(3)  /=iPOSN(3))errConn=1
        if(iLv<0 .and.                             iPOSN2(3)-1/=iPOSN(3)                            )errConn=1

        if(errConn/=0)then
           print*,'octNb6 iPOS error index=', index," iLv=",iLv,rank
           if(iLv/=-1)print*,"iFLG=",p0%iFLG,"prtiFLG=",p0%octPrt%iFLG,rank
           print*,"p0%iPos=",iPOSN,"NbiPos=", iPOSN2,rank
           print*,"octNb6 octN=",p0%octNb6%octN,"octType=",p0%octType,"Nb6 octType=",p0%octNb6%octType,rank
           if(iLv/=-1)print*,"octPrt=",p0%octPrt%octN,"Nb6 octPrt=",p0%octNb6%octPrt%octN,"Csort=",p0%Csort,"Nb6 Csort=",p0%octNb6%Csort
           print*,"---------------------------------------------------------------------",rank
           errFlag=1
        endif
     endif
!  endif

return
end subroutine Debug_one_oct_connection

!2011/4/25
!for debug
subroutine change_Mn2CPU
  use param
  use message_passing_interface
  implicit none
  integer(kind=4)::i

  print *,"before Mn2CPU=",Mn2CPU
  do i=1,nprocs
     if(i==1)then
        Mn2CPU(i)=Mn2CPU(i)+1
   !     Mn2CPU(i+1)=Mn2CPU(i+1)+1
!!$!!     elseif(i/=nprocs-1)then
!!$!!        Mn2CPU(i)=Mn2CPU(i)+4
     endif
  enddo
  print *,"after Mn2CPU=",Mn2CPU
  if(rank==0)then
     MinMn=0
  else
     MinMn=Mn2CPU(rank)+1
  endif
  MaxMn=Mn2CPU(rank+1)

end subroutine change_Mn2CPU
!2012/4/25
!2011/6/22
!2011/4/21
!debug subroutine
subroutine output_octs_and_family(sID,eID,outputLvMax,MeshType,MeshState,OutType)
  use param
  use message_passing_interface
  implicit none
  integer i
  integer,intent(in)::sID,eID,outputLvMax
  integer,intent(in)::MeshType,MeshState,OutType
  integer(kind=4)::errFlag
  character(LEN=128)::filename

  !MeshType=0 target Oct
  !MeshType=1 target GOct
  !MeshState=0 no GMesh
  !MeshState=1 Mesh and GMesh are exist
  !OutType=0 print error message only
  !OutType=1 full print 
  
  errFlag=0

  print *,"Debug octs ...",rank


  if(MeshType==1)then
    if(sID<MinID(3,-1))then
        print *,"Print out index err:",minID(3,-1),"-",maxID(3,-1)," rank:",rank
        print *,"sID-eID:",sID,eID
        stop
     endif
  elseif(MeshType==0)then
     if(eID>MaxID(1,-1))then
        print *,"Print out index err:",minID(1,-1),"-",maxID(1,-1)," rank:",rank
        print *,"sID-eID:",sID,eID
        stop
     endif
  else
     print *,"MeshType error",MeshType,rank
     stop
  endif

  if(MeshState/=0 .and. MeshState/=1)then
     print *,"MeshState error",MeshState,rank
  endif

  !for output
  write(filename,2311)version,rank
!  open(98,file=filename,form="formatted")
!  close(98)
  do i=sID,eID
     call output_oct_and_family(i,outputLvMax,-1,MeshState,OutType,errFlag)
  enddo

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(errFlag/=0)then
     print *,"error is detected istep=",istep,rank
     stop
  endif

  print *,"... Debug complete",rank
 
  2311 format("./output/",A,"_DOP",i3.3,".dat")
end subroutine output_octs_and_family

subroutine full_oct_check
  use param
  use message_passing_interface
  implicit none
  integer(kind=4)::errFlag,i
  character(LEN=128)::filename
  errFlag=0
  if(debugMode>0)print *,"full oct check start",rank
  write(filename,2311)version,rank
!  open(98,file=filename,form="formatted")
!  close(98)
  do i=MinID(1,-1),MaxID(1,-1)
     call output_oct_and_family(i,LvMax,-1,1,0,errFlag)
  enddo
  do i=MinID(3,-1),MaxID(3,-1)
     call output_oct_and_family(i,LvMax,-1,1,0,errFlag)
  enddo
  if(MinID(2,-1)<MaxID(2,-1))then
     do i=MinID(2,-1),MaxID(2,-1)
        call output_oct_and_family(i,LvMax,-1,1,0,errFlag)
     enddo
  endif
  if(MinID(4,-1)<MaxID(4,-1))then
     do i=MinID(4,-1),MaxID(4,-1)
        call output_oct_and_family(i,LvMax,-1,1,0,errFlag)
     enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(errFlag/=0)then
     print *,"error is detected istep=",istep,rank
     stop
  endif

  if(debugMode>0)print *,"full oct check complete",rank
  2311 format("./output/",A,"_DOP",i3.3,".dat")
end subroutine full_oct_check

recursive subroutine output_oct_and_family(index,LvFlag,nowLv,MeshState,OutType,errFlag)
  use init_mesh_size
  use oct_set
  use message_passing_interface
  implicit none

  integer i
  integer(kind=4),intent(in)::index,LvFlag,nowLv
  integer(kind=4),intent(in)::MeshState,OutType
  integer(kind=4)::errFlag,errFlag_cnct,errFlag_oct
  integer(kind=4)::defiPos(3),Csort
  type(oct),pointer::p0
  type(prtcl),pointer::pp
  character(LEN=128)::filename

  !MeshType=0 target Oct
  !MeshType=1 target GOct
  !MeshState=0 no GMesh
  !MeshState=1 Mesh and GMesh are exist
  !OutType=0 print error message only
  !OutType=1 full print 

  nullify(pp)

  errFlag_cnct=0
  errFlag_oct=0

  p0 => Mesh(index)
  if(p0%iFLG(1)<-3)return

!  defiPos(1)=(p0%iPos(1)-sConst)/intNxt+1
!  defiPos(2)=(p0%iPos(2)-sConst)/intNxt+1
!  defiPos(3)=(p0%iPos(3)-sConst)/intNxt+1

  defiPos = (p0%iPos + sConst(p0%octLv))/intNxt(p0%octLv)

  if(outType==1)then
     print *,""
     print *,"********************************************************"
     print '(A,I3,A,I7,A,I3,A,I3)',"  proc:",rank," index:",index," Lv:",p0%OctLv," octType:",p0%octType
     print *,"--------------------------------------------------------"
     print '(A,I3,A,I3,A,I3,A)',"iPos(",p0%iPos(1),",",p0%iPos(2),",",p0%iPos(3),")"
     print '(A,F8.5,A,F8.5,A,F8.5,A)',"rPos(",p0%rPos(1),",",p0%rPos(2),",",p0%rPos(3),")"
     print '(A,I3,A,I3,A,I3,A)',"Def_iPos(",defiPos(1),",",defiPos(2),",",defiPos(3),")"
     print '(A,F8.5,A,F8.5,A,F8.5,A)',"F1-3(",p0%F(1),",",p0%F(2),",",p0%F(3),")"
     print '(A,I2,A,I4,A,I4)',"Level:",p0%octLv," OctN:",p0%octN," MrtN:",p0%MrtN
     print *,"iFLG=",p0%iFLG
     print *,"octP=",p0%octP
     print *,"********************************************************"
  endif
  !check iFLG
  if(p0%iFLG(1)<=-4)then
     print *,"there is oct of iFLG(1)=-4 between MinID and MaxID",rank
     print *,"index=",index,"iFLG=",p0%iFLG,"rank=",rank
     errFlag=1
     errFlag_oct=1
  endif

  !check octLv and iLv
  if(p0%octLv/=nowLv)then
     print *,"Error! octLv /= nowLv [check MinMaxID] octLv:",p0%octLv,"nowLv=",nowLv,rank
     errFlag=1
     errFlag_oct=1
  endif

  !check representative particle
  if(p0%octLv>=0)then
     if(.not. associated(p0%ptcl))then
        print *,"Error! this oct doesn't have representative particle.",rank
        print *,"octN=",p0%octN,"octLv=",p0%octLv,rank
        errFlag=1
        errFlag_oct=1
        !stop
     else
        !check representative particle
        if(p0%ptcl%isort>=0)then
           print *,"representative particle isort error octN=",p0%octN,"octLv=",p0%octLv,rank
           print *,"octType=",p0%octType,"rep ptcl isort=",p0%ptcl%isort,rank
           errFlag=1
           errFlag_oct=1
        endif

        if(p0%ptcl%ioct/=p0%octN)then
           print *,"representative particle index error octN=",p0%octN,rank
           print *,"octType=",p0%octType,"rep ptcl ioct=",p0%ptcl%ioct,rank
           errFlag=1
           errFlag_oct=1
        endif
     endif
  endif
  !check particle connection
  if(p0%iFLG(1)>-4 .and. p0%iFLG(1)<=0)then
     !check particle linear list
     if(p0%octP>0)then
        i=0
        if(associated(p0%ptcl%prtnxt))then
           pp=>p0%ptcl%prtnxt
        else
           print *,"Error! num of particle =",i," octP=",p0%octP,rank
           errFlag=1
           errFlag_oct=1
           !stop
        endif
        do i=1,p0%octP
           if(pp%ioct/=p0%octN)then
              print *,"Error! ioct=",pp%ioct,"inum=",pp%inum," octN=",p0%octN,"isort=",pp%isort,rank
              print *,"index=",index,"iPos=",p0%iPos,"octType=",p0%octType,"octLv=",p0%octLv,rank
          
              !stop
              errFlag=1
              errFlag_oct=1
           endif
           if(associated(pp%prtnxt))then
              pp=>pp%prtnxt
           else
              if(i/=p0%octP)then
                 print *,"Error! num of particle =",i," octP=",p0%octP,rank
                 !stop
                 errFlag=1
                 errFlag_oct=1
              endif
           endif
        enddo
     endif
  endif
  
  !check index and octN
  if(index/=p0%octN)then
     print *,"index/=octN",index,p0%octN,rank
     errFlag=1
     errFlag_oct=1
     !stop
  endif

  !check neighbor connection
  if(outType==1)then
     print *,"******* Connection *******"

 

     if(associated(p0%octNb1))then 
        print '(A,3I5,A,$)',"[",p0%octNb1%iPos,"   ]"
     else
        print '(A,$)',"[    Null    ]"
     endif
     if(associated(p0%octNb2))then 
        print '(A,3I5,A,$)',"[",p0%octNb2%iPos,"   ]"
     else
        print '(A,$)',"[    Null    ]"
     endif
     if(associated(p0%octNb3))then 
        print '(A,3I5,A,$)',"[",p0%octNb3%iPos,"   ]"
     else
        print '(A,$)',"[    Null    ]"
     endif
     if(associated(p0%octNb4))then 
        print '(A,3I5,A,$)',"[",p0%octNb4%iPos,"   ]"
     else
        print '(A,$)',"[    Null    ]"
     endif
     if(associated(p0%octNb5))then
        print '(A,3I5,A,$)',"[",p0%octNb5%iPos,"   ]"
     else
        print '(A,$)',"[    Null    ]"
     endif
     if(associated(p0%octNb6))then
        print '(A,3I5,A,$)',"[",p0%octNb6%iPos,"   ]"
     else
        print '(A,$)',"[    Null    ]"
     endif

     print *," "

     if(associated(p0%octNb1))then 
        print '(A,$)',"[     -x     ]"
     else
        print '(A,$)',"[    Null    ]"
     endif
     if(associated(p0%octNb2))then 
        print '(A,$)',"[     +x     ]"
     else
        print '(A,$)',"[    Null    ]"
     endif
     if(associated(p0%octNb3))then 
        print '(A,$)',"[     -y     ]"
     else
        print '(A,$)',"[    Null    ]"
     endif
     if(associated(p0%octNb4))then 
        print '(A,$)',"[     +y     ]"
     else
        print '(A,$)',"[    Null    ]"
     endif
     if(associated(p0%octNb5))then
        print '(A,$)',"[     -z     ]"
     else
        print '(A,$)',"[    Null    ]"
     endif
     if(associated(p0%octNb6))then
        print '(A,$)',"[     +z     ]"
     else
        print '(A,$)',"[    Null    ]"
     endif

     print *," " 
  endif

  !oct Neighbor check (Oct only)
  if(MeshState==1)then
     !if(p0%octLv==-1 .and. rank==0)print *,"Lv-1 checking"
     call Debug_one_oct_connection(index,nowLv,MeshState,errFlag_cnct)
     if(errFlag_cnct/=0)errFlag_oct=1
  endif

  !check parent connection
  if(p0%octLv>=0)then
     if(associated(p0%octPrt))then
        if(outType==1)then
           print *,"******* Parent *******"
           print *,"OctPrt -> ",p0%octPrt%octN
        endif

        !check parent iFLG
        if(p0%octPrt%iFLG(1)<=0)then
           print *,"Parent iFLG error PrtiFLG=",p0%octPrt%iFLG,"iFLG=",p0%iFLG,rank
           print *,"iLv=",p0%octLv,rank
           errFlag=1
           errFlag_oct=1
        endif

        if(p0%octLv>0 .and. p0%iFLG(3)/=0 .and. p0%octPrt%iFLG(3)==0)then
           print *,"Parent iFLG3 error PrtiFLG=",p0%octPrt%iFLG,"iFLG=",p0%iFLG,rank
           print *,"iLv=",p0%octLv,rank
           errFlag=1
           errFlag_oct=1
        endif
     
     else
        print *,"Error this oct has no parent."
        !stop
        errFlag=1
        errFlag_oct=1
     endif

     !parent%MrtN check
     if(p0%octPrt%Mrtn/=p0%MrtN)then
        print *,"Parent's MrtN is not equal to this octN=",p0%octN,rank
        print *,"Parent MrtN=",p0%octPrt%octN,"child MrtN=",p0%MrtN,rank
        errFlag=1
        errFlag_oct=1
     endif

     !child->parent->child check
     Csort=p0%Csort
     if(Csort==1)then
        if(p0%octPrt%octCh1%octN/=index)then
           print *,"Parent is not connect to this oct octN=",p0%octN,rank
           print *,"Parent octN=",p0%octPrt%octN,"Ch1 octN",p0%octPrt%octCh1%octN,"index=",index,rank
           errFlag=1
           errFlag_oct=1
        endif
     endif
     if(Csort==2)then
        if(p0%octPrt%octCh2%octN/=index)then
           print *,"Parent is not connect to this oct octN=",p0%octN,rank
           print *,"Parent octN=",p0%octPrt%octN,"Ch1 octN",p0%octPrt%octCh2%octN,"index=",index,rank
           errFlag=1
           errFlag_oct=1
        endif
     endif
     if(Csort==3)then
        if(p0%octPrt%octCh3%octN/=index)then
           print *,"Parent is not connect to this oct octN=",p0%octN,rank
           print *,"Parent octN=",p0%octPrt%octN,"Ch1 octN",p0%octPrt%octCh3%octN,"index=",index,rank
           errFlag=1
           errFlag_oct=1
        endif
     endif
     if(Csort==4)then
        if(p0%octPrt%octCh4%octN/=index)then
           print *,"Parent is not connect to this oct octN=",p0%octN,rank
           print *,"Parent octN=",p0%octPrt%octN,"Ch1 octN",p0%octPrt%octCh4%octN,"index=",index,rank
           errFlag=1
           errFlag_oct=1
        endif
     endif
     if(Csort==5)then
        if(p0%octPrt%octCh5%octN/=index)then
           print *,"Parent is not connect to this oct octN=",p0%octN,rank
           print *,"Parent octN=",p0%octPrt%octN,"Ch1 octN",p0%octPrt%octCh5%octN,"index=",index,rank
           errFlag=1
           errFlag_oct=1
        endif
     endif
     if(Csort==6)then
        if(p0%octPrt%octCh6%octN/=index)then
           print *,"Parent is not connect to this oct octN=",p0%octN,rank
           print *,"Parent octN=",p0%octPrt%octN,"Ch1 octN",p0%octPrt%octCh6%octN,"index=",index,rank
           errFlag=1
           errFlag_oct=1
        endif
     endif
     if(Csort==7)then
        if(p0%octPrt%octCh7%octN/=index)then
           print *,"Parent is not connect to this oct octN=",p0%octN,rank
           print *,"Parent octN=",p0%octPrt%octN,"Ch1 octN",p0%octPrt%octCh7%octN,"index=",index,rank
           errFlag=1
           errFlag_oct=1
        endif
     endif
     if(Csort==8)then
        if(p0%octPrt%octCh8%octN/=index)then
           print *,"Parent is not connect to this oct octN=",p0%octN,rank
           print *,"Parent octN=",p0%octPrt%octN,"Ch1 octN",p0%octPrt%octCh8%octN,"index=",index,rank
           errFlag=1
           errFlag_oct=1
        endif
     endif
  endif


  !check children connection
  if(p0%iFLG(1)>0)then
     if(outType==1)print *,"******* Children *******"
     if(associated(p0%octCh1) .and. associated(p0%octCh2) .and. &
        associated(p0%octCh3) .and. associated(p0%octCh4) .and. &
        associated(p0%octCh5) .and. associated(p0%octCh6) .and. &
        associated(p0%octCh7) .and. associated(p0%octCh8) )then

        if(outType==1)then
           print '(A,I4,A,3I3,A,I3,A,I2)',"Ch1->",p0%octCh1%octN," iPos[",p0%octCh1%iPos,"] octP=",&
                p0%octCh1%octP,", Lv=",p0%octCh1%octLv 
           print '(A,I4,A,3I3,A,I3,A,I2)',"Ch2->",p0%octCh2%octN," iPos[",p0%octCh2%iPos,"] octP=",&
                p0%octCh2%octP,", Lv=",p0%octCh2%octLv 
           print '(A,I4,A,3I3,A,I3,A,I2)',"Ch3->",p0%octCh3%octN," iPos[",p0%octCh3%iPos,"] octP=",&
                p0%octCh3%octP,", Lv=",p0%octCh3%octLv 
           print '(A,I4,A,3I3,A,I3,A,I2)',"Ch4->",p0%octCh4%octN," iPos[",p0%octCh4%iPos,"] octP=",&
                p0%octCh4%octP,", Lv=",p0%octCh4%octLv 
           print '(A,I4,A,3I3,A,I3,A,I2)',"Ch5->",p0%octCh5%octN," iPos[",p0%octCh5%iPos,"] octP=",&
                p0%octCh5%octP,", Lv=",p0%octCh5%octLv 
           print '(A,I4,A,3I3,A,I3,A,I2)',"Ch6->",p0%octCh6%octN," iPos[",p0%octCh6%iPos,"] octP=",&
                p0%octCh6%octP,", Lv=",p0%octCh6%octLv 
           print '(A,I4,A,3I3,A,I3,A,I2)',"Ch7->",p0%octCh7%octN," iPos[",p0%octCh7%iPos,"] octP=",&
                p0%octCh7%octP,", Lv=",p0%octCh7%octLv 
           print '(A,I4,A,3I3,A,I3,A,I2)',"Ch8->",p0%octCh8%octN," iPos[",p0%octCh8%iPos,"] octP=",&
                p0%octCh8%octP,", Lv=",p0%octCh8%octLv 
        endif

        !check child iflg
        if(p0%iFLG(1)<4)then
           if(p0%octCh1%iFLG(1)>=0 .or. p0%octCh5%iFLG(1)>=0 .or. &
              p0%octCh2%iFLG(1)>=0 .or. p0%octCh6%iFLG(1)>=0 .or. &
              p0%octCh3%iFLG(1)>=0 .or. p0%octCh7%iFLG(1)>=0 .or. &
              p0%octCh4%iFLG(1)>=0 .or. p0%octCh8%iFLG(1)>=0 )then
              print *,"ID30=",minID(3,0),maxID(3,0),"ID31=",minID(3,1),maxID(3,1),rank
              print *,"ID40=",minID(4,0),maxID(4,0),"ID41=",minID(4,1),maxID(4,1),rank
              print *,"illegal flag connection with children rank=",rank
              print *,"index=",index,"iFLG=",p0%iFLG,"iLv=",p0%octLv,"type=",p0%octType,"prc_bndry=",p0%prc_bndry,rank
              print *,"Ch1 iFLG=",p0%octCh1%iFLG,"iLv=",p0%octCh1%octLv,"octN=",p0%octCh1%octN,&
                   "Ch5 iFLG=",p0%octCh5%iFLG,"iLv=",p0%octCh5%octLv,"octN=",p0%octCh5%octN,rank
              print *,"Ch2 iFLG=",p0%octCh2%iFLG,"iLv=",p0%octCh2%octLv,"octN=",p0%octCh2%octN,&
                   "Ch6 iFLG=",p0%octCh6%iFLG,"iLv=",p0%octCh6%octLv,"octN=",p0%octCh6%octN,rank
              print *,"Ch3 iFLG=",p0%octCh3%iFLG,"iLv=",p0%octCh3%octLv,"octN=",p0%octCh3%octN,&
                   "Ch7 iFLG=",p0%octCh7%iFLG,"iLv=",p0%octCh7%octLv,"octN=",p0%octCh7%octN,rank
              print *,"Ch4 iFLG=",p0%octCh4%iFLG,"iLv=",p0%octCh4%octLv,"octN=",p0%octCh4%octN,&
                   "Ch8 iFLG=",p0%octCh8%iFLG,"iLv=",p0%octCh8%octLv,"octN=",p0%octCh8%octN,rank
              
              errFlag=1
              errFlag_oct=1
           endif
        endif

        !parent->child->parent pointer check
        if(p0%octCh1%octType/=p0%octType .and. p0%octLv/=-1)then 
           print *,"OctCh1 is not same type",rank
           print *,"octCh1 octType=",p0%octCh1%octType,"Prt octType",p0%octType,rank
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh2%octType>p0%octType .and. p0%octLv/=-1)then 
           print *,"OctCh2 is not same type",rank
           print *,"octCh2 octType=",p0%octCh2%octType,"Prt octType",p0%octType,rank
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh3%octType>p0%octType .and. p0%octLv/=-1)then 
           print *,"OctCh3 is not same type",rank
           print *,"octCh3 octType=",p0%octCh3%octType,"Prt octType",p0%octType,rank
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh4%octType>p0%octType .and. p0%octLv/=-1)then 
           print *,"OctCh4 is not same type",rank
           print *,"octCh4 octType=",p0%octCh4%octType,"Prt octType",p0%octType,rank
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh5%octType>p0%octType .and. p0%octLv/=-1)then 
           print *,"OctCh5 is not same type",rank
           print *,"octCh5 octType=",p0%octCh5%octType,"Prt octType",p0%octType,rank
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh6%octType>p0%octType .and. p0%octLv/=-1)then 
           print *,"OctCh6 is not same type",rank
           print *,"octCh6 octType=",p0%octCh6%octType,"Prt octType",p0%octType,rank
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh7%octType>p0%octType .and. p0%octLv/=-1)then 
           print *,"OctCh7 is not same type",rank
           print *,"octCh7 octType=",p0%octCh7%octType,"Prt octType",p0%octType,rank
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh8%octType>p0%octType .and. p0%octLv/=-1)then 
           print *,"OctCh8 is not same type",rank
           print *,"octCh8 octType=",p0%octCh8%octType,"Prt octType",p0%octType,rank
           errFlag=1
           errFlag_oct=1
        endif

        if(p0%octCh1%octPrt%octN/=index)then 
           print *,"OctCh1 is not connect to this oct octN=",p0%octN,"index=",index,rank
           print *,"octCh1 octN=",p0%octCh1%octN,"Prt octN",p0%octCh1%octPrt%octN,&
                "myiFLG=",p0%iFLG,"octCh1 iFLG=",p0%iFLG,rank
           !stop
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh2%octPrt%octN/=index)then 
           print *,"OctCh2 is not connect to this oct octN=",p0%octN,"index=",index,rank
           print *,"octCh3 octN=",p0%octCh2%octN,"Prt octN",p0%octCh2%octPrt%octN,&
                "myiFLG=",p0%iFLG,"octCh2 iFLG=",p0%iFLG,rank
           !stop
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh3%octPrt%octN/=index)then
           print *,"OctCh3 is not connect to this oct octN=",p0%octN,"index=",index,rank
           print *,"octCh4 octN=",p0%octCh3%octN,"Prt octN",p0%octCh3%octPrt%octN,&
                "myiFLG=",p0%iFLG,"octCh3 iFLG=",p0%iFLG,rank
           !stop
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh4%octPrt%octN/=index)then
           print *,"OctCh4 is not connect to this oct octN=",p0%octN,"index=",index,rank
           print *,"octCh4 octN=",p0%octCh4%octN,"Prt octN",p0%octCh4%octPrt%octN,&
                "myiFLG=",p0%iFLG,"octCh4 iFLG=",p0%iFLG,rank
           !stop
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh5%octPrt%octN/=index)then
           print *,"OctCh5 is not connect to this oct octN=",p0%octN,"index=",index,rank
           print *,"octCh5 octN=",p0%octCh5%octN,"Prt octN",p0%octCh5%octPrt%octN,&
                "myiFLG=",p0%iFLG,"octCh5 iFLG=",p0%iFLG,rank
           !stop
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh6%octPrt%octN/=index)then
           print *,"OctCh6 is not connect to this oct octN=",p0%octN,"index=",index,rank
           print *,"octCh6 octN=",p0%octCh6%octN,"Prt octN",p0%octCh6%octPrt%octN,&
                "myiFLG=",p0%iFLG,"octCh6 iFLG=",p0%iFLG,rank
           !stop
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh7%octPrt%octN/=index)then
           print *,"OctCh7 is not connect to this oct octN=",p0%octN,"index=",index,rank
           print *,"octCh7 octN=",p0%octCh7%octN,"Prt octN",p0%octCh7%octPrt%octN,&
                "myiFLG=",p0%iFLG,"octCh7 iFLG=",p0%iFLG,rank
           !stop
           errFlag=1
           errFlag_oct=1
        endif
        if(p0%octCh8%octPrt%octN/=index)then
           print *,"OctCh8 is not connect to this oct octN=",p0%octN,"index=",index,rank
           print *,"octCh8 octN=",p0%octCh8%octN,"Prt octN",p0%octCh8%octPrt%octN,&
                "myiFLG=",p0%iFLG,"octCh8 iFLG=",p0%iFLG,rank
           !stop
           errFlag=1
           errFlag_oct=1
        endif
     else
        print *,"This oct doesnt have 8 children!",rank
        print *,"index=",index,"octType=",p0%octType,"octLv=",p0%octLv,rank
        print *,"MrtN=",p0%MrtN,"iFLG=",p0%iFLG,rank

        if(.not. associated(p0%octCh1))print *,"no Ch1",rank
        if(.not. associated(p0%octCh2))print *,"no Ch2",rank
        if(.not. associated(p0%octCh3))print *,"no Ch3",rank
        if(.not. associated(p0%octCh4))print *,"no Ch4",rank
        if(.not. associated(p0%octCh5))print *,"no Ch5",rank
        if(.not. associated(p0%octCh6))print *,"no Ch6",rank
        if(.not. associated(p0%octCh7))print *,"no Ch7",rank
        if(.not. associated(p0%octCh8))print *,"no Ch8",rank

     endif
  endif

  if(outType==1)print *,""

  if(errFlag==1 .or. errFlag_cnct==1)then
     errFlag=1
  endif

  if(errFlag_oct/=0 .and. outType/=1)then
     !call output_oct_and_family(p0%octN,LvFlag,nowLv,MeshState,1,errFlag)
     write(filename,2311)version,rank
!     open(98,file=filename,status="old",form="formatted",position="append")
!     write(98,'(3i7.1)')p0%iPOS
!     close(98)
  endif

  if(p0%iFLG(1)>0 .and. LvFlag>nowLv)then
     call output_oct_and_family(p0%octCh1%octN,LvFlag,nowLv+1,MeshState,outType,errFlag)
     call output_oct_and_family(p0%octCh2%octN,LvFlag,nowLv+1,MeshState,outType,errFlag)
     call output_oct_and_family(p0%octCh3%octN,LvFlag,nowLv+1,MeshState,outType,errFlag)
     call output_oct_and_family(p0%octCh4%octN,LvFlag,nowLv+1,MeshState,outType,errFlag)
     call output_oct_and_family(p0%octCh5%octN,LvFlag,nowLv+1,MeshState,outType,errFlag)
     call output_oct_and_family(p0%octCh6%octN,LvFlag,nowLv+1,MeshState,outType,errFlag)
     call output_oct_and_family(p0%octCh7%octN,LvFlag,nowLv+1,MeshState,outType,errFlag)
     call output_oct_and_family(p0%octCh8%octN,LvFlag,nowLv+1,MeshState,outType,errFlag)
  endif

2311 format("./output/",A,"_DOP",i3.3,".dat")

endsubroutine output_oct_and_family

subroutine output_oct_map(index,SX,SY,SZ,outtype)
  use init_mesh_size
  use oct_set
  use message_passing_interface
  implicit none
  integer(kind=4),intent(in)::SX,SY,SZ,outtype,index
  integer(kind=4)::i,j,k,x,y,z
  integer(kind=4)::hSX,hSY,hSZ
  integer(kind=4),dimension(SX,SY,SZ)::outData
  type(oct),pointer::p0,p1,p2,p3

  hSX=SX/2+1
  hSY=SY/2+1
  hSZ=SZ/2+1

  print *,"SX =", SX,"SY =", SY,"SZ =",SZ,rank
  print *,"hSX=",hSX,"hSY=",hSY,"hSZ=",hSZ,rank

  !SX,SY,SZ must be odd number
  if(mod(SX,2)==0 .or. mod(SY,2)==0 .or. mod(SZ,2)==0)then
     print *,"SX,SY,SZ must be odd number SX=",SX,"SY=",SY,"SZ=",SZ,rank
     return
  endif

  outData=-1000
  p0=>Mesh(index)
  
  !1: -x -y -z
  p1=>p0
  do z=1,hSZ
     p2=>p1
     do y=1,hSY
        p3=>p2
        do x=1,hSX
           if(associated(p3%octNb1))then
              if(outType==0)outData(hSX-x+1,hSY-y+1,hSZ-z+1)=p3%iFLG(1)
              if(outType==1)outData(hSX-x+1,hSY-y+1,hSZ-z+1)=p3%MrtN
              if(outType==2)outData(hSX-x+1,hSY-y+1,hSZ-z+1)=p3%octType
              if(p3%octN/=p3%octNb1%octN)then
                 p3=>p3%octNb1
              else
                 exit
              endif
           endif
        enddo
        if(associated(p2%octNb3))then
           if(p2%octN/=p2%octNb3%octN)p2=>p2%octNb3
        endif
     enddo
     if(associated(p1%octNb5))then
        if(p1%octN/=p1%octNb5%octN)p1=>p1%octNb5
     endif
  enddo
  !2: +x -y -z
  p1=>p0
  do z=1,hSZ
     p2=>p1
     do y=1,hSY
        p3=>p2
        do x=1,hSX
           if(associated(p3%octNb2))then
              if(outType==0)outData(hSX+x-1,hSY-y+1,hSZ-z+1)=p3%iFLG(1)
              if(outType==1)outData(hSX+x-1,hSY-y+1,hSZ-z+1)=p3%MrtN
              if(outType==2)outData(hSX+x-1,hSY-y+1,hSZ-z+1)=p3%octType
              if(p3%octN/=p3%octNb2%octN)then
                 p3=>p3%octNb2
              else
                 exit
              endif
           endif
        enddo
        if(associated(p2%octNb3))then
           if(p2%octN/=p2%octNb3%octN)p2=>p2%octNb3
        endif
     enddo
     if(associated(p1%octNb5))then
        if(p1%octN/=p1%octNb5%octN)p1=>p1%octNb5
     endif
  enddo
  !3: -x +y -z
  p1=>p0
  do z=1,hSZ
     p2=>p1
     do y=1,hSY
        p3=>p2
        do x=1,hSX
           if(associated(p3%octNb1))then
              if(outType==0)outData(hSX-x+1,hSY+y-1,hSZ-z+1)=p3%iFLG(1)
              if(outType==1)outData(hSX-x+1,hSY+y-1,hSZ-z+1)=p3%MrtN
              if(outType==2)outData(hSX-x+1,hSY+y-1,hSZ-z+1)=p3%octType
              if(p3%octN/=p3%octNb1%octN)then
                 p3=>p3%octNb1
              else
                 exit
              endif
           endif
        enddo
        if(associated(p2%octNb4))then
           if(p2%octN/=p2%octNb4%octN)p2=>p2%octNb4
        endif
     enddo
     if(associated(p1%octNb5))then
        if(p1%octN/=p1%octNb5%octN)p1=>p1%octNb5
     endif
  enddo
  !4: +x +y -z
  p1=>p0
  do z=1,hSZ
     p2=>p1
     do y=1,hSY
        p3=>p2
        do x=1,hSX
           if(associated(p3%octNb2))then
              if(outType==0)outData(hSX+x-1,hSY+y-1,hSZ-z+1)=p3%iFLG(1)
              if(outType==1)outData(hSX+x-1,hSY+y-1,hSZ-z+1)=p3%MrtN
              if(outType==2)outData(hSX+x-1,hSY+y-1,hSZ-z+1)=p3%octType
              if(p3%octN/=p3%octNb2%octN)then
                 p3=>p3%octNb2
              else
                 exit
              endif
           endif
        enddo
        if(associated(p2%octNb4))then
           if(p2%octN/=p2%octNb4%octN)p2=>p2%octNb4
        endif
     enddo
     if(associated(p1%octNb5))then
        if(p1%octN/=p1%octNb5%octN)p1=>p1%octNb5
     endif
  enddo  
  !5: -x -y +z
  p1=>p0
  do z=1,hSZ
     p2=>p1
     do y=1,hSY
        p3=>p2
        do x=1,hSX
           if(associated(p3%octNb1))then
              if(outType==0)outData(hSX-x+1,hSY-y+1,hSZ+z-1)=p3%iFLG(1)
              if(outType==1)outData(hSX-x+1,hSY-y+1,hSZ+z-1)=p3%MrtN
              if(outType==2)outData(hSX-x+1,hSY-y+1,hSZ+z-1)=p3%octType
              if(p3%octN/=p3%octNb1%octN)then
                 p3=>p3%octNb1
              else
                 exit
              endif
           endif
        enddo
        if(associated(p2%octNb3))then
           if(p2%octN/=p2%octNb3%octN)p2=>p2%octNb3
        endif
     enddo
     if(associated(p1%octNb6))then
        if(p1%octN/=p1%octNb6%octN)p1=>p1%octNb6
     endif
  enddo
  !6: +x -y +z
  p1=>p0
  do z=1,hSZ
     p2=>p1
     do y=1,hSY
        p3=>p2
        do x=1,hSX
           if(associated(p3%octNb2))then
              if(outType==0)outData(hSX+x-1,hSY-y+1,hSZ+z-1)=p3%iFLG(1)
              if(outType==1)outData(hSX+x-1,hSY-y+1,hSZ+z-1)=p3%MrtN
              if(outType==2)outData(hSX+x-1,hSY-y+1,hSZ+z-1)=p3%octType
              if(p3%octN/=p3%octNb2%octN)then
                 p3=>p3%octNb2
              else
                 exit
              endif
           endif
        enddo
        if(associated(p2%octNb3))then
           if(p2%octN/=p2%octNb3%octN)p2=>p2%octNb3
        endif
     enddo
     if(associated(p1%octNb6))then
        if(p1%octN/=p1%octNb6%octN)p1=>p1%octNb6
     endif
  enddo
  !7: -x +y +z
  p1=>p0
  do z=1,hSZ
     p2=>p1
     do y=1,hSY
        p3=>p2
        do x=1,hSX
           if(associated(p3%octNb1))then
              if(outType==0)outData(hSX-x+1,hSY+y-1,hSZ+z-1)=p3%iFLG(1)
              if(outType==1)outData(hSX-x+1,hSY+y-1,hSZ+z-1)=p3%MrtN
              if(outType==2)outData(hSX-x+1,hSY+y-1,hSZ+z-1)=p3%octType
              if(p3%octN/=p3%octNb1%octN)then
                 p3=>p3%octNb1
              else
                 exit
              endif
           endif
        enddo
        if(associated(p2%octNb4))then
           if(p2%octN/=p2%octNb4%octN)p2=>p2%octNb4
        endif
     enddo
     if(associated(p1%octNb6))then
        if(p1%octN/=p1%octNb6%octN)p1=>p1%octNb6
     endif
  enddo
  !8: +x +y +z
  p1=>p0
  do z=1,hSZ
     p2=>p1
     do y=1,hSY
        p3=>p2
        do x=1,hSX
           if(associated(p3%octNb2))then
              if(outType==0)outData(hSX+x-1,hSY+y-1,hSZ+z-1)=p3%iFLG(1)
              if(outType==1)outData(hSX+x-1,hSY+y-1,hSZ+z-1)=p3%MrtN
              if(outType==2)outData(hSX+x-1,hSY+y-1,hSZ+z-1)=p3%octType
              if(p3%octN/=p3%octNb2%octN)then
                 p3=>p3%octNb2
              else
                 exit
              endif
           endif
        enddo
        if(associated(p2%octNb4))then
           if(p2%octN/=p2%octNb4%octN)p2=>p2%octNb4
        endif
     enddo
     if(associated(p1%octNb6))then
        if(p1%octN/=p1%octNb6%octN)p1=>p1%octNb6
     endif
  enddo
  
  !output
  print *,"output oct map : centerPOS=",p0%iPOS,rank
  print *,"============================================",rank
  do j=1,SY
     print *,"centerPOS=",p0%iPOS(1:2),p0%iPOS(3)-hSY+j,rank
     print *,"---------------------------------------------",rank
     do k=1,SZ
        do i=1,SX
           if(outData(i,j,k)==-1000)then
              print '(A,$)'," nu"
           else
              if(outType==0)print '(I3,$)',outData(i,j,k)
              if(outType==1)print '(I5,$)',outData(i,j,k)
              if(outType==2)print '(I2,$)',outData(i,j,k)
           endif
        enddo
        print *,rank
     enddo
     print *,"---------------------------------------------",rank
  enddo
  print *,"============================================",rank

end subroutine output_oct_map

subroutine debug_oct_pack
  use param
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::i,j

  do j=0,nprocs-1
     if(j==rank)then
        if(sendSize(1)>0)then
           do i=1,sendSize(1)
              print *,"Left Send i,octP,iFLG,FR",i,sbufL_FI(i*2),sbufL_FI(i*2-1),sbufL_FR(1+(i-1)*12),rank
           enddo
        endif
        if(sendSize(3)>0)then
           do i=1,sendSize(3)
              print *,"Right Send i,octP,iFLG,FR",i,sbufR_FI(i*2),sbufR_FI(i*2-1),sbufR_FR(1+(i-1)*12),rank
           enddo
        endif

        if(recvSize(1)>0)then
           do i=1,recvSize(1)
              print *,"Left Recv i,octP,iFLG,FR",i,rbufL_FI(i*2),rbufL_FI(i*2-1),rbufL_FR(1+(i-1)*12),rank
           enddo
        endif
        if(recvSize(3)>0)then
           do i=1,recvSize(3)
              print *,"RightRecv i,octP,iFLG,FR",i,rbufR_FI(i*2),rbufR_FI(i*2-1),rbufR_FR(1+(i-1)*12),rank
           enddo
        endif
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo

end subroutine debug_oct_pack

subroutine debug_set_n_gcells
  use init_mesh_size
  use message_passing_interface
  use param
  use oct_set
  implicit none
  integer(kind=4)::index,i,indexP,iLv
  real(kind=8) :: R(9)
  type(oct),pointer::p0
  type(prtcl),pointer::PrtList

  !this subroutine is used for debugging [refresh_Particle] on 2Procs
  !and override n_gcell_proc and iGMesh_arr
  
  R(:)=0d0
  i=0
  do index=MinID(3,0),MaxID(3,0)
     p0=>Mesh(index)

     iLv=p0%octLv
    
     if(iLv/=-1)then
         !print *,"index=",index,"iLv=",iLv,rank
        indexP=MaxIP(iLv)
        indexP=indexP+1
        if(.not. associated(p0%ptcl))then
           nullify(PrtList)
           R(1:3) = p0%rPos(:)
           Pesh(indexP,iLv)=prtcl(R,-1,p0%octN,indexP,PrtList,0)
           p0%ptcl=>Pesh(indexP,iLv)
           indexP=indexP+1
           !print *,"add representative P into ",index,rank
        endif
        !pp=>p0%ptcl%prtnxt
        PrtList=>p0%ptcl%prtnxt
        R(1:3) = p0%rPos(:)
        Pesh(indexP,iLv)=prtcl(R,p0%MrtN+1,p0%octN,indexP,Prtlist,0)
        p0%ptcl%prtnxt=>Pesh(indexP,iLv)
        MaxIP(iLv)=indexP
        p0%octP=1
        cycle
     endif
     i=i+1
     iGMesh_arr(i)=index
  enddo

  do index=MinID(3,-1),MaxID(3,-1)
     p0=>Mesh(index)
     
     p0%octP=8
  enddo

!!$!  if(rank==0)then
!!$!     n_gcells_proc(1)=0
!!$!     n_gcells_proc(2)=64
!!$!
!!$!  else
!!$!     n_gcells_proc(2)=0
!!$!     n_gcells_proc(1)=64
!!$!  endif
  
!!$!  print *,"Parameter Changed : n_gcells_proc(1)=",n_gcells_proc(1)
!!$!  print *,"Parameter Changed : n_gcells_proc(2)=",n_gcells_proc(2)

end subroutine debug_set_n_gcells

subroutine debug_output_received_gcell
  use init_mesh_size
  use oct_set
  use message_passing_interface
  use param
  implicit none
  integer(kind=4)::index
  type(oct),pointer::p0
  type(prtcl),pointer::pp

  do index=MinID(1,0),MaxID(1,0)
     p0=>Mesh(index)
     !print *,"index=",index,"octP=",p0%octP
     if(p0%octP>0)then
        pp=>p0%ptcl%prtnxt
        print *,"p0%octN=",p0%octN,"p0%MrtN=",p0%MrtN,"pp%isort=",pp%isort
        if(p0%MrtN+1/=pp%isort)then
           print *,"error?"
           stop
        endif
     else
        
     endif
  enddo

end subroutine debug_output_received_gcell

!***********************************************************************************************
!----------------------------- Output For Display Routines ------------------------------------
!***********************************************************************************************
!2011/09/12 yagi added
subroutine Display_GMesh_All
  use param
  implicit none
  integer(kind=4):: iLv
  
  do iLv=-1,LvMax
     call Display_GMesh(iLv)
  end do
end subroutine Display_GMesh_All

subroutine Display_Grid_All
  use param
  implicit none
  integer(kind=4):: iLv

  if(LvMax>=0)then
     do iLv=0,LvMax
        call Display_Grid(iLv,-1,-1,NZ*NZR/2)
     enddo
  endif
end subroutine Display_Grid_All
!additional part end

subroutine Display_Grid(iLv,cutx,cuty,cutz)
  use oct_set
  use param
  use const
  use init_mesh_size
  use message_passing_interface
!-----------------------------------------
  implicit none
  integer(kind=4),intent(in)::iLv,cutX,cutY,cutZ
  integer(kind=4):: index, iPOSN(3), Nint,i
  integer(kind=4) :: cutoff,octCut
  integer(kind=4) :: iFLG
  integer(kind=4) :: filenum
  character(LEN=128) :: fo_name  
  type(oct), pointer :: p0
!-----------------------------------------
  if(MaxID(1,iLv)==MinID(1,iLv) .and. MaxID(3,iLv)==MinID(3,iLv))return

  filenum = 130+rank

  octCut=OoctIntvl
  cutoff=8**(iLv+1)

  write(fo_name,556) version,istep,rank,iLv
  open (filenum, file=fo_name, form=  'formatted')

  do i=1,4
     if(MinID(i,iLv)<MaxID(i,iLv))then
        do index=MinID(i,iLv), MaxID(i,iLv)
           !
           if(octCut==1 .and. mod(index,cutoff)/=0)cycle
           !

           p0 => Mesh(index)
           if(p0%iFLG(1)>0 .and. p0%iFLG(1)<=6)then

              Nint = 2**(LvMax-iLv)
              iPOSN = (p0%iPOS+Nint)/(2*Nint)
              iFLG=p0%iFLG(1)
              if(p0%iFLG(1)>0)then
                 if(p0%iFLG(1)>=4)then
                    if(p0%octCh1%iFLG(1)<0 .or. p0%octCh2%iFLG(1)<0 .or.&
                         p0%octCh3%iFLG(1)<0 .or. p0%octCh4%iFLG(1)<0 .or.&
                         p0%octCh5%iFLG(1)<0 .or. p0%octCh6%iFLG(1)<0 .or.&
                         p0%octCh7%iFLG(1)<0 .or. p0%octCh8%iFLG(1)<0 )then
                       iFLG=10 !error marker
                    endif
                 endif
                 if(p0%iFLG(1)>0 .and. p0%iFLG(1)<4)then
                    if(p0%octCh1%iFLG(1)/=p0%iFLG(1)-4 .or. p0%octCh2%iFLG(1)/=p0%iFLG(1)-4 .or.&
                         p0%octCh3%iFLG(1)/=p0%iFLG(1)-4 .or. p0%octCh4%iFLG(1)/=p0%iFLG(1)-4 .or.&
                         p0%octCh5%iFLG(1)/=p0%iFLG(1)-4 .or. p0%octCh6%iFLG(1)/=p0%iFLG(1)-4 .or.&
                         p0%octCh7%iFLG(1)/=p0%iFLG(1)-4 .or. p0%octCh8%iFLG(1)/=p0%iFLG(1)-4 &
                         )then
                       iFLG=10 !error marker
                    endif
                 endif
                 if(p0%octPrt%iFLG(1)<=0)then
                    iFLG=-10
                 endif
              endif

              if(octCut==2 .and. cutX>0 .and. iPOSN(1)/=cutX)cycle
              if(octCut==2 .and. cutY>0 .and. iPOSN(2)/=cutY)cycle
              if(octCut==2 .and. cutZ>0 .and. iPOSN(3)/=cutZ)cycle

              write(filenum,'(3F15.10,3I3)')p0%rPos,iFLG
           endif
        enddo
     endif
  enddo


  close(filenum)

556    format('./output/Data/',A,'_Grid',I5.5,'_P',I3.3,'_Lv',I2.2,'.dat') 
  
end subroutine Display_Grid

subroutine Display_GMesh(iLv)
  use oct_set
  use param
  use const
  use init_mesh_size
  use message_passing_interface
!-----------------------------------------
  implicit none
  integer(kind=4):: index, iPOSN(3), iPOSH(3), iLv,Nint,i
  integer(kind=4) :: ptcl_loops
  integer(kind=4) :: Nb(6)
  integer(kind=4) :: filenum
  real(kind=4):: rPOS(3)

  character(LEN=32) :: fo_name  
  type(oct), pointer :: p0
!-----------------------------------------
  
  filenum = 130+rank
  write(fo_name,556) './output/Data/GMesh_',rank,istep,(iLv+1)
  open (filenum, file=fo_name   , form=  'formatted')
  
  do i=3,4
     do index=MinID(i,iLv), MaxID(i,iLv)
        if(mod(index,8**(iLv+1))/=0)cycle
        Nb(:)=0
        p0 => Mesh(index)
        iPOSH = p0 % iPOS
        Nint = 2**(LvMax-iLv)
        iPOSN = (iPOSH+Nint)/(2*Nint)
        rPOS = p0 % rPOS
        ptcl_loops = p0% octP

        if( associated(p0%octNb1)) Nb(1) =1
        if( associated(p0%octNb2)) Nb(2) =1
        if( associated(p0%octNb3)) Nb(3) =1
        if( associated(p0%octNb4)) Nb(4) =1
        if( associated(p0%octNb5)) Nb(5) =1
        if( associated(p0%octNb6)) Nb(6) =1

        write(filenum,*)rPOS,p0%iFLG,p0%Csort
     enddo
  enddo
     !8 for octNb1, 14 for iFLG(1), 18 for prc_bndry
  close(filenum)

556    format(A,I2.2,'_',I5.5,'_',I2.2,'.dat') 
  
end subroutine Display_GMesh

!***********************************************************************
      subroutine out_particle(iLv)
! +-------------------------------------------------------------------+
! |                                                                   |
! |  output field data                                                |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
!-----------------------------------------
        implicit none
!        integer(kind=4)   :: istep
        integer(kind=4)   :: index,indexP,iLv
        character(LEN=32) :: fo_name
!------------------------------------------
        type(prtcl), pointer :: pp
!------------------------------------------
!
        write(fo_name,555) 'prIONG',iLv,istep
        open (93, file=fo_name   , form=  'formatted')
        write(fo_name,555) 'prELEG',iLv,istep
        open (94, file=fo_name   , form=  'formatted')
!
        indexP = maxIP(iLv)
        do index=1,indexP
          pp => Pesh(index,iLv)
          if(pp%Isort.eq.2) write(93,*) index,pp%R(1),pp%R(2) !,pp%Ioct,pp%Inum
          if(pp%Isort.eq.1) write(94,*) index,pp%R(1),pp%R(2) !,pp%Ioct,pp%Inum
        enddo
!
        close(93)
        close(94)
!---------------------------------------------------------
!
555    format(A,I2.2,'_',I6.6,'.datPL2') 
!
        return
      end subroutine out_particle
!
!***********************************************************************
      subroutine out_particle2(iLv)
! +-------------------------------------------------------------------+
! |                                                                   |
! |  output field data                                                |
! |                                                                   |
! +-------------------------------------------------------------------+
!***********************************************************************
        use oct_set
        use param
        use const
        use init_mesh_size
!-----------------------------------------
        implicit none
        integer(kind=4)   :: i,minindex,maxindex
        integer(kind=4)   :: index,iLv
        character(LEN=32) :: fo_name
!------------------------------------------
        type(  oct), pointer :: p0
        type(prtcl), pointer :: pp
!------------------------------------------
!
        write(fo_name,555) 'prIONx',iLv,istep
        open (93, file=fo_name   , form=  'formatted')
        write(fo_name,555) 'prELEx',iLv,istep
        open (94, file=fo_name   , form=  'formatted')
!
        maxindex=maxID(1,iLv)
        minindex=minID(1,iLv)

        do index=minindex,maxindex
           p0 => Mesh(index)
           pp => p0%ptcl%prtnxt
           do i=1,p0%octP
              if(pp%isort.eq.1) write(93,*) index,pp%R(1),pp%R(2)
              if(pp%isort.eq.2) write(94,*) index,pp%R(1),pp%R(2)
              pp => pp%prtnxt
           end do
        end do
!
        close(93)
        close(94)
!---------------------------------------------------------
!
555    format(A,I2.2,'_',I6.6,'.datPL2') 
!
        return
      end subroutine out_particle2
!
!-----------------------------------------------
      subroutine output_field(iLv)
!-----------------------------------------------
        use param
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        implicit none
        integer(kind=4),parameter        :: NXB2=NXB*4,NYB2=NYB*4
        integer(kind=4)                  :: iLv
        integer(kind=4)                  :: i,j,index,ix,iy
        integer(kind=4)                  :: minindex,maxindex
        integer(kind=4)                  :: iFLG1,iFLG0,NYBh
        real(kind=8)                     :: di
!-- output data
        real(kind=8)                     :: x,y,r,xhalf,yhalf
        integer(kind=4),dimension(NXB2,NYB2):: octdata1
        real(kind=8)   ,dimension(NXB2,NYB2)::  Exdata1, Eydata1, Ezdata1 &
                                              , Bxdata1, Bydata1, Bzdata1 &
                                              , Redata1, Ridata1          &
                                              ,Jexdata1,Jeydata1,Jezdata1 !&
!                                              ,Jixdata1,Jiydata1,Jizdata1
        real(kind=8)   ,dimension(NXB2,NYB2):: DPY,DPX,SCz,SDz,SEz
        real(kind=8)   ,dimension(NXB2) :: PBY0,DBY0,BBY0,CRE0,CRI0,JEZ0
        character(LEN=32) :: fo_name
! -- pointers --
        type(oct), pointer   :: p0
      
        if(iLv.gt.2) return
        maxindex=maxID(1,iLv)
        minindex=minID(1,iLv)
!
        octdata1=0
         Exdata1=0.d0 ;  Eydata1=0.d0 ;  Ezdata1=0.d0
         Bxdata1=0.d0 ;  Bydata1=0.d0 ;  Bzdata1=0.d0
         Redata1=0.d0 ;  Ridata1=0.d0 ;     
        Jexdata1=0.d0 ; Jeydata1=0.d0 ; Jezdata1=0.d0        

        do index=minindex,maxindex
           p0 => Mesh(index)
           di = dble(2**iLv)
           ix = int(p0%rPOS(1)*di/(dx(1)))+1
           iy = int(p0%rPOS(2)*di/(dx(2)))+1
           iFLG0=0
           iFLG1=p0%iFLG(1)
           if((iFLG1.ge.0).and.(iFLG1.le.3)) iFLG0=2

           octdata1(ix,iy)=octdata1(ix,iy)+iFLG0
           Exdata1( ix,iy)=Exdata1( ix,iy)+p0%F( 1)
           Eydata1( ix,iy)=Eydata1( ix,iy)+p0%F( 2)
           Ezdata1( ix,iy)=Ezdata1( ix,iy)+p0%F( 3)
           Bxdata1( ix,iy)=Bxdata1( ix,iy)+p0%F( 4)
           Bydata1( ix,iy)=Bydata1( ix,iy)+p0%F( 5)
           Bzdata1( ix,iy)=Bzdata1( ix,iy)+p0%F( 6)
           Redata1( ix,iy)=Redata1( ix,iy)+p0%Z( 1)
           Ridata1( ix,iy)=Ridata1( ix,iy)+p0%Z( 2)
           Jexdata1(ix,iy)=Jexdata1(ix,iy)+p0%F(10)
           Jeydata1(ix,iy)=Jeydata1(ix,iy)+p0%F(11)
           Jezdata1(ix,iy)=Jezdata1(ix,iy)+p0%F(12)
        end do

        xhalf=R_lim(1,1)*0.5d0
        yhalf=R_lim(1,2)*0.5d0
        do j=1,NYB2
           do i=1,NXB2
              x=dble(i-1)*dx(1)*0.25d0-xhalf
              y=dble(j-1)*dx(2)*0.25d0-yhalf
              r=x**2+y**2
              if(r.gt.1.d-8) then 
                 DPX(i,j)=-(x*y*dx(1)*dx(2)*BD)/(r**2)
                 DPY(i,j)=((x**2-y**2)*dx(1)*dx(2)*BD)/(r**2)/2.d0
              else
                 DPX(i,j)=0.d0
                 DPY(i,j)=0.d0
              end if
           end do
        end do
!


        call copy_func0(octdata1,NXB2,NYB2,iLv)
!
        call copy_func(  Exdata1,NXB2,NYB2,iLv)
        call copy_func(  Eydata1,NXB2,NYB2,iLv)
        call copy_func(  Ezdata1,NXB2,NYB2,iLv)
        call copy_func(  Bxdata1,NXB2,NYB2,iLv)
        call copy_func(  Bydata1,NXB2,NYB2,iLv)
        call copy_func(  Bzdata1,NXB2,NYB2,iLv)
!
        call copy_func( Jexdata1,NXB2,NYB2,iLv)
        call copy_func( Jeydata1,NXB2,NYB2,iLv)
        call copy_func( Jezdata1,NXB2,NYB2,iLv)
!
        call copy_func(  Redata1,NXB2,NYB2,iLv)
        call copy_func(  Ridata1,NXB2,NYB2,iLv)
!

        di = dble(NZB)*dble(2**(iLv))
        do j=1,NYB2
           do i=1,NXB2
              octdata1(i,j)=int(dble(octdata1(i,j))/di)
              Exdata1( i,j)=Exdata1( i,j)/di
              Eydata1( i,j)=Eydata1( i,j)/di
              Ezdata1( i,j)=Ezdata1( i,j)/di
              Bxdata1( i,j)=Bxdata1( i,j)/di
              Bydata1( i,j)=Bydata1( i,j)/di
              Bzdata1( i,j)=Bzdata1( i,j)/di
              Redata1( i,j)=Redata1( i,j)/di
              Ridata1( i,j)=Ridata1( i,j)/di
              Jexdata1(i,j)=Jexdata1(i,j)/di
              Jeydata1(i,j)=Jeydata1(i,j)/di
              Jezdata1(i,j)=Jezdata1(i,j)/di
           end do
        end do

!
        SCz=0.d0
        SDz=0.d0
        SEz=0.d0
        do i=2,NXB2
           SCz(i,1)=SCz(i-1,1)- Bydata1(i-1,1)
           SDz(i,1)=SDz(i-1,1)- DPY(i-1,1)
           SEz(i,1)=SEz(i-1,1)-(Bydata1(i-1,1)+DPY(i-1,1))
        end do
        do j=2,NYB2
           do i=1,NXB2
              SCz(i,j)=SCz(i,j-1)+ Bxdata1(i,j-1)
              SDz(i,j)=SDz(i,j-1)+ DPX(i,j-1)
              SEz(i,j)=SEz(i,j-1)+(Bxdata1(i,j-1)+DPX(i,j-1))
           end do
        end do
!

        NYBh=NYB2/2
        do i=1,NXB2
           PBY0(i)=dble(i-1)*dx(1)*0.5d0
           DBY0(i)=(     DPY(i,NYBh-1)+     DPY(i,NYBh)+     DPY(i,NYBh+1)+     DPY(i,NYBh+2))/4.d0
           BBY0(i)=( Bydata1(i,NYBh-1)+ Bydata1(i,NYBh)+ Bydata1(i,NYBh+1)+ Bydata1(i,NYBh+2))/4.d0
           CRE0(i)=( Redata1(i,NYBh-1)+ Redata1(i,NYBh)+ Redata1(i,NYBh+1)+ Redata1(i,NYBh+2))/4.d0
           CRI0(i)=( Ridata1(i,NYBh-1)+ Ridata1(i,NYBh)+ Ridata1(i,NYBh+1)+ Ridata1(i,NYBh+2))/4.d0
           JEZ0(i)=(Jezdata1(i,NYBh-1)+Jezdata1(i,NYBh)+Jezdata1(i,NYBh+1)+Jezdata1(i,NYBh+2))/4.d0
        end do

        write(fo_name,555) 'center',iLv,istep
        open (92, file=fo_name   , form=  'formatted')
        do i=1,NXB2
           write(92,*) PBY0(i),DBY0(i),BBY0(i),CRE0(i),CRI0(i),JEZ0(i)
        end do
        close(92)
        write(fo_name,555) 'plane_',iLv,istep
        open (93, file=fo_name   , form=  'unformatted')
        write(93)  octdata1                  &
                 , Exdata1, Eydata1, Ezdata1 &
                 , Bxdata1, Bydata1, Bzdata1 &
                 , Redata1, Ridata1          &
                 ,Jexdata1,Jeydata1,Jezdata1 &
                 , SCz,SDz,SEz
        close(93)

        write(fo_name,555) 'AVS_P',iLv,istep
        open (94, file=fo_name   , form=  'formatted')
        write(fo_name,555) 'AVS_c',iLv,istep
        open (95, file=fo_name   , form=  'formatted')
        write(fo_name,555) 'AVS_R',iLv,istep
        open (98, file=fo_name   , form=  'formatted')
        do j=1,NYB2/2
           do i=1,NXB2/2
              if(octdata1(i,j).gt.0) then 
                 write(94,*)  sngl(SDz(2*i-1,2*j-1)+SDz(2*i,2*j-1)+SDz(2*i-1,2*j)+SDz(2*i,2*j))/4.0  & 
                             ,sngl(SCz(2*i-1,2*j-1)+SCz(2*i,2*j-1)+SCz(2*i-1,2*j)+SCz(2*i,2*j))/4.0  &
                             ,sngl(SEz(2*i-1,2*j-1)+SEz(2*i,2*j-1)+SEz(2*i-1,2*j)+SEz(2*i,2*j))/4.0
                 write(95,*) sngl(dble(i)*dx(1)*0.5d0),sngl(dble(j)*dx(2)*0.5d0)
!                write(96,*)  sngl(Exdata1(i,j)), sngl(Eydata1(i,j)), sngl(Ezdata1(i,j))
!                write(97,*)  sngl(Bxdata1(i,j)), sngl(Bydata1(i,j)), sngl(Bzdata1(i,j))
                 write(98,*)  sngl(Redata1(2*i-1,2*j-1)+Redata1(2*i-1,2*j  )  &
                                  +Redata1(2*i  ,2*j-1)+Redata1(2*i  ,2*j  ))/4.0, &
                              sngl(Ridata1(2*i-1,2*j-1)+Ridata1(2*i-1,2*j  ) &
                                  +Ridata1(2*i  ,2*j-1)+Ridata1(2*i  ,2*j  ))/4.0 
!                write(99,*) sngl(Jexdata1(i,j)),sngl(Jeydata1(i,j)),sngl(Jezdata1(i,j))
              end if
           end do
        end do
        close(94)
        close(95)
        close(98)

555    format(A,I2.2,'_',I6.6,'.datPL2') 

        return 
        end subroutine output_field
!
!-----------------------------------------------
      subroutine output_wave(iLv)
!-----------------------------------------------
        use param
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        use message_passing_interface
        implicit none
        integer(kind=4)               :: iLv
        integer(kind=4)               :: i,index,ix,minindex,maxindex
        real(kind=8)                  :: di
        real(kind=8),dimension(NXB*LvMax2/2+1) :: debugE0,debugP0,debugT0,debugJ0
        integer(kind=4),dimension(NXB*LvMax2/2+1) :: debugI0
        character(LEN=32) :: fo_name
! -- pointers --
        type(oct), pointer   :: p0
       
        maxindex=maxID(1,iLv)
        minindex=minID(1,iLv)


        debugE0 = 0.d0
        debugP0 = 0.d0
        debugT0 = 0.d0
        debugJ0 = 0.d0
        debugI0 = 0
     
        di = dx(1)/dble((2**LvMax))

        do index=minindex,maxindex
           p0 => Mesh(index)
           ix = int((p0%rPOS(1)-R_lim(0,1))/dble(di))+1
           if((p0%iFLG(1).ge.-1).and.(p0%iFLG(1).le.4)) then
!           if( p0 % prc_bndry==0) then
              debugI0( ix)=debugI0( ix)+1
              debugE0( ix)=debugE0( ix)+dble(p0%iFLG(1))       
              debugT0( ix)=debugT0( ix)+p0%F(1)      
              debugJ0( ix)=debugJ0( ix)+p0%F(10)    
           endif

        end do


!print*, debugT0
        write(fo_name,555) 'wave_',rank,iLv,istep
        open (90+rank, file=fo_name   , form=  'formatted')

!        do i=1, NXB
        do i=1,NXB*(2**LvMax)
           if(debugI0(i).gt.0) then 
              write(90+rank,*) i,real(debugE0(i)/dble(debugI0(i))) &
                           ,real(debugT0(i)/dble(debugI0(i))) &
                           ,real(debugJ0(i)/dble(debugI0(i)))
           end if
        end do
        close(90+rank)

        return 

555    format(A,I2.2,'_', I2.2,'_',I6.6,'.dat') 

        end subroutine output_wave
!-----------------------------------------------
      subroutine output_fourier
!-----------------------------------------------
        use mpi
        use param
        use oct_set
        use particle_set
        use init_mesh_size
        use const
        implicit none
        integer(kind=4)               :: iLv
        integer(kind=4)               :: i,index,ix,minindex,maxindex
        real(kind=8)                  :: di
        real(kind=8),dimension(NXB)   :: debugP0,debugT0
        integer(kind=4),dimension(NXB):: debugI0
! -- pointers --
        type(oct), pointer   :: p0
     

        iLv=0
        maxindex=maxID(1,iLv)
        minindex=minID(1,iLv)

        fourE=0.d0

        debugT0 = 0.d0
        debugP0 = 0.d0
        debugI0 = 0
     
        di = dx(1)

        do index=minindex,maxindex
           p0 => Mesh(index)
           ix = int((p0%rPOS(1)-R_lim(0,1))/dble(di))+1
              debugI0( ix)=debugI0( ix)+1
              debugT0( ix)=debugT0( ix)+p0%F(1)           !/dble(NYB*NZB)
        end do

        do i=1,NXB
           if(debugI0(i).gt.0) then 
              debugP0(i)=debugT0(i)/dble(debugI0(i))
           end if
        end do

        debugT0=debugP0
        call realft(debugT0,NXB,1)
        
        fourE(istep0)=sqrt(debugT0(5)**2+debugT0(6)**2)

!        write(fo_name,555) 'fourier_',iLv,istep
!        open (92, file=fo_name   , form=  'formatted')
!
!        do i=1,NXB 
!              write(92,*) i,debugP0(i),debugT0(i)
!        end do

!        close(92)

        return 
        end subroutine output_fourier

!-------------------------------------------
        subroutine copy_func(A,ii,jj,iLv)
!-------------------------------------------
        use param
        implicit none
        integer(kind=4)               :: i,j,ii,jj,iLv,iii,jjj
        real(kind=8),dimension(ii,jj) :: A,B

        iii=ii/(2**(Lvmax-iLv+1))
        jjj=jj/(2**(Lvmax-iLv+1))

        if(iLv.eq.1) then 
           do j=1,jjj
              do i=1,iii
                 B(2*i-1,2*j-1)=A(i,j)
                 B(  2*i,2*j-1)=A(i,j)
                 B(2*i-1,  2*j)=A(i,j)
                 B(  2*i,  2*j)=A(i,j)
              end do
           end do
        else if(iLv.eq.0) then 
           do j=1,jjj
              do i=1,iii
                 B(4*i-3,4*j-3)=A(i,j) ; B(4*i-2,4*j-3)=A(i,j)
                 B(4*i-1,4*j-3)=A(i,j) ; B(4*i  ,4*j-3)=A(i,j)
                 B(4*i-3,4*j-2)=A(i,j) ; B(4*i-2,4*j-2)=A(i,j)
                 B(4*i-1,4*j-2)=A(i,j) ; B(4*i  ,4*j-2)=A(i,j)
                 B(4*i-3,4*j-1)=A(i,j) ; B(4*i-2,4*j-1)=A(i,j)
                 B(4*i-1,4*j-1)=A(i,j) ; B(4*i  ,4*j-1)=A(i,j)
                 B(4*i-3,4*j  )=A(i,j) ; B(4*i-2,4*j  )=A(i,j)
                 B(4*i-1,4*j  )=A(i,j) ; B(4*i  ,4*j  )=A(i,j)
              end do
           end do
        else
           do j=1,jjj
              do i=1,iii
                 B(i,j)=A(i,j)
              end do
           end do
        end if
        A = B
        return 
        end subroutine copy_func
!
!-------------------------------------------
        subroutine copy_func0(A,ii,jj,iLv)
!-------------------------------------------
        use param
        implicit none
        integer(kind=4)                  :: i,j,ii,jj,iLv,iii,jjj
        integer(kind=4),dimension(ii,jj) :: A,B

        iii=ii/(2**(Lvmax-iLv+1))
        jjj=jj/(2**(Lvmax-iLv+1))

        if(iLv.eq.1) then 
           do j=1,jjj
              do i=1,iii
                 B(2*i-1,2*j-1)=A(i,j)
                 B(  2*i,2*j-1)=A(i,j)
                 B(2*i-1,  2*j)=A(i,j)
                 B(  2*i,  2*j)=A(i,j)
              end do
           end do
        else if(iLv.eq.0) then 
           do j=1,jjj
              do i=1,iii
                 B(4*i-3,4*j-3)=A(i,j) ; B(4*i-2,4*j-3)=A(i,j)
                 B(4*i-1,4*j-3)=A(i,j) ; B(4*i  ,4*j-3)=A(i,j)
                 B(4*i-3,4*j-2)=A(i,j) ; B(4*i-2,4*j-2)=A(i,j)
                 B(4*i-1,4*j-2)=A(i,j) ; B(4*i  ,4*j-2)=A(i,j)
                 B(4*i-3,4*j-1)=A(i,j) ; B(4*i-2,4*j-1)=A(i,j)
                 B(4*i-1,4*j-1)=A(i,j) ; B(4*i  ,4*j-1)=A(i,j)
                 B(4*i-3,4*j  )=A(i,j) ; B(4*i-2,4*j  )=A(i,j)
                 B(4*i-1,4*j  )=A(i,j) ; B(4*i  ,4*j  )=A(i,j)
              end do
           end do
        else
           do j=1,jjj
              do i=1,iii
                 B(i,j)=A(i,j)
              end do
           end do
        end if
        A = B
        return 
        end subroutine copy_func0


!*********************************************************************
!2011/12/05 yagi
subroutine output_params(POIntvl)
  use oct_set
  use param
  use message_passing_interface
  use init_mesh_size
  implicit none
  integer(kind=4)::iLv
  integer(kind=4),intent(in)::POIntvl
  integer(kind=4)::Pint
  if(debugMode>0)print *,"start output_params",rank
  
  Pint=POIntvl
  if(POIntvl<=0)Pint=PoutIntvl

  if(istep==0)then
     if(VHistOutput/=0)then
        !call output_histogram(HistMax,HistMin,HistDatanum,"Vx",0,0)
        !call output_histogram(HistMax,HistMin,HistDatanum,"Vy",0,0)
        !call output_histogram(HistMax,HistMin,HistDatanum,"Vz",0,0)
     endif
     if(PHistOutput/=0)then
        call output_histogram(R_lim(1,1),R_lim(0,1),HistDatanum,"Px",0,0)
        call output_histogram(R_lim(1,2),R_lim(0,2),HistDatanum,"Py",0,0)
        call output_histogram(R_lim(1,3),R_lim(0,3),HistDatanum,"Pz",0,0)
     endif
  endif
  if(VHistOutput/=0)then
!    call output_histogram2(HistMax,HistMin,HistDatanum,"Vx",0,0)
!    call output_histogram2(HistMax,HistMin,HistDatanum,"Vy",0,0)
!    call output_histogram2(HistMax,HistMin,HistDatanum,"Vz",0,0)
    call output_histogram(HistMax,HistMin,HistDatanum,"Vx",0,0)
    call output_histogram(HistMax,HistMin,HistDatanum,"Vy",0,0)
    call output_histogram(HistMax,HistMin,HistDatanum,"Vz",0,0)
  endif
!  if(Doutput/=0)call output_field_snap(1,"D",iLv,Goctout,OoctIntvl,outputQuality)
  do iLv=0,LvMax
     if(    Routput/=0)call output_field_snap(   1,"R",iLv,Goctout,OoctIntvl,outputQuality)
     if(    Eoutput/=0)call output_field_snap(   1,"E",iLv,Goctout,OoctIntvl,outputQuality)  
     if(    Joutput/=0)call output_field_snap(   1,"J",iLv,Goctout,OoctIntvl,outputQuality)
     if(    Boutput/=0)call output_field_snap(   1,"B",iLv,Goctout,OoctIntvl,outputQuality)
     if(    Zoutput/=0)call output_field_snap(   1,"Z",iLv,Goctout,OoctIntvl,outputQuality)
     if(    Poutput/=0)call output_field_snap(Pint,"P",iLv,Goctout,OoctIntvl,outputQuality)
     if(BackBoutput/=0)call output_field_snap(   1,"D",iLv,Goctout,OoctIntvl,outputQuality)
     if( eleJoutput/=0)call output_field_snap(   1,"X",iLv,Goctout,OoctIntvl,outputQuality)
     if( ionJoutput/=0)call output_field_snap(   1,"Y",iLv,Goctout,OoctIntvl,outputQuality)
     if(    Voutput/=0)call output_field_snap(   1,"V",iLv,Goctout,OoctIntvl,outputQuality)
     if(  ExBoutput/=0)call output_field_snap(   1,"A",iLv,Goctout,OoctIntvl,outputQuality)
     if(gradBoutput/=0)call output_field_snap(   1,"G",iLv,Goctout,OoctIntvl,outputQuality)
     if(Gridoutput /=0)call Display_Grid_All
     if(GMeshoutput/=0)call Display_GMesh_All
  enddo
  
  if(debugMode>0)print *,"exit output_params",rank
end subroutine output_params

subroutine output_traces
   use oct_set
  use param
  use message_passing_interface
  use init_mesh_size
  implicit none
  if(PTraceOutput/=0)call output_trace("P")
  if(NTraceOutput/=0)call output_trace("N")
  if(LTraceOutput/=0)call output_trace("L")
  if(XTraceOutput/=0)call output_trace("X")
  if(ETraceOutput/=0)call output_trace("E")
  if(RTraceOutput/=0)call output_trace("R")
  if(TTraceOutput/=0)call output_trace("T")
end subroutine output_traces

subroutine output_details
  use oct_set
  use param
  use message_passing_interface
  use init_mesh_size
  implicit none
 
  if(IDDetail/=0)call output_MinMaxID
  call full_oct_check
  if(ptclDetail/=0)call output_particle_num
  if(iFLGDetail/=0)call output_iFLG_num(1)

end subroutine output_details

!2011/08/11 added by yagi
subroutine output_field_snap(skipNumber,outType,iLv,outOct,octCut,outQuality)
  use oct_set
  use param
  use const
  use message_passing_interface
  use init_mesh_size
  implicit none
  integer(kind=4),intent(in)::skipNumber,iLv,outOct,octCut,outQuality
  integer(kind=4)::index,i,output_count,output_count2
  integer(kind=4)::iMax
  real(kind=8)::maxrPos(3),minrPos(3)
  real(kind=8)::mPos(3),pPos(3),oF(6),Bx,By,Bz,Ex,Ey,Ez,Bs
  character(LEN=1),intent(in)::outType
  character(LEN=128)::filename1,filename2,filename3,filename4
  type(oct),pointer::p0
  type(prtcl),pointer::pp

  output_count=0
  output_count2=0

  do i=1,3
    mPos(i)=center(i)-0.5d0*dx(i)
    pPos(i)=center(i)+0.5d0*dx(i)
    oF(i)=0.d0
  end do
  if(outOct/=0)then
     iMax=4
  else
     iMax=2
  endif

  if(outType/="E" .and. outType/="B" .and. outType/="J" .and. &
     outType/="P" .and. outType/="Z" .and. outType/="R" .and. &
     outType/="D" .and. outType/="X" .and. outType/="Y" .and. &
     outType/="V" .and. outType/="A" .and. outType/="G" )then
        print *,"debug outputFile creating. type=",outType,rank
  endif
  
  if(outType=="P".or.outType=="G" )then
    write(filename1,2012)version,outType,"e",istep,rank,iLv
    write(filename2,2012)version,outType,"i",istep,rank,iLv
    open(98,file=filename1,status='replace',form='formatted')
    open(99,file=filename2,status='replace',form='formatted')
  else if(outType=="V")then
    write(filename1,2012)version,outType,"e",istep,rank,iLv
    write(filename2,2012)version,outType,"i",istep,rank,iLv
    write(filename3,2012)version,"U","e",istep,rank,iLv
    write(filename4,2012)version,"U","i",istep,rank,iLv
    open(98,file=filename1,status='replace',form='formatted')
    open(99,file=filename2,status='replace',form='formatted')
    open(97,file=filename3,status='replace',form='formatted')
    open(96,file=filename4,status='replace',form='formatted')
  else
    write(filename1,2011)version,outType,istep,rank,iLv
    open(98,file=filename1,status='replace',form='formatted')
  end if

  if(outType=="E")then
     do i=1,iMax
        if(MinID(i,iLv)>=MaxID(i,iLv))cycle
        do index=MinID(i,iLv),MaxID(i,iLv)
           p0=>Mesh(index)
           if(p0%iFLG(1)>3 .or. p0%iFLG(1)<0)cycle
           oF(1)=p0%F(1)+FIMF(1)
           oF(2)=p0%F(2)+FIMF(2)
           oF(3)=p0%F(3)+FIMF(3)
           if     (outQuality==0)then
              write(98,*)p0%rPos,oF(1:3)
           else if(outQuality==1)then
              write(98,'(3f12.6,3f14.7)')p0%rPos,oF(1:3)
           endif
        enddo
     enddo
  elseif(outType=="B")then
     do i=1,iMax
        if(MinID(i,iLv)>=MaxID(i,iLv))cycle
        do index=MinID(i,iLv),MaxID(i,iLv)
           p0=>Mesh(index)
           if((p0%iFLG(1)<0).or.(p0%iFLG(1)>3))cycle
           oF(1)=p0%F(4)+FIMF(4)
           oF(2)=p0%F(5)+FIMF(5)
           oF(3)=p0%F(6)+FIMF(6)
           if     (outQuality==0)then
              write(98,*)p0%rPos,oF(1:3)
           else if(outQuality==1)then
              write(98,'(3f12.6,3f14.7)')p0%rPos,oF(1:3)
           endif
        enddo
     enddo
  elseif(outType=="D")then
     do i=1,iMax
        if(MinID(i,iLv)>=MaxID(i,iLv))cycle
        do index=MinID(i,iLv),MaxID(i,iLv)
            p0=>Mesh(index)
            if((p0%iFLG(1)<0).or.(p0%iFLG(1)>3))cycle
            oF(1)=p0%D(1)+p0%F(4)+FIMF(4)
            oF(2)=p0%D(2)+p0%F(5)+FIMF(5)
            oF(3)=p0%D(3)+p0%F(6)+FIMF(6)
           if     (outQuality==0)then
              write(98,*)p0%iPOS,oF(1:3)
           else if(outQuality==1)then
              write(98,'(3f12.6,3f14.7)')p0%rPos,oF(1:3)
           endif
        enddo
     end do
  elseif(outType=="J")then
     do i=1,iMax
        if(MinID(i,iLv)>=MaxID(i,iLv))cycle
        do index=MinID(i,iLv),MaxID(i,iLv)
           p0=>Mesh(index)
           if((p0%iFLG(1)<0).or.(p0%iFLG(1)>3))cycle
           if     (outQuality==0)then
              write(98,*)p0%rPos,p0%F(7:12)
           else if(outQuality==1)then
              write(98,'(3f12.6,6f14.7)')p0%rPos,p0%F(7:12)
           endif
        enddo
     enddo
  elseif(outType=="X")then !electron J
     do i=1,iMax
        if(MinID(i,iLv)>=MaxID(i,iLv))cycle
        do index=MinID(i,iLv),MaxID(i,iLv)
           p0=>Mesh(index)
           if((p0%iFLG(1)<0).or.(p0%iFLG(1)>3))cycle
           if     (outQuality==0)then
              write(98,*)p0%rPos,p0%F(13:15)
           else if(outQuality==1)then
              write(98,'(3f12.6,6f14.7)')p0%rPos,p0%F(13:15)
           endif
        enddo
     enddo
  elseif(outType=="Y")then !ion J
     do i=1,iMax
        if(MinID(i,iLv)>=MaxID(i,iLv))cycle
        do index=MinID(i,iLv),MaxID(i,iLv)
           p0=>Mesh(index)
           if((p0%iFLG(1)<0).or.(p0%iFLG(1)>3))cycle
           if     (outQuality==0)then
              write(98,*)p0%rPos,p0%F(16:18)
           else if(outQuality==1)then
              write(98,'(3f12.6,6f14.7)')p0%rPos,p0%F(16:18)
           endif
        enddo
     enddo
  elseif(outType=="Z")then
     do index=MinID(1,iLv),MaxID(1,iLv)
        if(MinID(1,iLv) .ge. MaxID(1,iLv))cycle
        p0=>Mesh(index)
        if(p0%iFLG(1)>3 .or. p0%iFLG(1)<0)cycle
        if     (outQuality==0)then
           write(98,*)p0%rPos,p0%Z(1:2)
        else if(outQuality==1)then
           write(98,'(3f12.6,2f14.7)')p0%rPOS,p0%Z(1:2)
        end if
     enddo
  elseif(outType=="A")then
     do i=1,iMax
        if(MinID(i,iLv)>=MaxID(i,iLv))cycle
        do index=MinID(i,iLv),MaxID(i,iLv)
           p0=>Mesh(index)
           if((p0%iFLG(1)<0).or.(p0%iFLG(1)>3))cycle
           Ex=p0%F(1)+FIMF(1)
           Ey=p0%F(2)+FIMF(2)
           Ez=p0%F(3)+FIMF(3)
           Bx=p0%F(4)+FIMF(4)+p0%D(1)
           By=p0%F(5)+FIMF(5)+p0%D(2)
           Bz=p0%F(6)+FIMF(6)+p0%D(3)
           Bs=Bx**2+By**2+Bz**2
           oF(1)=(Ey*Bz-Ez*By)/Bs
           oF(2)=(Ez*Bx-Ex*Bz)/Bs
           oF(3)=(Ex*By-Ey*Bx)/Bs
           if     (outQuality==0)then
              write(98,*)p0%rPos,oF(1:3)
           else if(outQuality==1)then
              write(98,'(3f12.6,3f14.7)')p0%rPos,oF(1:3)
           endif
        enddo
     enddo
  elseif(outType=="G")then !gradB drift
     do i=1,iMax
        if(MinID(i,iLv)>=MaxID(i,iLv))cycle
        do index=MinID(i,iLv),MaxID(i,iLv)
           p0=>Mesh(index)
           if((p0%iFLG(1)<0).or.(p0%iFLG(1)>3))cycle
           call calculate_gradB(index)
           if     (outQuality==0)then
              write(98,*)p0%rPos,p0%O(7:9)
              write(99,*)p0%rPos,p0%O(10:12)
           else if(outQuality==1)then
              write(98,'(3f12.6,6f14.7)')p0%rPos,p0%O(7:9)
              write(98,'(3f12.6,6f14.7)')p0%rPos,p0%O(7:9)
              write(99,'(3f12.6,6f14.7)')p0%rPos,p0%O(10:12)
           endif
        enddo
     enddo
  elseif(outType=="V")then  ! calculation velocity
     do index=MinID(1,iLv),MaxID(1,iLv)
        call calculate_V(index)
     end do
!     do index=MinID(1,iLv),MaxID(1,iLv)
!        call LPF(index)
!     end do
     do index=MinID(1,iLv),MaxID(1,iLv)
        if(MinID(1,iLv) .ge. MaxID(1,iLv))cycle
        p0=>Mesh(index)
        if(p0%iFLG(1)>3 .or. p0%iFLG(1)<0)cycle
!        call calculate_V2(index)
        if     (outQuality==0)then
           write(98,*)p0%rPos,p0%O(1:3)
           write(99,*)p0%rPos,p0%O(4:6)
           write(97,*)p0%rPos,p0%C(4:6)
           write(96,*)p0%rPos,p0%C(7:9)
!           write(97,*)p0%rPOS,p0%O(13:15)
!           write(96,*)p0%rPOS,p0%O(16)
        else if(outQuality==1)then
           write(98,'(3f12.6,3f14.7)')p0%rPOS,p0%O(1:3)
           write(99,'(3f12.6,3f14.7)')p0%rPOS,p0%O(4:6)
           write(97,'(3f12.6,3f14.7)')p0%rPOS,p0%C(4:6)
           write(96,'(3f12.6,3f14.7)')p0%rPOS,p0%C(7:9)
!           write(97,'(3f12.6,3f14.7)')p0%rPOS,p0%O(13:15)
!           write(96,'(3f12.6,3f14.7)')p0%rPOS,p0%O(2),p0%O(16),p0%C(2)
        end if
     enddo
  elseif(outType=="P")then!output Particle
     if(MaxID(1,iLv)>MinID(1,iLv))then
        do index=MinID(1,iLv),MaxID(1,iLv)
           p0=>Mesh(index)
           if(p0%iFLG(1)<0 .or. p0%iFLG(1)>3)cycle
           if(p0%octP==0)cycle
           
!           if(p0%rPOS(3)<mPos(3))cycle
!           if(p0%rPOS(3)>pPos(3))cycle
!           if(p0%rPOS(2)<mPos(2))cycle
!           if(p0%rPOS(2)>pPos(2))cycle
!           if(p0%rPOS(1)<mPos(1))cycle
!           if(p0%rPOS(1)>pPos(1))cycle

           pp=>p0%ptcl
           do i=1,p0%octP
              pp=>pp%prtnxt
              !if(pp%isort==2)output_count2=output_count2+1
              !if(pp%isort/=1)cycle
              output_count=output_count+1
              if(mod(output_count,skipNumber)/=0)cycle        
              
              if     (outQuality==0)then
                 if(pp%isort==1)then
                    write(98,*)pp%R(1:6)
                 else if(pp%isort==2)then
                    write(99,*)pp%R(1:6)
                 end if
              else if(outQuality==1)then
                 if(pp%isort==1)then
                    write(98,'(6f12.8)')pp%R(1:6)
                 else if(pp%isort==2)then
                    write(99,'(6f12.8)')pp%R(1:6)
                 end if
              endif
           enddo
        enddo
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierr)
     if(debugMode>=1)then
        print *,"In output_field_snap...pnum isort1=",output_count,"isort2=",output_count2,rank
        print *,"In output_field_snap...maxIP=",maxIP,"Pesh size=",size(Pesh)/(LvMax+2),rank
     endif
  elseif(outType=="R")then
     !make header file
     if(rank==0)then
        open(99,file="./output/Data/techHeader.dat",form="formatted")
        !set dx , minrPos , maxrPos ,param dimension, nprocs ,LvMax
        write(99,'(A)')"dx, param dimension, nprocs, LvMax, Nstep, DataIntvl, version"
        write(99,'(3f12.6,3i3,2i7,A)')dx(1:3),1,nprocs,LvMax,Nstep,DataIntvl,version
        close(99)
     endif
     maxrPos=0.d0
     minrPos=1.0d20
     output_count=0
     !get minmaxRPos
     if(minID(1,iLv)<maxID(1,iLv))then
        do index=MinID(1,iLv),MaxID(1,iLv)
           p0=>Mesh(index)
           maxrPos(1)=max(maxrPos(1),p0%rPos(1))
           maxrPos(2)=max(maxrPos(2),p0%rPos(2))
           maxrPos(3)=max(maxrPos(3),p0%rPos(3))
           minrPos(1)=min(minrPos(1),p0%rPos(1))
           minrPos(2)=min(minrPos(2),p0%rPos(2))
           minrPos(3)=min(minrPos(3),p0%rPos(3))
           if(p0%iFLG(1)>=0 .and. p0%iFLG(1)<=3)output_count=output_count+1
        enddo
        !set region size
        write(98,'(9f12.6,i9)')minrPos(1:3),maxrPos(1:3),R_lim_local(0,1:3),output_count
     endif
     if(minID(1,iLv)<maxID(1,iLv))then
        do index=MinID(1,iLv),MaxID(1,iLv)
           p0=>Mesh(index)
           if(p0%iFLG(1)<0 .or. p0%iFLG(1)>3)cycle
           write(98,'(4f12.6)')p0%rPos(1:3),p0%Z(1)
        enddo
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierr) 
  end if

  close(98)
  if(debugMode>=1)print *,"output_field_snap completed istep=",istep,rank
  return

  2011 format("./output/Data/",A,"_",A,i5.5,"_P",i3.3,"_Lv",i2.2,".dat")
  2012 format("./output/Data/",A,"_",A,A,i5.5,"_P",i3.3,"_Lv",i2.2,".dat")
end subroutine output_field_snap

subroutine output_particle_num
  use oct_set
  use param
  use message_passing_interface
  use init_mesh_size
  implicit none
  integer(kind=4)::index,i,j,k,iFLG,Isort,iLv,ii
  integer(kind=4),dimension(-4:6,-4:IonSorts,-1:LvMax)::Apnum,Gpnum
  type(oct),pointer::p0
  type(prtcl),pointer::pp
  
  !print *,"debug output_particle_num",rank

  Apnum=0
  Gpnum=0
  nullify(pp)
  
  !get particle num in all oct
  do iLv=0,LvMax
     do ii=1,2
        if(MinID(ii,iLv)>=MaxID(ii,iLv))cycle
        do index=MinID(ii,iLv),MaxID(ii,iLv)
           p0=>Mesh(index)
           !if(p0%iFLG(1)<=-4 .or. index>MaxID(2,LvMax))cycle

           iFLG=p0%iFLG(1)
     !----------------
           if(iFLG>6)iFLG=6
     !----------------

           if(associated(p0%ptcl))then
              pp=>p0%ptcl
              Isort=pp%isort
              Apnum(iFLG,Isort,iLv)=Apnum(iFLG,Isort,iLv)+1
        
              if(Isort/=-1)then
                 print *,"rep isort error index=",index,p0%octType,"rep isort=",pp%isort,rank
              endif
           else
              if(iFLG>-4)print *,"rep error index=",index,p0%octType,p0%iFLG,rank
           endif
  
           if(p0%octP>0)then
              do i=1,p0%octP
                 if(.not. associated(pp%prtnxt))then
                    print *,"AM p link error",index,rank
      
                    stop
                 endif
                 pp=>pp%prtnxt
                 Isort=pp%isort
                 if(Isort>4)then
!!$              print *,"illegal isort is detected",Isort,index,rank
!!$              stop
                    Isort=4
                 endif
                 Apnum(iFLG,Isort,iLv)=Apnum(iFLG,Isort,iLv)+1

              enddo
           endif
        enddo
     enddo
  enddo
  
  do iLv=0,LvMax
     do ii=3,4
        do index=MinID(ii,iLv),MaxID(ii,iLv)
           p0=>Mesh(index)
           !if(index>GMaxID(2,LvMaxG) .or. index<=GMaxID(1,-1))cycle
           iFLG=p0%iFLG(1)
           !------------------
           if(iFLG>6)iFLG=6
           !------------------

           if(associated(p0%ptcl))then
              pp=>p0%ptcl
              Isort=pp%isort
              Gpnum(iFLG,Isort,iLv)=Gpnum(iFLG,Isort,iLv)+1
           else
              if(iFLG>-4)print *,"rep Gerror index=",index,p0%octType,rank
           endif

           if(p0%octP>0)then 
              !print *,"this goct has ptcl",p0%octN,p0%octP,p0%iFLG,rank

              do i=1,p0%octP
                 if(.not. associated(pp%prtnxt))then
                    print *,"GM link error",index,rank

                    stop
                 endif
                 pp=>pp%prtnxt
                 Isort=pp%isort

                 if(Isort>4)then
                    print *,"illegal isort is detected in G",Isort,index,rank
                    stop
                 endif

                 Gpnum(iFLG,Isort,iLv)=Gpnum(iFLG,Isort,iLv)+1
              enddo
           endif
        enddo
     enddo
  enddo
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if(rank==0)print *,"------------------------------------"
  do k=0,nprocs-1
     if(rank==k)then
    
        print *,"Particle List in Mesh",rank
        do i=-4,4
           print '(A,I2,$)',"iFLG:",i
           do j=-1,IonSorts
              do iLv=0,LvMax
                 print '(A,I2,A,I2,A,I6,$)',"|Lv:",iLv," P(",j,"):",Apnum(i,j,iLv)
              enddo
           enddo
           print '(A)',"|"
        enddo

        print *,"Particle List in GMesh",rank
        do i=-4,4
           print '(A,I2,$)',"iFLG",i
           do j=-1,IonSorts
              do iLv=0,Lvmax
                 print '(A,I2,A,I2,A,I6,$)',"|Lv:",iLv," P(",j,"):",Gpnum(i,j,iLv)
              enddo
           enddo
           print '(A)',"|"
        enddo
     
     endif
    
     call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo
  if(rank==0)print *,"------------------------------------"
  !print *,"debug output_particle_num end",rank

end subroutine output_particle_num

subroutine output_iFLG_num(flag_num)
  use oct_set
  use param
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::index,temp,i,j,iLv
  integer(kind=4),intent(in)::flag_num
  integer(kind=4),dimension(-1:LvMax,-4:8)::iFLG
  type(oct),pointer::p0

  if(flag_num<0 .or. flag_num>3)then
     print *,"flag_num error:",flag_num,"in output_iFLG_num",rank
     return
  endif

  iFLG=0
  do iLv=-1,LvMax
     do i=1,2
        if(MinID(i,iLv)<MaxID(i,iLv))then
        
           do index=MinID(i,iLv),MaxID(i,iLv)
              p0=>Mesh(index)
              if(p0%iFLG(1)/=4 .and. iLv==-1)print *,"BMEsh iFLG error:",index,p0%iFLG,rank
              temp=p0%iFLG(flag_num)
              iFLG(iLv,temp)=iFLG(iLv,temp)+1
           enddo
         
        endif
     enddo
  enddo

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  do j=0,nprocs-1
     if(rank==j)then
        print *,"-----------------------------"
        print *,"Mesh iFLG count",rank
        do iLv=-1,LvMax
           print '(A,I1,$)',"Lv=",iLv
           do i=-4,6
              print '(A,I2,A,I7,$)'," F(",i,"):",iFLG(iLv,i)
              !if(mod(i,4)==0)print *,""
           enddo
           print *,""
        enddo
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo
  
  iFLG=0
  do iLv=-1,LvMax
     do i=3,4
        if(MinID(i,iLv)<MaxID(i,iLv))then
           do index=MinID(i,iLv),MaxID(i,iLv)
              p0=>Mesh(index)
              temp=p0%iFLG(flag_num)
              iFLG(iLv,temp)=iFLG(iLv,temp)+1
           enddo
        endif
     enddo
  enddo

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  do j=0,nprocs-1
     if(rank==j)then
        print *,"-----------------------------"
        print *,"GMesh iFLG count rank:",rank
        do iLv=-1,LvMax
           print '(A,I2,$)',"Lv=",iLv
           do i=-4,6
              print '(A,I2,A,I7,$)'," F(",i,"):",iFLG(iLv,i)
              !if(mod(i,4)==0)print *,""
           enddo
           print *,""
        enddo
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo
  if(rank==nprocs-1)print *,"-----------------------------"
  call mpi_barrier(MPI_COMM_WORLD,ierr)
end subroutine output_iFLG_num

subroutine output_MinMaxID
  use oct_set
  use param
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::i,iLv,temp,j
  !output Mesh activity ratio
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  do i=0,nprocs-1
     if(i==rank)then
        print *       ,"========================================="
        print '(A,I3)',"      Mesh MinMaxID Check Rank:",rank
        print '(A,I7,A,I7)',"    MeshSize=",MeshBound," GMeshSize=",size(Mesh)-MeshBound
        print *       ,"========================================="
        temp=0
        do iLv=-1,LvMax
           do j=1,2
              print "(A,I1,A,I2,A,I7,A,$)","MinID(",j,",",iLv,")=",MinID(j,iLv)," | "
              print "(A,I1,A,I2,A,I7,A,$)","MaxID(",j,",",iLv,")=",MaxID(j,iLv)," | "
              temp=temp+1
              if(temp>=2)then
                 temp=0
                 print *,""
              endif
           enddo
        enddo
        print *       ,"========================================="
        temp=0
        
        do iLv=-1,LvMax
           do j=3,4
              print "(A,I1,A,I2,A,I7,A,$)","GMinID(",j,",",iLv,")=",MinID(j,iLv)," | "
              print "(A,I1,A,I2,A,I7,A,$)","GMaxID(",j,",",iLv,")=",MaxID(j,iLv)," | "
              temp=temp+1
              if(temp>=2)then
                 temp=0
                 print *,""
              endif
           enddo
        enddo
        print *,"========================================="
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo

end subroutine output_MinMaxID

subroutine output_trace(outType)
  use oct_set
  use time_evolution
  use param
  use message_passing_interface
  use init_mesh_size
  use const
  implicit none
  integer(kind=4)::index,i,iLv,output_count,LvSize,output_count_i,output_count_e
  integer(kind=4),dimension(3)::iPOSH,iPOSN
  integer(kind=4),dimension(:),allocatable::all_pnum,all_pnum_i,all_pnum_e
  integer(kind=4),dimension(4,-1:LvMax)::myMeshNum,allMeshNum

  real(kind=8)::kinePsum,kineP,EneEsum,EneE,EneB,EneBsum,JpX,JpY,JpZ,AllEne
  real(kind=8),dimension(0:LvMax)::Octvol
  real(kind=8),dimension(:),allocatable::AllProcEne,AllProcE,AllProcP

  character(LEN=1),intent(in)::outType
  character(LEN=128)::filename
  type(oct),pointer::p0
  type(prtcl),pointer::pp
  
  if(outType/="N" .and. outType/="X" .and. outType/="J" .and. &
     outType/="P" .and. outType/="E" .and. outType/="L" .and. &
     outType/="R" .and. outType/="T")then
     if(outType/="I" .and. outType/="H")then
        print *,"outType is wrong ...",outType,rank
        stop
     else
        print *,"debug outputFile creating. type=",outType,rank
     endif
  endif

  !----Time----
  if(outType=="T")then
     call mpi_barrier(MPI_COMM_WORLD,ierr)
     if(rank==0)then
        write(filename,2513)version,outType,nprocs
        if((InitialOutput/=0 .and. istep==0) .or. (InitialOutput==0 .and. istep==1))then
           open(97,file=filename,status='replace',form='formatted',position='append')
        else
           open(97,file=filename,status='old'    ,form='formatted',position='append')
        endif

        write(97,*)MPI_Wtime()-wtime
        print *,"Elapse time of this simulation is ",MPI_Wtime()-wtime
        print *,"Elapse time of this step is ",MPI_Wtime()-steptime
        close(97)
     endif
     
     steptime=MPI_Wtime()
  endif

  !----E----
  if(outType=="E")then
     LvSize=NYB*(2**0)/2
 
     do i=0,nprocs-1
        if(rank==i)then
           write(filename,2510)version,outType

           if(((InitialOutput/=0 .and. istep==0) .or. (InitialOutput==0 .and. istep==1)) .and. rank==0)then
              open(96,file=filename,status='replace',form='formatted',position='append')
           else
              open(96,file=filename,status='old'    ,form='formatted',position='append')
           endif

           do index=MinID(1,0),MaxID(1,0)
              p0=>Mesh(index)
              !if(p0%iFLG(1)>6 .or. p0%iFLG(1)<-3)cycle

              iPOSH = p0%iPos
              iPOSN = (iPOSH+sConst(0))/intNxt(0)
              if(iPOSN(2)/=LvSize .or. iPOSN(3)/=LvSize)cycle
              write(96,*)p0%rPOS(1),time,p0%F(1:6)
           enddo

           if(rank==nprocs-1)write(96,*)" "
           close(96)
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)
     enddo
  endif

  !----P----
  if(outType=="P")then
     output_count=0
     do iLv=0,LvMax
        write(filename,2512)version,outType,iLv,rank
        
 
        if((InitialOutput/=0 .and. istep==0) .or. (InitialOutput==0 .and. istep==1))then
           open(97,file=filename,status='replace',form='formatted',position='append')
        else
           open(97,file=filename,status='old'    ,form='formatted',position='append')
        endif
  

        do index=MinID(1,iLv),MaxID(1,iLv)
           p0=>Mesh(index)
           if(p0%octP==0)cycle

           pp=>p0%ptcl
           do i=1,p0%octP
              pp=>pp%prtnxt
              if(pp%isort/=1)cycle
              output_count=output_count+1
              write(97,*)pp%R(1:3)
           enddo
        enddo

        write(97,*)""
        close(97)
     enddo
   
  endif
  !----L----
  if(outType=="L")then

     call collect_ptcl_loops(output_count)!dummy variable

     output_count=0
     allocate(all_pnum(1:nprocs))
     all_pnum=0

     write(filename,2513)version,outType,nprocs

     if(rank==0)then
        if((InitialOutput/=0 .and. Ltrace_num==0) .or. (InitialOutput==0 .and. Ltrace_num==0))then
           open(98,file=filename,status='replace',form='formatted',position='append')
        else
           open(98,file=filename,status='old'    ,form='formatted',position='append')
        endif
     endif
     do iLv=0,LvMax
        if(MaxID(1,iLv)>MinID(1,iLv))then
           do index=MinID(1,iLv),MaxID(1,iLv)
              p0=>Mesh(index)
              
              if(p0%iFLG(1)<0 .or. p0%iFLG(1)>4)cycle
            
              output_count=output_count+p0%ptcl_loops
           enddo
        endif
     enddo

     !print *,"output_count=",output_count
     all_pnum=0
     call mpi_gather(output_count,1,MPI_INTEGER,all_pnum,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     if(rank==0)then
        print *,"Summation of particle loops =",sum(all_pnum)
        write(98,'(i5.5,A,i9.9,$)')istep," ",sum(all_pnum)
        if(nprocs<10)then
           do i=1,nprocs
              write(98,'(A,i9.9,$)')" ",all_pnum(i)
           enddo
        else
           write(98,'(A,i9.9,$)')" ",maxval(all_pnum)
        endif
        write(98,*)
     endif

     Ltrace_num=Ltrace_num+1
     close(98)
     deallocate(all_pnum)

  endif
  !----R----
  if(outType=="R")then
     output_count=0
     myMeshNum=0
     allMeshNum=0

     write(filename,2511)version,outType,rank

     if(rank==0)then
        if((InitialOutput/=0 .and. istep==0) .or. (InitialOutput==0 .and. istep==1))then
           open(98,file=filename,status='replace',form='formatted',position='append')
        else
           open(98,file=filename,status='old'    ,form='formatted',position='append')
        endif
     endif
     do iLv=-1,LvMax
        do i=1,4
           if(MaxID(i,iLv)>MinID(i,iLv))then
              myMeshNum(i,iLv)=MaxID(i,iLv)-MinID(i,iLv)+1
           else
              myMeshNum(i,iLv)=0
           endif
        enddo
     enddo

     call mpi_reduce(myMeshNum,allMeshNum,(4*(LvMax+2)),MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     if(rank==0)then
        print *,"********* Summation of cell number *********"
        do iLv=-1,LvMax
           do i=1,4
              print '(A,i1,A,i3,A,i7,$)'," i=",i," iLv=",iLv," cell num=",allMeshNum(i,iLv)
           enddo
           print *,""
        enddo
        print *,"*******************************************"
     endif

     call mpi_barrier(MPI_COMM_WORLD,ierr)
     close(98)

  endif
  !----N----
  if(outType=="N")then
     output_count_i=0
     output_count_e=0
     allocate(all_pnum_i(1:nprocs))
     allocate(all_pnum_e(1:nprocs))
     all_pnum_i=0
     all_pnum_e=0

     write(filename,2513)version,outType,nprocs

     if(rank==0)then
        if((InitialOutput/=0 .and. Ntrace_num==0) .or. (InitialOutput==0 .and. Ntrace_num==0))then
           open(98,file=filename,status='replace',form='formatted',position='append')
        else
           open(98,file=filename,status='old'    ,form='formatted',position='append')
        endif
     endif

     do iLv=0,LvMax
        if(MaxID(1,iLv)>MinID(1,iLv))then
           do index=MinID(1,iLv),MaxID(1,iLv)
              p0=>Mesh(index)
              if(p0%octP==0)cycle
              !if(p0%iFLG(1)<0 .or. p0%iFLG(1)>4)cycle
              pp=>p0%ptcl
              do i=1,p0%octP
                 pp=>pp%prtnxt
                 if(pp%isort<=0)cycle
                 if (pp%isort==1) then
                   output_count_e=output_count_e+1
                 else if (pp%isort==2) then
                   output_count_i=output_count_i+1
                 endif
              enddo
           enddo
        endif
     enddo

     !print *,"output_count=",output_count_i+output_count_e
     all_pnum_i=0
     all_pnum_e=0
     call mpi_gather(output_count_i,1,MPI_INTEGER,all_pnum_i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call mpi_gather(output_count_e,1,MPI_INTEGER,all_pnum_e,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     if(rank==0)then
        print *,"Summation of particle number =",sum(all_pnum_i)
        write(98,'(i5.5,A,i9.9,A,i9.9,$)')istep," ",sum(all_pnum_i)," ",sum(all_pnum_e)
        if(nprocs<10)then
           do i=1,nprocs
              write(98,'(A,i9.9,A,i9.9,$)')" ",all_pnum_i(i)," ",all_pnum_e(i)
           enddo
        endif
        write(98,*)
     endif

     close(98)
     deallocate(all_pnum_i)
     deallocate(all_pnum_e)
     Ntrace_num=Ntrace_num+1

  endif

  !----X----
  if(outType=="X")then

     do i=0,LvMax
        Octvol(i)=(dx(1)*(HALF**i))*(dx(2)*(HALF**i))*(dx(3)*(HALF**i))
     enddo

     output_count=0
     allocate(AllProcEne(1:nprocs))
     allocate(AllProcE(1:nprocs))
     allocate(AllProcP(1:nprocs))
 
     write(filename,2510)version,outType

     if(rank==0)then
        if((InitialOutput/=0 .and. istep==0) .or. (InitialOutput==0 .and. istep==1))then
           open(99,file=filename,status='replace',form='formatted',position='append')
        else
           open(99,file=filename,status='old'    ,form='formatted',position='append')
        endif
     endif

     EneEsum=0.d0
     EneBsum=0.d0
     kinePsum=0.d0
     JpX=0.d0
     JpY=0.d0
     JpZ=0.d0

     do iLv=0,LvMax
        do index=MinID(1,iLv),MaxID(1,iLv)
           p0=>Mesh(index)
           if(p0%iFLG(1)>3 .or. p0%iFLG(1)<0)cycle


           kineP=0.0d0
           EneE=(p0%F(1)**2+p0%F(2)**2+p0%F(3)**2)*Octvol(p0%octLv)
           EneB=(p0%F(4)**2+p0%F(5)**2+p0%F(6)**2)*Octvol(p0%octLv)

           JpX=JpX+p0%F(7)*(0.5d0**iLv)
           JpY=JpY+p0%F(8)*(0.5d0**iLv)
           JpZ=JpZ+p0%F(9)*(0.5d0**iLv)

           EneEsum=EneEsum+EneE
           EneBsum=EneBsum+EneB

           if(p0%octP==0)cycle
           !if(p0%iFLG(1)>=4 .or. p0%iFLG(1)<0)cycle
           pp=>p0%ptcl
           do i=1,p0%octP
              pp=>pp%prtnxt
              if(pp%isort<=0)cycle
              kineP=kineP+Rmass(pp%isort)*(pp%R(4)**2+pp%R(5)**2+pp%R(6)**2)
           enddo

           kinePsum=kinePsum+kineP
      
        enddo
     enddo
     
     JpX=JpX/(NXR*NYR*NZR*NXB*NYB*NZB)
     JpY=JpY/(NXR*NYR*NZR*NXB*NYB*NZB)
     JpZ=JpZ/(NXR*NYR*NZR*NXB*NYB*NZB)

     kinePsum=kinePsum/(omega**2)/dble(npart_per_cell)*dx(1)*dx(2)*dx(3) !change 2013/12/27

     AllEne=kinePsum+EneEsum+EneBsum
     if(debugMode>=1)print *,"kineP=",kinePsum,"Eene=",EneEsum,"Bene=",EneBsum,rank
     !print *,"Jp=",JpX,JpY,JpZ,rank
     EneEsum=EneEsum+EneBsum
     !print *,"output_count=",output_count
     call mpi_gather(AllEne,1,MPI_DOUBLE_PRECISION,AllProcEne,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call mpi_gather(EneEsum,1,MPI_DOUBLE_PRECISION,ALLProcE,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call mpi_gather(kinePsum,1,MPI_DOUBLE_PRECISION,ALLProcP,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     if(rank==0)then
        print *,"Summation of energy =",sum(ALLProcEne)
        write(99,'(I5,3f12.6,3f12.6)')istep,sum(ALLProcEne),sum(ALLProcE),sum(ALLProcP),JpX,JpY,JpZ
        close(99)
     endif
    
     deallocate(AllProcEne)
     deallocate(AllProcE)
     deallocate(AllProcP)
  endif
 
  2510 format("./output/",A,"_Trace_",A,".dat")
  2511 format("./output/",A,"_Trace_",A,"_P",i3.3,".dat")
  2512 format("./output/",A,"_Trace_",A,"_Lv",i3.3,"_P",i3.3,".dat")
  2513 format("./output/",A,"_Trace",A,"_Proc",i6.6,".dat")
end subroutine output_trace

!**********************************************************************
!2012/02/23 yagi
subroutine output_histogram(Hmax,Hmin,HDnum,tgt,iLvMax,iLvMin)
  use oct_set
  use particle_set
  use init_mesh_size
  use const
  use param
  use message_passing_interface
  implicit none
  integer(kind=4),intent(in)::HDnum,iLvMax,iLvMin
  real(kind=8),intent(in)::Hmax,Hmin
  character(LEN=2),intent(in)::tgt
  integer(kind=4)::tgt_num
  character(LEN=128)::filename
  type(oct),pointer::p0
  type(prtcl),pointer::pp
  integer(kind=4)::iLv,maxindex,minindex,index,i,octP
  integer(kind=4),dimension(HDnum)::Histogram,Htemp
  real(kind=8)::HistD,HistDrev
  real(kind=8)::hist_d
  real(kind=8)::dataMax,dataMin
  integer(kind=8)::dist

  if(iLvMax>LvMax)then
     print *,"Argument error iLvMax",iLvMax
     return
  elseif(iLvMin<0)then
     print *,"Argument error iLvMin",iLvMin
     return
  endif

  if(tgt/="Vx" .and. tgt/="Vy" .and. tgt/="Vz" .and. &
     tgt/="Px" .and. tgt/="Py" .and. tgt/="Pz")then
     print *,"Unknown output type ",tgt
     return
  endif

  if(tgt=="Vx")tgt_num=4
  if(tgt=="Vy")tgt_num=5
  if(tgt=="Vz")tgt_num=6
  if(tgt=="Px")tgt_num=1
  if(tgt=="Py")tgt_num=2
  if(tgt=="Pz")tgt_num=3

  if(Hmax<=Hmin)then
     print *,"Illegal histogram range Hmax=",Hmax,"Hmin=",Hmin
     return
  endif

  HistD=(Hmax-Hmin)/dble(HDnum)
  HistDrev=1.d0/HistD
  Histogram=0
  Htemp=0
  dataMax=-999
  dataMin=999

  if(tgt=="Vx" .or. tgt=="Vy" .or. tgt=="Vz" .or. tgt=="Px" .or. tgt=="Py" .or. tgt=="Pz")then
     do iLv=iLvMin,iLvMax
        minindex=MinID(1,iLv)
        maxindex=MaxID(1,iLv)
        if(minindex>=maxindex)cycle
        do index=minindex,maxindex
           p0=>Mesh(index)
           if(p0%iFLG(1)<0 .or. p0%iFLG(1)>3)cycle
           if(p0%octP==0)cycle
           octP=p0%octP
           pp=>p0%ptcl
           do i=1,octP
              pp=>pp%prtnxt
              if(pp%isort/=1)cycle
              hist_d=(pp%R(tgt_num)-Hmin)*HistDrev
              dist=int(hist_d)
              if(dist<1)dist=1
              if(dist>HDnum)dist=HDnum
              Histogram(dist)=Histogram(dist)+1
              dataMax=max(dataMax,pp%R(tgt_num))
              dataMin=min(dataMin,pp%R(tgt_num))
           end do
        end do
     enddo
  endif

!sum data using mpi
  call mpi_reduce(Histogram,Htemp,HDnum,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  Histogram=Htemp
  if(rank==0)then
     write(filename,2513)version,tgt,istep
     open(77,file=filename,form='formatted')
     do i=1,HDnum
        write(77,*)Hmin+HistD*dble(i),Histogram(i)    
     enddo
     close(77)
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if(debugMode>0)print *,"output_histogram ",tgt," completed ... max=",dataMax,"min=",dataMin
  return 

  2513 format("./output/",A,"_Hist_",A,i5.5,".dat")
end subroutine output_histogram

!***********************************************************************
! +-------------------------------------------------------------------+
! |     create octNb connections finer Mesh cells
! |     that MPI-communicate with GMesh                               |
! +-------------------------------------------------------------------+
!***********************************************************************
!   
!   The following is a conceptual diagram for this connection
!
!            |------------|---------------|---------------|---------------|
!   iLv=-1   |     G      |       B       |       B       |       B       |
!            |            |               |               |               |
!   prc_bndry|     0      |       0       |       0       |       0       |
!
!            |------------|---------------|---------------|---------------|
!   iLv=0    |     G      |       M       |       M       |       M       |
!            |            |               |               |               |
!   prc_bndry|     1      |       1       |       1       |       0       |
!
!            |------------|---------------|---------------|---------------|
!   iLv=1    |     G      |   M   |   M   |   M   |   M   |
!            |            |       |       |       |       |
!   prc_bndry|     1      |   1   |   1   |   1   |   2   |
!
!                       <=a=>   <=b=>   <=c=>   <=b=>   <=d=>
!
!   As shown in this diagram, there are 4 types (from (a) to (d) of octNb 
!                   connections in connect_oct_for_make_GMesh.
!   (a) for GMesh and Mesh
!   (b) for internal connection
!   (c) for connecting with adjacent Mesh cells
!   (d) for connecting at the edge, where there (might) not be children. In this
!       case, connect to itself 
!

subroutine connect_G_n_M
  use oct_set
  use param
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)    :: iLv,i,index
  integer(kind=4)    :: ich
  type(oct), pointer :: p0

!---------------- Connect Nb pointers b/w GMesh and Mesh -------------------------
!---------------------------- Extracted from connect_oct -------------------------- 


!Connecting octNBs within GMesh and GMesh => Mesh (But not Mesh => GMesh!)
  iLv=-1
  !do iLv = -1, LvMaxG-1
!!$omp parallel do private(index)
  do index=MinID(3,iLv),maxID(3,iLv)
     p0 => Mesh(index)
     ! - for octNb1 -
     if(associated(p0%octNb1))then  ! This means that the cell is not at the edge.
        p0%octCh1%octNb1 => p0%octNb1%octCh2     
        p0%octCh3%octNb1 => p0%octNb1%octCh4   
        p0%octCh5%octNb1 => p0%octNb1%octCh6
        p0%octCh7%octNb1 => p0%octNb1%octCh8
     else 
        p0%octCh1%octNb1 => p0%octCh1 !These corresponds to the edges
        p0%octCh3%octNb1 => p0%octCh3 !In case there are no neighbour, connect to themself
        p0%octCh5%octNb1 => p0%octCh5   
        p0%octCh7%octNb1 => p0%octCh7 
     endif
     ! - for octNb2 -                
     if(associated(p0%octNb2))then
        p0%octCh2%octNb2 => p0%octNb2%octCh1
        p0%octCh4%octNb2 => p0%octNb2%octCh3
        p0%octCh6%octNb2 => p0%octNb2%octCh5
        p0%octCh8%octNb2 => p0%octNb2%octCh7
     else
        p0%octCh2%octNb2 => p0%octCh2
        p0%octCh4%octNb2 => p0%octCh4
        p0%octCh6%octNb2 => p0%octCh6
        p0%octCh8%octNb2 => p0%octCh8
     endif
     ! - for octNb3 -
     if(associated(p0%octNb3))then
        p0%octCh1%octNb3 => p0%octNb3%octCh3
        p0%octCh2%octNb3 => p0%octNb3%octCh4
        p0%octCh5%octNb3 => p0%octNb3%octCh7
        p0%octCh6%octNb3 => p0%octNb3%octCh8
     else
        p0%octCh1%octNb3 => p0%octCh1 
        p0%octCh2%octNb3 => p0%octCh2
        p0%octCh5%octNb3 => p0%octCh5
        p0%octCh6%octNb3 => p0%octCh6
     endif
     ! - for octNb4 -
     if(associated(p0%octNb4))then
        p0%octCh3%octNb4 => p0%octNb4%octCh1
        p0%octCh4%octNb4 => p0%octNb4%octCh2
        p0%octCh7%octNb4 => p0%octNb4%octCh5
        p0%octCh8%octNb4 => p0%octNb4%octCh6
     else
        p0%octCh3%octNb4 => p0%octCh3
        p0%octCh4%octNb4 => p0%octCh4
        p0%octCh7%octNb4 => p0%octCh7
        p0%octCh8%octNb4 => p0%octCh8
     endif
     ! - for octNb5 -
     if(associated(p0%octNb5))then
        p0%octCh1%octNb5 => p0%octNb5%octCh5
        p0%octCh2%octNb5 => p0%octNb5%octCh6
        p0%octCh3%octNb5 => p0%octNb5%octCh7
        p0%octCh4%octNb5 => p0%octNb5%octCh8
     else
        p0%octCh1%octNb5 => p0%octCh1 
        p0%octCh2%octNb5 => p0%octCh2
        p0%octCh3%octNb5 => p0%octCh3
        p0%octCh4%octNb5 => p0%octCh4
     endif
     ! - for octNb6 -
     if(associated(p0%octNb6))then
        p0%octCh5%octNb6 => p0%octNb6%octCh1
        p0%octCh6%octNb6 => p0%octNb6%octCh2
        p0%octCh7%octNb6 => p0%octNb6%octCh3
        p0%octCh8%octNb6 => p0%octNb6%octCh4
     else
        p0%octCh5%octNb6 => p0%octCh5
        p0%octCh6%octNb6 => p0%octCh6
        p0%octCh7%octNb6 => p0%octCh7
        p0%octCh8%octNb6 => p0%octCh8
     endif
     !-internal connection
     p0%octCh1%octNb2 => p0       %octCh2
     p0%octCh1%octNb4 => p0       %octCh3
     p0%octCh1%octNb6 => p0       %octCh5

     p0%octCh2%octNb1 => p0       %octCh1
     p0%octCh2%octNb4 => p0       %octCh4
     p0%octCh2%octNb6 => p0       %octCh6

     p0%octCh3%octNb2 => p0       %octCh4
     p0%octCh3%octNb3 => p0       %octCh1
     p0%octCh3%octNb6 => p0       %octCh7

     p0%octCh4%octNb1 => p0       %octCh3
     p0%octCh4%octNb3 => p0       %octCh2
     p0%octCh4%octNb6 => p0       %octCh8

     p0%octCh5%octNb2 => p0       %octCh6
     p0%octCh5%octNb4 => p0       %octCh7
     p0%octCh5%octNb5 => p0       %octCh1

     p0%octCh6%octNb1 => p0       %octCh5           
     p0%octCh6%octNb4 => p0       %octCh8
     p0%octCh6%octNb5 => p0       %octCh2

     p0%octCh7%octNb2 => p0       %octCh8
     p0%octCh7%octNb3 => p0       %octCh5
     p0%octCh7%octNb5 => p0       %octCh3

     p0%octCh8%octNb1 => p0       %octCh7
     p0%octCh8%octNb3 => p0       %octCh6
     p0%octCh8%octNb5 => p0       %octCh4
  enddo
!!$omp end parallel do
        !end do

!-----------------------------Next, Connect Mesh=> GMesh-------------------------
!Connecting Mesh => GMesh. Connecting through BMesh.

  !if(debugMode>=3)print*,'rank,MinID,MaxID',rank,MinID(1,0),MaxID(1,0),MinID(2,0),MaxID(2,0)
  do i=1,2
     !do iLv = 0,LvMax
     iLv=0
        if(MinID(i,iLv) >= MaxID(i,iLv)) cycle
        do index=minID(i,iLv),MaxID(i,iLv)


           p0 => Mesh(index)
           ich = p0 % Csort
           
           select case(ich)
           case(1)           
              p0%octNb1 => p0%octPrt%octNb1%octCh2  
           case(3)
              p0%octNb1 => p0%octPrt%octNb1%octCh4  
           case(5)
              p0%octNb1 => p0%octPrt%octNb1%octCh6  
           case(7)
              p0%octNb1 => p0%octPrt%octNb1%octCh8  
           end select
           !              endif
           !             if(.not. associated(p0%octNb2).and. p0%octPrt%octNb2%octType == 1)then
           select case(ich)
           case(2)

              p0%octNb2 => p0%octPrt%octNb2%octCh1  
           case(4)
              p0%octNb2 => p0%octPrt%octNb2%octCh3  
           case(6)
              p0%octNb2 => p0%octPrt%octNb2%octCh5  
           case(8)
              p0%octNb2 => p0%octPrt%octNb2%octCh7  
           end select
!              endif
!             if(.not. associated(p0%octNb3).and. p0%octPrt%octNb3%octType == 1)then
           select case(ich)
           case(1)
              p0%octNb3 => p0%octPrt%octNb3%octCh3  
           case(2)
              p0%octNb3 => p0%octPrt%octNb3%octCh4  
           case(5)
              p0%octNb3 => p0%octPrt%octNb3%octCh7  
           case(6)
              p0%octNb3 => p0%octPrt%octNb3%octCh8  
           end select
           !              endif
           !             if(.not. associated(p0%octNb4).and. p0%octPrt%octNb4%octType == 1)then
           select case(ich)
           case(3)
              p0%octNb4 => p0%octPrt%octNb4%octCh1  
           case(4)
              p0%octNb4 => p0%octPrt%octNb4%octCh2  
           case(7)
              p0%octNb4 => p0%octPrt%octNb4%octCh5  
           case(8)
              p0%octNb4 => p0%octPrt%octNb4%octCh6  
           end select
!              endif
!             if(.not. associated(p0%octNb5).and. p0% octPrt% octNb5%octType == 1)then
           select case(ich)
           case(1)
              p0%octNb5 => p0%octPrt%octNb5%octCh5  
           case(2)
              p0%octNb5 => p0%octPrt%octNb5%octCh6  
           case(3)
              p0%octNb5 => p0%octPrt%octNb5%octCh7  
           case(4)
              p0%octNb5 => p0%octPrt%octNb5%octCh8  
           end select
           !              endif
           !             if(.not. associated(p0%octNb6).and. p0%octPrt% octNb6%octType == 1)then
           select case(ich)
           case(5)
              p0%octNb6 => p0%octPrt%octNb6%octCh1  
           case(6)
              p0%octNb6 => p0%octPrt%octNb6%octCh2  
           case(7)
              p0%octNb6 => p0%octPrt%octNb6%octCh3  
           case(8)
              p0%octNb6 => p0%octPrt%octNb6%octCh4  
           end select
!              endif
        enddo
     !enddo
  enddo
  !if(debugMode>=1)print *,"connect_G_n_M comp"
  return
     
end subroutine connect_G_n_M

subroutine detach_G_n_M
  ! This subroutine disconnects octNb connection B/W GMesh(iLv>=0) and Mesh 
  use oct_set
  use param
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)    :: iLv,index
  type(oct), pointer :: p0

  do iLv=0,LvMax
     do index = MinID(3,iLv),MaxID(3,iLv)
        p0 => Mesh(index)
        if(p0%iFLG(1)<=-4)cycle
        if(associated(p0%octNb1))then
           if( p0%octNb1%octType==0) then !if the octType is
              nullify(p0%octNb1%octNb2)
              nullify(p0%octNb1)
           endif
        endif
        if(associated(p0%octNb2))then
           if( p0%octNb2%octType==0) then
              nullify(p0%octNb2%octNb1)
              nullify(p0%octNb2)
           endif
        endif
        if(associated(p0%octNb3))then
           if( p0%octNb3%octType==0) then
              nullify(p0%octNb3%octNb4)
              nullify(p0%octNb3)
           endif
        endif
        if(associated(p0%octNb4))then
           if( p0%octNb4%octType==0) then
              nullify(p0%octNb4%octNb3)
              nullify(p0%octNb4)
           endif
        endif
        if(associated(p0%octNb5))then
           if( p0%octNb5%octType==0) then
              nullify(p0%octNb5%octNb6)
              nullify(p0%octNb5)
           endif
        endif
        if(associated(p0%octNb6))then
           if( p0%octNb6%octType==0) then
              nullify(p0%octNb6%octNb5)
              nullify(p0%octNb6)
           endif
        endif
     enddo
  enddo

  return
end subroutine detach_G_n_M


!------------------------------------------------
!new routine for initialization with Morton
!------------------------------------------------
subroutine reset_all_mesh
  use init_mesh_size
  use param
  use oct_set
  use particle_set
  use message_passing_interface
  implicit none
  integer(kind=4)::iLv,index
  type(oct)  ,pointer::p0
  type(prtcl),pointer::pp

  do index=1,MaxID(4,LvMax)
     p0=>Mesh(index)
     p0%octP=0
     p0%octLv=0
     p0%iFLG=0
     p0%Csort=0
     p0%F=0.d0
     p0%C=0.d0
     p0%Z=0.d0
     p0%D=0.d0
     p0%O=0.d0
     p0%iC=0
     p0%ptcl_loops=0
     p0%prc_bndry=0
     p0%MrtN=0
     p0%octType=-1
     nullify(p0%octPrt)
     nullify(p0%octNb1)
     nullify(p0%octNb2)
     nullify(p0%octNb3)
     nullify(p0%octNb4)
     nullify(p0%octNb5)
     nullify(p0%octNb6)
     nullify(p0%octCh1)
     nullify(p0%octCh2)
     nullify(p0%octCh3)
     nullify(p0%octCh4)
     nullify(p0%octCh5)
     nullify(p0%octCh6)
     nullify(p0%octCh7)
     nullify(p0%octCh8)
     nullify(p0%Psort )
     nullify(p0%ptcl  )
     nullify(p0%ptclA )
  enddo

  MinID=0
  MaxID=0

  do iLv=0,LvMax
     do index=1,MaxIP(iLv)
        pp=>Pesh(index,iLv)
        pp%isort=0
        pp%ioct=0
     enddo
     MaxIP(iLv)=0
  enddo

  deallocate(iBMesh_arr)
  deallocate(iGMesh_arr)
  deallocate(GMn2octN)

end subroutine reset_all_mesh

subroutine removeBoundaryB(iLv)
    use param
    implicit none
    integer(kind=4)::iLv
    
    
    if(boundaryFlag(1)==1) then
      call removeBoundaryBx2(iLv)
    endif
    if(boundaryFlag(2)==1) then
      call removeBoundaryBy2(iLv)
    endif
    if(boundaryFlag(3)==1) then
      call removeBoundaryBz2(iLv)
    endif
end subroutine removeBoundaryB

subroutine removeBoundaryBCIP
    use param
    implicit none
    if(boundaryFlag(1)==1) then
      call removeBoundaryBxCIP
    endif
    if(boundaryFlag(2)==1) then
      call removeBoundaryByCIP
    endif
    if(boundaryFlag(3)==1) then
      call removeBoundaryBzCIP
    endif
end subroutine removeBoundaryBCIP

subroutine removeBoundaryBx(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::Nint2,index,iLv,ixx
  real(kind=8)::wx,wxh
  type(oct),pointer::p0
  Nint2=2**LvMax

  !$omp parallel do private(index,p0,ixx,wx,wxh)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     ixx=int(p0%iPOS(1)+Nint2)/int(2*Nint2) !<-- for x direction
     if((ixx.ge.NXB*NXR-2).or.(ixx.le.3)) then !<-- notice!
        p0%F(4)=0.0d0
        p0%F(5)=0.0d0
        p0%F(6)=0.0d0
     else 
        if(ixx<=wall+4) then
           wx=dble((wall+4-ixx))/dble(wall)*0.5d0
           wxh=(dble((wall+4-ixx))+0.5d0)/dble(wall)*0.5d0
           ! -- for x direction --
           p0%F(4)=p0%F(4)*(dcos(wx *PI))
           p0%F(5)=p0%F(5)*(dcos(wx *PI))
           p0%F(6)=p0%F(6)*(dcos(wx *PI))
        endif
        if(ixx>=NXB*NXR-3-wall) then !<-- notice
           wx=dble(wall-(NXB*NXR-3-ixx))/dble(wall)*0.5d0 !<-- notice
           wxh=(dble(wall-(NXB*NXR-3-ixx))-0.5d0)/dble(wall)*0.5d0 !<-- notice
           ! -- for x direction --
           p0%F(4)=p0%F(4)*(dcos(wx *PI))
           p0%F(5)=p0%F(5)*(dcos(wx *PI))
           p0%F(6)=p0%F(6)*(dcos(wx *PI))
        endif
     endif
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryBx

subroutine removeBoundaryBy(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::Nint2,index,iLv,iyy
  real(kind=8)::wy,wyh
  type(oct),pointer::p0
  Nint2=2**LvMax


  !$omp parallel do private(index,p0,iyy,wy,wyh)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     iyy=int(p0%iPOS(2)+Nint2)/int(2*Nint2) !<-- for y direction
     if((iyy.ge.NYB*NYR-2))then!.or.(iyy.le.3)) then !<-- notice!
        p0%F(4)=0.0d0
        p0%F(5)=0.0d0
        p0%F(6)=0.0d0
     else 
        if(iyy<=wall+3) then
           wy=dble((wall+4-iyy))/dble(wall)*0.5d0
           wyh=(dble((wall+4-iyy))+0.5d0)/dble(wall)*0.5d0
           ! -- for x direction --
           p0%F(4)=p0%F(4)*(dcos(wy *PI))
           p0%F(5)=p0%F(5)*(dcos(wy *PI))
           p0%F(6)=p0%F(6)*(dcos(wy *PI))
        endif
        if(iyy>=NYB*NYR-2-wall) then !<-- notice
           wy=dble(wall-(NYB*NYR-3-iyy))/dble(wall)*0.5d0 !<-- notice
           wyh=(dble(wall-(NYB*NYR-3-iyy))-0.5d0)/dble(wall)*0.5d0 !<-- notice
           ! -- for x direction --
           p0%F(4)=p0%F(4)*(dcos(wy *PI))
           p0%F(5)=p0%F(5)*(dcos(wy *PI))
           p0%F(6)=p0%F(6)*(dcos(wy *PI))
        endif
     endif
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryBy

subroutine removeBoundaryBz(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::Nint2,index,iLv,izz
  real(kind=8)::wz,wzh
  type(oct),pointer::p0
  Nint2=2**LvMax


  !$omp parallel do private(index,p0,izz,wz,wzh)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     izz=int(p0%iPOS(3)+Nint2)/int(2*Nint2) !<-- for z direction
     if((izz.ge.NZB*NZR-2))then!.or.(izz.le.3)) then !<-- notice!
        p0%F(4)=0.0d0
        p0%F(5)=0.0d0
        p0%F(6)=0.0d0
     else 
        if(izz<=wall+3) then
           wz=dble((wall+4-izz))/dble(wall)*0.5d0
           wzh=(dble((wall+4-izz))+0.5d0)/dble(wall)*0.5d0
           ! -- for x direction --
           p0%F(4)=p0%F(4)*(dcos(wz *PI))
           p0%F(5)=p0%F(5)*(dcos(wz *PI))
           p0%F(6)=p0%F(6)*(dcos(wz *PI))
        endif
        if(izz>=NZB*NZR-2-wall) then !<-- notice
           wz=dble(wall-(NZB*NZR-3-izz))/dble(wall)*0.5d0 !<-- notice
           wzh=(dble(wall-(NZB*NZR-3-izz))-0.5d0)/dble(wall)*0.5d0 !<-- notice
           ! -- for x direction --
           p0%F(4)=p0%F(4)*(dcos(wz *PI))
           p0%F(5)=p0%F(5)*(dcos(wz *PI))
           p0%F(6)=p0%F(6)*(dcos(wz *PI))
        endif
     endif
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryBz

subroutine removeBoundaryBx2(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::index,iLv
  real(kind=8)::wx,wxh
  type(oct),pointer::p0
  
  wx=dble(wall)*dx(1)
  wxh=dble(NXB*NXR-wall)*dx(1)
  !$omp parallel do private(index,p0)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if(p0%rPos(1)<=wx)then
        p0%F(4)=p0%F(4)*dsin(p0%rPos(1)*0.5d0*PI/wx)
        p0%F(5)=p0%F(5)*dsin(p0%rPos(1)*0.5d0*PI/wx)
        p0%F(6)=p0%F(6)*dsin(p0%rPos(1)*0.5d0*PI/wx)
     else if(p0%rPos(1)>=wxh)then
        p0%F(4)=p0%F(4)*dcos((p0%rPos(1)-wxh)*0.5d0*PI/wx)
        p0%F(5)=p0%F(5)*dcos((p0%rPos(1)-wxh)*0.5d0*PI/wx)
        p0%F(6)=p0%F(6)*dcos((p0%rPos(1)-wxh)*0.5d0*PI/wx)
     end if
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryBx2

subroutine removeBoundaryBy2(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::index,iLv
  real(kind=8)::wy,wyh
  type(oct),pointer::p0

  wy=dble(wall)*dx(2)
  wyh=dble(NYB*NYR-wall)*dx(2)
  !$omp parallel do private(index,p0)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if(p0%rPos(2)<=wy)then
        p0%F(4)=p0%F(4)*dsin(p0%rPos(2)*0.5d0*PI/wy)
        p0%F(5)=p0%F(5)*dsin(p0%rPos(2)*0.5d0*PI/wy)
        p0%F(6)=p0%F(6)*dsin(p0%rPos(2)*0.5d0*PI/wy)
     else if(p0%rPos(2)>=wyh)then
        p0%F(4)=p0%F(4)*dcos((p0%rPos(2)-wyh)*0.5d0*PI/wy)
        p0%F(5)=p0%F(5)*dcos((p0%rPos(2)-wyh)*0.5d0*PI/wy)
        p0%F(6)=p0%F(6)*dcos((p0%rPos(2)-wyh)*0.5d0*PI/wy)
     end if
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryBy2

subroutine removeBoundaryBz2(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::index,iLv
  real(kind=8)::wz,wzh
  type(oct),pointer::p0

  wz=dble(wall)*dx(3)
  wzh=dble(NZB*NZR-wall)*dx(3)
  !$omp parallel do private(index,p0)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if(p0%rPos(3)<=wz)then
        p0%F(4)=p0%F(4)*dsin(p0%rPos(3)*0.5d0*PI/wz)
        p0%F(5)=p0%F(5)*dsin(p0%rPos(3)*0.5d0*PI/wz)
        p0%F(6)=p0%F(6)*dsin(p0%rPos(3)*0.5d0*PI/wz)
     else if(p0%rPos(3)>=wzh)then
        p0%F(4)=p0%F(4)*dcos((p0%rPos(3)-wzh)*0.5d0*PI/wz)
        p0%F(5)=p0%F(5)*dcos((p0%rPos(3)-wzh)*0.5d0*PI/wz)
        p0%F(6)=p0%F(6)*dcos((p0%rPos(3)-wzh)*0.5d0*PI/wz)
     end if
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryBz2

subroutine removeBoundaryBxCIP
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::index,iLv
  type(oct),pointer::p0
  do iLv=0,LvMax
    if (MinID(1,iLv).ge.MaxID(1,iLv)) cycle
    do index=MinID(1,iLv),MaxID(1,iLv)
      p0 => Mesh(index)
      if((p0%rPOS(1).ge.R_lim(1,1)-dx(1)) .or. (p0%rPOS(1).le.dx(1)))then
        p0%F(4)=0.0d0
        p0%F(5)=0.0d0
        p0%F(6)=0.0d0
        p0%G(10:18)=0.0d0
      endif
    enddo
  enddo
  return
end subroutine removeBoundaryBxCIP

subroutine removeBoundaryByCIP
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::index,iLv
  type(oct),pointer::p0
  do iLv=0,LvMax
    if (MinID(1,iLv).ge.MaxID(1,iLv)) cycle
    do index=MinID(1,iLv),MaxID(1,iLv)
      p0 => Mesh(index)
      if((p0%rPOS(2).ge.R_lim(1,2)-dx(2)) .or. (p0%rPOS(2).le.dx(2)))then
        p0%F(4)=0.0d0
        p0%F(5)=0.0d0
        p0%F(6)=0.0d0
        p0%G(10:18)=0.0d0
      endif
    enddo
  enddo
  return
end subroutine removeBoundaryByCIP

subroutine removeBoundaryBzCIP
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::index,iLv
  type(oct),pointer::p0
  do iLv=0,LvMax
    if (MinID(1,iLv).ge.MaxID(1,iLv)) cycle
    do index=MinID(1,iLv),MaxID(1,iLv)
      p0 => Mesh(index)
      if((p0%rPOS(3).ge.R_lim(1,3)-dx(3)) .or. (p0%rPOS(3).le.dx(3)))then
        p0%F(4)=0.0d0
        p0%F(5)=0.0d0
        p0%F(6)=0.0d0
        p0%G(10:18)=0.0d0
      endif
    enddo
  enddo
  return
end subroutine removeBoundaryBzCIP

subroutine removeBoundaryE(iLv)
    use param
    implicit none
    integer(kind=4)::iLv
    
    if(boundaryFlag(1)==1) then
      call removeBoundaryEx2(iLv)
    endif
    if(boundaryFlag(2)==1) then
      call removeBoundaryEy2(iLv)
    endif
    if(boundaryFlag(3)==1) then
      call removeBoundaryEz2(iLv)
    endif
end subroutine removeBoundaryE

subroutine removeBoundaryECIP
    use param
    implicit none
    if(boundaryFlag(1)==1) then
      call removeBoundaryExCIP
    endif
    if(boundaryFlag(2)==1) then
      call removeBoundaryEyCIP
    endif
    if(boundaryFlag(3)==1) then
      call removeBoundaryEzCIP
    endif
end subroutine removeBoundaryECIP

subroutine removeBoundaryEx(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::Nint2,index,iLv,ixx
  real(kind=8)::wx,wxh
  type(oct),pointer::p0
  Nint2=2**LvMax

  !$omp parallel do private(index,p0,ixx,wx,wxh)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     ixx=int(p0%iPOS(1)+Nint2)/int(2*Nint2) !<-- for x direction
     if((ixx.ge.NXB*NXR-2).or.(ixx.le.3)) then !<-- notice!
        p0%F(1)=0.0d0
        p0%F(2)=0.0d0
        p0%F(3)=0.0d0
     else 
        if(ixx<=wall+4) then
           wx=dble((wall+4-ixx))/dble(wall)*0.5d0
           wxh=(dble((wall+4-ixx))+0.5d0)/dble(wall)*0.5d0
           ! -- for x direction --
           p0%F(1)=p0%F(1)*(dcos(wx *PI))
           p0%F(2)=p0%F(2)*(dcos(wx *PI))
           p0%F(3)=p0%F(3)*(dcos(wx *PI))
        endif
        if(ixx>=NXB*NXR-3-wall) then !<-- notice
           wx=dble(wall-(NXB*NXR-3-ixx))/dble(wall)*0.5d0 !<-- notice
           wxh=(dble(wall-(NXB*NXR-3-ixx))-0.5d0)/dble(wall)*0.5d0 !<-- notice
           ! -- for x direction --
           p0%F(1)=p0%F(1)*(dcos(wx *PI))
           p0%F(2)=p0%F(2)*(dcos(wx *PI))
           p0%F(3)=p0%F(3)*(dcos(wx *PI))
        endif
     endif
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryEx

subroutine removeBoundaryEy(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::Nint2,index,iLv,iyy
  real(kind=8)::wy,wyh
  type(oct),pointer::p0
  Nint2=2**LvMax

  !$omp parallel do private(index,p0,iyy,wy,wyh)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     iyy=int(p0%iPOS(2)+Nint2)/int(2*Nint2) !<-- for y direction
     if((iyy.ge.NYB*NYR-2))then!.or.(iyy.le.3)) then !<-- notice!
        p0%F(1)=0.0d0
        p0%F(2)=0.0d0
        p0%F(3)=0.0d0
     else 
        if(iyy<=wall+3) then
           wy=dble((wall+4-iyy))/dble(wall)*0.5d0
           wyh=(dble((wall+4-iyy))+0.5d0)/dble(wall)*0.5d0
           ! -- for x direction --
           p0%F(1)=p0%F(1)*(dcos(wy *PI))
           p0%F(2)=p0%F(2)*(dcos(wy *PI))
           p0%F(3)=p0%F(3)*(dcos(wy *PI))
        endif
        if(iyy>=NYB*NYR-2-wall) then !<-- notice
           wy=dble(wall-(NYB*NYR-3-iyy))/dble(wall)*0.5d0 !<-- notice
           wyh=(dble(wall-(NYB*NYR-3-iyy))-0.5d0)/dble(wall)*0.5d0 !<-- notice
           ! -- for x direction --
           p0%F(1)=p0%F(1)*(dcos(wy *PI))
           p0%F(2)=p0%F(2)*(dcos(wy *PI))
           p0%F(3)=p0%F(3)*(dcos(wy *PI))
        endif
     endif
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryEy

subroutine removeBoundaryEz(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::Nint2,index,iLv,izz
  real(kind=8)::wz,wzh
  type(oct),pointer::p0
  Nint2=2**LvMax

  !$omp parallel do private(index,p0,izz,wz,wzh)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     izz=int(p0%iPOS(3)+Nint2)/int(2*Nint2) !<-- for z direction
     if((izz.ge.NZB*NZR-2))then!.or.(izz.le.3)) then !<-- notice!
        p0%F(1)=0.0d0
        p0%F(2)=0.0d0
        p0%F(3)=0.0d0
     else 
        if(izz<=wall+3) then
           wz=dble((wall+4-izz))/dble(wall)*0.5d0
           wzh=(dble((wall+4-izz))+0.5d0)/dble(wall)*0.5d0
           ! -- for x direction --
           p0%F(1)=p0%F(1)*(dcos(wz *PI))
           p0%F(2)=p0%F(2)*(dcos(wz *PI))
           p0%F(3)=p0%F(3)*(dcos(wz *PI))
        endif
        if(izz>=NZB*NZR-2-wall) then !<-- notice
           wz=dble(wall-(NZB*NZR-3-izz))/dble(wall)*0.5d0 !<-- notice
           wzh=(dble(wall-(NZB*NZR-3-izz))-0.5d0)/dble(wall)*0.5d0 !<-- notice
           ! -- for x direction --
           p0%F(1)=p0%F(1)*(dcos(wz *PI))
           p0%F(2)=p0%F(2)*(dcos(wz *PI))
           p0%F(3)=p0%F(3)*(dcos(wz *PI))
        endif
     endif
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryEz

subroutine removeBoundaryEx2(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::index,iLv
  real(kind=8)::wx,wxh
  type(oct),pointer::p0
  
  wx=dble(wall)*dx(1)
  wxh=dble(NXB*NXR-wall)*dx(1)
  !$omp parallel do private(index,p0)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if(p0%rPos(1)<=wx)then
        p0%F(1)=p0%F(1)*dsin(p0%rPos(1)*0.5d0*PI/wx)
        p0%F(2)=p0%F(2)*dsin(p0%rPos(1)*0.5d0*PI/wx)
        p0%F(3)=p0%F(3)*dsin(p0%rPos(1)*0.5d0*PI/wx)
     else if(p0%rPos(1)>=wxh)then
        p0%F(1)=p0%F(1)*dcos((p0%rPos(1)-wxh)*0.5d0*PI/wx)
        p0%F(2)=p0%F(2)*dcos((p0%rPos(1)-wxh)*0.5d0*PI/wx)
        p0%F(3)=p0%F(3)*dcos((p0%rPos(1)-wxh)*0.5d0*PI/wx)
     end if
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryEx2

subroutine removeBoundaryEy2(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::index,iLv
  real(kind=8)::wy,wyh
  type(oct),pointer::p0

  wy=dble(wall)*dx(2)
  wyh=dble(NYB*NYR-wall)*dx(2)
  !$omp parallel do private(index,p0)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if(p0%rPos(2)<=wy)then
        p0%F(1)=p0%F(1)*dsin(p0%rPos(2)*0.5d0*PI/wy)
        p0%F(2)=p0%F(2)*dsin(p0%rPos(2)*0.5d0*PI/wy)
        p0%F(3)=p0%F(3)*dsin(p0%rPos(2)*0.5d0*PI/wy)
     else if(p0%rPos(2)>=wyh)then
        p0%F(1)=p0%F(1)*dcos((p0%rPos(2)-wyh)*0.5d0*PI/wy)
        p0%F(2)=p0%F(2)*dcos((p0%rPos(2)-wyh)*0.5d0*PI/wy)
        p0%F(3)=p0%F(3)*dcos((p0%rPos(2)-wyh)*0.5d0*PI/wy)
     end if
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryEy2

subroutine removeBoundaryEz2(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none
  integer(kind=4)::index,iLv
  real(kind=8)::wz,wzh
  type(oct),pointer::p0
  
  wz=dble(wall)*dx(3)
  wzh=dble(NZB*NZR-wall)*dx(3)
  !$omp parallel do private(index,p0)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if(p0%rPos(3)<=wz)then
        p0%F(1)=p0%F(1)*dsin(p0%rPos(3)*0.5d0*PI/wz)
        p0%F(2)=p0%F(2)*dsin(p0%rPos(3)*0.5d0*PI/wz)
        p0%F(3)=p0%F(3)*dsin(p0%rPos(3)*0.5d0*PI/wz)
     else if(p0%rPos(3)>=wzh)then
        p0%F(1)=p0%F(1)*dcos((p0%rPos(3)-wzh)*0.5d0*PI/wz)
        p0%F(2)=p0%F(2)*dcos((p0%rPos(3)-wzh)*0.5d0*PI/wz)
        p0%F(3)=p0%F(3)*dcos((p0%rPos(3)-wzh)*0.5d0*PI/wz)
     end if
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryEz2

subroutine removeBoundaryExCIP
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::index,iLv
  type(oct),pointer::p0
  do iLv=0,LvMax
    if (MinID(1,iLv).ge.MaxID(1,iLv)) cycle
    do index=MinID(1,iLv),MaxID(1,iLv)
      p0 => Mesh(index)
      if((p0%rPOS(1).ge.R_lim(1,1)-dx(1)) .or. (p0%rPOS(1).le.dx(1)))then
        p0%F(1)=0.0d0
        p0%F(2)=0.0d0
        p0%F(3)=0.0d0
        p0%G(1:9)=0.0d0
      endif
    enddo
  enddo
  return
end subroutine removeBoundaryExCIP

subroutine removeBoundaryEyCIP
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::index,iLv
  type(oct),pointer::p0
  do iLv=0,LvMax
    if (MinID(1,iLv).ge.MaxID(1,iLv)) cycle
    do index=MinID(1,iLv),MaxID(1,iLv)
      p0 => Mesh(index)
      if((p0%rPOS(2).ge.R_lim(1,2)-dx(2)) .or. (p0%rPOS(2).le.dx(2)))then
        p0%F(1)=0.0d0
        p0%F(2)=0.0d0
        p0%F(3)=0.0d0
        p0%G(1:9)=0.0d0
      endif
    enddo
  enddo
  return
end subroutine removeBoundaryEyCIP

subroutine removeBoundaryEzCIP
  use param
  use oct_set
  use const
  use init_mesh_size
  use message_passing_interface
  implicit none

  integer(kind=4)::index,iLv
  type(oct),pointer::p0
  do iLv=0,LvMax
    if (MinID(1,iLv).ge.MaxID(1,iLv)) cycle
    do index=MinID(1,iLv),MaxID(1,iLv)
      p0 => Mesh(index)
      if((p0%rPOS(3).ge.R_lim(1,3)-dx(3)) .or. (p0%rPOS(3).le.dx(3)))then
        p0%F(1)=0.0d0
        p0%F(2)=0.0d0
        p0%F(3)=0.0d0
        p0%G(1:9)=0.0d0
      endif
    enddo
  enddo
  return
end subroutine removeBoundaryEzCIP

subroutine removeBoundaryJ(iLv)
  use param
  use oct_set
  use const
  use init_mesh_size
  implicit none

  integer(kind=4)::Nint2,index,iLv,ixx
  real(kind=8)::wx,wxh
  type(oct),pointer::p0
  Nint2=2**LvMax
  !$omp parallel do private(index,p0,ixx,wx,wxh)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     ixx=int(p0%iPOS(1)+Nint2)/int(2*Nint2) !<-- for x direction
     if((ixx.le.3).or.(ixx.ge.NXB*NXR-2)) then !<-- notice!
        p0%F(7)=0.0d0
        p0%F(8)=0.0d0
        p0%F(9)=0.0d0
        p0%F(10)=0.0d0
        p0%F(11)=0.0d0
        p0%F(12)=0.0d0
     else 
        if(ixx<=wall+3) then
           wx=dble((wall+4-ixx))/dble(wall)*0.5d0
           wxh=(dble((wall+4-ixx))+0.5d0)/dble(wall)*0.5d0
           ! -- for x direction --
           p0%F(7)=p0%F(7)*(dcos(wx *PI))
           p0%F(8)=p0%F(8)*(dcos(wx *PI))
           p0%F(9)=p0%F(9)*(dcos(wx *PI))
           p0%F(10)=p0%F(10)*(dcos(wx *PI))
           p0%F(11)=p0%F(11)*(dcos(wx *PI))
           p0%F(12)=p0%F(12)*(dcos(wx *PI))
        endif
        if(ixx>=NXB*NXR-2-wall) then !<-- notice
           wx=dble(wall-(NXB*NXR-3-ixx))/dble(wall)*0.5d0 !<-- notice
           wxh=(dble(wall-(NXB*NXR-3-ixx))-0.5d0)/dble(wall)*0.5d0 !<-- notice
           ! -- for x direction --
           p0%F(7)=p0%F(7)*(dcos(wx *PI))
           p0%F(8)=p0%F(8)*(dcos(wx *PI))
           p0%F(9)=p0%F(9)*(dcos(wx *PI))
           p0%F(10)=p0%F(10)*(dcos(wx *PI))
           p0%F(11)=p0%F(11)*(dcos(wx *PI))
           p0%F(12)=p0%F(12)*(dcos(wx *PI))
        endif
     endif
  enddo
  !$omp end parallel do
  return
end subroutine removeBoundaryJ

subroutine makeInjectionOctTable(size)
    use param
    use oct_set
    use init_mesh_size
    use message_passing_interface
    implicit none
    
    integer(kind=4)::i,index,size
    type(oct),pointer::p0
    i=1
    do index=MinID(1,-1),MaxID(1,-1)
        p0=>Mesh(index)
        if(p0%iPos(1)==2**(LvMax+1))then
            InjectionTable(i)=p0%octCh1%octN
            i=i+1
            InjectionTable(i)=p0%octCh3%octN
            i=i+1
            InjectionTable(i)=p0%octCh5%octN
            i=i+1
            InjectionTable(i)=p0%octCh7%octN
            i=i+1
        end if
    end do
!    print *,'tableNum =',i, rank
    size=i-1
end subroutine makeInjectionOctTable

subroutine particle_injctT_xs
  use param
  use oct_set
  use particle_set
  use init_mesh_size
  use const
  use init_condition
  use message_passing_interface
  use injectionOct
  implicit none
  integer(kind=4)::ii,index
  real(kind=8)::X0,Y0,Z0,X1,Y1,Z1,angleT,angleP,ratio
  real(kind=8),dimension(Niall)::PEinitX,PEinitY,PEinitZ
  real(kind=8),dimension(Niall)::PQinitX,PQinitY,PQinitZ
  real(kind=8),dimension(Niall)::PPinitX,PPinitY,PPinitZ
  integer(kind=4)::Ninjct,IPinit
  type(injectionStructure)::is
  ! -- inject number of particles --
  Ninjct=npart_per_cell*NYB*NZB
  angleT=0.5d0*PI
  angleP=0.0d0
  ratio =1.0d0
  !
  X0=0.0d0; X1=1.0d0
  Y0=0.0d0; Y1=1.0d0
  Z0=0.0d0; Z1=1.0d0
  call prand(PQinitX,Niall,Ninjct,X0,X1)
  call prand(PQinitY,Niall,Ninjct,Y0,Y1)
  call prand(PQinitZ,Niall,Ninjct,Z0,Z1)
  call set_thermal_velocity(TVELE,PEinitX,Niall,Ninjct)
  call set_thermal_velocity(TVELE,PEinitY,Niall,Ninjct)
  call set_thermal_velocity(TVELE,PEinitZ,Niall,Ninjct)
  call set_thermal_velocity(TVION,PPinitX,Niall,Ninjct)
  call set_thermal_velocity(TVION,PPinitY,Niall,Ninjct)
  call set_thermal_velocity(TVION,PPinitZ,Niall,Ninjct)
  ii=0
  
  ! Background Particle
  if(BackgroundParticle==1)then
    do index=1,tableSize
        IPinit=InjectionTable(index)
        ! electron
        call injectParticleEveryStep(IPinit,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,1,angleT,angleP,ratio)
    end do
    ii=0
    do index=1,tableSize
        IPinit=InjectionTable(index)
        ! ion
        call injectParticleEveryStep(IPinit,ii,PQinitX,PQinitY,PQinitZ,PPinitX,PPinitY,PPinitZ,2,angleT,angleP,ratio)
    end do
  end if
  
  !Beam particle injection
  if(BeamParticle==1)then
    do index=1,arraySize
        is = InjectionTableArray(index)
        if (is%iSort.eq.1) then
            call injectParticleEveryStep(is%octNum,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,is%iSort,is%theta,is%phi,is%ratio)
        else
            call injectParticleEveryStep(is%octNum,ii,PQinitX,PQinitY,PQinitZ,PPinitX,PPinitY,PPinitZ,is%iSort,is%theta,is%phi,is%ratio)
        endif
    enddo
  end if
  return 
end subroutine particle_injctT_xs

recursive subroutine injectParticleEveryStep(index,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,particle,angleT,angleP,ratio)
  use param
  use oct_set
  use particle_set
  use init_mesh_size
  use const
  use init_condition
  use message_passing_interface
  use Vparam
  implicit none
  integer(kind=4)::index, iLv, num, particle, j
  integer(kind=4),intent(inout)::ii
  real(kind=8)::angleT, angleP, ratio,sLv
  real(kind=8),dimension(Niall)::PEinitX,PEinitY,PEinitZ
  real(kind=8),dimension(Niall)::PQinitX,PQinitY,PQinitZ
  real(kind=8),dimension(1:9)::R
  real(kind=8),dimension(1:3)::dxs
  ! -- pointers --
  type(oct),pointer::p0
  type(prtcl),pointer::PrtList
  !===========Nakano==============
  if(beamParticle == 1)then 
	p0 => VMesh(index)
  else
    p0 => Mesh(index)
  end if
  !==============================
  if (p0%iFLG(1) .ge. 4) then
    call injectParticleEveryStep(p0%octCh1%octN,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,particle,angleT,angleP,ratio)
    call injectParticleEveryStep(p0%octCh2%octN,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,particle,angleT,angleP,ratio)
    call injectParticleEveryStep(p0%octCh3%octN,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,particle,angleT,angleP,ratio)
    call injectParticleEveryStep(p0%octCh4%octN,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,particle,angleT,angleP,ratio)
    call injectParticleEveryStep(p0%octCh5%octN,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,particle,angleT,angleP,ratio)
    call injectParticleEveryStep(p0%octCh6%octN,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,particle,angleT,angleP,ratio)
    call injectParticleEveryStep(p0%octCh7%octN,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,particle,angleT,angleP,ratio)
    call injectParticleEveryStep(p0%octCh8%octN,ii,PQinitX,PQinitY,PQinitZ,PEinitX,PEinitY,PEinitZ,particle,angleT,angleP,ratio)
  else if (p0%iFLG(1) .ge. 0 .and. p0%iFLG(1) .le. 3) then
    iLv = p0%octLv
    sLv=0.5d0**iLv
    dxs(1)=dx(1)*sLv
    dxs(2)=dx(2)*sLv
    dxs(3)=dx(3)*sLv
    num = 0
    j=MaxIP(iLv)+1
    do while (dble(num) .lt. dble(npart_per_cell)*ratio*(sLv**3.0d0))
      ii=ii+1         ! <-- storage index
      R(1)=PQinitX(ii)*dxs(1)+(p0%rPOS(1)-0.5d0*dxs(1))
      R(2)=PQinitY(ii)*dxs(2)+(p0%rPOS(2)-0.5d0*dxs(2))
      R(3)=PQinitZ(ii)*dxs(3)+(p0%rPOS(3)-0.5d0*dxs(3))
      R(4)=PEinitX(ii)+injectionFlow*dsin(angleT)*dcos(angleP)
      R(5)=PEinitY(ii)+injectionFlow*dsin(angleT)*dsin(angleP)
      R(6)=PEinitZ(ii)+injectionFlow*dcos(angleT)
      R(7)=R(4); R(8)=R(5); R(9)=R(6)
     ! -- strage to the list structure for particles --
      PrtList => p0%ptcl%prtnxt
      Pesh(j,iLv)=prtcl(R,particle,p0%octN,j,PrtList,beamParticle)
      p0%ptcl%prtnxt => Pesh(j,iLv)
      p0%octP=p0%octP+1
	  
	 ! print *,"octP",p0%octP
	  
      num = num + 1
      j=j+1           ! <-- Pesh index
    end do
    MaxIP(iLv)=MaxIP(iLv)+num
  end if
  return
end subroutine injectParticleEveryStep

subroutine removeParticleForInjection
  use param
  use injectionOct
  implicit none
  integer(kind=4)::iyz
  integer(kind=4)::index
  type(injectionStructure)::is
  
  !Background Particle remove
  if(backgroundParticle==1)then
    do iyz=1,TableSize
        index=InjectionTable(iyz)
        call removeParticleForInjectionEveryOct(index)
    enddo
  end if
  !beam Particle remove
  if(beamParticle==1)then
    do iyz=1,arraySize
        is = InjectionTableArray(iyz)
        call removeParticleForInjectionEveryOct(is%octNum)
    enddo
  end if
  if(satellite==1)then
    do iyz=1,satelliteArraySize
        index = satelliteTableArray(iyz)
        call removeParticleForInjectionEveryOct(index)
    enddo
  end if
  return
end subroutine removeParticleForInjection

recursive subroutine removeParticleForInjectionEveryOct(IPinit)
  use param
  use oct_set
  use particle_set
  use init_mesh_size
  use Vparam
  implicit none
  integer(kind=4)::i
  integer(kind=4)::indexP
  integer(kind=4)::IPinit
  ! -- pointers --
  type(oct),pointer::p0,p1
  type(prtcl),pointer::pptc
!==============Nakano=============
  if(beamParticle == 1)then
  p0 => VMesh(IPinit)
  else 
  p0 => Mesh(IPinit)
  end if
  !===============================
  if (p0%iFLG(1) .ge. 4) then
    call removeParticleForInjectionEveryOct(p0%octCh1%octN)
    call removeParticleForInjectionEveryOct(p0%octCh2%octN)
    call removeParticleForInjectionEveryOct(p0%octCh3%octN)
    call removeParticleForInjectionEveryOct(p0%octCh4%octN)
    call removeParticleForInjectionEveryOct(p0%octCh5%octN)
    call removeParticleForInjectionEveryOct(p0%octCh6%octN)
    call removeParticleForInjectionEveryOct(p0%octCh7%octN)
    call removeParticleForInjectionEveryOct(p0%octCh8%octN)
  else if (p0%octP .ne. 0 .and. p0%iFLG(1) .ge. 0 .and. p0%iFLG(1) .le. 3) then
    indexP = p0%octP
    pptc => p0%ptcl%prtnxt
    do i=1,indexP
       pptc%isort = 0
       pptc => pptc%prtNxt
    end do
    p0%octP = 0
  endif
 
  return
end subroutine removeParticleForInjectionEveryOct
  
subroutine prand(PEinit,Niall,Ninjct,X0,X1)
  use mtmod
  implicit none
  ! global variabeles
  integer(kind=4)::Niall,Ninjct
  real(kind=8)::PEinit(1:Niall),X0,X1
  ! local variables
  integer(kind=4)::i
  do i=1,Ninjct
     PEinit(i)=(X1-X0)*grnd()+X0
  enddo
  return
end subroutine prand

subroutine set_thermal_velocity(Vthermal,PEinitX,Niall,Ninjct)
  use mtmod
  implicit none
  ! -- global variables --
  integer(kind=4)::Niall,Ninjct
  real(kind=8)::Vthermal,PEinitX(1:Niall)
  ! -- local variables --
  integer(kind=4)::i
  real(kind=8)::al,pi2,b1,b2
  pi2=datan(1.0d0)*8.0d0
  do i=1,Ninjct
     b1=grnd()
     b2=grnd()
     do while(b1<1.0d-100)
        b1=grnd()
     enddo
     al=Vthermal*dsqrt(-2.0d0*dlog(b1))
     PEinitX(i)=al*dcos(b2*pi2)
  enddo
  return
end subroutine set_thermal_velocity

subroutine setAxis(random, output, ptclNum, start, end, endPos, startPos)
  use param
  implicit none
  real(kind=8)::startPos,endPos,width
  integer(kind=4)::ptclNum,start,end,i
  real(kind=8),dimension(ptclNum*NX*NY*NZ*NXR*NYR*NZR)::random
  real(kind=8),dimension(ptclNum)::output
  
  width = endPos - startPos
  do i=start,end
    output(i-start+1) = width*random(i)+startPos
  enddo
end subroutine setAxis

subroutine setRandom(seed,output,num,start,end)
    implicit none
    integer(kind=4)::seed,num,i,s
    integer,allocatable::new(:)
    real(kind=8)::start,end,temp
    real(kind=8),dimension(num)::output
    call random_seed(size=s)
    allocate (new(s))
    new(:) = seed

    temp = end - start
    call random_seed(put=new)
    call random_number(output)
    do i=1,num
        output(i) = output(i)*temp+start
    end do
end subroutine setRandom

subroutine non_advection(iLv,iFs,iFe,dstep)
  use param
  use oct_set
  use init_mesh_size
  implicit none
  integer(kind=4)::iLv,iFs,iFe
  integer(kind=4)::index
  real(kind=8)::dstep
  type(oct),pointer::p0
  if(MaxID(1,iLv)<=MinID(1,iLv)) return
  charge = 1 / omega / dble(npart_per_cell)
  ! -- non-advection term --
  !$omp parallel do private(index,p0)
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        p0%F(1)=p0%F(1)-charge*dstep*dt*(0.5d0**(iLv))*p0%F(7)
        p0%F(2)=p0%F(2)-charge*dstep*dt*(0.5d0**(iLv))*p0%F(8)
        p0%F(3)=p0%F(3)-charge*dstep*dt*(0.5d0**(iLv))*p0%F(9)
        p0%G(1)=p0%G(1)-charge*dstep*dt &
             *0.5d0*(p0%octNb2%F(7)-p0%octNb1%F(7))/dx(1)
        p0%G(2)=p0%G(2)-charge*dstep*dt &
             *0.5d0*(p0%octNb2%F(8)-p0%octNb1%F(8))/dx(1)
        p0%G(3)=p0%G(3)-charge*dstep*dt &
             *0.5d0*(p0%octNb2%F(9)-p0%octNb1%F(9))/dx(1)
        p0%G(4)=p0%G(4)-charge*dstep*dt &
             *0.5d0*(p0%octNb4%F(7)-p0%octNb3%F(7))/dx(2)
        p0%G(5)=p0%G(5)-charge*dstep*dt &
             *0.5d0*(p0%octNb4%F(8)-p0%octNb3%F(8))/dx(2)
        p0%G(6)=p0%G(6)-charge*dstep*dt &
             *0.5d0*(p0%octNb4%F(9)-p0%octNb3%F(9))/dx(2)
        p0%G(7)=p0%G(7)-charge*dstep*dt &
             *0.5d0*(p0%octNb6%F(7)-p0%octNb5%F(7))/dx(3)
        p0%G(8)=p0%G(8)-charge*dstep*dt &
             *0.5d0*(p0%octNb6%F(8)-p0%octNb5%F(8))/dx(3)
        p0%G(9)=p0%G(9)-charge*dstep*dt &
             *0.5d0*(p0%octNb6%F(9)-p0%octNb5%F(9))/dx(3)
     endif
  enddo
  !$omp end parallel do
  return
end subroutine non_advection

subroutine advection_x(iLv,iFs,iFe,dstep)
  use param
  use oct_set
  use init_mesh_size
  implicit none
  integer(kind=4)::iLv,iFs,iFe
  integer(kind=4)::index
  real(kind=8)::dstep
  real(kind=8)::us,PV,PVup,VIA,VIAup,nPIV,nVIA,dts,dxs,fi,fiup,fid,fn
  type(oct),pointer::p0
  if(MaxID(1,iLv)<=MinID(1,iLv)) return
  dts=dt*(0.5d0**iLv); dxs=dx(1)*(0.5d0**iLv)
  ! -- for x direction --
  ! -- 1st step --
  ! --- calc. numerical flux for Ey,Bz at real oct ---
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        ! Ey 1st
        us= 1.0d0; PV=p0%F(2); PVup=p0%octNb1%F(2)
        VIA=p0%G(2); VIAup=p0%octNb1%G(2)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 1)=nPIV; p0%C( 3)=nVIA
        us=-1.0d0; PV=p0%F(2); PVup=p0%octNb2%F(2)
        VIA=p0%G(2); VIAup=p0%octNb2%G(2)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 2)=nPIV; p0%C( 4)=nVIA
        ! Ey 2nd
        fi=p0%G(5); fiup=p0%octNb2%G(5); fid=p0%octNb1%G(5)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 5)=fn
        fi=p0%G(8); fiup=p0%octNb2%G(8); fid=p0%octNb1%G(8)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 6)=fn
        ! Bz 1st
        us= 1.0d0; PV=p0%F(6); PVup=p0%octNb1%F(6)
        VIA=p0%G(12); VIAup=p0%octNb1%G(12)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 9)=nPIV; p0%C(11)=nVIA
        us=-1.0d0; PV=p0%F(6); PVup=p0%octNb2%F(6)
        VIA=p0%G(12); VIAup=p0%octNb2%G(12)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C(10)=nPIV; p0%C(12)=nVIA
        ! Bz 2nd
        fi=p0%G(15); fiup=p0%octNb2%G(15); fid=p0%octNb1%G(15)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(13)=fn
        fi=p0%G(18); fiup=p0%octNb2%G(18); fid=p0%octNb1%G(18)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(14)=fn
     endif
  enddo
  ! -- reconstruct --
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        p0%F( 2)=0.5d0*((p0%C(1)+p0%C(2))+(p0%C( 9)-p0%C(10)))
        p0%G( 2)=0.5d0*((p0%C(3)+p0%C(4))+(p0%C(11)-p0%C(12)))
        p0%G( 5)=p0%C(5)
        p0%G( 8)=p0%C(6)
        p0%F( 6)=0.5d0*((p0%C(1)-p0%C(2))+(p0%C( 9)+p0%C(10)))
        p0%G(12)=0.5d0*((p0%C(3)-p0%C(4))+(p0%C(11)+p0%C(12)))
        p0%G(15)=p0%C(13)
        p0%G(18)=p0%C(14)
     endif
  enddo
  ! -- 2nd step --
  ! --- calc. numerical flux for Ez,By at real oct ---
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        ! Ez 1st
        us= 1.0d0; PV=p0%F(3); PVup=p0%octNb1%F(3)
        VIA=p0%G(3); VIAup=p0%octNb1%G(3)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 1)=nPIV; p0%C( 3)=nVIA
        us=-1.0d0; PV=p0%F(3); PVup=p0%octNb2%F(3)
        VIA=p0%G(3); VIAup=p0%octNb2%G(3)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 2)=nPIV; p0%C( 4)=nVIA
        ! Ez 2nd
        fi=p0%G(6); fiup=p0%octNb2%G(6); fid=p0%octNb1%G(6)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 5)=fn
        fi=p0%G(9); fiup=p0%octNb2%G(9); fid=p0%octNb1%G(9)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 6)=fn
        ! By 1st
        us= 1.0d0; PV=p0%F(5); PVup=p0%octNb1%F(5)
        VIA=p0%G(11); VIAup=p0%octNb1%G(11)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 9)=nPIV; p0%C(11)=nVIA
        us=-1.0d0; PV=p0%F(5); PVup=p0%octNb2%F(5)
        VIA=p0%G(11); VIAup=p0%octNb2%G(11)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C(10)=nPIV; p0%C(12)=nVIA
        ! By 2nd
        fi=p0%G(14); fiup=p0%octNb2%G(14); fid=p0%octNb1%G(14)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(13)=fn
        fi=p0%G(17); fiup=p0%octNb2%G(17); fid=p0%octNb1%G(17)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(14)=fn
     endif
  enddo
  ! -- reconstruct --
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        p0%F( 3)=0.5d0*(( p0%C(1)+p0%C(2))-(p0%C( 9)-p0%C(10)))
        p0%G( 3)=0.5d0*(( p0%C(3)+p0%C(4))-(p0%C(11)-p0%C(12)))
        p0%G( 6)=p0%C(5)
        p0%G( 9)=p0%C(6)
        p0%F( 5)=0.5d0*((-p0%C(1)+p0%C(2))+(p0%C( 9)+p0%C(10)))
        p0%G(11)=0.5d0*((-p0%C(3)+p0%C(4))+(p0%C(11)+p0%C(12)))
        p0%G(14)=p0%C(13)
        p0%G(17)=p0%C(14)
     endif
  enddo
  return
end subroutine advection_x

subroutine advection_y(iLv,iFs,iFe,dstep)
  use param
  use oct_set
  use init_mesh_size
  implicit none
  integer(kind=4)::iLv,iFs,iFe
  integer(kind=4)::index
  real(kind=8)::dstep
  real(kind=8)::us,PV,PVup,VIA,VIAup,nPIV,nVIA,dts,dxs,fi,fiup,fid,fn
  type(oct),pointer::p0
  if(MaxID(1,iLv)<=MinID(1,iLv)) return
  dts=dt*(0.5d0**iLv); dxs=dx(2)*(0.5d0**iLv)
  ! -- for y direction --
  ! -- 1st step --
  ! --- calc. numerical flux for Ez,Bx at real oct ---
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        ! Ez 1st
        us= 1.0d0; PV=p0%F(3); PVup=p0%octNb3%F(3)
        VIA=p0%G(6); VIAup=p0%octNb3%G(6)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 1)=nPIV; p0%C( 3)=nVIA
        us=-1.0d0; PV=p0%F(3); PVup=p0%octNb4%F(3)
        VIA=p0%G(6); VIAup=p0%octNb4%G(6)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 2)=nPIV; p0%C( 4)=nVIA
        ! Ez 2nd
        fi=p0%G(3); fiup=p0%octNb4%G(3); fid=p0%octNb3%G(3)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 5)=fn
        fi=p0%G(9); fiup=p0%octNb4%G(9); fid=p0%octNb3%G(9)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 6)=fn
        ! Bx 1st
        us= 1.0d0; PV=p0%F(4); PVup=p0%octNb3%F(4)
        VIA=p0%G(13); VIAup=p0%octNb3%G(13)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 9)=nPIV; p0%C(11)=nVIA
        us=-1.0d0; PV=p0%F(4); PVup=p0%octNb4%F(4)
        VIA=p0%G(13); VIAup=p0%octNb4%G(13)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C(10)=nPIV; p0%C(12)=nVIA
        ! Bx 2nd
        fi=p0%G(10); fiup=p0%octNb4%G(10); fid=p0%octNb3%G(10)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(13)=fn
        fi=p0%G(16); fiup=p0%octNb4%G(16); fid=p0%octNb3%G(16)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(14)=fn
     endif
  enddo
  ! -- reconstruct --
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        p0%F( 3)=0.5d0*((p0%C(1)+p0%C(2))+(p0%C( 9)-p0%C(10)))
        p0%G( 6)=0.5d0*((p0%C(3)+p0%C(4))+(p0%C(11)-p0%C(12)))
        p0%G( 3)=p0%C(5)
        p0%G( 9)=p0%C(6)
        p0%F( 4)=0.5d0*((p0%C(1)-p0%C(2))+(p0%C( 9)+p0%C(10)))
        p0%G(13)=0.5d0*((p0%C(3)-p0%C(4))+(p0%C(11)+p0%C(12)))
        p0%G(10)=p0%C(13)
        p0%G(16)=p0%C(14)
     endif
  enddo
  ! -- 2nd step --
  ! --- calc. numerical flux for Ex,Bz at real oct ---
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        ! Ex 1st
        us= 1.0d0; PV=p0%F(1); PVup=p0%octNb3%F(1)
        VIA=p0%G(4); VIAup=p0%octNb3%G(4)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 1)=nPIV; p0%C( 3)=nVIA
        us=-1.0d0; PV=p0%F(1); PVup=p0%octNb4%F(1)
        VIA=p0%G(4); VIAup=p0%octNb4%G(4)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 2)=nPIV; p0%C( 4)=nVIA
        ! Ex 2nd
        fi=p0%G(1); fiup=p0%octNb4%G(1); fid=p0%octNb3%G(1)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 5)=fn
        fi=p0%G(7); fiup=p0%octNb4%G(7); fid=p0%octNb3%G(7)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 6)=fn
        ! Bz 1st
        us= 1.0d0; PV=p0%F(6); PVup=p0%octNb3%F(6)
        VIA=p0%G(15); VIAup=p0%octNb3%G(15)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 9)=nPIV; p0%C(11)=nVIA
        us=-1.0d0; PV=p0%F(6); PVup=p0%octNb4%F(6)
        VIA=p0%G(15); VIAup=p0%octNb4%G(15)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C(10)=nPIV; p0%C(12)=nVIA
        ! Bz 2nd
        fi=p0%G(12); fiup=p0%octNb4%G(12); fid=p0%octNb3%G(12)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(13)=fn
        fi=p0%G(18); fiup=p0%octNb4%G(18); fid=p0%octNb3%G(18)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(14)=fn
     endif
  enddo
  ! -- reconstruct --
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        p0%F( 1)=0.5d0*((p0%C(1)+p0%C(2))-(p0%C( 9)-p0%C(10)))
        p0%G( 4)=0.5d0*((p0%C(3)+p0%C(4))-(p0%C(11)-p0%C(12)))
        p0%G( 1)=p0%C(5)
        p0%G( 7)=p0%C(6)
        p0%F( 6)=0.5d0*((-p0%C(1)+p0%C(2))+(p0%C( 9)+p0%C(10)))
        p0%G(15)=0.5d0*((-p0%C(3)+p0%C(4))+(p0%C(11)+p0%C(12)))
        p0%G(12)=p0%C(13)
        p0%G(18)=p0%C(14)
     endif
  enddo
  return
end subroutine advection_y

subroutine advection_z(iLv,iFs,iFe,dstep)
  use param
  use oct_set
  use init_mesh_size
  implicit none
  integer(kind=4)::iLv,iFs,iFe
  integer(kind=4)::index
  real(kind=8)::dstep
  real(kind=8)::us,PV,PVup,VIA,VIAup,nPIV,nVIA,dts,dxs,fi,fiup,fid,fn
  type(oct),pointer::p0
  if(MaxID(1,iLv)<=MinID(1,iLv)) return
  dts=dt*(0.5d0**iLv); dxs=dx(2)*(0.5d0**iLv)
  ! -- for y direction --
  ! -- 1st step --
  ! --- calc. numerical flux for Ez,Bx at real oct ---
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        ! Ex 1st
        us= 1.0d0; PV=p0%F(1); PVup=p0%octNb5%F(1)
        VIA=p0%G(7); VIAup=p0%octNb5%G(7)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 1)=nPIV; p0%C( 3)=nVIA
        us=-1.0d0; PV=p0%F(1); PVup=p0%octNb6%F(1)
        VIA=p0%G(7); VIAup=p0%octNb6%G(7)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 2)=nPIV; p0%C( 4)=nVIA
        ! Ex 2nd
        fi=p0%G(1); fiup=p0%octNb6%G(1); fid=p0%octNb5%G(1)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 5)=fn
        fi=p0%G(4); fiup=p0%octNb6%G(4); fid=p0%octNb5%G(4)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 6)=fn
        ! By 1st
        us= 1.0d0; PV=p0%F(5); PVup=p0%octNb5%F(5)
        VIA=p0%G(17); VIAup=p0%octNb5%G(17)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 9)=nPIV; p0%C(11)=nVIA
        us=-1.0d0; PV=p0%F(5); PVup=p0%octNb6%F(5)
        VIA=p0%G(17); VIAup=p0%octNb6%G(17)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C(10)=nPIV; p0%C(12)=nVIA
        ! By 2nd
        fi=p0%G(11); fiup=p0%octNb6%G(11); fid=p0%octNb5%G(11)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(13)=fn
        fi=p0%G(14); fiup=p0%octNb6%G(14); fid=p0%octNb5%G(14)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(14)=fn
     endif
  enddo
  ! -- reconstruct --
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        p0%F( 1)=0.5d0*((p0%C(1)+p0%C(2))+(p0%C( 9)-p0%C(10)))
        p0%G( 7)=0.5d0*((p0%C(3)+p0%C(4))+(p0%C(11)-p0%C(12)))
        p0%G( 1)=p0%C(5)
        p0%G( 4)=p0%C(6)
        p0%F( 5)=0.5d0*((p0%C(1)-p0%C(2))+(p0%C( 9)+p0%C(10)))
        p0%G(17)=0.5d0*((p0%C(3)-p0%C(4))+(p0%C(11)+p0%C(12)))
        p0%G(11)=p0%C(13)
        p0%G(14)=p0%C(14)
     endif
  enddo
  ! -- 2nd step --
  ! --- calc. numerical flux for Ex,Bz at real oct ---
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        ! Ey 1st
        us= 1.0d0; PV=p0%F(2); PVup=p0%octNb5%F(2)
        VIA=p0%G(8); VIAup=p0%octNb5%G(8)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 1)=nPIV; p0%C( 3)=nVIA
        us=-1.0d0; PV=p0%F(2); PVup=p0%octNb6%F(2)
        VIA=p0%G(8); VIAup=p0%octNb6%G(8)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 2)=nPIV; p0%C( 4)=nVIA
        ! Ey 2nd
        fi=p0%G(2); fiup=p0%octNb6%G(2); fid=p0%octNb5%G(2)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 5)=fn
        fi=p0%G(5); fiup=p0%octNb6%G(5); fid=p0%octNb5%G(5)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C( 6)=fn
        ! Bx 1st
        us= 1.0d0; PV=p0%F(4); PVup=p0%octNb5%F(4)
        VIA=p0%G(16); VIAup=p0%octNb5%G(16)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C( 9)=nPIV; p0%C(11)=nVIA
        us=-1.0d0; PV=p0%F(4); PVup=p0%octNb6%F(4)
        VIA=p0%G(16); VIAup=p0%octNb6%G(16)
        call CIP1D(PV,PVup,VIA,VIAup,us,nPIV,nVIA,dstep,dts,dxs)
        p0%C(10)=nPIV; p0%C(12)=nVIA
        ! Bx 2nd
        fi=p0%G(10); fiup=p0%octNb6%G(10); fid=p0%octNb5%G(10)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(13)=fn
        fi=p0%G(13); fiup=p0%octNb6%G(13); fid=p0%octNb5%G(13)
        call linear(fi,fiup,fid,fn,dstep,dts,dxs)
        p0%C(14)=fn
     endif
  enddo
  ! -- reconstruct --
  do index=MinID(1,iLv),MaxID(1,iLv)
     p0 => Mesh(index)
     if((p0%iFLG(1)>=iFs).and.(p0%iFLG(1)<=iFe)) then
        p0%F( 2)=0.5d0*((p0%C(1)+p0%C(2))-(p0%C( 9)-p0%C(10)))
        p0%G( 8)=0.5d0*((p0%C(3)+p0%C(4))-(p0%C(11)-p0%C(12)))
        p0%G( 2)=p0%C(5)
        p0%G( 5)=p0%C(6)
        p0%F( 4)=0.5d0*((-p0%C(1)+p0%C(2))+(p0%C( 9)+p0%C(10)))
        p0%G(16)=0.5d0*((-p0%C(3)+p0%C(4))+(p0%C(11)+p0%C(12)))
        p0%G(10)=p0%C(13)
        p0%G(13)=p0%C(14)
     endif
  enddo
  return
end subroutine advection_z

subroutine CIP1D(PV,PVup,VIA,VIAup,us,nPV,nVIA,dstep,dts,dxs)
  implicit none
  real(kind=8),intent(in)::PV,PVup,VIA,VIAup,us,dstep,dts,dxs
  real(kind=8)::nPV,nVIA
  real(kind=8)::xii,D,ai,bi
  real(kind=8)::Si,BB,alp,ci
  xii=-us*dstep*dts
  D=-dsign(1.0d0,us)*dxs
!  ai=(VIA+VIAup)/D**2+(2.0d0*(PV-PVup))/D**3
!  bi=(3.0d0*(PVup-PV))/D**2-(2.0d0*VIA+VIAup)/D
!  nPV=ai*xii**3+bi*xii**2+VIA*xii+PV
!  nVIA=3.0d0*ai*xii**2+2.0d0*bi*xii+VIA
  Si=(PVup-PV)/D
  if((VIAup-Si)==0.0d0) then
     BB=-1.0d0/D
  else
     BB=(dabs((Si-VIA)/(VIAup-Si+1.0d-10))-1.0d0)/D
  endif
  if((Si-VIA)/(VIAup-Si)>=0.0d0) then
     alp=1.0d0
  else
     alp=0.0d0
  endif
  ai=(VIA-Si+(VIAup-Si)*(1.0d0+alp*BB*D))/D**2
  bi=Si*alp*BB+(Si-VIA)/D-ai*D
  ci=VIA+PV*alp*BB
  nPV=(ai*xii**3+bi*xii**2+ci*xii+PV)/(1.0d0+alp*BB*xii)
  nVIA=(3.0d0*ai*xii**2+2.0d0*bi*xii+ci)/(1.0d0+alp*BB*xii) &
       -(ai*xii**3+bi*xii**2+ci*xii+PV)*alp*BB/(1.0d0+alp*BB*xii)**2
  return
end subroutine CIP1D

subroutine linear(fi,fiup,fid,fn,dstep,dts,dxs)
  implicit none
  real(kind=8),intent(in)::fi,fiup,fid,dstep,dts,dxs
  real(kind=8)::xii,fn
  xii=dstep*dts
  fn=fi+(fiup-fid)*xii/dxs*0.50d0
  return
end subroutine linear

subroutine makeSatelliteTableArray(size)
    use param
    use oct_set
    use init_mesh_size
    use message_passing_interface
    use injectionOct
    implicit none

    integer(kind=4)::i=1,index,size
    type(oct),pointer::p0
    do index=MinID(1,0),MaxID(1,0)
        p0=>Mesh(index)
        if(p0%iPos(1).le.76.and.p0%iPos(1).ge.60&
        &.and.p0%iPos(2).ge.348.and.p0%iPos(2).le.676&
        &.and.p0%iPos(3).ge.348.and.p0%iPos(3).le.676)then
            satelliteTableArray(i)=p0%octN
            i=i+1
        end if
    end do
    size=i-1
end subroutine makeSatelliteTableArray

subroutine setup_parameter
  use param
  use oct_set
  use particle_set
  use init_mesh_size
  use const
  use init_condition
  use mtmod

  implicit none
  integer(kind=4)::Isort
  real(kind=8)   ::rlol
  real(kind=8)   ::md,Ni,mu0
  real(kind=8)   ::Btheta,Bphi,BIMF

! --Rep particle 
  Isort = 0
  Rcharge(Isort) =  0.d0
  Rmass(Isort)   =  0.d0
! -- Electrons, particles of sort #1 --
  isort = 1
  Rcharge(isort) = -1.d0
  Rmass(isort)   =  1.d0
  do isort=2,IonSorts
! -- Ions,  particles of ion species #2 --
     if(isort==2) then
        if(Model==3)then
           Rcharge(isort) = 1.d0
           Rmass(isort)   = 1.d0
        else
           Rcharge(isort) = 1.d0
           Rmass(isort)   = 10000.d0
        end if
     endif
  enddo

  mu0   = 1.0d0
  Ni    = 1.0d0
  rlol  = 8.0d0
  Bphi  =(  0.d0/360.d0)*(2.0d0*PI)
  Btheta=(  0.d0/360.d0)*(2.0d0*PI)
  
  flow(1)=injectionFlow
  flow(2)=0.d0
  flow(3)=0.d0
!  flow(1)=injectionFlow*dsin(Ftheta)*dcos(Fphi)
!  flow(2)=injectionFlow*dsin(Ftheta)*dsin(Fphi)
!  flow(3)=injectionFlow*dcos(Ftheta)

  BIMF=omegace*Rmass(1)/Rcharge(1)
  FIMF(4)=BIMF*dsin(Btheta)*dcos(Bphi)
  FIMF(5)=BIMF*dsin(Btheta)*dsin(Bphi)
  FIMF(6)=BIMF*dcos(Btheta)
  FIMF(1)=-(flow(2)*FIMF(6)-flow(3)*FIMF(5))
  FIMF(2)=-(flow(3)*FIMF(4)-flow(1)*FIMF(6))
  FIMF(3)=-(flow(1)*FIMF(5)-flow(2)*FIMF(4))
  
  TVION=TVELE/(dsqrt(Rmass(2)))
  md=(4.d0 * PI * Rmass(2)**2 * injectionFlow) / (rlol**3 * Rcharge(2)**3 * mu0**2 * Ni)
  mx=0.d0
  my=0.d0
  mz=md
  
end subroutine setup_parameter

subroutine output_histogram2(Hmax,Hmin,HDnum,tgt,iLvMax,iLvMin)
  use oct_set
  use particle_set
  use init_mesh_size
  use const
  use param
  use message_passing_interface
  implicit none
  integer(kind=4),intent(in)::HDnum,iLvMax,iLvMin
  real(kind=8),intent(in)::Hmax,Hmin
  character(LEN=2),intent(in)::tgt
  integer(kind=4)::tgt_num
  character(LEN=128)::filename
  type(oct),pointer::p0
  type(prtcl),pointer::pp
  integer(kind=4)::iLv,maxindex,minindex,index,i,octP,IPinit
  integer(kind=4),dimension(HDnum)::Histogram,Htemp
  real(kind=8)::HistD,HistDrev
  real(kind=8)::hist_d
  real(kind=8)::dataMax,dataMin
  integer(kind=8)::dist

  if(iLvMax>LvMax)then
     print *,"Argument error iLvMax",iLvMax
     return
  elseif(iLvMin<0)then
     print *,"Argument error iLvMin",iLvMin
     return
  endif

  if(tgt/="Vx" .and. tgt/="Vy" .and. tgt/="Vz" .and. &
     tgt/="Px" .and. tgt/="Py" .and. tgt/="Pz")then
     print *,"Unknown output type ",tgt
     return
  endif

  if(tgt=="Vx")tgt_num=4
  if(tgt=="Vy")tgt_num=5
  if(tgt=="Vz")tgt_num=6
  if(tgt=="Px")tgt_num=1
  if(tgt=="Py")tgt_num=2
  if(tgt=="Pz")tgt_num=3


  if(Hmax<=Hmin)then
     print *,"Illegal histogram range Hmax=",Hmax,"Hmin=",Hmin
     return
  endif

  HistD=(Hmax-Hmin)/dble(HDnum)
  HistDrev=1.d0/HistD
  Histogram=0
  Htemp=0
  dataMax=-999
  dataMin=999

  if(tgt=="Vx" .or. tgt=="Vy" .or. tgt=="Vz" .or. tgt=="Px" .or. tgt=="Py" .or. tgt=="Pz")then
     do iLv=iLvMin,iLvMax
        minindex=MinID(1,iLv)
        maxindex=MaxID(1,iLv)
        if(minindex>=maxindex)cycle
        do index=1,tableSize
           IPinit=InjectionTable(index)
           p0=>Mesh(IPinit)
           if(p0%iFLG(1)<0 .or. p0%iFLG(1)>3)cycle
           if(p0%octP==0)cycle
           octP=p0%octP
           pp=>p0%ptcl
           do i=1,octP
              pp=>pp%prtnxt
              if(pp%isort/=1)cycle
              hist_d=(pp%R(tgt_num)-Hmin)*HistDrev
              dist=int(hist_d)
              if(dist<1)dist=1
              if(dist>HDnum)dist=HDnum
              Histogram(dist)=Histogram(dist)+1
              dataMax=max(dataMax,pp%R(tgt_num))
              dataMin=min(dataMin,pp%R(tgt_num))
           end do
        end do
     enddo
  endif

!sum data using mpi
  call mpi_reduce(Histogram,Htemp,HDnum,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  Histogram=Htemp
  if(rank==0)then
     write(filename,2513)version,tgt,istep
     open(77,file=filename,form='formatted')
     do i=1,HDnum
        write(77,*)Hmin+HistD*dble(i),Histogram(i)    
     enddo
     close(77)
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  return 

  2513 format("./output/",A,"_Hist_injection_",A,i5.5,".dat")
end subroutine output_histogram2

subroutine calculate_V(index)
    use oct_set
    use param
    use const
    use message_passing_interface
    use init_mesh_size
    implicit none
    integer(kind=4),intent(in)::index
    integer(kind=4)::Ecount,Icount,i
    real(kind=8)::Ve(1:6),Vi(1:6),oF(1:6)
    type(oct),pointer::p0
    type(prtcl),pointer::pp
        
    p0=>Mesh(index)
    pp=>p0%ptcl
    Ve(:)=0.d0
    Vi(:)=0.d0
    oF(:)=0.d0
    p0%O(1:6)=0.d0
    p0%C(4:9)=0.d0
    Ecount=0
    Icount=0
    do i=1,p0%octP
        pp=>pp%prtnxt
        if(pp%Isort==1)then
            Ve(1)=Ve(1)+pp%R(4)
            Ve(2)=Ve(2)+pp%R(5)
            Ve(3)=Ve(3)+pp%R(6)
            Ecount=Ecount+1
        else if(pp%Isort==2)then
            Vi(1)=Vi(1)+pp%R(4)
            Vi(2)=Vi(2)+pp%R(5)
            Vi(3)=Vi(3)+pp%R(6)
            Icount=Icount+1
        end if
    end do
    if(Ecount.ne.0)then
        p0%O(1)=Ve(1)/Ecount
        p0%O(2)=Ve(2)/Ecount
        p0%O(3)=Ve(3)/Ecount
    end if
    if(Icount.ne.0)then
        p0%O(4)=Vi(1)/Icount
        p0%O(5)=Vi(2)/Icount
        p0%O(6)=Vi(3)/Icount
    end if
    pp=>p0%ptcl
    do i=1,p0%octP
        pp=>pp%prtnxt
        if(pp%Isort==1)then
            Ve(4)=Ve(4)+(pp%R(4)-p0%O(1))**2
            Ve(5)=Ve(5)+(pp%R(5)-p0%O(2))**2
            Ve(6)=Ve(6)+(pp%R(6)-p0%O(3))**2
        else if(pp%Isort==2)then
            Vi(4)=Vi(4)+(pp%R(4)-p0%O(4))**2
            Vi(5)=Vi(5)+(pp%R(5)-p0%O(5))**2
            Vi(6)=Vi(6)+(pp%R(6)-p0%O(6))**2
        end if
    end do
    if(Ecount.ne.0)then
        p0%C(4)=dsqrt(Ve(4)/Ecount)
        p0%C(5)=dsqrt(Ve(5)/Ecount)
        p0%C(6)=dsqrt(Ve(6)/Ecount)
    end if
    if(Icount.ne.0)then
        p0%C(7)=dsqrt(Vi(4)/Icount)
        p0%C(8)=dsqrt(Vi(5)/Icount)
        p0%C(9)=dsqrt(Vi(6)/Icount)
    end if
end subroutine calculate_V
      
subroutine calculate_gradB(index)
    use oct_set
    use param
    use const
    use message_passing_interface
    use init_mesh_size
    implicit none
    integer(kind=4)::index
    real(kind=8)::Bpx,Bpy,Bpz,Bmx,Bmy,Bmz,Bm,Bp,Bx,By,Bz,Fx,Fy,Fz,tempF
    real(kind=8)::BB,Ve,Vi,te,ti,dxx
    type(oct),pointer::p0
  
    p0=>Mesh(index)

    Bpx=p0%octNb2%F(4)+p0%octNb2%D(1)+FIMF(4)
    Bpy=p0%octNb2%F(5)+p0%octNb2%D(2)+FIMF(5)
    Bpz=p0%octNb2%F(6)+p0%octNb2%D(3)+FIMF(6)
    Bp=dsqrt(Bpx**2+Bpy**2+Bpz**2)
    Bmx=p0%octNb1%F(4)+p0%octNb1%D(1)+FIMF(4)
    Bmy=p0%octNb1%F(5)+p0%octNb1%D(2)+FIMF(5)
    Bmz=p0%octNb1%F(6)+p0%octNb1%D(3)+FIMF(6)
    Bm=dsqrt(Bmx**2+Bmy**2+Bmz**2)
    tempF=Bp-Bm
    dxx=p0%octNb2%rPos(1)-p0%octNb1%rPos(1)
    Fx=tempF/dxx
    
    Bpx=p0%octNb4%F(4)+p0%octNb4%D(1)+FIMF(4)
    Bpy=p0%octNb4%F(5)+p0%octNb4%D(2)+FIMF(5)
    Bpz=p0%octNb4%F(6)+p0%octNb4%D(3)+FIMF(6)
    Bp=dsqrt(Bpx**2+Bpy**2+Bpz**2)
    Bmx=p0%octNb3%F(4)+p0%octNb3%D(1)+FIMF(4)
    Bmy=p0%octNb3%F(5)+p0%octNb3%D(2)+FIMF(5)
    Bmz=p0%octNb3%F(6)+p0%octNb3%D(3)+FIMF(6)
    Bm=dsqrt(Bmx**2+Bmy**2+Bmz**2)
    tempF=Bp-Bm
    dxx=p0%octNb4%rPos(2)-p0%octNb3%rPos(2)
    Fy=tempF/dxx
    
    Bpx=p0%octNb6%F(4)+p0%octNb6%D(1)+FIMF(4)
    Bpy=p0%octNb6%F(5)+p0%octNb6%D(2)+FIMF(5)
    Bpz=p0%octNb6%F(6)+p0%octNb6%D(3)+FIMF(6)
    Bp=dsqrt(Bpx**2+Bpy**2+Bpz**2)
    Bmx=p0%octNb5%F(4)+p0%octNb5%D(1)+FIMF(4)
    Bmy=p0%octNb5%F(5)+p0%octNb5%D(2)+FIMF(5)
    Bmz=p0%octNb5%F(6)+p0%octNb5%D(3)+FIMF(6)
    Bm=dsqrt(Bmx**2+Bmy**2+Bmz**2)
    tempF=Bp-Bm
    dxx=p0%octNb6%rPos(3)-p0%octNb5%rPos(3)
    Fz=tempF/dxx

    Bx=p0%F(4)+p0%D(1)+FIMF(4)
    By=p0%F(5)+p0%D(2)+FIMF(5)
    Bz=p0%F(6)+p0%D(3)+FIMF(6)
    BB=dsqrt(Bx**2+By**2+Bz**2)

    Ve=dsqrt(p0%O(1)**2+p0%O(2)**2)
    Vi=dsqrt(p0%O(4)**2+p0%O(5)**2)
    
    te=Rmass(1)*Ve**2/(2*Rcharge(1)*BB**3)
    ti=Rmass(2)*Vi**2/(2*Rcharge(2)*BB**3)
    
    p0%O(7) =te*(By*Fz-Bz*Fy)
    p0%O(8) =te*(Bz*Fx-Bx*Fz)
    p0%O(9) =te*(Bx*Fy-By*Fx)
    p0%O(10)=ti*(By*Fz-Bz*Fy)
    p0%O(11)=ti*(Bz*Fx-Bx*Fz)
    p0%O(12)=ti*(Bx*Fy-By*Fx)
end subroutine calculate_gradB

subroutine calculate_V2(index)
    use oct_set
    use param
    use const
    use message_passing_interface
    use init_mesh_size

    implicit none
    integer(kind=4)::index
    real(kind=8)::mq,mqn,B,BB,EE,VV,Vxx,Zxx
    type(oct),pointer::p0
    
    p0=>Mesh(index)
    p0%O(13)=0.d0
    p0%O(14)=0.d0
    p0%O(15)=0.d0
    p0%O(16)=0.d0
    
    B  = p0%F(6)+p0%D(3)+FIMF(6)
    EE = p0%F(1)+FIMF(1)
    mq = Rmass(1)/Rcharge(1)
    if(B==0.d0)then
        BB = 0.d0
    else
        BB = 1.d0/B
    end if
    if(p0%C(2)==0)then
        mqn = 0.d0
    else
        mqn = Rmass(1)/(Rcharge(1)*dble(p0%C(2)))
    end if
    Vxx=(p0%octNb2%C(1)-p0%octNb1%C(1))/(2*dx(1))
    Zxx=(p0%octNb2%C(2)-p0%octNb1%C(2))/(2*dx(1))
    
    p0%O(13)=-EE*BB
    p0%O(14)= mq*BB*p0%C(1)*Vxx
    p0%O(15)= 0.5d0*mqn*TVELE**2*BB*Zxx
    p0%O(16)= p0%O(13)+p0%O(14)+p0%O(15)

end subroutine calculate_V2

subroutine LPF(index)
    use oct_set
    use param
    use const
    use message_passing_interface
    use init_mesh_size

    implicit none
    integer(kind=4)::index
    type(oct),pointer::p0
    
    p0=>Mesh(index)
    p0%C(1)=0.d0
    p0%C(2)=0.d0
    p0%C(1)=0.2d0*(p0%octNb1%octNb1%O(1)+p0%octNb1%O(1)+p0%O(1)+p0%octNb2%O(1)+p0%octNb2%octNb2%O(1))
    p0%C(2)=0.2d0*(p0%octNb1%octNb1%Z(1)+p0%octNb1%Z(1)+p0%Z(1)+p0%octNb2%Z(1)+p0%octNb2%octNb2%Z(1))
end subroutine LPF


subroutine deleteFieldInSatellite
    use param
    use oct_set
    use init_mesh_size
    use injectionOct
    implicit none
    
    integer(kind=4)::index
    
    do index=1,satelliteArraySize
      call deleteFieldInSatelliteEveryOct(satelliteTableArray(index))
    end do
end subroutine deleteFieldInSatellite

recursive subroutine deleteFieldInSatelliteEveryOct(index)
    use param
    use oct_set
    use init_mesh_size
    implicit none
    
    integer(kind=4)::index
    type(oct),pointer::p0
    
    p0=>Mesh(index)
    if (p0%iFLG(1) .ge. 4) then
        call deleteFieldInSatelliteEveryOct(p0%octCh1%octN)
        call deleteFieldInSatelliteEveryOct(p0%octCh2%octN)
        call deleteFieldInSatelliteEveryOct(p0%octCh3%octN)
        call deleteFieldInSatelliteEveryOct(p0%octCh4%octN)
        call deleteFieldInSatelliteEveryOct(p0%octCh5%octN)
        call deleteFieldInSatelliteEveryOct(p0%octCh6%octN)
        call deleteFieldInSatelliteEveryOct(p0%octCh7%octN)
        call deleteFieldInSatelliteEveryOct(p0%octCh8%octN)
    else
        p0%F(1) = 0.0d0
        p0%F(2) = 0.0d0
        p0%F(3) = 0.0d0
        p0%F(4) = 0.0d0
        p0%F(5) = 0.0d0
        p0%F(6) = 0.0d0
        p0%G(1:18) = 0.0d0
    end if
end subroutine deleteFieldInSatelliteEveryOct

subroutine eraseElectronTable
    use param
    use oct_set
    use init_mesh_size
    use injectionOct
    implicit none
    
    integer(kind=4)::index
    
    do index=1,arraySize
      if (InjectionTableArray(index)%iSort.eq.1) then
        InjectionTableArray(index)%iSort = 0
      endif
    enddo
end subroutine eraseElectronTable
	
subroutine add_Voct
		use param
		use init_mesh_size
        use oct_set
        use particle_set
        use Vparam
        implicit none
		integer(kind=4)::index,indexP
		integer(kind=4)::octN,octLv,octP,Csort
		integer(kind=4) :: Isort
		integer(kind=4)::iPos(3)
		real(kind=8)   ::rPos(3)
		integer(kind=4)::iFLG(3)
		integer(kind=4)  :: MrtN      
		real(kind=8)    :: F(18),C(6*iomp0),Z(IonSorts),R(9),G(1:18),D(1:3),O(1:16)
        integer(kind=4)  :: iC(iomp0) 
        integer(kind=4) ::  octType, prc_bndry, ptcl_loops
        type(oct) :: p0
        type(prtcl), pointer :: PrtList
		type(prtcl), pointer :: ptcl    ! representative particle (0th)
        type(prtcl), pointer :: ptclA   ! representative particle (0th) ! We don't really use this pointer
		
		indexP = MaxIP(0)
		
		do index=1,VMeshSize
			octN=index	
			octLV=0
			Csort=0
			Isort = -1
			iPos=0
			rPos=0.d0
			MrtN=0
			F = 0.d0
			R = 0.d0
			C = 0.d0
			Z = 0.d0
			G = 0.d0
			D=0.d0
			O=0.d0
			octP= 0
			iC  = 0
			octType = 0
			prc_bndry = 0
			ptcl_loops=0
			iFLG=0
			nullify(PrtList)
			nullify(p0%octPrt)
			nullify(p0%octNb1) ; nullify(p0%octNb2)
			nullify(p0%octNb3) ; nullify(p0%octNb4)
			nullify(p0%octNb5) ; nullify(p0%octNb6)
			nullify(p0%octCh1) ; nullify(p0%octCh2)
			nullify(p0%octCh3) ; nullify(p0%octCh4)
			nullify(p0%octCh5) ; nullify(p0%octCh6)
			nullify(p0%octCh7) ; nullify(p0%octCh8)
			nullify(p0%Psort)
			nullify(p0%ptclA)
          			
			call get_VInverseIndex(octN,iPos)
			rPos(1)=dble(iPos(1)-1)*dx(1)+0.5d0*dx(1)
			rPos(2)=dble(iPos(2)-1)*dx(2)+0.5d0*dx(2)
			rPos(3)=dble(iPos(3)-1)*dx(3)+0.5d0*dx(3)
            R(1)=rPOS(1)
            R(2)=rPOS(2)
            R(3)=rPOS(3)
			indexP = indexP + index
			Pesh(indexP,0) =                  &
             prtcl(R,Isort,octN,indexP,PrtList,0)
            p0%ptcl => Pesh(indexP,0)				 
			
			VMesh(index)=oct(octN,octLv,octP,Csort,iFLG,     &
                             iPOS,rPOS,                       &
                             MrtN,iC,                    &
                             prc_bndry,                  &
                             octType,                      &
                             F,C,Z,G,D,O,ptcl_loops,    & 
                             p0%octPrt,                     &
                             p0%octNb1, p0%octNb2,        &
                             p0%octNb3, p0%octNb4,        &
                             p0%octNb5, p0%octNb6,        &
                             p0%octCh1, p0%octCh2,        &
                             p0%octCh3, p0%octCh4,        &
                             p0%octCh5, p0%octCh6,        &
                             p0%octCh7, p0%octCh8,        &
                             p0%Psort ,                     &
                             p0%ptcl  , p0%ptclA)

		end do
    MaxIP(0) = indexP

end subroutine add_Voct

subroutine connect_octNb
	use Vparam
	use oct_set
	implicit none
	type(oct),pointer :: p0	
	integer(kind=4)::index,Nindex
	integer(kind=4)::iPos(3)
	
	do index=1,VMeshSize
	p0 => VMesh(index)
     iPOS=p0%iPOS
     !x_length is not bigger than 2 
     p0%octNb1=>p0
     
     p0%octNb2=>p0

     if(iPOS(2)==1)then
        call get_VIndex(iPOS(1)  ,       VNY,iPOS(3)  ,Nindex)
     else
        call get_VIndex(iPOS(1)  ,iPOS(2)-1,iPOS(3)  ,Nindex)
     end if
     p0%octNb3=>VMesh(Nindex)

     if(iPOS(2)==VNY)then
        call get_VIndex(iPOS(1)  ,        1,iPOS(3)  ,Nindex)
     else
        call get_VIndex(iPOS(1)  ,iPOS(2)+1,iPOS(3)  ,Nindex)
     end if
     p0%octNb4=>VMesh(Nindex)
     
     if(iPOS(3)==1)then
        call get_VIndex(iPOS(1)  ,iPOS(2)  ,       VNZ,Nindex)
     else
        call get_VIndex(iPOS(1)  ,iPOS(2)  ,iPOS(3)-1,Nindex)
     end if
     p0%octNb5=>VMesh(Nindex)
     
     if(iPOS(3)==VNZ)then
        call get_VIndex(iPOS(1)  ,iPOS(2)  ,        1,Nindex)
     else
        call get_VIndex(iPOS(1)  ,iPOS(2)  ,iPOS(3)+1,Nindex)
     end if
     p0%octNb6=>VMesh(Nindex)
  end do
  
end subroutine connect_octNb

subroutine get_VInverseIndex(MrtN,iPOS)
  use Vparam
  implicit none
  integer(kind=4),intent( in)::MrtN
  integer(kind=4),intent(out)::iPOS(3)
  integer(kind=4)::index
  index=MrtN-1
  iPOS(1)=mod(index,VNX)+1
  iPOS(2)=mod(int(index/VNX),VNY)+1
  iPOS(3)=mod(int(index/(VNX*VNY)),VNZ)+1
  return
end subroutine get_VInverseIndex

subroutine get_VIndex(iPOSX,iPOSY,iPOSZ,index)
  use Vparam
  implicit none
  integer(kind=4),intent( in)::iPOSX,iPOSY,iPOSZ
  integer(kind=4),intent(out)::index
  
  index=iPOSX + (iPOSY-1)*VNX + (iPOSZ-1)*VNX*VNY
  return
end subroutine get_VIndex

subroutine makeVoct
  use init_mesh_size
  use Vparam
  use particle_set
  use oct_set
  implicit none
  allocate(VMesh(VMeshSize))
  call add_Voct
  call connect_octNb
  return 
end subroutine makeVoct

subroutine move_Vparticle
        use message_passing_interface
        use oct_set
        use particle_set
        use init_mesh_size
        use const
		use Vparam
        use param
        implicit none
        real(kind=8)    :: dvlevel, deltaT,dvlevel3
        integer(kind=4)      :: index
		integer(kind=4) :: iLv=0
        integer(kind=4)      :: Nin
        real(kind=4)         :: Rinin
        real(kind=8)         :: R1(10)
        integer(kind=4)      :: k,kp,indexP,particle
        type(oct), pointer   :: p0,p1
        type(prtcl), pointer :: pptc
        real(kind=8)         :: R(23),R0(-2:20,3)
!------------------------
!!$        include 'omp_lib.h'
!
! == BEGIN: ============================================================

        if((maxID(1,iLv).eq.minID(1,iLv)) .and. (maxID(3,iLv).eq.minID(3,iLv))) return 

!
!=======================================================
        R1=0.d0
! -- Set spacing parameters --
        dvlevel = 0.5d0**(iLv)
        dvlevel3=(dvlevel**3)
        deltaT  = dvlevel * dt
        Nin   = 2**(LvMax-iLv)
        Rinin = 1.d0/real(2*Nin)
        PICnumber(1) = npart_per_cell
        PICnumber(2) = npart_per_cell
        R1( 1) = 0.5d0*QM*deltaT*Rcharge(1)/Rmass(1)*omega
        R1( 2) = Rcharge(1)/dvlevel3*(-1.d0)/dble(PICnumber(1))
        R1( 3) = ONE/(dvlevel*dx(1))
        R1( 4) = 0.5d0*QM*deltaT*Rcharge(2)/Rmass(2)*omega
        R1( 5) = Rcharge(2)/dvlevel3*(-1.d0)/dble(PICnumber(2))
        R1( 6) = ONE/(dvlevel*dx(2))
        R1( 7)=  Rinin
        R1( 8)=  deltaT
        R1( 9)=-dx(1)*dx(2)*BD
        R1(10) = ONE/(dvlevel*dx(3))

          R0(-1,1)=R1( 1) ; R0( 0,1)=R1( 2) ; R0( 1,1)=R1( 3)
          R0(-1,2)=R1( 4) ; R0( 0,2)=R1( 5) ; R0( 1,2)=R1( 6)
          R0(-2,3)=R1( 7) ; R0(-1,3)=R1( 8) ; R0( 0,3)=R1( 9) ; R0( 1,3)=R1(10)        


!!$omp parallel do private(index   &
!!$omp&        ,R0,indexP,pptc,kp,R,p1 & 
!!$omp&        ,particle) &
!!$omp&             shared(R1)        
!!!$omp&           schedule(dynamic,2)

		do index = 1 , VMeshSize
			p1 => VMesh(index)
				indexP =  p1%octP 
				if(indexP==0) cycle
				
				pptc   => p1%ptcl%prtnxt
                do kp=1,indexP
                      R( 1)= pptc%R(1)
                      R( 2)= pptc%R(2)
                      R( 3)= pptc%R(3)
                      R( 4)= pptc%R(4)
                      R( 5)= pptc%R(5)
                      R( 6)= pptc%R(6)
                      R( 7)= 0.d0
                      R( 8)= 0.d0
                      R( 9)= 0.d0
                      R(10)= 0.d0
                      R(11)= 0.d0
                      R(12)= 0.d0
                      R(13)= R(4)
                      R(14)= R(5)
                      R(15)= R(6)
                      R(16)= R0(-1,pptc%Isort)
                      R(23)= R0( 0,pptc%Isort)
                      particle= pptc%Isort

! =====================================================
!    Calculation of particle movement
! =====================================================                    
                     
                      R( 7)=R( 7)*R(16)
                      R( 8)=R( 8)*R(16)
                      R( 9)=R( 9)*R(16)
                      R(10)=R(10)*R(16)
                      R(11)=R(11)*R(16)
                      R(12)=R(12)*R(16)

! ---  u^- = u^{n-1/2} + (qdt/2m)*E  ----
                      R(17)=R(13)+R(7)
                      R(18)=R(14)+R(8)
                      R(19)=R(15)+R(9)
!
! ---  gamma = sqrt{ 1 + |u^-|^2 }  ---
                      R(20)=1.d0/SQRT( ONE + R(17)*R(17) + R(18)*R(18) + R(19)*R(19) )
!
! ---  h = (qdt/2m)*(1/gamma)*B  ---

                      R(10)=R(10)*R(20)
                      R(11)=R(11)*R(20)
                      R(12)=R(12)*R(20)
!
! ---  coefficients for u^{n+1/2}  ---
                      R(20) = R(10)*R(10) + R(11)*R(11) + R(12)*R(12)
                      R(21) = (ONE - R(20))
                      R(22) = ONE/(ONE + R(20))
                      R(20) = R(10)*R(17) + R(11)*R(18) + R(12)*R(19) 
!
! ---  update velocities
!      u^{n+1/2} = (1-h^2)/(1+h^2)*u^-
!                + 2/(1+h^2)*( (h,u^-) + [u^- x h] ) + (qdt/2m)*E
                      R(13) = R(7) + R(22)*(R(21)*R(17) + 2.d0*(R(20)*R(10) + R(18)*R(12)-R(19)*R(11)) ) 
                      R(14) = R(8) + R(22)*(R(21)*R(18) + 2.d0*(R(20)*R(11) + R(19)*R(10)-R(17)*R(12)) ) 
                      R(15) = R(9) + R(22)*(R(21)*R(19) + 2.d0*(R(20)*R(12) + R(17)*R(11)-R(18)*R(10)) )

!
! ---  update coordinates
!      x^{n+1} = x^{n} + u^{n+1/2}*dt/gamma^{n+1/2}
                      R(20)= R0(-1,3)/SQRT( ONE + R(13)*R(13) + R(14)*R(14) + R(15)*R(15) )
                      R(1) = R(1) + R(13)*R(20)
                      R(2) = R(2) + R(14)*R(20)
                      R(3) = R(3) + R(15)*R(20) 

!======================important==========================
                      pptc%R(1) = R(1)
                      pptc%R(2) = R(2)
                      pptc%R(3) = R(3)
                      pptc%R(4) = R(13)
                      pptc%R(5) = R(14)
                      pptc%R(6) = R(15)
                      pptc%R(7) = (R(13)-R(4))*0.5d0+R(13)
                      pptc%R(8) = (R(14)-R(5))*0.5d0+R(14)
                      pptc%R(9) = (R(15)-R(6))*0.5d0+R(15)
                      pptc => pptc%prtNxt
!=========================================================
                   end do ! particle loop
            end do  ! oct loop   
       return
       
end subroutine move_Vparticle


subroutine makeInjectionOctTableArray(size)
    use param
    use oct_set
    use init_mesh_size
    use message_passing_interface
    use injectionOct
    use Vparam
    implicit none
    integer(kind=4)::i=1,index,size
    type(oct),pointer::p0,p1
    call makeVoct
    do index=MinID(1,0),MaxID(1,0)
        p0=>Mesh(index)
        if(i .le. VMeshSize)then
            p1=>VMesh(i)
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==468.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==468.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==468.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==468.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==468)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==468)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==468)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==468)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==468.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==468.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==468.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==468.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==476.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==484.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==492.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==556)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==500.and.p0%iPos(3)==556)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==556)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==508.and.p0%iPos(3)==556)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==468)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==468)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==468)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==468)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==476)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==484)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==492)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==556.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==556.and.p0%iPos(3)==500)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==556.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==556.and.p0%iPos(3)==508)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==556)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==516.and.p0%iPos(3)==556)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==556)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==524.and.p0%iPos(3)==556)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==532.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==540.and.p0%iPos(3)==548)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==532)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==548.and.p0%iPos(3)==540)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==556.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==556.and.p0%iPos(3)==516)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==556.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 2, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
        if(p0%iPos(1)==76.and.p0%iPos(2)==556.and.p0%iPos(3)==524)then
            InjectionTableArray(i)=injectionStructure(p1%octN, 1, 0d0, 1.5707963267949d0, 1d0)
            p1%rPos = p0%rPos
            p1%octNb2 => p0%octNb2
            i = i + 1
        end if
    end do
    size = i - 1
end subroutine makeInjectionOctTableArray
