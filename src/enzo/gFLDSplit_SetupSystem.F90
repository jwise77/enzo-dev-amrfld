#include "fortran.def"
#include "phys_const.def"
!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================


function limiter(E1, E2, k1, k2, nUn, lUn, tUn, dxi)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       February 2013
!  modified:   
!
!  PURPOSE: Computes the flux limiter at a given face
!=======================================================================
  implicit none
  real*8, intent(in) :: E1, E2, k1, k2, nUn, lUn, tUn, dxi
  real*8 :: limiter, Eavg, kap, R, Emin, Rmin, Dmax
  
  ! set limiter bounds
!  Rmin = 1.d-20/lUn
  Rmin = 1.d0/lUn
!  Rmin = 1.d-20
  Emin = 1.d-30
  Dmax = 0.0021565d0 * c_light * lUn  +  0.0167231d0 * lUn * lUn / tUn
!  Dmax = 2.0539e-3 * c_light * lUn
!  Dmax = 1.d300

  ! compute limiter
  Eavg = max((E1 + E2)*0.5d0, Emin)
  kap = 2.d0*k1*k2/(k1+k2)*nUn        ! harmonic average
!!$  kap = (k1 + k2)*0.5d0*nUn           ! arithmetic average
!!$  kap = sqrt(k1)*sqrt(k2)*nUn         ! geometric average
  R = max(dxi*abs(E1 - E2)/Eavg, Rmin)
  limiter = min(c_light/sqrt(9.d0*kap*kap + R*R), Dmax)

end function limiter
!=======================================================================



subroutine gFLDSplit_SetupSystem(matentries, rhsentries, rhsnorm, E0,   &
     E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,     &
     theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn, rank, dx, dy, dz,     &
     BCXl, BCXr, BCYl, BCYr, BCZl, BCZr, x0s, x0e, x1s, x1e, x2s, x2e,  &
     Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, xlface, xrface,    &
     ylface, yrface, zlface, zrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!  modified:   
!
!  PURPOSE: Computes the array of matrix stencil elements and vector of 
!           rhs entries for the Grey FLD radiation problem,
!              d_t E - Div(D(E)*Grad(E)) = -adot/a*E - c*kappa*E + eta + src
!           where D(E) is a nonlinear flux-limiter 
!           depending on E0 (time lagged).  We define the values
!              R_i = |Grad(E)_i|/E,
!           The '_i' subscript implies the gradient in the ith 
!           direction; these quantities are all required at cell faces, 
!           as that is the location of the divergence calculations.
!           With these components, we allow any of the following three 
!           forms of the limiter, 
!             [Levermore-Pomraning, 1981],
!                 D_i(E) = c/kappa/R_i*[coth(R_i)-1/R_i],
!             [rational approx. to above, Levermore-Pomraning, 1981],
!                 D_i(E) = c/kappa*(2+R_i)/(6+3*R_i+R_i**2),
!             [Reynolds approximation to LP],
!                 D_i(E) = 2/pi*c*atan(R_i*pi/6/kappa)/R_i
!             [Zeus form of rational approx. to LP],
!                 D_i(E) = c*(2*kappa+R_i)/(6*kappa*kappa+3*kappa*R_i+R_i**2)
!           where we have the [small] parameter
!              kappa = absorption coefficient.
!           As the stencil has {7,5,3} non-zero elements per matrix row 
!           (depending on whether the problem is 3D, 2D or 1D), we 
!           set these entries over the computational domain, with the 
!           proper adjustments due to the choice of limiter.
!
!           We in fact solve a scaled version of the equation.  Since 
!           the values of E are in fact in normalized units 
!           (E_true = E*rUn), we must scale src by rUn to achieve the
!           correct equation.  Moreover, we do not solve the equation 
!           directly, and instead solve for a correction to the current
!           state such that the corrected solution satisfies the above 
!           equation.  This helps with enforcement of boundary conditions,
!           since they may be directly placed onto the current state, 
!           and the correction need only refrain from interfering.
!
!  INPUTS:
!     E0         - Grey radiation energy density (prev time step)
!     E          - Grey radiation energy density (current guess)
!     Temp       - gas temperature for black-body radiation
!     Temp0      - gas temperature (prev time step)
!     kappa      - opacity array
!     src        - spatially-dependent radiation source
!     dt         - time step size
!     a,a0       - cosmological expansion factor (new and old time steps)
!     adot,adot0 - da/dt and da0/dt
!     ESpectrum  - flag denoting what type of radiation field we have:
!                  (-1=>monochromatic)
!     theta      - overall implicitness parameter
!     *Un,*Un0   - variable scaling constants (new and old time steps)
!     rank       - 1, 2 or 3; the dimensionality of the problem
!     dx,dy,dz   - mesh spacing in each direction
!     BC*        - boundary condition type in each direction, face
!                     0->periodic
!                     1->Dirichlet
!                     2->Neumann
!     x*{s,e}    - start/end indices of linear solver domain; 
!                  typically 1:Nx for standard dims, but Dirichlet 
!                  BCs may move these to 0:Nx, 1:Nx+1, etc.
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     *{l,r}face - integer flag denoting whether direction/face 
!                  is external to the domain (0->int, 1->ext)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     matentries - array of stencil values over the active domain.  
!                  Since the stencil has 7 nonzero entries, and as 
!                  this array should not include ghost cells, it 
!                  has dimensions (7,x0s:x0e,x1s:x1e,x2s:x2e).
!     rhsentries - array of rhs values over the active domain.  As 
!                  this array should not include ghost cells, it 
!                  has dimensions (x0s:x0e,x1s:x1e,x2s:x2e)
!     rhsnorm    - 2-norm of rhs array
!     ier        - success/failure flag (0->failure, 1->success)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
  implicit none
  
  !--------------
  ! argument declarations
  integer, intent(in)  :: rank, ESpectrum
  integer, intent(in)  :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer, intent(in)  :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  integer, intent(in)  :: BCZl, BCZr, x2s, x2e, Nz, NGzl, NGzr, zlface, zrface
  REALSUB, intent(in)  :: a, a0, adot, adot0
  real,    intent(in)  :: dt, theta, dx, dy, dz
  real,    intent(in)  :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn
  real,    intent(in)  :: E0(*), E(*), Temp(*), Temp0(*), kappa(*), src(*)
  real*8,  intent(out) :: matentries(*)
  real*8,  intent(out) :: rhsentries(*)
  real,    intent(out) :: rhsnorm
  integer, intent(out) :: ier

  !=======================================================================
  

  ! call the apprpriate routine based on rank
  if (rank == 3) then

     call gFLDSplit_SetupSystem3D(matentries, rhsentries, rhsnorm, E0,    &
          E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,  &
          theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn, dx, dy, dz, BCXl,  &
          BCXr, BCYl, BCYr, BCZl, BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, &
          Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, xlface, xrface,     &
          ylface, yrface, zlface, zrface, ier)

  elseif (rank == 2) then

     call gFLDSplit_SetupSystem2D(matentries, rhsentries, rhsnorm, E0,    &
          E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,  &
          theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn, dx, dy, BCXl,      &
          BCXr, BCYl, BCYr, x0s, x0e, x1s, x1e, Nx, Ny, NGxl, NGxr, NGyl, &
          NGyr, xlface, xrface, ylface, yrface, ier)

  elseif (rank == 1) then

     call gFLDSplit_SetupSystem1D(matentries, rhsentries, rhsnorm, E0,    &
          E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,  &
          theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn, dx, BCXl, BCXr,    &
          x0s, x0e, Nx, NGxl, NGxr, xlface, xrface, ier)

  else
     write(0,*) 'gFLDSplit_SetupSystem error: illegal rank =',rank
  end if

end subroutine gFLDSplit_SetupSystem
!=======================================================================






subroutine gFLDSplit_SetupSystem3D(matentries, rhsentries, rhsnorm, E0, &
     E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,     &
     theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn, dx, dy, dz, BCXl,     &
     BCXr, BCYl, BCYr, BCZl, BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx,    &
     Ny, Nz, NGxl, NGxr, NGyl, NGyr,  NGzl, NGzr, xlface, xrface,       &
     ylface, yrface, zlface, zrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       June 2009
!  modified:   
!
!  PURPOSE: 3D version of the routine
!=======================================================================
  implicit none
  
  !--------------
  ! argument declarations
  integer,  intent(in) :: ESpectrum
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer,  intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  integer,  intent(in) :: BCZl, BCZr, x2s, x2e, Nz, NGzl, NGzr, zlface, zrface
  REALSUB,  intent(in) :: a, a0, adot, adot0
  real,     intent(in) :: dt, theta, dx, dy, dz
  real,     intent(in) :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) &
       :: E0, E, kappa, src, Temp, Temp0
  real*8,  intent(out) :: matentries(7,x0s:x0e,x1s:x1e,x2s:x2e)
  real*8,  intent(out) :: rhsentries(x0s:x0e,x1s:x1e,x2s:x2e)
  real,    intent(out) :: rhsnorm
  integer, intent(out) :: ier

  !--------------
  ! locals
  integer :: i, j, k
  real*8  :: dtfac, dtfac0, kap, kap0, eta, eta0, StBz
  real*8  :: dxi, dxi0, dyi, dyi0, dzi, dzi0, afac, afac0
  real*8  :: D_xl, D0_xl, D_xr, D0_xr, E0d_xl, E0d_xr, Ed_xl, Ed_xr
  real*8  :: D_yl, D0_yl, D_yr, D0_yr, E0d_yl, E0d_yr, Ed_yl, Ed_yr
  real*8  :: D_zl, D0_zl, D_zr, D0_zr, E0d_zl, E0d_zr, Ed_zl, Ed_zr
  real*8, external :: limiter

!=======================================================================

  ! initialize outputs to zero, flag to success
  matentries = 0.d0
  rhsentries = 0.d0
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta         ! time step conversion factor
  dtfac0 = dt*(1.d0-theta)  ! time step conversion factor
  if (ESpectrum == -1) then
     afac  = 0.d0
     afac0 = 0.d0
  else
     afac  = adot/a         ! expansion factor (new time)
     afac0 = adot0/a0       ! expansion factor (old time)
  endif
  dxi   = 1.d0/dx/lUn
  dyi   = 1.d0/dy/lUn
  dzi   = 1.d0/dz/lUn
  dxi0  = 1.d0/dx/lUn0
  dyi0  = 1.d0/dy/lUn0
  dzi0  = 1.d0/dz/lUn0
  StBz  = 5.6704d-5         ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]

  ! iterate over the active domain
  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1

           !--------------
           ! z-direction, lower face
           E0d_zl = E0(i,j,k) - E0(i,j,k-1)
           Ed_zl  = E(i,j,k) - E(i,j,k-1)
           D0_zl  = limiter(E0(i,j,k), E0(i,j,k-1), kappa(i,j,k), &
                            kappa(i,j,k-1), nUn0, lUn0, tUn, dzi)
           D_zl   = limiter(E(i,j,k), E(i,j,k-1), kappa(i,j,k), &
                            kappa(i,j,k-1), nUn, lUn, tUn, dzi)

           !--------------
           ! y-direction, lower face
           E0d_yl = E0(i,j,k) - E0(i,j-1,k)
           Ed_yl  = E(i,j,k) - E(i,j-1,k)
           D0_yl  = limiter(E0(i,j,k), E0(i,j-1,k), kappa(i,j,k), &
                            kappa(i,j-1,k), nUn0, lUn0, tUn, dyi)
           D_yl   = limiter(E(i,j,k), E(i,j-1,k), kappa(i,j,k), &
                            kappa(i,j-1,k), nUn, lUn, tUn, dyi)

           !--------------
           ! x-direction, lower face
           E0d_xl = E0(i,j,k) - E0(i-1,j,k)
           Ed_xl  = E(i,j,k) - E(i-1,j,k)
           D0_xl  = limiter(E0(i,j,k), E0(i-1,j,k), kappa(i,j,k), &
                            kappa(i-1,j,k), nUn0, lUn0, tUn, dxi)
           D_xl   = limiter(E(i,j,k), E(i-1,j,k), kappa(i,j,k), &
                            kappa(i-1,j,k), nUn, lUn, tUn, dxi)

           !--------------
           ! x-direction, upper face
           E0d_xr = E0(i+1,j,k) - E0(i,j,k)
           Ed_xr  = E(i+1,j,k) - E(i,j,k)
           D0_xr  = limiter(E0(i,j,k), E0(i+1,j,k), kappa(i,j,k), &
                            kappa(i+1,j,k), nUn0, lUn0, tUn, dxi)
           D_xr   = limiter(E(i,j,k), E(i+1,j,k), kappa(i,j,k), &
                            kappa(i+1,j,k), nUn, lUn, tUn, dxi)

           !--------------
           ! y-direction, upper face
           E0d_yr = E0(i,j+1,k) - E0(i,j,k)
           Ed_yr  = E(i,j+1,k) - E(i,j,k)
           D0_yr  = limiter(E0(i,j,k), E0(i,j+1,k), kappa(i,j,k), &
                            kappa(i,j+1,k), nUn0, lUn0, tUn, dyi)
           D_yr   = limiter(E(i,j,k), E(i,j+1,k), kappa(i,j,k), &
                            kappa(i,j+1,k), nUn, lUn, tUn, dyi)

           !--------------
           ! z-direction, upper face
           E0d_zr = E0(i,j,k+1) - E0(i,j,k)
           Ed_zr  = E(i,j,k+1) - E(i,j,k)
           D0_zr  = limiter(E0(i,j,k), E0(i,j,k+1), kappa(i,j,k), &
                            kappa(i,j,k+1), nUn0, lUn0, tUn, dzi)
           D_zr   = limiter(E(i,j,k), E(i,j,k+1), kappa(i,j,k), &
                            kappa(i,j,k+1), nUn, lUn, tUn, dzi)

           ! opacity values in this cell
           kap  = kappa(i,j,k)*nUn
           kap0 = kappa(i,j,k)**nUn0

           ! black-body radiation in this cell (if applicable; otherwise Temp=0)
           eta  = 4.d0*kap*StBz/rUn*Temp(i,j,k)**4
           eta0 = 4.d0*kap0*StBz/rUn0*Temp0(i,j,k)**4

           ! set the matrix entries.  Note: the diffusive component 
           ! need not be rescaled, since scaling and chain rule cancel 
           matentries(1,i,j,k) = -dtfac*dzi*dzi*D_zl         ! z-left
           matentries(2,i,j,k) = -dtfac*dyi*dyi*D_yl         ! y-left
           matentries(3,i,j,k) = -dtfac*dxi*dxi*D_xl         ! x-left
           matentries(4,i,j,k) = 1.d0 &
                + dtfac*(afac + c_light*kap + dxi*dxi*(D_xl+D_xr)   &   ! self
                + dyi*dyi*(D_yl+D_yr) + dzi*dzi*(D_zl+D_zr))
           matentries(5,i,j,k) = -dtfac*dxi*dxi*D_xr         ! x-right
           matentries(6,i,j,k) = -dtfac*dyi*dyi*D_yr         ! y-right
           matentries(7,i,j,k) = -dtfac*dzi*dzi*D_zr         ! z-right

           ! set the rhs entries
           rhsentries(i,j,k) = ( (dtfac/rUn + dtfac0/rUn0)*src(i,j,k)           &
                               + dtfac*eta + dtfac0*eta0                        &
                               + (1.d0 - dtfac0*(afac0+c_light*kap0))*E0(i,j,k) &
                               + dtfac0*dxi0*dxi0*(D0_xr*E0d_xr-D0_xl*E0d_xl)   &
                               + dtfac0*dyi0*dyi0*(D0_yr*E0d_yr-D0_yl*E0d_yl)   &
                               + dtfac0*dzi0*dzi0*(D0_zr*E0d_zr-D0_zl*E0d_zl)   &
                               - (1.d0 + dtfac*(afac+c_light*kap))*E(i,j,k)     &
                               + dtfac*dxi*dxi*(D_xr*Ed_xr-D_xl*Ed_xl)          &
                               + dtfac*dyi*dyi*(D_yr*Ed_yr-D_yl*Ed_yl)          &
                               + dtfac*dzi*dzi*(D_zr*Ed_zr-D_zl*Ed_zl) )

!!$           if (i==126 .and. j==62 .and. k==64) then
!!$              print *,'Evals =',E(i,j,k-1),E(i,j-1,k),E(i-1,j,k),E(i,j,k),E(i+1,j,k),E(i,j+1,k),E(i,j,k+1)
!!$              print *,'kappas =',kappa(i,j,k-1),kappa(i,j-1,k),kappa(i-1,j,k),kappa(i,j,k),kappa(i+1,j,k),kappa(i,j+1,k),kappa(i,j,k+1)
!!$              print *,'nUn =',nUn,',  nUn0 =',nUn0,',  rUn =',rUn,',  rUn0 =',rUn0
!!$              print *,'dxi =',dxi,',  dyi =',dyi,',  dzi =',dzi
!!$              print *,'dtfac =',dtfac,',  dtfac0 =',dtfac0
!!$              print *,'eta =',eta,',  eta0 =',eta0
!!$              print *,'src =',src(i,j,k)
!!$              print *,'D_* =',D_zl,D_yl,D_xl,D_xr,D_yr,D_zr
!!$           end if

        enddo
     enddo
  enddo


  ! update matrix/rhs based on boundary conditions/location
  !    z-left face
  if (zlface == 1) then
     ! Dirichlet
     if (BCZl==1) then
        k = x2s
        do j=1,Ny
           do i=1,Nx
              matentries(1,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCZl==2) then
        k = x2s
        do j=1,Ny
           do i=1,Nx
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(1,i,j,k)
              matentries(1,i,j,k) = 0.d0
           enddo
        enddo
     endif
  end if

  !    y-left face
  if (ylface == 1) then
     ! Dirichlet
     if (BCYl==1) then
        j = x1s
        do k=1,Nz
           do i=1,Nx
              matentries(2,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCYl==2) then
        j = x1s
        do k=1,Nz
           do i=1,Nx
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(2,i,j,k)
              matentries(2,i,j,k) = 0.d0
           enddo
        enddo
     endif
  end if

  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        do k=1,Nz
           do j=1,Ny
              matentries(3,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCXl==2) then
        i = x0s
        do k=1,Nz
           do j=1,Ny
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(3,i,j,k)
              matentries(3,i,j,k) = 0.d0
           enddo
        enddo
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        do k=1,Nz
           do j=1,Ny
              matentries(5,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCXr==2) then
        i = x0e
        do k=1,Nz
           do j=1,Ny
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(5,i,j,k)
              matentries(5,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  !    y-right face
  if (yrface==1) then
     ! Dirichlet
     if (BCYr==1) then
        j = x1e
        do k=1,Nz
           do i=1,Nx
              matentries(6,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCYr==2) then
        j = x1e
        do k=1,Nz
           do i=1,Nx
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(6,i,j,k)
              matentries(6,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  !    z-right face
  if (zrface==1) then
     ! Dirichlet
     if (BCZr==1) then
        k = x2e
        do j=1,Ny
           do i=1,Nx
              matentries(7,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCZr==2) then
        k = x2e
        do j=1,Ny
           do i=1,Nx
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(7,i,j,k)
              matentries(7,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  rhsnorm = sum(rhsentries*rhsentries)

  return

end subroutine gFLDSplit_SetupSystem3D
!=======================================================================






subroutine gFLDSplit_SetupSystem2D(matentries, rhsentries, rhsnorm, E0,   &
     E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,       &
     theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn, dx, dy, BCXl, BCXr,     &
     BCYl, BCYr, x0s, x0e, x1s, x1e, Nx, Ny, NGxl, NGxr, NGyl, NGyr,      &
     xlface, xrface, ylface, yrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       July 2009
!  modified:   
!
!  PURPOSE: 2D version of the routine
!=======================================================================
  implicit none
  
  !--------------
  ! argument declarations
  integer,  intent(in) :: ESpectrum
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer,  intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  REALSUB,  intent(in) :: a, a0, adot, adot0
  real,     intent(in) :: dt, theta, dx, dy
  real,     intent(in) :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), intent(in) &
       :: E0, E, src, kappa, Temp, Temp0
  real*8,   intent(out) :: matentries(5,x0s:x0e,x1s:x1e)
  real*8,   intent(out) :: rhsentries(x0s:x0e,x1s:x1e)
  real,     intent(out) :: rhsnorm
  integer,  intent(out) :: ier

  !--------------
  ! locals
  integer :: i, j
  real*8  :: dtfac, dtfac0, kap, kap0, StBz, eta, eta0
  real*8  :: dxi, dxi0, dyi, dyi0, afac, afac0
  real*8  :: D_xl, D0_xl, D_xr, D0_xr, E0d_xl, E0d_xr, Ed_xl, Ed_xr
  real*8  :: D_yl, D0_yl, D_yr, D0_yr, E0d_yl, E0d_yr, Ed_yl, Ed_yr
  real*8, external :: limiter

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  matentries = 0.d0
  rhsentries = 0.d0
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta         ! time step conversion factor
  dtfac0 = dt*(1.d0-theta)  ! time step conversion factor
  if (ESpectrum == -1) then
     afac  = 0.d0
     afac0 = 0.d0
  else
     afac  = adot/a         ! expansion factor (new time)
     afac0 = adot0/a0       ! expansion factor (old time)
  endif
  dxi   = 1.d0/dx/lUn
  dyi   = 1.d0/dy/lUn
  dxi0  = 1.d0/dx/lUn0
  dyi0  = 1.d0/dy/lUn0
  StBz  = 5.6704d-5         ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]

  ! iterate over the active domain
  do j=1,Ny,1
     do i=1,Nx,1

        !--------------
        ! y-direction, lower face
        E0d_yl = E0(i,j) - E0(i,j-1)
        Ed_yl  = E(i,j) - E(i,j-1)
        D0_yl  = limiter(E0(i,j), E0(i,j-1), kappa(i,j), &
                         kappa(i,j-1), nUn0, lUn0, tUn, dyi)
        D_yl   = limiter(E(i,j), E(i,j-1), kappa(i,j), &
                         kappa(i,j-1), nUn, lUn, tUn, dyi)

        !--------------
        ! x-direction, lower face
        E0d_xl = E0(i,j) - E0(i-1,j)
        Ed_xl  = E(i,j) - E(i-1,j)
        D0_xl  = limiter(E0(i,j), E0(i-1,j), kappa(i,j), &
                         kappa(i-1,j), nUn0, lUn0, tUn, dxi)
        D_xl   = limiter(E(i,j), E(i-1,j), kappa(i,j), &
                         kappa(i-1,j), nUn, lUn, tUn, dxi)

        !--------------
        ! x-direction, upper face
        E0d_xr = E0(i+1,j) - E0(i,j)
        Ed_xr  = E(i+1,j) - E(i,j)
        D0_xr  = limiter(E0(i,j), E0(i+1,j), kappa(i,j), &
                         kappa(i+1,j), nUn0, lUn0, tUn, dxi)
        D_xr   = limiter(E(i,j), E(i+1,j), kappa(i,j), &
                         kappa(i+1,j), nUn, lUn, tUn, dxi)

        !--------------
        ! y-direction, upper face
        E0d_yr = E0(i,j+1) - E0(i,j)
        Ed_yr  = E(i,j+1) - E(i,j)
        D0_yr  = limiter(E0(i,j), E0(i,j+1), kappa(i,j), &
                         kappa(i,j+1), nUn0, lUn0, tUn, dyi)
        D_yr   = limiter(E(i,j), E(i,j+1), kappa(i,j), &
                         kappa(i,j+1), nUn, lUn, tUn, dyi)

        ! opacity values in this cell
        kap  = kappa(i,j)*nUn
        kap0 = kappa(i,j)*nUn0

        ! black-body radiation in this cell (if applicable)
        eta  = 4.d0*kap*StBz/rUn*Temp(i,j)**4
        eta0 = 4.d0*kap0*StBz/rUn0*Temp0(i,j)**4

        ! set the matrix entries.  Note: the diffusive component 
        ! need not be rescaled, since scaling and chain rule cancel 
        matentries(1,i,j) = -dtfac*dyi*dyi*D_yl         ! y-left
        matentries(2,i,j) = -dtfac*dxi*dxi*D_xl         ! x-left
        matentries(3,i,j) = 1.d0 + dtfac*(afac + c_light*kap  &  ! self
              + dxi*dxi*(D_xl+D_xr)+dyi*dyi*(D_yl+D_yr))
        matentries(4,i,j) = -dtfac*dxi*dxi*D_xr         ! x-right
        matentries(5,i,j) = -dtfac*dyi*dyi*D_yr         ! y-right

        ! set the rhs entries
        rhsentries(i,j) = ( (dtfac/rUn + dtfac0/rUn0)*src(i,j)            &
                          + dtfac*eta + dtfac0*eta0                       &
                          + (1.d0 - dtfac0*(afac0+c_light*kap0))*E0(i,j)  &
                          + dtfac0*dxi0*dxi0*(D0_xr*E0d_xr-D0_xl*E0d_xl)  &
                          + dtfac0*dyi0*dyi0*(D0_yr*E0d_yr-D0_yl*E0d_yl)  &
                          - (1.d0 + dtfac*(afac+c_light*kap))*E(i,j)      &
                          + dtfac*dxi*dxi*(D_xr*Ed_xr-D_xl*Ed_xl)         &
                          + dtfac*dyi*dyi*(D_yr*Ed_yr-D_yl*Ed_yl) )
     enddo
  enddo
  
  ! update matrix/rhs based on boundary conditions/location
  !    y-left face
  if (ylface == 1) then
     ! Dirichlet
     if (BCYl==1) then
        j = x1s
        do i=1,Nx
           matentries(1,i,j) = 0.d0
        enddo
        ! Neumann
     else if (BCYl==2) then
        j = x1s
        do i=1,Nx
           matentries(3,i,j) = matentries(3,i,j) + matentries(1,i,j)
           matentries(1,i,j) = 0.d0
        enddo
     endif
  end if

  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        do j=1,Ny
           matentries(2,i,j) = 0.d0
        enddo
        ! Neumann
     else if (BCXl==2) then
        i = x0s
        do j=1,Ny
           matentries(3,i,j) = matentries(3,i,j) + matentries(2,i,j)
           matentries(2,i,j) = 0.d0
        enddo
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        do j=1,Ny
           matentries(4,i,j) = 0.d0
        enddo
        ! Neumann
     else if (BCXr==2) then
        i = x0e
        do j=1,Ny
           matentries(3,i,j) = matentries(3,i,j) + matentries(4,i,j)
           matentries(4,i,j) = 0.d0
        enddo
     endif
  endif

  !    y-right face
  if (yrface==1) then
     ! Dirichlet
     if (BCYr==1) then
        j = x1e
        do i=1,Nx
           matentries(5,i,j) = 0.d0
        enddo
        ! Neumann
     else if (BCYr==2) then
        j = x1e
        do i=1,Nx
           matentries(3,i,j) = matentries(3,i,j) + matentries(5,i,j)
           matentries(5,i,j) = 0.d0
        enddo
     endif
  endif

  rhsnorm = sum(rhsentries*rhsentries)

  return
end subroutine gFLDSplit_SetupSystem2D
!=======================================================================






subroutine gFLDSplit_SetupSystem1D(matentries, rhsentries, rhsnorm, E0, &
     E, Temp, Temp0, kappa, src, dt, a, a0, adot, adot0, ESpectrum,     &
     theta, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn, dx, BCXl, BCXr, x0s,  &
     x0e, Nx, NGxl, NGxr, xlface, xrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       June 2009
!  modified:   
!
!  PURPOSE: 1D version of the routine
!=======================================================================
  implicit none
  
  !--------------
  ! argument declarations
  integer,  intent(in) :: ESpectrum
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  REALSUB,  intent(in) :: a, a0, adot, adot0
  real,     intent(in) :: dt, theta, dx
  real,     intent(in) :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, tUn
  real, dimension(1-NGxl:Nx+NGxr), intent(in) :: E0, E, src, kappa, Temp, Temp0
  real*8,   intent(out) :: matentries(3,x0s:x0e)
  real*8,   intent(out) :: rhsentries(x0s:x0e)
  real,     intent(out) :: rhsnorm
  integer,  intent(out) :: ier

  !--------------
  ! locals
  integer :: i
  real*8  :: dtfac, dtfac0, kap, kap0, StBz, eta, eta0
  real*8  :: dxi, dxi0, afac, afac0
  real*8  :: D_xl, D0_xl, D_xr, D0_xr, E0d_xl, E0d_xr, Ed_xl, Ed_xr
  real*8, external :: limiter

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  matentries = 0.d0
  rhsentries = 0.d0
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta         ! time step conversion factor
  dtfac0 = dt*(1.d0-theta)  ! time step conversion factor
  if (ESpectrum == -1) then
     afac  = 0.d0
     afac0 = 0.d0
  else
     afac  = adot/a         ! expansion factor (new time)
     afac0 = adot0/a0       ! expansion factor (old time)
  endif
  dxi   = 1.d0/dx/lUn
  dxi0  = 1.d0/dx/lUn0
  StBz  = 5.6704d-5         ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]

  ! iterate over the active domain
  do i=1,Nx,1

     !--------------
     ! x-direction, lower face
     E0d_xl = E0(i) - E0(i-1)
     Ed_xl  = E(i) - E(i-1)
     D0_xl  = limiter(E0(i), E0(i-1), kappa(i), kappa(i-1), nUn0, lUn0, tUn, dxi)
     D_xl   = limiter(E(i), E(i-1), kappa(i), kappa(i-1), nUn, lUn, tUn, dxi)

     !--------------
     ! x-direction, upper face
     E0d_xr = E0(i+1) - E0(i)
     Ed_xr  = E(i+1) - E(i)
     D0_xr  = limiter(E0(i), E0(i+1), kappa(i), kappa(i+1), nUn0, lUn0, tUn, dxi)
     D_xr   = limiter(E(i), E(i+1), kappa(i), kappa(i+1), nUn, lUn, tUn, dxi)

     ! opacity values in this cell
     kap  = kappa(i)*nUn
     kap0 = kappa(i)*nUn0

     ! black-body radiation in this cell (if applicable)
     eta  = 4.d0*kap*StBz/rUn*Temp(i)**4
     eta0 = 4.d0*kap0*StBz/rUn0*Temp0(i)**4

     ! set the matrix entries.  Note: the diffusive component 
     ! need not be rescaled, since scaling and chain rule cancel 
     matentries(1,i) = -dtfac*dxi*dxi*D_xl                  ! x-left
     matentries(2,i) = 1.d0 + dtfac*(afac + c_light*kap  &  ! self
                            + dxi*dxi*(D_xl+D_xr))
     matentries(3,i) = -dtfac*dxi*dxi*D_xr                  ! x-right

     ! set the rhs entries
     rhsentries(i) = ( (dtfac/rUn + dtfac0/rUn0)*src(i)              &
                     + dtfac*eta + dtfac0*eta0                       &
                     + (1.d0 - dtfac0*(afac0+c_light*kap0))*E0(i)    &
                     + dtfac0*dxi0*dxi0*(D0_xr*E0d_xr-D0_xl*E0d_xl)  &
                     - (1.d0 + dtfac*(afac+c_light*kap))*E(i)        &
                     + dtfac*dxi*dxi*(D_xr*Ed_xr-D_xl*Ed_xl) )

  enddo

  ! update matrix/rhs based on boundary conditions/location
  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        matentries(1,i) = 0.d0
     ! Neumann
     else if (BCXl==2) then
        i = x0s
        matentries(2,i) = matentries(2,i) + matentries(1,i)
        matentries(1,i) = 0.d0
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        matentries(3,i) = 0.d0
     ! Neumann
     else if (BCXr==2) then
        i = x0e
        matentries(2,i) = matentries(2,i) + matentries(3,i)
        matentries(3,i) = 0.d0
     endif
  endif

  rhsnorm = sum(rhsentries*rhsentries)

  return
end subroutine gFLDSplit_SetupSystem1D
!=======================================================================
