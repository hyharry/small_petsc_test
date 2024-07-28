!
!  This example shows how to avoid Fortran line lengths larger than 132 characters.
!  It avoids used of certain macros such as PetscCallA() and PetscCheckA() that
!  generate very long lines
!
!  We recommend starting from src/snes/tutorials/ex9f90.F90 instead of this example
!  because that does not have the restricted formatting that makes this version
!  more difficult to read
!
!  Description: This example solves a nonlinear system in parallel with SNES.
!  We solve the  Bratu (SFI - solid fuel ignition) problem in a 2D rectangular
!  domain, using distributed arrays (DMDAs) to partition the parallel grid.
!  The command line options include:
!    -par <param>, where <param> indicates the nonlinearity of the problem
!       problem SFI:  <parameter> = Bratu parameter (0 <= par <= 6.81)
!
!  --------------------------------------------------------------------------
!
!  Solid Fuel Ignition (SFI) problem.  This problem is modeled by
!  the partial differential equation
!
!          -Laplacian u - lambda*exp(u) = 0,  0 < x,y < 1,
!
!  with boundary conditions
!
!           u = 0  for  x = 0, x = 1, y = 0, y = 1.
!
!  A finite difference approximation with the usual 5-point stencil
!  is used to discretize the boundary value problem to obtain a nonlinear
!  system of equations.
!
!  --------------------------------------------------------------------------
      module ex9fmodule
      use petscsnes
      use petscdmda
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
      PetscInt xs,xe,xm,gxs,gxe,gxm
      PetscInt ys,ye,ym,gys,gye,gym
      PetscInt mx,my
      PetscMPIInt rank,size
      PetscReal lambda

      contains

      function psi(x, y) result(retval)
            PetscReal :: x, y
            PetscReal :: retval
            PetscReal, parameter :: r0 = 0.9

            PetscReal, parameter :: psi0 = sqrt(1.0 - r0**2)
            PetscReal, parameter :: dpsi0 = -r0 / psi0

            PetscReal :: r

            r = x**2 + y**2

            if ( r <= r0 ) then
                  retval = 1.0 - r
            else 
                  retval = psi0 + dpsi0 * (r - r0)
            end if
            
      end function psi

      function u_exact_f(x, y) result(retval)
            PetscReal :: x, y
            PetscReal :: retval
            PetscReal :: afree = 0.697965148223374, A = 0.680259411891719, B = 0.471519893402112
            PetscReal :: r 

            r = sqrt(x**2 + y**2)
            if ( r <= afree ) then
                  retval = psi(x, y)
            else 
                  retval = -A * log(r) + B
            end if

      end function u_exact_f

      end module ex9fmodule

      program main
      use ex9fmodule
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                   Variable declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     snes        - nonlinear solver
!     x, r        - solution, residual vectors
!     its         - iterations for convergence
!
!  See additional variable declarations in the file ex9f.h
!
      SNES           snes
      DM             da, da_after
      Vec            u, u_exact
      DMDALocalInfo :: info(DMDA_LOCAL_INFO_SIZE)
      PetscReal      error1, errorinf
      PetscErrorCode ierr

      PetscReal two,zero,one

      character(len=50) :: outputString

!  Note: Any user-defined Fortran routines (such as FormJacobianLocal)
!  MUST be declared as external.

      ! Yi: todo change it!
      external FormBounds
      external FormExactSolution
      external FormFunctionLocal,FormJacobianLocal

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Initialize program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      two = 2.0
      zero = 0.0
      one = 1.0

      call PetscInitialize(ierr)
      CHKERRA(ierr)
      call DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,5,5, &
                        PETSC_DECIDE,PETSC_DECIDE,1,1, &
                        PETSC_NULL_INTEGER_ARRAY,PETSC_NULL_INTEGER_ARRAY,da,ierr)
      CHKERRA(ierr)
      call DMSetFromOptions(da,ierr)
      CHKERRA(ierr)
      call DMSetUp(da,ierr)
      CHKERRA(ierr)
      call DMDASetUniformCoordinates(da, -two, two, -two, two, zero, one,ierr)
      CHKERRA(ierr)

      call DMCreateGlobalVector(da,u,ierr)
      CHKERRA(ierr)
      call VecSet(u, zero, ierr)
      CHKERRA(ierr)

      call SNESCreate(PETSC_COMM_WORLD,snes,ierr)
      CHKERRA(ierr)
      call SNESSetDM(snes,da,ierr)
      CHKERRA(ierr)
      call SNESSetType(snes, SNESVINEWTONRSLS, ierr)
      CHKERRA(ierr)
      call SNESVISetComputeVariableBounds(snes, FormBounds, ierr)
      CHKERRA(ierr)
      call DMDASNESSetFunctionLocal(da, INSERT_VALUES, FormFunctionLocal, PETSC_NULL_SNES, ierr)
      CHKERRA(ierr)
      call DMDASNESSetJacobianLocal(da, FormJacobianLocal, PETSC_NULL_SNES, ierr)
      CHKERRA(ierr)
      call SNESSetFromOptions(snes, ierr)
      CHKERRA(ierr)

      call DMDAGetInfo(da,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_ENUM,PETSC_NULL_ENUM, &
                       PETSC_NULL_ENUM,PETSC_NULL_ENUM,ierr)
      CHKERRA(ierr)
      call DMDAGetCorners(da,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)
      CHKERRA(ierr)
      call DMDAGetGhostCorners(da,gxs,gys,PETSC_NULL_INTEGER,gxm,gym,PETSC_NULL_INTEGER,ierr)
      CHKERRA(ierr)

      xs  = xs+1
      ys  = ys+1
      gxs = gxs+1
      gys = gys+1

      ye  = ys+ym-1
      xe  = xs+xm-1
      gye = gys+gym-1
      gxe = gxs+gxm-1

      call SNESSolve(snes,PETSC_NULL_VEC,u,ierr)
      CHKERRA(ierr)
      call VecDestroy(u,ierr)
      CHKERRA(ierr)
      call DMDestroy(da,ierr)
      CHKERRA(ierr)

      call SNESGetDM(snes, da_after, ierr)
      CHKERRA(ierr)
      call SNESGetSolution(snes, u, ierr) ! do not destroy u 
      CHKERRA(ierr)
      call DMDAGetLocalInfo(da_after, info, ierr)
      CHKERRA(ierr)
      call VecDuplicate(u, u_exact, ierr)
      CHKERRA(ierr)
      call FormExactSolution(info,da_after, u_exact, ierr) ! TODO: which da, da_after?
      CHKERRA(ierr)
      call VecAXPY(u, -one, u_exact, ierr) ! u <-- u - u_exact
      CHKERRA(ierr)
      call VecNorm(u, NORM_1, error1, ierr)
      CHKERRA(ierr)
      error1 = error1 / mx * my
      call VecNorm(u, NORM_INFINITY, errorinf, ierr)
      CHKERRA(ierr)
      write(outputString,'(A,I0,A,I0,A,A,1PE10.3,A,A,1PE10.3)') 'errors on ', mx, ' x ', my, &
       ' grid:  av |u-uexact|  = ', error1, ',  |u-uexact|_inf = ', errorinf
      call PetscPrintf(PETSC_COMM_WORLD,outputString, ierr)
      CHKERRA(ierr)
      call VecDestroy(u_exact, ierr)
      CHKERRA(ierr)
      call SNESDestroy(snes, ierr)
      CHKERRA(ierr)
      call DMDestroy(da, ierr)
      CHKERRA(ierr)
      call PetscFinalize(ierr)
      CHKERRA(ierr)



      
      end

!================ Yi add ==============
!====================================== 
      subroutine FormExactSolution(info,da,u,ierr)
            use ex9fmodule
            implicit none
            
            DMDALocalInfo info(DMDA_LOCAL_INFO_SIZE)
            Vec :: u
            DM  :: da

            PetscInt  :: i, j
            PetscReal :: dx, dy, x, y
            PetscScalar, pointer :: au(:,:)
            PetscErrorCode :: ierr

            xs     = info(DMDA_LOCAL_INFO_XS)+1
            xe     = xs+info(DMDA_LOCAL_INFO_XM)-1
            ys     = info(DMDA_LOCAL_INFO_YS)+1
            ye     = ys+info(DMDA_LOCAL_INFO_YM)-1
            mx     = info(DMDA_LOCAL_INFO_MX)
            my     = info(DMDA_LOCAL_INFO_MY)
            
            dx = 4.0 / (mx - 1)
            dy = 4.0 / (my - 1)
            call DMDAVecGetArrayF90(da, u, au, ierr)
            CHKERRQ(ierr)
            do j = ys,ye
              y = -2.0 + j * dy
              do i = xs,xe
                x = -2.0 + i * dx
                au(i,j) = u_exact_f(x, y)
              end do
            end do
            call DMDAVecRestoreArrayF90(da, u, au, ierr)
            CHKERRQ(ierr)
      
            
      end subroutine FormExactSolution

      subroutine FormBounds(snes,Xl,Xu,ierr)
            use ex9fmodule
            implicit none
            
            SNES :: snes
            Vec :: Xl, Xu
            DM :: da
            PetscErrorCode :: ierr

            PetscInt  :: i, j
            PetscReal :: dx, dy, x, y
            PetscScalar, pointer :: aXl(:,:)

            DMDALocalInfo :: info(DMDA_LOCAL_INFO_SIZE)

            call SNESGetDM(snes, da, ierr)
            CHKERRQ(ierr)
            call DMDAGetLocalInfo(da, info, ierr)
            CHKERRQ(ierr)

            xs     = info(DMDA_LOCAL_INFO_XS)+1
            xe     = xs+info(DMDA_LOCAL_INFO_XM)-1
            ys     = info(DMDA_LOCAL_INFO_YS)+1
            ye     = ys+info(DMDA_LOCAL_INFO_YM)-1
            mx     = info(DMDA_LOCAL_INFO_MX)
            my     = info(DMDA_LOCAL_INFO_MY)
            
            dx = 4.0 / (mx - 1)
            dy = 4.0 / (my - 1)

            call DMDAVecGetArrayF90(da, Xl, aXl, ierr)
            CHKERRQ(ierr)
            do j = ys,ye
              y = -2.0 + j * dy
              do i = xs,xe
                x = -2.0 + i * dx
                aXl(i,j) = psi(x, y)
              end do
            end do
            call DMDAVecRestoreArrayF90(da, Xl, aXl, ierr)
            CHKERRQ(ierr)
            call VecSet(Xu, PETSC_INFINITY, ierr)
            CHKERRQ(ierr)
            
      end subroutine FormBounds

      subroutine FormFunctionLocal(info, au, af, user, ierr)
            use ex9fmodule
            implicit none

            DMDALocalInfo :: info(DMDA_LOCAL_INFO_SIZE)
            PetscScalar :: au(:,:), af(:,:) ! Yi: how to get size of au and af??
            PetscInt :: user 
            PetscErrorCode :: ierr

            PetscReal :: twelve

            PetscInt  i, j
            PetscReal dx, dy, x, y, ue, un, us, uw

            twelve = 12.0
            ! Yi: this piece of info should be embedded in a module and func?
            xs     = info(DMDA_LOCAL_INFO_XS)+1
            xe     = xs+info(DMDA_LOCAL_INFO_XM)-1
            ys     = info(DMDA_LOCAL_INFO_YS)+1
            ye     = ys+info(DMDA_LOCAL_INFO_YM)-1
            mx     = info(DMDA_LOCAL_INFO_MX)
            my     = info(DMDA_LOCAL_INFO_MY)
            
            dx = 4.0 / (mx - 1)
            dy = 4.0 / (my - 1)

            do j = ys,ye
              y = -2.0 + j * dy
              do i = xs,xe
                x = -2.0 + i * dx
                if ( i==1 .or. j==1 .or. i==mx .or. j==my ) then
                  af(i,j) = 4.0 * (au(i,j) - u_exact_f(x, y))
                else
                  if ( i - 1 == 1 ) then
                        uw = u_exact_f(x-dx, y)
                  else
                        uw = au(i-1,j)
                  end if
                  if ( i + 1 == mx ) then
                        ue = u_exact_f(x+dx, y)
                  else
                        ue = au(i+1,j)
                  end if
                  if ( j - 1 == 1 ) then
                        us = u_exact_f(x,y-dy)
                  else
                        us = au(i,j-1)
                  end if
                  if ( j + 1 == my ) then
                        un = u_exact_f(x,y+dy)
                  else
                        un = au(i,j+1)
                  end if
                  af(i,j) = -(dy / dx) * (uw - 2.0 * au(i,j) + ue) - (dx / dy) * (us - 2.0 * au(i,j) + un)
                end if
              end do
            end do
            call PetscLogFlops(twelve * ym * xm, ierr)
            CHKERRQ(ierr)
            
      end subroutine FormFunctionLocal

      subroutine FormJacobianLocal(info, au, A, jac, user, ierr)
            use ex9fmodule
            implicit none

            DMDALocalInfo :: info(DMDA_LOCAL_INFO_SIZE)
            PetscScalar :: au(:,:) ! Yi: how to get size of au and af??
            Mat :: A, jac
            PetscInt :: user 
            PetscErrorCode :: ierr

            PetscReal :: two

            PetscInt   i, j, n
            MatStencil col(4,5), row(4)
            PetscReal  v(5), dx, dy, oxx, oyy
            
            two = 2.0
            ! Yi: this piece of info should be embedded in a module and func?
            xs     = info(DMDA_LOCAL_INFO_XS)+1
            xe     = xs+info(DMDA_LOCAL_INFO_XM)-1
            ys     = info(DMDA_LOCAL_INFO_YS)+1
            ye     = ys+info(DMDA_LOCAL_INFO_YM)-1
            mx     = info(DMDA_LOCAL_INFO_MX)
            my     = info(DMDA_LOCAL_INFO_MY)
            
            dx = 4.0 / (mx - 1)
            dy = 4.0 / (my - 1)
            oxx = dy / dx
            oyy = dx / dy

            do j = ys,ye
              do i = xs,xe
                  row(MatStencil_i) = i
                  row(MatStencil_j) = j
                  if ( i==1 .or. j==1 .or. i==mx .or. j==my ) then
                        v(1) = 4.0
                        call MatSetValuesStencil(jac,1,row,1,row,v,INSERT_VALUES,ierr)
                        CHKERRQ(ierr)
                  else
                        v(1) = 2.0 * (oxx + oyy)
                        col(MatStencil_i,1) = i
                        col(MatStencil_j,1) = j
                        n = 1
                        if ( i-1 > 1 ) then
                              v(2) = -oxx
                              col(MatStencil_i,2) = i-1
                              col(MatStencil_j,2) = j
                              n = n+1
                        end if
                        if ( i+1 < mx ) then
                              v(3) = -oxx
                              col(MatStencil_i,3) = i+1
                              col(MatStencil_j,3) = j
                              n = n+1
                        end if
                        if ( j-1 > 1 ) then
                              v(4) = -oyy
                              col(MatStencil_i,4) = i
                              col(MatStencil_j,4) = j-1
                              n = n+1
                        end if
                        if ( j+1 < my ) then
                              v(5) = -oyy
                              col(MatStencil_i,5) = i
                              col(MatStencil_j,5) = j+1
                              n = n+1
                        end if
                        call MatSetValuesStencil(jac,1,row,n,col,v,INSERT_VALUES,ierr)
                        CHKERRQ(ierr)
                  end if
              end do
            end do

            call MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY, ierr)
            CHKERRQ(ierr)
            call MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY, ierr)
            CHKERRQ(ierr)
            if (A /= jac) then
              call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
              CHKERRQ(ierr)
              call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
              CHKERRQ(ierr)
             end if
            call PetscLogFlops(two* ym * xm, ierr)
            CHKERRQ(ierr)
            
      end subroutine FormJacobianLocal

!================ Yi end ==============
!======================================
