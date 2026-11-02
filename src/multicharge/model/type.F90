! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> @file multicharge/model/type.f90
!> Provides a general base class for the charge models

#ifndef IK
#define IK i4
#endif

!> General charge model
module multicharge_model_type

   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_cutoff, only: get_lattice_points
   use mctc_env, only: error_type, fatal_error, wp, ik => IK
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matinv_3x3
   use mctc_ncoord, only: ncoord_type
   use mctc_wignerseitz, only: wignerseitz_cell
   use multicharge_blas, only: gemv, symv, gemm
   use multicharge_ewald, only: reduce_basis_2d
   use multicharge_lapack, only: sytrf, sytrs, sytri
   use multicharge_model_cache, only: model_cache, cache_container
   implicit none
   private

   public :: mchrg_model_type, get_dir_trans, get_rec_trans, check_structure
   public :: init_wignerseitz_cell, get_periodic_axes, cross_product_3d

   !> Abstract multicharge model type
   type, abstract :: mchrg_model_type
      !> Electronegativity
      real(wp), allocatable :: chi(:)
      !> Charge width
      real(wp), allocatable :: rad(:)
      !> Chemical hardness
      real(wp), allocatable :: eta(:)
      !> CN scaling factor for electronegativity
      real(wp), allocatable :: kcnchi(:)
      !> Local charge scaling factor for electronegativity
      real(wp), allocatable :: kqchi(:)
      !> Local charge scaling factor for chemical hardness
      real(wp), allocatable :: kqeta(:)
      !> CN scaling factor for charge width
      real(wp), allocatable :: kcnrad
      !> Coordination number
      class(ncoord_type), allocatable :: ncoord
      !> Electronegativity weighted CN for local charge
      class(ncoord_type), allocatable :: ncoord_en
   contains
      !> Solve linear equations for the charge model
      procedure :: solve
      !> Contract the second charge response with fixed atomic weights
      procedure :: add_charge_hessian
      !> Model-specific curvature u^T (x_ab - A_ab q)
      procedure :: add_model_hessian
      !> Calculate local charges from electronegativity weighted CN
      procedure :: local_charge
      !> Update cache
      procedure(update), deferred :: update
      !> Calculate right-hand side (electronegativity)
      procedure(get_xvec), deferred :: get_xvec
      !> Calculate xvec Gradients
      procedure(get_xvec_derivs), deferred :: get_xvec_derivs
      !> Calculate Coulomb matrix
      procedure(get_coulomb_matrix), deferred :: get_coulomb_matrix
      !> Calculate Coulomb matrix derivatives
      procedure(get_coulomb_derivs), deferred :: get_coulomb_derivs
   end type mchrg_model_type

   abstract interface
      subroutine update(self, mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(in) :: cn(:)
         real(wp), intent(in), optional :: qloc(:)
         real(wp), intent(in), optional :: dcndr(:, :, :)
         real(wp), intent(in), optional :: dcndL(:, :, :)
         real(wp), intent(in), optional :: dqlocdr(:, :, :)
         real(wp), intent(in), optional :: dqlocdL(:, :, :)
      end subroutine update

      subroutine get_coulomb_matrix(self, mol, cache, amat)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(out) :: amat(:, :)
      end subroutine get_coulomb_matrix

      subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(in) :: qvec(:)
         real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
      end subroutine get_coulomb_derivs

      subroutine get_xvec(self, mol, cache, xvec)
         import :: mchrg_model_type, cache_container, structure_type, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(out) :: xvec(:)
      end subroutine get_xvec

      subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(out), contiguous :: dxdr(:, :, :)
         real(wp), intent(out), contiguous :: dxdL(:, :, :)
      end subroutine get_xvec_derivs
   end interface

   real(wp), parameter :: twopi = 2 * pi

contains

!> Compute 3D cross product of two vectors
pure function cross_product_3d(a, b) result(c)
   !> First vector
   real(wp), intent(in) :: a(3)
   !> Second vector
   real(wp), intent(in) :: b(3)
   !> Cross product a x b
   real(wp) :: c(3)

   c(1) = a(2) * b(3) - a(3) * b(2)
   c(2) = a(3) * b(1) - a(1) * b(3)
   c(3) = a(1) * b(2) - a(2) * b(1)
end function cross_product_3d

!> Extract the number and indices of active periodic axes
pure subroutine get_periodic_axes(periodic, ndim, p_axes)
   !> Periodicity mask
   logical, intent(in) :: periodic(:)
   !> Number of periodic dimensions (0, 1, 2, or 3)
   integer, intent(out) :: ndim
   !> Indices of active periodic axes in 1..ndim
   integer, intent(out) :: p_axes(3)

   integer :: i

   ndim = 0
   p_axes = 0
   do i = 1, min(3, size(periodic))
      if (periodic(i)) then
         ndim = ndim + 1
         p_axes(ndim) = i
      end if
   end do
end subroutine get_periodic_axes

!> Validate periodic subspace geometry for a structure
subroutine check_structure(mol, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ndim, p_axes(3)
   real(wp) :: a1(3), a2(3), a3(3), len1, len2, area, vol
   real(wp), parameter :: tol_deg = 1.0e-10_wp

   if (.not. allocated(mol%periodic)) return
   call get_periodic_axes(mol%periodic, ndim, p_axes)

   select case (ndim)
   case default
      return
   case (1)
      a1 = mol%lattice(:, p_axes(1))
      if (.not. all(ieee_is_finite(a1))) then
         call fatal_error(error, "1D periodic lattice vector contains non-finite values")
         return
      end if
      len1 = norm2(a1)
      if (len1 <= tol_deg) then
         call fatal_error(error, "Degenerate 1D periodic lattice vector")
         return
      end if
   case (2)
      a1 = mol%lattice(:, p_axes(1))
      a2 = mol%lattice(:, p_axes(2))
      if (.not. all(ieee_is_finite(a1)) .or. .not. all(ieee_is_finite(a2))) then
         call fatal_error(error, "2D periodic lattice vectors contain non-finite values")
         return
      end if
      len1 = norm2(a1)
      len2 = norm2(a2)
      area = norm2(cross_product_3d(a1, a2))
      if (len1 <= tol_deg .or. len2 <= tol_deg .or. area <= tol_deg * len1 * len2) then
         call fatal_error(error, "Degenerate 2D periodic lattice vectors")
         return
      end if
   case (3)
      a1 = mol%lattice(:, 1)
      a2 = mol%lattice(:, 2)
      a3 = mol%lattice(:, 3)
      if (.not. all(ieee_is_finite(mol%lattice))) then
         call fatal_error(error, "3D periodic lattice vectors contain non-finite values")
         return
      end if
      len1 = norm2(a1)
      len2 = norm2(a2)
      vol = abs(dot_product(a1, cross_product_3d(a2, a3)))
      if (len1 <= tol_deg .or. len2 <= tol_deg .or. norm2(a3) <= tol_deg &
         & .or. vol <= tol_deg * len1 * len2 * norm2(a3)) then
         call fatal_error(error, "Degenerate 3D periodic lattice vectors")
         return
      end if
   end select
end subroutine check_structure

!> Bound multi-cell integer coordinate offsets of atom pairs in the periodic basis
pure subroutine get_cell_shift_rep(mol, ndim, p_axes, u1, u2, rep_extra)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of periodic dimensions
   integer, intent(in) :: ndim
   !> Indices of active periodic axes
   integer, intent(in) :: p_axes(3)
   !> First reduced periodic basis vector (used for 1D/2D)
   real(wp), intent(in) :: u1(3)
   !> Second reduced periodic basis vector (used for 2D)
   real(wp), intent(in) :: u2(3)
   !> Extra integer translation shells per axis
   integer, intent(out) :: rep_extra(3)

   integer :: iat, jat
   real(wp) :: rij(3), l11, g11, g22, g12, det, b1, b2, s1, s2
   real(wp), parameter :: tol_deg = 1.0e-10_wp
   real(wp), parameter :: tol_cell = 1.0_wp + 1.0e-6_wp

   rep_extra = 0
   if (mol%nat <= 1) return

   select case (ndim)
   case (1)
      l11 = dot_product(u1, u1)
      if (l11 <= tol_deg) return
      do iat = 1, mol%nat
         do jat = 1, iat - 1
            rij = mol%xyz(:, iat) - mol%xyz(:, jat)
            s1 = abs(dot_product(rij, u1) / l11)
            if (s1 > tol_cell) then
               rep_extra(p_axes(1)) = max(rep_extra(p_axes(1)), ceiling(s1 - 0.5_wp))
            end if
         end do
      end do
   case (2)
      g11 = dot_product(u1, u1)
      g22 = dot_product(u2, u2)
      g12 = dot_product(u1, u2)
      det = g11 * g22 - g12 * g12
      if (g11 <= tol_deg .or. g22 <= tol_deg .or. det <= tol_deg * g11 * g22) return
      do iat = 1, mol%nat
         do jat = 1, iat - 1
            rij = mol%xyz(:, iat) - mol%xyz(:, jat)
            b1 = dot_product(rij, u1)
            b2 = dot_product(rij, u2)
            s1 = abs((g22 * b1 - g12 * b2) / det)
            s2 = abs((g11 * b2 - g12 * b1) / det)
            if (s1 > tol_cell) then
               rep_extra(p_axes(1)) = max(rep_extra(p_axes(1)), ceiling(s1 - 0.5_wp))
            end if
            if (s2 > tol_cell) then
               rep_extra(p_axes(2)) = max(rep_extra(p_axes(2)), ceiling(s2 - 0.5_wp))
            end if
         end do
      end do
   end select

end subroutine get_cell_shift_rep

!> Generate direct lattice translation vectors strictly in the periodic subspace
subroutine get_dir_trans(mol, trans, cutoff, wsc_wrapped)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   !> Optional lattice-vector cutoff
   real(wp), intent(in), optional :: cutoff

   !> Whether pair displacements are already wrapped into the Wigner-Seitz cell
   logical, intent(in), optional :: wsc_wrapped

   integer :: ndim, p_axes(3), rep(3), rep_extra(3)
   real(wp) :: eff_lat(3, 3), a1(3), a2(3), u1(3), u2(3), len1, len2, area
   real(wp), parameter :: tol_deg = 1.0e-10_wp
   logical :: need_shift_rep

   if (.not. allocated(mol%periodic)) then
      allocate(trans(3, 1), source=0.0_wp)
      return
   end if

   need_shift_rep = .true.
   if (present(wsc_wrapped)) need_shift_rep = .not. wsc_wrapped

   call get_periodic_axes(mol%periodic, ndim, p_axes)
   rep = 0
   rep_extra = 0
   eff_lat = 0.0_wp

   select case (ndim)
   case default
      allocate(trans(3, 1), source=0.0_wp)
      return
   case (1)
      a1 = mol%lattice(:, p_axes(1))
      eff_lat(:, p_axes(1)) = a1
      if (need_shift_rep) call get_cell_shift_rep(mol, ndim, p_axes, a1, a1, rep_extra)
      if (present(cutoff)) then
         len1 = norm2(a1)
         if (len1 > tol_deg) then
            rep(p_axes(1)) = ceiling(cutoff / len1) + rep_extra(p_axes(1))
         end if
      else
         rep(p_axes(1)) = 2 + rep_extra(p_axes(1))
      end if
   case (2)
      a1 = mol%lattice(:, p_axes(1))
      a2 = mol%lattice(:, p_axes(2))
      call reduce_basis_2d(a1, a2, u1, u2)
      eff_lat(:, p_axes(1)) = u1
      eff_lat(:, p_axes(2)) = u2
      if (need_shift_rep) call get_cell_shift_rep(mol, ndim, p_axes, u1, u2, rep_extra)
      if (present(cutoff)) then
         len1 = norm2(u1)
         len2 = norm2(u2)
         area = norm2(cross_product_3d(u1, u2))
         if (len1 > tol_deg .and. len2 > tol_deg .and. area > tol_deg * len1 * len2) then
            rep(p_axes(1)) = ceiling(cutoff * len2 / area) + rep_extra(p_axes(1))
            rep(p_axes(2)) = ceiling(cutoff * len1 / area) + rep_extra(p_axes(2))
         end if
      else
         rep(p_axes(1)) = 2 + rep_extra(p_axes(1))
         rep(p_axes(2)) = 2 + rep_extra(p_axes(2))
      end if
   case (3)
      eff_lat = mol%lattice
      if (present(cutoff)) then
         call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)
         return
      else
         rep = [2, 2, 2]
      end if
   end select

   call get_lattice_points(eff_lat, rep, .true., trans)

end subroutine get_dir_trans

!> Construct a Wigner-Seitz cell using only periodic lattice translations
subroutine init_wignerseitz_cell(self, mol)
   !> Wigner-Seitz cell instance
   type(wignerseitz_cell), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), parameter :: thr = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: tol_lowd = 1.0e-8_wp
   integer :: iat, jat, ntr, nimg, ndim, p_axes(3)
   integer, allocatable :: tridx(:)
   real(wp) :: vec(3)
   real(wp), allocatable :: trans(:, :)
   logical :: use_lowd_tol

   use_lowd_tol = .false.
   if (allocated(mol%periodic)) then
      call get_periodic_axes(mol%periodic, ndim, p_axes)
      use_lowd_tol = (ndim == 1 .or. ndim == 2)
   end if

   call get_dir_trans(mol, trans, thr)
   ntr = size(trans, 2)
   allocate(self%nimg(mol%nat, mol%nat), self%tridx(ntr, mol%nat, mol%nat), &
      & tridx(ntr))

   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(mol, trans, self, use_lowd_tol) private(iat, jat, vec, nimg, tridx)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         if (use_lowd_tol) then
            call get_wsc_pairs(nimg, trans, vec, tridx, tol=tol_lowd)
         else
            call get_wsc_pairs(nimg, trans, vec, tridx)
         end if
         self%nimg(jat, iat) = nimg
         self%tridx(:, jat, iat) = tridx
      end do
   end do

   if (mol%nat > 0) self%nimg_max = maxval(self%nimg)
   call move_alloc(trans, self%trans)

end subroutine init_wignerseitz_cell

!> Find periodic images with minimum distance for Wigner-Seitz cell construction.
!> Safe for 1D/2D where all candidate images may be equidistant.
pure subroutine get_wsc_pairs(iws, trans, rij, list, tol)
   !> Number of equivalent minimum-distance images
   integer, intent(out) :: iws
   !> Lattice translation vectors
   real(wp), intent(in) :: trans(:, :)
   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)
   !> Indices of equivalent minimum-distance images in trans
   integer, intent(out) :: list(:)
   !> Optional relative squared-distance tie tolerance (defaults to legacy 0.3 Bohr^2)
   real(wp), intent(in), optional :: tol

   integer :: itr, img, pos
   integer, allocatable :: index(:)
   logical, allocatable :: mask(:)
   real(wp), allocatable :: dist(:)
   real(wp) :: vec(3), r2, tol_eff
   real(wp), parameter :: thr = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: default_tol = 0.3_wp

   list = 0
   iws = 0
   allocate(dist(size(trans, 2)), mask(size(trans, 2)), index(size(trans, 2)))

   img = 0
   do itr = 1, size(trans, 2)
      vec(:) = rij - trans(:, itr)
      r2 = dot_product(vec, vec)
      if (r2 < thr) cycle
      img = img + 1
      index(img) = itr
      dist(img) = r2
   end do

   if (img == 0) return

   mask = .true.
   pos = minloc(dist(:img), dim=1, mask=mask(:img))
   mask(pos) = .false.
   r2 = dist(pos)
   if (present(tol)) then
      tol_eff = tol * max(1.0_wp, r2)
   else
      tol_eff = default_tol
   end if
   iws = 1
   list(iws) = index(pos)

   do
      if (iws >= img) exit
      pos = minloc(dist(:img), dim=1, mask=mask(:img))
      if (pos <= 0) exit
      if (dist(pos) - r2 > tol_eff) exit
      mask(pos) = .false.
      iws = iws + 1
      list(iws) = index(pos)
   end do

end subroutine get_wsc_pairs

subroutine get_rec_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = [2, 2, 2]
   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi * transpose(matinv_3x3(lattice))
   call get_lattice_points(rec_lat, rep, .false., trans)

end subroutine get_rec_trans

subroutine solve(self, mol, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
   & energy, gradient, sigma, qvec, dqdr, dqdL, response)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Coordination number
   real(wp), intent(in), contiguous :: cn(:)
   !> Local atomic partial charges
   real(wp), intent(in), contiguous :: qloc(:)
   !> Optional derivative of the coordination number w.r.t. atomic positions
   real(wp), intent(in), contiguous, optional :: dcndr(:, :, :)
   !> Optional derivative of the coordination number w.r.t. lattice vectors
   real(wp), intent(in), contiguous, optional :: dcndL(:, :, :)
   !> Optional derivative of the local atomic partial charges w.r.t. atomic positions
   real(wp), intent(in), contiguous, optional :: dqlocdr(:, :, :)
   !> Optional derivative of the local atomic partial charges w.r.t. lattice vectors
   real(wp), intent(in), contiguous, optional :: dqlocdL(:, :, :)
   !> Optional atomic partial charges result
   real(wp), intent(out), contiguous, optional :: qvec(:)
   !> Optional electrostatic energy result
   real(wp), intent(inout), contiguous, optional :: energy(:)
   !> Optional gradient for electrostatic energy
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)
   !> Optional stress tensor for electrostatic energy
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)
   !> Optional derivative of the atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)
   !> Optional derivative of the atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   !> Retain data for a subsequent weighted second-response calculation
   type(cache_container), optional, target, intent(inout) :: response

   integer :: ic, jc, iat, ndim
   logical :: grad, cpq, dcn
   integer(ik) :: info
   integer(ik), allocatable :: ipiv(:)

   ! Variables for solving ES equation
   real(wp), allocatable :: xvec(:), vrhs(:), amat(:, :)
   real(wp), allocatable :: ainv(:, :), jmat(:, :)
   ! Gradients
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: dxdr(:, :, :), dxdL(:, :, :)
   type(cache_container), target :: local_cache
   type(cache_container), pointer :: cache

   ! Calculate gradient if the respective arrays are present
   dcn = present(dcndr) .and. present(dcndL)
   grad = present(gradient) .and. present(sigma) .and. dcn
   cpq = ((present(dqdr) .and. present(dqdL)) .or. present(response)) .and. dcn

   call check_structure(mol, error)
   if (allocated(error)) return
   if (.not. allocated(self%ncoord_en) .and. allocated(mol%periodic)) then
      if ((count(mol%periodic) == 1 .or. count(mol%periodic) == 2) &
         & .and. abs(mol%charge) > 1.0e-10_wp) then
         call fatal_error(error, &
            & "Charged 1D/2D periodic unit cells are not supported by EEQ; use EEQBC instead")
         return
      end if
   end if

   cache => local_cache
   if (present(response)) then
      cache => response
      if (allocated(cache%qvec)) deallocate(cache%qvec)
      if (allocated(cache%ainv)) deallocate(cache%ainv)
      if (allocated(cache%dqdr)) deallocate(cache%dqdr)
      if (allocated(cache%raw)) deallocate(cache%raw)
      if (allocated(cache%error)) deallocate(cache%error)
      if (.not. dcn) then
         call fatal_error(error, "Retaining charge response requires CN derivatives")
         return
      end if
   end if

   ! Update cache
   call self%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   if (allocated(cache%error)) then
      call move_alloc(cache%error, error)
      return
   end if

   ! Setup the Coulomb matrix
   ndim = mol%nat + 1
   allocate(amat(ndim, ndim))
   call self%get_coulomb_matrix(mol, cache, amat)
   if (allocated(cache%error)) then
      call move_alloc(cache%error, error)
      return
   end if

   ! Get RHS of ES equation
   allocate(xvec(ndim))
   call self%get_xvec(mol, cache, xvec)

   vrhs = xvec
   ainv = amat

   ! Factorize the Coulomb matrix
   allocate(ipiv(ndim))
   call sytrf(ainv, ipiv, info=info, uplo="l")
   if (info /= 0) then
      call fatal_error(error, "Bunch-Kaufman factorization failed.")
      return
   end if

   if (cpq) then
      ! Inverted matrix is needed for coupled-perturbed equations
      call sytri(ainv, ipiv, info=info, uplo="l")
      if (info /= 0) then
         call fatal_error(error, "Inversion of factorized matrix failed.")
         return
      end if
      ! Solve the linear system
      call symv(ainv, xvec, vrhs, uplo="l")
      do ic = 1, ndim
         do jc = ic + 1, ndim
            ainv(ic, jc) = ainv(jc, ic)
         end do
      end do
   else
      ! Solve the linear system
      call sytrs(ainv, vrhs, ipiv, info=info, uplo="l")
      if (info /= 0) then
         call fatal_error(error, "Solution of linear system failed.")
         return
      end if

   end if

   if (present(qvec)) then
      qvec(:) = vrhs(:mol%nat)
   end if

   if (present(energy)) then
      ! Extract only the Coulomb matrix without the constraints
      allocate(jmat(mol%nat, mol%nat))
      jmat = amat(:mol%nat, :mol%nat)
      call symv(jmat, vrhs(:mol%nat), xvec(:mol%nat), &
         & alpha=0.5_wp, beta=-1.0_wp, uplo="l")
      energy(:) = energy(:) + vrhs(:mol%nat) * xvec(:mol%nat)
   end if

   ! Allocate and get amat derivatives
   if (grad .or. cpq) then
      allocate(dadr(3, mol%nat, ndim), dadL(3, 3, ndim), atrace(3, mol%nat))
      allocate(dxdr(3, mol%nat, ndim), dxdL(3, 3, ndim))
      call self%get_xvec_derivs(mol, cache, dxdr, dxdL)
      call self%get_coulomb_derivs(mol, cache, vrhs, dadr, dadL, atrace)
      if (allocated(cache%error)) then
         call move_alloc(cache%error, error)
         return
      end if
      do iat = 1, mol%nat
         dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
      end do
   end if

   if (grad) then
      gradient = 0.0_wp
      call gemv(dadr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=0.5_wp)
      call gemv(dxdr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=-1.0_wp)
      call gemv(dadL, vrhs, sigma, beta=1.0_wp, alpha=0.5_wp)
      call gemv(dxdL, vrhs, sigma, beta=1.0_wp, alpha=-1.0_wp)
   end if

   if (cpq) then
      do iat = 1, mol%nat
         dadr(:, :, iat) = -dxdr(:, :, iat) + dadr(:, :, iat)
         dadL(:, :, iat) = -dxdL(:, :, iat) + dadL(:, :, iat)
      end do

      if (present(dqdL)) call gemm(dadL, ainv(:, :mol%nat), dqdL, alpha=-1.0_wp)
      if (present(response)) then
         allocate(cache%dqdr(3, mol%nat, mol%nat))
         call gemm(dadr, ainv(:, :mol%nat), cache%dqdr, alpha=-1.0_wp)
         if (present(dqdr)) dqdr(:, :, :) = cache%dqdr
         call move_alloc(vrhs, cache%qvec)
         call move_alloc(ainv, cache%ainv)
      else
         call gemm(dadr, ainv(:, :mol%nat), dqdr, alpha=-1.0_wp)
      end if
   end if
end subroutine solve

!> Add the Cartesian charge Hessian contracted with fixed weights and lattice.
!> Ay=x, A^T u=(weights,0): H_ab=u^T(x_ab-A_ab y-A_a y_b-A_b y_a).
subroutine add_charge_hessian(self, mol, error, weights, hessian, cache)
   use multicharge_model_hessian, only: add_symmetric_product
   class(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(in) :: weights(:)
   real(wp), intent(inout) :: hessian(:, :)
   type(cache_container), target, intent(inout) :: cache

   integer :: nat, iat
   real(wp), allocatable :: adj(:), rhs(:), dadL(:, :,:), atrace(:, :)
   real(wp), allocatable, target :: dadr(:, :, :)
   real(wp), pointer :: da(:, :), dq(:, :)

   nat = mol%nat
   call check_structure(mol, error)
   if (allocated(error)) return
   if (.not. allocated(self%ncoord_en) .and. allocated(mol%periodic)) then
      if ((count(mol%periodic) == 1 .or. count(mol%periodic) == 2) &
         & .and. abs(mol%charge) > 1.0e-10_wp) then
         call fatal_error(error, &
            & "Charged 1D/2D periodic unit cells are not supported by EEQ; use EEQBC instead")
         return
      end if
   end if
   if (size(weights) /= nat .or. any(shape(hessian) /= 3*nat)) then
      call fatal_error(error, "Invalid dimensions for charge Hessian contraction")
      return
   end if
   if (.not. allocated(cache%qvec) .or. .not. allocated(cache%ainv) &
      & .or. .not. allocated(cache%dqdr) .or. .not. allocated(cache%raw)) then
      call fatal_error(error, "Charge Hessian requires a response cache from get_charges")
      return
   end if
   if (size(cache%qvec) /= nat+1 .or. any(shape(cache%ainv) /= nat+1) &
      & .or. any(shape(cache%dqdr) /= [3, nat, nat])) then
      call fatal_error(error, "Charge response cache has incompatible dimensions")
      return
   end if

   allocate(adj(nat+1), rhs(nat+1))
   rhs(:nat) = weights
   rhs(nat+1) = 0.0_wp
   call symv(cache%ainv, rhs, adj, uplo="l")
   call self%add_model_hessian(mol, error, cache, adj, hessian)
   if (allocated(error)) return

   ! Reuse the first derivative action with the adjoint in place of the charges.
   allocate(dadr(3, nat, nat+1), dadL(3, 3, nat+1), atrace(3, nat))
   call self%get_coulomb_derivs(mol, cache, adj, dadr, dadL, atrace)
   if (allocated(cache%error)) then
      call move_alloc(cache%error, error)
      return
   end if
   do iat = 1, nat
      dadr(:, iat, iat) = dadr(:, iat, iat) + atrace(:, iat)
   end do
   da(1:3*nat, 1:nat+1) => dadr
   dq(1:3*nat, 1:nat) => cache%dqdr
   call add_symmetric_product(da(:, :nat), dq, hessian, alpha=-1.0_wp)
end subroutine add_charge_hessian


!> Concrete models override this curvature contraction when supported.
subroutine add_model_hessian(self, mol, error, cache, adj, hessian)
   class(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: adj(:)
   real(wp), intent(inout) :: hessian(:, :)

   call fatal_error(error, "Analytical charge Hessian is not implemented for this model")
end subroutine add_model_hessian


subroutine local_charge(self, mol, trans, qloc, dqlocdr, dqlocdL)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Local atomic partial charges
   real(wp), intent(out) :: qloc(:)
   !> Optional derivative of local atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), optional :: dqlocdr(:, :, :)
   !> Optional derivative of local atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out), optional :: dqlocdL(:, :, :)

   qloc = 0.0_wp
   if (present(dqlocdr) .and. present(dqlocdL)) then
      dqlocdr = 0.0_wp
      dqlocdL = 0.0_wp
   end if
   ! Get the electronegativity weighted CN for local charge
   ! Derivatives depend only in this CN
   if (allocated(self%ncoord_en)) then
      call self%ncoord_en%get_coordination_number(mol, trans, qloc, dqlocdr, dqlocdL)
   end if

   ! Distribute the total charge equally
   qloc = qloc + mol%charge / real(mol%nat, wp)

end subroutine local_charge

end module multicharge_model_type
