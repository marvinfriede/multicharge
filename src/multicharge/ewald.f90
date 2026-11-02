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

module multicharge_ewald
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matdet_3x3, matinv_3x3
   implicit none
   private

   public :: get_alpha, get_alpha_2d, erfc_exp_2d
   public :: get_ewald_dir_trans_2d, get_ewald_rec_trans_2d

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   abstract interface
      !> Returns the max. value of a term in the reciprocal space part of the Ewald
      !> summation for a given vector length.
      pure function get_rec_term_gen(gg, alpha, vol) result(gTerm)
         import :: wp

         !> Length of the reciprocal space vector
         real(wp), intent(in) :: gg

         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha

         !> Volume (or area/length) of the real space unit cell
         real(wp), intent(in) :: vol

         !> Reciprocal term
         real(wp) :: gTerm

      end function get_rec_term_gen
   end interface

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


subroutine get_alpha(lattice, alpha, periodic)
   !> Lattice vectors (columns)
   real(wp), intent(in) :: lattice(:, :)
   !> Ewald splitting parameter
   real(wp), intent(out) :: alpha
   !> Optional periodicity mask
   logical, intent(in), optional :: periodic(:)

   integer :: ndim, p_axes(3)
   real(wp) :: vol, rec_lat(3, 3)

   if (present(periodic)) then
      call get_periodic_axes(periodic, ndim, p_axes)
      select case (ndim)
      case (0)
         alpha = 0.25_wp
         return
      case (1)
         call get_alpha_1d_subspace(lattice(:, p_axes(1)), alpha)
         return
      case (2)
         call get_alpha_2d(lattice(:, p_axes(1)), lattice(:, p_axes(2)), alpha)
         return
      case (3)
         continue
      end select
   end if

   vol = abs(matdet_3x3(lattice))
   rec_lat = twopi * transpose(matinv_3x3(lattice))

   call search_alpha(lattice, rec_lat, vol, eps, alpha, get_rec_term_3d)

end subroutine get_alpha


!> Select optimal Ewald alpha for a 1D periodic wire from its periodic vector
subroutine get_alpha_1d_subspace(a1, alpha)
   !> Periodic lattice vector
   real(wp), intent(in) :: a1(3)
   !> Optimal Ewald splitting parameter
   real(wp), intent(out) :: alpha

   real(wp) :: len1, p_lat(3, 1), rec_lat(3, 1)

   len1 = norm2(a1)
   p_lat(:, 1) = a1
   rec_lat(:, 1) = (twopi / (len1 * len1)) * a1

   call search_alpha(p_lat, rec_lat, len1, eps, alpha, get_rec_term_1d)

end subroutine get_alpha_1d_subspace


!> Select optimal Ewald alpha for a 2D periodic slab from its two periodic vectors
subroutine get_alpha_2d(a1, a2, alpha)
   !> First periodic lattice vector
   real(wp), intent(in) :: a1(3)
   !> Second periodic lattice vector
   real(wp), intent(in) :: a2(3)
   !> Optimal Ewald splitting parameter
   real(wp), intent(out) :: alpha

   real(wp) :: area, m11, m12, m22, p_lat(3, 2), rec_lat(3, 2)

   area = norm2(cross_product_3d(a1, a2))
   m11 = dot_product(a1, a1)
   m12 = dot_product(a1, a2)
   m22 = dot_product(a2, a2)

   p_lat(:, 1) = a1
   p_lat(:, 2) = a2
   rec_lat(:, 1) = (twopi / (area * area)) * (m22 * a1 - m12 * a2)
   rec_lat(:, 2) = (twopi / (area * area)) * (-m12 * a1 + m11 * a2)

   call search_alpha(p_lat, rec_lat, area, eps, alpha, get_rec_term_2d)

end subroutine get_alpha_2d


!> Overflow-safe evaluation of exp(G * s) * erfc(G / (2 * alpha) + alpha * s)
pure function erfc_exp_2d(g_norm, s, alpha) result(val)
   !> Norm of the 2D reciprocal vector (G > 0)
   real(wp), intent(in) :: g_norm
   !> Signed out-of-plane coordinate (+z or -z)
   real(wp), intent(in) :: s
   !> Ewald splitting parameter (alpha > 0)
   real(wp), intent(in) :: alpha
   !> Value of exp(G * s) * erfc(G / (2 * alpha) + alpha * s)
   real(wp) :: val

   real(wp) :: u, exponent

   u = 0.5_wp * g_norm / alpha + alpha * s
   if (u > 0.0_wp) then
      exponent = -0.25_wp * (g_norm / alpha)**2 - (alpha * s)**2
      val = exp(exponent) * erfc_scaled(u)
   else
      val = exp(g_norm * s) * erfc(u)
   end if

end function erfc_exp_2d


!> Generate converged direct-space lattice translations for 2D Ewald
pure subroutine get_ewald_dir_trans_2d(a1, a2, alpha, gam_min, trans)
   !> First periodic lattice vector
   real(wp), intent(in) :: a1(3)
   !> Second periodic lattice vector
   real(wp), intent(in) :: a2(3)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Minimum Gaussian exponent min_ij(gamma_ij)
   real(wp), intent(in) :: gam_min
   !> Direct-space lattice translations (3, ntrans)
   real(wp), allocatable, intent(out) :: trans(:, :)

   real(wp) :: kappa, rcut, area, len1, len2
   integer :: rep1, rep2, ntr, i1, i2, idx

   kappa = min(alpha, gam_min)
   rcut = 6.0_wp / kappa
   len1 = norm2(a1)
   len2 = norm2(a2)
   area = norm2(cross_product_3d(a1, a2))

   rep1 = min(20, max(2, ceiling(rcut * len2 / area)))
   rep2 = min(20, max(2, ceiling(rcut * len1 / area)))

   ntr = (2 * rep1 + 1) * (2 * rep2 + 1)
   allocate(trans(3, ntr))
   idx = 0
   do i1 = -rep1, rep1
      do i2 = -rep2, rep2
         idx = idx + 1
         trans(:, idx) = real(i1, wp) * a1 + real(i2, wp) * a2
      end do
   end do

end subroutine get_ewald_dir_trans_2d


!> Generate converged nonzero reciprocal-space lattice vectors for 2D Ewald
pure subroutine get_ewald_rec_trans_2d(a1, a2, alpha, trans, area, n_hat)
   !> First periodic lattice vector
   real(wp), intent(in) :: a1(3)
   !> Second periodic lattice vector
   real(wp), intent(in) :: a2(3)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Nonzero reciprocal lattice vectors (3, nrec)
   real(wp), allocatable, intent(out) :: trans(:, :)
   !> Physical area of the 2D unit cell
   real(wp), intent(out) :: area
   !> Unit normal to the 2D periodic plane
   real(wp), intent(out) :: n_hat(3)

   real(wp) :: n_vec(3), m11, m12, m22, b1(3), b2(3), gcut
   integer :: rep1, rep2, ntr, i1, i2, idx

   n_vec = cross_product_3d(a1, a2)
   area = norm2(n_vec)
   n_hat = n_vec / area

   m11 = dot_product(a1, a1)
   m12 = dot_product(a1, a2)
   m22 = dot_product(a2, a2)
   b1 = (twopi / (area * area)) * (m22 * a1 - m12 * a2)
   b2 = (twopi / (area * area)) * (-m12 * a1 + m11 * a2)

   gcut = 12.0_wp * alpha
   rep1 = min(20, max(2, ceiling(gcut * sqrt(m11) / twopi)))
   rep2 = min(20, max(2, ceiling(gcut * sqrt(m22) / twopi)))

   ntr = (2 * rep1 + 1) * (2 * rep2 + 1) - 1
   allocate(trans(3, ntr))
   idx = 0
   do i1 = -rep1, rep1
      do i2 = -rep2, rep2
         if (i1 == 0 .and. i2 == 0) cycle
         idx = idx + 1
         trans(:, idx) = real(i1, wp) * b1 + real(i2, wp) * b2
      end do
   end do

end subroutine get_ewald_rec_trans_2d


!> Get optimal alpha-parameter for the Ewald summation by finding alpha, where
!> decline of real and reciprocal part of Ewald are equal.
subroutine search_alpha(lattice, rec_lat, volume, tolerance, alpha, get_rec_term)
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)
   !> Reciprocal vectors
   real(wp), intent(in) :: rec_lat(:, :)
   !> Volume (or area/length) of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance for difference in real and rec. part
   real(wp), intent(in) :: tolerance
   !> Optimal alpha
   real(wp), intent(out) :: alpha
   !> Reciprocal term evaluator
   procedure(get_rec_term_gen) :: get_rec_term

   real(wp) :: alpl, alpr, rlen, dlen, diff
   real(wp), parameter :: alpha0 = 1.0e-8_wp
   integer, parameter :: niter = 30
   integer :: ibs, stat

   rlen = sqrt(minval(sum(rec_lat(:, :)**2, dim=1)))
   dlen = sqrt(minval(sum(lattice(:, :)**2, dim=1)))

   stat = 0
   alpha = alpha0
   diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
   end do
   if (alpha > huge(1.0_wp)) then
      stat = 1
   else if (alpha == alpha0) then
      stat = 2
   end if

   if (stat == 0) then
      alpl = 0.5_wp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_wp))
         alpha = 2.0_wp * alpha
         diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         stat = 3
      end if
   end if

   if (stat == 0) then
      alpr = alpha
      alpha = (alpl + alpr) * 0.5_wp
      ibs = 0
      diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
      do while (abs(diff) > tolerance .and. ibs <= niter)
         if (diff < 0) then
            alpl = alpha
         else
            alpr = alpha
         end if
         alpha = (alpl + alpr) * 0.5_wp
         diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
         ibs = ibs + 1
      end do
      if (ibs > niter) then
         stat = 4
      end if
   end if

   if (stat /= 0) then
      alpha = 0.25_wp
   end if

end subroutine search_alpha

!> Returns the difference in the decrease of the real and reciprocal parts of the
!> Ewald sum. In order to make the real space part shorter than the reciprocal
!> space part, the values are taken at different distances for the real and the
!> reciprocal space parts.
pure function rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume) result(diff)

   !> Parameter for the Ewald summation
   real(wp), intent(in) :: alpha

   !> Procedure pointer to reciprocal routine
   procedure(get_rec_term_gen) :: get_rec_term

   !> Length of the shortest reciprocal space vector in the sum
   real(wp), intent(in) :: rlen

   !> Length of the shortest real space vector in the sum
   real(wp), intent(in) :: dlen

   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume

   !> Difference between changes in the two terms
   real(wp) :: diff

   diff = ((get_rec_term(4*rlen, alpha, volume) &
      & - get_rec_term(5*rlen, alpha, volume))) &
      & - (get_dir_term(2*dlen, alpha) - get_dir_term(3*dlen, alpha))

end function rec_dir_diff

!> Returns the max. value of a term in the real space part of the Ewald summation
!> for a given vector length.
pure function get_dir_term(rr, alpha) result(dval)

   !> Length of the real space vector
   real(wp), intent(in) :: rr

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Real space term
   real(wp) :: dval

   dval = erfc(alpha*rr)/rr

end function get_dir_term


!> Returns the max. value of a term in the 3D reciprocal space part of the Ewald
!> summation for a given vector length.
pure function get_rec_term_3d(gg, alpha, vol) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol

   !> Reciprocal term
   real(wp) :: rval

   rval = 4.0_wp*pi*(exp(-0.25_wp*gg*gg/(alpha**2))/(vol*gg*gg))

end function get_rec_term_3d


!> Returns the max. value of a term in the 2D reciprocal space part of the Ewald
!> summation for a given vector length (at z = 0).
pure function get_rec_term_2d(gg, alpha, area) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Area of the 2D real space unit cell
   real(wp), intent(in) :: area

   !> Reciprocal term
   real(wp) :: rval

   rval = twopi * erfc(0.5_wp * gg / alpha) / (area * gg)

end function get_rec_term_2d


!> Returns the max. value of a term in the 1D reciprocal space part of the Ewald
!> summation for a given vector length (at rho = 0).
pure function get_rec_term_1d(gg, alpha, len1) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Length of the 1D real space unit cell
   real(wp), intent(in) :: len1

   !> Reciprocal term
   real(wp) :: rval
   real(wp) :: u

   u = 0.25_wp * gg * gg / (alpha * alpha)
   rval = 2.0_wp * exp(-u) / (len1 * (u + 1.0_wp))

end function get_rec_term_1d

end module multicharge_ewald
