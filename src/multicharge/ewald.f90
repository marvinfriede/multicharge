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
   public :: get_ewald_dir_trans_1d, get_ewald_rec_trans_1d
   public :: fast_f012_1d, fast_k012_1d

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: euler_gamma = 0.57721566490153286060_wp

   !> Positive nodes of 16-point Gauss-Legendre quadrature on [-1, 1]
   real(wp), parameter :: gl16_x(8) = [ &
      & 0.095012509837637440185_wp, &
      & 0.281603550779258913230_wp, &
      & 0.458016777657227386342_wp, &
      & 0.617876244402643748447_wp, &
      & 0.755404408355003033895_wp, &
      & 0.865631202387831743880_wp, &
      & 0.944575023073232576078_wp, &
      & 0.989400934991649932596_wp]

   !> Weights of 16-point Gauss-Legendre quadrature on [-1, 1]
   real(wp), parameter :: gl16_w(8) = [ &
      & 0.189450610455068496285_wp, &
      & 0.182603415044923588867_wp, &
      & 0.169156519395002538189_wp, &
      & 0.149595988816576732081_wp, &
      & 0.124628971255533872052_wp, &
      & 0.095158511682492784810_wp, &
      & 0.062253523938647892863_wp, &
      & 0.027152459411754094852_wp]

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


!> Generate converged direct-space lattice translations for 1D Ewald
pure subroutine get_ewald_dir_trans_1d(a1, alpha, gam_min, trans)
   !> Periodic lattice vector
   real(wp), intent(in) :: a1(3)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Minimum Gaussian exponent min_ij(gamma_ij)
   real(wp), intent(in) :: gam_min
   !> Direct-space lattice translations (3, ntrans)
   real(wp), allocatable, intent(out) :: trans(:, :)

   real(wp) :: kappa, rcut, len1
   integer :: rep1, ntr, i1, idx

   kappa = min(alpha, gam_min)
   rcut = 6.0_wp / kappa
   len1 = norm2(a1)
   rep1 = min(40, max(2, ceiling(rcut / len1)))

   ntr = 2 * rep1 + 1
   allocate(trans(3, ntr))
   idx = 0
   do i1 = -rep1, rep1
      idx = idx + 1
      trans(:, idx) = real(i1, wp) * a1
   end do

end subroutine get_ewald_dir_trans_1d


!> Generate converged nonzero reciprocal-space lattice vectors for 1D Ewald
pure subroutine get_ewald_rec_trans_1d(a1, alpha, trans, len1, u_hat)
   !> Periodic lattice vector
   real(wp), intent(in) :: a1(3)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Nonzero reciprocal lattice vectors (3, nrec)
   real(wp), allocatable, intent(out) :: trans(:, :)
   !> Physical length of the 1D unit cell
   real(wp), intent(out) :: len1
   !> Unit vector along the 1D periodic axis
   real(wp), intent(out) :: u_hat(3)

   real(wp) :: b1(3), gcut
   integer :: rep1, ntr, i1, idx

   len1 = norm2(a1)
   u_hat = a1 / len1
   b1 = (twopi / len1) * u_hat

   gcut = 12.0_wp * alpha
   rep1 = min(40, max(2, ceiling(gcut * len1 / twopi)))

   ntr = 2 * rep1
   allocate(trans(3, ntr))
   idx = 0
   do i1 = -rep1, rep1
      if (i1 == 0) cycle
      idx = idx + 1
      trans(:, idx) = real(i1, wp) * b1
   end do

end subroutine get_ewald_rec_trans_1d


!> Evaluate 1D zero-mode functions f0(u), f1(u), f2(u) for u = alpha^2 * rho^2 >= 0
pure subroutine fast_f012_1d(u, f0, f1, f2)
   !> Non-negative scalar argument u = alpha^2 * rho^2
   real(wp), intent(in) :: u
   !> Zero-mode potential function f0(u) = gamma_E + E1(u) + ln(u)
   real(wp), intent(out) :: f0
   !> First derivative f1(u) = f0'(u) = (1 - exp(-u)) / u
   real(wp), intent(out) :: f1
   !> Second derivative f2(u) = f0''(u) = ((u + 1)*exp(-u) - 1) / u^2
   real(wp), intent(out) :: f2

   real(wp) :: term, cf, e1_val, exp_u
   integer :: n, k

   if (u <= 0.0_wp) then
      f0 = 0.0_wp
      f1 = 1.0_wp
      f2 = -0.5_wp
      return
   end if

   if (u <= 4.0_wp) then
      f0 = 0.0_wp
      term = u
      do n = 1, 35
         f0 = f0 + term / real(n, wp)
         term = -term * u / real(n + 1, wp)
      end do
   else
      cf = 0.0_wp
      do k = 40, 1, -1
         cf = -real(k * k, wp) / (u + real(2 * k + 1, wp) + cf)
      end do
      e1_val = exp(-u) / (u + 1.0_wp + cf)
      f0 = euler_gamma + e1_val + log(u)
   end if

   if (u <= 0.5_wp) then
      f1 = 0.0_wp
      f2 = 0.0_wp
      term = 1.0_wp
      do n = 0, 19
         f1 = f1 + term / real(n + 1, wp)
         f2 = f2 - term / real(n + 2, wp)
         term = -term * u / real(n + 1, wp)
      end do
   else
      exp_u = exp(-u)
      f1 = (1.0_wp - exp_u) / u
      f2 = ((u + 1.0_wp) * exp_u - 1.0_wp) / (u * u)
   end if

end subroutine fast_f012_1d


!> Evaluate incomplete Bessel/exponential integrals K0(x, y), K1(x, y), K2(x, y)
!> using composite 16-point Gauss-Legendre quadrature in w = ln(t).
pure subroutine fast_k012_1d(x, y, k0, k1, k2)
   !> Positive reciprocal argument x = G^2 / (4 * alpha^2)
   real(wp), intent(in) :: x
   !> Non-negative radial argument y = alpha^2 * rho^2
   real(wp), intent(in) :: y
   !> Integral K0(x, y)
   real(wp), intent(out) :: k0
   !> Integral K1(x, y) = -dK0/dy
   real(wp), intent(out) :: k1
   !> Integral K2(x, y) = d^2K0/dy^2
   real(wp), intent(out) :: k2

   real(wp) :: w0, psi_min, w_max, s0, s1, s2

   if (y > x) then
      w0 = 0.5_wp * log(y / x)
      psi_min = 2.0_wp * sqrt(x * y)
   else
      w0 = 0.0_wp
      psi_min = x + y
   end if

   if (psi_min > 700.0_wp) then
      k0 = 0.0_wp
      k1 = 0.0_wp
      k2 = 0.0_wp
      return
   end if

   w_max = log((psi_min + 42.0_wp) / x)
   if (w_max <= 0.0_wp) then
      k0 = 0.0_wp
      k1 = 0.0_wp
      k2 = 0.0_wp
      return
   end if

   if (w0 > 0.0_wp .and. w0 < w_max) then
      call integrate_k012_interval(x, y, 0.0_wp, w0, k0, k1, k2)
      call integrate_k012_interval(x, y, w0, w_max, s0, s1, s2)
      k0 = k0 + s0
      k1 = k1 + s1
      k2 = k2 + s2
   else
      call integrate_k012_interval(x, y, 0.0_wp, w_max, k0, k1, k2)
   end if

end subroutine fast_k012_1d


!> Composite 16-point Gauss-Legendre quadrature of exp(-x*exp(w) - y*exp(-w)) on [a_start, b_end]
pure subroutine integrate_k012_interval(x, y, a_start, b_end, s0, s1, s2)
   !> Positive reciprocal argument x
   real(wp), intent(in) :: x
   !> Non-negative radial argument y
   real(wp), intent(in) :: y
   !> Lower integration bound in w = ln(t)
   real(wp), intent(in) :: a_start
   !> Upper integration bound in w = ln(t)
   real(wp), intent(in) :: b_end
   !> Accumulated K0 contribution
   real(wp), intent(out) :: s0
   !> Accumulated K1 contribution
   real(wp), intent(out) :: s1
   !> Accumulated K2 contribution
   real(wp), intent(out) :: s2

   integer :: n_sub, isub, q
   real(wp) :: h, a_sub, wt, w, ew, emw, val

   n_sub = max(1, ceiling((b_end - a_start) / 0.5_wp))
   h = (b_end - a_start) / real(n_sub, wp)
   s0 = 0.0_wp
   s1 = 0.0_wp
   s2 = 0.0_wp

   do isub = 0, n_sub - 1
      a_sub = a_start + real(isub, wp) * h
      do q = 1, 8
         wt = 0.5_wp * h * gl16_w(q)

         w = a_sub + 0.5_wp * h * (1.0_wp - gl16_x(q))
         ew = exp(w)
         emw = 1.0_wp / ew
         val = wt * exp(-x * ew - y * emw)
         s0 = s0 + val
         s1 = s1 + val * emw
         s2 = s2 + val * emw * emw

         w = a_sub + 0.5_wp * h * (1.0_wp + gl16_x(q))
         ew = exp(w)
         emw = 1.0_wp / ew
         val = wt * exp(-x * ew - y * emw)
         s0 = s0 + val
         s1 = s1 + val * emw
         s2 = s2 + val * emw * emw
      end do
   end do

end subroutine integrate_k012_interval


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
