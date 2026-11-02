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

!> @file multicharge/model/eeq.f90
!> Provides implementation of the electronegativity equilibration model (EEQ)

!> Electronegativity equlibration charge model published in
!>
!> E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth
!> and S. Grimme, *J. Chem. Phys.*, **2019**, 150, 154122.
!> DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
module multicharge_model_eeq
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matdet_3x3
   use mctc_ncoord, only: new_ncoord, cn_count
   use mctc_wignerseitz, only: wignerseitz_cell
   use multicharge_blas, only: gemm
   use multicharge_ewald, only: get_alpha, erfc_exp_2d, &
      & get_ewald_dir_trans_3d, get_ewald_rec_trans_3d, &
      & get_ewald_dir_trans_2d, get_ewald_rec_trans_2d, &
      & get_ewald_dir_trans_1d, get_ewald_rec_trans_1d, &
      & fast_f012_1d, fast_k012_1d
   use multicharge_model_cache, only: cache_container, model_cache
   use multicharge_model_hessian, only: coulomb_kernel, radial_hessian, add_pair_hessian, add_pair_block
   use multicharge_model_type, only: mchrg_model_type, get_dir_trans, &
      & init_wignerseitz_cell, get_periodic_axes
   implicit none
   private

   public :: eeq_model, new_eeq_model

   type, extends(model_cache), public :: eeq_cache
   end type eeq_cache

   type, extends(mchrg_model_type) :: eeq_model
   contains
      !> Update and allocate cache
      procedure :: update
      !> Calculate Coulomb matrix
      procedure :: get_coulomb_matrix
      !> Calculate derivatives of Coulomb matrix
      procedure :: get_coulomb_derivs
      !> Calculate right-hand side (electronegativity vector)
      procedure :: get_xvec
      !> Calculate EN vector derivatives
      procedure :: get_xvec_derivs
      !> Cartesian second-response curvature at fixed lattice
      procedure :: add_model_hessian
   end type eeq_model

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

contains

subroutine new_eeq_model(self, mol, error, chi, rad, eta, kcnchi, &
   & cutoff, cn_exp, rcov, cn_max)
   !> Electronegativity equilibration model
   type(eeq_model), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Electronegativity
   real(wp), intent(in) :: chi(:)
   !> Exponent gaussian charge
   real(wp), intent(in) :: rad(:)
   !> Chemical hardness
   real(wp), intent(in) :: eta(:)
   !> CN scaling factor for electronegativity
   real(wp), intent(in) :: kcnchi(:)
   !> Cutoff radius for coordination number
   real(wp), intent(in), optional :: cutoff
   !> Steepness of the CN counting function
   real(wp), intent(in), optional :: cn_exp
   !> Covalent radii for CN
   real(wp), intent(in), optional :: rcov(:)
   !> Maximum CN cutoff for CN
   real(wp), intent(in), optional :: cn_max

   self%chi = chi
   self%rad = rad
   self%eta = eta
   self%kcnchi = kcnchi

   call new_ncoord(self%ncoord, mol, cn_count%erf, error, &
      & cutoff=cutoff, kcn=cn_exp, rcov=rcov, cut=cn_max)

end subroutine new_eeq_model

subroutine update(self, mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in), optional :: qloc(:)
   real(wp), intent(in), optional :: dcndr(:, :, :)
   real(wp), intent(in), optional :: dcndL(:, :, :)
   real(wp), intent(in), optional :: dqlocdr(:, :, :)
   real(wp), intent(in), optional :: dqlocdL(:, :, :)

   type(eeq_cache), pointer :: ptr

   call taint(cache, ptr)

   ! Refer CN arrays in cache
   ptr%cn = cn
   if (present(dcndr) .and. present(dcndL)) then
      ptr%dcndr = dcndr
      ptr%dcndL = dcndL
   end if

   if (any(mol%periodic)) then
      ! Create WSC
      call init_wignerseitz_cell(ptr%wsc, mol)
      call get_alpha(mol%lattice, ptr%alpha, mol%periodic)
   end if

end subroutine update

subroutine get_xvec(self, mol, cache, xvec)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(out) :: xvec(:)
   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp
   real(wp) :: tmp

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, self, xvec, ptr) private(iat, izp, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      tmp = self%kcnchi(izp) / sqrt(ptr%cn(iat) + reg)
      xvec(iat) = -self%chi(izp) + tmp * ptr%cn(iat)
   end do
   xvec(mol%nat + 1) = mol%charge

end subroutine get_xvec

subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(out), contiguous :: dxdr(:, :, :)
   real(wp), intent(out), contiguous :: dxdL(:, :, :)
   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp
   real(wp) :: tmp

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   dxdr(:, :, :) = 0.0_wp
   dxdL(:, :, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, self, ptr, dxdr, dxdL) &
   !$omp private(iat, izp, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      tmp = self%kcnchi(izp) * (0.5_wp*ptr%cn(iat) + reg) / (ptr%cn(iat) + reg)**1.5_wp
      dxdr(:, :, iat) = tmp * ptr%dcndr(:, :, iat)
      dxdL(:, :, iat) = tmp * ptr%dcndL(:, :, iat)
   end do

end subroutine get_xvec_derivs

!> Compute maximum interatomic pair displacement within the Wigner-Seitz cell
pure function get_max_wsc_disp(mol, wsc) result(max_rij)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wigner-Seitz cell instance
   type(wignerseitz_cell), intent(in) :: wsc
   !> Maximum pair displacement norm
   real(wp) :: max_rij

   integer :: iat, jat, img
   real(wp) :: vec(3)

   max_rij = 0.0_wp
   if (.not. allocated(wsc%nimg) .or. .not. allocated(wsc%tridx) &
      & .or. .not. allocated(wsc%trans)) return

   do iat = 1, mol%nat
      do jat = 1, iat - 1
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) &
               & + wsc%trans(:, wsc%tridx(img, jat, iat))
            max_rij = max(max_rij, norm2(vec))
         end do
      end do
   end do
end function get_max_wsc_disp

!> u^T(x_ab-A_ab q), with u and q held fixed.
subroutine add_model_hessian(self, mol, error, cache, adj, hessian)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: adj(:)
   real(wp), intent(inout) :: hessian(:, :)

   type(eeq_cache), pointer :: ptr
   integer :: iat, jat, izp, jzp, nat, ipair, npair, ii
   real(wp), parameter :: reg = 1.0e-14_wp
   real(wp) :: c, vec(3), r, gam, val, dr, drr, weight
   real(wp), allocatable :: dEdcn(:), dcndr(:, :), work(:, :), trans(:, :)
   real(wp), allocatable :: diagonal(:, :, :)

   call view(cache, ptr)
   if (.not. associated(ptr)) then
      call fatal_error(error, "Incompatible EEQ response cache")
      return
   end if
   nat = mol%nat
   allocate(dEdcn(nat), work(3*nat, nat))
   dcndr = reshape(ptr%dcndr, [3*nat, nat])
   do iat = 1, nat
      izp = mol%id(iat)
      c = ptr%cn(iat)
      dEdcn(iat) = adj(iat)*self%kcnchi(izp)*(0.5_wp*c+reg)/(c+reg)**1.5_wp
      weight = -adj(iat)*self%kcnchi(izp)*(0.25_wp*c+reg)/(c+reg)**2.5_wp
      work(:, iat) = weight*dcndr(:, iat)
   end do
   call get_dir_trans(mol, trans, self%ncoord%cutoff)
   call self%ncoord%add_coordination_number_hessian(mol, trans, dEdcn, hessian)
   call gemm(work, dcndr, hessian, transb="t", beta=1.0_wp)

   allocate(diagonal(3, 3, nat), source=0.0_wp)
   select case (count(mol%periodic))
   case (3)
      call add_periodic_hessian(self, mol, ptr%wsc, ptr%alpha, cache%qvec, adj, hessian, diagonal, error)
      if (allocated(error)) return
   case (2)
      call add_periodic_hessian_2d(self, mol, ptr%wsc, ptr%alpha, cache%qvec, adj, hessian, diagonal, error)
      if (allocated(error)) return
   case (1)
      call add_periodic_hessian_1d(self, mol, ptr%wsc, ptr%alpha, cache%qvec, adj, hessian, diagonal, error)
      if (allocated(error)) return
   case default
      npair = nat*(nat+1)/2
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(self, mol, cache, adj, hessian, npair) &
      !$omp private(iat, jat, izp, jzp, vec, r, gam, val, dr, drr, weight) &
      !$omp reduction(+:diagonal)
      do ipair = 1, npair
         iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
         if (iat*(iat+1)/2 < ipair) iat = iat+1
         jat = ipair-iat*(iat-1)/2
         if (iat == jat) cycle
         izp = mol%id(iat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, iat)-mol%xyz(:, jat)
         r = norm2(vec)
         gam = 1.0_wp/sqrt(self%rad(izp)**2+self%rad(jzp)**2)
         call coulomb_kernel(r, gam, val, dr, drr)
         weight = -(adj(iat)*cache%qvec(jat)+adj(jat)*cache%qvec(iat))
         call add_pair_hessian(iat, jat, vec, weight*dr, weight*drr, hessian, diagonal)
      end do
      !$omp end parallel do
   end select

   !$omp parallel do default(none) schedule(static) shared(nat, hessian, diagonal) private(ii)
   do iat = 1, nat
      ii = 3*(iat-1)
      hessian(ii+1:ii+3, ii+1:ii+3) = hessian(ii+1:ii+3, ii+1:ii+3)+diagonal(:, :, iat)
   end do
   !$omp end parallel do
end subroutine add_model_hessian

!> Curvature of the same finite Ewald sums and WS average used for A.
!> The lattice, Ewald parameter, and self-image distances are fixed here.
subroutine add_periodic_hessian(self, mol, wsc, alpha, qvec, adj, hessian, diagonal, error)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha, qvec(:), adj(:)
   real(wp), intent(inout) :: hessian(:, :)
   real(wp), intent(inout) :: diagonal(:, :, :)
   type(error_type), allocatable, intent(out), optional :: error
   integer :: iat, jat, izp, jzp, img, itr, ipair, npair
   real(wp) :: rij(3), vec(3), box(3, 3), gam, gam_min, max_rij, r, weight, wsw
   real(wp) :: val, dr, drr, ar, arr
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rec_weight(:)

   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
   max_rij = get_max_wsc_disp(mol, wsc)
   call get_ewald_dir_trans_3d(mol%lattice, alpha, gam_min, dtrans, max_rij=max_rij, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_ewald_rec_trans_3d(mol%lattice, alpha, rtrans, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_rec_weights(abs(matdet_3x3(mol%lattice)), alpha, rtrans, rec_weight)

   ! Each pair owns its off-diagonal blocks. Only the diagonal blocks are
   ! reduced, avoiding a full (3N, 3N) Hessian copy for every OpenMP thread.
   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, wsc, alpha, qvec, adj, hessian, dtrans, rtrans, rec_weight, npair) &
   !$omp private(iat, jat, izp, jzp, img, itr, rij, vec, box, gam, r, weight, wsw, val, dr, drr, ar, arr) &
   !$omp reduction(+:diagonal)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      if (iat == jat) cycle
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp/sqrt(self%rad(izp)**2+self%rad(jzp)**2)
      weight = -(adj(iat)*qvec(jat)+adj(jat)*qvec(iat))
      wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
      box = 0.0_wp
      do img = 1, wsc%nimg(jat, iat)
         rij = mol%xyz(:, iat)-mol%xyz(:, jat)-wsc%trans(:, wsc%tridx(img, jat, iat))
         do itr = 1, size(dtrans, 2)
            vec = rij+dtrans(:, itr)
            r = norm2(vec)
            call coulomb_kernel(r, gam, val, dr, drr)
            call coulomb_kernel(r, alpha, val, ar, arr)
            box = box+wsw*radial_hessian(vec, dr-ar, drr-arr)
         end do
         do itr = 1, size(rtrans, 2)
            vec = rtrans(:, itr)
            box = box-wsw*rec_weight(itr)*cos(dot_product(rij, vec)) &
               & *spread(vec, 1, 3)*spread(vec, 2, 3)
         end do
      end do
      call add_pair_block(iat, jat, weight*box, hessian, diagonal)
   end do
   !$omp end parallel do
end subroutine add_periodic_hessian


!> Curvature of the 2D Ewald sum and WS average at fixed lattice and alpha.
subroutine add_periodic_hessian_2d(self, mol, wsc, alpha, qvec, adj, hessian, diagonal, error)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha, qvec(:), adj(:)
   real(wp), intent(inout) :: hessian(:, :)
   real(wp), intent(inout) :: diagonal(:, :, :)
   type(error_type), allocatable, intent(out), optional :: error

   integer :: iat, jat, izp, jzp, img, itr, ipair, npair, ndim, p_axes(3)
   real(wp) :: rij(3), vec(3), g_vec(3), box(3, 3), gam, gam_min, max_rij, r, weight, wsw
   real(wp) :: val, dr, drr, ar, arr, area, n_hat(3), a1(3), a2(3)
   real(wp) :: z, g2, g_norm, gv, cos_gv, sin_gv, ep, em, fp, fm, gauss, pref, gauss_z
   real(wp) :: gg_mat(3, 3), nn_mat(3, 3), ng_sym(3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   call get_periodic_axes(mol%periodic, ndim, p_axes)
   a1 = mol%lattice(:, p_axes(1))
   a2 = mol%lattice(:, p_axes(2))
   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
   max_rij = get_max_wsc_disp(mol, wsc)
   call get_ewald_dir_trans_2d(a1, a2, alpha, gam_min, dtrans, max_rij=max_rij, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_ewald_rec_trans_2d(a1, a2, alpha, rtrans, area, n_hat, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   nn_mat = spread(n_hat, 2, 3) * spread(n_hat, 1, 3)

   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, wsc, alpha, qvec, adj, hessian, dtrans, rtrans, area, n_hat, nn_mat, npair) &
   !$omp private(iat, jat, izp, jzp, img, itr, rij, vec, g_vec, box, gam, r, weight, wsw) &
   !$omp private(val, dr, drr, ar, arr, z, g2, g_norm, gv, cos_gv, sin_gv, ep, em, fp, fm) &
   !$omp private(gauss, pref, gauss_z, gg_mat, ng_sym) &
   !$omp reduction(+:diagonal)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      if (iat == jat) cycle
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp/sqrt(self%rad(izp)**2+self%rad(jzp)**2)
      weight = -(adj(iat)*qvec(jat)+adj(jat)*qvec(iat))
      wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
      box = 0.0_wp
      do img = 1, wsc%nimg(jat, iat)
         rij = mol%xyz(:, iat)-mol%xyz(:, jat)-wsc%trans(:, wsc%tridx(img, jat, iat))
         do itr = 1, size(dtrans, 2)
            vec = rij+dtrans(:, itr)
            r = norm2(vec)
            call coulomb_kernel(r, gam, val, dr, drr)
            call coulomb_kernel(r, alpha, val, ar, arr)
            box = box+wsw*radial_hessian(vec, dr-ar, drr-arr)
         end do
         z = dot_product(rij, n_hat)
         do itr = 1, size(rtrans, 2)
            g_vec = rtrans(:, itr)
            g2 = dot_product(g_vec, g_vec)
            if (g2 < eps) cycle
            g_norm = sqrt(g2)
            gv = dot_product(g_vec, rij)
            cos_gv = cos(gv)
            sin_gv = sin(gv)
            ep = erfc_exp_2d(g_norm, +z, alpha)
            em = erfc_exp_2d(g_norm, -z, alpha)
            fp = ep + em
            fm = ep - em
            gauss = exp(-0.25_wp * g2 / (alpha * alpha) - (alpha * z)**2)
            pref = pi / (area * g_norm)
            gg_mat = spread(g_vec, 2, 3) * spread(g_vec, 1, 3)
            ng_sym = spread(g_vec, 2, 3) * spread(n_hat, 1, 3) &
               & + spread(n_hat, 2, 3) * spread(g_vec, 1, 3)
            box = box + wsw * (-pref * cos_gv * fp * gg_mat &
               & - (pi / area) * sin_gv * fm * ng_sym &
               & + (pi / area) * cos_gv * (g_norm * fp - 4.0_wp * alpha * gauss / sqrtpi) * nn_mat)
         end do
         gauss_z = exp(-(alpha * z)**2)
         box = box - wsw * (4.0_wp * sqrtpi * alpha / area) * gauss_z * nn_mat
      end do
      call add_pair_block(iat, jat, weight*box, hessian, diagonal)
   end do
   !$omp end parallel do
end subroutine add_periodic_hessian_2d


!> Curvature of the 1D Ewald sum and WS average used for A at fixed lattice.
subroutine add_periodic_hessian_1d(self, mol, wsc, alpha, qvec, adj, hessian, diagonal, error)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha, qvec(:), adj(:)
   real(wp), intent(inout) :: hessian(:, :), diagonal(:, :, :)
   type(error_type), allocatable, intent(out), optional :: error

   integer :: iat, jat, izp, jzp, img, itr, ipair, npair, idir
   real(wp) :: gam, gam_min, max_rij, wsw, rij(3), vec(3), r, weight, box(3, 3)
   real(wp) :: val_g, dr_g, drr_g, val_a, dr_a, drr_a
   real(wp) :: a1(3), len1, u_hat(3), r_perp(3), rho2, alp2, x, y
   real(wp) :: g_vec(3), g2, gv, cos_gv, sin_gv
   real(wp) :: k0, k1, k2, f0, f1, f2
   real(wp) :: p_par(3, 3), p_perp(3, 3), gg_mat(3, 3), rg_sym(3, 3), rr_mat(3, 3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   do idir = 1, 3
      if (mol%periodic(idir)) then
         a1 = mol%lattice(:, idir)
         exit
      end if
   end do
   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
   max_rij = get_max_wsc_disp(mol, wsc)
   call get_ewald_dir_trans_1d(a1, alpha, gam_min, dtrans, max_rij=max_rij, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_ewald_rec_trans_1d(a1, alpha, rtrans, len1, u_hat, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   p_par = spread(u_hat, 2, 3) * spread(u_hat, 1, 3)
   p_perp = unity - p_par
   alp2 = alpha * alpha

   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, wsc, alpha, alp2, len1, p_perp, u_hat, dtrans, rtrans, qvec, adj, hessian, npair) &
   !$omp private(iat, jat, izp, jzp, img, itr, gam, wsw, rij, vec, r, weight, box, &
   !$omp& val_g, dr_g, drr_g, val_a, dr_a, drr_a, r_perp, rho2, x, y, &
   !$omp& g_vec, g2, gv, cos_gv, sin_gv, k0, k1, k2, f0, f1, f2, gg_mat, rg_sym, rr_mat) &
   !$omp reduction(+:diagonal)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      if (iat == jat) cycle
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp/sqrt(self%rad(izp)**2+self%rad(jzp)**2)
      weight = -(adj(iat)*qvec(jat)+adj(jat)*qvec(iat))
      wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
      box = 0.0_wp
      do img = 1, wsc%nimg(jat, iat)
         rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + wsc%trans(:, wsc%tridx(img, jat, iat)))
         do itr = 1, size(dtrans, 2)
            vec = rij - dtrans(:, itr)
            r = norm2(vec)
            call coulomb_kernel(r, gam, val_g, dr_g, drr_g)
            call coulomb_kernel(r, alpha, val_a, dr_a, drr_a)
            box = box + wsw * radial_hessian(vec, dr_g - dr_a, drr_g - drr_a)
         end do
         r_perp = rij - dot_product(rij, u_hat) * u_hat
         rho2 = dot_product(r_perp, r_perp)
         y = alp2 * rho2
         rr_mat = spread(r_perp, 2, 3) * spread(r_perp, 1, 3)
         do itr = 1, size(rtrans, 2)
            g_vec = rtrans(:, itr)
            g2 = dot_product(g_vec, g_vec)
            if (g2 < eps) cycle
            x = 0.25_wp * g2 / alp2
            gv = dot_product(g_vec, rij)
            cos_gv = cos(gv)
            sin_gv = sin(gv)
            call fast_k012_1d(x, y, k0, k1, k2)
            gg_mat = spread(g_vec, 2, 3) * spread(g_vec, 1, 3)
            rg_sym = spread(g_vec, 2, 3) * spread(r_perp, 1, 3) &
               & + spread(r_perp, 2, 3) * spread(g_vec, 1, 3)
            box = box + (wsw / len1) * (-cos_gv * k0 * gg_mat &
               & + 2.0_wp * alp2 * sin_gv * k1 * rg_sym &
               & + cos_gv * (-2.0_wp * alp2 * k1 * p_perp + 4.0_wp * alp2 * alp2 * k2 * rr_mat))
         end do
         call fast_f012_1d(y, f0, f1, f2)
         box = box + (wsw / len1) * (-2.0_wp * alp2 * f1 * p_perp - 4.0_wp * alp2 * alp2 * f2 * rr_mat)
      end do
      call add_pair_block(iat, jat, weight*box, hessian, diagonal)
   end do
   !$omp end parallel do
end subroutine add_periodic_hessian_1d


subroutine get_coulomb_matrix(self, mol, cache, amat)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(out) :: amat(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   select case (count(mol%periodic))
   case (3)
      call get_amat_3d(self, mol, ptr%wsc, ptr%alpha, amat, cache%error)
   case (2)
      call get_amat_2d(self, mol, ptr%wsc, ptr%alpha, amat, cache%error)
   case (1)
      call get_amat_1d(self, mol, ptr%wsc, ptr%alpha, amat, cache%error)
   case default
      call get_amat_0d(self, mol, amat)
   end select
end subroutine get_coulomb_matrix

subroutine get_amat_0d(self, mol, amat)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, ipair, npair
   real(wp) :: vec(3), r2, gam, tmp

   ! Every matrix element belongs to one pair; no thread-local matrix is needed.
   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat, mol, self, npair) &
   !$omp private(iat, jat, izp, jzp, gam, vec, r2, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      izp = mol%id(iat)
      if (iat == jat) then
         tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      else
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         if (r2 * gam < 1.0e-4_wp) then
            tmp = (2.0_wp * sqrt(gam) / sqrtpi) * (1.0_wp - (r2 * gam) * (1.0_wp / 3.0_wp &
               & - (r2 * gam) * (0.1_wp - (r2 * gam) / 42.0_wp)))
         else
            tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
         end if
      end if
      amat(jat, iat) = tmp
      amat(iat, jat) = tmp
   end do
   !$omp end parallel do

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
end subroutine get_amat_0d

subroutine get_amat_3d(self, mol, wsc, alpha, amat, error)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(out) :: amat(:, :)
   type(error_type), allocatable, intent(out), optional :: error

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: vec(3), gam, gam_min, max_rij, wsw, dtmp, rtmp, tmp
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rec_weight(:)

   amat = 0.0_wp
   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
   max_rij = get_max_wsc_disp(mol, wsc)
   call get_ewald_dir_trans_3d(mol%lattice, alpha, gam_min, dtrans, max_rij=max_rij, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_ewald_rec_trans_3d(mol%lattice, alpha, rtrans, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_rec_weights(abs(matdet_3x3(mol%lattice)), alpha, rtrans, rec_weight)

   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat, mol, self, wsc, dtrans, rtrans, rec_weight, alpha, npair) &
   !$omp private(iat, jat, izp, jzp, img, gam, wsw, vec, dtmp, rtmp, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      if (iat == jat) then
         vec = 0.0_wp
         call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp, skip_zero=.true.)
         call get_amat_rec_3d(vec, rtrans, rec_weight, rtmp)
         tmp = dtmp + rtmp + self%eta(izp) + sqrt2pi / self%rad(izp) - 2.0_wp * alpha / sqrtpi
      else
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         tmp = 0.0_wp
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp, skip_zero=.false.)
            call get_amat_rec_3d(vec, rtrans, rec_weight, rtmp)
            tmp = tmp + (dtmp + rtmp) * wsw
         end do
      end if
      amat(jat, iat) = tmp
      amat(iat, jat) = tmp
   end do
   !$omp end parallel do

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
end subroutine get_amat_3d

subroutine get_amat_dir_3d(rij, gam, alp, trans, amat, skip_zero)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat
   logical, intent(in), optional :: skip_zero

   integer :: itr
   real(wp) :: vec(3), r1, r2, tmp, kmax
   logical :: exclude_self

   amat = 0.0_wp
   exclude_self = .true.
   if (present(skip_zero)) exclude_self = skip_zero
   kmax = max(gam, alp)

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (exclude_self .and. r1 < eps) cycle
      if (kmax * r1 < 1.0e-2_wp) then
         r2 = r1 * r1
         tmp = (2.0_wp / sqrtpi) * ((gam - alp) &
            & - ((gam**3 - alp**3) / 3.0_wp) * r2 &
            & + (0.1_wp * (gam**5 - alp**5)) * r2 * r2 &
            & - ((gam**7 - alp**7) / 42.0_wp) * r2 * r2 * r2)
      else
         tmp = erf(gam * r1) / r1 - erf(alp * r1) / r1
      end if
      amat = amat + tmp
   end do

end subroutine get_amat_dir_3d

subroutine get_amat_rec_3d(rij, trans, weight, amat)
   real(wp), intent(in) :: rij(3), trans(:, :), weight(:)
   real(wp), intent(out) :: amat
   integer :: itr

   amat = 0.0_wp
   do itr = 1, size(trans, 2)
      amat = amat + weight(itr) * cos(dot_product(rij, trans(:, itr)))
   end do
end subroutine get_amat_rec_3d

!> Evaluate 2D Ewald Coulomb matrix for a slab
subroutine get_amat_2d(self, mol, wsc, alpha, amat, error)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(out) :: amat(:, :)
   type(error_type), allocatable, intent(out), optional :: error

   integer :: iat, jat, izp, jzp, img, ipair, npair, ndim, p_axes(3)
   real(wp) :: vec(3), gam, gam_min, max_rij, wsw, dtmp, rtmp, tmp, area, n_hat(3), a1(3), a2(3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   amat = 0.0_wp
   call get_periodic_axes(mol%periodic, ndim, p_axes)
   a1 = mol%lattice(:, p_axes(1))
   a2 = mol%lattice(:, p_axes(2))
   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
   max_rij = get_max_wsc_disp(mol, wsc)
   call get_ewald_dir_trans_2d(a1, a2, alpha, gam_min, dtrans, max_rij=max_rij, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_ewald_rec_trans_2d(a1, a2, alpha, rtrans, area, n_hat, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if

   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat, mol, self, wsc, dtrans, rtrans, alpha, area, n_hat, npair) &
   !$omp private(iat, jat, izp, jzp, img, gam, wsw, vec, dtmp, rtmp, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      if (iat == jat) then
         vec = 0.0_wp
         call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp, skip_zero=.true.)
         call get_amat_rec_2d(vec, alpha, area, n_hat, rtrans, rtmp)
         tmp = dtmp + rtmp + self%eta(izp) + sqrt2pi / self%rad(izp) - 2.0_wp * alpha / sqrtpi
      else
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         tmp = 0.0_wp
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp, skip_zero=.false.)
            call get_amat_rec_2d(vec, alpha, area, n_hat, rtrans, rtmp)
            tmp = tmp + (dtmp + rtmp) * wsw
         end do
      end if
      amat(jat, iat) = tmp
      amat(iat, jat) = tmp
   end do
   !$omp end parallel do

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
end subroutine get_amat_2d

!> Reciprocal-space and G=0 sheet contribution for 2D Ewald
pure subroutine get_amat_rec_2d(rij, alpha, area, n_hat, trans, amat)
   real(wp), intent(in) :: rij(3), alpha, area, n_hat(3), trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: z, g_vec(3), g_norm, gv, ep, em

   z = dot_product(rij, n_hat)
   amat = 0.0_wp
   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g_norm = norm2(g_vec)
      if (g_norm < eps) cycle
      gv = dot_product(rij, g_vec)
      ep = erfc_exp_2d(g_norm, +z, alpha)
      em = erfc_exp_2d(g_norm, -z, alpha)
      amat = amat + (pi / (area * g_norm)) * cos(gv) * (ep + em)
   end do
   amat = amat - (2.0_wp * pi / area) &
      & * (z * erf(alpha * z) + exp(-(alpha * z)**2) / (alpha * sqrtpi))
end subroutine get_amat_rec_2d

subroutine get_amat_1d(self, mol, wsc, alpha, amat, error)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(out) :: amat(:, :)
   type(error_type), allocatable, intent(out), optional :: error

   integer :: iat, jat, izp, jzp, img, ipair, npair, idir
   real(wp) :: vec(3), gam, gam_min, max_rij, wsw, dtmp, rtmp, tmp
   real(wp) :: a1(3), len1, u_hat(3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   amat = 0.0_wp
   do idir = 1, 3
      if (mol%periodic(idir)) then
         a1 = mol%lattice(:, idir)
         exit
      end if
   end do
   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
   max_rij = get_max_wsc_disp(mol, wsc)
   call get_ewald_dir_trans_1d(a1, alpha, gam_min, dtrans, max_rij=max_rij, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_ewald_rec_trans_1d(a1, alpha, rtrans, len1, u_hat, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if

   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat, mol, self, wsc, dtrans, rtrans, alpha, len1, u_hat, npair) &
   !$omp private(iat, jat, izp, jzp, img, gam, wsw, vec, dtmp, rtmp, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      if (iat == jat) then
         vec = 0.0_wp
         call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp, skip_zero=.true.)
         call get_amat_rec_1d(vec, alpha, len1, u_hat, rtrans, rtmp)
         tmp = dtmp + rtmp + self%eta(izp) + sqrt2pi / self%rad(izp) - 2.0_wp * alpha / sqrtpi
      else
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         tmp = 0.0_wp
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp, skip_zero=.false.)
            call get_amat_rec_1d(vec, alpha, len1, u_hat, rtrans, rtmp)
            tmp = tmp + (dtmp + rtmp) * wsw
         end do
      end if
      amat(jat, iat) = tmp
      amat(iat, jat) = tmp
   end do
   !$omp end parallel do

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
end subroutine get_amat_1d

subroutine get_amat_rec_1d(rij, alpha, len1, u_hat, trans, amat)
   real(wp), intent(in) :: rij(3), alpha, len1, u_hat(3), trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: r_perp(3), rho2, alp2, x, y, g_vec(3), g2, gv
   real(wp) :: k0, k1, k2, f0, f1, f2

   amat = 0.0_wp
   r_perp = rij - dot_product(rij, u_hat) * u_hat
   rho2 = dot_product(r_perp, r_perp)
   alp2 = alpha * alpha
   y = alp2 * rho2

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g2 = dot_product(g_vec, g_vec)
      if (g2 < eps) cycle
      x = 0.25_wp * g2 / alp2
      gv = dot_product(rij, g_vec)
      call fast_k012_1d(x, y, k0, k1, k2)
      amat = amat + cos(gv) * k0 / len1
   end do
   call fast_f012_1d(y, f0, f1, f2)
   amat = amat + (-f0 + 2.0_wp * log(alpha)) / len1
end subroutine get_amat_rec_1d

subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   select case (count(mol%periodic))
   case (3)
      call get_damat_3d(self, mol, ptr%wsc, ptr%alpha, qvec, dadr, dadL, atrace, cache%error)
   case (2)
      call get_damat_2d(self, mol, ptr%wsc, ptr%alpha, qvec, dadr, dadL, atrace, cache%error)
   case (1)
      call get_damat_1d(self, mol, ptr%wsc, ptr%alpha, qvec, dadr, dadL, atrace, cache%error)
   case default
      call get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
   end select
end subroutine get_coulomb_derivs

subroutine get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

   integer :: iat, jat, izp, jzp, ipair, npair
   real(wp) :: vec(3), r2, gam, arg, dtmp, dG(3), dS(3, 3)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   npair = mol%nat*(mol%nat+1)/2
   ! Pair derivatives have unique owners. Only atom-wise traces and strain
   ! derivatives need reduction, with O(N) storage per thread.
   !$omp parallel default(none) shared(dadr, mol, self, qvec, npair) &
   !$omp private(iat, jat, izp, jzp, ipair, gam, r2, vec, dG, dS, dtmp, arg) &
   !$omp reduction(+:atrace, dadL)
   !$omp do schedule(static)
   do iat = 1, size(dadr, 3)
      dadr(:, :, iat) = 0.0_wp
   end do
   !$omp end do
   !$omp do schedule(runtime)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      if (iat == jat) cycle
      izp = mol%id(iat)
      jzp = mol%id(jat)
      vec = mol%xyz(:, jat) - mol%xyz(:, iat)
      r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      arg = gam * gam * r2
      if (arg < 1.0e-4_wp) then
         dtmp = (2.0_wp * gam**3 / sqrtpi) * (-2.0_wp / 3.0_wp &
            & + arg * (0.4_wp - arg * (1.0_wp / 7.0_wp - arg / 27.0_wp)))
      else
         dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) / (r2 * sqrt(r2))
      end if
      dG = dtmp * vec
      dS = spread(dG, 1, 3) * spread(vec, 2, 3)
      atrace(:, iat) = atrace(:, iat) - dG * qvec(jat)
      atrace(:, jat) = atrace(:, jat) + dG * qvec(iat)
      dadr(:, iat, jat) = -dG * qvec(iat)
      dadr(:, jat, iat) = +dG * qvec(jat)
      dadL(:, :, jat) = dadL(:, :, jat) + dS * qvec(iat)
      dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(jat)
   end do
   !$omp end do
   !$omp end parallel
end subroutine get_damat_0d

subroutine get_damat_3d(self, mol, wsc, alpha, qvec, dadr, dadL, atrace, error)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha, qvec(:)
   real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   type(error_type), allocatable, intent(out), optional :: error

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: gam, gam_min, max_rij, wsw, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rec_weight(:), rec_strain(:, :, :)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
   max_rij = get_max_wsc_disp(mol, wsc)
   call get_ewald_dir_trans_3d(mol%lattice, alpha, gam_min, dtrans, max_rij=max_rij, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_ewald_rec_trans_3d(mol%lattice, alpha, rtrans, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_rec_weights(abs(matdet_3x3(mol%lattice)), alpha, rtrans, rec_weight, rec_strain)

   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel default(none) &
   !$omp shared(mol, self, wsc, alpha, dtrans, rtrans, rec_weight, rec_strain, qvec, dadr, npair) &
   !$omp private(iat, jat, izp, jzp, img, ipair, gam, wsw, vec, dG, dS, dGr, dSr, dGd, dSd) &
   !$omp reduction(+:atrace, dadL)
   !$omp do schedule(static)
   do iat = 1, size(dadr, 3)
      dadr(:, :, iat) = 0.0_wp
   end do
   !$omp end do
   !$omp do schedule(runtime)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      dG(:) = 0.0_wp
      dS(:, :) = 0.0_wp
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      if (iat == jat) then
         vec = 0.0_wp
         call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd, skip_zero=.true.)
         call get_damat_rec_3d(vec, rtrans, rec_weight, rec_strain, dGr, dSr)
         dS = dSd + dSr
         dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(iat)
      else
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd, skip_zero=.false.)
            call get_damat_rec_3d(vec, rtrans, rec_weight, rec_strain, dGr, dSr)
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
         end do
         atrace(:, iat) = atrace(:, iat) - dG * qvec(jat)
         atrace(:, jat) = atrace(:, jat) + dG * qvec(iat)
         dadr(:, iat, jat) = -dG * qvec(iat)
         dadr(:, jat, iat) = +dG * qvec(jat)
         dadL(:, :, jat) = dadL(:, :, jat) + dS * qvec(iat)
         dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(jat)
      end if
   end do
   !$omp end do
   !$omp end parallel
end subroutine get_damat_3d

subroutine get_damat_dir_3d(rij, gam, alp, trans, dg, ds, skip_zero)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)
   logical, intent(in), optional :: skip_zero

   integer :: itr
   real(wp) :: vec(3), r1, r2, gtmp, atmp, gam2, alp2, kmax, dtmp
   logical :: exclude_self

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   exclude_self = .true.
   if (present(skip_zero)) exclude_self = skip_zero
   kmax = max(gam, alp)
   gam2 = gam * gam
   alp2 = alp * alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (exclude_self .and. r1 < eps) cycle
      r2 = r1 * r1
      if (kmax * r1 < 1.0e-2_wp) then
         dtmp = (2.0_wp / sqrtpi) * (-2.0_wp * (gam**3 - alp**3) / 3.0_wp &
            & + 0.4_wp * (gam**5 - alp**5) * r2 &
            & - ((gam**7 - alp**7) / 7.0_wp) * r2 * r2)
      else
         gtmp = +2 * gam * exp(-r2 * gam2) / (sqrtpi * r2) - erf(r1 * gam) / (r2 * r1)
         atmp = -2 * alp * exp(-r2 * alp2) / (sqrtpi * r2) + erf(r1 * alp) / (r2 * r1)
         dtmp = gtmp + atmp
      end if
      dg(:) = dg + dtmp * vec
      ds(:, :) = ds + dtmp * spread(vec, 1, 3) * spread(vec, 2, 3)
   end do

end subroutine get_damat_dir_3d

subroutine get_damat_rec_3d(rij, trans, weight, strain, dg, ds)
   real(wp), intent(in) :: rij(3), trans(:, :), weight(:), strain(:, :, :)
   real(wp), intent(out) :: dg(3), ds(3, 3)
   integer :: itr
   real(wp) :: gv

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp
   do itr = 1, size(trans, 2)
      gv = dot_product(rij, trans(:, itr))
      dg(:) = dg - weight(itr) * sin(gv) * trans(:, itr)
      ds(:, :) = ds + cos(gv) * strain(:, :, itr)
   end do
end subroutine get_damat_rec_3d

subroutine get_damat_2d(self, mol, wsc, alpha, qvec, dadr, dadL, atrace, error)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha, qvec(:)
   real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   type(error_type), allocatable, intent(out), optional :: error

   integer :: iat, jat, izp, jzp, img, ipair, npair, idir, ip
   real(wp) :: gam, gam_min, max_rij, wsw, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp) :: a1(3), a2(3), area, n_hat(3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   ip = 0
   do idir = 1, 3
      if (mol%periodic(idir)) then
         ip = ip + 1
         if (ip == 1) a1 = mol%lattice(:, idir)
         if (ip == 2) a2 = mol%lattice(:, idir)
      end if
   end do
   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
   max_rij = get_max_wsc_disp(mol, wsc)
   call get_ewald_dir_trans_2d(a1, a2, alpha, gam_min, dtrans, max_rij=max_rij, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_ewald_rec_trans_2d(a1, a2, alpha, rtrans, area, n_hat, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if

   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel default(none) &
   !$omp shared(mol, self, wsc, alpha, area, n_hat, dtrans, rtrans, qvec, dadr, npair) &
   !$omp private(iat, jat, izp, jzp, img, ipair, gam, wsw, vec, dG, dS, dGr, dSr, dGd, dSd) &
   !$omp reduction(+:atrace, dadL)
   !$omp do schedule(static)
   do iat = 1, size(dadr, 3)
      dadr(:, :, iat) = 0.0_wp
   end do
   !$omp end do
   !$omp do schedule(runtime)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      dG(:) = 0.0_wp
      dS(:, :) = 0.0_wp
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      if (iat == jat) then
         vec = 0.0_wp
         call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd, skip_zero=.true.)
         call get_damat_rec_2d(vec, alpha, area, n_hat, rtrans, dGr, dSr)
         dS = dSd + dSr
         dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(iat)
      else
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd, skip_zero=.false.)
            call get_damat_rec_2d(vec, alpha, area, n_hat, rtrans, dGr, dSr)
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
         end do
         atrace(:, iat) = atrace(:, iat) - dG * qvec(jat)
         atrace(:, jat) = atrace(:, jat) + dG * qvec(iat)
         dadr(:, iat, jat) = -dG * qvec(iat)
         dadr(:, jat, iat) = +dG * qvec(jat)
         dadL(:, :, jat) = dadL(:, :, jat) + dS * qvec(iat)
         dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(jat)
      end if
   end do
   !$omp end do
   !$omp end parallel
end subroutine get_damat_2d

subroutine get_damat_rec_2d(rij, alpha, area, n_hat, trans, dg, ds)
   real(wp), intent(in) :: rij(3), alpha, area, n_hat(3), trans(:, :)
   real(wp), intent(out) :: dg(3), ds(3, 3)

   integer :: itr
   real(wp) :: z, r_perp(3), g_vec(3), g_norm, gv, cos_gv, sin_gv
   real(wp) :: ep, em, f_gz, df_dg, df_dz, pref, phi_0, dphi0_dz
   real(wp) :: i_par(3, 3), nn(3, 3), outer_gg(3, 3), outer_rg(3, 3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   z = dot_product(rij, n_hat)
   r_perp = z * n_hat
   nn = spread(n_hat, 1, 3) * spread(n_hat, 2, 3)
   i_par = unity - nn

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g_norm = norm2(g_vec)
      if (g_norm < eps) cycle
      gv = dot_product(rij, g_vec)
      cos_gv = cos(gv)
      sin_gv = sin(gv)
      ep = erfc_exp_2d(g_norm, +z, alpha)
      em = erfc_exp_2d(g_norm, -z, alpha)
      f_gz = ep + em
      df_dg = z * (ep - em) - (2.0_wp / (alpha * sqrtpi)) &
         & * exp(-g_norm**2 / (4.0_wp * alpha**2) - alpha**2 * z**2)
      df_dz = g_norm * (ep - em)
      pref = pi / (area * g_norm)

      dg(:) = dg + pref * (-sin_gv * f_gz * g_vec + cos_gv * df_dz * n_hat)

      outer_gg = spread(g_vec, 1, 3) * spread(g_vec, 2, 3)
      outer_rg = spread(r_perp, 1, 3) * spread(g_vec, 2, 3) &
         & + spread(g_vec, 1, 3) * spread(r_perp, 2, 3)
      ds(:, :) = ds + pref * cos_gv * (-f_gz * i_par &
         & + ((f_gz / g_norm - df_dg) / g_norm) * outer_gg + z * df_dz * nn) &
         & - pref * sin_gv * f_gz * outer_rg
   end do

   phi_0 = -(2.0_wp * pi / area) &
      & * (z * erf(alpha * z) + exp(-(alpha * z)**2) / (alpha * sqrtpi))
   dphi0_dz = -(2.0_wp * pi / area) * erf(alpha * z)
   dg(:) = dg + dphi0_dz * n_hat
   ds(:, :) = ds - phi_0 * i_par + dphi0_dz * z * nn
end subroutine get_damat_rec_2d

subroutine get_damat_1d(self, mol, wsc, alpha, qvec, dadr, dadL, atrace, error)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: alpha, qvec(:)
   real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   type(error_type), allocatable, intent(out), optional :: error

   integer :: iat, jat, izp, jzp, img, ipair, npair, idir
   real(wp) :: gam, gam_min, max_rij, wsw, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp) :: a1(3), len1, u_hat(3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   do idir = 1, 3
      if (mol%periodic(idir)) then
         a1 = mol%lattice(:, idir)
         exit
      end if
   end do
   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
   max_rij = get_max_wsc_disp(mol, wsc)
   call get_ewald_dir_trans_1d(a1, alpha, gam_min, dtrans, max_rij=max_rij, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if
   call get_ewald_rec_trans_1d(a1, alpha, rtrans, len1, u_hat, error=error)
   if (present(error)) then
      if (allocated(error)) return
   end if

   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel default(none) &
   !$omp shared(mol, self, wsc, alpha, len1, u_hat, dtrans, rtrans, qvec, dadr, npair) &
   !$omp private(iat, jat, izp, jzp, img, ipair, gam, wsw, vec, dG, dS, dGr, dSr, dGd, dSd) &
   !$omp reduction(+:atrace, dadL)
   !$omp do schedule(static)
   do iat = 1, size(dadr, 3)
      dadr(:, :, iat) = 0.0_wp
   end do
   !$omp end do
   !$omp do schedule(runtime)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      dG(:) = 0.0_wp
      dS(:, :) = 0.0_wp
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      if (iat == jat) then
         vec = 0.0_wp
         call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd, skip_zero=.true.)
         call get_damat_rec_1d(vec, alpha, len1, u_hat, rtrans, dGr, dSr)
         dS = dSd + dSr
         dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(iat)
      else
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd, skip_zero=.false.)
            call get_damat_rec_1d(vec, alpha, len1, u_hat, rtrans, dGr, dSr)
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
         end do
         atrace(:, iat) = atrace(:, iat) - dG * qvec(jat)
         atrace(:, jat) = atrace(:, jat) + dG * qvec(iat)
         dadr(:, iat, jat) = -dG * qvec(iat)
         dadr(:, jat, iat) = +dG * qvec(jat)
         dadL(:, :, jat) = dadL(:, :, jat) + dS * qvec(iat)
         dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(jat)
      end if
   end do
   !$omp end do
   !$omp end parallel
end subroutine get_damat_1d

subroutine get_damat_rec_1d(rij, alpha, len1, u_hat, trans, dg, ds)
   real(wp), intent(in) :: rij(3), alpha, len1, u_hat(3), trans(:, :)
   real(wp), intent(out) :: dg(3), ds(3, 3)

   integer :: itr
   real(wp) :: r_perp(3), rho2, alp2, x, y, g_vec(3), g2, gv, cos_gv, sin_gv
   real(wp) :: k0, k1, k2, f0, f1, f2, exp_xy
   real(wp) :: p_par(3, 3), outer_rr(3, 3), outer_rg(3, 3)

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   r_perp = rij - dot_product(rij, u_hat) * u_hat
   rho2 = dot_product(r_perp, r_perp)
   alp2 = alpha * alpha
   y = alp2 * rho2
   p_par = spread(u_hat, 1, 3) * spread(u_hat, 2, 3)
   outer_rr = spread(r_perp, 1, 3) * spread(r_perp, 2, 3)

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g2 = dot_product(g_vec, g_vec)
      if (g2 < eps) cycle
      x = 0.25_wp * g2 / alp2
      gv = dot_product(rij, g_vec)
      cos_gv = cos(gv)
      sin_gv = sin(gv)
      call fast_k012_1d(x, y, k0, k1, k2)
      exp_xy = exp(-x - y)

      dg(:) = dg - (sin_gv * k0 / len1) * g_vec &
         & - (2.0_wp * alp2 * cos_gv * k1 / len1) * r_perp

      outer_rg = spread(g_vec, 1, 3) * spread(r_perp, 2, 3) &
         & + spread(r_perp, 1, 3) * spread(g_vec, 2, 3)
      ds(:, :) = ds + (cos_gv * (2.0_wp * exp_xy + 2.0_wp * y * k1 - k0) / len1) * p_par &
         & - (2.0_wp * alp2 * cos_gv * k1 / len1) * outer_rr &
         & - (sin_gv * k0 / len1) * outer_rg
   end do

   call fast_f012_1d(y, f0, f1, f2)
   dg(:) = dg - (2.0_wp * alp2 * f1 / len1) * r_perp
   ds(:, :) = ds + ((f0 - 2.0_wp * log(alpha)) / len1) * p_par &
      & - (2.0_wp * alp2 * f1 / len1) * outer_rr
end subroutine get_damat_rec_1d

!> Reciprocal-space factors independent of the atom pair.
!> Prepare once for each matrix/response evaluation, including strain when needed.
subroutine get_rec_weights(vol, alpha, trans, weight, strain)
   real(wp), intent(in) :: vol, alpha, trans(:, :)
   real(wp), allocatable, intent(out) :: weight(:)
   real(wp), allocatable, intent(out), optional :: strain(:, :, :)
   integer :: itr
   real(wp) :: g2, fac, vec(3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])

   fac = 4 * pi / vol
   allocate(weight(size(trans, 2)), source=0.0_wp)
   if (present(strain)) allocate(strain(3, 3, size(trans, 2)), source=0.0_wp)
   do itr = 1, size(trans, 2)
      vec = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      weight(itr) = fac * exp(-0.25_wp * g2 / alpha**2) / g2
      if (present(strain)) then
         strain(:, :, itr) = weight(itr) * ((2.0_wp / g2 + 0.5_wp / alpha**2) &
            & * spread(vec, 1, 3) * spread(vec, 2, 3) - unity)
      end if
   end do
end subroutine get_rec_weights

!> Inspect cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the cache
   type(cache_container), target, intent(inout) :: cache
   !> Reference to the cache
   type(eeq_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not. allocated(cache%raw)) then
      block
         type(eeq_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the cache
   type(cache_container), target, intent(inout) :: cache
   !> Reference to the cache
   type(eeq_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(eeq_cache)
      ptr => target
   end select
end subroutine view

end module multicharge_model_eeq
