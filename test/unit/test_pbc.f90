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

module test_pbc
   use mctc_cutoff, only: get_lattice_points
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge, only: get_charges, add_charge_hessian
   use multicharge_ewald, only: get_ewald_dir_trans_1d, get_ewald_rec_trans_1d, &
      & get_ewald_dir_trans_2d, get_ewald_rec_trans_2d, &
      & get_ewald_dir_trans_3d, get_ewald_rec_trans_3d, &
      & fast_f012_1d, fast_k012_1d, erfc_exp_2d
   use multicharge_model, only: mchrg_model_type
   use multicharge_model_cache, only: cache_container
   use multicharge_model_eeq, only: eeq_model, eeq_cache
   use multicharge_model_eeqbc, only: eeqbc_model, eeqbc_cache
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   implicit none
   private

   public :: collect_pbc

   real(wp), parameter :: thr = 1000 * epsilon(1.0_wp)
   real(wp), parameter :: thr2 = 10 * sqrt(epsilon(1.0_wp))

contains

!> Collect all exported unit tests
subroutine collect_pbc(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("eeq-charges-cyanamide", test_eeq_q_cyanamide), &
      & new_unittest("eeq-energy-formamide", test_eeq_e_formamide), &
      & new_unittest("eeq-dbdr-co2", test_eeq_dbdr_co2), &
      & new_unittest("eeq-dbdL-co2", test_eeq_dbdL_co2), &
      & new_unittest("eeq-dadr-ice", test_eeq_dadr_ice), &
      & new_unittest("eeq-dadL-ice", test_eeq_dadL_ice), &
      & new_unittest("eeq-gradient-co2", test_eeq_g_co2), &
      & new_unittest("eeq-sigma-ice", test_eeq_s_ice), &
      & new_unittest("eeq-dqdr-urea", test_eeq_dqdr_urea), &
      & new_unittest("eeq-dqdL-oxacb", test_eeq_dqdL_oxacb), &
      & new_unittest("eeqbc-dbdr-co2", test_eeqbc_dbdr_co2), &
      & new_unittest("eeqbc-dbdL-co2", test_eeqbc_dbdL_co2), &
      & new_unittest("eeqbc-dadr-ice", test_eeqbc_dadr_ice), &
      & new_unittest("eeqbc-response-znoh", test_eeqbc_response_znoh), &
      & new_unittest("eeqbc-dadL-ice", test_eeqbc_dadL_ice), &
      & new_unittest("eeqbc-gradient-co2", test_eeqbc_g_co2), &
      & new_unittest("eeqbc-sigma-ice", test_eeqbc_s_ice), &
      & new_unittest("eeqbc-dqdr-urea", test_eeqbc_dqdr_urea), &
      & new_unittest("eeqbc-dqdL-oxacb", test_eeqbc_dqdL_oxacb), &
      & new_unittest("eeqbc-1d-all-masks", test_eeqbc_1d_all_masks), &
      & new_unittest("eeqbc-2d-all-masks", test_eeqbc_2d_all_masks), &
      & new_unittest("eeq-1d-all-masks", test_eeq_1d_all_masks), &
      & new_unittest("eeq-2d-all-masks", test_eeq_2d_all_masks), &
      & new_unittest("ewald-convergence-1d-2d", test_ewald_convergence_1d_2d), &
      & new_unittest("pbc-wsc-invariance-1d-2d", test_pbc_wsc_invariance_1d_2d), &
      & new_unittest("eeqbc-convergence-1d-2d", test_eeqbc_convergence_1d_2d), &
      & new_unittest("eeqbc-water-helium-358", test_eeqbc_water_helium_358), &
      & new_unittest("eeq-3d-slab-285", test_eeq_3d_slab_285), &
      & new_unittest("special-functions-1d-2d", test_special_functions_1d_2d), &
      & new_unittest("degenerate-lattice-error", test_degenerate_lattice_error), &
      & new_unittest("structure-and-contract-validation", test_structure_and_contract_validation), &
      & new_unittest("charged-lowd-contracts", test_charged_lowd_contracts) &
      & ]

end subroutine collect_pbc

!> Compact skew cell probes nonlocal widths and their periodic self images.
subroutine test_eeqbc_response_znoh(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call new(mol, [30, 8, 1], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.1_wp, 0.2_wp, &
      & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3]), charge=1.0_wp, &
      & lattice=reshape([5.1_wp, 0.0_wp, 0.0_wp, 0.4_wp, 6.2_wp, 0.0_wp, &
      & 0.2_wp, 0.7_wp, 7.1_wp], [3, 3]), periodic=[.true., .true., .true.])
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)
   if (allocated(error)) return
   call test_dadL(error, mol, model)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)
end subroutine test_eeqbc_response_znoh

subroutine gen_test(error, mol, model, qref, eref)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Reference charges
   real(wp), intent(in), optional :: qref(:)

   !> Reference energies
   real(wp), intent(in), optional :: eref(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), allocatable :: cn(:), qloc(:), trans(:, :)
   real(wp), allocatable :: energy(:)
   real(wp), allocatable :: qvec(:)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), qloc(mol%nat))

   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)

   if (present(eref)) then
      allocate(energy(mol%nat))
      energy(:) = 0.0_wp
   end if
   if (present(qref)) then
      allocate(qvec(mol%nat))
   end if

   call model%solve(mol, error, cn, qloc, energy=energy, qvec=qvec)
   if (allocated(error)) return

   if (present(qref)) then
      if (any(abs(qvec - qref) > thr)) then
         call test_failed(error, "Partial charges do not match")
         print"(a)", "Charges:"
         print"(3es21.14)", qvec
         print'("---")'
         print"(3es21.14)", qref
         print'("---")'
         print"(3es21.14)", qvec - qref
      end if
   end if
   if (allocated(error)) return

   if (present(eref)) then
      if (any(abs(energy - eref) > thr)) then
         call test_failed(error, "Energies do not match")
         print"(a)", "Energy:"
         print"(3es21.14)", energy
         print'("---")'
         print"(3es21.14)", eref
         print'("---")'
         print"(3es21.14)", energy - eref
      end if
   end if

end subroutine gen_test

subroutine test_numgrad(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: numgrad(:, :)
   real(wp) :: er, el

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & energy(mol%nat), gradient(3, mol%nat), sigma(3, 3), numgrad(3, mol%nat))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         energy(:) = 0.0_wp
         er = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp * (er - el) / step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   energy(:) = 0.0_wp
   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, energy, gradient, sigma)
   if (allocated(error)) return

   if (any(abs(gradient(:, :) - numgrad(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
      print"(a)", "gradient:"
      print"(3es21.14)", gradient
      print"(a)", "numgrad:"
      print"(3es21.14)", numgrad
      print"(a)", "diff:"
      print"(3es21.14)", gradient - numgrad
   end if

end subroutine test_numgrad

subroutine test_numsigma(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic, jc
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :)
   real(wp) :: er, el, eps(3, 3), numsigma(3, 3), sigma(3, 3), lattice(3, 3)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & energy(mol%nat), gradient(3, mol%nat), xyz(3, mol%nat))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         lattr(:, :) = trans
         numsigma(jc, ic) = 0.5_wp * (er - el) / step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   energy(:) = 0.0_wp
   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, energy, gradient, sigma)
   if (allocated(error)) return

   if (any(abs(sigma(:, :) - numsigma(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
      print"(a)", "sigma:"
      print"(3es21.14)", sigma
      print"(a)", "numgrad:"
      print"(3es21.14)", numsigma
      print"(a)", "diff:"
      print"(3es21.14)", sigma(:, :) - numsigma(:, :)
   end if

end subroutine test_numsigma

subroutine test_dbdr(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dbdr(:, :, :), dbdL(:, :, :)
   real(wp), allocatable :: numgrad(:, :, :), xvecr(:), xvecl(:)
   type(cache_container), allocatable :: cache
   allocate(cache)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & xvecr(mol%nat + 1), xvecl(mol%nat + 1), numgrad(3, mol%nat, mol%nat + 1), &
      & dbdr(3, mol%nat, mol%nat + 1), dbdL(3, 3, mol%nat + 1))

   numgrad = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         ! Right-hand side
         xvecr(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if
         call model%get_xvec(mol, cache, xvecr)

         ! Left-hand side
         xvecl(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if
         call model%get_xvec(mol, cache, xvecl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat, :) = 0.5_wp * (xvecr(:) - xvecl(:)) / step
      end do
   end do lp
   if (allocated(error)) return

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   if (allocated(cache%error)) then
      call move_alloc(cache%error, error)
      return
   end if
   call model%get_xvec(mol, cache, xvecl) ! need to call this for xtmp in cache (eeqbc)
   call model%get_xvec_derivs(mol, cache, dbdr, dbdL)

   if (any(abs(dbdr(:, :, :) - numgrad(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of the b vector does not match")
      print"(a)", "dbdr:"
      print"(3es21.14)", dbdr
      print"(a)", "numgrad:"
      print"(3es21.14)", numgrad
      print"(a)", "diff:"
      print"(3es21.14)", dbdr - numgrad
   end if

end subroutine test_dbdr

subroutine test_dbdL(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic, jc
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dbdr(:, :, :), dbdL(:, :, :)
   real(wp), allocatable :: numsigma(:, :, :), xvecr(:), xvecl(:)
   real(wp), allocatable :: xyz(:, :), lattr(:, :), trans(:, :)
   real(wp) :: lattice(3, 3)
   real(wp) :: eps(3, 3)
   type(cache_container), allocatable :: cache
   allocate(cache)

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & xvecr(mol%nat + 1), xvecl(mol%nat + 1), numsigma(3, 3, mol%nat + 1), &
      & dbdr(3, mol%nat, mol%nat + 1), dbdL(3, 3, mol%nat + 1), xyz(3, mol%nat))

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   numsigma = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         ! Right-hand side
         xvecr(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, cache, cn, qloc)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if
         call model%get_xvec(mol, cache, xvecr)

         ! Left-hand side
         xvecl(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, cache, cn, qloc)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if
         call model%get_xvec(mol, cache, xvecl)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         lattr(:, :) = trans
         do iat = 1, mol%nat
            numsigma(jc, ic, iat) = 0.5_wp * (xvecr(iat) - xvecl(iat)) / step
         end do
      end do
   end do lp
   if (allocated(error)) return

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   if (allocated(cache%error)) then
      call move_alloc(cache%error, error)
      return
   end if
   call model%get_xvec(mol, cache, xvecl) ! need to call this for xtmp in cache (eeqbc)
   call model%get_xvec_derivs(mol, cache, dbdr, dbdL)

   if (any(abs(dbdL(:, :, :) - numsigma(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of the b vector does not match")
      print"(a)", "dbdL:"
      print"(3es21.14)", dbdL
      print"(a)", "numsigma:"
      print"(3es21.14)", numsigma
      print"(a)", "diff:"
      print"(3es21.14)", dbdL - numsigma
   end if

end subroutine test_dbdL

subroutine test_dadr(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic, jat, kat
   real(wp) :: thr2_local
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: qloc(:)
   real(wp), allocatable :: trans(:, :)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: qvec(:), numgrad(:, :, :), amatr(:, :), amatl(:, :), numtrace(:, :)
   type(cache_container), allocatable :: cache
   allocate(cache)

   allocate(cn(mol%nat), qloc(mol%nat), amatr(mol%nat + 1, mol%nat + 1), amatl(mol%nat + 1, mol%nat + 1), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), dqlocdr(3, mol%nat, mol%nat), &
      & dqlocdL(3, 3, mol%nat), dadr(3, mol%nat, mol%nat + 1), dadL(3, 3, mol%nat + 1), &
      & atrace(3, mol%nat), numtrace(3, mol%nat), numgrad(3, mol%nat, mol%nat + 1), qvec(mol%nat))

   ! Set tolerance higher if testing eeqbc model
   select type(model)
   type is(eeqbc_model)
      thr2_local = 3.0_wp * thr2
   class default
      thr2_local = thr2
   end select

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   ! Obtain the vector of charges
   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)
   call model%solve(mol, error, cn, qloc, qvec=qvec)
   if (allocated(error)) return

   numgrad = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         ! Right-hand side
         amatr(:, :) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if
         call model%get_coulomb_matrix(mol, cache, amatr)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if

         ! Left-hand side
         amatl(:, :) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if
         call model%get_coulomb_matrix(mol, cache, amatl)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         do kat = 1, mol%nat
            do jat = 1, mol%nat
               ! Numerical gradient of the A matrix
               numgrad(ic, iat, kat) = 0.5_wp * qvec(jat) * (amatr(kat, jat) - amatl(kat, jat)) / step &
                  & + numgrad(ic, iat, kat)
            end do
         end do
      end do
   end do lp
   if (allocated(error)) return

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   if (allocated(cache%error)) then
      call move_alloc(cache%error, error)
      return
   end if
   call model%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)
   if (allocated(cache%error)) then
      call move_alloc(cache%error, error)
      return
   end if

   ! Add trace of the A matrix
   do iat = 1, mol%nat
      dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
   end do

   ! higher tolerance for numerical gradient
   if (any(abs(dadr(:, :, :) - numgrad(:, :, :)) > thr2_local)) then
      call test_failed(error, "Derivative of the A matrix does not match")
      print"(a)", "dadr:"
      print"(3es21.14)", dadr
      print"(a)", "numgrad:"
      print"(3es21.14)", numgrad
      print"(a)", "diff:"
      print"(3es21.14)", dadr - numgrad
   end if

end subroutine test_dadr

subroutine test_dadL(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic, jc, iat
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: xyz(:, :), lattr(:, :), trans(:, :)
   real(wp), allocatable :: qvec(:), numsigma(:, :, :), amatr(:, :), amatl(:, :)
   real(wp) :: lattice(3, 3)
   real(wp) :: eps(3, 3)
   type(cache_container), allocatable :: cache
   allocate(cache)

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & amatr(mol%nat + 1, mol%nat + 1), amatl(mol%nat + 1, mol%nat + 1), &
      & dadr(3, mol%nat, mol%nat + 1), dadL(3, 3, mol%nat + 1), atrace(3, mol%nat), &
      & numsigma(3, 3, mol%nat + 1), qvec(mol%nat), xyz(3, mol%nat))

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)
   call model%solve(mol, error, cn, qloc, qvec=qvec)
   if (allocated(error)) return

   numsigma = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         amatr(:, :) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, cache, cn, qloc)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if
         call model%get_coulomb_matrix(mol, cache, amatr)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if

         amatl(:, :) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, cache, cn, qloc)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if
         call model%get_coulomb_matrix(mol, cache, amatl)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            exit lp
         end if

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         lattr(:, :) = trans
         do iat = 1, mol%nat
            ! Numerical sigma of the a matrix
            numsigma(jc, ic, :) = 0.5_wp * qvec(iat) * (amatr(iat, :) - amatl(iat, :)) / step + numsigma(jc, ic, :)
         end do
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   if (allocated(cache%error)) then
      call move_alloc(cache%error, error)
      return
   end if

   call model%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)
   if (allocated(cache%error)) then
      call move_alloc(cache%error, error)
      return
   end if

   if (any(abs(dadL(:, :, :) - numsigma(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of the A matrix does not match")
      print"(a)", "dadL:"
      print"(3es21.14)", dadL
      print"(a)", "numsigma:"
      print"(3es21.14)", numsigma
      print"(a)", "diff:"
      print"(3es21.14)", dadL - numsigma
   end if

end subroutine test_dadL

subroutine test_numdqdr(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: ql(:), qr(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & ql(mol%nat), qr(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         qr = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=qr)
         if (allocated(error)) exit lp

         ql = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=ql)
         if (allocated(error)) exit lp

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numdr(ic, iat, :) = 0.5_wp * (qr - ql) / step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call model%solve(mol, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, dqdr=dqdr, dqdL=dqdL)
   if (allocated(error)) return

   if (any(abs(dqdr(:, :, :) - numdr(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of charges does not match")
      print"(a)", "dqdr:"
      print"(3es21.14)", dqdr
      print"(a)", "numdr:"
      print"(3es21.14)", numdr
      print"(a)", "diff:"
      print"(3es21.14)", dqdr - numdr
   end if

end subroutine test_numdqdr

subroutine test_numdqdL(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic, jc
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: qr(:), ql(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :), numdL(:, :, :)
   real(wp) :: eps(3, 3), lattice(3, 3)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & qr(mol%nat), ql(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & xyz(3, mol%nat), numdL(3, 3, mol%nat))

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         qr = 0.0_wp
         ql = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%solve(mol, error, cn, qloc, qvec=qr)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%solve(mol, error, cn, qloc, qvec=ql)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         lattr(:, :) = trans
         numdL(jc, ic, :) = 0.5_wp * (qr - ql) / step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, dqdr=dqdr, dqdL=dqdL)
   if (allocated(error)) return

   if (any(abs(dqdL(:, :, :) - numdL(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of charges does not match")
      print"(a)", "dqdL:"
      print"(3es21.14)", dqdL
      print"(a)", "numdL:"
      print"(3es21.14)", numdL
      print"(a)", "diff:"
      print"(3es21.14)", dqdL - numdL
   end if

end subroutine test_numdqdL

subroutine test_eeq_q_cyanamide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   ! Converged 3D Ewald reference charges (differs by <= 2.3e-6 from legacy rep=[2,2,2]).
   real(wp), parameter :: ref(40) = [&
      & 3.47273696479021E-1_wp, 3.47272705640270E-1_wp, 3.47276492903749E-1_wp, &
      & 3.47266922178364E-1_wp, 3.47288439460215E-1_wp, 3.47296528581096E-1_wp, &
      & 3.47299555433948E-1_wp, 3.47299047247046E-1_wp, 3.45885448217334E-1_wp, &
      & 3.45873984082476E-1_wp, 3.45887979518336E-1_wp, 3.45877189683493E-1_wp, &
      & 3.45902103370702E-1_wp, 3.45901899732038E-1_wp, 3.45900074770712E-1_wp, &
      & 3.45903507923248E-1_wp, 3.58994427905596E-1_wp, 3.58985684782718E-1_wp, &
      & 3.58998551488410E-1_wp, 3.58988652945645E-1_wp, 3.58988114752443E-1_wp, &
      & 3.58995533629608E-1_wp, 3.58994816857993E-1_wp, 3.58999419102342E-1_wp, &
      &-5.86636852376099E-1_wp, -5.86601030236713E-1_wp, -5.86624269538219E-1_wp, &
      &-5.86579192001313E-1_wp, -5.86653588719833E-1_wp, -5.86668553721606E-1_wp, &
      &-5.86666673919381E-1_wp, -5.86674064898699E-1_wp, -4.65527015747556E-1_wp, &
      &-4.65525637708548E-1_wp, -4.65544375954499E-1_wp, -4.65538830945934E-1_wp, &
      &-4.65527785578432E-1_wp, -4.65531681382053E-1_wp, -4.65525608429865E-1_wp, &
      &-4.65525615528052E-1_wp]

   call get_structure(mol, "X23", "cyanamide")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

end subroutine test_eeq_q_cyanamide

subroutine test_eeq_e_formamide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   ! Converged 3D Ewald reference energies (differs by <= 1.35e-4 from legacy rep=[2,2,2]).
   real(wp), parameter :: ref(24) = [&
      & 4.01820317641649E-1_wp, 4.01826970293715E-1_wp, 4.01778618929473E-1_wp, &
      & 4.01801620737465E-1_wp, 4.30133215169715E-1_wp, 4.30106241841601E-1_wp, &
      & 4.30086201257563E-1_wp, 4.30109482177696E-1_wp, 1.90391275024120E-1_wp, &
      & 1.90392904370207E-1_wp, 1.90395221896778E-1_wp, 1.90398527506007E-1_wp, &
      & 3.33798778194695E-1_wp, 3.33779134873330E-1_wp, 3.33806406791091E-1_wp, &
      & 3.33823483957305E-1_wp, -8.35501801552073E-1_wp, -8.35481140201098E-1_wp, &
      &-8.35408885793255E-1_wp, -8.35455735270739E-1_wp, -6.31069087329026E-1_wp, &
      &-6.31049654240391E-1_wp, -6.31075181828644E-1_wp, -6.31071707101098E-1_wp]

   call get_structure(mol, "X23", "formamide")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, eref=ref)

end subroutine test_eeq_e_formamide

subroutine test_eeq_dbdr_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model)

end subroutine test_eeq_dbdr_co2

subroutine test_eeq_dbdL_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdL(error, mol, model)

end subroutine test_eeq_dbdL_co2

subroutine test_eeq_dadr_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)

end subroutine test_eeq_dadr_ice

subroutine test_eeq_dadL_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model)

end subroutine test_eeq_dadL_ice

subroutine test_eeq_g_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)

end subroutine test_eeq_g_co2

subroutine test_eeq_s_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)

end subroutine test_eeq_s_ice

subroutine test_eeq_dqdr_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "urea")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)

end subroutine test_eeq_dqdr_urea

subroutine test_eeq_dqdL_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "oxacb")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model)

end subroutine test_eeq_dqdL_oxacb

subroutine test_eeqbc_dbdr_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model)

end subroutine test_eeqbc_dbdr_co2

subroutine test_eeqbc_dbdL_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdL(error, mol, model)

end subroutine test_eeqbc_dbdL_co2

subroutine test_eeqbc_dadr_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)

end subroutine test_eeqbc_dadr_ice

subroutine test_eeqbc_dadL_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model)

end subroutine test_eeqbc_dadL_ice

subroutine test_eeqbc_g_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)

end subroutine test_eeqbc_g_co2

subroutine test_eeqbc_s_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)

end subroutine test_eeqbc_s_ice

subroutine test_eeqbc_dqdr_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "urea")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)

end subroutine test_eeqbc_dqdr_urea

subroutine test_eeqbc_dqdL_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "oxacb")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model)

end subroutine test_eeqbc_dqdL_oxacb

subroutine test_eeqbc_1d_all_masks(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol, mol_big
   class(mchrg_model_type), allocatable :: model
   real(wp) :: xyz(3, 3), a_wire(3), lat1(3, 3), lat2(3, 3)
   real(wp) :: q1(3), q2(3), q_ref(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3)
   real(wp) :: dqdL1(3, 3, 3), dqdL2(3, 3, 3)
   logical :: masks(3, 3)
   integer :: imask, iq
   real(wp) :: charges(2)

   xyz = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 2.4_wp, 0.3_wp, -0.4_wp, &
      & 1.1_wp, 1.8_wp, 0.7_wp], [3, 3])
   a_wire = [5.3_wp, 0.6_wp, -0.5_wp]
   masks = reshape([ &
      & .true., .false., .false., &
      & .false., .true., .false., &
      & .false., .false., .true.], [3, 3])
   charges = [0.0_wp, 1.0_wp]

   do iq = 1, 2
      do imask = 1, 3
         lat1 = 0.0_wp
         lat1(:, imask) = a_wire
         ! Give nonperiodic columns arbitrary values in lat1 vs 10x larger in lat2
         lat1(:, mod(imask, 3) + 1) = [0.0_wp, 20.0_wp, 0.0_wp]
         lat1(:, mod(imask + 1, 3) + 1) = [0.0_wp, 0.0_wp, 20.0_wp]

         lat2 = 0.0_wp
         lat2(:, imask) = a_wire
         lat2(:, mod(imask, 3) + 1) = [15.0_wp, 200.0_wp, -30.0_wp]
         lat2(:, mod(imask + 1, 3) + 1) = [-40.0_wp, 10.0_wp, 250.0_wp]

         call new(mol, [30, 8, 1], xyz, charge=charges(iq), &
            & lattice=lat1, periodic=masks(:, imask))
         call new(mol_big, [30, 8, 1], xyz, charge=charges(iq), &
            & lattice=lat2, periodic=masks(:, imask))

         call new_eeqbc2025_model(mol, model, error)
         if (allocated(error)) return

         call get_charges(model, mol, error, q1, dqdr1, dqdL1)
         if (allocated(error)) return
         call get_charges(model, mol_big, error, q2, dqdr2, dqdL2)
         if (allocated(error)) return

         if (abs(sum(q1) - charges(iq)) > thr) then
            call test_failed(error, "1D EEQBC total charge constraint violated")
            return
         end if
         if (maxval(abs(q1 - q2)) > thr .or. maxval(abs(dqdr1 - dqdr2)) > thr &
            & .or. maxval(abs(dqdL1 - dqdL2)) > thr) then
            call test_failed(error, "1D EEQBC depends on nonperiodic lattice vectors")
            return
         end if
         if (imask == 1) then
            q_ref = q1
         else if (maxval(abs(q1 - q_ref)) > thr) then
            call test_failed(error, "1D EEQBC charges depend on periodic mask index")
            return
         end if

         ! Full derivative validation across all 1D masks: TFF, FTF (neutral) and FFT (charged)
         if ((iq == 1 .and. imask <= 2) .or. (iq == 2 .and. imask == 3)) then
            call test_dbdr(error, mol, model)
            if (allocated(error)) return
            call test_dbdL(error, mol, model)
            if (allocated(error)) return
            call test_dadr(error, mol, model)
            if (allocated(error)) return
            call test_dadL(error, mol, model)
            if (allocated(error)) return
            call test_numgrad(error, mol, model)
            if (allocated(error)) return
            call test_numsigma(error, mol, model)
            if (allocated(error)) return
            call test_numdqdr(error, mol, model)
            if (allocated(error)) return
            call test_numdqdL(error, mol, model)
            if (allocated(error)) return
         end if
      end do
   end do
end subroutine test_eeqbc_1d_all_masks

subroutine test_eeqbc_2d_all_masks(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol, mol_big
   class(mchrg_model_type), allocatable :: model
   real(wp) :: xyz(3, 3), a1(3), a2(3), lat1(3, 3), lat2(3, 3)
   real(wp) :: q1(3), q2(3), q_ref(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3)
   real(wp) :: dqdL1(3, 3, 3), dqdL2(3, 3, 3)
   logical :: masks(3, 3)
   integer :: imask, iq, p1, p2, np
   real(wp) :: charges(2)

   xyz = reshape([ &
      & 0.0_wp, 0.0_wp, -0.3_wp, &
      & 2.5_wp, 0.2_wp, 0.5_wp, &
      & 0.6_wp, 2.1_wp, 0.4_wp], [3, 3])
   a1 = [5.2_wp, 0.3_wp, 0.1_wp]
   a2 = [0.5_wp, 6.1_wp, -0.2_wp]
   masks = reshape([ &
      & .true., .true., .false., &
      & .true., .false., .true., &
      & .false., .true., .true.], [3, 3])
   charges = [0.0_wp, -1.0_wp]

   do iq = 1, 2
      do imask = 1, 3
         select case (imask)
         case (1); p1 = 1; p2 = 2; np = 3
         case (2); p1 = 1; p2 = 3; np = 2
         case default; p1 = 2; p2 = 3; np = 1
         end select
         lat1 = 0.0_wp
         lat1(:, p1) = a1
         lat1(:, p2) = a2
         lat1(:, np) = [0.0_wp, 0.0_wp, 20.0_wp]

         lat2 = 0.0_wp
         lat2(:, p1) = a1
         lat2(:, p2) = a2
         lat2(:, np) = [35.0_wp, -70.0_wp, 200.0_wp]

         call new(mol, [30, 8, 1], xyz, charge=charges(iq), &
            & lattice=lat1, periodic=masks(:, imask))
         call new(mol_big, [30, 8, 1], xyz, charge=charges(iq), &
            & lattice=lat2, periodic=masks(:, imask))

         call new_eeqbc2025_model(mol, model, error)
         if (allocated(error)) return

         call get_charges(model, mol, error, q1, dqdr1, dqdL1)
         if (allocated(error)) return
         call get_charges(model, mol_big, error, q2, dqdr2, dqdL2)
         if (allocated(error)) return

         if (abs(sum(q1) - charges(iq)) > thr) then
            call test_failed(error, "2D EEQBC total charge constraint violated")
            return
         end if
         if (maxval(abs(q1 - q2)) > thr .or. maxval(abs(dqdr1 - dqdr2)) > thr &
            & .or. maxval(abs(dqdL1 - dqdL2)) > thr) then
            call test_failed(error, "2D EEQBC depends on nonperiodic lattice vector")
            return
         end if
         if (imask == 1) then
            q_ref = q1
         else if (maxval(abs(q1 - q_ref)) > thr) then
            call test_failed(error, "2D EEQBC charges depend on periodic mask index")
            return
         end if

         ! Full derivative validation across all 2D masks: TTF, TFT (neutral) and FTT (charged)
         if ((iq == 1 .and. imask <= 2) .or. (iq == 2 .and. imask == 3)) then
            call test_dbdr(error, mol, model)
            if (allocated(error)) return
            call test_dbdL(error, mol, model)
            if (allocated(error)) return
            call test_dadr(error, mol, model)
            if (allocated(error)) return
            call test_dadL(error, mol, model)
            if (allocated(error)) return
            call test_numgrad(error, mol, model)
            if (allocated(error)) return
            call test_numsigma(error, mol, model)
            if (allocated(error)) return
            call test_numdqdr(error, mol, model)
            if (allocated(error)) return
            call test_numdqdL(error, mol, model)
            if (allocated(error)) return
         end if
      end do
   end do
end subroutine test_eeqbc_2d_all_masks

subroutine test_degenerate_lattice_error(error)
   type(error_type), allocatable, intent(out) :: error
   type(error_type), allocatable :: expected
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp) :: q(2), lat(3, 3)

   ! 1D zero periodic vector
   lat = 0.0_wp
   lat(:, 2) = [0.0_wp, 10.0_wp, 0.0_wp]
   call new(mol, [8, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp, 0.0_wp, 0.0_wp], [3, 2]))
   mol%lattice = lat
   mol%periodic = [.true., .false., .false.]
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call get_charges(model, mol, expected, q)
   if (.not. allocated(expected)) then
      call test_failed(error, "Degenerate 1D lattice was not rejected")
      return
   end if
   deallocate(expected)

   ! 2D collinear periodic vectors
   lat = 0.0_wp
   lat(:, 1) = [5.0_wp, 0.0_wp, 0.0_wp]
   lat(:, 2) = [10.0_wp, 0.0_wp, 0.0_wp]
   lat(:, 3) = [0.0_wp, 0.0_wp, 10.0_wp]
   call new(mol, [8, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp, 0.0_wp, 0.0_wp], [3, 2]))
   mol%lattice = lat
   mol%periodic = [.true., .true., .false.]
   call get_charges(model, mol, expected, q)
   if (.not. allocated(expected)) then
      call test_failed(error, "Degenerate 2D collinear lattice was not rejected")
      return
   end if
   deallocate(expected)
end subroutine test_degenerate_lattice_error

subroutine test_eeq_2d_all_masks(error)
   type(error_type), allocatable, intent(out) :: error
   type(error_type), allocatable :: expected
   type(structure_type) :: mol, mol_big
   class(mchrg_model_type), allocatable :: model
   real(wp) :: xyz(3, 3), a1(3), a2(3), lat1(3, 3), lat2(3, 3)
   real(wp) :: q1(3), q2(3), q_ref(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3)
   real(wp) :: dqdL1(3, 3, 3), dqdL2(3, 3, 3)
   logical :: masks(3, 3)
   integer :: imask, p1, p2, np

   xyz = reshape([ &
      & 0.0_wp, 0.0_wp, -0.3_wp, &
      & 2.5_wp, 0.2_wp, 0.5_wp, &
      & 0.6_wp, 2.1_wp, 0.4_wp], [3, 3])
   a1 = [5.2_wp, 0.3_wp, 0.1_wp]
   a2 = [0.5_wp, 6.1_wp, -0.2_wp]
   masks = reshape([ &
      & .true., .true., .false., &
      & .true., .false., .true., &
      & .false., .true., .true.], [3, 3])

   do imask = 1, 3
      select case (imask)
      case (1); p1 = 1; p2 = 2; np = 3
      case (2); p1 = 1; p2 = 3; np = 2
      case default; p1 = 2; p2 = 3; np = 1
      end select
      lat1 = 0.0_wp
      lat1(:, p1) = a1
      lat1(:, p2) = a2
      lat1(:, np) = [0.0_wp, 0.0_wp, 20.0_wp]

      lat2 = 0.0_wp
      lat2(:, p1) = a1
      lat2(:, p2) = a2
      lat2(:, np) = [35.0_wp, -70.0_wp, 200.0_wp]

      call new(mol, [30, 8, 1], xyz, charge=0.0_wp, &
         & lattice=lat1, periodic=masks(:, imask))
      call new(mol_big, [30, 8, 1], xyz, charge=0.0_wp, &
         & lattice=lat2, periodic=masks(:, imask))

      call new_eeq2019_model(mol, model, error)
      if (allocated(error)) return

      call get_charges(model, mol, error, q1, dqdr1, dqdL1)
      if (allocated(error)) return
      call get_charges(model, mol_big, error, q2, dqdr2, dqdL2)
      if (allocated(error)) return

      if (abs(sum(q1)) > thr) then
         call test_failed(error, "2D EEQ total charge constraint violated")
         return
      end if
      if (maxval(abs(q1 - q2)) > thr .or. maxval(abs(dqdr1 - dqdr2)) > thr &
         & .or. maxval(abs(dqdL1 - dqdL2)) > thr) then
         call test_failed(error, "2D EEQ depends on nonperiodic lattice vector")
         return
      end if
      if (imask == 1) then
         q_ref = q1
      else if (maxval(abs(q1 - q_ref)) > thr) then
         call test_failed(error, "2D EEQ charges depend on periodic mask index")
         return
      end if

      ! Full derivative validation for all 2D masks (TTF, TFT, FTT)
      call test_dbdr(error, mol, model)
      if (allocated(error)) return
      call test_dbdL(error, mol, model)
      if (allocated(error)) return
      call test_dadr(error, mol, model)
      if (allocated(error)) return
      call test_dadL(error, mol, model)
      if (allocated(error)) return
      call test_numgrad(error, mol, model)
      if (allocated(error)) return
      call test_numsigma(error, mol, model)
      if (allocated(error)) return
      call test_numdqdr(error, mol, model)
      if (allocated(error)) return
      call test_numdqdL(error, mol, model)
      if (allocated(error)) return
   end do

   ! Charged 2D EEQ cell must be explicitly rejected
   mol%charge = 1.0_wp
   call get_charges(model, mol, expected, q1)
   if (.not. allocated(expected)) then
      call test_failed(error, "Charged 2D EEQ cell was not rejected")
      return
   end if
   deallocate(expected)
end subroutine test_eeq_2d_all_masks

subroutine test_eeq_1d_all_masks(error)
   type(error_type), allocatable, intent(out) :: error
   type(error_type), allocatable :: expected
   type(structure_type) :: mol, mol_big
   class(mchrg_model_type), allocatable :: model
   real(wp) :: xyz(3, 3), a_wire(3), lat1(3, 3), lat2(3, 3)
   real(wp) :: q1(3), q2(3), q_ref(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3)
   real(wp) :: dqdL1(3, 3, 3), dqdL2(3, 3, 3)
   logical :: masks(3, 3)
   integer :: imask

   ! Includes both an on-axis atom pair (atoms 1 and 2 have rho = 0 along a_wire)
   ! and off-axis atom pairs (atom 3 has rho > 0)
   a_wire = [5.3_wp, 0.6_wp, -0.5_wp]
   xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
   xyz(:, 2) = 0.42_wp * a_wire
   xyz(:, 3) = [1.1_wp, 1.8_wp, 0.7_wp]
   masks = reshape([ &
      & .true., .false., .false., &
      & .false., .true., .false., &
      & .false., .false., .true.], [3, 3])

   do imask = 1, 3
      lat1 = 0.0_wp
      lat1(:, imask) = a_wire
      lat1(:, mod(imask, 3) + 1) = [0.0_wp, 20.0_wp, 0.0_wp]
      lat1(:, mod(imask + 1, 3) + 1) = [0.0_wp, 0.0_wp, 20.0_wp]

      lat2 = 0.0_wp
      lat2(:, imask) = a_wire
      lat2(:, mod(imask, 3) + 1) = [15.0_wp, 200.0_wp, -30.0_wp]
      lat2(:, mod(imask + 1, 3) + 1) = [-40.0_wp, 10.0_wp, 250.0_wp]

      call new(mol, [30, 8, 1], xyz, charge=0.0_wp, &
         & lattice=lat1, periodic=masks(:, imask))
      call new(mol_big, [30, 8, 1], xyz, charge=0.0_wp, &
         & lattice=lat2, periodic=masks(:, imask))

      call new_eeq2019_model(mol, model, error)
      if (allocated(error)) return

      call get_charges(model, mol, error, q1, dqdr1, dqdL1)
      if (allocated(error)) return
      call get_charges(model, mol_big, error, q2, dqdr2, dqdL2)
      if (allocated(error)) return

      if (abs(sum(q1)) > thr) then
         call test_failed(error, "1D EEQ total charge constraint violated")
         return
      end if
      if (maxval(abs(q1 - q2)) > thr .or. maxval(abs(dqdr1 - dqdr2)) > thr &
         & .or. maxval(abs(dqdL1 - dqdL2)) > thr) then
         call test_failed(error, "1D EEQ depends on nonperiodic lattice vectors")
         return
      end if
      if (imask == 1) then
         q_ref = q1
      else if (maxval(abs(q1 - q_ref)) > thr) then
         call test_failed(error, "1D EEQ charges depend on periodic mask index")
         return
      end if

      ! Full derivative validation for all 1D masks (TFF, FTF, FFT)
      call test_dbdr(error, mol, model)
      if (allocated(error)) return
      call test_dbdL(error, mol, model)
      if (allocated(error)) return
      call test_dadr(error, mol, model)
      if (allocated(error)) return
      call test_dadL(error, mol, model)
      if (allocated(error)) return
      call test_numgrad(error, mol, model)
      if (allocated(error)) return
      call test_numsigma(error, mol, model)
      if (allocated(error)) return
      call test_numdqdr(error, mol, model)
      if (allocated(error)) return
      call test_numdqdL(error, mol, model)
      if (allocated(error)) return
   end do

   ! Charged 1D EEQ cell must be explicitly rejected
   mol%charge = 1.0_wp
   call get_charges(model, mol, expected, q1)
   if (.not. allocated(expected)) then
      call test_failed(error, "Charged 1D EEQ cell was not rejected")
      return
   end if
   deallocate(expected)
end subroutine test_eeq_1d_all_masks


!> Verify uncapped Ewald convergence, resource-limit errors, and alpha invariance
!> across adverse 1D and 2D geometries.
subroutine test_ewald_convergence_1d_2d(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(error_type), allocatable :: expected
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: cache
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), cn(:), qloc(:), trans(:, :)
   real(wp) :: a1(3), a2(3), u_hat(3), n_hat(3), len1, area
   real(wp) :: amat_ref(4, 4), amat_cur(4, 4)
   real(wp) :: dadr_ref(3, 3, 4), dadL_ref(3, 3, 4), atrace_ref(3, 3)
   real(wp) :: dadr_cur(3, 3, 4), dadL_cur(3, 3, 4), atrace_cur(3, 3)
   real(wp) :: hess_ref(9, 9), hess_cur(9, 9)
   real(wp) :: qvec(4), adj(4), alphas(3)
   integer :: ia, ig

   ! 1. Verify removal of hard caps (>40 in 1D, >20 per axis in 2D)
   a1 = [2.0_wp, 0.0_wp, 0.0_wp]
   call get_ewald_dir_trans_1d(a1, 0.06_wp, 0.06_wp, dtrans, error)
   if (allocated(error)) return
   if (size(dtrans, 2) <= 81) then
      call test_failed(error, "1D direct-space Ewald remains capped at 40 shells")
      return
   end if

   a1 = [45.0_wp, 0.0_wp, 0.0_wp]
   call get_ewald_rec_trans_1d(a1, 0.50_wp, rtrans, len1, u_hat, error)
   if (allocated(error)) return
   if (size(rtrans, 2) <= 80) then
      call test_failed(error, "1D reciprocal-space Ewald remains capped at 40 shells")
      return
   end if

   a1 = [3.0_wp, 0.0_wp, 0.0_wp]
   a2 = [0.0_wp, 3.0_wp, 0.0_wp]
   call get_ewald_dir_trans_2d(a1, a2, 0.08_wp, 0.08_wp, dtrans, error)
   if (allocated(error)) return
   if (size(dtrans, 2) <= 1681) then
      call test_failed(error, "2D direct-space Ewald remains capped at 20 shells")
      return
   end if

   a1 = [30.0_wp, 0.0_wp, 0.0_wp]
   a2 = [0.0_wp, 30.0_wp, 0.0_wp]
   call get_ewald_rec_trans_2d(a1, a2, 0.45_wp, rtrans, area, n_hat, error)
   if (allocated(error)) return
   if (size(rtrans, 2) <= 1680) then
      call test_failed(error, "2D reciprocal-space Ewald remains capped at 20 shells")
      return
   end if

   ! 2. Verify explicit resource-limit error reporting (no silent truncation)
   call get_ewald_dir_trans_1d([2.0_wp, 0.0_wp, 0.0_wp], 0.06_wp, 0.06_wp, dtrans, &
      & expected, max_trans=10)
   if (.not. allocated(expected)) then
      call test_failed(error, "1D direct-space resource limit was not enforced")
      return
   end if

   call get_ewald_rec_trans_1d([45.0_wp, 0.0_wp, 0.0_wp], 0.50_wp, rtrans, len1, u_hat, &
      & expected, max_trans=10)
   if (.not. allocated(expected)) then
      call test_failed(error, "1D reciprocal-space resource limit was not enforced")
      return
   end if

   call get_ewald_dir_trans_2d([3.0_wp, 0.0_wp, 0.0_wp], [0.0_wp, 3.0_wp, 0.0_wp], &
      & 0.08_wp, 0.08_wp, dtrans, expected, max_trans=10)
   if (.not. allocated(expected)) then
      call test_failed(error, "2D direct-space resource limit was not enforced")
      return
   end if

   call get_ewald_rec_trans_2d([30.0_wp, 0.0_wp, 0.0_wp], [0.0_wp, 30.0_wp, 0.0_wp], &
      & 0.45_wp, rtrans, area, n_hat, expected, max_trans=10)
   if (.not. allocated(expected)) then
      call test_failed(error, "2D reciprocal-space resource limit was not enforced")
      return
   end if

   call get_ewald_dir_trans_3d(reshape([3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 3.0_wp], [3, 3]), 0.08_wp, 0.08_wp, dtrans, expected, max_trans=10)
   if (.not. allocated(expected)) then
      call test_failed(error, "3D direct-space resource limit was not enforced")
      return
   end if

   call get_ewald_rec_trans_3d(reshape([30.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 30.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 30.0_wp], [3, 3]), 0.45_wp, rtrans, expected, max_trans=10)
   if (.not. allocated(expected)) then
      call test_failed(error, "3D reciprocal-space resource limit was not enforced")
      return
   end if

   ! Public get_charges API must propagate resource-limit errors cleanly when
   ! an extreme aspect-ratio 2D cell requires > 200,000 translations
   call new(mol, [8, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp, 0.5_wp, 0.0_wp], [3, 2]), &
      & charge=0.0_wp, lattice=reshape([ &
      & 0.01_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 5000.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 20.0_wp], [3, 3]), periodic=[.true., .true., .false.])
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call get_charges(model, mol, expected, qvec(1:2))
   if (.not. allocated(expected)) then
      call test_failed(error, "Public get_charges did not propagate Ewald resource-limit error")
      return
   end if
   deallocate(expected, model)

   ! 3. Alpha-invariance across adverse 1D and 2D geometries:
   ! ig = 1: Short-period 1D wire with broad Gaussians
   ! ig = 2: Strongly skewed 2D lattice a1=(10,0,0), a2=(9,2,0)
   ! ig = 3: High-aspect-ratio 2D lattice (3.6 x 18.0 bohr)
   ! ig = 4: Large 2D periodic cell (16.0 x 18.0 bohr)
   alphas = [0.16_wp, 0.26_wp, 0.40_wp]
   qvec = [0.35_wp, -0.55_wp, 0.20_wp, 0.05_wp]
   adj = [-0.40_wp, 0.25_wp, 0.15_wp, 0.0_wp]
   allocate(cn(3), qloc(3))

   do ig = 1, 4
      select case (ig)
      case (1)
         call new(mol, [8, 1, 6], reshape([ &
            & 0.0_wp, 0.0_wp, 0.0_wp, &
            & 1.1_wp, 0.8_wp, -0.3_wp, &
            & 0.5_wp, -0.9_wp, 0.6_wp], [3, 3]), charge=0.0_wp, &
            & lattice=reshape([ &
            & 2.6_wp, 0.0_wp, 0.0_wp, &
            & 0.0_wp, 20.0_wp, 0.0_wp, &
            & 0.0_wp, 0.0_wp, 20.0_wp], [3, 3]), &
            & periodic=[.true., .false., .false.])
      case (2)
         call new(mol, [8, 1, 1], reshape([ &
            & 0.1_wp, 0.1_wp, -0.2_wp, &
            & 1.2_wp, 0.5_wp, 0.3_wp, &
            & -0.4_wp, 0.9_wp, 0.1_wp], [3, 3]), charge=0.0_wp, &
            & lattice=reshape([ &
            & 10.0_wp, 0.0_wp, 0.0_wp, &
            & 9.0_wp, 2.0_wp, 0.0_wp, &
            & 0.0_wp, 0.0_wp, 20.0_wp], [3, 3]), &
            & periodic=[.true., .true., .false.])
      case (3)
         call new(mol, [8, 1, 1], reshape([ &
            & 0.0_wp, 0.0_wp, -0.2_wp, &
            & 1.4_wp, 1.5_wp, 0.4_wp, &
            & 0.7_wp, 4.2_wp, -0.1_wp], [3, 3]), charge=0.0_wp, &
            & lattice=reshape([ &
            & 3.6_wp, 0.0_wp, 0.0_wp, &
            & 0.0_wp, 18.0_wp, 0.0_wp, &
            & 0.0_wp, 0.0_wp, 20.0_wp], [3, 3]), &
            & periodic=[.true., .true., .false.])
      case (4)
         call new(mol, [8, 1, 1], reshape([ &
            & 0.0_wp, 0.0_wp, -0.2_wp, &
            & 1.8_wp, 0.4_wp, 0.3_wp, &
            & -0.5_wp, 1.7_wp, 0.2_wp], [3, 3]), charge=0.0_wp, &
            & lattice=reshape([ &
            & 16.0_wp, 0.0_wp, 0.0_wp, &
            & 0.0_wp, 18.0_wp, 0.0_wp, &
            & 0.0_wp, 0.0_wp, 20.0_wp], [3, 3]), &
            & periodic=[.true., .true., .false.])
      case default
         continue
      end select

      call new_eeq2019_model(mol, model, error)
      if (allocated(error)) return
      if (ig == 1) then
         ! Broaden Gaussian charge widths on the short-period 1D wire
         model%rad = 2.2_wp * model%rad
      end if

      call get_lattice_points(mol%periodic, mol%lattice, 25.0_wp, trans, error)
      if (allocated(error)) return
      call model%ncoord%get_coordination_number(mol, trans, cn)
      call model%local_charge(mol, trans, qloc)
      call model%update(mol, cache, cn, qloc)
      if (allocated(cache%error)) then
         call move_alloc(cache%error, error)
         return
      end if
      allocate(cache%qvec(4), source=qvec)
      select type (raw => cache%raw)
      type is (eeq_cache)
         allocate(raw%dcndr(3, 3, 3), source=0.0_wp)
      end select

      do ia = 1, size(alphas)
         select type (raw => cache%raw)
         type is (eeq_cache)
            raw%alpha = alphas(ia)
         end select
         call model%get_coulomb_matrix(mol, cache, amat_cur)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            return
         end if
         call model%get_coulomb_derivs(mol, cache, qvec, dadr_cur, dadL_cur, atrace_cur)
         if (allocated(cache%error)) then
            call move_alloc(cache%error, error)
            return
         end if
         hess_cur = 0.0_wp
         call model%add_model_hessian(mol, error, cache, adj, hess_cur)
         if (allocated(error)) return

         if (ia == 1) then
            amat_ref = amat_cur
            dadr_ref = dadr_cur
            dadL_ref = dadL_cur
            atrace_ref = atrace_cur
            hess_ref = hess_cur
         else
            if (maxval(abs(amat_cur - amat_ref)) > 1.0e-11_wp) then
               call test_failed(error, "Coulomb matrix depends on Ewald alpha")
               return
            end if
            if (maxval(abs(dadr_cur - dadr_ref)) > 1.0e-11_wp &
               & .or. maxval(abs(dadL_cur - dadL_ref)) > 1.0e-11_wp &
               & .or. maxval(abs(atrace_cur - atrace_ref)) > 1.0e-11_wp) then
               call test_failed(error, "Coulomb derivatives depend on Ewald alpha")
               return
            end if
            if (maxval(abs(hess_cur - hess_ref)) > 1.0e-10_wp) then
               call test_failed(error, "Coulomb Hessian depends on Ewald alpha")
               return
            end if
         end if
      end do
      deallocate(cache%qvec, cache%raw)
   end do

end subroutine test_ewald_convergence_1d_2d

!> Verify Wigner-Seitz invariance under multi-cell shifts, unimodular basis transforms,
!> and smoothness across Wigner-Seitz cell boundaries in 1D and 2D.
subroutine test_pbc_wsc_invariance_1d_2d(error)
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol_ref, mol_test, mol_m, mol_0, mol_p
   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: cache
   real(wp) :: a1(3), a2(3), lat_ref(3, 3), lat_test(3, 3), xyz_ref(3, 3), xyz_test(3, 3)
   real(wp) :: q_ref(3), dqdr_ref(3, 3, 3), dqdL_ref(3, 3, 3), hess_ref(9, 9)
   real(wp) :: q_test(3), dqdr_test(3, 3, 3), dqdL_test(3, 3, 3), hess_test(9, 9)
   real(wp) :: q_m(2), dqdr_m(3, 2, 2), dqdL_m(3, 3, 2), hess_m(6, 6)
   real(wp) :: q_0(2), dqdr_0(3, 2, 2), dqdL_0(3, 3, 2), hess_0(6, 6)
   real(wp) :: q_p(2), dqdr_p(3, 2, 2), dqdL_p(3, 3, 2), hess_p(6, 6)
   real(wp) :: cn2(2), qloc2(2), dcndr2(3, 2, 2), dcndL2(3, 3, 2), dqlocdr2(3, 2, 2), dqlocdL2(3, 3, 2)
   real(wp), allocatable :: trans(:, :)
   real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
   real(wp), parameter :: angle_65 = 65.0_wp * pi / 180.0_wp
   real(wp), parameter :: weights3(3) = [0.6_wp, -0.4_wp, 0.25_wp]
   real(wp), parameter :: weights2(2) = [0.7_wp, -0.3_wp]

   ! 1. 1D wire rotated in 3D space with zero nonperiodic columns and multi-cell shifts (|n| > 1)
   a1 = [3.6_wp, 4.8_wp, 0.0_wp]
   lat_ref = 0.0_wp
   lat_ref(:, 3) = a1
   xyz_ref = reshape([ &
      & 0.0_wp,  0.0_wp,  0.1_wp, &
      & 1.1_wp,  1.6_wp, -0.2_wp, &
      & 2.1_wp,  2.7_wp,  0.3_wp], [3, 3])
   xyz_test = xyz_ref
   xyz_test(:, 2) = xyz_test(:, 2) + 4.0_wp * a1
   xyz_test(:, 3) = xyz_test(:, 3) - 3.0_wp * a1

   call new(mol_ref, [8, 1, 1], xyz_ref, charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.false., .false., .true.])
   call new(mol_test, [8, 1, 1], xyz_test, charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.false., .false., .true.])

   call new_eeq2019_model(mol_ref, model, error)
   if (allocated(error)) return
   call get_charges(model, mol_ref, error, q_ref, dqdr_ref, dqdL_ref, cache=cache)
   if (allocated(error)) return
   hess_ref = 0.0_wp
   call add_charge_hessian(model, mol_ref, error, weights3, hess_ref, cache)
   if (allocated(error)) return

   call get_charges(model, mol_test, error, q_test, dqdr_test, dqdL_test, cache=cache)
   if (allocated(error)) return
   hess_test = 0.0_wp
   call add_charge_hessian(model, mol_test, error, weights3, hess_test, cache)
   if (allocated(error)) return

   if (maxval(abs(q_test - q_ref)) > 1.0e-12_wp &
      & .or. maxval(abs(dqdr_test - dqdr_ref)) > 1.0e-11_wp &
      & .or. maxval(abs(dqdL_test - dqdL_ref)) > 1.0e-11_wp &
      & .or. maxval(abs(hess_test - hess_ref)) > 1.0e-10_wp) then
      call test_failed(error, "1D EEQ not invariant under multi-cell atom translation")
      return
   end if

   ! 2. 2D slab in rotated 3D subspace with zero nonperiodic column,
   !    skewed unimodular basis transformation (a2' = a2 + 5*a1), and multi-cell shifts
   a1 = [5.5_wp, 0.0_wp, 0.0_wp]
   a2 = [0.5_wp, 6.0_wp, 2.0_wp]
   lat_ref = 0.0_wp
   lat_ref(:, 1) = a1
   lat_ref(:, 3) = a2
   lat_test = 0.0_wp
   lat_test(:, 1) = a1
   lat_test(:, 3) = a2 + 5.0_wp * a1

   xyz_ref = reshape([ &
      & 0.1_wp, 0.2_wp, -0.1_wp, &
      & 1.7_wp, 0.6_wp,  0.3_wp, &
      & 0.3_wp, 1.8_wp,  0.5_wp], [3, 3])
   xyz_test = xyz_ref
   xyz_test(:, 2) = xyz_test(:, 2) + 3.0_wp * a1 - 2.0_wp * a2
   xyz_test(:, 3) = xyz_test(:, 3) - 4.0_wp * a1 + 3.0_wp * a2

   call new(mol_ref, [8, 1, 1], xyz_ref, charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.true., .false., .true.])
   call new(mol_test, [8, 1, 1], xyz_test, charge=0.0_wp, &
      & lattice=lat_test, periodic=[.true., .false., .true.])

   call get_charges(model, mol_ref, error, q_ref, dqdr_ref, dqdL_ref, cache=cache)
   if (allocated(error)) return
   hess_ref = 0.0_wp
   call add_charge_hessian(model, mol_ref, error, weights3, hess_ref, cache)
   if (allocated(error)) return

   call get_charges(model, mol_test, error, q_test, dqdr_test, dqdL_test, cache=cache)
   if (allocated(error)) return
   hess_test = 0.0_wp
   call add_charge_hessian(model, mol_test, error, weights3, hess_test, cache)
   if (allocated(error)) return

   if (maxval(abs(q_test - q_ref)) > 1.0e-11_wp &
      & .or. maxval(abs(dqdr_test - dqdr_ref)) > 1.0e-10_wp &
      & .or. maxval(abs(dqdL_test - dqdL_ref)) > 1.0e-10_wp &
      & .or. maxval(abs(hess_test - hess_ref)) > 1.0e-9_wp) then
      call test_failed(error, &
         & "2D EEQ not invariant under unimodular basis transform and multi-cell shift")
      return
   end if

   ! 3. Smoothness across a 2D Wigner-Seitz cell boundary (x = L/2 - 1e-6, L/2, L/2 + 1e-6)
   lat_ref = 0.0_wp
   lat_ref(:, 1) = [6.0_wp, 0.0_wp, 0.0_wp]
   lat_ref(:, 2) = [0.0_wp, 6.0_wp, 0.0_wp]
   call new(mol_m, [8, 1], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 3.0_wp - 1.0e-6_wp, 1.2_wp, 0.4_wp], [3, 2]), charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.true., .true., .false.])
   call new(mol_0, [8, 1], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 3.0_wp, 1.2_wp, 0.4_wp], [3, 2]), charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.true., .true., .false.])
   call new(mol_p, [8, 1], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 3.0_wp + 1.0e-6_wp, 1.2_wp, 0.4_wp], [3, 2]), charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.true., .true., .false.])

   deallocate(model)
   call new_eeq2019_model(mol_0, model, error)
   if (allocated(error)) return

   call get_charges(model, mol_m, error, q_m, dqdr_m, dqdL_m, cache=cache)
   if (allocated(error)) return
   hess_m = 0.0_wp
   call add_charge_hessian(model, mol_m, error, weights2, hess_m, cache)
   if (allocated(error)) return

   call get_charges(model, mol_0, error, q_0, dqdr_0, dqdL_0, cache=cache)
   if (allocated(error)) return
   hess_0 = 0.0_wp
   call add_charge_hessian(model, mol_0, error, weights2, hess_0, cache)
   if (allocated(error)) return

   call get_charges(model, mol_p, error, q_p, dqdr_p, dqdL_p, cache=cache)
   if (allocated(error)) return
   hess_p = 0.0_wp
   call add_charge_hessian(model, mol_p, error, weights2, hess_p, cache)
   if (allocated(error)) return

   if (maxval(abs(0.5_wp * (q_m + q_p) - q_0)) > 1.0e-10_wp &
      & .or. maxval(abs(0.5_wp * (dqdr_m + dqdr_p) - dqdr_0)) > 1.0e-9_wp &
      & .or. maxval(abs(0.5_wp * (dqdL_m + dqdL_p) - dqdL_0)) > 1.0e-9_wp &
      & .or. maxval(abs(0.5_wp * (hess_m + hess_p) - hess_0)) > 1.0e-8_wp) then
      call test_failed(error, "2D EEQ exhibits discontinuity across Wigner-Seitz boundary")
      return
   end if

   ! 4. 65-degree 2D counterexample (unwrapped 1.39*(a1+a2) vs wrapped 0.39*(a1+a2))
   !    tested through both public get_charges and direct model%solve for EEQ and EEQBC
   lat_ref = 0.0_wp
   lat_ref(:, 1) = [5.0_wp, 0.0_wp, 0.0_wp]
   lat_ref(:, 2) = [5.0_wp * cos(angle_65), 5.0_wp * sin(angle_65), 0.0_wp]
   call new(mol_0, [6, 8], reshape([ &
      & [0.0_wp, 0.0_wp, 0.0_wp], &
      & 0.39_wp * (lat_ref(:, 1) + lat_ref(:, 2))], [3, 2]), charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.true., .true., .false.])
   call new(mol_p, [6, 8], reshape([ &
      & [0.0_wp, 0.0_wp, 0.0_wp], &
      & 1.39_wp * (lat_ref(:, 1) + lat_ref(:, 2))], [3, 2]), charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.true., .true., .false.])

   deallocate(model)
   call new_eeqbc2025_model(mol_0, model, error)
   if (allocated(error)) return
   call get_charges(model, mol_0, error, q_0, dqdr_0, dqdL_0, cache=cache)
   if (allocated(error)) return
   hess_0 = 0.0_wp
   call add_charge_hessian(model, mol_0, error, weights2, hess_0, cache)
   if (allocated(error)) return

   ! Direct model%solve on unwrapped 65-degree counterexample (mol_p) using CN/qloc from mol_0
   call get_lattice_points(mol_0%periodic, mol_0%lattice, model%ncoord%cutoff, trans, error)
   if (allocated(error)) return
   call model%ncoord%get_coordination_number(mol_0, trans, cn2, dcndr2, dcndL2)
   call model%local_charge(mol_0, trans, qloc2, dqlocdr2, dqlocdL2)
   call model%solve(mol_p, error, cn2, qloc2, dcndr2, dcndL2, dqlocdr2, dqlocdL2, &
      & qvec=q_p, dqdr=dqdr_p, dqdL=dqdL_p, response=cache)
   if (allocated(error)) return
   hess_p = 0.0_wp
   call model%add_charge_hessian(mol_p, error, weights2, hess_p, cache)
   if (allocated(error)) return

   if (maxval(abs(q_p - q_0)) > 1.0e-12_wp &
      & .or. maxval(abs(dqdr_p - dqdr_0)) > 1.0e-11_wp &
      & .or. maxval(abs(dqdL_p - dqdL_0)) > 1.0e-11_wp &
      & .or. maxval(abs(hess_p - hess_0)) > 1.0e-10_wp) then
      call test_failed(error, "65-degree 2D counterexample failed in direct model%solve")
      return
   end if

   ! 5. 3D periodic multi-cell shift invariance via get_charges and add_charge_hessian
   lat_ref = reshape([ &
      & 5.4_wp, 0.0_wp, 0.0_wp, &
      & 0.5_wp, 6.1_wp, 0.0_wp, &
      & 0.3_wp, 0.6_wp, 6.8_wp], [3, 3])
   xyz_ref = reshape([ &
      & 0.2_wp, 0.3_wp, 0.4_wp, &
      & 2.1_wp, 1.2_wp, 0.9_wp, &
      & 0.8_wp, 2.5_wp, 1.7_wp], [3, 3])
   xyz_test = xyz_ref
   xyz_test(:, 2) = xyz_test(:, 2) + 3.0_wp * lat_ref(:, 1) - 2.0_wp * lat_ref(:, 2)
   xyz_test(:, 3) = xyz_test(:, 3) - 2.0_wp * lat_ref(:, 1) + 4.0_wp * lat_ref(:, 3)

   call new(mol_ref, [30, 8, 1], xyz_ref, charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.true., .true., .true.])
   call new(mol_test, [30, 8, 1], xyz_test, charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.true., .true., .true.])

   deallocate(model)
   call new_eeqbc2025_model(mol_ref, model, error)
   if (allocated(error)) return
   call get_charges(model, mol_ref, error, q_ref, dqdr_ref, dqdL_ref, cache=cache)
   if (allocated(error)) return
   hess_ref = 0.0_wp
   call add_charge_hessian(model, mol_ref, error, weights3, hess_ref, cache)
   if (allocated(error)) return

   call get_charges(model, mol_test, error, q_test, dqdr_test, dqdL_test, cache=cache)
   if (allocated(error)) return
   hess_test = 0.0_wp
   call add_charge_hessian(model, mol_test, error, weights3, hess_test, cache)
   if (allocated(error)) return

   if (maxval(abs(q_test - q_ref)) > 1.0e-12_wp &
      & .or. maxval(abs(dqdr_test - dqdr_ref)) > 1.0e-11_wp &
      & .or. maxval(abs(dqdL_test - dqdL_ref)) > 1.0e-11_wp &
      & .or. maxval(abs(hess_test - hess_ref)) > 1.0e-10_wp) then
      call test_failed(error, "3D get_charges/add_charge_hessian not invariant under multi-cell shift")
      return
   end if

   ! 6. 3D skewed-cell counterexample (a1=(5,0,0), a2=(25,1,0), a3=(0,0,5), r2=0.5*a2)
   !    vs LLL-reduced basis (a1'=(5,0,0), a2'=(0,1,0), a3'=(0,0,5), r2'=r2-2*a1=(2.5,0.5,0))
   !    across get_charges, add_charge_hessian, model%solve, and model%add_charge_hessian
   lat_ref = reshape([ &
      & 5.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 5.0_wp], [3, 3])
   lat_test = reshape([ &
      &  5.0_wp, 0.0_wp, 0.0_wp, &
      & 25.0_wp, 1.0_wp, 0.0_wp, &
      &  0.0_wp, 0.0_wp, 5.0_wp], [3, 3])
   call new(mol_0, [6, 8], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 2.5_wp, 0.5_wp, 0.0_wp], [3, 2]), charge=0.0_wp, &
      & lattice=lat_ref, periodic=[.true., .true., .true.])
   call new(mol_p, [6, 8], reshape([ &
      &  0.0_wp, 0.0_wp, 0.0_wp, &
      & 12.5_wp, 0.5_wp, 0.0_wp], [3, 2]), charge=0.0_wp, &
      & lattice=lat_test, periodic=[.true., .true., .true.])

   deallocate(model)
   call new_eeq2019_model(mol_0, model, error)
   if (allocated(error)) return
   call get_charges(model, mol_0, error, q_0, dqdr_0, dqdL_0, cache=cache)
   if (allocated(error)) return
   hess_0 = 0.0_wp
   call add_charge_hessian(model, mol_0, error, weights2, hess_0, cache)
   if (allocated(error)) return

   call get_charges(model, mol_p, error, q_p, dqdr_p, dqdL_p, cache=cache)
   if (allocated(error)) return
   hess_p = 0.0_wp
   call model%add_charge_hessian(mol_p, error, weights2, hess_p, cache)
   if (allocated(error)) return

   if (maxval(abs(q_p - q_0)) > 1.0e-11_wp &
      & .or. maxval(abs(dqdr_p - dqdr_0)) > 1.0e-10_wp &
      & .or. maxval(abs(dqdL_p - dqdL_0)) > 1.0e-10_wp &
      & .or. maxval(abs(hess_p - hess_0)) > 1.0e-9_wp) then
      call test_failed(error, "3D skewed-cell counterexample failed for EEQ")
      return
   end if

   deallocate(model)
   call new_eeqbc2025_model(mol_0, model, error)
   if (allocated(error)) return
   call get_charges(model, mol_0, error, q_0, dqdr_0, dqdL_0, cache=cache)
   if (allocated(error)) return
   hess_0 = 0.0_wp
   call add_charge_hessian(model, mol_0, error, weights2, hess_0, cache)
   if (allocated(error)) return

   call get_charges(model, mol_p, error, q_p, dqdr_p, dqdL_p, cache=cache)
   if (allocated(error)) return
   hess_p = 0.0_wp
   call model%add_charge_hessian(mol_p, error, weights2, hess_p, cache)
   if (allocated(error)) return

   if (maxval(abs(q_p - q_0)) > 1.0e-11_wp &
      & .or. maxval(abs(dqdr_p - dqdr_0)) > 1.0e-10_wp &
      & .or. maxval(abs(dqdL_p - dqdL_0)) > 1.0e-10_wp &
      & .or. maxval(abs(hess_p - hess_0)) > 1.0e-9_wp) then
      call test_failed(error, "3D skewed-cell counterexample failed for EEQBC")
      return
   end if

end subroutine test_pbc_wsc_invariance_1d_2d

!> Verify 1D and 2D EEQBC direct-lattice summation convergence across screening
!> parameters, lattice lengths, skewed cells, and neutral/charged systems, as well
!> as multi-cell shift and unimodular basis invariance.
subroutine test_eeqbc_convergence_1d_2d(error)
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol, mol_test
   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: cache_cur, cache_ref
   real(wp) :: cmat_cur(4, 4), cmat_ref(4, 4), amat_cur(4, 4), amat_ref(4, 4)
   real(wp) :: q_cur(3), q_ref(3), dqdr_cur(3, 3, 3), dqdr_ref(3, 3, 3)
   real(wp) :: dqdL_cur(3, 3, 3), dqdL_ref(3, 3, 3), hess_cur(9, 9), hess_ref(9, 9)
   real(wp) :: a1(3), a2(3), lat_ref(3, 3), lat_test(3, 3), xyz_ref(3, 3), xyz_test(3, 3)
   real(wp), parameter :: weights(3) = [0.6_wp, -0.4_wp, 0.25_wp]
   integer :: icase

   ! Test convergence against an independently over-converged direct sum (cutoff = 120 Bohr)
   ! across 4 cases covering short/long cells, small/large kbc, skewed 2D, neutral/charged:
   do icase = 1, 4
      select case (icase)
      case (1)
         ! Short-period neutral 1D cell with small kbc = 0.40 (long-range capacitance tail)
         call new(mol, [8, 1, 6], reshape([ &
            & 0.0_wp,  0.0_wp,  0.0_wp, &
            & 1.2_wp,  0.7_wp, -0.3_wp, &
            & 0.5_wp, -0.8_wp,  0.6_wp], [3, 3]), charge=0.0_wp, &
            & lattice=reshape([ &
            & 3.2_wp, 0.0_wp, 0.0_wp, &
            & 0.0_wp, 20.0_wp, 0.0_wp, &
            & 0.0_wp, 0.0_wp, 20.0_wp], [3, 3]), &
            & periodic=[.true., .false., .false.])
      case (2)
         ! Longer charged (+1) 1D cell with larger kbc = 1.10
         call new(mol, [30, 8, 1], reshape([ &
            & 0.0_wp, 0.0_wp,  0.0_wp, &
            & 2.4_wp, 0.3_wp, -0.4_wp, &
            & 1.1_wp, 1.8_wp,  0.7_wp], [3, 3]), charge=1.0_wp, &
            & lattice=reshape([ &
            & 12.0_wp, 0.0_wp, 0.0_wp, &
            & 0.0_wp, 20.0_wp, 0.0_wp, &
            & 0.0_wp, 0.0_wp, 20.0_wp], [3, 3]), &
            & periodic=[.true., .false., .false.])
      case (3)
         ! Strongly skewed neutral 2D cell with default kbc = 0.65
         call new(mol, [8, 1, 1], reshape([ &
            & 0.1_wp, 0.1_wp, -0.2_wp, &
            & 1.4_wp, 0.5_wp,  0.3_wp, &
            &-0.3_wp, 1.1_wp,  0.1_wp], [3, 3]), charge=0.0_wp, &
            & lattice=reshape([ &
            & 8.0_wp, 0.0_wp, 0.0_wp, &
            & 7.2_wp, 2.2_wp, 0.0_wp, &
            & 0.0_wp, 0.0_wp, 20.0_wp], [3, 3]), &
            & periodic=[.true., .true., .false.])
      case (4)
         ! Charged (-1) 2D cell with small kbc = 0.45
         call new(mol, [30, 8, 1], reshape([ &
            & 0.0_wp, 0.0_wp, -0.3_wp, &
            & 2.3_wp, 0.2_wp,  0.5_wp, &
            & 0.6_wp, 1.9_wp,  0.4_wp], [3, 3]), charge=-1.0_wp, &
            & lattice=reshape([ &
            & 4.5_wp, 0.0_wp, 0.0_wp, &
            & 0.4_wp, 5.2_wp, 0.0_wp, &
            & 0.0_wp, 0.0_wp, 20.0_wp], [3, 3]), &
            & periodic=[.true., .true., .false.])
      case default
         continue
      end select

      if (allocated(model)) deallocate(model)
      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      select type (m => model)
      type is (eeqbc_model)
         select case (icase)
         case (1); m%kbc = 0.40_wp
         case (2); m%kbc = 1.10_wp
         case (3); m%kbc = 0.65_wp
         case (4); m%kbc = 0.45_wp
         case default; m%kbc = 0.65_wp
         end select
      end select

      ! Default cutoff evaluation
      call get_charges(model, mol, error, q_cur, dqdr_cur, dqdL_cur, cache=cache_cur)
      if (allocated(error)) return
      call model%get_coulomb_matrix(mol, cache_cur, amat_cur)
      select type (raw => cache_cur%raw)
      type is (eeqbc_cache)
         cmat_cur = raw%cmat
      end select
      hess_cur = 0.0_wp
      call add_charge_hessian(model, mol, error, weights, hess_cur, cache_cur)
      if (allocated(error)) return

      ! Over-converged reference evaluation (cutoff = 120 Bohr)
      model%ncoord%cutoff = 120.0_wp
      model%ncoord_en%cutoff = 120.0_wp
      call get_charges(model, mol, error, q_ref, dqdr_ref, dqdL_ref, cache=cache_ref)
      if (allocated(error)) return
      call model%get_coulomb_matrix(mol, cache_ref, amat_ref)
      select type (raw => cache_ref%raw)
      type is (eeqbc_cache)
         cmat_ref = raw%cmat
      end select
      hess_ref = 0.0_wp
      call add_charge_hessian(model, mol, error, weights, hess_ref, cache=cache_ref)
      if (allocated(error)) return

      if (maxval(abs(cmat_cur - cmat_ref)) > 1.0e-13_wp) then
         call test_failed(error, "1D/2D EEQBC capacitance matrix cmat is not converged")
         return
      end if
      if (maxval(abs(amat_cur - amat_ref)) > 1.0e-13_wp) then
         call test_failed(error, "1D/2D EEQBC Coulomb matrix amat is not converged")
         return
      end if
      if (maxval(abs(q_cur - q_ref)) > 1.0e-13_wp &
         & .or. maxval(abs(dqdr_cur - dqdr_ref)) > 1.0e-12_wp &
         & .or. maxval(abs(dqdL_cur - dqdL_ref)) > 1.0e-12_wp &
         & .or. maxval(abs(hess_cur - hess_ref)) > 1.0e-11_wp) then
         call test_failed(error, "1D/2D EEQBC charges or derivatives are not cutoff-converged")
         return
      end if
   end do

   ! Verify 2D EEQBC invariance under unimodular basis transformation and multi-cell shifts
   a1 = [5.4_wp, 0.2_wp, 0.0_wp]
   a2 = [0.6_wp, 5.8_wp, 1.2_wp]
   lat_ref = 0.0_wp
   lat_ref(:, 1) = a1
   lat_ref(:, 2) = a2
   lat_test = 0.0_wp
   lat_test(:, 1) = a1
   lat_test(:, 2) = a2 + 4.0_wp * a1

   xyz_ref = reshape([ &
      & 0.0_wp, 0.0_wp, -0.2_wp, &
      & 2.1_wp, 0.5_wp,  0.3_wp, &
      & 0.5_wp, 1.9_wp,  0.4_wp], [3, 3])
   xyz_test = xyz_ref
   xyz_test(:, 2) = xyz_test(:, 2) + 3.0_wp * a1 - 2.0_wp * a2
   xyz_test(:, 3) = xyz_test(:, 3) - 3.0_wp * a1 + 2.0_wp * a2

   call new(mol, [30, 8, 1], xyz_ref, charge=-1.0_wp, &
      & lattice=lat_ref, periodic=[.true., .true., .false.])
   call new(mol_test, [30, 8, 1], xyz_test, charge=-1.0_wp, &
      & lattice=lat_test, periodic=[.true., .true., .false.])

   deallocate(model)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return

   call get_charges(model, mol, error, q_ref, dqdr_ref, dqdL_ref, cache=cache_ref)
   if (allocated(error)) return
   hess_ref = 0.0_wp
   call add_charge_hessian(model, mol, error, weights, hess_ref, cache=cache_ref)
   if (allocated(error)) return

   call get_charges(model, mol_test, error, q_cur, dqdr_cur, dqdL_cur, cache=cache_cur)
   if (allocated(error)) return
   hess_cur = 0.0_wp
   call add_charge_hessian(model, mol_test, error, weights, hess_cur, cache=cache_cur)
   if (allocated(error)) return

   if (maxval(abs(q_cur - q_ref)) > 1.0e-12_wp &
      & .or. maxval(abs(dqdr_cur - dqdr_ref)) > 1.0e-11_wp &
      & .or. maxval(abs(dqdL_cur - dqdL_ref)) > 1.0e-11_wp &
      & .or. maxval(abs(hess_cur - hess_ref)) > 1.0e-10_wp) then
      call test_failed(error, &
         & "2D EEQBC not invariant under unimodular basis transform and multi-cell shift")
      return
   end if

end subroutine test_eeqbc_convergence_1d_2d

!> Regression and characterization test for DFT-D4 issue #358: asymptotic EEQBC
!> decoupling across 0D, 1D, and 2D for:
!> 1. Neutral water-helium complex across separations up to 1000 A (~1889.73 Bohr),
!> 2. Two identical neutral molecules (H2O...H2O) with zero net charge transfer,
!> 3. Two dissimilar neutral molecules (H2O...NH3) with analytical Lagrange multiplier check,
!> 4. Charged fragment combinations (H2O...NH4+, Q = +1), and
!> 5. Invariance under fragment reordering and rigid translation.
subroutine test_eeqbc_water_helium_358(error)
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: aa2au = 1.8897261246257702_wp
   real(wp), parameter :: z_seps(5) = [20.0_wp, 35.0_wp, 60.0_wp, 200.0_wp, 1000.0_wp * aa2au]
   real(wp), parameter :: q_he_0d_ref = 0.07309211_wp
   type(structure_type) :: mol_h2o, mol_he, mol_comb, mol_nh3, mol_nh4
   type(structure_type) :: mol_dimer, mol_het, mol_het_rev, mol_het_shift, mol_chg
   class(mchrg_model_type), allocatable :: model_eeq, model_eeqbc
   type(cache_container) :: cache_eeq, cache_eeqbc, cache_h2o, cache_he, cache_nh3, cache_nh4
   real(wp) :: xyz_h2o(3, 3), xyz_he(3, 1), xyz_comb(3, 4), lat(3, 3)
   real(wp) :: xyz_nh3(3, 4), xyz_nh4(3, 5), xyz_dimer(3, 6), xyz_het(3, 7), xyz_het_rev(3, 7)
   real(wp) :: xyz_het_shift(3, 7), xyz_chg(3, 8), shift_vec(3)
   real(wp) :: q_h2o_bc(3), dqdr_h2o_bc(3, 3, 3), dqdL_h2o_bc(3, 3, 3), q_he_bc(1)
   real(wp) :: q_nh3_bc(4), q_nh4_bc(5)
   real(wp) :: amat_h2o(4, 4), cmat_h2o(4, 4), amat_he(2, 2), cmat_he(2, 2)
   real(wp) :: amat_nh3(5, 5), cmat_nh3(5, 5), amat_nh4(6, 6), cmat_nh4(6, 6)
   real(wp) :: res_h2o(3), res_he, res_nh3(4), res_nh4(5)
   real(wp) :: lam_h2o, lam_he, lam_nh3, lam_nh4, lam_comb
   real(wp) :: q_eeq(4, 5), dqdr_eeq(3, 4, 4), dqdL_eeq(3, 3, 4)
   real(wp) :: q_bc(4, 5), dqdr_bc(3, 4, 4), dqdL_bc(3, 3, 4)
   real(wp) :: amat_eeq(5, 5), amat_bc(5, 5), cmat_bc(5, 5)
   real(wp) :: q_dimer(6), dqdr_dimer(3, 6, 6), dqdL_dimer(3, 3, 6)
   real(wp) :: q_het(7), dqdr_het(3, 7, 7), dqdL_het(3, 3, 7), amat_het(8, 8), hess_het(21, 21)
   real(wp) :: q_het_rev(7), q_het_shift(7), weights7(7)
   real(wp) :: q_chg(8), q_chg_ref(8), dqdr_chg(3, 8, 8), dqdL_chg(3, 3, 8), amat_chg(9, 9)
   logical :: pmask(3)
   integer :: idim, ir

   xyz_h2o = reshape([ &
      &  0.0000000000_wp,  0.0000000000_wp,  0.2214080000_wp, &
      &  0.0000000000_wp,  1.4309280000_wp, -0.8856320000_wp, &
      &  0.0000000000_wp, -1.4309280000_wp, -0.8856320000_wp], [3, 3])
   xyz_he = reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1])
   xyz_nh3 = reshape([ &
      &  0.0000000000_wp,  0.0000000000_wp,  0.2190000000_wp, &
      &  0.0000000000_wp,  1.7720000000_wp, -0.5110000000_wp, &
      &  1.5345970000_wp, -0.8860000000_wp, -0.5110000000_wp, &
      & -1.5345970000_wp, -0.8860000000_wp, -0.5110000000_wp], [3, 4])
   xyz_nh4 = reshape([ &
      &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, &
      &  1.1200000000_wp,  1.1200000000_wp,  1.1200000000_wp, &
      &  1.1200000000_wp, -1.1200000000_wp, -1.1200000000_wp, &
      & -1.1200000000_wp,  1.1200000000_wp, -1.1200000000_wp, &
      & -1.1200000000_wp, -1.1200000000_wp,  1.1200000000_wp], [3, 5])
   xyz_comb(:, 1:3) = xyz_h2o
   weights7 = [0.4_wp, -0.2_wp, 0.1_wp, -0.3_wp, 0.25_wp, -0.15_wp, 0.35_wp]
   shift_vec = [1.7_wp, -2.3_wp, 4.1_wp]

   lat = reshape([ &
      & 12.0_wp,  0.0_wp,    0.0_wp, &
      &  0.0_wp, 12.0_wp,    0.0_wp, &
      &  0.0_wp,  0.0_wp, 4000.0_wp], [3, 3])

   do idim = 0, 2
      select case (idim)
      case (0); pmask = [.false., .false., .false.]
      case (1); pmask = [.true.,  .false., .false.]
      case (2); pmask = [.true.,  .true.,  .false.]
      case default; pmask = [.false., .false., .false.]
      end select

      ! 1. Isolated H2O, He, NH3, and NH4+ reference fragments
      if (idim == 0) then
         call new(mol_h2o, [8, 1, 1], xyz_h2o, charge=0.0_wp)
         call new(mol_he, [2], xyz_he, charge=0.0_wp)
         call new(mol_nh3, [7, 1, 1, 1], xyz_nh3, charge=0.0_wp)
         call new(mol_nh4, [7, 1, 1, 1, 1], xyz_nh4, charge=1.0_wp)
      else
         call new(mol_h2o, [8, 1, 1], xyz_h2o, charge=0.0_wp, lattice=lat, periodic=pmask)
         call new(mol_he, [2], xyz_he, charge=0.0_wp, lattice=lat, periodic=pmask)
         call new(mol_nh3, [7, 1, 1, 1], xyz_nh3, charge=0.0_wp, lattice=lat, periodic=pmask)
         call new(mol_nh4, [7, 1, 1, 1, 1], xyz_nh4, charge=1.0_wp, lattice=lat, periodic=pmask)
      end if

      if (allocated(model_eeqbc)) deallocate(model_eeqbc)
      call new_eeqbc2025_model(mol_h2o, model_eeqbc, error)
      if (allocated(error)) return
      call get_charges(model_eeqbc, mol_h2o, error, q_h2o_bc, dqdr_h2o_bc, dqdL_h2o_bc, &
         & cache=cache_h2o)
      if (allocated(error)) return
      call model_eeqbc%get_coulomb_matrix(mol_h2o, cache_h2o, amat_h2o)
      cmat_h2o = 0.0_wp
      select type (raw => cache_h2o%raw)
      type is (eeqbc_cache)
         cmat_h2o = raw%cmat
      end select
      lam_h2o = cache_h2o%qvec(4)

      deallocate(model_eeqbc)
      call new_eeqbc2025_model(mol_he, model_eeqbc, error)
      if (allocated(error)) return
      call get_charges(model_eeqbc, mol_he, error, q_he_bc, cache=cache_he)
      if (allocated(error)) return
      call model_eeqbc%get_coulomb_matrix(mol_he, cache_he, amat_he)
      cmat_he = 0.0_wp
      select type (raw => cache_he%raw)
      type is (eeqbc_cache)
         cmat_he = raw%cmat
      end select
      lam_he = cache_he%qvec(2)
      if (abs(q_he_bc(1)) > 1.0e-15_wp) then
         call test_failed(error, "Isolated neutral He atom has non-zero charge")
         return
      end if

      deallocate(model_eeqbc)
      call new_eeqbc2025_model(mol_nh3, model_eeqbc, error)
      if (allocated(error)) return
      call get_charges(model_eeqbc, mol_nh3, error, q_nh3_bc, cache=cache_nh3)
      if (allocated(error)) return
      call model_eeqbc%get_coulomb_matrix(mol_nh3, cache_nh3, amat_nh3)
      cmat_nh3 = 0.0_wp
      select type (raw => cache_nh3%raw)
      type is (eeqbc_cache)
         cmat_nh3 = raw%cmat
      end select
      lam_nh3 = cache_nh3%qvec(5)

      deallocate(model_eeqbc)
      call new_eeqbc2025_model(mol_nh4, model_eeqbc, error)
      if (allocated(error)) return
      call get_charges(model_eeqbc, mol_nh4, error, q_nh4_bc, cache=cache_nh4)
      if (allocated(error)) return
      call model_eeqbc%get_coulomb_matrix(mol_nh4, cache_nh4, amat_nh4)
      cmat_nh4 = 0.0_wp
      select type (raw => cache_nh4%raw)
      type is (eeqbc_cache)
         cmat_nh4 = raw%cmat
      end select
      lam_nh4 = cache_nh4%qvec(6)

      ! 2. Combined H2O...He system across progressive separations up to 1000 A
      do ir = 1, size(z_seps)
         xyz_comb(:, 4) = [0.0_wp, 0.0_wp, z_seps(ir)]
         if (idim == 0) then
            call new(mol_comb, [8, 1, 1, 2], xyz_comb, charge=0.0_wp)
         else
            call new(mol_comb, [8, 1, 1, 2], xyz_comb, charge=0.0_wp, lattice=lat, periodic=pmask)
         end if

         if (ir == 1) then
            if (allocated(model_eeq)) deallocate(model_eeq)
            if (allocated(model_eeqbc)) deallocate(model_eeqbc)
            call new_eeq2019_model(mol_comb, model_eeq, error)
            if (allocated(error)) return
            call new_eeqbc2025_model(mol_comb, model_eeqbc, error)
            if (allocated(error)) return
         end if

         call get_charges(model_eeq, mol_comb, error, q_eeq(:, ir), dqdr_eeq, dqdL_eeq, cache=cache_eeq)
         if (allocated(error)) return
         call model_eeq%get_coulomb_matrix(mol_comb, cache_eeq, amat_eeq)

         call get_charges(model_eeqbc, mol_comb, error, q_bc(:, ir), dqdr_bc, dqdL_bc, cache=cache_eeqbc)
         if (allocated(error)) return
         call model_eeqbc%get_coulomb_matrix(mol_comb, cache_eeqbc, amat_bc)
         cmat_bc = 0.0_wp
         select type (raw => cache_eeqbc%raw)
         type is (eeqbc_cache)
            cmat_bc = raw%cmat
         end select
         lam_comb = cache_eeqbc%qvec(5)

         ! Conventional EEQ exhibits artificial long-range charge transfer even at 1000 A in 0D
         if (idim == 0 .and. abs(q_eeq(4, ir)) < 1.0e-2_wp) then
            call test_failed(error, "Expected residual asymptotic charge transfer in 0D conventional EEQ")
            return
         end if

         ! For R >= 35 Bohr (including 1000 A), EEQBC screens out inter-fragment coupling,
         ! decouples the H2O and He diagonal blocks to match the isolated fragments exactly,
         ! and suppresses all distance dependence and spatial derivatives on He
         if (ir >= 2) then
            if (maxval(abs(cmat_bc(1:3, 4))) > 1.0e-14_wp &
               & .or. maxval(abs(amat_bc(1:3, 4))) > 1.0e-14_wp) then
               call test_failed(error, "EEQBC failed to screen long-range H2O-He coupling")
               return
            end if
            if (maxval(abs(cmat_bc(1:3, 1:3) - cmat_h2o(1:3, 1:3))) > 1.0e-14_wp &
               & .or. maxval(abs(amat_bc(1:3, 1:3) - amat_h2o(1:3, 1:3))) > 1.0e-14_wp &
               & .or. abs(cmat_bc(4, 4) - cmat_he(1, 1)) > 1.0e-14_wp &
               & .or. abs(amat_bc(4, 4) - amat_he(1, 1)) > 1.0e-14_wp) then
               call test_failed(error, "EEQBC failed to recover isolated H2O and He matrix blocks")
               return
            end if
            res_h2o = matmul(amat_bc(1:3, 1:3), q_bc(1:3, ir) - q_h2o_bc) + (lam_comb - lam_h2o)
            res_he = amat_bc(4, 4) * (q_bc(4, ir) - q_he_bc(1)) + (lam_comb - lam_he)
            if (maxval(abs(res_h2o)) > 1.0e-13_wp &
               & .or. abs(res_he) > 1.0e-13_wp) then
               call test_failed(error, "EEQBC fragment blocks inconsistent with isolated fragments")
               return
            end if
            if (maxval(abs(q_bc(:, ir) - q_bc(:, 2))) > 1.0e-13_wp &
               & .or. maxval(abs(dqdr_bc(:, 4, :))) > 1.0e-14_wp) then
               call test_failed(error, "EEQBC failed to suppress asymptotic charge response on He")
               return
            end if
            if (idim == 0 .and. (abs(q_bc(4, ir) + lam_comb) > 1.0e-14_wp &
               & .or. abs(q_bc(4, ir) - q_he_0d_ref) > 1.0e-7_wp)) then
               call test_failed(error, "0D H2O...He asymptotic He charge deviates from dftd4 EEQBC reference")
               return
            end if
         end if
      end do
      if (abs(q_eeq(4, 5) - q_eeq(4, 1)) < 1.0e-4_wp) then
         call test_failed(error, "Expected distance-dependent 1/R coupling in unscreened EEQ")
         return
      end if

      ! 3. Two identical neutral molecules (H2O...H2O) at 1000 A:
      !    Equal isolated Lagrange multipliers => strictly zero asymptotic charge transfer
      xyz_dimer(:, 1:3) = xyz_h2o
      xyz_dimer(1:2, 4:6) = xyz_h2o(1:2, :)
      xyz_dimer(3, 4:6) = xyz_h2o(3, :) + z_seps(5)
      if (idim == 0) then
         call new(mol_dimer, [8, 1, 1, 8, 1, 1], xyz_dimer, charge=0.0_wp)
      else
         call new(mol_dimer, [8, 1, 1, 8, 1, 1], xyz_dimer, charge=0.0_wp, lattice=lat, periodic=pmask)
      end if
      deallocate(model_eeqbc)
      call new_eeqbc2025_model(mol_dimer, model_eeqbc, error)
      if (allocated(error)) return
      call get_charges(model_eeqbc, mol_dimer, error, q_dimer, dqdr_dimer, dqdL_dimer, &
         & cache=cache_eeqbc)
      if (allocated(error)) return
      if (abs(sum(q_dimer(1:3))) > 1.0e-13_wp &
         & .or. abs(sum(q_dimer(4:6))) > 1.0e-13_wp &
         & .or. maxval(abs(q_dimer(1:3) - q_h2o_bc)) > 1.0e-13_wp &
         & .or. maxval(abs(q_dimer(4:6) - q_h2o_bc)) > 1.0e-13_wp &
         & .or. abs(cache_eeqbc%qvec(7) - lam_h2o) > 1.0e-13_wp &
         & .or. maxval(abs(sum(dqdr_dimer(:, 1:3, :), dim=2))) > 1.0e-13_wp &
         & .or. maxval(abs(sum(dqdr_dimer(:, 4:6, :), dim=2))) > 1.0e-13_wp) then
         call test_failed(error, "Identical neutral H2O dimer at 1000 A failed zero-transfer check")
         return
      end if

      ! 4. Two dissimilar neutral molecules (H2O...NH3) at 1000 A:
      !    Verify exact block-diagonal decoupling, vanishing rigid-translation inter-fragment
      !    derivatives, analytical Lagrange-multiplier identity, and reordering/shift invariance
      xyz_het(:, 1:3) = xyz_h2o
      xyz_het(1:2, 4:7) = xyz_nh3(1:2, :)
      xyz_het(3, 4:7) = xyz_nh3(3, :) + z_seps(5)
      xyz_het_rev(:, 1:4) = xyz_het(:, 4:7)
      xyz_het_rev(:, 5:7) = xyz_het(:, 1:3)
      do ir = 1, 7
         xyz_het_shift(:, ir) = xyz_het(:, ir) + shift_vec
      end do
      if (idim == 0) then
         call new(mol_het, [8, 1, 1, 7, 1, 1, 1], xyz_het, charge=0.0_wp)
         call new(mol_het_rev, [7, 1, 1, 1, 8, 1, 1], xyz_het_rev, charge=0.0_wp)
         call new(mol_het_shift, [8, 1, 1, 7, 1, 1, 1], xyz_het_shift, charge=0.0_wp)
      else
         call new(mol_het, [8, 1, 1, 7, 1, 1, 1], xyz_het, charge=0.0_wp, lattice=lat, periodic=pmask)
         call new(mol_het_rev, [7, 1, 1, 1, 8, 1, 1], xyz_het_rev, charge=0.0_wp, &
            & lattice=lat, periodic=pmask)
         call new(mol_het_shift, [8, 1, 1, 7, 1, 1, 1], xyz_het_shift, charge=0.0_wp, &
            & lattice=lat, periodic=pmask)
      end if
      deallocate(model_eeqbc)
      call new_eeqbc2025_model(mol_het, model_eeqbc, error)
      if (allocated(error)) return
      call get_charges(model_eeqbc, mol_het, error, q_het, dqdr_het, dqdL_het, cache=cache_eeqbc)
      if (allocated(error)) return
      call model_eeqbc%get_coulomb_matrix(mol_het, cache_eeqbc, amat_het)
      hess_het = 0.0_wp
      call add_charge_hessian(model_eeqbc, mol_het, error, weights7, hess_het, cache_eeqbc)
      if (allocated(error)) return

      lam_comb = cache_eeqbc%qvec(8)
      res_h2o = matmul(amat_h2o(1:3, 1:3), q_het(1:3) - q_h2o_bc) + (lam_comb - lam_h2o)
      res_nh3 = matmul(amat_nh3(1:4, 1:4), q_het(4:7) - q_nh3_bc) + (lam_comb - lam_nh3)
      if (maxval(abs(amat_het(1:3, 4:7))) > 1.0e-14_wp &
         & .or. maxval(abs(res_h2o)) > 1.0e-13_wp &
         & .or. maxval(abs(res_nh3)) > 1.0e-13_wp &
         & .or. maxval(abs(sum(dqdr_het(:, 1:3, :), dim=2))) > 1.0e-13_wp &
         & .or. maxval(abs(sum(dqdr_het(:, 4:7, :), dim=2))) > 1.0e-13_wp) then
         call test_failed(error, "Dissimilar neutral H2O...NH3 complex failed decoupling check")
         return
      end if

      call get_charges(model_eeqbc, mol_het_shift, error, q_het_shift)
      if (allocated(error)) return
      deallocate(model_eeqbc)
      call new_eeqbc2025_model(mol_het_rev, model_eeqbc, error)
      if (allocated(error)) return
      call get_charges(model_eeqbc, mol_het_rev, error, q_het_rev)
      if (allocated(error)) return
      if (maxval(abs(q_het_shift - q_het)) > 1.0e-13_wp &
         & .or. maxval(abs(q_het_rev(1:4) - q_het(4:7))) > 1.0e-13_wp &
         & .or. maxval(abs(q_het_rev(5:7) - q_het(1:3))) > 1.0e-13_wp) then
         call test_failed(error, "EEQBC asymptotic charges not invariant to shift or fragment reordering")
         return
      end if

      ! 5. Charged fragment combination (H2O...NH4+, Q = +1) at R = 60 Bohr vs 1000 A:
      !    Verify exact off-diagonal block decoupling (A_AB = 0, C_AB = 0), distance
      !    independence beyond cutoff, charge conservation, and zero rigid-shift response
      xyz_chg(:, 1:3) = xyz_h2o
      xyz_chg(1:2, 4:8) = xyz_nh4(1:2, :)
      xyz_chg(3, 4:8) = xyz_nh4(3, :) + z_seps(3)
      if (idim == 0) then
         call new(mol_chg, [8, 1, 1, 7, 1, 1, 1, 1], xyz_chg, charge=1.0_wp)
      else
         call new(mol_chg, [8, 1, 1, 7, 1, 1, 1, 1], xyz_chg, charge=1.0_wp, &
            & lattice=lat, periodic=pmask)
      end if
      deallocate(model_eeqbc)
      call new_eeqbc2025_model(mol_chg, model_eeqbc, error)
      if (allocated(error)) return
      call get_charges(model_eeqbc, mol_chg, error, q_chg_ref, cache=cache_eeqbc)
      if (allocated(error)) return
      lam_nh4 = cache_eeqbc%qvec(9)

      xyz_chg(3, 4:8) = xyz_nh4(3, :) + z_seps(5)
      if (idim == 0) then
         call new(mol_chg, [8, 1, 1, 7, 1, 1, 1, 1], xyz_chg, charge=1.0_wp)
      else
         call new(mol_chg, [8, 1, 1, 7, 1, 1, 1, 1], xyz_chg, charge=1.0_wp, &
            & lattice=lat, periodic=pmask)
      end if
      call get_charges(model_eeqbc, mol_chg, error, q_chg, dqdr_chg, dqdL_chg, &
         & cache=cache_eeqbc)
      if (allocated(error)) return
      call model_eeqbc%get_coulomb_matrix(mol_chg, cache_eeqbc, amat_chg)
      lam_comb = cache_eeqbc%qvec(9)
      if (abs(sum(q_chg) - 1.0_wp) > 1.0e-13_wp &
         & .or. maxval(abs(amat_chg(1:3, 4:8))) > 1.0e-14_wp &
         & .or. maxval(abs(q_chg - q_chg_ref)) > 1.0e-13_wp &
         & .or. abs(lam_comb - lam_nh4) > 1.0e-13_wp &
         & .or. maxval(abs(sum(dqdr_chg(:, 1:3, :), dim=2))) > 1.0e-13_wp &
         & .or. maxval(abs(sum(dqdr_chg(:, 4:8, :), dim=2))) > 1.0e-13_wp) then
         call test_failed(error, "Charged H2O...NH4+ complex at 1000 A failed decoupling check")
         return
      end if
   end do

end subroutine test_eeqbc_water_helium_358

!> Regression test for DFT-D4 issue #285: 3D periodic slab with large vacuum gap
!> and in-plane lattice scan across exact hexagonal symmetry.
subroutine test_eeq_3d_slab_285(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol2d, mol3d
   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: cache
   real(wp), parameter :: aa2au = 1.8897261246257702_wp
   real(wp), parameter :: pi = 3.14159265358979323846264338327950288_wp
   real(wp), parameter :: q2d_hex_ref(3) = [ &
      & -1.06091408_wp, 0.54516560_wp, 0.51574848_wp]
   real(wp), parameter :: q2d_ahex_ref(3) = [ &
      & -1.06069500_wp, 0.54509370_wp, 0.51560131_wp]
   real(wp), parameter :: q2d_orth_ref(3) = [ &
      & -1.06073102_wp, 0.54653780_wp, 0.51419322_wp]
   real(wp) :: xyz(3, 3), lat(3, 3), cn(3), qloc(3)
   real(wp), allocatable :: trans(:, :)
   real(wp) :: q2d(3), q_z25(3), q_z50(3), q_z100(3)
   real(wp) :: amat0(4, 4), amat_s(4, 4), alpha0, scales(2), vol, g0_shift
   real(wp) :: q_scan(3, 13), e_scan(3, 13), etot(13), x_ang, d2q, d2e
   integer :: ix, is

   xyz = reshape([ &
      & -0.3910542064024872_wp,  0.3976645623910747_wp, 2.977929077294636_wp, &
      & -0.3542054272959850_wp,  1.0567259314317760_wp, 2.610987477729398_wp, &
      & -0.9445722089116435_wp, -0.0304014548863456_wp, 3.370602014356701_wp], [3, 3]) * aa2au

   lat = 0.0_wp
   lat(:, 1) = [2.794646093662685_wp, 0.0_wp, 0.0_wp] * aa2au
   lat(:, 2) = [-1.397323046831344_wp, 2.420234511698831_wp, 0.0_wp] * aa2au
   lat(:, 3) = [0.0_wp, 0.0_wp, 25.0_wp] * aa2au

   call new(mol2d, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .false.])
   call new_eeq2019_model(mol2d, model, error)
   if (allocated(error)) return
   call get_charges(model, mol2d, error, q2d)
   if (allocated(error)) return
   if (maxval(abs(q2d - q2d_hex_ref)) > 1.0e-7_wp) then
      call test_failed(error, "2D hexagonal #285 slab charges deviate from reference")
      return
   end if

   ! 1. Verify large-vacuum 3D supercell convergence toward true 2D EEQ
   call new(mol3d, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .true.])
   call get_charges(model, mol3d, error, q_z25, cache=cache)
   if (allocated(error)) return

   ! Check 3D alpha-invariance on the slab cell (accounting for the omitted G=0
   ! constant -pi / (alpha^2 * V) in raw 3D A_ij, which shifts into lambda)
   call model%get_coulomb_matrix(mol3d, cache, amat0)
   select type (raw => cache%raw)
   type is (eeq_cache)
      vol = lat(1, 1) * lat(2, 2) * lat(3, 3)
      alpha0 = raw%alpha
      scales = [0.75_wp, 1.3_wp]
      do is = 1, size(scales)
         raw%alpha = alpha0 * scales(is)
         call model%get_coulomb_matrix(mol3d, cache, amat_s)
         g0_shift = (pi / vol) * (1.0_wp / raw%alpha**2 - 1.0_wp / alpha0**2)
         if (maxval(abs((amat_s(1:3, 1:3) - amat0(1:3, 1:3)) - g0_shift)) > 1.0e-12_wp) then
            call test_failed(error, "3D Ewald matrix is not invariant to alpha on slab cell")
            return
         end if
      end do
      raw%alpha = alpha0
   end select

   ! Check 3D analytical derivatives on the hexagonal slab cell
   call test_numgrad(error, mol3d, model)
   if (allocated(error)) return
   call test_numsigma(error, mol3d, model)
   if (allocated(error)) return
   call test_numdqdr(error, mol3d, model)
   if (allocated(error)) return
   call test_numdqdL(error, mol3d, model)
   if (allocated(error)) return

   lat(:, 3) = [0.0_wp, 0.0_wp, 50.0_wp] * aa2au
   call new(mol3d, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .true.])
   call get_charges(model, mol3d, error, q_z50)
   if (allocated(error)) return

   lat(:, 3) = [0.0_wp, 0.0_wp, 100.0_wp] * aa2au
   call new(mol3d, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .true.])
   call get_charges(model, mol3d, error, q_z100)
   if (allocated(error)) return

   if (maxval(abs(q_z25 - q2d)) > 6.0e-5_wp &
      & .or. maxval(abs(q_z50 - q2d)) > 3.0e-5_wp &
      & .or. maxval(abs(q_z100 - q2d)) > 1.5e-5_wp) then
      call test_failed(error, "3D slab charges failed to converge toward 2D EEQ limit")
      return
   end if

   ! Yeh-Berkowitz 1/Lz dipole-sheet scaling: Lz * (q_3D(Lz) - q_2D) is constant
   if (maxval(abs(50.0_wp * (q_z50 - q2d) - 100.0_wp * (q_z100 - q2d))) > 1.0e-5_wp) then
      call test_failed(error, "3D-to-2D slab difference does not follow 1/Lz asymptotic scaling")
      return
   end if

   ! 2. Verify smoothness across hexagonal symmetry (X = 2.7940..2.7952 A at Z = 100 A)
   do ix = 1, 13
      x_ang = 2.0_wp + real(7939 + ix, wp) / 10000.0_wp
      lat(:, 1) = [x_ang, 0.0_wp, 0.0_wp] * aa2au
      call new(mol3d, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .true.])
      call get_lattice_points(mol3d%periodic, mol3d%lattice, model%ncoord%cutoff, trans, error)
      if (allocated(error)) return
      call model%ncoord%get_coordination_number(mol3d, trans, cn)
      call model%local_charge(mol3d, trans, qloc)
      e_scan(:, ix) = 0.0_wp
      call model%solve(mol3d, error, cn, qloc, energy=e_scan(:, ix), qvec=q_scan(:, ix))
      if (allocated(error)) return
      etot(ix) = sum(e_scan(:, ix))
   end do

   d2q = maxval(abs(q_scan(:, 3:13) - 2.0_wp * q_scan(:, 2:12) + q_scan(:, 1:11)))
   d2e = maxval(abs(etot(3:13) - 2.0_wp * etot(2:12) + etot(1:11)))
   if (d2q > 1.0e-8_wp .or. d2e > 1.0e-8_wp) then
      call test_failed(error, "Discontinuity detected across hexagonal symmetry in 3D slab X-scan")
      return
   end if

   ! 3. Verify almost_hexagonal and orthorhombic #285 slab geometries at Z = 100 A
   lat(:, 1) = [2.8_wp, 0.0_wp, 0.0_wp] * aa2au
   lat(:, 2) = [-1.397323046831344_wp, 2.420234511698831_wp, 0.0_wp] * aa2au
   lat(:, 3) = [0.0_wp, 0.0_wp, 100.0_wp] * aa2au
   call new(mol2d, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .false.])
   call new(mol3d, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .true.])
   call get_charges(model, mol2d, error, q2d)
   if (allocated(error)) return
   call get_charges(model, mol3d, error, q_z100)
   if (allocated(error)) return
   if (maxval(abs(q2d - q2d_ahex_ref)) > 1.0e-7_wp &
      & .or. maxval(abs(q_z100 - q2d)) > 1.5e-5_wp) then
      call test_failed(error, "almost_hexagonal #285 slab charges deviate from reference")
      return
   end if

   lat(:, 1) = [2.8_wp, 0.0_wp, 0.0_wp] * aa2au
   lat(:, 2) = [0.0_wp, 2.8_wp, 0.0_wp] * aa2au
   lat(:, 3) = [0.0_wp, 0.0_wp, 100.0_wp] * aa2au
   call new(mol2d, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .false.])
   call new(mol3d, [8, 1, 1], xyz, lattice=lat, periodic=[.true., .true., .true.])
   call get_charges(model, mol2d, error, q2d)
   if (allocated(error)) return
   call get_charges(model, mol3d, error, q_z100)
   if (allocated(error)) return
   if (maxval(abs(q2d - q2d_orth_ref)) > 1.0e-7_wp &
      & .or. maxval(abs(q_z100 - q2d)) > 1.5e-5_wp) then
      call test_failed(error, "orthorhombic #285 slab charges deviate from reference")
      return
   end if

end subroutine test_eeq_3d_slab_285

!> Verify 1D and 2D special-function accuracy (fast_f012_1d, fast_k012_1d, erfc_exp_2d)
!> against 50-digit mpmath reference values across small, moderate, and large arguments.
subroutine test_special_functions_1d_2d(error)
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: u_pts(6) = [ &
      & 0.0_wp, 1.0e-8_wp, 0.5_wp, 4.0_wp, 4.000001_wp, 25.0_wp]
   real(wp), parameter :: f0_ref(6) = [ &
      & 0.0_wp, 9.9999999750000000556e-9_wp, 0.44384207911774836294_wp, &
      & 1.9672893784312723859_wp, 1.9672896238523342756_wp, 3.7960914897702684998_wp]
   real(wp), parameter :: f1_ref(6) = [ &
      & 1.0_wp, 0.99999999500000001667_wp, 0.78693868057473315279_wp, &
      & 0.24542109027781645493_wp, 0.24542103350146551229_wp, 0.039999999999444482245_wp]
   real(wp), parameter :: f2_ref(6) = [ &
      & -0.5_wp, -0.49999999666666667917_wp, -0.36081604172419945838_wp, &
      & -0.056776362847270568658_wp, -0.056776339038005506227_wp, -0.0015999999994222615352_wp]

   real(wp), parameter :: x_pts(5) = [1.0e-8_wp, 1.0e-4_wp, 0.5_wp, 2.0_wp, 10.0_wp]
   real(wp), parameter :: y_pts(5) = [0.0_wp, 1.0e-2_wp, 0.5_wp, 0.1_wp, 25.0_wp]
   real(wp), parameter :: k0_ref(5) = [ &
      & 17.843465089050832587_wp, 8.6232592773087082697_wp, 0.42102443824070833334_wp, &
      & 0.045293668272415987242_wp, 8.1625283072617411727e-15_wp]
   real(wp), parameter :: k1_ref(5) = [ &
      & 0.99999981156534915949_wp, 0.99405430461170725993_wp, 0.23402778902579225314_wp, &
      & 0.034642560210863181845_wp, 5.2280069047533075135e-15_wp]
   real(wp), parameter :: k2_ref(5) = [ &
      & 0.49999999000000096717_wp, 0.49657968235608294261_wp, 0.15332113394940819643_wp, &
      & 0.027734685027132461108_wp, 3.4489111320542408121e-15_wp]

   real(wp) :: f0, f1, f2, k0, k1, k2, v_2d
   integer :: i

   do i = 1, size(u_pts)
      call fast_f012_1d(u_pts(i), f0, f1, f2)
      if (abs(f0 - f0_ref(i)) > 1.0e-14_wp * max(1.0_wp, abs(f0_ref(i))) &
         & .or. abs(f1 - f1_ref(i)) > 1.0e-14_wp * max(1.0_wp, abs(f1_ref(i))) &
         & .or. abs(f2 - f2_ref(i)) > 1.0e-14_wp * max(1.0_wp, abs(f2_ref(i)))) then
         call test_failed(error, "fast_f012_1d deviates from 50-digit reference")
         return
      end if
   end do

   do i = 1, size(x_pts)
      call fast_k012_1d(x_pts(i), y_pts(i), k0, k1, k2)
      if (abs(k0 - k0_ref(i)) > 2.0e-14_wp * max(1.0e-15_wp, abs(k0_ref(i))) &
         & .or. abs(k1 - k1_ref(i)) > 2.0e-14_wp * max(1.0e-15_wp, abs(k1_ref(i))) &
         & .or. abs(k2 - k2_ref(i)) > 2.0e-14_wp * max(1.0e-15_wp, abs(k2_ref(i)))) then
         call test_failed(error, "fast_k012_1d deviates from 50-digit reference")
         return
      end if
   end do

   v_2d = erfc_exp_2d(2.0_wp, -1.0_wp, 0.35_wp)
   if (abs(v_2d - 0.00005300645451857405009_wp) > 1.0e-14_wp * 0.00005300645451857405009_wp) then
      call test_failed(error, "erfc_exp_2d (u > 0) deviates from 50-digit reference")
      return
   end if

   v_2d = erfc_exp_2d(0.5_wp, -5.0_wp, 0.35_wp)
   if (abs(v_2d - 0.15243208859081749745_wp) > 1.0e-14_wp * 0.15243208859081749745_wp) then
      call test_failed(error, "erfc_exp_2d (u < 0) deviates from 50-digit reference")
      return
   end if

   v_2d = erfc_exp_2d(5.0_wp, 20.0_wp, 0.35_wp)
   if (abs(v_2d - 1.4504329021329886148e-45_wp) > 5.0e-14_wp * 1.4504329021329886148e-45_wp) then
      call test_failed(error, "erfc_exp_2d (large positive Gs) deviates from 50-digit reference")
      return
   end if

end subroutine test_special_functions_1d_2d

!> Verify structural, model, and array-shape contract validation across public and
!> type-bound APIs without relying on logical short-circuit evaluation.
subroutine test_structure_and_contract_validation(error)
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol, mol_bad
   class(mchrg_model_type), allocatable :: model
   type(eeq_model) :: uninit_eeq
   type(eeqbc_model) :: uninit_eeqbc
   type(cache_container) :: empty_cache
   type(error_type), allocatable :: expected
   real(wp) :: qvec(2), qshort(1), cn(2), qloc(2), hess(6, 6), weights(2)
   real(wp) :: nan_val, inf_val

   nan_val = ieee_value(1.0_wp, ieee_quiet_nan)
   inf_val = ieee_value(1.0_wp, ieee_positive_inf)
   cn = 1.0_wp
   qloc = 0.0_wp
   weights = 1.0_wp

   call new(mol, [8, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp], [3, 2]), &
      & charge=0.0_wp, lattice=reshape([ &
      & 10.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 10.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 10.0_wp], [3, 3]), periodic=[.true., .false., .false.])
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return

   ! 1. nat <= 0
   mol_bad = mol
   mol_bad%nat = 0
   call get_charges(model, mol_bad, expected, qvec)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject mol%nat <= 0")
      return
   end if
   deallocate(expected)

   ! 2. Unallocated / undersized xyz
   mol_bad = mol
   deallocate(mol_bad%xyz)
   call get_charges(model, mol_bad, expected, qvec)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject unallocated mol%xyz")
      return
   end if
   deallocate(expected)

   ! 3. Out-of-bounds species id
   mol_bad = mol
   mol_bad%id(1) = 99
   call get_charges(model, mol_bad, expected, qvec)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject out-of-bounds mol%id")
      return
   end if
   deallocate(expected)

   ! 4. Nonfinite coordinates, charge, or lattice
   mol_bad = mol
   mol_bad%xyz(1, 1) = nan_val
   call get_charges(model, mol_bad, expected, qvec)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject NaN in mol%xyz")
      return
   end if
   deallocate(expected)

   mol_bad = mol
   mol_bad%charge = inf_val
   call get_charges(model, mol_bad, expected, qvec)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject Inf in mol%charge")
      return
   end if
   deallocate(expected)

   mol_bad = mol
   mol_bad%lattice(1, 1) = nan_val
   call get_charges(model, mol_bad, expected, qvec)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject NaN in mol%lattice")
      return
   end if
   deallocate(expected)

   ! 5. Periodic structure with undersized periodic mask or unallocated lattice
   mol_bad = mol
   deallocate(mol_bad%periodic)
   allocate(mol_bad%periodic(2), source=[.true., .false.])
   call get_charges(model, mol_bad, expected, qvec)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject periodic mask of length < 3")
      return
   end if
   deallocate(expected)

   mol_bad = mol
   deallocate(mol_bad%lattice)
   call get_charges(model, mol_bad, expected, qvec)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject periodic structure with unallocated lattice")
      return
   end if
   deallocate(expected)

   ! 6. Uninitialized models and undersized output arrays
   call uninit_eeq%solve(mol, expected, cn, qloc, qvec=qvec)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject uninitialized eeq_model")
      return
   end if
   deallocate(expected)

   call uninit_eeqbc%solve(mol, expected, cn, qloc, qvec=qvec)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject uninitialized eeqbc_model")
      return
   end if
   deallocate(expected)

   call get_charges(model, mol, expected, qshort)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject undersized qvec in get_charges")
      return
   end if
   deallocate(expected)

   call model%solve(mol, expected, cn, qloc, qvec=qshort)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject undersized qvec in model%solve")
      return
   end if
   deallocate(expected)

   ! 7. Calling model%add_charge_hessian with an empty cache must fail cleanly
   hess = 0.0_wp
   call model%add_charge_hessian(mol, expected, weights, hess, empty_cache)
   if (.not. allocated(expected)) then
      call test_failed(error, "Failed to reject empty cache in model%add_charge_hessian")
      return
   end if
   deallocate(expected)

end subroutine test_structure_and_contract_validation

!> Verify charged 1D and 2D contracts: EEQ rejects non-neutral 1D/2D systems across
!> all entry points, whereas EEQBC supports charged 1D/2D systems with exact charge
!> conservation and analytical derivatives.
subroutine test_charged_lowd_contracts(error)
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model_eeq, model_eeqbc
   type(cache_container) :: cache
   type(error_type), allocatable :: expected
   real(wp) :: xyz(3, 3), lat(3, 3), qvec(3), dqdr(3, 3, 3), dqdL(3, 3, 3)
   real(wp) :: cn(3), qloc(3), hess(9, 9), weights(3), charges(2)
   logical :: masks(3, 6)
   integer :: im, iq

   xyz = reshape([ &
      & 0.0_wp,  0.0_wp,  0.0_wp, &
      & 1.4_wp,  0.6_wp, -0.3_wp, &
      & 0.5_wp, -1.1_wp,  0.8_wp], [3, 3])
   lat = reshape([ &
      & 8.5_wp, 0.0_wp, 0.0_wp, &
      & 0.4_wp, 9.0_wp, 0.0_wp, &
      & 0.2_wp, 0.5_wp, 9.5_wp], [3, 3])
   masks(:, 1) = [.true.,  .false., .false.]
   masks(:, 2) = [.false., .true.,  .false.]
   masks(:, 3) = [.false., .false., .true.]
   masks(:, 4) = [.true.,  .true.,  .false.]
   masks(:, 5) = [.true.,  .false., .true.]
   masks(:, 6) = [.false., .true.,  .true.]
   charges = [1.0_wp, -0.5_wp]
   weights = [0.5_wp, -0.3_wp, 0.2_wp]
   cn = 1.0_wp
   qloc = 0.0_wp

   call new(mol, [8, 1, 6], xyz, charge=0.0_wp, lattice=lat, periodic=masks(:, 1))
   call new_eeq2019_model(mol, model_eeq, error)
   if (allocated(error)) return
   call new_eeqbc2025_model(mol, model_eeqbc, error)
   if (allocated(error)) return

   if (model_eeq%supports_charged_lowd()) then
      call test_failed(error, "eeq_model%supports_charged_lowd() must be .false.")
      return
   end if
   if (.not. model_eeqbc%supports_charged_lowd()) then
      call test_failed(error, "eeqbc_model%supports_charged_lowd() must be .true.")
      return
   end if

   do im = 1, 6
      do iq = 1, size(charges)
         call new(mol, [8, 1, 6], xyz, charge=charges(iq), lattice=lat, periodic=masks(:, im))

         ! EEQ must reject charged 1D and 2D across get_charges, solve, and add_charge_hessian
         call get_charges(model_eeq, mol, expected, qvec)
         if (.not. allocated(expected)) then
            call test_failed(error, "EEQ get_charges did not reject charged 1D/2D structure")
            return
         end if
         deallocate(expected)

         call model_eeq%solve(mol, expected, cn, qloc, qvec=qvec)
         if (.not. allocated(expected)) then
            call test_failed(error, "EEQ solve did not reject charged 1D/2D structure")
            return
         end if
         deallocate(expected)

         hess = 0.0_wp
         call add_charge_hessian(model_eeq, mol, expected, weights, hess, cache)
         if (.not. allocated(expected)) then
            call test_failed(error, "EEQ add_charge_hessian did not reject charged 1D/2D structure")
            return
         end if
         deallocate(expected)

         ! EEQBC must support charged 1D and 2D and conserve total charge
         call get_charges(model_eeqbc, mol, error, qvec, dqdr, dqdL, cache=cache)
         if (allocated(error)) return
         if (abs(sum(qvec) - charges(iq)) > 1.0e-12_wp) then
            call test_failed(error, "Charged 1D/2D EEQBC did not conserve total charge")
            return
         end if
         hess = 0.0_wp
         call add_charge_hessian(model_eeqbc, mol, error, weights, hess, cache)
         if (allocated(error)) return
      end do
   end do

   ! Finite-difference derivative verification for charged 1D and 2D EEQBC
   call new(mol, [8, 1, 6], xyz, charge=1.0_wp, lattice=lat, periodic=[.true., .false., .false.])
   call test_numgrad(error, mol, model_eeqbc)
   if (allocated(error)) return
   call test_numsigma(error, mol, model_eeqbc)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model_eeqbc)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model_eeqbc)
   if (allocated(error)) return

   call new(mol, [8, 1, 6], xyz, charge=-0.5_wp, lattice=lat, periodic=[.true., .true., .false.])
   call test_numgrad(error, mol, model_eeqbc)
   if (allocated(error)) return
   call test_numsigma(error, mol, model_eeqbc)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model_eeqbc)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model_eeqbc)
   if (allocated(error)) return

end subroutine test_charged_lowd_contracts

end module test_pbc
