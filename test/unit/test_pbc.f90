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
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge, only: get_charges, add_charge_hessian
   use multicharge_ewald, only: get_ewald_dir_trans_1d, get_ewald_rec_trans_1d, &
      & get_ewald_dir_trans_2d, get_ewald_rec_trans_2d
   use multicharge_model, only: mchrg_model_type
   use multicharge_model_cache, only: cache_container
   use multicharge_model_eeq, only: eeq_cache
   use multicharge_model_eeqbc, only: eeqbc_model, eeqbc_cache
   use multicharge_model_type, only: get_dir_trans
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
      & new_unittest("degenerate-lattice-error", test_degenerate_lattice_error) &
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

   call get_dir_trans(mol, trans, cutoff)

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

   call get_dir_trans(mol, trans, cutoff)

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

   call get_dir_trans(mol, trans, cutoff)

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

   call get_dir_trans(mol, trans, cutoff)

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
         call model%get_xvec(mol, cache, xvecr)

         ! Left-hand side
         xvecl(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat, :) = 0.5_wp * (xvecr(:) - xvecl(:)) / step
      end do
   end do lp

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
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

   call get_dir_trans(mol, trans, cutoff)

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

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
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

   call get_dir_trans(mol, trans, cutoff)

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
         call model%get_coulomb_matrix(mol, cache, amatr)

         ! Left-hand side
         amatl(:, :) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatl)

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

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   call model%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)

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

   call get_dir_trans(mol, trans, cutoff)

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
         call model%get_coulomb_matrix(mol, cache, amatr)
         if (allocated(error)) exit lp

         amatl(:, :) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatl)
         if (allocated(error)) exit lp

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

   call model%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)
   if (allocated(error)) return

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

   call get_dir_trans(mol, trans, cutoff)

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

   call get_dir_trans(mol, trans, cutoff)

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
   real(wp), parameter :: ref(40) = [&
      & 3.47274007973765E-1_wp, 3.47273017259661E-1_wp, 3.47276806231462E-1_wp, &
      & 3.47267235121650E-1_wp, 3.47288751410308E-1_wp, 3.47296840650318E-1_wp, &
      & 3.47299868243747E-1_wp, 3.47299359754750E-1_wp, 3.45885710186856E-1_wp, &
      & 3.45874246032015E-1_wp, 3.45888242047875E-1_wp, 3.45877451600398E-1_wp, &
      & 3.45902365123333E-1_wp, 3.45902162041418E-1_wp, 3.45900336974539E-1_wp, &
      & 3.45903770121699E-1_wp, 3.58996736955583E-1_wp, 3.58987991603151E-1_wp, &
      & 3.59000859720791E-1_wp, 3.58990960169615E-1_wp, 3.58990422147650E-1_wp, &
      & 3.58997842028106E-1_wp, 3.58997124263551E-1_wp, 3.59001728635141E-1_wp, &
      &-5.86637659146611E-1_wp, -5.86601835343945E-1_wp, -5.86625078062709E-1_wp, &
      &-5.86579999312624E-1_wp, -5.86654394920376E-1_wp, -5.86669360667441E-1_wp, &
      &-5.86667481028248E-1_wp, -5.86674873278027E-1_wp, -4.65529092171338E-1_wp, &
      &-4.65527713392387E-1_wp, -4.65546450764424E-1_wp, -4.65540905286785E-1_wp, &
      &-4.65529860758785E-1_wp, -4.65533757413217E-1_wp, -4.65527683275370E-1_wp, &
      &-4.65527691475100E-1_wp]

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
   ! Updated for the mctc-lib Wigner-Seitz image tolerance in the finite Ewald sum.
   real(wp), parameter :: ref(24) = [&
      & 4.01878374787191E-1_wp, 4.01885037899597E-1_wp, 4.01836688390911E-1_wp, &
      & 4.01859670984182E-1_wp, 4.30132953804280E-1_wp, 4.30105983229582E-1_wp, &
      & 4.30085942904623E-1_wp, 4.30109219580760E-1_wp, 1.90400205083764E-1_wp, &
      & 1.90401834703106E-1_wp, 1.90404142905890E-1_wp, 1.90407459232625E-1_wp, &
      & 3.33868772086876E-1_wp, 3.33849119518316E-1_wp, 3.33876414218980E-1_wp, &
      & 3.33893467423880E-1_wp, -8.35636285533486E-1_wp, -8.35615627485463E-1_wp, &
      &-8.35543385304371E-1_wp, -8.35590189215555E-1_wp, -6.31079470991667E-1_wp, &
      &-6.31060026823144E-1_wp, -6.31085555821538E-1_wp, -6.31082095933280E-1_wp]

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

         ! Full derivative validation for TFF (neutral) and FTF (charged)
         if ((iq == 1 .and. imask == 1) .or. (iq == 2 .and. imask == 2)) then
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

         ! Full derivative validation for TTF (neutral) and FTT (charged)
         if ((iq == 1 .and. imask == 1) .or. (iq == 2 .and. imask == 3)) then
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

      ! Full derivative validation for TTF and TFT
      if (imask <= 2) then
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

      ! Full derivative validation for TFF and FTF
      if (imask <= 2) then
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
   call get_ewald_dir_trans_1d(a1, 0.06_wp, 0.06_wp, dtrans, error=error)
   if (allocated(error)) return
   if (size(dtrans, 2) <= 81) then
      call test_failed(error, "1D direct-space Ewald remains capped at 40 shells")
      return
   end if

   a1 = [45.0_wp, 0.0_wp, 0.0_wp]
   call get_ewald_rec_trans_1d(a1, 0.50_wp, rtrans, len1, u_hat, error=error)
   if (allocated(error)) return
   if (size(rtrans, 2) <= 80) then
      call test_failed(error, "1D reciprocal-space Ewald remains capped at 40 shells")
      return
   end if

   a1 = [3.0_wp, 0.0_wp, 0.0_wp]
   a2 = [0.0_wp, 3.0_wp, 0.0_wp]
   call get_ewald_dir_trans_2d(a1, a2, 0.08_wp, 0.08_wp, dtrans, error=error)
   if (allocated(error)) return
   if (size(dtrans, 2) <= 1681) then
      call test_failed(error, "2D direct-space Ewald remains capped at 20 shells")
      return
   end if

   a1 = [30.0_wp, 0.0_wp, 0.0_wp]
   a2 = [0.0_wp, 30.0_wp, 0.0_wp]
   call get_ewald_rec_trans_2d(a1, a2, 0.45_wp, rtrans, area, n_hat, error=error)
   if (allocated(error)) return
   if (size(rtrans, 2) <= 1680) then
      call test_failed(error, "2D reciprocal-space Ewald remains capped at 20 shells")
      return
   end if

   ! 2. Verify explicit resource-limit error reporting (no silent truncation)
   call get_ewald_dir_trans_1d([2.0_wp, 0.0_wp, 0.0_wp], 0.06_wp, 0.06_wp, dtrans, &
      & max_trans=10, error=expected)
   if (.not. allocated(expected)) then
      call test_failed(error, "1D direct-space resource limit was not enforced")
      return
   end if
   deallocate(expected)

   call get_ewald_rec_trans_1d([45.0_wp, 0.0_wp, 0.0_wp], 0.50_wp, rtrans, len1, u_hat, &
      & max_trans=10, error=expected)
   if (.not. allocated(expected)) then
      call test_failed(error, "1D reciprocal-space resource limit was not enforced")
      return
   end if
   deallocate(expected)

   call get_ewald_dir_trans_2d([3.0_wp, 0.0_wp, 0.0_wp], [0.0_wp, 3.0_wp, 0.0_wp], &
      & 0.08_wp, 0.08_wp, dtrans, max_trans=10, error=expected)
   if (.not. allocated(expected)) then
      call test_failed(error, "2D direct-space resource limit was not enforced")
      return
   end if
   deallocate(expected)

   call get_ewald_rec_trans_2d([30.0_wp, 0.0_wp, 0.0_wp], [0.0_wp, 30.0_wp, 0.0_wp], &
      & 0.45_wp, rtrans, area, n_hat, max_trans=10, error=expected)
   if (.not. allocated(expected)) then
      call test_failed(error, "2D reciprocal-space resource limit was not enforced")
      return
   end if
   deallocate(expected)

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
      end select

      call new_eeq2019_model(mol, model, error)
      if (allocated(error)) return
      if (ig == 1) then
         ! Broaden Gaussian charge widths on the short-period 1D wire
         model%rad = 2.2_wp * model%rad
      end if

      call get_dir_trans(mol, trans, 25.0_wp)
      call model%ncoord%get_coordination_number(mol, trans, cn)
      call model%local_charge(mol, trans, qloc)
      call model%update(mol, cache, cn, qloc)
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

!> Regression test for DFT-D4 issue #358: neutral water-helium complex at large
!> separation in 0D, 1D, and 2D. Verifies that EEQBC screens out the long-range
!> inter-fragment capacitance and Coulomb coupling, suppressing distance-dependent
!> asymptotic charge transfer and charge gradients in contrast to unscreened EEQ.
subroutine test_eeqbc_water_helium_358(error)
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol_r1, mol_r2
   class(mchrg_model_type), allocatable :: model_eeq, model_eeqbc
   type(cache_container) :: cache_eeq, cache_eeqbc
   real(wp) :: xyz_r1(3, 4), xyz_r2(3, 4), lat(3, 3)
   real(wp) :: q_eeq1(4), q_eeq2(4), dqdr_eeq2(3, 4, 4), dqdL_eeq2(3, 3, 4)
   real(wp) :: q_bc1(4), q_bc2(4), dqdr_bc2(3, 4, 4), dqdL_bc2(3, 3, 4)
   real(wp) :: amat_eeq(5, 5), amat_bc(5, 5), cmat_bc(5, 5)
   logical :: pmask(3)
   integer :: idim

   xyz_r1 = reshape([ &
      &  0.0000000000_wp,  0.0000000000_wp,  0.2214080000_wp, &
      &  0.0000000000_wp,  1.4309280000_wp, -0.8856320000_wp, &
      &  0.0000000000_wp, -1.4309280000_wp, -0.8856320000_wp, &
      &  0.0000000000_wp,  0.0000000000_wp, 35.0000000000_wp], [3, 4])
   xyz_r2 = xyz_r1
   xyz_r2(:, 4) = [0.0_wp, 0.0_wp, 60.0_wp]

   lat = reshape([ &
      & 12.0_wp,  0.0_wp,   0.0_wp, &
      &  0.0_wp, 12.0_wp,   0.0_wp, &
      &  0.0_wp,  0.0_wp, 140.0_wp], [3, 3])

   do idim = 0, 2
      select case (idim)
      case (0); pmask = [.false., .false., .false.]
      case (1); pmask = [.true.,  .false., .false.]
      case (2); pmask = [.true.,  .true.,  .false.]
      end select

      if (idim == 0) then
         call new(mol_r1, [8, 1, 1, 2], xyz_r1, charge=0.0_wp)
         call new(mol_r2, [8, 1, 1, 2], xyz_r2, charge=0.0_wp)
      else
         call new(mol_r1, [8, 1, 1, 2], xyz_r1, charge=0.0_wp, &
            & lattice=lat, periodic=pmask)
         call new(mol_r2, [8, 1, 1, 2], xyz_r2, charge=0.0_wp, &
            & lattice=lat, periodic=pmask)
      end if

      if (allocated(model_eeq)) deallocate(model_eeq)
      if (allocated(model_eeqbc)) deallocate(model_eeqbc)
      call new_eeq2019_model(mol_r2, model_eeq, error)
      if (allocated(error)) return
      call new_eeqbc2025_model(mol_r2, model_eeqbc, error)
      if (allocated(error)) return

      call get_charges(model_eeq, mol_r1, error, q_eeq1)
      if (allocated(error)) return
      call get_charges(model_eeq, mol_r2, error, q_eeq2, dqdr_eeq2, dqdL_eeq2, cache=cache_eeq)
      if (allocated(error)) return
      call model_eeq%get_coulomb_matrix(mol_r2, cache_eeq, amat_eeq)

      call get_charges(model_eeqbc, mol_r1, error, q_bc1)
      if (allocated(error)) return
      call get_charges(model_eeqbc, mol_r2, error, q_bc2, dqdr_bc2, dqdL_bc2, cache=cache_eeqbc)
      if (allocated(error)) return
      call model_eeqbc%get_coulomb_matrix(mol_r2, cache_eeqbc, amat_bc)
      select type (raw => cache_eeqbc%raw)
      type is (eeqbc_cache)
         cmat_bc = raw%cmat
      end select

      ! Unscreened EEQ has long-range 1/R coupling and distance-dependent charge transfer
      if (maxval(abs(amat_eeq(1:3, 4))) < 1.0e-2_wp &
         & .or. abs(q_eeq2(4) - q_eeq1(4)) < 1.0e-4_wp) then
         call test_failed(error, "Expected long-range 1/R coupling in unscreened EEQ")
         return
      end if

      ! EEQBC screens out inter-fragment capacitance and Coulomb coupling at large R,
      ! suppressing distance-dependent charge variation and He charge response to zero
      if (maxval(abs(cmat_bc(1:3, 4))) > 1.0e-14_wp &
         & .or. maxval(abs(amat_bc(1:3, 4))) > 1.0e-14_wp) then
         call test_failed(error, "EEQBC failed to screen long-range H2O-He coupling")
         return
      end if
      if (maxval(abs(q_bc2 - q_bc1)) > 1.0e-12_wp &
         & .or. maxval(abs(dqdr_bc2(:, 4, :))) > 1.0e-14_wp) then
         call test_failed(error, "EEQBC failed to suppress asymptotic charge response on He")
         return
      end if
   end do

end subroutine test_eeqbc_water_helium_358

end module test_pbc
