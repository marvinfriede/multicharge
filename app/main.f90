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

program main
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, read_structure, filetype, get_filetype
   use multicharge, only : mchrg_model_type, new_eeq2019_model, &
      & write_ascii_model, write_ascii_properties, write_ascii_results, &
      & get_coordination_number, get_covalent_rad, get_lattice_points, &
      & get_multicharge_version
   implicit none
   character(len=*), parameter :: prog_name = "multicharge"

   character(len=:), allocatable :: input
   integer, allocatable :: input_format
   type(error_type), allocatable :: error
   type(structure_type) :: mol
   type(mchrg_model_type) :: model
   logical :: grad
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), rcov(:), trans(:, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :), hessian(:, :, :, :)
   real(wp), allocatable :: qvec(:), dqdr(:, :, :), dqdL(:, :, :)

   integer :: iat, ix
   real(wp), parameter :: step = 1.0e-4_wp
   type(structure_type) :: displ
   real(wp), allocatable :: el(:), er(:)
   real(wp), allocatable :: gl(:, :), gr(:, :), sl(:, :), sr(:, :)
   real(wp), allocatable :: qvec_l(:), dqdr_l(:, :, :), dqdL_l(:, :, :)
   real(wp), allocatable :: qvec_r(:), dqdr_r(:, :, :), dqdL_r(:, :, :)

   real(wp), allocatable :: dcndr_l(:, :, :), dcndL_l(:, :, :)
   real(wp), allocatable :: dcndr_r(:, :, :), dcndL_r(:, :, :)

   call get_arguments(input, input_format, grad, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (input == "-") then
      if (.not.allocated(input_format)) input_format = filetype%xyz
      call read_structure(mol, input_unit, input_format, error)
   else
      call read_structure(mol, input, error, input_format)
   end if
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

  ! mol%charge = -1.0_wp

   call new_eeq2019_model(mol, model)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   call write_ascii_model(output_unit, mol, model)

   allocate(cn(mol%nat))
   if (grad) then
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   end if

   allocate(hessian(mol%nat, 3, mol%nat, 3))

   rcov = get_covalent_rad(mol%num)
   call get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL, cut=cn_max)

   allocate(energy(mol%nat), qvec(mol%nat))
   energy(:) = 0.0_wp
   if (grad) then
      allocate(gradient(3, mol%nat), sigma(3, 3))
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp

      allocate(dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat))
      dqdr(:, :, :) = 0.0_wp
      dqdL(:, :, :) = 0.0_wp
   end if

   call model%solve(mol, cn, dcndr, dcndL, energy, gradient, sigma, &
      & qvec, dqdr, dqdL)

   call write_ascii_properties(output_unit, mol, model, cn, qvec)
   call write_ascii_results(output_unit, mol, energy, gradient, sigma)

   write (*, *) "energy"
   write(*, '(SP,es23.16e2,",")') energy
   write (*, *) ""

      write (*, *) "qvec"
   write(*, '(SP,es23.16e2,",")') qvec
   write (*, *) ""

   write (*, *) "grad"
   write(*, '(*(6x,SP,"[",3(es23.16e2, "":, ","), "],", /))', advance='no') gradient
   write (*, *) ",]"

   write (*, *) "dqdr"
   write(*, '(SP,es23.16e2,",")') dqdr
   write (*, *) ""



   hessian(:, :, :, :) = 0.0_wp
   
   allocate(dcndr_l(3, mol%nat, mol%nat), dcndL_l(3, 3, mol%nat))
   allocate(dcndr_r(3, mol%nat, mol%nat), dcndL_r(3, 3, mol%nat))

   
   allocate(gl(3, mol%nat), gr(3, mol%nat), sl(3, 3), sr(3, 3))


   allocate(el(mol%nat), er(mol%nat))

   
   allocate(qvec_l(mol%nat), qvec_r(mol%nat))



   allocate(dqdr_l(3, mol%nat, mol%nat), dqdL_l(3, 3, mol%nat))
   allocate(dqdr_r(3, mol%nat, mol%nat), dqdL_r(3, 3, mol%nat))

   return

   displ = mol
   do iat = 1, mol%nat
      do ix = 1, 3
         displ%xyz(ix, iat) = mol%xyz(ix, iat) + step
          
          dqdr_l(:, :, :) = 0.0_wp
          dqdr_r(:, :, :) = 0.0_wp
          dqdL_l(:, :, :) = 0.0_wp
          dqdL_r(:, :, :) = 0.0_wp
           qvec_l(:) = 0.0_wp
          qvec_r(:) = 0.0_wp

          el(:) = 0.0_wp
          er(:) = 0.0_wp

              gl(:, :) = 0.0_wp
          gr(:, :) = 0.0_wp
          sl(:, :) = 0.0_wp
          sr(:, :) = 0.0_wp
          
          dcndr_l(:, :, :) = 0.0_wp
          dcndr_r(:, :, :) = 0.0_wp
          dcndL_l(:, :, :) = 0.0_wp
          dcndL_r(:, :, :) = 0.0_wp

         call get_lattice_points(displ%periodic, displ%lattice, cutoff, trans)
         call get_coordination_number(displ, trans, cutoff, rcov, cn, dcndr_l, dcndL_l, cut=cn_max)
         call model%solve(displ, cn, dcndr_l, dcndL_l, el, gl, sl, &
          & qvec_l, dqdr_l, dqdL_l)


         displ%xyz(ix, iat) = mol%xyz(ix, iat) - step
         call get_lattice_points(displ%periodic, displ%lattice, cutoff, trans)
         call get_coordination_number(displ, trans, cutoff, rcov, cn, dcndr_r, dcndL_r, cut=cn_max)
         call model%solve(displ, cn, dcndr_r, dcndL_r, er, gr, sr, &
          & qvec_r, dqdr_r, dqdL_r)

         displ%xyz(ix, iat) = mol%xyz(ix, iat)
         hessian(:, :, iat, ix) = (transpose(gl) - transpose(gr)) / (2 * step)
      end do
   end do


   write(*, '(SP,es23.16e2,",")') hessian

contains


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      "Electronegativity equilibration model for atomic charges and", &
      "higher multipole moments", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "-i, --input <format>", "Hint for the format of the input file", &
      "--grad", "Evaluate molecular gradient and virial", &
      "--version", "Print program version and exit", &
      "--help", "Show this help message"

   write(unit, '(a)')

end subroutine help


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_multicharge_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


subroutine get_arguments(input, input_format, grad, error)

   !> Input file name
   character(len=:), allocatable :: input

   !> Input file format
   integer, allocatable, intent(out) :: input_format

   !> Evaluate gradient
   logical, intent(out) :: grad

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   character(len=:), allocatable :: arg

   grad = .false.
   iarg = 0
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      select case(arg)
      case("-help", "--help")
         call help(output_unit)
         stop
      case("-version", "--version")
         call version(output_unit)
         stop
      case default
         if (.not.allocated(input)) then
            call move_alloc(arg, input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case("-i", "-input", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         input_format = get_filetype("."//arg)
      case("-grad", "--grad")
         grad = .true.
      end select
   end do

   if (.not.allocated(input)) then
      if (.not.allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_arguments

end program main
