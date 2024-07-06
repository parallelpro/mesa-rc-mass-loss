! *****************************************************
!  p-mode asteroseismology for solor-like oscillators. 
!  save pi modes and p modes.
!  for mesa-r24.03.1 and gyre-7.1 (the version shipped with mesa)
! *****************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use gyre_lib

      implicit none

      ! self-defined variables for the whole module to use
      character(len=6) :: model_number
      character(len=80) :: summary_filename
      integer :: summary_unit
      logical :: flag_gyre, flag_flash

      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         s% other_wind => sc_wind

         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      
      subroutine sc_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         use star_lib, only: star_ptr
         integer, intent(in) :: id
         type (star_info), pointer :: s
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         ! x_ctrl(1) = my_special_param
         w = 2 * 1d-13 * s% x_ctrl(1) * (Lsurf / Lsun) * &
                        (Rsurf / Rsun) * & 
                        (Msurf / Msun)**(-1.) * &
                        (Tsurf / 4000.)**(3.5) * &
                        (1 + (Msun*Rsun**(-2.)) / (4300*Msurf*(Rsurf**(-2.))) )
      end subroutine sc_wind
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         ! if ((s% power_he_burn / (s%power_he_burn + s%power_h_burn)) > 0.001) then
         !    ! stop when center H1 mass drops to specified level
         !    extras_check_model = terminate
         !    write(*, *) 'have reached flash.'
         !    return
         ! end if

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 3
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         integer :: k 
         real(dp) :: vals(n), dr
         logical :: entered_g_mode_cavity
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         names(1) = 'flag_gyre'
         names(2) = 'flag_flash'
         names(3) = 'delta_Pg_outward'

         if (flag_gyre) then
            vals(1) = 1.0
         else
            vals(1) = 0.0
         endif

         if (flag_flash) then
            vals(2) = 1.0
         else
            vals(2) = 0.0
         endif

         vals(3) = 0
         entered_g_mode_cavity = .false.
         do k = s% nz, 2, -1
            dr = s% rmid(k-1) - s% rmid(k)
            if (s% brunt_N2(k) > 0) then
               entered_g_mode_cavity = .true.
               vals(3) = vals(3) + sqrt(s% brunt_N2(k))/s% r(k)*dr
            else
               if (entered_g_mode_cavity) exit
            end if
         end do

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         character(len=50), dimension(2) :: strings_to_replace(4), replace_with_strings(4)
         character(len=50) :: freq_min_radial, freq_max_radial, freq_min_nonradial, freq_max_nonradial
         real(dp) :: width !, obs_Teff_min, obs_Teff_max, obs_L_min, obs_L_max

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! flag_gyre controls whether we run gyre for this step
         ! 1) MS and post-TAMS !.and. ((s%L_He/s%L_nuc_burn_total)<=0.0001)
         ! 3) profile_interval not too large
         flag_gyre = ((s%center_h1 < 0.5) .or. &
         ( ( (s%L_nuc_burn_total/s%L_surf)>=0.999 ) )  ) &
            .and. (s%xtra_old(7) >= s%xtra_old(6)) &
            .and. (s%xtra_old(7) > 0)


         if (flag_gyre) then

            ! 1. Calculate pi modes

            ! 1.1 Modify gyre_template_pi.in and spits out gyre.in
            ! needs to be changed at each step
            ! dynamically change the frequency range
            width = exp(0.9638 * log(s%nu_max) - 1.7145) ! according to a fitting relation
            if (s%nu_max < 2.) then 
               write(freq_min_radial, '(F50.6)')  s%delta_nu * 0.1 ! s%nu_max - width*3
               write(freq_min_nonradial, '(F50.6)')  s%delta_nu * 0.1 ! s%nu_max - width*3
            else
               write(freq_min_radial, '(F50.6)')  s%nu_max - width*3
               write(freq_min_nonradial, '(F50.6)')  s%nu_max - width*3
            endif 

            write(freq_max_radial, '(F50.6)')  s%nu_max + width*3 ! s%nu_max + width*3
            write(freq_max_nonradial, '(F50.6)') s%nu_max + width*3 ! s%nu_max + width*3

            strings_to_replace = (/ string_pad_left('{freq_min_radial}', 50), &
                                    string_pad_left('{freq_max_radial}', 50), &
                                    string_pad_left('{freq_min_nonradial}', 50), & 
                                    string_pad_left('{freq_max_nonradial}', 50) /)
                                 
            replace_with_strings = (/ freq_min_radial, freq_max_radial, freq_min_nonradial, freq_max_nonradial /)

            ! 1.2 Initialize GYRE
            ! see $MESA_DIR/gyre/make/gyre_lib.f90
            call find_and_replace_string_in_file(strings_to_replace, replace_with_strings, 'gyre_template_pi.in', 'gyre.in')
            call gyre_init('gyre.in')

            ! 1.3 Initialize summary file
            ! assign unit numbers
            summary_unit = 100
            write(model_number, '(I0.6)') s%model_number
            summary_filename = trim( 'LOGS/' // trim(adjustl(s%star_history_name)) // '.model' // model_number // '.txt' )
            open(summary_unit, file=summary_filename, status="replace", action="write")
            write(summary_unit, '(5A20, 5A40)') 'id', 'l', 'n_p', 'n_g', 'n_pg', 'Re(freq)', 'Re(freq)_corr', 'E_norm', 'E_p', 'E_g'

            ! 1.4 Run GYRE
            call run_gyre(id, ierr)

            ! 1.5 close file
            close(summary_unit)

         endif

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! Record variables that needs for saving profiles
         ! use a custom frequency, save 50 profiles per 1 dex in Dnu, 
         ! or 0.2 profiles per 1K in Teff, or 50 profiles per 1dex in log10(L/Lsun)
         ! d(N_m)/d(log10(Dnu)) = d(N_m)/d(N_p) * d(N_p)/d(log10(Dnu)) = gyre_interval * 200/1dex
         ! d(N_m)/d(Teff) = d(N_m)/d(N_p) * d(N_p)/d(Teff) = gyre_interval * 0.1/1K
         ! d(N_m)/d(log10(L)) = d(N_m)/d(N_p) * d(N_p)/d(log10(L)) = gyre_interval * 100/1dex
         s%xtra(1) = log10(s%delta_nu)
         s%xtra(2) = s%Teff
         s%xtra(3) = log10(s%L_phot)
         s%xtra(4) = s%center_h1 
         s%xtra(5) = s%center_he4
         if ((s%center_he4 > 0.001) .and. (s%center_h1 < 0.001)) then 
            s%xtra(6) = max( min(int( 1./abs(s%xtra(1)-s%xtra_old(1)) / 100. ), &
                                 int( 1./abs(s%xtra(2)-s%xtra_old(2)) / 0.005 ), &
                                 int( 1./abs(s%xtra(3)-s%xtra_old(3)) / 100. ), &
                                 int( 1./abs(s%xtra(5)-s%xtra_old(5)) / 200. ) ) &
                           , 1) ! gyre_interval
         else 
            s%xtra(6) = max( min(int( 1./abs(s%xtra(1)-s%xtra_old(1)) / 100. ), &
                                 int( 1./abs(s%xtra(2)-s%xtra_old(2)) / 0.005 ), &
                                 int( 1./abs(s%xtra(3)-s%xtra_old(3)) / 100. ) ) &
                           , 1) ! gyre_interval
         endif

         if (flag_gyre) then
            ! print *, 'Modifying xtra(5) ', flag_gyre
            s%xtra(7) = 1 ! actual_gyre_interval
         else 
            s%xtra(7) = s%xtra_old(7) + 1 ! actual_gyre_interval
         endif

         flag_flash = ((s% power_he_burn / (s%power_he_burn + s%power_h_burn)) > 0.001) .and. (s%center_he4 > 0.93)

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      
      subroutine run_gyre(id, ierr)

         integer, intent(in)  :: id
         integer, intent(out) :: ierr

         real(dp), allocatable :: global_data(:)
         real(dp), allocatable :: point_data(:,:)
         integer               :: ipar(0)
         real(dp)              :: rpar(0)

         ! Set constants
         call gyre_set_constant('G_GRAVITY', standard_cgrav)
         call gyre_set_constant('C_LIGHT', clight)
         call gyre_set_constant('A_RADIATION', crad)
         call gyre_set_constant('M_SUN', Msun)
         call gyre_set_constant('R_SUN', Rsun)
         call gyre_set_constant('L_SUN', Lsun)
         call gyre_set_constant('GYRE_DIR', trim(mesa_dir)//'/gyre/gyre')

         ! print *, 'inside run_gyre'
         ! Pass model data to GYRE
         call star_get_pulse_data(id, 'GYRE', .false., .true., .true., &
              global_data, point_data, ierr)
              ! see $MESA_DIR/star/private/pulse.f90 and  $MESA_DIR/star/private/pulse_gyre.f90 
              ! The function gets structure from MESA and save in global_data and point_data
              ! Inputs of subroutine star_get_pulse_data
              ! id, data_format, add_center_point, keep_surface_point, add_atmosphere, 
              ! global_data: global_data(1) = m_outer. global_data(2) = r_outer, global_data(3) = s%L(1)
              ! point_data: a 2d array for stellar structure
              ! ierr
         if (ierr /= 0) then
            print *,'Failed when calling star_get_pulse_data'
            return
         end if

         ! Calculate relavent variables from point_data
         ! see $MESA_DIR/gyre/make/gyre_lib.f90 and $MESA_DIR/gyre/make/gyre_mesa_file.f90
         call gyre_set_model(global_data, point_data, 101)


         ! Run GYRE to get modes
         ! see $MESA_DIR/gyre/make/gyre_lib.f90
         ! the first input is the l-degree.
         ! users can pass process_mode as a subroutine to do some processing (outputing results in our case) 
         ! after each mode if found. 
         call gyre_get_modes(0, process_mode_radial, ipar, rpar)
         call gyre_get_modes(1, process_mode_nonradial, ipar, rpar)
         call gyre_get_modes(2, process_mode_nonradial, ipar, rpar)
         call gyre_get_modes(3, process_mode_nonradial, ipar, rpar)
         
         ! Deallocate memories
         ! see $MESA_DIR/gyre/make/gyre_lib.f90
         call gyre_final()

      contains

         subroutine process_mode_radial (md, ipar, rpar, retcode)
         
            type(mode_t), intent(in) :: md
            integer, intent(inout)   :: ipar(:)
            real(dp), intent(inout)  :: rpar(:)
            integer, intent(out)     :: retcode

            ! Save results to the summary file
            write(summary_unit, '(5I20, 5ES40.20)') md%id, md%l, md%n_p, md%n_g, md%n_pg, real(md%freq('UHZ')), real(md%freq('UHZ')), md%E_norm(), md%E_p(), md%E_g()

            retcode = 0
         end subroutine process_mode_radial

         subroutine process_mode_nonradial (md, ipar, rpar, retcode)
         
            type(mode_t), intent(in) :: md
            integer, intent(inout)   :: ipar(:)
            real(dp), intent(inout)  :: rpar(:)
            integer, intent(out)     :: retcode

            integer :: j
            real(dp) :: lambda, integral1, integral2, Rpipi
            real(wp), allocatable :: x(:), rho(:), N2(:), xi_r(:), xi_h(:)
            logical, allocatable :: mask(:)

            type (star_info), pointer :: s
            call star_ptr(id, s, ierr)

            ! We will output two types of modes
            ! the uncorrected modes, and the first-order corrected modes (Ong & Basu 2020)
   
            ! Pull structure data and mode data
            allocate(x(size(point_data, 2)))
            allocate(rho(size(point_data, 2)))
            allocate(N2(size(point_data, 2)))
            allocate(xi_r(size(point_data, 2)))
            allocate(xi_h(size(point_data, 2)))
            allocate(mask(size(point_data, 2)))  ! allocate(mask(size( s% nz)))

            lambda = md%l * (md%l+1)

            ! see also $MESA_DIR/gyre/gyre/src/model/gyre_mesa_file.fpp
            x = point_data(1,:) / (10.0**s%log_surface_radius * Rsun) !point_data(1, size(point_data, 2)) ! fractional radius
            do j = 1, size(point_data, 2)
               xi_r(j) = md%xi_r(j)
               xi_h(j) = md%xi_h(j)
            enddo

            rho = point_data(6,:) 
            N2 = point_data(8,:) 

            mask = (N2 > 0) .and. (x < 0.9)

            ! performing the integration only on the required indices
            integral1 = trapz(apply_mask(x**2.0 * rho * xi_r**2.0 * N2, mask), apply_mask(x, mask))
            integral2 = trapz(x**2.0 * rho * (xi_r**2.0 + lambda * xi_h**2.0), x)
            
            Rpipi = integral1 / integral2 / (2.0 * 3.1415926535/1e6)**2.0 ! from angular frequency (s^-1) to frequency (muHz)

            ! Save results to the summary file
            write(summary_unit, '(5I20, 5ES40.20)') md%id, md%l, md%n_p, md%n_g, md%n_pg, real(md%freq('UHZ')),  (real(md%freq('UHZ'))**2.0 + Rpipi)**0.5, md%E_norm(), md%E_p(), md%E_g()

            retcode = 0
         end subroutine process_mode_nonradial
      end subroutine run_gyre


      subroutine find_and_replace_string_in_file(strings_to_replace, replace_with_strings, input_filename, output_filename)
         
         character(len=50), intent(in) :: strings_to_replace(:), replace_with_strings(:)
         character(len=*), intent(in) :: input_filename, output_filename
         character(500) :: input_line, temp_line
         integer :: iunit, ounit, ios, pos, i

         ! Ensure the two arrays have the same size
         if (size(strings_to_replace) /= size(replace_with_strings)) then
            stop "Error: The size of stringToReplace and replacementString arrays do not match!"
         end if

         ! Assign unit numbers
         iunit = 10
         ounit = 20

         ! Open input and output files
         open(iunit, file=input_filename, status='old', action='read')
         open(ounit, file=output_filename, status='replace', action='write')

         ! Loop through each line of the input file
         do
            read(iunit, '(A)', iostat=ios) input_line
            if (ios /= 0) exit  ! Exit the loop if we've reached the end of the file or on error

            ! Replace substrings for each pair of strings_to_replace and replace_with_strings
            do i = 1, size(strings_to_replace)
               do
                  pos = index(input_line, trim(adjustl(strings_to_replace(i))))
                  if (pos == 0) exit
                  temp_line = input_line(1:pos-1) // trim(adjustl(replace_with_strings(i))) // &
                           input_line(pos+len(trim(adjustl(strings_to_replace(i)))):)
                  input_line = temp_line
               end do
            end do

            ! Write the modified line to the output file
            write(ounit, '(A)') trim(input_line)
         end do

         ! Close files
         close(iunit)
         close(ounit)
    
      end subroutine find_and_replace_string_in_file

      function string_pad_left(input_string, desired_length) result(padded_string)
         character(len=*), intent(in) :: input_string
         integer, intent(in) :: desired_length
         character(len=desired_length) :: padded_string
         padded_string = repeat(' ', desired_length - len(input_string)) // input_string
      end function string_pad_left

      function string_pad_right(input_string, desired_length) result(padded_string)
         character(len=*), intent(in) :: input_string
         integer, intent(in) :: desired_length
         character(len=desired_length) :: padded_string
         padded_string = input_string // repeat(' ', desired_length - len(input_string))
      end function string_pad_right

      function trapz(y, x) result(integral)
         real(dp), dimension(:), intent(in) :: y, x
         real(dp) :: integral
         integer :: n

         ! integral = 0.0
         ! do i = 1, size(x) - 1
         !    integral = integral + 0.5_WP * (x(i+1) - x(i)) * (y(i+1) + y(i))
         ! end do
         n = size(x)
         integral = sum(0.5_WP*(y(2:) + y(:n-1))*(x(2:) - x(:n-1)))

      end function trapz

      function apply_mask(input_array, mask) result(output_array)
         real(dp), dimension(:), intent(in) :: input_array
         logical, dimension(:), intent(in) :: mask
         real(dp), dimension(:), allocatable :: output_array

         integer :: i, counter

         ! Check for size mismatch
         if (size(input_array) /= size(mask)) then
            print *, "Error: Size mismatch between input array and mask."
            return
         end if

         counter = count(mask)
         allocate(output_array(counter))

         counter = 0
         do i = 1, size(input_array)
            if (mask(i)) then
                  counter = counter + 1
                  output_array(counter) = input_array(i)
            end if
         end do
      end function apply_mask

      end module run_star_extras