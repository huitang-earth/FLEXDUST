!**********************************************************************
! Copyright 2016,2017,2018                                            *
! Christine Groot Zwaaftink                                           *
! cgz@nilu.no                                                         *
!                                                                     *
! modified by Christian Maurer                                        *
! christian.maurer@zamg.ac.at                                         *
! 201807                                                              *
!                                                                     *
! This file is part of FLEXDUST.                                      *
!                                                                     *
! FLEXDUST is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXDUST is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXDUST.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************
subroutine writeRELEASEfile(filename, typeSizeDistr, particlesPerTonDust, Junge_index, emission, &
    timestep, juldate, lon_left, nx_out, lat_bottom, ny_out, dxdy_degr, &
    release_dxdy_step, start_date_day, releaseDays, totalEmission, totalParticles, &
    minMassWrite)

    use par_mod
    use com_mod

    implicit none
    !***********************************************************************
    character(*):: filename
    real        :: particlesPerTonDust, dxdy_degr, lon_left, lat_bottom, combEmission, totalEmission, releaseDays
    integer     :: timestep, i, typeSizeDistr, numberSpecies, Junge_index, xx, yy
    integer     :: nx_out, ny_out, ix, iy, date_day, date_hour, date_day_orig, date_hour_orig
    integer     :: date_day_ran, date_hour_ran, date_day_next_ran, date_hour_next_ran
    integer     :: date_day_next_orig, date_hour_next_orig, nr, ii, release_dxdy_step
    integer     :: date_day_next, date_hour_next
    integer     :: start_date_day, startSpecies, current_species, nrPart, stepPart, totalParticles, estPart
    real        :: emission(0:nx_out - 1, 0:ny_out - 1)
    real        :: cum_fract, rand_class, randomNumber, minMassWrite
    real(kind = dp)                         :: juldate, mass_species, lat, lon, start_ran
    logical                                 :: fileRELEASE
    integer, parameter                      :: releaseunit = 333, aboveSfc = 1
    real, parameter                         :: bottomHeight_dustProf = 1.0
    real, parameter                         :: topHeight_dustProf = 100.0
    real(8), dimension(:), allocatable      :: inisize
    real(8), dimension(:,:), allocatable    :: inifrac
    INTEGER :: ALLOC_ERR
    !***********************************************************************

    !First adjust settings to chosen size distribution
    !***********************************************************************
    if (typeSizeDistr .eq. 1) then !type defined in dust_mod.f90
        numberSpecies = 15
        startSpecies = 201
    elseif (typeSizeDistr .eq. 2) then
        numberSpecies = 10
        startSpecies = 301
    elseif (typeSizeDistr .eq. 3) then
        numberSpecies = 10
        startSpecies = 401
    endif

    !allocate size & fraction arrays
    !***********************************************************************
    allocate(inisize(0:numberSpecies - 1), STAT = ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"

    allocate(inifrac(0:numberSpecies - 1, 0:2), STAT = ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    !***********************************************************************

    !Check if the RELEASE file exists
    !************************************************************************
    inquire( file = filename, exist = fileRELEASE)

    !Open release file
    !************************************************************************
    open(unit = releaseunit, file = filename, access = 'append', action = 'write')
    !************************************************************************

    totalEmission = 0

    !Write header in release file if it did not exist yet
    !************************************************************************
    if (.not.fileRELEASE)then

        !Fill 11 lines with something
        !**************************************************************
        do i = 0, 3
            write(releaseunit, *) "****** FLEXPART DUST RELEASES HEADER ******"
        end do
        !**************************************************************
        write(releaseunit, *) "***Some information of dust simulation, NOT read in FLEXPART***"
        write(releaseunit, *) "Date Start:", start_date_day, ", nr of Days:", releaseDays
        write(releaseunit, *) "Domain:", lat_bottom, lat_bottom + ny_out * dxdy_degr, lon_left, lon_left + nx_out * dxdy_degr
        write(releaseunit, *) "Resolution of emission calculation:", dxdy_degr
        write(releaseunit, *) "Resolution of releasefile:", release_dxdy_step * dxdy_degr
        write(releaseunit, *) "Type of size distribution, nr species:", typeSizeDistr, numberSpecies
        write(releaseunit, *) "***Some information of dust simulation, NOT read in FLEXPART***"

        !Write names of species
        !**************************************************************
        write(releaseunit, '(I3)') numberSpecies
        do i = 0, numberSpecies - 1
            write(releaseunit, '(I5)') i + startSpecies
        end do
        !**************************************************************

        !Fill 1 line with something
        !**************************************************************
        write(releaseunit, *) "-----------------------------"
        !**************************************************************
        
        !Preparing for namelist format not working but stopped implementing because size increases
        !write(releaseunit, *) '&RELEASES_CTRL'
        ! write(releaseunit, '(A6XI3A)') 'NSPEC=',numberSpecies,','
        ! do i = 1, numberSpecies-1
        !     if (i.eq.1)then
        !         write (7, "(*(A12G0,:,','))") 'SPECNUM_REL=',i + startSpecies
        !     else
        !         write (7, "(*(G0,:,','))") i + startSpecies
        !     endif
        ! end do        
        ! write(releaseunit, *) '/'



    endif
    !************************************************************************

    !Get information on particle size distribution
    !************************************************************************
    call getSizeDistribution(inisize, inifrac, typeSizeDistr, numberSpecies)
    !************************************************************************

    !Get date for this timestep and end of emission period
    !************************************************************************
    call caldate(juldate-real(real(timestep-1)/24.), date_day_orig, date_hour_orig)
    call caldate(juldate + real(real(1)/24.), date_day_next_orig, date_hour_next_orig)
    !print*, 'RELEASE PERIOD: ', juldate,  date_day_orig, date_hour_orig, date_day_next_orig, date_hour_next_orig

    !************************************************************************

    !Loop through grid and write data to release file
    !************************************************************************
    nr = 0
    stepPart = 0
    do iy = 0, ny_out - 1, release_dxdy_step
        do ix = 0, nx_out - 1, release_dxdy_step

            lat = lat_bottom + real(iy) * dxdy_degr
            lon = lon_left + real(ix) * dxdy_degr

            !reset
            combEmission = 0 

            !gather release of several grid points
            do yy = iy, min(iy + release_dxdy_step - 1, ny_out - 1)
                do xx = ix, min(ix + release_dxdy_step - 1, nx_out - 1)
                    combEmission = combEmission + emission(xx, yy)
                end do
            end do
            
            if (combEmission .gt. minMassWrite) then!emission will be written this time step, set emission to 0
                do yy = iy, min(iy + release_dxdy_step - 1, ny_out - 1)
                    do xx = ix, min(ix + release_dxdy_step - 1, nx_out - 1)
                        emission(xx, yy) = 0;
                    end do
                end do
                totalEmission = totalEmission + combEmission
            endif

            !If there was any suspension, write to release file
            !*************************************************************
            if (combEmission .gt. minMassWrite)then

                !Nr particles indication
                estPart = max(1, min(1000, int(particlesPerTonDust * combEmission/1000.)))

                call flush(5)
                if (estPart > 25)then
                    
                    !Loop through complete size distribution if there are many particles emitted (otherwise pick one size class)
                    !*********************************************************
                    do i = 0, numberSpecies - 1
                        !Calculate mass for this particle size
                        !*************************************************
                        mass_species = inifrac(i, Junge_index) * combEmission;
                        !Determine number of particles
                        !*************************************************
                        if (inisize(i) .lt. 15.0)then
                            nrPart = max(1, min(2000, int(particlesPerTonDust * mass_species/1000.)))
                        else
                            nrPart = max(1, min(1000, int(particlesPerTonDust * mass_species/1000.)))
                        endif
                        stepPart = stepPart + nrPart
                        
                        !If less than 8 particles, choose random release time within release period that covers 10% of release period
                        !This is done because depending on settings in FLEXPART it may start all particles at once at the beginning of
                        !a release if there are few particles.
                        if (nrPart .lt. 8)then
                            call RANDOM_NUMBER(randomNumber)
                            start_ran=juldate-real(real(timestep-1)/24.) + (randomNumber * real(timestep))/24. !Random start time within release period
                            call caldate(start_ran, date_day_ran, date_hour_ran) !Start of release in right format                           
                            call caldate(MIN(start_ran+real(real(0.10*timestep))/24.,juldate + real(real(1)/24.)), & !End of random release time may not succeed end of release period
                            date_day_next_ran, date_hour_next_ran)
                            
                            date_day = date_day_ran
                            date_hour = date_hour_ran
                            date_day_next = date_day_next_ran
                            date_hour_next = date_hour_next_ran

                        else
                            date_day = date_day_orig
                            date_hour = date_hour_orig
                            date_day_next = date_day_next_orig
                            date_hour_next = date_hour_next_orig
                        endif

                        !write release info
                        !*************************************************
                        write(releaseunit, '(I8, I07.6)') date_day, date_hour
                        write(releaseunit, '(I8, I07.6)') date_day_next, date_hour_next
                        write(releaseunit, '(F9.4)') lon
                        write(releaseunit, '(F9.4)') lat
                        write(releaseunit, '(F9.4)') lon + real(dxdy_degr * release_dxdy_step)
                        write(releaseunit, '(F9.4)') lat + real(dxdy_degr * release_dxdy_step)
                        write(releaseunit, '(I9)') aboveSfc
                        write(releaseunit, '(F10.3)') bottomHeight_dustProf
                        write(releaseunit, '(F10.3)') topHeight_dustProf
                        write(releaseunit, '(I9)') nrPart

                        !write mass at this point for all species
                        !*************************************************
                        do ii = 0, numberSpecies - 1
                            if (ii .eq. i) then
                                write(releaseunit, '(E10.4)') mass_species
                            else
                                write(releaseunit, '(E10.4)') 0.0
                            endif
                        end do

                        !Finish with name of this particular release
                        !*************************************************
                        write(releaseunit, '(A14,I8,I06.6,A1,I0.4,A8,I3)') 'Dust_location_', &
                        date_day, date_hour, '_', nr, '_species', startSpecies + i
                        write(releaseunit, *) '-----------------------------'
                        nr = nr + 1
                    end do

                else
                   
                    !Not much mass release, gather in 1 particle size, pick random size from size distribution
                    !*************************************************
                    call RANDOM_NUMBER(rand_class)
                    current_species = 0
                    cum_fract = inifrac(current_species, Junge_index)
                    do while (rand_class .gt. cum_fract .and. current_species .lt.9)
                        !print*, rand_class, current_species, cum_fract
                        current_species = current_species + 1
                        cum_fract = cum_fract + inifrac(current_species, Junge_index)
                    end do

                    mass_species = combEmission

                    !Get number of particles
                    !*************************************************
                    if (inisize(current_species) .lt. 15.0)then
                        nrPart = max(1, min(2000, int(particlesPerTonDust * mass_species/1000.)))
                    else
                        nrPart = max(1, min(1000, int(particlesPerTonDust * mass_species/1000.)))
                    endif
                    stepPart = stepPart + nrPart

                    !If less than 8 particles, choose random release time within release period
                    if (nrPart .lt. 8)then
                        call RANDOM_NUMBER(randomNumber)
                        start_ran=juldate-real(real(timestep-1)/24.) + (randomNumber * real(timestep))/24. !Random start time within release period
                        call caldate(start_ran, date_day_ran, date_hour_ran) !Start of release in right format                           
                        call caldate(MIN(start_ran+real(real(0.10*timestep))/24.,juldate + real(real(1)/24.)), & !End of random release time may not succeed end of release period
                        date_day_next_ran, date_hour_next_ran)
                     
                        date_day = date_day_ran
                        date_hour = date_hour_ran
                        date_day_next = date_day_next_ran
                        date_hour_next = date_hour_next_ran
                    else
                        date_day = date_day_orig
                        date_hour = date_hour_orig
                        date_day_next = date_day_next_orig
                        date_hour_next = date_hour_next_orig
                    endif

                    !write release info
                    !*************************************************
                    write(releaseunit, '(I8, I07.6)') date_day, date_hour
                    write(releaseunit, '(I8, I07.6)') date_day_next, date_hour_next
                    write(releaseunit, '(F9.4)') lon
                    write(releaseunit, '(F9.4)') lat
                    write(releaseunit, '(F9.4)') lon + real(dxdy_degr * release_dxdy_step)
                    write(releaseunit, '(F9.4)') lat + real(dxdy_degr * release_dxdy_step)
                    write(releaseunit, '(I9)') aboveSfc
                    write(releaseunit, '(F10.3)') bottomHeight_dustProf
                    write(releaseunit, '(F10.3)') topHeight_dustProf
                    write(releaseunit, '(I9)') nrPart

                    !write mass at this point for all species
                    !*************************************************
                    do ii = 0, numberSpecies - 1
                        if (ii .eq. current_species) then
                            write(releaseunit, '(E10.4)') combEmission
                        else
                            write(releaseunit, '(E10.4)')0.0
                        endif
                    end do

                    !Finish with name of this particular release
                    !*************************************************
                    write(releaseunit, '(A14,I8,I06.6,A1,I0.4,A8,I3)') 'Dust_location_', &
                    date_day, date_hour, '_', nr, '_species', startSpecies + current_species
                    write(releaseunit, *) '-----------------------------'
                    nr = nr + 1

                endif
            endif
            !********************************************************
        end do
    end do
    !************************************************************************

    !Finished writing for now, close the release file
    !************************************************************************
    close(releaseunit)
    !************************************************************************
    totalParticles = totalParticles + stepPart
end subroutine
