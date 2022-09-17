!**********************************************************************
! Copyright 1998,1999,2000,2001,2002,2005,2007,2008,2009,2010         *
! Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa,             *
! Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
!                                                                     *
! This file is part of FLEXPART.                                      *
!                                                                     *
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

subroutine readcommand

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the current model run.  *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     18 May 1996                                                            *
  !     HSO, 1 July 2014                                                       *
  !     Added optional namelist input                                          *
  !     Ove Haugvaldstad 2020 adopted for FLEXDUST                             *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  !                                                                            *
  ! Variables:                                                                 *
  ! unitcommand          unit connected to file COMMAND                        *
  ! start_date_day       start date of simulation                              *
  ! start_date_hour      start time of simulation                              *
  ! time_step            output time step                                      *
  ! releaseDays          Number of days to simulate emission flux              *
  ! output_directory     path to directory where output files shoul be stored  *
  ! dx_dy_out            output gridbox size                                   *
  ! release_dxdy_step    Interval of x and y in FLEXPART release file          *
  ! nx_lat_out           number of gridboxes in x direction                    *
  ! nx_lon_out           number of gridboxes in y direction                    *
  !                                                                            *
  !                                                                            *
  !                                                                            *
  !                                                                            *
  !                                                                            *
  !                                                                            *
  !                                                                            *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use dust_mod   
  implicit none

  real(kind=dp) :: juldate
  character(len=50) :: line
  logical :: old
  integer :: readerror
  
  old = .false.

  namelist /command_flexdust/ &
  start_date_day, &
  start_date_hour, &
  time_step, &
  releaseDays, &
  output_directory, &
  lat_bottom, &
  lon_left, &
  dx_dy_out, &
  release_dxdy_step, &
  release ,&
  ny_lat_out, &
  nx_lon_out, &
  nc_file_name, &
  summary_file_name

    
  ! Presetting namelist command
  start_date_hour=000000
  start_date_day=20190301
  time_step	  = 3
  releaseDays	  = 90
  output_directory  ='/cluster/work/users/ovewh/'
  lat_bottom        = 30    
  lon_left          = 60    
  dx_dy_out         = 0.1   
  release_dxdy_step = 1   
  ny_lat_out        = 320   
  nx_lon_out        = 680
  release= 'RELEASES_FLEXDUST'
  summary_file_name='Summary.txt'
  nc_file_name='FLEXDUST_out.nc'   
  ! Reading namelist
  open(unitcommand, file='../COMMAND_FLEXDUST', status='old', form='formatted', err=999)
  read(unitcommand,command_flexdust,iostat=readerror)
  close(unitcommand)

  if (time_step.eq.0) then
    write(*,*) ' #### FLEXDUST MODEL ERROR! TIME STEP MUST    #### '
    write(*,*) ' #### NOT BE ZERO                             #### '
    write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
    stop
  endif

  return

999   write(*,*) ' #### FLEXDUST MODEL ERROR! FILE "COMMAND"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  stop

1000   write(*,*) ' #### FLEXDUST MODEL ERROR! FILE "COMMAND"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  stop
end subroutine readcommand
