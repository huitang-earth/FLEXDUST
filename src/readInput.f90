!**********************************************************************
! Copyright 2016,2017                                                 *
! Christine Groot Zwaaftink                                           *
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

subroutine readlanduse_IGBP(inputFile, landinvent_loc, nx_lu, ny_lu)

  !*****************************************************************************
  !                                                                            *
  !      Reads the landuse inventory into memory and relates it to Leaf Area   *
  !      Index and roughness length.                                           *
  !                                                                            *
  !      AUTHOR: Andreas Stohl, 10 January 1994                                *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! i                       loop indices                                       *
  ! landinvent(1200,600,13) area fractions of 13 landuse categories            *
  ! LENGTH(numpath)         length of the path names                           *
  ! PATH(numpath)           contains the path names                            *
  ! unitland                unit connected with landuse inventory              *
  !                                                                            *
  ! -----                                                                      *
  ! Sabine Eckhardt, Dec 06 - new landuse inventary                            *
  ! after                                                                      *
  ! Belward, A.S., Estes, J.E., and Kline, K.D., 1999,                         *
  ! The IGBP-DIS 1-Km Land-Cover Data Set DISCover:                            *
  ! A Project Overview: Photogrammetric Engineering and Remote Sensing,        *
  ! v. 65, no. 9, p. 1013-1020                                                 *
  !                                                                            *
  ! LANDUSE CATEGORIES:                                                        *
  !                                                                            *
  ! 1	Urban land                                                             *
  ! 2	Agricultural land					               *
  ! 3	Range land					                       *
  ! 4	Deciduous forest			                               *
  ! 5	Coniferous forest		                                       *
  ! 6	Mixed forest including wetland                                         *
  ! 7	water, both salt and fresh                                             *
  ! 8	barren land mostly desert                                              *
  ! 9	nonforested wetland	                                               *
  ! 10 	mixed agricultural and range land		                       *
  ! 11 	rocky open areas with low growing shrubs                               *
  ! 12  ice                                                                    *
  ! 13  rainforest                                                             *
  !                                                                            *
  !*****************************************************************************
  use par_mod
  use com_mod
  implicit none

  integer :: ix,jy,i,k,lu_cat,lu_perc, nx_lu, ny_lu
  integer(kind=1) :: ilr
  integer(kind=1) :: ilr_buffer(2160000)
  integer :: il,irecread
  real :: rlr, r2lr
  integer(kind=1) :: landinvent_loc(0:nx_lu-1,0:ny_lu-1,6)
  character(*):: inputFile
 ! integer, parameter :: unitland=11, unitsurfdata=12


  ! Read landuse inventory
  !***********************
  ! The landuse information is saved in a compressed format and written
  ! out by records of the length of 1 BYTE. Each grid cell consists of 3
  ! Bytes, which include 3 landuse categories (val 1-13 and 16 percentage
  ! categories) So one half byte is used to store the Landusecat the other
  ! for the percentageclass in 6.25 steps (100/6.25=16)
  ! e.g.
  ! 4 3  percentage 4 = 4*6.25 => 25% landuse class 3
  ! 2 1  percentage 2 = 2*6.25 => 13% landuse class 1
  ! 1 12 percentage 1 = 1*6.26 => 6.25% landuse class 12
  open(unitland,file=inputFile,status='old', &
       form='UNFORMATTED', convert='little_endian')
!  print*,unitland
  read (unitland) (ilr_buffer(i),i=1,3*nx_lu*ny_lu)
  close(unitland)
  irecread=1
  do ix=0,nx_lu-1
    do jy=0,ny_lu-1
  ! the 3 most abundant landuse categories in the inventory
  ! first half byte contains the landuse class
  ! second half byte contains the respective percentage
      do k=1,3
  ! 1 byte is read
        ilr=ilr_buffer(irecread)
  !      ilr=0
        irecread=irecread+1
  ! as only signed integer values exist an unsigned value is constructed
        if (ilr.lt.0) then
           il=ilr+256
        else
           il=ilr
        endif
  ! dividing by 16 has the effect to get rid of the right half of the byte
  ! so just the left half remains, this corresponds to a shift right of 4
  ! bits
        rlr=real(il)/16.
        lu_cat=int(rlr)
  ! the left half of the byte is substracted from the whole in order to
  ! get only the right half of the byte
        r2lr=rlr-int(rlr)
  ! shift left by 4
        lu_perc=r2lr*16.
        landinvent_loc(ix,jy,k)=lu_cat
        landinvent_loc(ix,jy,k+3)=lu_perc
      end do
    end do
  end do
  !call writeGridInteger('Landinventory.txt', landinvent_loc(:,:,1), nx_lu-1, ny_lu-1)
    ! Read relation landuse,z0
  !*****************************

  open(unitsurfdata,file='surfdata.t', &
       status='old')

  do i=1,4
    read(unitsurfdata,*)
  end do
  do i=1,numclass
    read(unitsurfdata,'(45x,f15.3)') z0(i)
  end do
  close(unitsurfdata)
    !*****************************
  write(*,*) '(readInput) Finished reading landuse file: ', inputFile
end subroutine readlanduse_IGBP

subroutine readClayContent()
    use dust_mod
    implicit none
    integer             :: ii,jj
    integer,parameter   :: unitclayfile=14
    !******************************************    
    
    
    if (ISRIC_soilmaps) then
      !Read the Soilgrids/ISRIC clay file
      !******************************************
      open(unit=unitclayfile, file=clayFile, form = 'unformatted')
      read(unitclayfile) clayContent
      close(unitclayfile)
      !******************************************   
    else
    
      !Read the vector into grid (Clay file from global soil data task)
      !******************************************
      open(unit=unitclayfile, file=clayFile, action='read')
      do ii=ny_c-1, 0 , -1 !change order iy to get grid in same direction as wind field and land use
              do jj=0, nx_c-1
               read(unitclayfile,*) clayContent(jj,ii)
               !Fill clay content for missing soil types with default value
               if (clayContent(jj,ii).eq.-1) then
                clayContent(jj,ii)=5.
               endif
               !print*, 'Clay', ii,jj, clayContent(jj,ii)
          end do
      end do
      !******************************************   
      close(unitclayfile)
    endif

    
end subroutine readClayContent

subroutine readSandContent()
    use dust_mod
    implicit none
    integer             :: i,j
    integer,parameter   :: unitsandfile=15
    
    if (ISRIC_soilmaps) then
      !Read the Soilgrids/ISRIC clay file
      !******************************************
      open(unit=unitsandfile, file=sandFile, form = 'unformatted')
      read(unitsandfile)sandContent
      close(unitsandfile)
      !******************************************
    else
      !Read the vector into grid (Sand file from global soil data task)
      !******************************************
      open(unit=unitsandfile, file=sandFile, action='read')
      do i=ny_c-1, 0 , -1 !change order iy to get grid in same direction as wind field and land use
              do j=0, nx_c-1
               read(unitsandfile,*) sandContent(j,i)
               !Fill sand content for missing soil types with default value
               if (sandContent(j,i).eq.-1) then
                sandContent(j,i)=95.
               endif          
              end do
      end do
      !******************************************   
      close(unitsandfile)

    endif

end subroutine readSandContent

subroutine read_landuse_nest(filename, ny_file, nx_file, out_grid_turned, binary)
    implicit none
    integer,parameter   :: unitfile=16
    integer             :: ny_file, nx_file, i, ix, iy
    character(*)        :: filename
    logical             :: binary
    integer, dimension(0:ny_file-1, 0:nx_file-1)      :: tmp
    integer, dimension(0:nx_file-1, 0:ny_file-1)      :: out_grid_turned
    
    !Switch for binary file written by this code previously
    !******************************************************
    if(.not.binary)then
        
        !Read Ascii grid
        !******************************************
        open(unit=unitfile, file=filename, status='old')
        !dump first 6 lines
        do i=1,6
            read(unitfile,*)
        end do
    
        !read grid into temporary grid   
        do iy=0,ny_file-1
            read(unitfile,*) tmp(iy,:)
            !print*, tmp(iy,:)
        end do
        close(unitfile)
        !******************************************
        !call writeGridInteger('../Nest_landinvent_tmp.dat',tmp,nx_file, ny_file)

        !convert grid to be consistent with wind fields and other landuse
        !******************************************
        do iy=0,ny_file-1
            do ix=0,nx_file-1
                out_grid_turned(ix,iy)=tmp(ny_file-1-iy,ix)
            end do    
        end do
        !******************************************

        !write a binary file for quicker reading
        !******************************************
        open(unit=987, file='ASCII_landusenest.bin',form='unformatted')
        write(987) out_grid_turned
        print*, ' Prepared ASCII_landusenest.bin from ', filename, ',save and use in next simulation'
        close(987)
        !******************************************
    else
        !Open binary file (as was prepared by lines above)
        !******************************************
        open(unit = unitfile, file = filename, form = 'unformatted')
	read(unitfile) out_grid_turned
        close(unitfile)
        !******************************************
    endif
end subroutine  read_landuse_nest

subroutine readlanduse(inputFile, landinvent_loc, nx_lu, ny_lu, binaryFile)
  implicit none
  integer :: ix,iy,nx_lu, ny_lu
  integer(kind=1) :: landinvent_loc(0:nx_lu-1,0:ny_lu-1)
  integer(kind=1) :: tmp(0:nx_lu-1,0:ny_lu-1)
  character(*):: inputFile
  integer, parameter :: unitfile=101
  logical :: binaryFile
  
  !Switch if the landuse file has been used before and saved as binary file
  !************************************************************************
  if (.not.binaryFile)then
        open(unit = unitfile, file = inputFile, status = 'old')

        !read grid into temporary grid 
        !******************************************
        do iy = 0, ny_lu - 1
            read(unitfile, *) tmp(:, ny_lu - 1 - iy)
        end do
        close(unitfile)
        !*****************************************

        !write a binary file for quicker reading
        !******************************************
        open(unit = 987, file = 'ASCII_landuse_convert.bin', form = 'unformatted')
        write(987) tmp
        print*, 'Save landuse file ASCII_landuse_convert.bin for quicker reading'
        close(987)
        !******************************************
    else
        !Open bindary file (as was prepared in previous simulation by lines above)
        !******************************************
        print*, inputfile
        !open(unit = unitfile, file = inputFile, form = 'unformatted')
        open(unit = unitfile, file = inputFile, form ='unformatted',status='unknown',action='read',access='stream')
        read(unitfile) tmp
        close(unitfile)
        !******************************************
    endif
    
    !Copy temporary grid to landinventory
    !******************************************
    landinvent_loc(0:nx_lu - 1, 0:ny_lu - 1) = tmp
    !******************************************
    
end subroutine readlanduse

