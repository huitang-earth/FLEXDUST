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
subroutine getLanduseGlobal(landuseFileName, nx, ny, file_type,landuseFile_bin, landinvent)
implicit none

character(*)    :: landuseFileName
integer         :: nx, ny, file_type
integer(kind=1) :: landinvent(0:nx-1,0:ny-1)
logical         :: landuseFile_bin
!***********************************************************************
 
!Get the landuse inventory 
!***********************************************************************
if(file_type.eq.1)then
    print*, ' WARNING: THIS LANDUSE FILE SHOULD NOT BE USED, OLD ROUTINES MAY NOT FUNCTION ANY MORE'
    call readlanduse_IGBP(landuseFileName, landinvent, nx, ny)
else
    call readlanduse(landuseFileName, landinvent, nx, ny, landuseFile_bin)
endif
!***********************************************************************
end subroutine getLanduseGlobal 

subroutine getErodibility(erodibility, ix_wind, iy_wind, ix_min, ix_max, iy_min, iy_max)
    use par_mod
    use com_mod
    implicit none
    
    integer             :: ix_wind, iy_wind, ix_min, ix_max, iy_min, iy_max
    real                :: erodibility, z_max, z_min, z_loc
    real,dimension(4)   :: values_max, values_min
    
    !Determine max and min value, but take care of boundary
    if(ix_max.lt.ix_min)then
        if(iy_max.lt.iy_min)then
            !split in 4 areas
            values_max(1)=maxval(oro(0:ix_max, 0:iy_max))
            values_max(2)=maxval(oro(0:ix_max, iy_min:ny-1))
            values_max(3)=maxval(oro(ix_min:nx-1, iy_min:ny-1))
            values_max(4)=maxval(oro(ix_min:nx-1, 0:iy_max)) 
            
            values_min(1)=minval(oro(0:ix_max, 0:iy_max))
            values_min(2)=minval(oro(0:ix_max, iy_min:ny-1))
            values_min(3)=minval(oro(ix_min:nx-1, iy_min:ny-1))
            values_min(4)=minval(oro(ix_min:nx-1, 0:iy_max)) 
        else
            !split in 2 areas
            values_max(1)=maxval(oro(0:ix_max, iy_min:iy_max))
            values_max(2)=maxval(oro(ix_min:nx-1, iy_min:iy_max))
            values_max(3)=values_max(2)
            values_max(4)=values_max(2)
            
            values_min(1)=minval(oro(0:ix_max, iy_min:iy_max))
            values_min(2)=minval(oro(ix_min:nx-1, iy_min:iy_max))
            values_min(3)=values_min(2)
            values_min(4)=values_min(2)
            
        endif
    else
        if(iy_max.lt.iy_min)then
            !split in 2 areas
            values_max(1)=maxval(oro(ix_min:ix_max, 0:iy_max))
            values_max(2)=maxval(oro(ix_min:ix_max,iy_min:ny-1))
            values_max(3)=values_max(2)
            values_max(4)=values_max(2)
            
            values_min(1)=minval(oro(ix_min:ix_max, 0:iy_max))
            values_min(2)=minval(oro(ix_min:ix_max,iy_min:ny-1))
            values_min(3)=values_min(2)
            values_min(4)=values_min(2)
        else
            !don't split
            values_max(1)=maxval(oro(ix_min:ix_max, iy_min:iy_max))
            values_max(2)=values_max(1)
            values_max(3)=values_max(1)
            values_max(4)=values_max(1)
            
            values_min(1)=minval(oro(ix_min:ix_max, iy_min:iy_max))
            values_min(2)=values_min(1)
            values_min(3)=values_min(1)
            values_min(4)=values_min(1)
        endif
    endif
    z_max=maxval(values_max(:))
    z_min=minval(values_min(:))
 
    z_loc=oro(ix_wind,iy_wind) !Orography in unit m (conversion in readwind)
    
    erodibility=real((z_max-z_loc)/(z_max-z_min))**5.0
    
    if(erodibility.gt.1.0)then
        print*, 'WARNING, erosion>1:', ix_wind, iy_wind, ix_min, ix_max, iy_min, iy_max, z_min, z_max, z_loc, erodibility
    endif
    
end subroutine getErodibility

subroutine getSoilFromLU(soilFraction, landinventory_global, landinventory_n, inLU_n, ix_lu_n, iy_lu_n, ix_lu, iy_lu, &
    ix_wind, iy_wind)
    use dust_mod
    use com_mod
    implicit none
    integer :: ix, iy
    real    :: soilFraction(0:nx_lon_out - 1, 0:ny_lat_out - 1)
    logical :: inLU_n(0:nx_lon_out - 1, 0:ny_lat_out - 1), useVEG2010
    integer :: ix_lu_n(0:nx_lon_out, 0:numbnests_landuse)
    integer :: iy_lu_n(0:ny_lat_out, 0:numbnests_landuse)
    integer :: ix_lu(0:nx_lon_out), ix_wind(0:nx_lon_out)
    integer :: iy_lu(0:ny_lat_out), iy_wind(0:ny_lat_out)
    real :: vegetationECMWF(0:nx-1, 0:ny-1)
    real*8 :: gridSum
    real :: vegFrac
    integer :: count, i, j
    integer(kind = 1) :: landinventory_global(0:nx_landuse - 1, 0:ny_landuse - 1)
    integer :: landinventory_n(0:nx_landuse_n(1) - 1, 0:ny_landuse_n(1) - 1)
    logical :: city

    print*, 'Building soil fraction grid'
    
    useVEG2010=.false.
    
    !Get vegetation cover of ECMWF in 2010 if no current information is available
    !***************************************************************
    if (sum(sum(cvh(:,:, 1, 1), 2), 1) .lt.0.01 .or. sum(sum(cvl(:,:, 1, 1), 2), 1).lt.0.01)then
        print*, 'Using vegetation of different year'
        open(unit = 654, file = '../INPUT/Vegetation_ECMWF_2010', action = 'read', form = 'unformatted', &
        access = 'stream')
        read(654) vegetationECMWF(:,:)
        close(654)
        useVEG2010=.true.
    endif

    !call flush()
    !open(unit = 654, file = '../INPUT/Vegetation_ECMWF_2010_NASA', form = 'unformatted')
    !write(654) cvh(:,:, 1, 1)+cvl(:,:, 1, 1)
    !close(654)

    !**************************************************************
    soilFraction(:,:) = 0.D0
    !Loop through output grid
    !**************************************************************
    do iy = 0, ny_lat_out - 2
        do ix = 0, nx_lon_out - 2

            !reset for new point
            !**************************************************
            gridSum = 0.D0;
            vegFrac = 0.;
            count = 0;
            city = .false.;

            !Make a grid with soil fraction from global landuse > only valid for GLCNMO land cover data!!
            !**************************************************
            do j = iy_lu(iy), min(iy_lu(iy + 1), ny_landuse - 1)
                do i = ix_lu(ix), min(ix_lu(ix + 1), nx_landuse - 1)
                    if (landinventory_global(i, j) .eq. 17)then
                        !bare land - sand
                        gridSum = gridSum + 1.0

                    else if (landinventory_global(i, j) .eq. 16)then
                        !bare land - gravel/rock
                        gridSum = gridSum + 0.4 !Partly erodible, topography to identify rock

                    else if ((landinventory_global(i, j) .eq. 10) .or. (landinventory_global(i, j) .eq. 8) .or. (landinventory_global(i, j) .eq. 9) .or. (landinventory_global(i, j) .eq. 11) .or. (landinventory_global(i, j) .eq. 13) ) then
                        !sparse vegetation, get vegetation fraction from ECWMF but due to lower resolution only allow a vegetation cover between 10 and 90%
                        if (useVEG2010)then
                            !Not available in wind field, get from fixed file
                            vegFrac=vegetationECMWF(ix_wind(ix), iy_wind(iy))
                        else
                            vegFrac = cvh(ix_wind(ix), iy_wind(iy), 1, 1) + cvl(ix_wind(ix), iy_wind(iy), 1, 1)
                        endif

                        !restrict vegetation fraction so high resolution land use 'sparse vegetation' cannot be overwritten by low-resolution ECMWF vegetation
                        vegFrac = max(vegFrac, 0.30)
                        vegFrac = min(vegFrac, 0.90)
                        gridSum = gridSum + (1 - vegFrac)
                    else
                        !no soil fraction to be added
                        gridSum = gridSum
                    endif
                    count = count + 1
                    if (landinventory_global(i, j) .eq. 18)then
                        city = .true.
                    endif
                end do
            end do
            

            !Checked all land use grid points for single FLEXDUST output grid point > get mean soilFraction
            !***********************************************************
            if( gridSum .gt. 1.)then
            soilFraction(ix, iy) = gridSum/count
            endif
            !if(soilFraction(ix,iy).gt. 1e-8)then
            !    print*, 'soil fraction:',soilFraction(ix, iy), gridSum, count
            !    print*, landinventory_global(ix_lu(ix):ix_lu(ix + 1),iy_lu(iy):iy_lu(iy + 1))
            !endif
            !Overwrite the area where a nested landuse grid is available > only applicable for sandy desert data Iceland
            !************************************************************************************************
            if (inLU_n(ix, iy))then
                count = 0
                gridSum = 0
             !   print*, ix, iy, iy_lu_n(iy, 1), ix_lu_n(ix, 1)

                do j = iy_lu_n(iy, 1), min(iy_lu_n(iy + 1, 1), ny_landuse_n(1) - 1)
                    do i = ix_lu_n(ix, 1), min(ix_lu_n(ix + 1, 1), nx_landuse_n(1) - 1)
                            if (landinventory_n(i, j) .eq. landuse_n_bare)then
                                !Iceland:Sand & pumice
                                gridSum = gridSum + 0.95
                             elseif (any(landuse_n_altBare .eq. landinventory_n(i, j))) then
                                !Iceland:Sandy lava and sandy gravel. Assume that not all transportable
                                gridSum = gridSum + 0.70
                            else
                                !no soil fraction to be added
                                gridSum = gridSum
                            endif
                        !endif
                        count = count + 1
                    end do
                end do
                soilFraction(ix, iy) = gridSum/max(count, 1)
                !print*, 'Soil fraction: ', soilFraction(ix, iy) 

                if(city)then
                    !In the global landuse file there is urban landcover, which is not represented in Nytjaland. Correct the soil factor to 0.
                    soilFraction(ix, iy) = 0.
                endif
                
            endif
            !********************************************         
        end do
    end do
end subroutine getSoilFromLU