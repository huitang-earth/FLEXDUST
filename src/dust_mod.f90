!**********************************************************************
! Copyright 2016,2017                                                 *
! Christine Groot Zwaaftink                                           *
! cgz@nilu.no
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

!*******************************************************************************
!This module contains parameters of the flexdust model

module dust_mod
    
    !Input files/settings
    !***********************************************************************
    !output time frame
    integer :: start_date_day, start_date_hour, time_step	  
    real    :: releaseDays	      
    !***********************************************************************
    
    !output grid
    !***********************************************************************

    character(len=256)   :: output_directory
    real           :: lat_bottom
    real           :: lon_left 
    real           :: dx_dy_out   !resolution of emission calculation in degree, should be larger than resolution of global landuse file (15/3600)
    integer        :: release_dxdy_step = 1    !Interval of x and y in which release file should be written 
                                                            !(2 means that calculated emission of 4 grid cells with resolution dx_dy_out will be combined in 1 FLEXPART release)
    integer         :: ny_lat_out 
    integer         :: nx_lon_out
    !***********************************************************************
    
    !Output files
    !***********************************************************************
    character(len=256)     :: release
    character(len=64)     :: summary_file_name
    character(len=64)     :: nc_file_name
    character(len=256)    :: summary_file
    character(len=256)    :: nc_file_out
    !***********************************************************************
    !Input files/settings (HARD coded) recompile is require when parameters
    !are adjusted
    ! - Changed so that parameters such start date and duration of simulation
    !   can be adjusted using a namelist COMAND file instead.
    !***********************************************************************
    
    !Windfields and properties
    !***********************************************************************
  
    character(*), parameter    :: ECMWF_input='/cluster/projects/nn2806k/huit/flexpart/AVAILABLE_WINDFIELDS_EA'
    character(*), parameter    :: ECMWF_input_nest= '../AVAILABLE_WINDFIELDS' !FLEXPART AVAILABLE file for nested wind field
    integer, parameter         :: numberOfNests = 0 !Number of nested wind fields, if more than 1 remember to change path and length! NOT tested....
    integer, parameter         :: time_step_wind = 3 !time step wind fields in hours, default 3
    !***********************************************************************
   
    !properties of the clayContent & sand file
    !************************ src/dust_mod.f90***********************************************

    !for default version (old version based on soil types)
    !logical, parameter         :: ISRIC_soilmaps=.true.
    !character(*), parameter    :: clayFile= '../INPUT/Clay.srf'        !File with clay content
    !character(*), parameter    :: sandFile= '../INPUT/Sand.srf'        !File with sand content
    !real, parameter            :: dx_c= 0.0833, dy_c=0.0833         !Resolution sand/clay grids
    !real, parameter            :: xlon0_c= -180.00, ylat0_c=-56.50  !Lower left corner of clay/sand grids
    !integer,parameter          :: nx_c=4320 , ny_c=1686             !Size sand/clay grids
    !real, dimension(0:nx_c-1,0:ny_c-1) :: clayContent, sandContent  !Sand and clay should be equal grids!!!
    
    !for alternative ISRIC version (new version, 3 min resolution), not validated
    logical, parameter         :: ISRIC_soilmaps=.true.
    character(*), parameter    :: clayFile= '../INPUT/ISRIC_clay.bin'        !File with clay content
    character(*), parameter    :: sandFile= '../INPUT/ISRIC_sand.bin'        !File with sand content
    real, parameter            :: dx_c= 0.05, dy_c=0.05         !Resolution sand/clay grids
    real, parameter            :: xlon0_c= -180.0, ylat0_c=-55.98  !Lower left corner of clay/sand grids
    integer,parameter          :: nx_c=7200 , ny_c=2774             !Size sand/clay grids
    real, dimension(0:nx_c-1,0:ny_c-1) :: clayContent, sandContent  !Sand and clay should be equal grids!!!

    !***********************************************************************
    
    !global landuse file and properties
    !***********************************************************************
    integer, parameter         :: landuse_file_type=2               !1>same as flexpart 2>MODIS
    logical, parameter         :: landuse_binary=.true.             !Is the landuse file already converted from ASCII to binary file?
    character(*), parameter    :: landuse_file= '../INPUT/landcover_GLCNMO_l.bin'
    integer, parameter         :: nx_landuse=86400, ny_landuse=43200 !Size landuse file
    real, parameter            :: dxdy_degr_landuse=15./3600.        !Resolution landuse file
    !***********************************************************************
    
    !Nested landuse file and properties
    !***********************************************************************
    integer, parameter         :: numbnests_landuse=1   !Developed and tested for only 1 nested field (sandy deserts Iceland or Antarctic), 
                                                        !requires further changes in the source code if other fields are used! 
                                                        !(Adjust code for bare land and possibly soil fraction calculation.)
    
	!Iceland nest
	character(*), parameter    :: landuse_file_n(1)= '/home/christine/Dusty/INPUT/Iceland_bareLand.bin'
    logical, parameter         :: landuse_bin_fileType=.true. !true for ASCII already converted to binary with routine in readInput.f90
    integer, parameter         :: nx_landuse_n(1)=39984, ny_landuse_n(1)= 11177
    real, parameter            :: xlon0_n(1)=-24.708731621941, ylat0_n(1)=63.300009852139
    real, parameter            :: dxdy_landuse_n(1)=0.00028783278362664!degr
    integer, parameter         :: landuse_n_bare=3 !With what value is bare land indicated in the nested land use?
    integer, parameter         :: landuse_n_altBare(1:2)=(/4,5/) !Alternative value that has limited soil availability, check code for details integer, dimension(0:nx_landuse_n(1)-1,0:ny_landuse_n(1)-1):: landuse_n

	
    !Additional information on erosion classes for Iceland (Arnalds 2001)
    !***********************************************************************
    logical, parameter         :: applyClassErosion=.false.       !Change mobilization threshold depending on erosion class?
    character(*), parameter    :: erClass_file= '/home/christine/Dusty/INPUT/ErosionClass.bin'
    logical, parameter         :: erClass_bin_fileType=.true. !true for ASCII already converted to binary with routine in readInput.f90
    integer, parameter         :: nx_erC=39927, ny_erC=11177
    real, parameter            :: xlon0_erC=-24.692251367597, ylat0_erC=63.299997666465
    real, parameter            :: dxdy_erC= 0.00028783387369149!degr
    integer, dimension(0:nx_erC-1,0:ny_erC-1):: erClass
    !***********************************************************************
    
    !Switches output
    !***********************************************************************
    logical, parameter          :: RELEASEFILE=.true.       !Write a FLEXPART release file
    logical, parameter          :: writeGridEmission=.true. !For each output time step, write a grid with emission flux (kg m-2), 
                                                            !practical for splitting in regions and doing FLEXPART simulations with changing number of particles
    !***********************************************************************
    
    !Model parameters
    !***********************************************************************
    real, parameter             :: mobThreshold = 0.3        !Default mobilization threshold should be wind speed or friction velocity, depending on choice "emissionModel", default should be 0.3 for emissionModel 2
    real, parameter             :: particlesPerTonDust = 0.8 !Number of particles to be released per ton of dust, adjust with resolution
    integer, parameter          :: typeSizeDistr=3           !Use size distribution as in DustBowl-Sodemann et al.2015 (1), or similar to Kok 2011 (2 & 3) with many small particles in 3
    integer, parameter          :: Junge_index = 0           !only for typeSizeDistr 1
    real*8, parameter           :: scalingFactor = 4.8e-4    !Default value 4.8e-4 for emissionModel 2
    integer, parameter          :: emissionModel = 2         !Choose from several emission Models (1: HSO, 2:MB95, 3:Kok et al. 2014), default and tested: 2. Options 1 and 3 may currently not be up to date with other changes in the model
    real, parameter             :: snowLimit =0.02           !From which snow amount should mobilization not be possible?
    real, parameter             :: minMassWrite=11.0         !Minimum emission (kg) for which to write a release > change depending on wanted resolution
    !***********************************************************************
    
    !Switches model
    !***********************************************************************
    logical, parameter          :: OBSTACLES=.true.          !Influence on roughness length by obstables?
    logical, parameter          :: EROSION_TOPO=.true.       !Add erodibility depending on topography acc. to Ginoux et al., 2001?
    logical, parameter          :: PRECIP_BLOCK=.false.      !Block dust emission in case of precipitation?
    logical, parameter          :: SOILMOISTURE_DEP=.true.   !Should threshold friction velocity increase with soil moisture?
    logical, parameter          :: correctLSM_SNOW=.true.    !Use snow cover of nearby land-points if the point of interest is on a 'sea-point' in ECMWF but has soil fraction>0 
    logical, parameter          :: ustarVersion=.true.	     !Calculate emission based on friction velocity, or shear stress > implemented for testing but does not change results with current model setup!
    !***********************************************************************
    
end module dust_mod
