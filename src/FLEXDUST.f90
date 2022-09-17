!**********************************************************************
! Copyright 2016,2017                                                 *
! Christine Groot Zwaaftink                                           *
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
program FLEXDUST
    
    use par_mod
    use com_mod
    use dust_mod
    
    implicit none
    
    !Loop variables
    !*************************************************************************************************
    integer 			:: ix, iy
    real 				:: lat_out, lon_out
    !*************************************************************************************************
    
    !Model variables
    !*************************************************************************************************
    real                            :: frictVelThres, shearStressThres
    real(kind=dp)                   :: juldate, start_date, start_date_available, end_date
    integer 			            :: tot_sec, tot_sec_end, time_last_wind,nstop, tmp_day, tmp_hour, time_wind_field, step_nc
    integer                         :: totalParticles, ix_ll, ix_ur, iy_ll, iy_ur, dummy_int
    integer                         :: calculationStep
    real            				:: snow, tmp_tot_emission, totalEmission
    character(len=200)              :: grid_filename
    !character(len=200)              :: tmp
    character(len=100)              :: release_file
    real,dimension(:,:), allocatable    :: soilFraction(:,:)
    real,dimension(:,:), allocatable    :: gridarea(:,:)
    real,dimension(:,:), allocatable    :: emission_mass(:,:),emission_flux(:,:), soilMoisture(:,:), cum_emission(:,:)
    real,dimension(:,:), allocatable    :: emission_mass_step(:,:),emission_flux_step(:,:)
    real,dimension(:,:), allocatable    :: outputField(:,:), erodibility(:,:), erodibility_total(:,:)
    real,dimension(:,:), allocatable    :: sandMap(:,:), clayMap(:,:)
    real,dimension(:,:,:), allocatable  :: precipitation
    integer,dimension(:,:), allocatable :: inNestNr,  ix_wind_n, iy_wind_n, landcovertype
    integer,dimension(:), allocatable   :: ix_lu, iy_lu, ix_wind, iy_wind, ix_clay,  iy_clay, ix_erClass,  iy_erClass
    logical,dimension(:,:), allocatable :: inLU_n, inClayGrid,  inClassErosion
    integer,dimension(:,:),allocatable  :: ix_lu_n, iy_lu_n
    INTEGER                             :: ALLOC_ERR
    integer(kind=1)                     :: landinventory_global(0:nx_landuse-1,0:ny_landuse-1)
    integer            	        	    :: landinventory_n(0:nx_landuse_n(1)-1,0:ny_landuse_n(1)-1)
    integer, dimension(2)              :: timetesting
    real,dimension(:), allocatable   :: lons, lats
    integer :: emission_time
    !*************************************************************************************************
    
    !Set some variables used in the code to read files copied from FLEXPART
    !*************************************************************************************************
    call readcommand
    write(*,*) 'Reading Namelist Complete' 
    
    numbnests           = numberOfNests   
    ldirect				= 1
    path(4) 			= ECMWF_input
    path(6)             = ECMWF_input_nest  
    length(4)			= index(path(4),' ')-1
    length(6)			= index(path(6),' ')-1
    start_date			= juldate(start_date_day, start_date_hour)
    start_date_available            = start_date	
    end_date			= start_date+releaseDays
    ideltas				= nint((end_date-start_date)*86400.)
    nc_file_out         = trim(output_directory)//trim(nc_file_name)
    summary_file    = trim(output_directory)//trim(summary_file_name)
    !*************************************************************************************************

    !Save settings dust_mod to output
    !*************************************************************************************************
    call system('cp '//trim('dust_mod.f90')//' '//trim(output_directory))
    !*************************************************************************************************


    !Check settings of time step and releaseDays
    !*************************************************************************************************
    if (mod(releaseDays*24,real(time_step)) .gt. (1e-4)) then
        print*, mod(releaseDays,real(time_step/24.))
        print*, 'Time step and number of release days (',releaseDays,') do not match.',&
        'Please change releaseDays to ',&
        (releaseDays-mod(releaseDays,real(time_step/24.))), 'or', &
        (releaseDays-mod(releaseDays,real(time_step/24.)))+time_step/24.
        stop
    endif

    calculationStep=1 !Step through simulation in hourly intervals. 
    !Emission is calculated per wind field and copied to next hours until next wind field is available.
    !The hourly calculationStep allows for time_steps of releases independent of the interval of wind fields.
    !Note that it does not improve the emission estimates as long as there is no additional input available.
    !*************************************************************************************************
    
    !Initialize the output grids
    !*************************************************************************************************
    allocate(lons(0:nx_lon_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(lats(0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(emission_mass(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(emission_mass_step(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"

    allocate(cum_emission(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(emission_flux(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(emission_flux_step(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"

    allocate(soilFraction(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(gridarea(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"

    allocate(erodibility(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(erodibility_total(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"

    allocate(sandMap(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(clayMap(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"   
    
    allocate(soilMoisture(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(outputField(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(precipitation(0:nx_lon_out-1,0:ny_lat_out-1,1:40), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(landcovertype(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(inLU_n(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(inClayGrid(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(inClassErosion(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(inNestNr(0:nx_lon_out-1,0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(ix_lu(0:nx_lon_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(iy_lu(0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(ix_erClass(0:nx_lon_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(iy_erClass(0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(ix_wind(0:nx_lon_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(iy_wind(0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(ix_clay(0:nx_lon_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(iy_clay(0:ny_lat_out-1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(ix_lu_n(0:nx_lon_out,0:numbnests_landuse), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(iy_lu_n(0:ny_lat_out,0:numbnests_landuse), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(ix_wind_n(0:nx_lon_out,0:numbnests), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    
    allocate(iy_wind_n(0:ny_lat_out,0:numbnests), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) STOP "*** Not enough memory ***"
    !*************************************************************************************************
    
    release_file=trim(output_directory)//trim(release)
    write(*,*) 'Start FLEXDUST simulation'
    !print*, 'Scaling factor:', scalingFactor
    
    !Read the global and nested landuse file and determine soil fraction per grid point
    !*************************************************************************************************
    write(*,*) 'Get land cover and soil data'
    
    call getLanduseGlobal(landuse_file, nx_landuse, ny_landuse,landuse_file_type,landuse_binary,landinventory_global)
    
    if( numbnests_landuse.ge.1)then
        !read nested land use
        call read_landuse_nest(landuse_file_n(1), ny_landuse_n(1), nx_landuse_n(1), landinventory_n, landuse_bin_fileType)
    endif
    
    !Retrieve erosion classes Iceland
    if( applyClassErosion)then
        call read_landuse_nest(erClass_file, ny_erC, nx_erC, erClass, erClass_bin_fileType)
    else
        erClass(:,:)=0
        ix_erClass(:)=0
        iy_erClass(:)=0
    endif
    
    call readClayContent()
    call readSandContent()
    !*************************************************************************************************
    
    
    !Prepare to read wind fields (use FLEXPART routines)
    !*************************************************************************************************
    write(*,*) 'Get wind field info'
    call readavailable(start_date_available, start_date+releaseDays+1)
    call gridcheck()
    if(numbnests.gt.0) call gridcheck_nests()
    
    !Prepare links between grids
    !*************************************************************************************************
    call getGridPoints(ix_lu,iy_lu,inLU_n,ix_lu_n,iy_lu_n,ix_wind,iy_wind,inClayGrid,ix_clay,iy_clay,&
    inNestNr,ix_wind_n,iy_wind_n,inClassErosion, ix_erClass,iy_erClass)
    !*************************************************************************************************

    !Some initializations
    !*************************************************************************************************
    tot_sec=0
    time_wind_field=0
    tot_sec_end= releaseDays*24*3600
    totalEmission=0
    totalParticles=0
    cum_emission(:,:)=0
    precipitation(:,:,:)=0
    step_nc=0
    time_last_wind=0
    !*************************************************************************************************

    !Loop through complete time series of releaseDays in hourly steps
    !*************************************************************************************************
    do while(tot_sec.lt.tot_sec_end)
        write(*,*) 'Currently at ' ,tot_sec, ' (sec) out of total:', tot_sec_end!, ', last wind:', time_last_wind
        call flush()
        
        !Check if the wind field is up to date
        !****************************************************************
        if (tot_sec.eq.0 .or. tot_sec.ge.time_last_wind+time_step_wind*3600) then
                !Get wind field
            !*************************************************************
            call getfields(tot_sec,nstop, time_last_wind)
            print*, 'Finished reading wind, save time:', time_last_wind
            !*************************************************************
        endif
              
        !Store precipitation previous time steps
        !*************************************************************
        precipitation(:,:,1:39)=precipitation(:,:,2:40)
        !*************************************************************
        
        !In first time step determine soil fraction and check if some variables not standard in FLEXPART are available
        !*************************************************************
        if (tot_sec .eq. 0 )then
            call getSoilFromLU(soilFraction, landinventory_global, landinventory_n, inLU_n, ix_lu_n, iy_lu_n, ix_lu, iy_lu, &
            ix_wind, iy_wind)
            
            !Check if some read fields actually have values
            print* 
            print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print*, 'CHECK if the following values are non-zero'
            print*, 'If 0 missing variables in ECMWF fields; reconsider if you should change switches or get additional data!'
            print*, 'Total vegetation cover high/low:', sum(sum(cvh(:,:, 1, 1), 2), 1), sum(sum(cvl(:,:, 1, 1), 2), 1)
            print*, 'Total soil moisture:', sum(sum(svw(:,:, 1, 1), 2), 1)
            print*, 'Total precipitation:', sum(sum(lsprec(:,:, 1, 1), 2), 1), sum(sum(convprec(:,:, 1, 1), 2), 1)
            print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print* 
            !call writeGrid(output_directory//'Soil_fraction.dat',soilFraction,nx_lon_out, ny_lat_out)
        endif
        !*************************************************************
        
        
        if (time_last_wind.eq.tot_sec .or. PRECIP_BLOCK)then! Need to check how to efficiently run with PRECIPBLOCK!
            
            !Loop through release grid
            !*************************************************************
            do iy=0, ny_lat_out-1
                do ix=0,nx_lon_out-1

                    !Reset emission that was based on previous wind field
                    !********************************************************
                    emission_mass_step(ix,iy)=0.
                    emission_flux_step(ix,iy)=0.

                    !Calculate coordinates and catch coordinates outside grid
                    !********************************************************
                    if (tot_sec .eq. 0) then
                        lats(iy)=lat_bottom+iy*dx_dy_out
                        lons(ix)=lon_left+ix*dx_dy_out
                        
                        if(lons(ix).gt. 180.) lons(ix)= lons(ix)-360.
                        if(lons(ix).lt. -180.) lons(ix)= lons(ix)+360.

                        call calcGridArea(gridarea(ix,iy), lats(iy),dx_dy_out)
                        
                    endif
                    lat_out=lats(iy)
                    lon_out=lons(ix)
                    
                    if(lat_out.lt.-90 .or. lat_out.gt.90 .or. lon_out.lt.-180 .or. lon_out.gt.180)then
                        write(*,*) 'Incorrect latitude or longitude', lon_out, lat_out
                        stop
                    endif
                    !********************************************************
                    
                    !Add current precipitation
                    !********************************************************
                    precipitation(ix,iy,40)=lsprec(ix_wind(ix), iy_wind(iy),1,1)+convprec( ix_wind(ix), iy_wind(iy),1,1) !Large scale + convective precipitation in m
                    !********************************************************
                    if(tot_sec.eq.0)then
                        clayMap(ix,iy)=clayContent(ix_clay(ix),iy_clay(iy))
                        sandMap(ix,iy)=sandContent(ix_clay(ix),iy_clay(iy))
                    endif
                    !In first time step only determine erodibility at grid point if switched on
                    !********************************************************
                    if(tot_sec.eq.0 .and. EROSION_TOPO)then
                        
                        erodibility_total(ix,iy)=0
                        
                        !get lower left corner (ix_ll, iy_ll) of erosion area
                        call getGridPointWind(lat_out-2., lon_out-2,dummy_int, ix_ll, iy_ll, dummy_int, dummy_int)
                        !get upper right corner (ix_ur, iy_ur) of erosion area
                        call getGridPointWind(lat_out+2., lon_out+2.,dummy_int, ix_ur, iy_ur,dummy_int,dummy_int)
                        !scale erodibility in this area
                        call getErodibility(erodibility(ix,iy), ix_wind(ix), iy_wind(iy), ix_ll, ix_ur, iy_ll, iy_ur)
                        
                        erodibility_total(ix,iy)=erodibility_total(ix,iy)+erodibility(ix,iy)
                        
                        !get lower left corner (ix_ll, iy_ll) of erosion area
                        call getGridPointWind(lat_out-5., lon_out-5,dummy_int, ix_ll, iy_ll, dummy_int, dummy_int)
                        !get upper right corner (ix_ur, iy_ur) of erosion area
                        call getGridPointWind(lat_out+5., lon_out+5.,dummy_int, ix_ur, iy_ur,dummy_int,dummy_int)
                        !scale erodibility in this area
                        call getErodibility(erodibility(ix,iy), ix_wind(ix), iy_wind(iy), ix_ll, ix_ur, iy_ll, iy_ur)
                        
                        erodibility_total(ix,iy)=erodibility_total(ix,iy)+erodibility(ix,iy)
                        
                        !get lower left corner (ix_ll, iy_ll) of erosion area
                        call getGridPointWind(lat_out-10., lon_out-10,dummy_int, ix_ll, iy_ll, dummy_int, dummy_int)
                        !get upper right corner (ix_ur, iy_ur) of erosion area
                        call getGridPointWind(lat_out+10., lon_out+10.,dummy_int, ix_ur, iy_ur,dummy_int,dummy_int)
                        !scale erodibility in this area
                        call getErodibility(erodibility(ix,iy), ix_wind(ix), iy_wind(iy), ix_ll, ix_ur, iy_ll, iy_ur)
                        
                        erodibility_total(ix,iy)=erodibility_total(ix,iy)+erodibility(ix,iy)
                        
                        !get lower left corner (ix_ll, iy_ll) of erosion area
                        call getGridPointWind(lat_out-15., lon_out-15,dummy_int, ix_ll, iy_ll, dummy_int, dummy_int)
                        !get upper right corner (ix_ur, iy_ur) of erosion area
                        call getGridPointWind(lat_out+15., lon_out+15.,dummy_int, ix_ur, iy_ur,dummy_int,dummy_int)
                        !scale erodibility in this area
                        call getErodibility(erodibility(ix,iy), ix_wind(ix), iy_wind(iy), ix_ll, ix_ur, iy_ll, iy_ur)
                        
                        erodibility_total(ix,iy)=erodibility_total(ix,iy)+erodibility(ix,iy)
                        
                        erodibility_total(ix,iy)=erodibility_total(ix,iy)/4
                        
                        soilFraction(ix,iy)=soilFraction(ix,iy)*erodibility_total(ix,iy)

                        !Store soil fraction grid when finished > now only save in netcdf
                        !if(iy.eq.ny_lat_out-1 .and. ix.eq.nx_lon_out-1)then
                        !    print*,'Write a soil fraction grid for inspection'
                        !   call writeGrid(output_directory//'Soil_fraction_topo.dat',soilFraction,nx_lon_out, ny_lat_out)
                        !endif
                    endif
                    !********************************************************
                    snow=0.
                    !Only do further checks/calculations for points with bare soil
                    !********************************************************
                    !if(soilFraction(ix,iy).gt.1e-8)then
                    if(soilFraction(ix,iy).gt.1./(dx_dy_out/dxdy_degr_landuse))then !At least 1 landuse tile must be bare soil

                        !Get snow depth
                        !********************************************************
                        if(inNestNr(ix,iy).eq.0)then
                            snow=sd(ix_wind(ix),iy_wind(iy),1,1)
                        else
                            snow=sdn(ix_wind_n(ix,inNestNr(ix,iy)),iy_wind_n(iy,inNestNr(ix,iy)),&
                            1,1,inNestNr(ix,iy))
                        endif
                        !********************************************************                      
                        
                        !Check landseamask: if there is soil fraction at a 'ECMWF sea-point' than use mean snow cover of surrounding land-gridpoints (if available)
                        !********************************************************
                        if(correctLSM_SNOW)then
                            if(lsm(ix_wind(ix),iy_wind(iy)).eq.0. .and. snow.eq.0.) then
                                if(inNestNr(ix,iy).eq.0)then
                                    snow=sum(sd(ix_wind(max(ix-1,0):min(ix+1,nx_lon_out-1)), &
                                    iy_wind(max(iy-1,0):min(iy+1,ny_lat_out-1)),1,1)*&
                                    lsm(ix_wind(max(ix-1,0):min(ix+1,nx_lon_out-1)),iy_wind( max(iy-1,0):min(iy+1,ny_lat_out-1))))/&
                                    (sum(lsm(ix_wind(max(ix-1,0):min(ix+1,nx_lon_out-1)),&
                                    iy_wind( max(iy-1,0):min(iy+1,ny_lat_out-1))))+1.e-9)
                                endif
                                
                            endif
                        endif
                        !********************************************************
                        
                        !Erosion is only possible if there is less than snowLimit in dust_mod. 
                        !(Test:If the snow depth is greater than 7 m assume it's a glacier in ECMWF,
                        !but it's not according to the high-resolution land use file so allow erosion anyway)
                        !********************************************************
                        !if(snow.lt.snowLimit .or. snow.gt.7.)then 
                        if(snow.lt.snowLimit)then     
                            
                            emission_time=calculationStep !Emission is representative of how long period?
                            
                            if(emissionModel.eq.1)then
                                !simple model from Harald Sodemann, see Sodemann et al., 2015
                                !*****************************
                                call dustEmmission_HSO(emission_mass_step(ix,iy), soilFraction(ix,iy), &
                                lat_out, dxdy_degr_landuse, inNestNr(ix,iy), ix_wind(ix), iy_wind(iy), &
                                ix_wind_n(ix,inNestNr(ix,iy)), iy_wind_n(iy,inNestNr(ix,iy)), emission_time, &
                                scalingFactor, mobThreshold)
                                !*****************************
                            endif
                            
                            if(emissionModel.ge.2)then
                                call calcThreshold_tau(frictVelThres, shearStressThres,inClayGrid(ix,iy), clayContent(ix_clay(ix),&
                                iy_clay(iy)),sandContent(ix_clay(ix),iy_clay(iy)),inNestNr(ix,iy), ix_wind(ix), &
                                iy_wind(iy), ix_wind_n(ix,inNestNr(ix,iy)), iy_wind_n(iy,inNestNr(ix,iy)), &
                                inClassErosion(ix,iy),ix_erClass(ix),iy_erClass(iy), soilMoisture(ix,iy), &
                                soilFraction(ix,iy), precipitation(ix,iy,:))
                                
                                if(emissionModel.eq.2)then                                         
                                    call dustEmmission_MB95_tau(emission_mass_step(ix,iy), soilFraction(ix,iy), &
                                    lat_out, dx_dy_out, inNestNr(ix,iy),  ix_wind(ix), iy_wind(iy), &
                                    ix_wind_n(ix,inNestNr(ix,iy)), iy_wind_n(iy,inNestNr(ix,iy)), emission_time, &
                                    scalingFactor, frictVelThres, shearStressThres, inClayGrid(ix,iy), &
                                    clayContent(ix_clay(ix),iy_clay(iy)),&
                                    emission_flux_step(ix,iy), ustarVersion, gridarea(ix,iy))
                                endif
                                
                                if(emissionModel.eq.3)then     
                                    !Use model from Kok et al., 2014
                                    !*****************************
                                    call dustEmmission_Kok14(emission_mass_step(ix,iy), soilFraction(ix,iy), &
                                    lat_out, dx_dy_out, inNestNr(ix,iy),  ix_wind(ix), iy_wind(iy), &
                                    ix_wind_n(ix,inNestNr(ix,iy)), iy_wind_n(iy,inNestNr(ix,iy)),emission_time, scalingFactor, &
                                    frictVelThres, inClayGrid(ix,iy), clayContent(ix_clay(ix),iy_clay(iy)), & 
                                    emission_flux_step(ix,iy), gridarea(ix,iy))
                                endif
                                !*****************************
                            endif!emission model with varying threshold
                            
                        endif !snow
                    endif !soil
                end do !x
            end do!y

        else
            !No new wind field in this step. 
            !As long as no interpolation is implemented, no need to calculate emission again for this step.
        endif

        !Add emission of this hourly step to emission to be written in release/grid once necessary
        !*********************************
        emission_mass=emission_mass+emission_mass_step
        emission_flux=emission_flux+emission_flux_step
        !*********************************
        
        !Save emission in cumulative field of complete simulation
        !*********************************
        cum_emission=cum_emission+emission_flux_step
        !*********************************

        !Check if release and emission needs to be written. 
        !************************************************************************
        if ( mod((tot_sec/3600)+calculationStep,time_step).eq.0. ) then
            
            !Write the dust emission for this time step as a grid in a binary file
            !************************************************************************
            if (writeGridEmission)then
                !print*, 'cgz debug: writing grid at time: ',tot_sec
                !Old binary format:
                !grid_filename = output_directory // 'DustEmissionFlux_' // trim(tmp) // '.bin'
                !call writeGridBin(grid_filename, emission_flux, nx_lon_out, ny_lat_out)
                call caldate(start_date + real(tot_sec)/(3600. * 24.), tmp_day, tmp_hour)
                
                !Switch to netCDF output
                if(step_nc.eq.0)then
                    !Initialize netcdf file
                    call netCDF_prepareEmission(nc_file_out, lons, lats)
                    !Save bare soil fraction in netcdf out
                    call netCDF_write_grid(nc_file_out,"erodibility",erodibility_total)
                    call netCDF_write_grid(nc_file_out, "soil", soilFraction)
                    call netCDF_write_grid(nc_file_out, "area", gridarea)
                    call netCDF_write_grid(nc_file_out, "sand", sandMap)
                    call netCDF_write_grid(nc_file_out, "clay", clayMap)  
                endif
                timetesting(1)=tot_sec-(time_step-1)*3600
                timetesting(2)=tot_sec+3600
                call netCDF_writeEmission(nc_file_out, emission_flux,timetesting, step_nc)           
                step_nc= step_nc+1 
            endif
            !************************************************************************
            
            !Write the dust emission for this time step in a RELEASE file for FLEXPART
            !************************************************************************
            if(RELEASEFILE)then
                call writeRELEASEfile(release_file, typeSizeDistr, particlesPerTonDust,&
                Junge_index, emission_mass, time_step, &
                start_date+real(tot_sec)/real(24*3600), &
                lon_left, nx_lon_out,lat_bottom, ny_lat_out, dx_dy_out, &
                release_dxdy_step,start_date_day, releaseDays, tmp_tot_emission, totalParticles,&
                minMassWrite)
                !Track total emission in this simulation
                totalEmission=totalEmission+tmp_tot_emission
                write(*,*) 'Mass of emitted mineral dust up to this time step:', totalEmission, totalParticles
                write(*,*) 'Mass still in memory:', sum(sum(emission_mass(:,:),2),1)
            else
                !Set emission to 0
                totalEmission=totalEmission+sum(sum(emission_mass(:,:),2),1)
                emission_mass(:,:)=0.
                write(*,*) 'Mass of emitted mineral dust up to this time step:', totalEmission
            endif
            !************************************************************************
            
            !Reset emission and soilMoisture for next step
            !*************************************************************
            emission_flux(:,:)=0. !kg/m2
            soilMoisture(:,:)=0.
            !*************************************************************
        else
            !The interval of output time steps is lower than that of wind field time steps. 
            !Emission is calculated for each new wind field but not written to output yet.
        endif
        
        !Advance time step
        !************************************************************************
        tot_sec=tot_sec+real(calculationStep)*3600.
        !************************************************************************
        
        call flush()
    end do
    
    !*************************************************************************************************
    DEALLOCATE(emission_mass, STAT = ALLOC_ERR)
    
    !Save cumulative dust emission flux field
    !grid_filename = output_directory // 'CumEmissionFlux.bin'
    !call writeGridBin(grid_filename, cum_emission, nx_lon_out, ny_lat_out)
    call netCDF_write_grid(nc_file_out, "cum_emission", cum_emission)
    
    write(*,*) 'Finsihed FLEXDUST simulation, total mass of emitted mineral dust:', totalEmission, ' kg , # of particles:', &
    totalParticles
    
    call writeSummary(totalEmission)
    
end program FLEXDUST
