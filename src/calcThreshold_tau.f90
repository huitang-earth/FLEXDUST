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
subroutine calcThreshold_tau(frictVelThres, shearStressThres, inClay, clayCont_loc, sandCont_loc, &
    inNestNr, ix_wind, iy_wind, ix_wind_n, iy_wind_n, inErC, ix_erC, iy_erC, w, soilFract, precipMem)
    use com_mod
    use par_mod
    use dust_mod
    implicit none
    !*********************************************************************************
    integer :: ix_wind, iy_wind, ix_wind_n, iy_wind_n, inNestNr, ix_erC, iy_erC, ii, jj, count, gridSum
    real :: frictVelThres, shearStressThres, clayCont_loc, sandCont_loc, out_erClass
    real :: precipMem(1:40)
    real :: f_d, soilFract, tmp_1, tmp_2, z0s, z0m
    real :: w, w_acc, rho_air, Diam_p, precip
    real, parameter :: gamma_shao = 2.9e-4, An = 0.05  !An = 0.111
    real, parameter :: rho_p = 2500
    logical :: inClay, inErC
    !********************************************************************************* 

    !Get the current soil volumetric water content, air density and precipitation
    !******************************************************************************
    if (inNestNr .gt. 0)then
        w = svw(ix_wind, iy_wind, 1, 1) * 100 !Soil volumetric water content >not available in nested grids so far
        rho_air = rhon(ix_wind_n, iy_wind_n, 1, 1, inNestNr) !Air density
        precip = lsprecn(ix_wind_n, iy_wind_n, 1, 1, inNestNr) + convprecn(ix_wind_n, iy_wind_n, 1, 1, inNestNr) !Large scale + convective precipitation in m
        !soil_T = st(ix_wind, iy_wind, 1, 1) - 273.16
    else
        w = svw(ix_wind, iy_wind, 1, 1) * 100 !convert to percentage
        rho_air = rho(ix_wind, iy_wind, 1, 1)
        precip = lsprec(ix_wind, iy_wind, 1, 1) + convprec(ix_wind, iy_wind, 1, 1)
        !soil_T = st(ix_wind, iy_wind, 1, 1) - 273.16
    endif
    !save precip into precipMem
    precipMem(40) = precip

    !First check if there is information about the clay content available at this point
    !**********************************************************************************
    if (.not.inClay .or. clayCont_loc .le. 0)then

        !no information on clay amount, use clay sized particles (higher threshold)
        !******************************************************************************
        Diam_p = 10e-6
        shearStressThres = An**2 * ((rho_p - rho_air) * ga * Diam_p + &
        gamma_shao/Diam_p)
        frictVelThres = An * sqrt(((rho_p - rho_air)/rho_air) * ga * Diam_p + &
        gamma_shao/(rho_air * Diam_p))
        !*****************************************************************************

    else

        !first calculate threshold friction velocity depending on grain size
        !******************************************************************************
        if (sandCont_loc .gt. 0.)then
            !If sand is present, assume low threshold friction velocity, sandblasting can occur
            Diam_p = 76e-6
        else
            !Assume threshold friction velocity for clay sized particles
            Diam_p = 10e-6
        endif

        shearStressThres = An**2 * ((rho_p - rho_air) * ga * Diam_p + gamma_shao/Diam_p)
        frictVelThres = An * sqrt(((rho_p - rho_air)/rho_air) * ga * Diam_p + gamma_shao/(rho_air * Diam_p))
        !******************************************************************************
        
    endif

        !Change threshold friction velocity for certain regions in Iceland
        !******************************************************************************
        if (applyClassErosion)then
            if (inErC)then
                !Get mean erosion class
                count = 0
                gridSum = 0
                do ii = ix_erC, min(ix_erC + int(dx_dy_out/dxdy_erC), nx_erC - 1)
                    do jj = iy_erC, min(iy_erC + int(dx_dy_out/dxdy_erC), ny_erC - 1)
                        if (erClass(ii, jj) .ge. 3)then
                            gridSum = gridSum + erClass(ii, jj)
                            count = count + 1
                        endif
                    end do
                end do
                
                if (count .gt. 0) then
                    out_erClass = real(gridSum)/real(count)
                else
                    out_erClass = -999
                endif

                !New map erosion fields
                if (out_erClass .ge. 3 .and. out_erClass .lt. 3.7)then
                    !Erosion class has some value, adjust threshold friction velocity accordingly
                    frictVelThres = 0.70
                else if (out_erClass .ge. 3.7 .and. out_erClass .lt. 4.4)then
                    !Severe erosion
                    frictVelThres = 0.58
                else if (out_erClass .ge. 4.4)then
                    !extreme erosion
                    frictVelThres = 0.33
                else
                    frictVelThres = 0.87
                endif

            endif
            shearStressThres = frictVelThres * frictVelThres * rho_air
        endif
        !******************************************************************************

        !Block events in case of precipitation
        !******************************************************************************
        if (PRECIP_BLOCK) then
            !if(sum(precipMem(33:40))*3.0 .lt. 0.5 .or. sum(precipMem(1:39))*3.0 .lt. 0.1)then
            if (sum(precipMem(40:40)) .lt. 1.0 .or. sum(precipMem(1:39)) .lt. 0.1)then
                !No significant precipitation or very dry soil > allow erosion   
                frictVelThres = frictVelThres
                shearStressThres = shearStressThres
            else
                frictVelThres = 5.
                shearStressThres = 25.
            endif
        endif
        !******************************************************************************

        !If the water content exceeds the water content for capillary forces increase threshold
        !******************************************************************************
        if (SOILMOISTURE_DEP .and. inClay .and. clayCont_loc .gt. 0)then

            !check the critical soil volumetric water content
            !******************************************************************************
            !acc. to Kok et al 2012, following Fecan et al. 1999
            w_acc = 0.17 * clayCont_loc + 0.0014 * clayCont_loc * clayCont_loc !clayCont in percent!
            !******************************************************************************

            if (w .lt. w_acc)then
                !volumetric water content too low, no change of threshold shear stress
                frictVelThres = frictVelThres
                shearStressThres = shearStressThres
            else
                shearStressThres = shearStressThres * (1 + 1.21 * (w - w_acc)**0.68)
                frictVelThres = frictVelThres * sqrt(1 + 1.21 * (w - w_acc)**0.68)
            endif
        endif
        !******************************************************************************
    
    !If the region is vegetated, increase the threshold by f_d (Zender, 2003)
    !******************************************************************************
    if (OBSTACLES)then
        if (soilFract .lt. 0.2)then
            if ((applyClassErosion .and. .not.inErC) .or. .not.applyClassErosion)then
                z0m = 100e-6 !values z0m=100e-6 and z0s=33.3e-6 result in f_d=1.26
                z0s = 33.3e-6
                tmp_1 = log(z0m/z0s)
                tmp_2 = log(0.35 * ((0.1/z0s)**0.8))
                f_d = 1.0/(1.0 - tmp_1/tmp_2)
                !multiply with f_d squared because of conversion to stress
                shearStressThres = shearStressThres * f_d * f_d
                frictVelThres = frictVelThres * f_d
            endif
        endif
    endif
    !******************************************************************************


    !********************************************************************************** 
end subroutine
