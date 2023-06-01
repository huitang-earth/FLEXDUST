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

subroutine dustEmmission_Kok14(em_mass, soilF, lat,dxdy_degr, inNestNr, ix_wind, iy_wind, ix_wind_n, iy_wind_n, &
                                time_int, scalingFactor, mobilisationThreshold,inClay, f_clay_tmp, em_flux, &
                                gridarea)

    use par_mod
    use com_mod

    implicit none

    real     :: mfac, em_mass, em_flux, lat, mass, mobilisationThreshold, dxdy_degr
    real*8   :: scalingFactor
    real     :: soilF
    real     :: u_star_local, rho_air, flux
    integer  :: ix_wind, iy_wind, ix_wind_n, iy_wind_n, time_int, inNestNr
    real, parameter :: C_alpha=1.8!2.7 !(Kok et al. 2014, part 2)
    real, parameter :: C_d0=4.4e-5
    real, parameter :: C_e=2.0
    real, parameter :: rho_air_0=1.225
    real, parameter :: u_star_st_0=0.16
    real :: f_clay, f_clay_tmp, C_d, tmp, f_bare, u_star_st, u_star_t, gridarea
    logical ::inClay

        !Get friction velocity
        if(inNestNr.eq.0)then
            u_star_local=ustar(ix_wind,iy_wind,1,1)
            rho_air=rho(ix_wind,iy_wind,1,1)
        else
            u_star_local=ustarn(ix_wind_n,iy_wind_n,1,1,inNestNr)
            rho_air=rhon(ix_wind_n,iy_wind_n,1,1,inNestNr)
        endif
        
        !Get threshold
        u_star_t=mobilisationThreshold
        u_star_st=u_star_t*SQRT(real(rho_air/rho_air_0))
        C_d=C_d0*exp(-C_e*(u_star_st-u_star_st_0)/u_star_st_0)
        f_bare=soilF
        
        !Get clay content
        if(inClay)then
            f_clay=real(f_clay_tmp/100)
            if(f_clay.lt.0)then
                f_clay=0.0 !No information on clay content
            endif
        else
            f_clay=0.1 !Should never occur since clay is available at each grid point with soil
        endif
        
       
        !******************************************************************
        !Calculate emitted mass (em_mass) per time step if the threshold for saltation is exceeded
        !******************************************************************	
	if(u_star_local.ge.u_star_t)then	
		
		!multiplication factor accounting for time step, area and conversion to kg
		!**********************************************************
		mfac=3600.0*time_int*gridarea	
				
		!Emitted mass in kg when mfac is used (otherwise mass flux in kg m-2 s-1)
                tmp=C_alpha*(u_star_st-u_star_st_0)/u_star_st_0
                flux=scalingFactor*C_d*f_bare*f_clay*(rho_air*(u_star_local*u_star_local-u_star_t*u_star_t))/u_star_st &
                        *((u_star_local/u_star_t)**tmp)
                mass=flux*mfac
                       
               em_mass=em_mass+mass
               em_flux=em_flux+flux*3600.*time_int !kg/m2
               
	endif
        
end subroutine dustEmmission_Kok14
