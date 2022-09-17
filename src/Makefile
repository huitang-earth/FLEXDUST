SHELL = /bin/bash
MAIN = FLEXDUST

FC       = ifort
#INCPATH  = /xnilu_wrk/flex_wrk/bin64/grib_api/include
#INCPATH1 = /usr/include
#INCPATH_NCDF =  /homevip/flexpart/gcc-5.4.0/include
INCPATH1 = ${EBROOTNETCDFMINFORTRAN}/include
INCPATH2 = /cluster/software/ecCodes/2.9.2-intel-2018b/include
#/usr/include

#LIBPATH1 = /xnilu_wrk/flex_wrk/bin64/grib_api/lib/
#LIBPATH1 = /homevip/flexpart/gcc-5.4.0/lib
#/usr/lib/x86_64-linux-gnu/ 
#/homevip/flexpart/gcc-5.4.0/lib
#LIBPATH2 =   /usr/lib/x86_64-linux-gnu/
LIBPATH1 = /cluster/software/netCDF-Fortran/4.4.4-intel-2018b/lib
LIBPATH2 = /cluster/software/ecCodes/2.9.2-intel-2018b/lib
LIBS = -leccodes_f90 -leccodes -lm -ljasper -lnetcdff -DUSE_NCF

FFLAGS   = -I$(INCPATH1) -I$(INCPATH2) -O2 -g -cpp -m64  -mcmodel=medium -lnetcdff -traceback
#FFLAGS   =   -I$(INCPATH1) -I$(INCPATH_NCDF) -g  -fbacktrace -mcmodel=medium -fbounds-check -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
LDFLAGS  = $(FFLAGS)  -L$(LIBPATH1) -L$(LIBPATH2) $(LIBS)  #-Wl,-rpath,$(LIBPATH1) $(LIBPATH2) $(LIBS) 
#LDFLAGS  = $(FFLAGS)  -L$(LIBPATH1) -lgrib_api_f90 -lgrib_api -lm -ljasper -lnetcdff
#-lnetcdff
#-L$(LIBPATH2)


MODOBJS = \
par_mod.o \
com_mod.o \
cmapf_mod.o \
conv_mod.o \
dust_mod.o

OBJECTS = \
FLEXDUST.o \
subroutinesSoil.o \
readInput.o \
juldate.o \
caldate.o \
shift_field.o \
shift_field_0.o \
pbl_profile.o \
psim.o \
psih.o \
ew.o \
getfields.o \
readwind_nests.o \
calcpar.o \
calcpar_nests.o \
verttransform.o \
verttransform_nests.o \
richardson.o \
scalev.o \
getvdep.o \
getvdep_nests.o \
obukhov.o \
calcpv.o \
calcpv_nests.o \
getrb.o \
getrc.o \
raerod.o \
partdep.o \
qvsat.o \
readavailable.o \
readwind.o \
gridcheck.o \
dustEmission_HSO.o \
getSizeDistribution.o \
writeRELEASEfile.o \
subroutines.o \
dustEmission_MB95.o \
dustEmission_MB95_tau.o\
dustEmission_Kok14.o \
gridcheck_nests.o \
calcThreshold.o \
calcThreshold_tau.o \
readcommand.o



$(MAIN): $(MODOBJS) $(OBJECTS)
	$(FC) *.o -o $(MAIN) $(LDFLAGS)

$(OBJECTS): $(MODOBJS)

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm *.o *.mod

cleanall:
	rm *.o *.mod $(MAIN)
