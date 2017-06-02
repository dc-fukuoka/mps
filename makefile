fc       = ifort
fppflags = -cpp -D_DEBUG
fflags   = -g -O3 -mavx
openmp   = 
ldflags  =
libs     =
bin      = mps

all: $(bin)
$(bin): mps.F90
	$(fc) $(fppflags) $(fflags) $^ -o $@
clean:
	rm -f $(bin) *~ *.mod core.* *.gif
