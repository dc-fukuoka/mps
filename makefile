fc       = ifort
fppflags = -cpp -D_DEBUG
fflags   = -g -O3 -march=core-avx2
ldflags  =
libs     =
openmp   = -fopenmp
src      = mps.F90
bin      = mps

all: $(bin)
$(bin): $(src)
	$(fc) $(fppflags) $(fflags) $(openmp) $^ -o $@
TAGS: $(src)
	etags $^
tags: $(src)
	ctags $^
clean:
	rm -f $(bin) *~ *.mod core.* *.gif TAGS tags
