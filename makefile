fc       = ifort
fppflags = -cpp
fflags   = -g -O3 -mavx
openmp   = 
ldflags  =
libs     =
src      = mps.F90
bin      = mps

all: $(bin)
$(bin): $(src)
	$(fc) $(fppflags) $(fflags) $^ -o $@
TAGS: $(src)
	etags $^
tags: $(src)
	ctags $^
clean:
	rm -f $(bin) *~ *.mod core.* *.gif TAGS tags
