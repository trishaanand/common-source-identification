# The source files for Chapel
CHAPEL_SRC	= common-source-identification.chpl prnu.chpl \
		  read_jpg.chpl

# The source files for C
C_SRC		= ../c/read_jpg.c ../c/fastnoise.c ../c/grayscale.c \
		  ../c/prnu.c ../c/util.c ../c/wiener.c ../c/zeromean.c

# Each C file has a header file with the same name
C_HEADER	= $(C_SRC:.c=.h)

# The name of the program that we are making
PROGNAME	= common-source-identification

# The compiler and its flags.  We link to the jpeg, FFT, and math library
CC		= chpl
CFLAGS		= -ljpeg -lfftw3 -lm --fast

# The runtime flags.  With this flag we set environment variables that help the
# compiler to generate code for the right runtime.
RT_FLAGS	= CHPL_LAUNCHER=slurm-srun

# The following targets do not represent files and are therefore phony targets
.PHONY:		clean cleanslurm cleanoutput cleanall

# The program depends on the Chapel and C source
$(PROGNAME): $(CHAPEL_SRC) $(C_SRC) $(C_HEADER)
	$(RT_FLAGS) $(CC) -o $@ $(OBJ) $(CFLAGS) $(CHAPEL_SRC)

# Clean the generated binaries
clean: 
	rm -f $(PROGNAME) $(PROGNAME)_real

# Clean the Slurm output files
cleanslurm:
	rm -f *.out

# Clean the output matrix
cleanoutput:
	rm -f corrMatrix

# Clean everyting
cleanall: clean cleanslurm cleanoutput
