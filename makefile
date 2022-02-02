#~ all: shell
OBJ = real_field.o functions.o io_functions.o
#~ OBJ = src/*.o
COMPILER = mpiicpc -std=c++11
DEBUG = -g
PROFILING = -p
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
LIBS = -lhdf5 -lh5si
ATTEMPT = 0
# ATTEMPT ?= $(shell bash -c 'read -p "Target name: " pwd; echo $$pwd')

target: $(OBJ)
	$(COMPILER) $(DEBUG) $(OBJ) $(LIBS) -o ./rf_$(ATTEMPT)
	
real_field.o: real_field.cc real_field.h
	$(COMPILER) $(CFLAGS) real_field.cc 
	
functions.o: functions.cc real_field.h
	$(COMPILER) $(CFLAGS) functions.cc
 
io_functions.o: io_functions.cc real_field.h
	$(COMPILER) $(CFLAGS) io_functions.cc -lhdf5 -lh5si
	
run:
	@../binaries/rf_$(ATTEMPT) | tee std_out.d
	
clean:
	@rm *.o
	
tar:
	tar cvf shell_$(shell date +"%F_T%H-%M").tar real_field.cc real_field.h \
	functions.cc io.cc json.hpp parameters.json makefile velocity_field.d

	


