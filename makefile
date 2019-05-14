FORT     = ifort
NICE     = 
OPTS    =  -O3 -g -fPIC
#OPTS = -O0 -g -CB  -traceback -debug extended -debug-parameters all  -check -fpe0 -no-wrap-margin
LIB = -mkl
 

EXE     = cavity
OBJECTS =  prms.o mesh.o io.o BC.o timeadv.o trace.o TDMA.o SOR.o ADI.o postpro.o main.o
%.o: %.f90
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG) $<

$(EXE): $(OBJECTS)
	$(FORT) $(OPTS) $(NICE) -o $(EXE) $(OBJECTS) $(LIB)
prms.o: prms.f90 makefile
mesh.o: mesh.f90 prms.o makefile
trace.o: trace.f90 mesh.o prms.o makefile
TDMA.o: TDMA.f90 makefile
BC.o: BC.f90 mesh.o prms.o makefile
SOR.o: SOR.f90 mesh.o prms.o BC.o makefile
ADI.o: ADI.f90 mesh.o prms.o BC.o TDMA.o makefile
io.o: io.f90 mesh.o prms.o makefile
timeadv.o: timeadv.f90 mesh.o prms.o TDMA.o SOR.o ADI.o io.o makefile
postpro.o: postpro.f90 timeadv.o makefile
main.o: main.f90 prms.o mesh.o trace.o timeadv.o BC.o postpro.o makefile

clean:
	\rm *.o
	\rm $(EXE)
