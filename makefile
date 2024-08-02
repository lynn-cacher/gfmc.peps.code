PRES_DIR = .
WORK_DIR = ./console
DIRG = ${PRES_DIR}/general
DIRR = ${PRES_DIR}/randomc
DIRS = ${PRES_DIR}/special
EXE  = ${WORK_DIR}/exe
VPATH = ${DIRG}:${DIRR}:${DIRS}

CC       = mpiicpc
CPPFLAGS = -std=c++17 -qmkl -O3
CFLAGS   =
CXXFLAGS = $(CFLAGS)
COMPILE  = $(CC) $(CPPFLAGS) $(CXXFLAGS) -c

SRCS := $(wildcard ${DIRG}/*.cpp) $(wildcard ${DIRR}/*.cpp) $(wildcard ${DIRS}/*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))
DEPS := $(patsubst %.cpp,%.d,$(SRCS))

default: srun

srun: ${EXE}
	@echo '______________________________________________________________________________________________ sbatch _______'
	@echo ''
	@qsub ${WORK_DIR}/slurm.bscc.t6.txt
	@echo '_____________________________________________________________________________________________________________'
	@echo ''

${EXE}: compile $(DEPS) $(OBJS)
	@echo '_______________________________________________________________________________________________ link ________'
	@echo ''
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -o ${EXE} $(OBJS) -L /es01/paratera/sce4637/hylin/anaconda3/envs/pytorch/lib -lpython3.9 -lpthread
	

start:
	@echo '______________________________________________________________________________________________ start ________'
	@echo ''

explain:
	@echo '____________________________________________________________________________________________ explain ________'
	@echo ''
	@echo "the following information represents your prgram"
	@echo "final EXE name: $(EXE)"
	@echo "source files: $(SRCS)"
	@echo "object files: $(OBJS)"
	@echo "dependency files: $(DEPS)"

compile:
	@echo '____________________________________________________________________________________________ compile ________'
	@echo ''

%.d:%.cpp
	$(CC) -MM $(CPPFLAGS) -I /es01/paratera/sce4637/hylin/anaconda3/envs/pytorch/include/python3.9 $< > $@
	$(CC) -MM $(CPPFLAGS) -I /es01/paratera/sce4637/hylin/anaconda3/envs/pytorch/include/python3.9 $< | sed s/\\.o/\\.d/   >> $@

%.o:%.cpp
	$(COMPILE) -o  -I /es01/paratera/sce4637/hylin/anaconda3/envs/pytorch/include/python3.9 $@ $<

clean:
	-rm -rf $(OBJS) $(DEPS) $(EXE)

clear clear:
	-rm -rf bi/*
	-rm -rf io/*
	-rm -rf haha.*
	-rm -rf slurm-*

depend:$(DEPS)
	@echo "dependencies are now up-to-date."

-include $(DEPS)
