CC  := g++
LINKER     := $(CC)	      

CFLAGS     := -std=c++0x -g3 -ggdb3 -Wall  -O3  -march=native  -fdevirtualize -fexpensive-optimizations \
          -fthread-jumps  -falign-functions  -falign-jumps -fdelete-null-pointer-checks -ftree-vrp\
          -falign-loops  -falign-labels -fcaller-saves -fcrossjumping -fgcse -ftree-slp-vectorize\
					-ftree-pre -finline-functions -funswitch-loops -fgcse-after-reload -fvect-cost-model\
          -fstrict-aliasing -fstrict-overflow -ftree-switch-conversion\
          -ftree-tail-merge\
	        -fvect-cost-model\
					-fgcse-lm -fsched-interblock  -fsched-spec -fschedule-insns  -fschedule-insns2\
					-fcse-follow-jumps  -fcse-skip-blocks -finline-small-functions\
          -frerun-cse-after-loop -ftree-partial-pre -fipa-cp-clone\
					-foptimize-sibling-calls -fregmove -ftree-vectorize\
          -findirect-inlining -fipa-sra -fpartial-inlining -fpeephole2\
          -freorder-blocks\
	        -freorder-functions


UTIL := ./Solver/Energy_Equation.o ./Solver/Energy_Equation_Corrector.o\
	 ./Solver/Energy_Viscous.o ./Solver/Energy_Convection.o\
	 ./Initialize/Mesh_Generation/Cubic_Mesh.o\
	 ./Initialize/Alloc_Mem/Matrix_Allocator.o ./Initialize/Alloc_Mem/Free_Matrix.o\
	 ./Initialize/Print_3D.o  ./Initialize/Initial_Conditions/Initial_Zero.o \
	 ./Solver/Intermediate_Velocity_X.o ./Solver/Residuals/Velocity_Residual_X.o\
	 ./Solver/Intermediate_Velocity_Y.o ./Solver/Residuals/Velocity_Residual_Y.o\
	 ./Solver/Intermediate_Velocity_Z.o ./Solver/Residuals/Velocity_Residual_Z.o\
	 ./Solver/Poisson/Vector_Constructor.o\
	 ./Solver/Poisson/Right_Hand_Side_Poisson.o ./Solver/Poisson/Divergence_X.o\
   ./Solver/Poisson/Divergence_Y.o ./Solver/Poisson/Divergence_Z.o ./Solver/Velocity_Update_X.o\
	 ./Solver/Velocity_Update_Y.o ./Solver/Velocity_Update_Z.o ./Solver/Flux_Evaluation_Y.o\
	 ./Solver/Flux_Evaluation_X.o ./Solver/Flux_Evaluation_Z.o\
	 ./Initialize/Initial_Conditions/Initial_One.o\
	 ./Initialize/Print_3D_File.o ./Initialize/Alloc_Mem/Allocator.o ./Initialize/Alloc_Mem/DeAllocator.o\
	 ./Initialize/Boundary_Conditions/BC_Flux.o ./Initialize/Print_2D_Data.o\
	 ./Initialize/Boundary_Conditions/BC_Single.o ./Initialize/Print_3D_Plane.o\
	 ./Initialize/Print_3D_Binary.o ./Header_Files/Density_Calculator.o ./Initialize/Initial_Conditions/Initial_Cos.o\
	 ./Solver/Poisson/bcgc_solver.o ./Initialize/Print_3D_Single.o\
	 ./Solver/Residuals/Convection_Term.o ./Solver/Residuals/Viscous_Component_XX.o ./Solver/Residuals/Viscous_Component_XY.o\
	 ./Solver/Residuals/Viscous_Component_XZ.o ./Solver/Residuals/Viscous_Component_YX.o ./Solver/Residuals/Viscous_Component_YY.o\
	 ./Solver/Residuals/Viscous_Component_YZ.o ./Solver/Residuals/Viscous_Component_ZX.o ./Solver/Residuals/Viscous_Component_ZY.o\
	 ./Solver/Residuals/Viscous_Component_ZZ.o

TEST_OBJS := main.o


%.o:  %.cpp
	$(CC) $(CFLAGS)  -c $< -o $@


all:	
	make emain;


emain :   $(TEST_OBJS) $(UTIL) 
	  $(LINKER) $(CFLAGS)  $(TEST_OBJS) $(UTIL) \
	  -o $@ 




run:
	make all	
	./emain 
	mv -f Z* ./Data/Z/
	mv -f Y* ./Data/Y/
	mv -f X* ./Data/X/	




check-syntax:
	g++ -Wall -o  ${UTIL}


clean:
	rm -f *.o *~core *.x

cleanall:
	rm -f *.o *~core *.x  *.dat 
	rm -f ./Solver/*.o ./Solver/*~ 
	rm -f ./Initialize/Mesh_Generation/*.o ./Initialize/Mesh_Generation/*~
	rm -f ./Solver/Residuals/*.o ./Solver/Residuals/*~
	rm -f ./Solver/Poisson/*.o ./Solver/Poisson/*~	
	rm -f ./Initialize/*.o ./Initialize/*~ 
	rm -f ./Data/X/* ./Data/Y/* ./Data/Z/*
	clear 