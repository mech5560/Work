CC  := g++
LINKER     := $(CC)

CFLAGS     := -std=c++0x -g3 -ggdb3 -Wall -O3


UTIL := ./Solver/Energy_Equation.o ./Solver/Energy_Equation_Corrector.o\
	 ./Solver/Energy_Viscous.o ./Solver/Energy_Convection.o\
	 ./Initialize/Mesh_Generation/Cubic_Mesh.o \
	./Initialize/Mesh_Generation/Hyperbolic_Mesh.o\
	 ./Initialize/Alloc_Mem/Matrix_Allocator.o ./Initialize/Alloc_Mem/Free_Matrix.o\
	 ./Initialize/Print_3D.o ./Solver/Next_Step.o\
	 ./Initialize/Initial_Conditions/Initial_Zero.o \
	./Initialize/Initial_Conditions/Pertubation_Introducer.o\
	./Initialize/Initial_Conditions/Initial_Reader.o\
	./Initialize/Initial_Conditions/Initial_Conditions_Turbulence.o\
	./Initialize/Initial_Conditions/Turbulence_Reader.o\
	./Solver/Intermediate_Velocity_X.o ./Solver/Residuals/Velocity_Residual_X.o\
	 ./Solver/Intermediate_Velocity_Y.o ./Solver/Residuals/Velocity_Residual_Y.o\
	 ./Solver/Intermediate_Velocity_Z.o ./Solver/Residuals/Velocity_Residual_Z.o\
	 ./Solver/Poisson/Vector_Constructor.o\
	 ./Solver/Poisson/Right_Hand_Side_Poisson.o ./Solver/Poisson/Divergence_X.o\
   ./Solver/Poisson/Divergence_Y.o ./Solver/Poisson/Divergence_Z.o ./Solver/Velocity_Update_X.o\
	 ./Solver/Velocity_Update_Y.o ./Solver/Velocity_Update_Z.o ./Solver/Flux_Evaluation_Y.o\
	 ./Solver/Flux_Evaluation_X.o ./Solver/Flux_Evaluation_Z.o\
	 ./Initialize/Initial_Conditions/Initial_One.o\
	 ./Initialize/Initial_Conditions/Initial_Christos.o\
	 ./Initialize/Initial_Conditions/Initial_Brown_2.o\
	 ./Initialize/Boundary_Conditions/BC_Velocities.o\
	 ./Initialize/Print_2D_Matrix.o ./Initialize/Print_1D_Matrix.o\
	 ./Initialize/Print_2D_Matrix_Ghost.o\
	 ./Initialize/Alloc_Mem/Allocator.o\
	 ./Initialize/Alloc_Mem/DeAllocator.o\
	 ./Initialize/Boundary_Conditions/BC_Tilda.o ./Initialize/Print_2D_Curve.o\
	 ./Initialize/Boundary_Conditions/BC_Single.o\
	 ./Header_Files/Density_Calculator.o\
	 ./Solver/Poisson/bcgc_solver_Printing.o\
	 ./Solver/Residuals/Convection_Term.o\
	 ./Solver/Residuals/Viscous_Component_XX.o\
	 ./Solver/Residuals/Viscous_Component_XY.o\
	 ./Solver/Residuals/Viscous_Component_XZ.o\
	 ./Solver/Residuals/Viscous_Component_YX.o\
	 ./Solver/Residuals/Viscous_Component_YY.o\
	 ./Solver/Residuals/Viscous_Component_YZ.o\
	 ./Solver/Residuals/Viscous_Component_ZX.o\
	 ./Solver/Residuals/Viscous_Component_ZY.o\
	 ./Solver/Residuals/Viscous_Component_ZZ.o\
	 ./Solver/Residuals/Forcing_Term_Christos_X.o\
	 ./Solver/Residuals/Forcing_Term_Christos_Y.o

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




clear-data:
	rm *.dat
	rm *.pdf


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
