# *********************************************************
#                                                     
#                  PLUTO 4.2  Makefile  
#                                                     
# *********************************************************

 # Default target
all: pluto                             

ARCH         = Linux.cuda.gcc.mpi.defs
PLUTO_DIR    = /mnt/c/Users/gabri/OneDrive/Documents/Physics_Y4/pluto
SRC          = $(PLUTO_DIR)/Src
INCLUDE_DIRS = -I. -I$(SRC)
VPATH        = ./:$(SRC)/New:$(SRC):$(SRC)/Time_Stepping:$(SRC)/States

include $(PLUTO_DIR)/Config/$(ARCH)

# ---------------------------------------------------------
#         Set headers and object files 
# ---------------------------------------------------------

HEADERS = pluto.h prototypes.h structs.h definitions.h macros.h mod_defs.h \
		  plm_coeffs.h
CUDA_HEADERS = plutoCuda.cuh block_size.cuh macrosCuda.cuh structsCuda.cuh \
			   prototypesCuda.cuh inlineCuda.cuh

OBJ = adv_flux.o arrays.o boundary.o check_states.o  \
      cmd_line_opt.o entropy_switch.o  \
      flag_shock.o flatten.o get_nghost.o   \
      init.o int_bound_reset.o input_data.o \
      mappers3D.o mean_mol_weight.o \
      parse_file.o plm_coeffs.o rbox.o \
      set_indexes.o set_geometry.o set_output.o \
      tools.o var_names.o  

OBJ += bin_io.o colortable.o initialize.o jet_domain.o \
       main.o restart.o runtime_setup.o show_config.o  \
       set_image.o set_grid.o startup.o split_source.o \
       userdef_output.o write_data.o write_tab.o \
       write_img.o write_vtk.o 

CUDA_OBJ = toolsCuda.o memoryCuda.o mappers3DCuda.o boundaryCuda.o initCuda.o \
		   flagCuda.o flag_shockCuda.o adv_fluxCuda.o

include $(SRC)/Math_Tools/makefile

# ---------------------------------------------------------
#  Define macros by adding -D<name> where <name> has been
#  set to TRUE in the system configuration file (.defs) 
# ---------------------------------------------------------

ifeq ($(strip $(PARALLEL)), TRUE)
 CFLAGS += -I$(SRC)/Parallel -DPARALLEL
 NVFLAGS += -DPARALLEL_CUDA
 include $(SRC)/Parallel/makefile
 
 ifeq ($(strip $(USE_ASYNC_IO)), TRUE)
  CFLAGS += -DUSE_ASYNC_IO
 endif
endif

ifeq ($(strip $(USE_CUDA)), TRUE)
 CFLAGS += -DUSE_CUDA
 NVFLAGS += -DUSE_CUDA
endif

ifeq ($(strip $(USE_HDF5)), TRUE)
 CFLAGS += -DUSE_HDF5
 OBJ    += hdf5_io.o
endif
      
ifeq ($(strip $(USE_PNG)), TRUE)
 CFLAGS += -DUSE_PNG
endif

ifeq ($(strip $(DEBUG_CUDA)), TRUE)
	CFLAGS += -DDEBUG_CUDA
	NVFLAGS += -DDEBUG_CUDA
endif

-include local_make

# ---------------------------------------------------------
#   Additional_header_files_here   ! dont change this line
# ---------------------------------------------------------


# ---------------------------------------------------------
#   Additional_object_files_here   ! dont change this line
# ---------------------------------------------------------

OBJ += plm_states.o
OBJ += rk_step.o
OBJ += update_stage.o
OBJ += parabolic_flux.o
CUDA_OBJ += plm_statesCuda.o
CUDA_OBJ += rk_stepCuda.o
CUDA_OBJ += analysis_reduction.o
include $(SRC)/HD/makefile
include $(SRC)/Viscosity/makefile
include $(SRC)/EOS/Isothermal/makefile
include $(SRC)/Nbody/makefile

# ---------------------------------------------------------
#    PLUTO target rule
# ---------------------------------------------------------

ifeq ($(strip $(USE_CUDA)), TRUE)
pluto: $(OBJ) $(CUDA_OBJ)
	$(NVCC) $(GPU_ARCH) $(CUDA_LINK_FLAGS) -o $(CUDA_LINK_OBJ) $(CUDA_OBJ)
	$(CC) $(OBJ) $(CUDA_OBJ) $(CUDA_LINK_OBJ) $(LDFLAGS) -o $@
else
pluto: $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $@
endif

# ---------------------------------------------------------
#                    Pattern rules
# ---------------------------------------------------------

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE_DIRS) $<
	
ifeq ($(strip $(USE_CUDA)), TRUE)
%.o: %.cu
	$(NVCC) $(GPU_ARCH) $(NVFLAGS) $(INCLUDE_DIRS) $<
endif

.PHONY: clean
clean:
	@rm -f	*.o
	@echo make clean: done

# ---------------------------------------------------------
#          Dependencies for object files
# ---------------------------------------------------------

$(OBJ):  $(HEADERS) makefile
$(CUDA_OBJ): $(HEADERS) $(CUDA_HEADERS) makefile
