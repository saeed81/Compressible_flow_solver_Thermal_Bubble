#
include make.macro
#
LIB_SRC = src/par_kind.f90\
	  src/in_out_manager.f90\
          src/par_dcn.f90\
          src/step.f90\
	  src/phys_cst.f90\
	  src/iom.f90\
	  src/dcn_vctr.f90\
	  src/dcn_opt.f90\
	  src/io_ezcdf.f90\
	  src/dom_init.f90\
	  src/dcn_cycle.f90\
	  src/dcn.f90\
	  src/model.f90

LIB_OBJ = $(LIB_SRC:.f90=.o)
#
LIBDCN = lib/libdcn.a
#
LIB = -L./lib -ldcn
#
.SUFFIXES:
.SUFFIXES: .f90 .o
#
#
#main program
#
all: dcn.x
	@echo "                 "
	@echo DCN model is OK
	@echo "                 "

#
dcn.x:  $(LIBDCN)  
	$(FC) $(FF) -o dcn.x $(LIB) $(LIB_CDF) 
#
$(LIBDCN): $(LIB_OBJ)
	@echo ""
	@echo $(LIB_OBJ)
	@echo ""
	@mkdir -p lib
	ar -rv lib/libdcn.a  $(LIB_OBJ)
	ranlib lib/libdcn.a
	@echo ""
#
#
.f90.o: $(LIB_SRC)
	$(FC) -c $(FF) $< -o $*.o
#
#
clean:
	rm -f *.x *.out *.log *~ out src/*.o mod/* lib/* *.nc src/*.nc src/*.x obj/*.o *.PLT *.dat *.data
#
install:
	@cp -f *.x src/ && cd src && ./*.x
unistall:
	rm -rf bin