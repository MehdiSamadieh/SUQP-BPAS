CILKPP  = g++
CC = gcc
INCLUDES = ../../include/
COMPILEARG = -I$(INCLUDES) -c -g -fpic -fcilkplus -DLINUXINTEL64=1 -march=native
#-Wall
CCOMPILEARG = -I$(INCLUDES) -c  -g -fpic -DLINUXINTEL64=1 -march=native

vpath %.cpp ./ 
vpath %.c ./ 
vpath %.h $(INCLUDES)/Interval/ $(INCLUDES)/RationalNumberPolynomial/ $(INCLUDES)RingPolynomial/ $(INCLUDES)/DyadicRationalNumber $(INCLUDES)/

all: SMQP_Support-AA.o SMQP_Support_Test-AA.o SMQP_Support_Recursive-AA.o urpolynomial.o urpolynomial_taylorshift.o urpolynomial_realroot.o mrpolynomial.o mrpolynomial_mdd-altarr.o  urpolynomial-altarr.o SUQP_Support.o
#urpolynomial-altarr.o SUQP_Support.o

#all: SMQP_Support.o SMQP_Support_Test.o SMQP_Support_Recursive.o urpolynomial.o urpolynomial_taylorshift.o urpolynomial_realroot.o mrpolynomial.o mrpolynomial_mdd.o

altarray: SMQP_Support-AA.o SMQP_Support_Test-AA.o SMQP_Support_Recursive-AA.o mrpolynomial-altarr.o

%.o: %.cpp
	$(CILKPP) $(COMPILEARG) $<

%.o : %.c
	$(CC) $(CCOMPILEARG) $<

serial: COMPILEARG += -DSERIAL=1
serial: SMQP_Support.o SMQP_Support_Test.o SMQP_Support_Recursive.o mrpolynomial.o urpolynomial.o urpolynomial_taylorshift.o urpolynomial_realroot.o mrpolynomial_mdd.o

test:

clean:
	rm -rf *.o
