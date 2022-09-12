# Use: 
CILKPP  = g++
LIBARG  = -lmaplec -lbpas -lgmpxx -lgmp -lcilkrts -lmodpnLINUXINTEL64 -lmps -lpthread
COMPILEARG = -c -O2 -g -fcilkplus -DLINUXINTEL64=1 -march=native
TARGET  = test.bin


OBJECTFILES = ../../src/RationalNumberPolynomial/SUQP_Support.o ../MapleTestTool/MapleTestTool.o


#vpath %.o ./ ## ../../ /usr/local/lib

all: clean $(TARGET) test

$(TARGET): test.o ../../src/RationalNumberPolynomial/SUQP_Support.c
	$(CILKPP) -o $@ $< $(OBJECTFILES) $(LIBARG)

test.o: test.cpp 
		$(CILKPP) $(COMPILEARG) $<

debug: COMPILEARG += -DBPASDEBUG=1
debug: $(TARGET)

serial: COMPILEARG += -DSERIAL=1
serial: $(TARGET) test

test:
		./$(TARGET)

clean:
		rm -rf $(TARGET) *.out* *~ *.log *.o *.dat
