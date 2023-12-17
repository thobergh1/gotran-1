sources = demo.cpp
executable = demo
libfile = libdemo.so

CPPFLAGS = -Wall -Wno-unknown-pragmas -O3 -fopenmp -ffast-math -march=native -DHINT_OMP_SIMD -DVECTOR_LENGTH=4
LDLIBS = -lm -fopenmp -lstdc++


gotran_generator = python -m gotran gotran2c

all: $(executable) $(libfile)
$(executable): demo.o


$(libfile): CPPFLAGS += -fPIC
$(libfile): $(sources) hodgkin_huxley.hpp
	g++ $(CPPFLAGS) -shared -o $@ $< $(LDLIBS)

demo.o: hodgkin_huxley.hpp
hodgkin_huxley.hpp: hodgkin_huxley.ode
	$(gotran_generator) $< --solvers.explicit_euler.generate=1 --code.body.use_enum=1 --output=$@

clean:
	$(RM) $(libfile) $(executable) $(wildcard *.o) hodgkin_huxley.hpp