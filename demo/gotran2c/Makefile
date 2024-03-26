sources = demo.cpp
executable = demo
libfile = libdemo.so

# CPPFLAGS = -Wall -Wno-unknown-pragmas -O3 -fopenmp -ffast-math -march=native -DHINT_OMP_SIMD -DVECTOR_LENGTH=4
# CPPFLAGS = -Wall -O3 -fopenmp -ffast-math -march=native -DHINT_OMP_SIMD -DVECTOR_LENGTH=4
# CPPFLAGS = -Wall -O3 -fopenmp -ffast-math -march=native -DVECTOR_LENGTH=4
# CPPFLAGS = -Wall -O3 -fopenmp -ffast-math -march=native -DHINT_OMP_SIMD
CPPFLAGS = -Wall -O3 -fopenmp -ffast-math -march=native
# CPPFLAGS = -Wall -O3 -ffast-math -march=native


# LDLIBS = -lm -fopenmp -lstdc++ -g
# sudo sysctl kernel.perf_event_paranoid=-1
# sudo sysctl kernel.kptr_restrict=0

# perf record -e task-clock ./demo 10000
# perf report


LDLIBS = -lm -fopenmp -lstdc++

# gprof demo gmon.out > analysis.txt


# Working models:
#-----------------------------------------------
# MODEL_NAME = beeler_reuter_1977_version06
# MODEL_NAME = GPB
# MODEL_NAME = grandi_pasqualini_bers_2010
# MODEL_NAME = iyer_2004
# MODEL_NAME = JT21
# MODEL_NAME = maleckar_2009
# MODEL_NAME = maltsev_2009
# MODEL_NAME = pandit_et_al_2001_endo
# MODEL_NAME = Pandit_Hinch_Niederer_Model
# MODEL_NAME = shannon_2004
# MODEL_NAME = tentusscher_2004_mcell
MODEL_NAME = tentusscher_panfilov_2006_M_cell
# MODEL_NAME = terkildsen_2008
# MODEL_NAME = winslow_model_1999
#----------------------------------------------



# Working models, without lookup
#----------------------------------
# MODEL_NAME = FHN
# MODEL_NAME = hodgkin_huxley

# MODEL_NAME = irvine_model_1999

# MODEL_NAME = Hinch_et_al_2004
# MODEL_NAME = niederer_2006
# MODEL_NAME = Niederer_et_al_2006
# MODEL_NAME = rice_model_2008
#----------------------------------




MODEL_PATH = models


MODEL_ODE_FILE = $(MODEL_PATH)/$(MODEL_NAME).ode
MODEL_HPP_FILE = $(MODEL_NAME).hpp

candidates = ""

gotran_generator = python -m gotran gotran2c 

all: clean $(executable) $(libfile)
$(executable): demo.o



$(libfile): CPPFLAGS += -fPIC
$(libfile): $(sources) $(MODEL_HPP_FILE)
	gcc $(CPPFLAGS) -shared -o $@ $< $(LDLIBS)


demo.o: $(MODEL_HPP_FILE)
$(MODEL_HPP_FILE): $(MODEL_ODE_FILE)
	$(gotran_generator) $< $(candidates) --solvers.explicit_euler.generate=1 --solvers.rush_larsen.generate=1 --solvers.generalized_rush_larsen.generate=1 --code.body.use_enum=1 --output=$@


#$(gotran_generator) $< $(candidates) --solvers.generalized_rush_larsen.generate=1 --code.body.use_enum=1 --output=$@


clean:
	$(RM) $(libfile) $(executable) $(wildcard *.o) $(MODEL_HPP_FILE)

