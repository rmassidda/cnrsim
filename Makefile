CFLAGS = -Iinclude -Wall -O3 -g
LDFLAGS = -lhts -lm -ledlib -lz

VAROBJ = parse_frequency.c user_variation.c variator.c wrapper.c allele.c
ERROBJ = error_profiler.c translate_notation.c allele.c stats.c source.c model.c tandem.c
SIMOBJ = simulator.c stats.c source.c model.c tandem.c

variator: $(addprefix src/, ${VAROBJ})
	cc ${CFLAGS} -o $@ $^ ${LDFLAGS}

error: $(addprefix src/, ${ERROBJ})
	cc ${CFLAGS} -o $@ $^ ${LDFLAGS} 

simulator: $(addprefix src/, ${SIMOBJ})
	cc ${CFLAGS} -o $@ $^ ${LDFLAGS} 

.PHONY: clean all

all: variator error simulator

clean:
	-rm variator error simulator
