CFLAGS = -Iinclude -Wall -lhts -lm -O3 -ledlib -lz

VAROBJ = parse_frequency.c user_variation.c variator.c wrapper.c allele.c
ERROBJ = error_profiler.c translate_notation.c allele.c stats.c source.c

variator: $(addprefix src/, ${VAROBJ})
	cc ${CFLAGS} -o $@ $^

error: $(addprefix src/, ${ERROBJ})
	cc ${CFLAGS} -o $@ $^

.PHONY: clean
clean:
	-rm variator error
