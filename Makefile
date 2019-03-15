CFLAGS = -Iinclude -Wall -lhts -lm -g -ledlib

VAROBJ = fileManager.c parse_frequency.c user_variation.c variator.c wrapper.c allele.c
ERROBJ = align.c error_profiler.c fileManager.c translate_notation.c allele.c

variator: $(addprefix src/, ${VAROBJ})
	cc ${CFLAGS} -o $@ $^

error: $(addprefix src/, ${ERROBJ})
	cc ${CFLAGS} -o $@ $^

.PHONY: clean
clean:
	-rm variator error
