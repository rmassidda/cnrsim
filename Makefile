CFLAGS = -Iinclude -Wall -lhts -lm -g

VAROBJ = fileManager.c parse_frequency.c user_variation.c variator.c wrapper.c

variator: $(addprefix src/, ${VAROBJ})
	cc ${CFLAGS} -Wall -g -lhts -lm -o $@ $^

.PHONY: clean
clean:
	-rm variator
