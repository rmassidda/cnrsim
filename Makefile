CFLAGS = -Iinclude -Wall -lhts -lm -g

variator: src/*.c
	cc ${CFLAGS} -Wall -g -lhts -lm -o $@ $^

.PHONY: clean
clean:
	-rm variator
