exp = e00 e01 e02 e03 e04 e05

variator: variator.c fileManager.o parse_frequency.o
	cc -Wall -g -lhts -lm -o $@ $^

experiments: $(exp)

e00: experiment_00.c fileManager.o 
	cc -Wall -o $@ $^

e01: experiment_01.c
	cc -Wall -g -lhts -o $@ $< 

e02: experiment_02.c parse_frequency.o
	cc -Wall -g -lhts -o $@ $^ 
	         
e03: experiment_03.c fileManager.o
	cc -Wall -g -lhts -o $@ $^ 

e04: experiment_04.c fileManager.o
	cc -Wall -g -lhts -o $@ $^ 

e05: experiment_05.c fileManager.o
	cc -Wall -g -lhts -o $@ $^

fileManager.o: fileManager.h common.h 

parse_frequency.o: parse_frequency.h

.PHONY: clean
clean:
	-rm $(exp) variator *.o
