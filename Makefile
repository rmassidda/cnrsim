objects = e00 e01 e02 e03

experiments: $(objects)

e00: experiment_00.c fileManager.o 
	cc -Wall -o $@ $^

e01: experiment_01.c
	cc -Wall -g -lhts -o $@ $< 

e02: experiment_02.c
	cc -Wall -g -lhts -o $@ $< 
	         
e03: experiment_03.c fileManager.o
	cc -Wall -g -lhts -o $@ $^ 

fileManager.o: fileManager.h common.h 

.PHONY: clean
clean:
	-rm $(objects) *.o
