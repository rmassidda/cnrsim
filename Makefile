example: example.c fileManager.o 
	cc -Wall -o example example.c fileManager.o
	         
fileManager.o: lib/fileManager.h lib/common.h 
	cc -Wall -c lib/fileManager.c

clean:
	rm *.o example