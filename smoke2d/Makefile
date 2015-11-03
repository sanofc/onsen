.DEFAULT_GOAL := all 
all: clean build run
build: 
	clang -c *.c -framework GLUT -framework OpenGL -Wno-deprecated
	clang -o smoke2d *.o -framework GLUT -framework OpenGL -Wno-deprecated
run:
	./smoke2d
clean:
	rm -f smoke2d 
