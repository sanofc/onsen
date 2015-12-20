.DEFAULT_GOAL := all 
all: clean build
build: 
	clang++ -c *.cpp -framework GLUT -framework OpenGL -Wno-deprecated
	clang++ -o steam2d *.o -framework GLUT -framework OpenGL -Wno-deprecated
run:
	./steam2d
clean:
	rm -f steam2d 
