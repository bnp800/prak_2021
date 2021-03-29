pnum = 8
qnum = 15
target = 1
appname = main
write = 1
read = 2
all: main run

main: 
	mpic++ main.cpp -o $(appname) -std=c++11 

run: main
	mpirun -np $(pnum) $(appname) $(qnum) $(target) $(write)

test: main tester
	mpirun -np 8 $(appname) $(qnum) $(target) $(write)
	mpirun -np 1 $(appname) $(qnum) $(target) $(read)
	mpirun -np 2 $(appname) $(qnum) $(target) $(read)
	mpirun -np 4 $(appname) $(qnum) $(target) $(read)
	mpirun -np 8 $(appname) $(qnum) $(target) $(read)
	./test $(qnum)

tester:
	g++ test.cpp -o test

clean:
	rm $(appname)
	rm *.txt
