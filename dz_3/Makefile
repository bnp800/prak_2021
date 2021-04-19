pnum = 2
qnum = 24
e = 0.01
appname = main
write = 1
read = 2
mpirun = mpirun
n = -n

all: main run

main: 
	mpic++ main.cpp -o $(appname) -std=c++11 -fopenmp

run: main
	$(mpirun) $(n) $(pnum) $(appname) $(qnum) $(e) $(write)

test: main tester
	$(mpirun) $(n) 8 $(appname) $(qnum) $(e) $(write)
	$(mpirun) $(n) 1 $(appname) $(qnum) $(e) $(read)
	$(mpirun) $(n) 2 $(appname) $(qnum) $(e) $(read)
	$(mpirun) $(n) 4 $(appname) $(qnum) $(e) $(read)
	$(mpirun) $(n) 8 $(appname) $(qnum) $(e) $(read)
	./test $(qnum)

tester:
	g++ test.cpp -o test

clean:
	rm $(appname)
	rm *.txt
	rm *.err
	rm *.out
	rm test
