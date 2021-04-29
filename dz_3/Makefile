SHELL := /bin/bash
pnum = 4
qnum = 8
e = 0.01
appname = main
write = 1
read = 2
mpic++ = mpic++
mpirun = mpirun
n = -n
times = 60

all: main run

main: 
	$(mpic++) main.cpp -o $(appname) -std=c++11 -fopenmp

run: main
	$(mpirun) $(n) $(pnum) $(appname) $(qnum) $(e) $(write)

test: main tester
	$(mpirun) $(n) 4 $(appname) $(qnum) $(e) $(write)
	$(mpirun) $(n) 1 $(appname) $(qnum) $(e) $(read)
	$(mpirun) $(n) 2 $(appname) $(qnum) $(e) $(read)
	$(mpirun) $(n) 4 $(appname) $(qnum) $(e) $(read)
	$(mpirun) $(n) 8 $(appname) $(qnum) $(e) $(read)
	./test $(qnum)

tester:
	g++ test.cpp -o test

distribute: main run
	for ((i = 1; i <= ${times} ; i++)) do $(mpirun) $(n) $(pnum) $(appname) $(qnum) $(e) $(read); done

.PHONY: distribute

clean:
	rm $(appname)
	rm *.txt
	rm *.err
	rm *.out
	rm test
