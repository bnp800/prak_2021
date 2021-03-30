pnum = 8
qnum = 27
target = 27
appname = main
write = 1
read = 2
mpirun = mpisubmit.pl
n = -p

all: main run

main: 
	mpic++ main.cpp -o $(appname) -std=c++11 

run: main
	$(mpirun) $(n) $(pnum) $(appname) $(qnum) $(target) $(write)

test: main tester
	$(mpirun) $(n) 8 $(appname) $(qnum) $(target) $(write)
	$(mpirun) $(n) 1 $(appname) $(qnum) $(target) $(read)
	$(mpirun) $(n) 2 $(appname) $(qnum) $(target) $(read)
	$(mpirun) $(n) 4 $(appname) $(qnum) $(target) $(read)
	$(mpirun) $(n) 8 $(appname) $(qnum) $(target) $(read)
	./test $(qnum)

tester:
	g++ test.cpp -o test

clean:
	rm $(appname)
	rm *.txt
	rm *.err
	rm *.out
	rm test
