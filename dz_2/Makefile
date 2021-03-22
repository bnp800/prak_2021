pnum = 2
qnum = 20
target = 9
appname = main
all: main run

main: 
	mpic++ main.cpp -o $(appname)

run:
	mpirun -np $(pnum) $(appname) $(qnum) $(target)

test:
	mpirun -np 1 $(appname) $(qnum) $(target)
	mpirun -np 2 $(appname) $(qnum) $(target)
	mpirun -np 4 $(appname) $(qnum) $(target)
	mpirun -np 8 $(appname) $(qnum) $(target)

clean:
	rm $(appname)
	rm *.txt
