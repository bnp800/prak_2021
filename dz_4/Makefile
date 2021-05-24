qnum = 24
pnum = 4
target1 = 9
target2 = 7
mpirun = mpirun
mpic++ = mpic++
n = -n

generator: generator.cpp
	$(mpic++) generator.cpp -o generator

generate: generator
	$(mpirun) $(n) $(pnum) generator $(qnum) 

check: check.cpp
	$(mpic++) check.cpp -o check


hadamard: hadamard.cpp hadamard_test.cpp
	$(mpic++) hadamard.cpp -o hadamard
	$(mpic++) hadamard_test.cpp -o hadamard_test

test_hadamard: hadamard
	$(mpirun) $(n) $(pnum) hadamard $(qnum) $(target1)
	$(mpirun) $(n) $(pnum) hadamard_test $(qnum) $(target1)

check_hadamard: check generate test_hadamard
	$(mpirun) $(n) $(pnum) check $(qnum) H_res$(pnum).txt.t


n-hadamard: n-hadamard.cpp n-hadamard_test.cpp
	$(mpic++) n-hadamard.cpp -o n-hadamard
	$(mpic++) n-hadamard_test.cpp -o n-hadamard_test

test_n-hadamard: n-hadamard
	$(mpirun) $(n) $(pnum) n-hadamard $(qnum) 
	$(mpirun) $(n) $(pnum) n-hadamard_test $(qnum) 

check_n-hadamard: check generate test_n-hadamard 
	$(mpirun) $(n) $(pnum) check $(qnum) n-H_res$(pnum).txt.t


NOT: NOT.cpp NOT_test.cpp
	$(mpic++) NOT.cpp -o NOT
	$(mpic++) NOT_test.cpp -o NOT_test

test_NOT: NOT
	$(mpirun) $(n) $(pnum) NOT $(qnum) $(target1)
	$(mpirun) $(n) $(pnum) NOT_test $(qnum) $(target1)

check_NOT: check generate test_NOT
	$(mpirun) $(n) $(pnum) check $(qnum) NOT_res$(pnum).txt.t

CNOT: CNOT.cpp CNOT_test.cpp
	$(mpic++) CNOT.cpp -o CNOT
	$(mpic++) CNOT_test.cpp -o CNOT_test

test_CNOT: CNOT
	$(mpirun) $(n) $(pnum) CNOT $(qnum) $(target1) $(target2)
	$(mpirun) $(n) $(pnum) CNOT_test $(qnum) $(target1) $(target2)

check_CNOT: check generate test_CNOT
	$(mpirun) $(n) $(pnum) check $(qnum) CNOT_res$(pnum).txt.t

ROT: ROT.cpp ROT_test.cpp
	$(mpic++) ROT.cpp -o ROT
	$(mpic++) ROT_test.cpp -o ROT_test

test_ROT: ROT
	$(mpirun) $(n) $(pnum) ROT $(qnum) $(target1)
	$(mpirun) $(n) $(pnum) ROT_test $(qnum) $(target1)

check_ROT: check generate test_ROT
	$(mpirun) $(n) $(pnum) check $(qnum) ROT_res$(pnum).txt.t

CROT: CROT.cpp CROT_test.cpp
	$(mpic++) CROT.cpp -o CROT
	$(mpic++) CROT_test.cpp -o CROT_test

test_CROT: CROT
	$(mpirun) $(n) $(pnum) CROT $(qnum) $(target1) $(target2)
	$(mpirun) $(n) $(pnum) CROT_test $(qnum) $(target1) $(target2)

check_CROT: check generate test_CROT
	$(mpirun) $(n) $(pnum) check $(qnum) CROT_res$(pnum).txt.t

clean:
	rm *txt
	rm *.t
	rm check
	rm generator
	rm hadamard
	rm n-hadamard
	rm NOT
	rm CNOT
	rm ROT
	rm CROT
	rm *_test