Qnum = 20
Target = 1
list = 1 2 4 8
all: main run clean

main: 
	g++ main.cpp -o main -fopenmp

run:
	./main
clean:
	rm main
test: main 
	./main $(Qnum) 1 $(Target)
	./main $(Qnum) 2 $(Target)
	./main $(Qnum) 4 $(Target)
	./main $(Qnum) 8 $(Target)


	

