all:
	cc -std=c11 -I . cdol.c main.c -g -o cdol_test.bin

clean:
	rm -fv cdol_test.bin
	rm -rfv *.dSYM
