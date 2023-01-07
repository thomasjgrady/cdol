all:
	cc -std=c99 -I . cdol.c main.c -g3 -o cdol_test.bin

clean:
	rm -fv cdol_test.bin
	rm -rfv *.dSYM
