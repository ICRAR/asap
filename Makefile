# Makefile to make integration into IDE's easier as they don't (yet) support scons

all:
	@$(SCONS)

clean:
	@$(SCONS) -c
