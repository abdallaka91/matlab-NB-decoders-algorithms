#makefile
all: essai

#init.o
obj/init.o: init.c include/init.h include/tools.h include/struct.h  -lm
	gcc -O3 -Wall -Iinclude/ -c init.c -o obj/init.o

#tools.o
obj/tools.o: tools.c include/tools.h include/init.h include/struct.h -lm
	gcc -O3 -Wall -Iinclude/ -c tools.c -o obj/tools.o

#channel.o
obj/channel.o: channel.c include/channel.h include/init.h include/tools.h include/struct.h -lm
	gcc -O3 -Wall -Iinclude/ -c channel.c -o obj/channel.o
	
#codec.o
obj/polarcodec.o:  polarcodec.c include/polarcodec.h include/init.h include/tools.h include/struct.h -lm
	gcc -O3 -Wall -Iinclude/ -c polarcodec.c -o obj/polarcodec.o

#simulator.o
obj/simulator.o: simulator.c include/init.h include/tools.h include/channel.h include/struct.h include/polarcodec.h -lm
	gcc -O3 -Wall -Iinclude/ -c simulator.c -o obj/simulator.o


########## Generation de l'executable #######################
essai: obj/init.o obj/tools.o obj/channel.o obj/polarcodec.o obj/simulator.o 
	gcc  -O3 -Wall -lm -o essai obj/init.o obj/tools.o obj/channel.o obj/polarcodec.o obj/simulator.o -lm

#fin makefile
