#makefile
all: essai

#init.o
obj/init.o: init.c include/init.h include/tools.h include/struct.h -lm
	gcc -O3 -Wall -Iinclude/ -c init.c -o obj/init.o -lm

#tools.o
obj/tools.o: tools.c include/tools.h include/init.h include/struct.h -lm
	gcc -O3 -Wall -Iinclude/ -c tools.c -o obj/tools.o -lm

#channel.o
obj/channel.o: channel.c include/channel.h include/init.h include/tools.h include/struct.h -lm
	gcc -O3 -Wall -Iinclude/ -c channel.c -o obj/channel.o -lm
	
#codec.o
obj/polarcodec.o:  init.c include/init.h include/tools.h include/struct.h -lm
	gcc -O3 -Wall -Iinclude/ -c polarcodec.c -o obj/polarcodec.o -lm

#simulator.o
obj/simulator.o: simulator.c include/init.h include/tools.h include/channel.h include/struct.h include/polarcodec.h -lm
	gcc -O3 -Wall -Iinclude/ -c simulator.c -o obj/simulator.o -lm


########## Generation de l'executable #######################
essai: obj/init.o obj/tools.o obj/channel.o obj/polarcodec.o obj/simulator.o -lm
	gcc  -O3 -Wall -lm -o essai obj/init.o obj/tools.o obj/channel.o obj/polarcodec.o obj/simulator.o -lm

#fin makefile
