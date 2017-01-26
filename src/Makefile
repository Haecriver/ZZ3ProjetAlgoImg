#On purge la liste des suffixes utilis� pour les r�les implicites
.SUFFIXES:

#On ajoute simplements les extensions dont l'on a besoin
.SUFFIXES:.cpp .o

#Nom de l'executable
EXEC=mosuofe

#Liste des fichiers sources separes par des espaces
SOURCES=main.cpp
HEADERS=flotOpti.hpp contAct.hpp

#Liste des fichiers objets
OBJETS=$(SOURCES:%.cpp=%.o)

#Compilateur et options de compilation
CCPP=g++
CFLAGS=-Wall -ansi -pedantic -ffast-math -I /usr/X11R6/include -I ./CImg -O3


LFLAGS= -L . -L /usr/X11R6/lib  -lpthread -lX11 -lXext -Dcimg_use_xshm  -lm -lpng

#R�le explicite de construction de l'ex�utable
$(EXEC):$(OBJETS) $(HEADERS) Makefile
	$(CCPP) -o  $(EXEC) $(OBJETS) $(LFLAGS)
.cpp.o:
	$(CCPP) $(CFLAGS) -c $< -o $@

run: $(EXEC)
	./$(EXEC)

clean:
	rm $(OBJETS)
clear:
	rm $(EXEC)
depend:
	sed -e "/^#DEPENDANCIES/,$$ d" Makefile >dependances
	echo "#DEPENDANCIES" >> dependances
	$(CCPP) -MM $(SOURCES) >> dependances
	cat dependances >Makefile
	rm dependances

#DEPENDANCIES
main.o: main.cpp CImg.h contAct.hpp flotOpti.hpp
