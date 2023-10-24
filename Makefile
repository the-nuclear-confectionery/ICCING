CC = g++
DEBUG = -g -O0
CFLAGS = -std=c++17 -Wall -c -fopenmp $(DEBUG)
LFLAGS = -Wall -fopenmp $(DEBUG)
MALLOC_CHECK = 2

# Define GSL flags separately
GSL_CFLAGS = `gsl-config --cflags`
GSL_LIBS = `gsl-config --libs`

CPP_FILES = ecc.cpp event.cpp io.cpp functions.cpp main.cpp splitting.cpp correlation.cpp greenfunction.cpp mask.cpp
HEADER_FILES = ecc.h event.h io.h functions.h splitting.h global.h correlation.h greenfunction.h mask.h
OBJECT_FILES = $(CPP_FILES:.cpp=.o)

# Add dependencies for object files
%.o: %.cpp $(HEADER_FILES)
	$(CC) $(CFLAGS) $(GSL_CFLAGS) $< -o $@

iccing: $(OBJECT_FILES)
	$(CC) $(LFLAGS) $(GSL_LIBS) $(OBJECT_FILES) -o iccing

clean:
	rm -f $(OBJECT_FILES) iccing
	echo Clean done
