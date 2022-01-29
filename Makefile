FLAGS = -O3 -ftree-vectorize -fPIC -std=c++11 -fopenmp 
CXX = g++ $(FLAGS)
LIBNAME = src/grid_otflib.so

OBJ:= src/grid_otf.o
	
lib: $(OBJ)
	$(CXX) -shared -Wl,-install_name,$(LIBNAME) -o $(LIBNAME) $^

liblinux: $(OBJ)
	$(CXX) -shared -o $(LIBNAME) $^
	
$(OBJ): %.o : %.cpp
	$(CXX) -c -o $@ $< 
	
src/grid_otflib.so: $(OBJ)
	

clean:
	rm -rf src/*.o *.o 