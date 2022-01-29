FLAGS = -O3 -ftree-vectorize -fPIC -std=c++11 -fopenmp 
CXX = g++ $(FLAGS)

OBJ:= src/grid_otf.o
	
lib: src/grid_otflib.so
	
$(OBJ): %.o : %.cpp
	$(CXX) -c -o $@ $< 
	
src/grid_otflib.so: $(OBJ)
	$(CXX) -shared -Wl,-install_name,$@ -o $@ $^

clean:
	rm -rf src/*.o *.o 