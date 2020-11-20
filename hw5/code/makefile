OBJ = fitcurve.o main.o glad.o myapp.o 
OBJ += ./imgui/imgui_draw.o ./imgui/imgui_impl_glfw.o
OBJ += ./imgui/imgui_impl_opengl3.o ./imgui/imgui_widgets.o ./imgui/imgui.o
OBJ += ./imgui/implot_items.o ./imgui/implot.o 

CXXFLAGS = -I ./glad -I /usr/local/include
CXXFLAGS += -DIMGUI_IMPL_OPENGL_LOADER_GLAD -O2

LIBS = -lGL -lm -lrt -ldl -L/usr/local/lib -lglfw

hw5:$(OBJ)
	g++ $(OBJ) $(CXXFLAGS) $(LIBS) -o hw5

glad.o:glad.c
	gcc $(CXXFLAGS) -c -o $@ $<

%.o:%.cpp
	g++ $(CXXFLAGS) -c -o $@ $<

%.o:./imgui/%.cpp
	g++ $(CXXFLAGS) -c -o $@ $<

cleanOBJS:
	rm $(OBJ)

clean:
	rm hw5 $(OBJ)