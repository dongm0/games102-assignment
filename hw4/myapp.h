#include "imgui/imgui.h"
#include "imgui_impl/imgui_impl_glfw.h"
#include "imgui_impl/imgui_impl_opengl3.h"
#include "glad/include/glad/glad.h"
#include <GLFW/glfw3.h>
#include "implot/implot.h"
#include "implot/implot_internal.h"

#include "fitcurve.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>

class MyApp {
public:
    void run();
private:
    GLFWwindow *window;

    void initAll();
    void mainLoop();
    void cleanUp();
};