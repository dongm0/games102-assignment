#include "imgui/imgui.h"
#include "imgui_impl/imgui_impl_glfw.h"
#include "imgui_impl/imgui_impl_opengl3.h"
#include "glad/include/glad/glad.h"
#include <GLFW/glfw3.h>
#include "implot/implot.h"
#include "implot/implot_internal.h"

#include "implot/implot_demo.cpp"
#include "fitcurve.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>

static void glfw_error_callback(int error, const char* description) {
    std::cerr << "Glfw error"+std::to_string(error)+
                 ": "+std::string(description) << std::endl;
}

int main(int argc, char **argv) {
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;

    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); 

    GLFWwindow* window = glfwCreateWindow(1280, 720, "Hw1 by dongmo", NULL, NULL);
    if (window == nullptr)
        return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    if (gladLoadGL() == 0) {
        std::cerr << "Failed to init glad." << std::endl;
        return 1;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO &io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    ImVector<ImPlotPoint> data;
    NodeArr arr1, arr2, arr3, arr4, src;

    bool mainopen = true;
    bool draw = false;
    double mu = 1;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        {
            /*
            const float my_pi = 3.14159265;
            static float_t xs[1001], ys[1001];
            for (int i=0; i<1001; ++i) {
                xs[i] = cos(i*0.001f*my_pi);
                ys[i] = sin(i*0.001f*my_pi);//xs[i]*xs[i] - 2*xs[i] + 1;
            }
            if (ImPlot::BeginPlot("test", "x", "y")) {
                ImPlot::PlotLine("test func1", xs, ys, 1001);
                ImPlot::EndPlot();
            }
            */
            
            ImGui::SetNextWindowSize(ImVec2(1366,768), ImGuiCond_Appearing);
            ImGui::Begin("main", &mainopen);
            bool cal = false;
            if (ImPlot::BeginPlot("Fig", "x", "y", ImVec2(-3, 0))) {
                if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0)) {
                    ImPlotPoint pt = ImPlot::GetPlotMousePos();
                    data.push_back(pt);
                }
                if (data.size() > 0) {
                    ImPlot::PlotScatter("Raw Data", &data[0].x, &data[0].y, data.size(), 0, 2*sizeof(double));
                }
                else {
                    double xx = 0, yy = 0;
                    ImPlot::PlotScatter("Raw Data", &xx, &yy, 0);
                }
                
                if (draw == true && arr1.size>0) {
                    ImPlot::PlotLine("InterPolation_Poly", &arr1.xs[0], &arr1.ys[0], arr1.size, 0, sizeof(double));
                    ImPlot::PlotLine("InterPolation_Gaussian", &arr2.xs[0], &arr2.ys[0], arr2.size, 0, sizeof(double));
                }
                
                ImPlot::EndPlot();
            }
            if (ImGui::Button("Draw", ImVec2(50, 26))) {
                cal = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("Clear", ImVec2(50, 26))) {
                data.clear();
                draw = false;
            }
            ImGui::InputDouble("Input mu", &mu);
            
            if (cal == true && data.size()>0) {
                std::sort(data.begin(), data.end(), [](ImPlotPoint &p1, ImPlotPoint &p2){ return p1.x < p2.x; });
                src.xs.resize(data.size()), src.ys.resize(data.size());
                for (int i=0; i<data.size(); ++i) {
                    src.xs.at(i) = data[i].x;
                    src.ys.at(i) = data[i].y;
                }
                src.size = src.xs.size();
                InterPolation_1(src, arr1);
                InterPolation_2(src, arr2, 0.1);
                draw = true;
            }
            ImGui::End();
        
            /*
            bool kkkkk = true;
            ImPlot::ShowDemoWindow(&kkkkk);
            */
        }

        ImGui::Render();
        int _width, _height;
        glfwGetFramebufferSize(window, &_width, &_height);
        glViewport(0, 0, _width, _height);
        glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}