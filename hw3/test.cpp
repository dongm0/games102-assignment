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

    GLFWwindow* window = glfwCreateWindow(1200, 1000, "Hw3 by dongmo", NULL, NULL);
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
    std::vector<double> para1;
    std::vector<double> ones;

    bool mainopen = true;
    bool draw = false;
    double mu = 1;
    int order = 5;
    double lam = 1;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        {            
            ImGui::SetNextWindowSize(ImVec2(800,900), ImGuiCond_Appearing);
            ImGui::Begin("main", &mainopen);
            bool cal = false;
            if (ImPlot::BeginPlot("Fig", "x", "y", ImVec2(600, 600))) {
                if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0)) {
                    ImPlotPoint pt = ImPlot::GetPlotMousePos();
                    data.push_back(pt);
                }
                if (data.size() > 0) {
                    ImPlot::PlotScatter("Raw Data", &data[0].x, &data[0].y, data.size(), 0, 2*sizeof(double));
                }
                
                if (draw == true) {
                    ImPlot::PlotLine("InterPolation_Gaussian", &arr1.xs[0], &arr1.ys[0], arr1.size, 0, sizeof(double));
                    //ImPlot::PlotLine("InterPolation_Gaussian", &arr2.xs[0], &arr2.ys[0], arr2.size, 0, sizeof(double));
                    //ImPlot::PlotLine("Least_Square", &arr3.xs[0], &arr3.ys[0], arr3.size, 0, sizeof(double));
                    //ImPlot::PlotLine("Ridge_Regression", &arr4.xs[0], &arr4.ys[0], arr4.size, 0, sizeof(double));
                }
                
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Para", "t", "", ImVec2(600, 120), 0, 0, 6)) {
                if (draw == true) {
                    ImPlot::PlotScatter("", &para1[0], &ones[0], para1.size(), 0, sizeof(double));
                    ImPlot::PlotLine("", &para1[0], &ones[0], para1.size(), 0, sizeof(double));
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
            static int e = 0;
            ImGui::RadioButton("uniform", &e, 0); 
            ImGui::SameLine();
            ImGui::RadioButton("chordal", &e, 1); 
            ImGui::SameLine();
            ImGui::RadioButton("centripetal", &e, 2);
            ImGui::SameLine();
            ImGui::RadioButton("foley", &e, 3);
            
            if (cal == true && data.size()>0) {
                //std::sort(data.begin(), data.end(), [](ImPlotPoint &p1, ImPlotPoint &p2){ return p1.x < p2.x; });
                src.xs.resize(data.size()), src.ys.resize(data.size());
                for (int i=0; i<data.size(); ++i) {
                    src.xs.at(i) = data[i].x;
                    src.ys.at(i) = data[i].y;
                }
                src.size = src.xs.size();
                
                if (e == 0) {
                    para1 = Parametrization_uniform(src);
                }
                else if (e == 1) {
                    para1 = Parametrization_chordal(src);
                }
                else if (e == 2) {
                    para1 = Parametrization_centripetal(src);
                }
                else {
                    para1 = Parametrization_foley(src);
                }
                ones.assign(para1.size(), 0.5);
                NodeArr src_x(para1, src.xs);
                NodeArr src_y(para1, src.ys);
                //InterPolation_1(src, arr1);
                //InterPolation_2(src, arr2, mu);
                InterPolation_2(src_x, arr1);
                InterPolation_2(src_y, arr2);
                arr1.xs = arr1.ys, arr1.ys = arr2.ys;
                //Ridge_Regression(src, arr4, order, lam);
                draw = true;
                cal = false;
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