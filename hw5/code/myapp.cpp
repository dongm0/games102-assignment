#include "myapp.h"

static void glfw_error_callback(int error, const char* description) {
    std::cerr << "Glfw error"+std::to_string(error)+
                 ": "+std::string(description) << std::endl;
}

void MyApp::initAll() {
    glfwSetErrorCallback(glfw_error_callback);

    glfwInit();

    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); 

    window = glfwCreateWindow(1200, 1200, "Hw5 by dongmo", NULL, NULL);

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    gladLoadGL();

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO &io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);
}

void MyApp::cleanUp() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
}

void MyApp::mainLoop() {
    //预定变量
    ControlPointArray2D arr;
    bool edit = false;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        {            
            ImGui::SetNextWindowSize(ImVec2(1000,1000), ImGuiCond_Appearing);
            ImGui::Begin("main");

            if (ImPlot::BeginPlot("", "", "", ImVec2(900, 900))) {
                if (ImPlot::IsPlotHovered() and ImGui::IsMouseClicked(0) and !arr.closed()) {
                    ImPlotPoint pt = ImPlot::GetPlotMousePos();
                    arr.pushBack({pt.x, pt.y});
                }
                else if (ImPlot::IsPlotHovered() and ImGui::IsMouseDoubleClicked(0)) {
                    arr.makeClose();
                }
            }
            {
                auto cp = arr.getCtrlPoint();
                auto dp = arr.getDrawPoint();

                if (cp.size() > 0)
                    ImPlot::PlotLine("control points", &cp[0].x, &cp[0].y, cp.size(), 0, 2*sizeof(double));
                if (dp.size() > 0)
                    ImPlot::PlotLine("curve", &dp[0].x, &dp[0].y, dp.size(), 0, 2*sizeof(double));

            }
            ImPlot::EndPlot();

            ImGui::Checkbox("edit ctrl points", &edit);
            int it = 0;
            ImGui::SliderInt("curve iteration times", &it, 0, 10);
            arr.setItTime(it);
            float l = 0.12;
            ImGui::SliderFloat("curve iteration times", &l, 0, 0.125f);
            arr.setInterpolationPara(l);
            int me = 0;
            ImGui::RadioButton("2nd B-spline", &me, 0);ImGui::SameLine();
            ImGui::RadioButton("3rd B-spline", &me, 1);ImGui::SameLine();
            ImGui::RadioButton("Intrtpolation", &me, 2);
            arr.setSegmentation(Segmentation(me));

            ImGui::End();
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
}