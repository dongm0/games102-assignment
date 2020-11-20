#include <math.h>
#include <vector>

enum Segmentation {
    B_SPLINE_2ND,
    B_SPLINE_3RD,
    INTERPOLATION
};

struct mypoint {
    double x;
    double y;
    mypoint operator+(mypoint p1) {return{p1.x+x, p1.y+y};}
    mypoint operator-(mypoint p1) {return{x-p1.x, y-p1.y};}
    mypoint operator*(double a) {return{a*x, a*y};}
};



class ControlPointArray2D {
public:
    bool closed() {return m_closed;}
    int sizeC() {return m_ctrlPoints.size();}
    int sizeD() {return m_drawPoints.size();}
    int getClosePoint(mypoint p);
    std::vector<mypoint>& getDrawPoint() {return m_drawPoints;}
    std::vector<mypoint>& getCtrlPoint() {return m_ctrlPoints;}
    
    void makeClose() {
        if (sizeC() >= 2) 
            m_closed = true;
        updateDrawPoints();
    }
    void pushBack(mypoint p) {
        m_ctrlPoints.push_back(p);
        updateDrawPoints();
    }
    void setPos(int i, mypoint newp) {
        m_ctrlPoints.at(i) = newp;
        updateDrawPoints();
    }

    void setItTime(int time) {
        m_itTime = time;
        updateDrawPoints();
    }
    void setSegmentation(Segmentation m) {
        m_method = m;
        updateDrawPoints();
    }
    void setInterpolationPara(double l) {
        m_interpolationPara = l;
        updateDrawPoints();
    }

    void clear() {
        m_closed = false;
        m_ctrlPoints.clear();
        m_drawPoints.clear();
    }


private:
    void updateDrawPoints();
    void updateDrawPoints_1();
    void updateDrawPoints_2();
    void updateDrawPoints_3();

private:
    bool m_closed = false;
    std::vector<mypoint> m_ctrlPoints;
    std::vector<mypoint> m_drawPoints;

    double m_closePointDist = 0.02;

    //控制层
    int m_itTime = 0;
    Segmentation m_method = Segmentation::B_SPLINE_2ND;
    double m_interpolationPara = 0.05;

};
