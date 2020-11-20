#include "fitcurve.h"



int ControlPointArray2D::getClosePoint(mypoint p) {
    int res = -1;
    double minv = 1e10;
    for (int i = 0; i < sizeC(); ++i) {
        double dist = sqrt(pow(p.x - m_ctrlPoints[i].x, 2) + pow(p.y - m_ctrlPoints[i].y, 2));
        if (dist < minv and dist < m_closePointDist) {
            minv = dist;
            res = i;
        }
    }
    return res;
}

void ControlPointArray2D::updateDrawPoints() {
    if (m_method == Segmentation::B_SPLINE_2ND) {
        updateDrawPoints_1();
    }
    else if (m_method == Segmentation::B_SPLINE_3RD) {
        updateDrawPoints_2();
    }
    else if (m_method == Segmentation::INTERPOLATION) {
        updateDrawPoints_3();
    }
    else {}
}

void ControlPointArray2D::updateDrawPoints_1() {
    if (sizeC() <= 1)
        return;
    auto arrcpy = m_ctrlPoints;
    for (int i=0; i<m_itTime; ++i) {
        std::vector<mypoint> tmp(arrcpy.size()*2-2);
        for (int j=0; j<arrcpy.size()-1; ++j) {
            tmp[2*j] = arrcpy[j]*0.75 + arrcpy[j+1]*0.25;
            tmp[2*j+1] = arrcpy[j]*0.25 + arrcpy[j+1]*0.75;
        }
        if (m_closed) {
            tmp.push_back(arrcpy.back()*0.75 + arrcpy.front()*0.25);
            tmp.push_back(arrcpy.back()*0.25 + arrcpy.front()*0.75);
        }
        arrcpy = std::move(tmp);
    }
    m_drawPoints = arrcpy;
}

void ControlPointArray2D::updateDrawPoints_2() {
    if (sizeC() <= 2)
        return;
    auto arrcpy = m_ctrlPoints;
    for (int i=0; i<m_itTime; ++i) {
        std::vector<mypoint> tmp(arrcpy.size()*2-3);
        for (int j=0; j<arrcpy.size()-1; ++j) {
            tmp[2*j] = arrcpy[j]*0.5 + arrcpy[j+1]*0.5;
            if (j != arrcpy.size()-2)
                tmp[2*j+1] = arrcpy[j]*0.125 + arrcpy[j+1]*0.75 + arrcpy[j+2]*0.125;
        }
        if (m_closed) {
            tmp.push_back(arrcpy[arrcpy.size()-2]*0.125 + arrcpy[arrcpy.size()-1]*0.75 + arrcpy[0]*0.125);
            tmp.push_back(arrcpy[arrcpy.size()-1]*0.5 + arrcpy[0]*0.5);
            tmp.push_back(arrcpy[arrcpy.size()-1]*0.125 + arrcpy[0]*0.75 + arrcpy[1]*0.125);
        }
        arrcpy = std::move(tmp);
    }
    m_drawPoints = arrcpy;
}

void ControlPointArray2D::updateDrawPoints_3() {
    if (sizeC() <= 4)
        return;
    auto arrcpy = m_ctrlPoints;
    for (int i=0; i<m_itTime; ++i) {
        if (!m_closed) {
            std::vector<mypoint> tmp(arrcpy.size()*2-1);
            mypoint dumph, dumpt;
            dumph = arrcpy[0]*2 - arrcpy[1];
            dumpt = arrcpy[arrcpy.size()-1]*2 - arrcpy[arrcpy.size()-2];
            for (int j=0; j<arrcpy.size(); ++j) {
                tmp[j*2] = arrcpy[j];
                if (j != arrcpy.size()-1) {
                    if (j == 0) {
                        tmp[1] = ((arrcpy[0]+arrcpy[1])*(1+m_interpolationPara) - 
                                  (dumph+arrcpy[2])*m_interpolationPara)*(0.5);
                    }
                    else if (j == arrcpy.size()-2) {
                        tmp[j*2+1] = ((arrcpy[j]+arrcpy[j+1])*(1+m_interpolationPara) - 
                                      (arrcpy[j-1]+dumpt)*m_interpolationPara)*(0.5);
                    }
                    else {
                        tmp[j*2+1] = ((arrcpy[j]+arrcpy[j+1])*(1+m_interpolationPara) - 
                                      (arrcpy[j-1]+arrcpy[j+2])*m_interpolationPara)*(0.5);
                    }
                }
            }
            arrcpy = std::move(tmp);
        }
        else {
            std::vector<mypoint> tmp(arrcpy.size()*2);
            int len = arrcpy.size();
            for (int j=0; j<arrcpy.size(); ++j) {
                tmp[j*2] = arrcpy[j];
                tmp[j*2+1] = ((arrcpy[j]+arrcpy[(j+1+len)%len])*(1+m_interpolationPara) - 
                              (arrcpy[(j-1+len)%len]+arrcpy[(j+2+len)%len])*m_interpolationPara)*(0.5);
            }
            arrcpy = std::move(tmp);
        }
    }
    m_drawPoints = arrcpy;
}