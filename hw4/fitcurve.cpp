#include "fitcurve.h"

void ControlPointArray2D::ThreeOrder() {
    int st = 0, ed = 0;
    for (int i=(int)nodenum()-1; i>=0; --i) {
        if (fixed[i] == true) {
            st = i;
            break;
        }
    }
    auto xparm = ThreeOrderSample(xs.begin()+st, xs.end());
    auto yparm = ThreeOrderSample(ys.begin()+st, ys.end());

    for (int i=st; i<ed-1; ++i) {
        while (param.size()<=i) {
            param.push_back({});
        }
        param[i] = {xparm[i*4], xparm[i*4+1], xparm[i*4+2], xparm[i*4+3], 
                    yparm[i*4], yparm[i*4+1], yparm[i*4+2], yparm[i*4+3]};
    }

}

void ControlPointArray2D::calculateRange(int p) {
    
}