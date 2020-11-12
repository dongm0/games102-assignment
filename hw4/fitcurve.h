#include <math.h>
#include <vector>
#include <eigen3/Eigen/Dense>



struct controlPoint {
    double val;
    double ldiff;
    double rdiff;
    bool fixed_diff = false;
};

void ThreeOrderSample(std::vector<controlPoint> &points, std::vector<controlPoint>::iterator begin, std::vector<controlPoint>::iterator end) {
    using namespace Eigen;

    int pointnum = end-begin;
    if (pointnum < 2) return;
    int n = pointnum - 1;
    
    Matrix<double, Dynamic, Dynamic> A(n, n);
    VectorXd b(n), x(x);
    int disp = 1;
    for (auto it=begin; it!=end; ++it) {
        int num = it-begin;

        //lvalue
        A(num*4+disp, num*4) = pow(num, 3);
        A(num*4+disp, num*4+1) = pow(num, 2);
        A(num*4+disp, num*4+2) = num;
        A(num*4+disp, num*4+3) = 1;
        b(num*4+disp) = it->val;

        //rvalue
        A(num*4+1+disp, num*4) = pow(num+1, 3);
        A(num*4+1+disp, num*4+1) = pow(num+1, 2);
        A(num*4+1+disp, num*4+2) = num+1;
        A(num*4+1+disp, num*4+3) = 1;
        b(num*4+1+disp) = (it+1)->val;

        if (num != pointnum - 1) {
            A(num*4+2+disp, num*4) = 3*pow(num+1, 2);
            A(num*4+2+disp, num*4+1) = 2*(num+1);
            A(num*4+2+disp, num*4+2) = 1;
            A(num*4+2+disp, (num+1)*4) = -3*pow(num+1, 2);
            A(num*4+2+disp, (num+1)*4+1) = -2*(num+1);
            A(num*4+2+disp, (num+1)*4+2) = -1;

            A(num*4+3+disp, num*4) = 3*(num+1);
            A(num*4+3+disp, num*4+1) = 1;
            A(num*4+3+disp, (num+1)*4) = -3*(num+1);
            A(num*4+3+disp, (num+1)*4+1) = -1;
        }

    }
    if (begin->fixed_diff == true) {
        A(0, 2) = 1;
        b(0) = begin->rdiff;
    }
    else {
        A(0, 1) = 1;
    }
    if ((end-1)->fixed_diff == true) {
        A(n*4-1, (n-1)*4) = 3*pow(n, 2);
        A(n*4-1, (n-1)*4+1) = 2*n;
        A(n*4-1, (n-1)*4+2) = 1;
        b(n*4-1) = (end-1)->ldiff;
    }
    else {
        A(n*4-1, (n-1)*4) = 6*n;
        A(n*4-1, (n-1)*4+1) = 2;
    }
    x = A.colPivHouseholderQr().solve(b);
}

class ControlPoints2D {
public:

private:
    std::vector<controlPoint> xs;
    std::vector<controlPoint> ys;
    std::vector<bool> fixed;
}

struct NodeArr {
    NodeArr() {}
    NodeArr(std::vector<double> x, std::vector<double> y) : xs(x), ys(y), size(x.size()) {}
    unsigned int size;
    std::vector<double> xs;
    std::vector<double> ys;
};

void Normalize(std::vector<double> &arr) {
    double len = arr.back();
    for (auto &x : arr)
        x /= len;
}


std::vector<double> Parametrization_uniform(const NodeArr &src) {
    using namespace std;
    assert(src.size > 1);
    vector<double> res(src.size, 0);
    for (int i=0; i<res.size(); ++i) {
        res[i] = i;
    }
    Normalize(res);
    return res;
}

std::vector<double> Parametrization_chordal(const NodeArr &src) {
    using namespace std;
    assert(src.size > 1);
    vector<double> res(src.size, 0);
    for (int i=1; i<res.size(); ++i) {
        res[i] = res[i-1] + sqrt(pow(src.xs[i]-src.xs[i-1], 2) + 
                                pow(src.ys[i]-src.ys[i-1], 2));
    }
    Normalize(res);
    return res;
}

std::vector<double> Parametrization_centripetal(const NodeArr &src) {
    using namespace std;
    assert(src.size > 1);
    vector<double> res(src.size, 0);
    for (int i=1; i<res.size(); ++i) {
        res[i] = res[i-1] + sqrt(sqrt(pow(src.xs[i]-src.xs[i-1], 2) + 
                                pow(src.ys[i]-src.ys[i-1], 2)));
    }
    Normalize(res);
    return res;
}

double foley_dist(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2-x1, 2) + pow(y2-y1, 2));
}

double foley_alpha_hat(double x1, double y1, double x2, double y2, double x3, double y3) {
    double dx1 = x2-x1, dy1 = y2-y1, dx2 = x3-x2, dy2 = y3-y2;
    double dot = -dx1*dx2 - dy1*dy2;
    double l1 = foley_dist(x1, y1, x2, y2), l2 = foley_dist(x2, y2, x3, y3);
    double angle = acos(dot/(l1*l2));
    angle = std::min(4*atan(1)-angle, 2*atan(1));
    return angle;
}


std::vector<double> Parametrization_foley(const NodeArr &src) {
    using namespace std;
    assert(src.size > 3);
    vector<double> res(src.size, 0);
    for (int i=1; i<res.size(); ++i) {
        bool has_left = true, has_right = true;
        double angle_l, angle_r, length, length_l, length_r;
        if (i == 1) has_left = false;
        if (i == res.size()-1) has_right = false;
        length = foley_dist(src.xs[i-1], src.ys[i-1], src.xs[i], src.ys[i]);
        if (has_left) {
            angle_l = foley_alpha_hat(src.xs[i-2], src.ys[i-2], src.xs[i-1], src.ys[i-1], src.xs[i], src.ys[i]);
            length_l = foley_dist(src.xs[i-2], src.ys[i-2], src.xs[i-1], src.ys[i-1]);
        }
        if (has_right) {
            angle_r = foley_alpha_hat(src.xs[i-1], src.ys[i-1], src.xs[i], src.ys[i], src.xs[i+1], src.ys[i+1]);
            length_r = foley_dist(src.xs[i], src.ys[i], src.xs[i+1], src.ys[i+1]);
        }
        double base_len = length;
        double mul = 1;
        if (has_left) {
            mul += 1.5*angle_l*length_l/(length_l+length);
        }
        if (has_right) {
            mul += 1.5*angle_r*length_r/(length_r+length);
        }
        res[i] = res[i-1] + mul*base_len;
    }
    Normalize(res);
    return res;
}


// 多项式
void InterPolation_1(NodeArr &src, NodeArr &tar) {
    using namespace Eigen;
    const int num_pow = src.size;
    Matrix<double, Dynamic, Dynamic> A(num_pow, num_pow);
    VectorXd x(num_pow), b(num_pow);
    
    // build A and b
    for (int i=0; i<num_pow; ++i) {
        for (int j=0; j<num_pow; ++j) {
            A(i, j) = pow(src.xs[i], j);
        }
    }

    for (int i=0; i<num_pow; ++i) {
        b[i] = src.ys[i];
    }

    x = A.colPivHouseholderQr().solve(b);
    
    // calculate fit curve
    double st_region = src.xs[0], ed_region = src.xs[src.size-1];
    std::vector<double> resxs, resys;
    for (double xx=st_region; xx<ed_region; xx+=0.001) {
        double yy = 0;
        for (int i=0; i<num_pow; ++i) {
            yy += x(i)*std::pow(xx, i);
        }
        resxs.push_back(xx);
        resys.push_back(yy);
    }
    tar.xs.swap(resxs);
    tar.ys.swap(resys);
    tar.size = tar.xs.size();
    return;
}

// 高斯基函数
void InterPolation_2(NodeArr &src, NodeArr &tar) {
    using namespace Eigen;
    const int num_pow = src.size;
    Matrix<double, Dynamic, Dynamic> A(num_pow+1, num_pow+1);
    VectorXd x(num_pow+1), b(num_pow+1);
    double para1 = (src.xs.back()-src.xs.front())/(src.size-1);
    
    // build A and b
    for (int i=0; i<num_pow; ++i) {
        for (int j=0; j<num_pow; ++j) {
            A(i, j) = exp(-pow(src.xs[i]-src.xs[j], 2)/(2*pow(para1, 2)));
        }
        A(i, num_pow) = 1;
    }
    for (int j=0; j<num_pow; ++j) {
        A(num_pow, j) = 1;
    }

    for (int i=0; i<num_pow; ++i) {
        b[i] = src.ys[i];
    }
    b[num_pow] = 0;

    x = A.colPivHouseholderQr().solve(b);
    
    // calculate fit curve
    double st_region = src.xs[0], ed_region = src.xs[src.size-1];
    std::vector<double> resxs, resys;
    for (double xx=st_region; xx<ed_region; xx+=0.001) {
        double yy = 0;
        for (int i=0; i<num_pow; ++i) {
            yy += x(i)*exp(-pow(xx-src.xs[i], 2)/(2*para1*para1));
        }
        yy += x(num_pow);
        resxs.push_back(xx);
        resys.push_back(yy);
    }
    tar = NodeArr(resxs, resys);
}

void LeastSquare(NodeArr &src, NodeArr &tar, int order) {
    using namespace Eigen;
    int nodenum = src.size;
    Matrix<double, Dynamic, Dynamic> A(order, order);
    VectorXd x(order), b(order);
    
    // build A and b
    for (int i=0; i<order; ++i) {
        for (int j=0; j<order; ++j) {
            A(i, j) = 0;
            for (int k=0; k<nodenum; ++k) {
                A(i, j) += pow(src.xs[k], i+j);
            }
        }
    }

    for (int i=0; i<order; ++i) {
        b(i) = 0;
        for (int k=0; k<nodenum; ++k) {
            b(i) += pow(src.xs[k], i)*src.ys[k];
        }
    }

    x = A.colPivHouseholderQr().solve(b);
    
    // calculate fit curve
    double st_region = src.xs[0], ed_region = src.xs[src.size-1];
    std::vector<double> resxs, resys;
    for (double xx=st_region; xx<ed_region; xx+=0.001) {
        double yy = 0;
        for (int i=0; i<order; ++i) {
            yy += x(i)*pow(xx, i);
        }
        resxs.push_back(xx);
        resys.push_back(yy);
    }
    tar = NodeArr(resxs, resys);
}

void Ridge_Regression(NodeArr &src, NodeArr &tar, int order, double l) {
    using namespace Eigen;
    int nodenum = src.size;
    Matrix<double, Dynamic, Dynamic> A(order, order);
    VectorXd x(order), b(order);
    
    // build A and b
    for (int i=0; i<order; ++i) {
        for (int j=0; j<order; ++j) {
            A(i, j) = 0;
            for (int k=0; k<nodenum; ++k) {
                A(i, j) += pow(src.xs[k], i+j);
            }
            if (i == j) {
                A(i, j) += 2*l;
            }
        }
    }

    for (int i=0; i<order; ++i) {
        b(i) = 0;
        for (int k=0; k<nodenum; ++k) {
            b(i) += pow(src.xs[k], i)*src.ys[k];
        }
    }

    x = A.colPivHouseholderQr().solve(b);
    
    // calculate fit curve
    double st_region = src.xs[0], ed_region = src.xs[src.size-1];
    std::vector<double> resxs, resys;
    for (double xx=st_region; xx<ed_region; xx+=0.001) {
        double yy = 0;
        for (int i=0; i<order; ++i) {
            yy += x(i)*pow(xx, i);
        }
        resxs.push_back(xx);
        resys.push_back(yy);
    }
    tar = NodeArr(resxs, resys);
}