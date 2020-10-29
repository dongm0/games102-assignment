#include <math.h>
#include <vector>
#include <eigen3/Eigen/Dense>

struct NodeArr {
    NodeArr() {}
    NodeArr(std::vector<double> x, std::vector<double> y) : xs(x), ys(y), size(x.size()) {}
    unsigned int size;
    std::vector<double> xs;
    std::vector<double> ys;
};

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
void InterPolation_2(NodeArr &src, NodeArr &tar, double para1 /*分母的delta*/) {
    using namespace Eigen;
    const int num_pow = src.size;
    Matrix<double, Dynamic, Dynamic> A(num_pow+1, num_pow+1);
    VectorXd x(num_pow+1), b(num_pow+1);
    
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