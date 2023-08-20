#include "matrix.h"

namespace mtrx{
//A * x = B
Matrix<double> InverseSearch(const Matrix<double> &matA, const Matrix<double> &matB){
    return Inverse(matA) * matB;
}

Matrix<double> Cramer(const Matrix<double> &matA, const Matrix<double> &matB){
    int i, j;
    double detA = matA.Determinant();
    std::cout << "detA = " << detA << std::endl;
    if(std::abs(detA) < ZERO_EPS) throw std::exception();
    std::vector<double> dets(matA.GetNCols());
    Matrix<double> matR = matA;
    for(j = 0; j < matA.GetNCols(); ++j){
        for(i = 0; i < matA.GetNRows(); ++i){
            matR.Set(i, j, matB.Get(i, 0));
        }
        dets[j] = matR.Determinant();
        for(i = 0; i < matA.GetNRows(); ++i){
            matR.Set(i, j, matA.Get(i, j));
        }
    }
    Matrix<double> matX(matA.GetNRows(), 1);
    for(i = 0; i < matX.GetNRows(); ++i){
        matX.Set(i, 0, dets[i]/detA);
    }
    return matX;
}

void Print(const Matrix<double> &source){
    printf("( ");
    int i, j;
    double value;
    for(i = 0; i < source.GetNRows(); ++i){
        if(i != 0) printf("  ");
        for(j = 0; j < source.GetNCols(); ++j){
            value = source.Get(i, j);
            printf( "%12.6lf ",
                   (std::abs(value) < ZERO_EPS) ? std::abs(value) : value );
        }
        if(i != source.GetNRows()-1) printf("\n");
    }
    printf(")\n");
}

void Print(const Matrix<int> &source){
    printf("( ");
    int i, j;
    int value;
    for(i = 0; i < source.GetNRows(); ++i){
        if(i != 0) printf("  ");
        for(j = 0; j < source.GetNCols(); ++j){
            value = source.Get(i, j);
            printf( "%5d ", value );
        }
        if(i != source.GetNRows()-1) printf("\n");
    }
    printf(")\n");
}

Eigen Eigen_QR(const Matrix<double> &source, int max_iters){
    int n = source.GetNRows(), m = source.GetNCols();
    assert(n == m);
    int i, j, iter;
    QR_decomposition<double> qr;
    Eigen eig;
    double max;
    eig.eig_val = source;
    eig.eig_vec = Matrix<double>(n, m, [](int p, int q){return static_cast<double>(p == q);});
    for(iter = 0; iter < max_iters; ++iter){
        qr = QR_decompose(eig.eig_val);
        eig.eig_vec = eig.eig_vec*qr.Q;
        eig.eig_val = (qr.R)*(qr.Q);
        max = eig.eig_val.Get(1, 0);
        for(j = 0; j < m; ++j){
            for(i = j + 1; i < n; ++i){
                if(max < std::abs(eig.eig_val.Get(i, j)))
                    max = std::abs(eig.eig_val.Get(i, j));
            }
        }
        if(max < ZERO_EPS) break;
    }
    eig.iters = iter;
    return eig;
}

Eigen Eigen_Jacobi(const Matrix<double> &source, int max_iters){
    int n = source.GetNRows(), m = source.GetNCols();
    assert(n == m);
    if(n == 1){
        return {Matrix<double>(1, 1, [](int i, int j){return 1;}), Matrix<double>(1, 1, [source](int i, int j){return source.Get(0, 0);})};
    }
    int i, j, iter, i_max, j_max;
    Eigen res;
    double max, theta, sin_th, cos_th;
    res.eig_val = source;
    res.eig_vec = Matrix<double>(n, m, [](int p, int q){return static_cast<double>(p == q);});
    for(iter = 0; iter < max_iters; ++iter){
        Matrix<double> matR(n, n, [](int i, int j){return static_cast<double>(i == j);});
        max = std::abs(res.eig_val.Get(0, 1));
        i_max = 0;
        j_max = 1;
        for(i = 0; i < n - 1; ++i){
            for(j = i + 1; j < m; ++j){
                if(std::abs(res.eig_val.Get(i, j)) > max){
                    i_max = i;
                    j_max = j;
                    max = std::abs(res.eig_val.Get(i, j));
                }
            }
        }
        if(max < ZERO_EPS) break;
        theta = 0.5 * atan( 2*res.eig_val.Get(i_max, j_max)/
                           ( res.eig_val.Get(i_max, i_max) - res.eig_val.Get(j_max, j_max) ) );
        sin_th = sin(theta);
        cos_th = cos(theta);
        matR.Set(i_max, i_max, cos_th);
        matR.Set(j_max, j_max, cos_th);
        matR.Set(i_max, j_max, -sin_th);
        matR.Set(j_max, i_max, sin_th);
        res.eig_val = Transpose(matR) * res.eig_val * matR;
        res.eig_vec = res.eig_vec * matR;
    }
    if(iter == max_iters) std::cout << "MATRIX DIAG. ERROR!!!" << std::endl;
    std::vector<int> nums(n);
    for(i = 0; i < n; ++i){
        nums[i] = i;
    }
    std::sort(nums.begin(), nums.end(), [&res](int p, int q){
        return res.eig_val.Get(p, p) < res.eig_val.Get(q, q);
    });
    Matrix<double> vec_box(n, n);
    Matrix<double> val_box(n, n, [](int i, int j){return 0;});
    for(j = 0; j < n; ++j){
        val_box.Set(j, j, res.eig_val.Get(nums[j], nums[j]));
        for(i = 0; i < n; ++i){
            vec_box.Set(i, j, res.eig_vec.Get(i, nums[j]));
        }
    }
    return {vec_box, val_box, iter};
}

SingleEigen Eigen_Power(const Matrix<double> &source, int max_iters = POWER_MAX_ITERS,
                        const Matrix<double> &start = Matrix<double>()){
    int n = source.GetNRows(), m = source.GetNCols();
    assert(n == m);
    Matrix<double> res(n, 1);
    Matrix<double> prev;
    bool stop = false;
    double norm = 0;
    int i, iter;
    if(start.Data() == nullptr){
        double start_value = 1./sqrt(n);
        res.Fill( [start_value](int i, int j){return /*static_cast<int>(i==0)*/start_value;} );
    } else {
        assert(start.GetNRows() == n && start.GetNCols() == 1);
        res = start;
    }
    for(iter = 0; iter < max_iters; ++iter){
        prev = res;
        res = source * res;
        norm = 0;
        for(i = 0; i < n; ++i){
            norm += res.Get(i, 0)*res.Get(i, 0);
        }
        norm = sqrt(norm);
        res = res * (1./norm);
        stop = true;
        for(i = 0; i < n; ++i){
            if(std::abs(std::abs(prev.Get(i, 0)) - std::abs(res.Get(i, 0))) > ZERO_EPS){
                stop = false;
                break;
            }
        }
        if(stop) break;
    }
    return {res, (Transpose(res)*source*res).Get(0, 0)};
}

Matrix<double> Cholesky(const Matrix<double> &source){
    int n = source.GetNRows(), m = source.GetNCols();
    assert(source.GetNRows() == source.GetNCols());
    double *ch = new double[n*n];
    const double *a = source.Data();
    double sij, sii;
    int i, j, p;
    for(i = 0; i < n*n ; ++i){
        ch[i] = 0;
    }
    ch[0] = sqrt(a[0]);
    for(i = 1; i < n; ++i){
        ch[m*i]=a[m*i]/ch[0];
    }
    ch[m+1] = sqrt( a[m+1] - ch[m]*ch[m] );
    for(i = 2; i < n; ++i){
        for(j = 1; j <= i; ++j){
            if(i != j){
                sij = 0;
                for(p = 0; p < j; ++p){
                    sij += ch[m*j + p]*ch[m*i + p];
                }
                ch[m*i + j] = (a[m*i + j] - sij)/(ch[m*j + j]);
            }
            else{
                sii = 0;
                for(p = 0; p < i; ++p){
                    sii += ch[m*i + p]*ch[m*i + p];
                }
                ch[m*i + i]=sqrt(a[m*i + i] - sii);
            }
        }
    }
    Matrix<double> res(n, n, [ch, n](int i, int j){return ch[n*i + j];});
    a = nullptr;
    delete []ch;
    ch = nullptr;
    return res;
}

Eigen VarTask(const Matrix<double> &matH, const Matrix<double> &matS,
              Eigen (*df)(const Matrix<double>&, int), int max_iters){ // HC = SCe
    Matrix<double> matChi = SLTInverse(Cholesky(matS));
    Matrix<double> matChiT = Transpose(matChi);
    Matrix<double> matB = matChi * matH * matChiT;
    Eigen pre_res = df(matB, max_iters);
    return {matChiT*pre_res.eig_vec, pre_res.eig_val};
}

}

const std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs){
    std::vector<double> res(lhs.size());
    for(std::size_t i = 0; i < lhs.size(); ++i){
        res[i] = lhs[i] + rhs[i];
    }
    return res;
}

const std::vector<double> operator-(const std::vector<double> &lhs, const std::vector<double> &rhs){
    std::vector<double> res(lhs.size());
    for(std::size_t i = 0; i < lhs.size(); ++i){
        res[i] = lhs[i] - rhs[i];
    }
    return res;
}
