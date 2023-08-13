#pragma once

#include <cmath>
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <vector>

#define QR_MAX_ITERS 1000000
#define JACOBI_MAX_ITERS 1000000
#define POWER_MAX_ITERS 1000000
#define VAR_TASK_MAX_ITERS 1000000
#define GRADMIN_MAX_ITERS 1000000
#define GRAD_NULL_SEARCH_MAX_ITERS 1000000
#define GRAD_NULL_SEARCH_EL_ZERO 1e-6
#define GRADMIN_EPS 5e-8
#define ZERO_EPS 1e-12
#define HESSIAN_ZERO_EPS 1e-5

namespace mtrx{
	
const double PI = 3.1415926535;

template <typename T>
T GenRand(T min, T max){
	T res = min + static_cast<T>(rand())/( static_cast<T>(RAND_MAX/(max-min)));
	return res;
}
	
template <typename T>
class Matrix{
	public:
		Matrix();
		
		Matrix(int, int);
		
		template <typename Function>
		Matrix(int, int, Function);
		
		Matrix(const Matrix<T>&);
		
		~Matrix();
		
		void Set(int, int, T); // set element on i-th row, j-th column
		
		template <typename Function>
		void Fill(Function f);
		
		template <typename Function>
		void SymFill(Function f);
		
		void swap(Matrix<T>&) noexcept;
		
		void Transpose();
		
		int GetNRows() const;
		
		int GetNCols() const;
		
		Matrix<T> GetCol(int) const;
		
		Matrix<T>& operator=(const Matrix<T>&);
		
		Matrix<T> operator*(const Matrix<T>&) const;
		
		Matrix<T> operator+(const Matrix<T>&) const;
		
		Matrix<T> operator-(const Matrix<T>&) const;
		
		T Get(int, int) const;
		
		T Determinant() const;
		
		T Minor(int, int) const;
		
		const T* Data() const;
		
		template <typename Type>
		friend Matrix<Type> Inverse(const Matrix<Type>&);
		
		template <typename Type>
		friend Matrix<Type> GramSchmidt(const Matrix<Type>&);
		
	private:
		int n_rows_;
		int m_cols_;
		T *mat_ptr_ = nullptr;
		
		void Copy(const Matrix<T>&);
		Matrix<T> proj(const Matrix<T>&) const;
};

//A * x = B
Matrix<double> InverseSearch(const Matrix<double>&, const Matrix<double>&);

Matrix<double> Cramer(const Matrix<double>&, const Matrix<double>&);

template <typename T>
Matrix<T>::Matrix(){
	n_rows_ = 0;
	m_cols_ = 0;
	mat_ptr_ = nullptr;
}

template <typename T>
Matrix<T>::Matrix(int n, int m){
	assert(n > 0 && m > 0);
	n_rows_ = n;
	m_cols_ = m;
	mat_ptr_ = new T[n_rows_ * m_cols_];
}

template<typename T>
template<typename Function>
Matrix<T>::Matrix(int n, int m, Function f){
	assert(n > 0 && m > 0);
	n_rows_ = n;
	m_cols_ = m;
	mat_ptr_ = new T[n_rows_ * m_cols_];
	this->Fill(f);
}

template <typename T>
Matrix<T>::~Matrix(){
	delete []mat_ptr_;
	mat_ptr_ = nullptr;
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T> &other){
	this->Copy(other);
}

template <typename T>
void Matrix<T>::Set(int i, int j, T value){
	assert( i < n_rows_ && j < m_cols_ );
	mat_ptr_[m_cols_ * i + j] = value;
}

template <typename T>
T Matrix<T>::Get(int i, int j) const {
	assert( i < n_rows_ && j < m_cols_ );
	return mat_ptr_[m_cols_ * i + j];
}

template <typename T>
const T* Matrix<T>::Data() const {
	return mat_ptr_;
}

template <typename T>
Matrix<T> Matrix<T>::GetCol(int k) const {
	if(k >= m_cols_) throw std::out_of_range("Column number is out of range");
	return Matrix<T>(n_rows_, 1, [this, k](int i, int j){return this->mat_ptr_[this->m_cols_*i + k];});
}

template <typename T>
template <typename Function>
void Matrix<T>::Fill(Function f){
	int i, j;
	for(i = 0; i < n_rows_; ++i){
		for(j = 0; j < m_cols_; ++j){
			mat_ptr_[m_cols_ * i + j] = f(i, j);
		}
	}
}

template <typename T>
template <typename Function>
void Matrix<T>::SymFill(Function f){
	assert(n_rows_ == m_cols_);
	int i, j;
	for(i = 0; i < n_rows_; ++i){
		for(j = i; j < m_cols_; ++j){
			mat_ptr_[m_cols_ * i + j] = f(i, j);
			mat_ptr_[m_cols_ * j + i] = mat_ptr_[m_cols_ * i + j];
		}
	}
}

template <typename T>
void Matrix<T>::swap(Matrix<T> &other) noexcept{
	std::swap(this->mat_ptr_, other.mat_ptr_);
	std::swap(this->n_rows_, other.n_rows_);
	std::swap(this->m_cols_, other.m_cols_);
}

template <typename T>
void Matrix<T>::Copy(const Matrix<T> &other){
	Matrix<T> temp(other.n_rows_, other.m_cols_);
	temp.Fill( [&other](int i, int j){
		return other.mat_ptr_[other.m_cols_*i+j];
		} );
	this->swap(temp);
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &rhs){
	if( this == &rhs ) return *this;
	this->Copy(rhs);
	return *this;
}

template <typename T>
int Matrix<T>::GetNRows() const {
	return n_rows_;
}

template <typename T>
int Matrix<T>::GetNCols() const {
	return m_cols_;
}

template <typename T>
void Matrix<T>::Transpose(){
	Matrix<T> temp(m_cols_, n_rows_);
	temp.Fill( [this](int i, int j){
		return this->mat_ptr_[m_cols_*j+i];
		} );
	this->swap(temp);
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &rhs) const {
	assert(this->m_cols_ == rhs.n_rows_);
	Matrix<T> res(this->n_rows_, rhs.m_cols_);
	res.Fill( [this, &rhs](int i, int j){
			T s = 0;
			for(int k = 0; k < this->m_cols_; ++k){
				s += this->mat_ptr_[this->m_cols_*i+k] * rhs.mat_ptr_[rhs.m_cols_*k+j];
			}
			return s;
		} );
	return res;
}

template <typename T>
Matrix<T> operator*(T lhs, const Matrix<T> &rhs){
	Matrix<T> res(rhs.GetNRows(), rhs.GetNCols(), [lhs, &rhs](int i, int j){
			return rhs.Get(i, j)*(lhs);
		});
	return res;
}

template <typename T>
Matrix<T> operator*(const Matrix<T> &lhs, T rhs){
	Matrix<T> res(lhs.GetNRows(), lhs.GetNCols(), [&lhs, rhs](int i, int j){
			return lhs.Get(i, j)*(rhs);
		});
	return res;
}

template <typename T>
Matrix<T> operator*(int lhs, const Matrix<T> &rhs){
	Matrix<T> res(rhs.GetNRows(), rhs.GetNCols(), [lhs, &rhs](int i, int j){
			return rhs.Get(i, j)*static_cast<T>(lhs);
		});
	return res;
}

template <typename T>
Matrix<T> operator*(const Matrix<T> &lhs, int rhs){
	Matrix<T> res(lhs.GetNRows(), lhs.GetNCols(), [&lhs, rhs](int i, int j){
			return lhs.Get(i, j)*static_cast<T>(rhs);
		});
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &rhs) const{
	assert(this->n_rows_ == rhs.n_rows_ && this->m_cols_ == rhs.m_cols_);
	Matrix<T> res(this->n_rows_, this->m_cols_,
			[this, &rhs](int i, int j){return this->mat_ptr_[this->m_cols_*i + j] + rhs.mat_ptr_[rhs.m_cols_*i + j];});
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &rhs) const{
	assert(this->n_rows_ == rhs.n_rows_ && this->m_cols_ == rhs.m_cols_);
	Matrix<T> res(this->n_rows_, this->m_cols_,
			[this, &rhs](int i, int j){return this->mat_ptr_[this->m_cols_*i + j] - rhs.mat_ptr_[rhs.m_cols_*i + j];});
	return res;
}

template <typename T>
T Matrix<T>::Determinant() const {
	assert(n_rows_ == m_cols_);
	//printf("det-f called\n");
	if(n_rows_ == 1) return mat_ptr_[0];
	if(n_rows_ == 2) return mat_ptr_[0]*mat_ptr_[3] - mat_ptr_[1]*mat_ptr_[2];
	T res = 0;
	int k;
	for(k = 0; k < n_rows_; ++k){
		Matrix<T> matM(n_rows_ - 1, m_cols_ - 1);
		matM.Fill( [this, k](int i, int j){
				return this->mat_ptr_[this->m_cols_*(i + 1) + (j + static_cast<int>(j >= k))];
			} );
		/*matM.Print(std::cout);
		std::cout << ( (k%2 == 0) ? (1) : (-1) ) << std::endl;*/
		res += mat_ptr_[k]*matM.Determinant()*( (k%2 == 0) ? (1) : (-1) );
	}
	return res;
}

template <typename T>
T Matrix<T>::Minor(int k, int l) const {
	Matrix matM(n_rows_ - 1, m_cols_ - 1);
	matM.Fill( [this, k, l](int i, int j){
			return this->mat_ptr_[this->m_cols_*(i + static_cast<int>(i >= k)) + (j + static_cast<int>(j >= l))];
		} );
	return matM.Determinant();
}

template <typename T> //Usable for T = float, double
Matrix<T> Inverse(const Matrix<T> &source){
	assert(source.n_rows_ == source.m_cols_);
	T det = source.Determinant();
	if(std::abs(det) < ZERO_EPS) throw std::exception();
	Matrix<T> res(source.n_rows_, source.m_cols_);
	res.Fill( [&source, det](int i, int j){
			return source.Minor(j, i)*(( (i+j)%2 == 0 ) ? (1) : (-1))/det;
		} );
	return res;
}

void Print(const Matrix<double>&);

void Print(const Matrix<int>&);

template <typename T>
Matrix<T> Matrix<T>::proj(const Matrix<T> &base) const {
	assert(base.m_cols_ == 1 && this->m_cols_ == 1 && this->n_rows_ == base.n_rows_);
	T scalar_ab = 0, scalar_bb = 0, mult;
	for(int i = 0; i < base.n_rows_; ++i){
		scalar_bb += base.mat_ptr_[i]*base.mat_ptr_[i];
		scalar_ab += this->mat_ptr_[i]*base.mat_ptr_[i];
	}
	mult = scalar_ab / scalar_bb;
	return Matrix<T>(this->n_rows_, 1, [&base, mult](int i, int j){return base.mat_ptr_[i]*mult;});
}

template <typename T>
Matrix<T> GramSchmidt(const Matrix<T> &source){
	Matrix<T> *res_vecs = new Matrix<T>[source.m_cols_]; 
	Matrix<T> *source_vecs = new Matrix<T>[source.m_cols_];
	T *norms = new T[source.m_cols_]; 
	Matrix<T> res(source.n_rows_, source.m_cols_);
	int j, k;
	for(j = 0; j < source.m_cols_; ++j){
		source_vecs[j] = source.GetCol(j);
	}
	for(j = 0; j < source.m_cols_; ++j){
		res_vecs[j] = source_vecs[j];
		for(k = 0; k < j; ++k){
			res_vecs[j] = res_vecs[j] - source_vecs[j].proj(res_vecs[k]);
		}
	}
	for(j = 0; j < source.m_cols_; ++j){
		norms[j] = 0;
		for(k = 0; k < source.n_rows_; ++k){
			norms[j] += res_vecs[j].mat_ptr_[k]*res_vecs[j].mat_ptr_[k];
		}
		norms[j] = sqrt(norms[j]);
	}
	res.Fill( [&res_vecs, &norms](int i, int j){
			return res_vecs[j].mat_ptr_[i]/norms[j];
		} );
	delete []res_vecs;
	delete []source_vecs;
	delete []norms;
	return res;
}

template <typename T>
struct QR_decomposition{
	Matrix<T> Q;
	Matrix<T> R;
};

struct Eigen{
	Matrix<double> eig_vec;
	Matrix<double> eig_val; // <- diagonalized source matrix
	int iters = 0;
};

struct SingleEigen{
	Matrix<double> eig_vec;
	double eig_val;
};

template <typename T>
Matrix<T> Transpose(const Matrix<T> &source){
	Matrix<T> res = source;
	res.Transpose();
	return res;
}

template <typename T>
QR_decomposition<T> QR_decompose(const Matrix<T> &source){
	Matrix<T> matQ = GramSchmidt(source);
	Matrix<T> matR = Transpose(matQ)*source;
	return {matQ, matR};
}

Eigen Eigen_QR(const Matrix<double> &source, int max_iters = QR_MAX_ITERS);

Eigen Eigen_Jacobi(const Matrix<double> &source, int max_iters = JACOBI_MAX_ITERS);

SingleEigen Eigen_Power(const Matrix<double>&, int,
                        const Matrix<double>&);

template <typename T>
Matrix<T> SLTInverse(const Matrix<T> &source){ //Builds inverse matrix for square lower triangular matrix
	int n = source.GetNRows();
	assert(n == source.GetNCols());
	int i, j, k;
    T *res = new double[n*n];
    T *mi = new double[n*n];
    const T *a = source.Data();
    
    for(i = 0; i < n*n ; ++i){
        res[i] = a[i];
        mi[i] = 0;
    }
    for(i = 0; i < n; ++i){
        mi[n*i + i] = 1;
    }
    for(i = 0; i < n; ++i){
        for(j = 0; j < n; ++j){
            mi[n*i + j] = mi[n*i + j]/res[n*i + i];
        }
        for(j = 0; j < n; ++j){
            res[n*i + j] = res[n*i + j]/res[n*i + i];
        }
        for(j = 0; j < n; ++j){
            for(k = i + 1; k < n; ++k){
                mi[n*k + j] = mi[n*k + j] - mi[n*i + j]*res[n*k + i];
            }
        }
    }
    Matrix<T> matChi(n, n, [n, mi](int i, int j){return mi[n*i + j];});
    delete []mi;
    delete []res;
    return matChi;
}

Eigen VarTask(const Matrix<double> &matH, const Matrix<double> &matS,
              Eigen (*df)(const Matrix<double>&, int) = &(mtrx::Eigen_Jacobi), int max_iters = VAR_TASK_MAX_ITERS);

//Container must have:
//1. Copy constructor
//2. Method size()
//3. Operator =
//4. Operator []
//5. Constructor Container(std::size_t)
//6. Iterators, begin(), end()
//7. push_back()
//8. Operator +, -
 
template<typename Container, typename Function> 
double GradMin(Function f, Container &v, double eps = GRADMIN_EPS, double step = 0.5){
	double gradF, *dfdr, h = 1e-6, s, box, u, un;
	int iter = 0, i, n = static_cast<int>(v.size());
	dfdr = new double[n];
	Container v_res = v;
	s = 0;
	u = f(v);
	for(i = 0; i < n; ++i){
		box = v_res[i];
		v_res[i] = v[i] + h;
		dfdr[i] = (f(v_res) - u)/h;
		s += dfdr[i]*dfdr[i];
		v_res[i] = box;
	}
	gradF = sqrt(s);
	for(auto x : v){
		std::cout << x << " ";
	}
	std::cout << "U = " << u << " gradU = " << gradF << std::endl;
	if(gradF < eps) return u;
	for(iter = 0; iter < GRADMIN_MAX_ITERS;){
		for(i = 0; i < n; ++i){
			v_res[i] = v[i] - step*dfdr[i]/gradF;
		}
		un = f(v_res);
		if(un < u){
			u = un;
			v = v_res;
			s = 0;
			for(i = 0; i < n; ++i){
				box = v_res[i];
				v_res[i] = v[i] + h;
				dfdr[i] = (f(v_res) - u)/h;
				s += dfdr[i]*dfdr[i];
				v_res[i] = box;
			}
			gradF = sqrt(s);
			++iter;
			if(iter%1 == 0){
			for(auto x : v){
				std::cout << x << " ";
			}
			std::cout << "U = " << u << " gradU = " << gradF << " step = " << step << " | iter " << iter << std::endl;
			}
			if(gradF < eps) break;
		} else {
			step = 0.5 * step;
			if(step < 1e-7) break;
		}
	}
	delete []dfdr;
	return u;
}

template<typename Container, typename Function> 
double Derivative(Function f, const Container &vars, int index, double h = 1e-4){
	Container vars_ph = vars;
	Container vars_mh = vars;
	vars_ph[index] = vars[index] + h;
	vars_mh[index] = vars[index] - h;
	return (f(vars_ph) - f(vars_mh))/(2*h);
}

template<typename Container, typename Function> 
double SecondDerivative(Function f, const Container &vars, int i, int j, double h = 1e-4){
	Container vars_ph = vars;
	Container vars_mh = vars;
	vars_ph[i] = vars[i] + h;
	vars_mh[i] = vars[i] - h;
	return ( Derivative(f, vars_ph, j, h) - Derivative(f, vars_mh, j, h) )/(2*h);
}

template<typename Container, typename Function> 
Matrix<double> Hessian(Function f, const Container &vars, double h = 1e-4){
	double value;
	Matrix<double> res(static_cast<int>(vars.size()), static_cast<int>(vars.size()));
	//res.SymFill([f, &vars, h](int i, int j){return SecondDerivative(f, vars, i, j, h);});
	for(int i = 0; i < res.GetNRows(); ++i){
		for(int j = i; j < res.GetNCols(); ++j){
			value = SecondDerivative(f, vars, i, j, h);
			res.Set(i, j, value);
			res.Set(j, i, value);
		}
	}
	return res;
}

template<typename Container, typename Function> 
Container Gradient(Function f, const Container &vars, double h = 1e-7){
	Container res(vars.size());
	Container vars_ph = vars;
	double u = f(vars);
	int i, n = static_cast<int>(vars.size());
	for(i = 0; i < n; ++i){
		vars_ph[i] = vars[i] + h;
		res[i] = (f(vars_ph) - u)/h;
		vars_ph[i] = vars[i];
	}
	return res;
}

template<typename Container>
double Norm(const Container &v){
	double res = 0;
	for(std::size_t i = 0; i < v.size(); ++i){
		res += v[i]*v[i];
	}
	return sqrt(res);
}
			
template<typename Container>
double Norm2(const Container &v){
	double res = 0;
	for(std::size_t i = 0; i < v.size(); ++i){
		res += v[i]*v[i];
	}
	return res;
}

template<typename Container, typename Function> 
bool GradNullSearch(Function f, Container &x){
	int n = static_cast<int>(x.size());
	Matrix<double> matA(n, n);
	Matrix<double> matB(n, 1);
	Matrix<double> Xwave(n, 1);
	std::vector<double> derivatives(n);
	Container reserve_copy_x = x;
	double gradU = 0;
	matA = Hessian(f, x);
	int minus_count = 0;
	Matrix<double> diag_matA = (Eigen_Jacobi(matA)).eig_val;
	std::cout << "Hessians:" << std::endl;
	Print(matA);
	Print(diag_matA);
	for(int i = 0; i < n; ++i){
		if( diag_matA.Get(i, i) < 0 ) ++minus_count;
		if( std::abs(diag_matA.Get(i, i)) < HESSIAN_ZERO_EPS ) {
			std::cout << "Zero(s) in diag. hessian" << std::endl;
			return false;
		}
	}
	if(minus_count != 1) {
		std::cout << "Negative value is not even" << std::endl;
		return false;
	}
	for(int iter = 0; iter < GRAD_NULL_SEARCH_MAX_ITERS; ++iter){
		//matA = Hessian(f, x);
		gradU = 0;
		for(int i = 0; i < n; ++i){
			derivatives[i] = Derivative(f, x, i);
			matB.Set(i, 0, (-1)*derivatives[i]);
			gradU += derivatives[i]*derivatives[i];
		}
		gradU = sqrt(gradU);
		if(gradU < GRADMIN_EPS) break;
		//Xwave = InverseSearch(matA, matB);
		try{
			Xwave = Inverse(matA) * matB;
			//Xwave = Cramer(matA, matB);
		} catch(const std::exception &e){
				std::cout << "Unsolvable hessian" << std::endl;
				x = reserve_copy_x;
				return false;
			}
		for(int i = 0; i < n; ++i){
			std::cout << x[i] << " ";
			x[i] = Xwave.Get(i, 0) + x[i];
			//std::cout << x[i] << " ";
		}
		std::cout << "| gradU = " << gradU << std::endl;
		matA = Hessian(f, x);
		minus_count = 0;
		diag_matA = (Eigen_Jacobi(matA)).eig_val;
		std::cout << "Hessians:" << std::endl;
		Print(matA);
		Print(diag_matA);
		for(int i = 0; i < n; ++i){
			if( diag_matA.Get(i, i) < 0 ) ++minus_count;
			if( std::abs(diag_matA.Get(i, i)) < HESSIAN_ZERO_EPS ){
				x = reserve_copy_x;
				std::cout << "Zero(s) in diag. hessian" << std::endl;
				return false;
			}
		}
		if(minus_count != 1) {
			x = reserve_copy_x;
			std::cout << "Negative value is not even" << std::endl;
			return false;
		}
	}
	for(int i = 0; i < n; ++i){
		std::cout << x[i] << " ";
	}
	std::cout << "| gradU = " << gradU << std::endl;
	return true;
}


template<typename Container>
typename std::vector<Container>::iterator AddCenter(std::vector<Container> &line){
	std::size_t max_index = 0;
	double max_length = 0, length = 0;
	Container res(line[0].size());
	max_length = Norm(line[1] - line[0]);
	max_index = 0;
	for(std::size_t i = 1; i < line.size() - 1; ++i){
		length = Norm(line[i+1] - line[i]);
		if(length > max_length){
			max_length = length;
			max_index = i;
		}
	}
	for(std::size_t i = 0; i < res.size(); ++i){
		res[i] = (line[max_index][i] + line[max_index + 1][i])*0.5;
	}
	return line.insert(line.begin() + max_index + 1, res);
}

template<typename Container, typename Function> 
Container TSSearch(Function f, const Container &source1, const Container &source2){
	assert(source1.size() == source2.size());
	int iter;
	double step = 0.5, abs_gradF, u, un;
	typename std::vector<Container>::iterator center;
	Container center_res;
	Container gradF;
	std::vector<Container> way = {source1, source2};
	for(iter = 0; iter < GRADMIN_MAX_ITERS; ++iter){
		center = AddCenter(way);
		std::cout << "center: ( ";
		for(auto it = center->begin(); it != center->end(); ++it){
			std::cout << *it;
			if(it + 1 != center->end()) std::cout << " ";
		}
		std::cout << " )" << std::endl;
		if(GradNullSearch(f, *center)) {
			std::cout << "ts: ( ";
			for(auto it = center->begin(); it != center->end(); ++it){
				std::cout << *it;
				if(it + 1 != center->end()) std::cout << " ";
			}
			std::cout << " )" << std::endl;
			std::cout << "way: ";
			for(std::size_t i = 0; i < way.size(); ++i){
				std::cout << "( ";
				for(auto it = way[i].begin(); it != way[i].end(); ++it){
					std::cout << *it;
					if(it + 1 != way[i].end()) std::cout << " ";
				}
				std::cout << " )" << std::endl;
			}
			return *center;
		}
		gradF = Gradient(f, *center);
		abs_gradF = Norm(gradF);
		step = 0.5;
		std::cout << "Loc. minimizing..." << std::endl;
		for(;;){
			center_res = (*center);
			u = f(*center);
			for(std::size_t j = 0; j < (*center).size(); ++j){
				(*center)[j] -= step*gradF[j]/abs_gradF;
			}
			un = f(*center);
			if(un > u){
				*center = center_res;
				step = step*0.5;
			}
			else{
				u = un;
			}
			if(step < 1e-7) break;
		}
		std::cout << "Loc. minimizing finished" << std::endl;
		if(GradNullSearch(f, *center)) {
			std::cout << "ts: ( ";
			for(auto it = center->begin(); it != center->end(); ++it){
				std::cout << *it;
				if(it + 1 != center->end()) std::cout << " ";
			}
			std::cout << " )" << std::endl;
			std::cout << "way: ";
			for(std::size_t i = 0; i < way.size(); ++i){
				std::cout << "( ";
				for(auto it = way[i].begin(); it != way[i].end(); ++it){
					std::cout << *it;
					if(it + 1 != way[i].end()) std::cout << " ";
				}
				std::cout << " )" << std::endl;
			}
			return *center;
		}
		std::vector<std::size_t> nums = {0, way.size()-1};
		std::vector<double> dists;
		double dist = 0, min_dist = 0;
		std::size_t min_index;
		typename std::vector<Container> new_way;
		new_way.push_back(source1);
		for(std::size_t i = 0; i < way.size() - 1; ++i){
			min_dist = Norm(way[i] - way[i+1]);
			min_index = i + 1;
			for(std::size_t j = 0; j < way.size() - 1; ++j){
				if(find(nums.begin(), nums.end(), j) != nums.end()) continue;
				dist = Norm(way[i] - way[j]);
				if(dist < min_dist){
					min_dist = dist;
					min_index = j;
				}
			}
			nums.push_back(min_index);
			new_way.push_back(way[min_index]);
		}
		for(std::size_t i = 0; i < way.size(); ++i){
			std::cout << "( ";
			for(auto it = way[i].begin(); it != way[i].end(); ++it){
				std::cout << *it;
				if(it + 1 != way[i].end()) std::cout << " ";
			}
			std::cout << " )" << std::endl;
		}
	}
	return *center;
}
}

const std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs);

const std::vector<double> operator-(const std::vector<double> &lhs, const std::vector<double> &rhs);

