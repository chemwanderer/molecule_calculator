#pragma once

#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include <map>
#include "matrix.h"

#define SCF_MAX_ITERS 1000
#define SCF_EPS 1e-9

namespace mlc{
	
struct BasisFunction{
	int nuc_num;
	double a;
	double dx;
	double dy;
	double dz;
};

class IntegralBox{
	public:
		double Get(int i, int j, int p, int q) const {
			if(i > j) std::swap(i, j);
			if(p > q) std::swap(p, q);
			return integrals.at({{i, j}, {p, q}});
		}
		void Set(int i, int j, int p, int q, double value){
			if(i > j) std::swap(i, j);
			if(p > q) std::swap(p, q);
			integrals[{{i, j}, {p, q}}] = value;
		}
		void Clear(){
			integrals.clear();
		}
	private:
		std::map<std::pair<std::pair<int, int>, std::pair<int, int>>, double> integrals;
};

template <typename Container>
class Molecule{
	public:
		Molecule(int, const Container&, const std::vector<int>&);
		Molecule(const Molecule<Container>&);
		
		void SetCoord(int, int, double);
		
		void SetCoord(int, double);
		
		void SetCoords(Container&);
		
		double GetCoord(int, int) const;
		
		double GetCoord(int) const;
		
		int GetNNucs() const;
		
		void SetNAEls(int);
		
		void SetNBEls(int);
		
		int GetNEls() const;
		
		std::vector<int> GetCharges() const;
		
		Container GetCoords() const;
		
		void Print() const;
		
		double Dist(int, int) const;
		
		void AddBasisFunction(int, double, double, double, double);
		
		void AddS(int, double);
		
		void AddPx(int, double, double);
		void AddPy(int, double, double);
		void AddPz(int, double, double);
		
		double Vnn() const;
		
		double Ee_1e() const;
		
		double Ee_RHF(bool);
		double Ee_UHF(bool);
		
		template <typename C>
		friend double Ee_1e(Molecule<C>&);
		
		void Optimize(const std::string&, bool);
		
		bool TS_Refine(const std::string&, bool);
		
		std::pair<Molecule<Container>, Molecule<Container>> TS_Relax(const std::string&, bool);
		
		mtrx::Matrix<double> Hessian(const std::string&, bool);
		
		void Normalize_Coords();
		
	private:
		int N_els_ = 0;
		int N_els_a_ = 0;
		int N_els_b_ = 0;
		Container coords_;
		std::vector<int> charges_;
		//std::vector<Container> basis_;
		std::vector<BasisFunction> basis_ = {};
		IntegralBox integrals_;
		mtrx::Matrix<double> matC_;
		mtrx::Matrix<double> matCa_;
		mtrx::Matrix<double> matCb_;
		bool hasC_ = false;
		bool hasCab_ = false;
		
		int n_bas_ = 0;
		
		std::pair<mtrx::Matrix<double>, mtrx::Matrix<double>> BuildSH() const;
		mtrx::Matrix<double> BuildP_RHF(mtrx::Matrix<double>) const;
		mtrx::Matrix<double> BuildP_UHF(mtrx::Matrix<double>, char) const;
		
		std::pair<int, double> Choose(int) const;
		
		double Ee_calc(const std::string&, bool);
};

template<typename Container>
mtrx::Matrix<double> Molecule<Container>::Hessian(const std::string &method, bool use_prev_C){
	//this->Normalize_Coords(); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Container norm_coords;
	Container res;
	for(unsigned int i = 0; i < coords_.size(); ++i){
		if(i <= 2 || i == 4 || i == 5 || i == 8) continue;
		else norm_coords.push_back(coords_[i]);
	}
	for(std::size_t i = 0; i < norm_coords.size(); ++i){
	}
	return mtrx::Hessian([this, &res, &method, use_prev_C](const Container &v){
			double energy;
			int k = 0;
			res = this->coords_;
			for(std::size_t i = 0; i < coords_.size(); ++i){
				if(i <= 2 || i == 4 || i == 5 || i == 8) coords_[i] = 0;
				else {
					coords_[i] = v.at(k);
					++k;
				}
			}
			energy = this->Ee_calc(method, use_prev_C);
			this->coords_ = res;
			return energy;
		}, norm_coords);
}

template<typename Container>
double Molecule<Container>::Ee_calc(const std::string &method, bool use_prev_C){
	using namespace std::string_literals;
	if(method == "RHF"s) return this->Ee_RHF(use_prev_C);
	if(method == "UHF"s) return this->Ee_UHF(use_prev_C);
	return this->Ee_UHF(false);
}

template<typename Container>
void Molecule<Container>::SetNAEls(int N_els_a){
	N_els_a_ = N_els_a;
}

template<typename Container>
void Molecule<Container>::SetNBEls(int N_els_b){
	N_els_b_ = N_els_b;
}

template<typename Container>
void Molecule<Container>::SetCoord(int i, double value){
	this->coords_[i] = value;
}

template<typename Container>
void Molecule<Container>::SetCoords(Container &c){
	this->coords_ = c;
}

template<typename Container>
std::vector<int> Molecule<Container>::GetCharges() const {
	return this->charges_;
}

template<typename Container>
Container Molecule<Container>::GetCoords() const {
	return this->coords_;
}

template<typename Container>
int Molecule<Container>::GetNEls() const {
	return this->N_els_;
}

template<typename Container>
void Molecule<Container>::Normalize_Coords(){
	double dx = -coords_[0], dy = -coords_[1], dz = -coords_[2];
	double fi, x0, y0, z0;
	for(int i = 0; i < static_cast<int>(charges_.size()); ++i){
		coords_[3*i + 0] += dx;
		coords_[3*i + 1] += dy;
		coords_[3*i + 2] += dz;
	}
	if(charges_.size() == 1) return;
	fi = atan(coords_[5]/coords_[3]);
	for(int i = 0; i < static_cast<int>(charges_.size()); ++i){
		x0 = coords_[3*i + 0];
		z0 = coords_[3*i + 2];
		coords_[3*i + 0] = x0*cos(fi) + z0*sin(fi);
		coords_[3*i + 2] = -x0*sin(fi) + z0*cos(fi);
	}
	fi = atan(-coords_[4]/coords_[3]);
	for(int i = 0; i < static_cast<int>(charges_.size()); ++i){
		x0 = coords_[3*i + 0];
		y0 = coords_[3*i + 1];
		coords_[3*i + 0] = x0*cos(fi) - y0*sin(fi);
		coords_[3*i + 1] = x0*sin(fi) + y0*cos(fi);
	}
	if(charges_.size() == 2) return;
	fi = atan(-coords_[8]/coords_[7]);
	for(int i = 0; i < static_cast<int>(charges_.size()); ++i){
		y0 = coords_[3*i + 1];
		z0 = coords_[3*i + 2];
		coords_[3*i + 1] = y0*cos(fi) - z0*sin(fi);
		coords_[3*i + 2] = y0*sin(fi) + z0*cos(fi);
	}
}

template<typename Container>
mtrx::Matrix<double> Molecule<Container>::BuildP_RHF(mtrx::Matrix<double> matC) const {
	mtrx::Matrix<double> matP(n_bas_, n_bas_);
	int i, j, k;
	double p;
	for(i = 0; i < n_bas_; ++i){
		for(j = 0; j < n_bas_; ++j){
			p = 0;
			for(k = 0; k < N_els_/2; ++k){
				p += matC.Get(i, k)*matC.Get(j, k);
			}
			matP.Set(i, j, p);
		}
	}
	return matP;
}

template<typename Container>
mtrx::Matrix<double> Molecule<Container>::BuildP_UHF(mtrx::Matrix<double> matC, char g) const {
	mtrx::Matrix<double> matP(n_bas_, n_bas_);
	int i, j, k, N_els_g = (g == 'a') ? N_els_a_ : N_els_b_;
	double p;
	for(i = 0; i < n_bas_; ++i){
		for(j = 0; j < n_bas_; ++j){
			p = 0;
			for(k = 0; k < N_els_g; ++k){
				p += matC.Get(i, k)*matC.Get(j, k);
			}
			matP.Set(i, j, p);
		}
	}
	return matP;
}
/*
template<typename Container>
void Molecule<Container>::Optimize_RHF(){
	Container res;
	mtrx::GradMin([this, &res](const Container &v){
			//Molecule<Container> res = (*this);
			double energy;
			res = this->coords_;
			this->coords_ = v;
			energy = this->Ee_RHF();
			this->coords_ = res;
			return energy;
		}, this->coords_);
}*/

template<typename Container>
bool Molecule<Container>::TS_Refine(const std::string &method, bool use_prev_C){
	//this->Normalize_Coords(); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Container norm_coords;
	Container res;
	for(unsigned int i = 0; i < coords_.size(); ++i){
		if(i <= 2 || i == 4 || i == 5 || i == 8) continue;
		else norm_coords.push_back(coords_[i]);
	}
	for(std::size_t i = 0; i < norm_coords.size(); ++i){
	}
	bool result = mtrx::GradNullSearch([this, &res, &method, use_prev_C](const Container &v){
			double energy;
			int k = 0;
			res = this->coords_;
			for(std::size_t i = 0; i < coords_.size(); ++i){
				if(i <= 2 || i == 4 || i == 5 || i == 8) coords_[i] = 0;
				else {
					coords_[i] = v.at(k);
					++k;
				}
			}
			energy = this->Ee_calc(method, use_prev_C);
			this->coords_ = res;
			return energy;
		}, norm_coords);
		int k = 0;
		for(std::size_t i = 0; i < coords_.size(); ++i){
			if(i <= 2 || i == 4 || i == 5 || i == 8) coords_[i] = 0;
			else {
				coords_[i] = norm_coords[k];
				++k;
			}
		}
		//this->Print();
		return result;
}

template<typename Container>
void Molecule<Container>::Optimize(const std::string &method, bool use_prev_C){
	using namespace std::string_literals;
	//this->Normalize_Coords(); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Container norm_coords;
	Container res;
	
	for(std::size_t i = 0; i < coords_.size(); ++i){
		if(i <= 2 || i == 4 || i == 5 || i == 8) continue;
		else norm_coords.push_back(coords_[i]);
	}
	mtrx::GradMin([this, &res, &method, use_prev_C](const Container &v){
			//Molecule<Container> res = (*this);
			double energy;
			int k = 0;
			res = this->coords_;
			//this->coords_ = v;
			for(std::size_t i = 0; i < coords_.size(); ++i){
				if(i <= 2 || i == 4 || i == 5 || i == 8) coords_[i] = 0;
				else {
					coords_[i] = v[k];
					++k;
				}
			}
			energy = this->Ee_calc(method, use_prev_C);
			
			this->coords_ = res;
			return energy;
		}, norm_coords);
		int k = 0;
		for(std::size_t i = 0; i < coords_.size(); ++i){
			if(i <= 2 || i == 4 || i == 5 || i == 8) coords_[i] = 0;
			else {
				coords_[i] = norm_coords[k];
				++k;
			}
		}
}
/*
template<typename Container>
void Molecule<Container>::Optimize_UHF(){
	//this->Normalize_Coords(); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Container norm_coords;
	Container res;
	for(std::size_t i = 0; i < coords_.size(); ++i){
		if(i <= 2 || i == 4 || i == 5 || i == 8) continue;
		else norm_coords.push_back(coords_[i]);
	}
	mtrx::GradMin([this, &res](const Container &v){
			//Molecule<Container> res = (*this);
			double energy;
			int k = 0;
			res = this->coords_;
			//this->coords_ = v;
			for(std::size_t i = 0; i < coords_.size(); ++i){
				if(i <= 2 || i == 4 || i == 5 || i == 8) coords_[i] = 0;
				else {
					coords_[i] = v[k];
					++k;
				}
			}
			energy = this->Ee_UHF(false);
			this->coords_ = res;
			return energy;
		}, norm_coords);
		int k = 0;
		for(std::size_t i = 0; i < coords_.size(); ++i){
			if(i <= 2 || i == 4 || i == 5 || i == 8) coords_[i] = 0;
			else {
				coords_[i] = norm_coords[k];
				++k;
			}
		}
}*/

template <typename Container>
Molecule<Container> TSSearch(Molecule<Container> &m1, Molecule<Container> &m2){
	/*m1.Normalize_Coords();
	m2.Normalize_Coords();*/
	Container norm_coords1, norm_coords2, TS_coords;
	Container res(m1.GetNNucs()*3);
	for(int i = 0; i < m1.GetNNucs()*3; ++i){
		if(i <= 2 || i == 4 || i == 5 || i == 8) continue;
		else norm_coords1.push_back(m1.GetCoord(i));
	}
	for(int i = 0; i < m2.GetNNucs()*3; ++i){
		if(i <= 2 || i == 4 || i == 5 || i == 8) continue;
		else norm_coords2.push_back(m2.GetCoord(i));
	}
	TS_coords = mtrx::TSSearch([&m1, &res](const Container &v){
			double energy;
			int k = 0;
			res = m1.GetCoords();
			for(int i = 0; i < m1.GetNNucs()*3; ++i){
				if(i <= 2 || i == 4 || i == 5 || i == 8)  m1.SetCoord(i, 0);
				else {
					m1.SetCoord(i, v[k]);
					++k;
				}
			}
			energy = m1.Ee_UHF(false);
			m1.SetCoords(res);
			return energy;
		}, norm_coords1, norm_coords2);
		int k = 0;
		for(int i = 0; i < m1.GetNNucs()*3; ++i){
			if(i <= 2 || i == 4 || i == 5 || i == 8) res[i] = 0;
			else {
				res[i] = TS_coords[k];
				++k;
			}
		}
		return Molecule<Container>(m1.GetNEls(), res, m1.GetCharges());
}

template <typename Container>
std::pair<Molecule<Container>, Molecule<Container>> Molecule<Container>::TS_Relax(const std::string &method, bool use_prev_C){
	mtrx::Matrix<double> matHes = this->Hessian(method, use_prev_C);
	std::cout << "TS relaxation; start Hessian:" << std::endl;
	mtrx::Print(matHes);
	auto eig = mtrx::Eigen_Jacobi(matHes);
	std::cout << "start Hessian (diag):" << std::endl;
	mtrx::Print(eig.eig_val);
	/*Molecule<Container> m1 = (*this);
	Molecule<Container> m2 = (*this);*/
	Molecule<Container> m1(this->N_els_, this->coords_, this->charges_);
	Molecule<Container> m2(this->N_els_, this->coords_, this->charges_);
	m1.basis_ = this->basis_;
	m2.basis_ = this->basis_;
	m1.n_bas_ = this->n_bas_;
	m2.n_bas_ = this->n_bas_;
	m1.N_els_a_ = this->N_els_a_;
	m2.N_els_a_ = this->N_els_a_;
	m1.N_els_b_ = this->N_els_b_;
	m2.N_els_b_ = this->N_els_b_;
	int k = 0;
	for(int i = 0; i < m1.GetNNucs()*3; ++i){
		if(i <= 2 || i == 4 || i == 5 || i == 8){
			 //res[i] = 0;
			 m1.SetCoord(i, 0);
			 m2.SetCoord(i, 0);
		 }
		else {
			m1.SetCoord(i, m1.coords_[i] + 0.5*eig.eig_vec.Get(k, 0));
			m2.SetCoord(i, m2.coords_[i] - 0.5*eig.eig_vec.Get(k, 0));
			/*res[i] = TS_coords[k];*/
			++k;
		}
	}
	/*m1.Optimize(method, use_prev_C);
	m1.Print();*/
	m2.Optimize(method, use_prev_C);
	m2.Print();
	m1.Optimize(method, use_prev_C);
	m1.Print();
	return std::make_pair(m1, m2);
	//return {Molecule<Container>(m1.N_els_, m1.coords_, m1.charges_), Molecule<Container>(m2.N_els_, m2.coords_, m2.charges_)};
}

template <typename Container>
Molecule<Container>::Molecule(int N_els, const Container &coords, const std::vector<int> &charges){
	N_els_ = N_els;
	coords_ = coords;
	charges_ = charges;
}

template <typename Container>
Molecule<Container>::Molecule(const Molecule<Container> &other){
	this->N_els_ = other.N_els_;
	this->N_els_a_ = other.N_els_a_;
	this->N_els_b_ = other.N_els_b_;
	this->coords_ = other.coords_;
	this->charges_ = other.charges_;
	this->basis_ = other.basis_;
	/*this->matC_ = other.matC_;
	this->matCa_ = other.matCa_;
	this->matCb_ = other.matCb_;*/
	this->hasC_ = false;
	this->hasCab_ = false;
	this->n_bas_ = other.n_bas_;
}

template <typename Container>
double Molecule<Container>::Dist(int i, int j) const {
	Container res(3);
	res[0] = coords_[3*i + 0] - coords_[3*j + 0];
	res[1] = coords_[3*i + 1] - coords_[3*j + 1];
	res[2] = coords_[3*i + 2] - coords_[3*j + 2];
	return mtrx::Norm(res);
}

template <typename Container>
void Molecule<Container>::SetCoord(int i_nuc, int i_coord, double value){ 
	coords_[3*i_nuc + i_coord] = value;
}

template <typename Container>
double Molecule<Container>::GetCoord(int i_nuc, int i_coord) const {
	return coords_[3*i_nuc + i_coord];
}

template <typename Container>
double Molecule<Container>::GetCoord(int i_coord) const {
	return coords_[i_coord];
}

template <typename Container>
int Molecule<Container>::GetNNucs() const {
	return static_cast<int>(charges_.size());
}

template <typename Container>
void Molecule<Container>::AddBasisFunction(int i, double a, double dx, double dy, double dz){
	basis_.push_back({i, a, dx, dy, dz});
	++n_bas_;
}

template <typename Container>
void Molecule<Container>::AddS(int i, double a){
	basis_.push_back({i, a, 0, 0, 0});
	++n_bas_;
}

template <typename Container>
void Molecule<Container>::AddPx(int i, double a, double dx){
	basis_.push_back({i, a, dx, 0, 0});
	++n_bas_;
	basis_.push_back({i, a, -dx, 0, 0});
	++n_bas_;
}

template <typename Container>
void Molecule<Container>::AddPy(int i, double a, double dy){
	basis_.push_back({i, a, 0, dy, 0});
	++n_bas_;
	basis_.push_back({i, a, 0, -dy, 0});
	++n_bas_;
}

template <typename Container>
void Molecule<Container>::AddPz(int i, double a, double dz){
	basis_.push_back({i, a, 0, 0, dz});
	++n_bas_;
	basis_.push_back({i, a, 0, 0, -dz});
	++n_bas_;
}
/*
template <typename Container>
void Molecule<Container>::AddS(int i, double a){
	basis_[i].push_back(a);
	++n_bas_;
}*/

template <typename Container>
std::pair<mtrx::Matrix<double>, mtrx::Matrix<double>> Molecule<Container>::BuildSH() const {
	/*mtrx::Matrix<double> matS(n_bas_, n_bas_);
	mtrx::Matrix<double> matH(n_bas_, n_bas_);*/
	//std::pair<int, double> bas_f1, bas_f2;
	Container rp(3);
	Container r1(3);
	Container r2(3);
	int nuc1, nuc2;
	double a1, a2, su, AB2, CP2;
	int i, j, k;
	mtrx::Matrix<double> matS(n_bas_, n_bas_);
	mtrx::Matrix<double> matH(n_bas_, n_bas_);
	
	for(i = 0; i < n_bas_; ++i){
		nuc1 = basis_[i].nuc_num;
		a1 = basis_[i].a; 
		
		r1[0] = coords_[3*nuc1+0] + basis_[i].dx;
		r1[1] = coords_[3*nuc1+1] + basis_[i].dy;
		r1[2] = coords_[3*nuc1+2] + basis_[i].dz;
		
		for(j = i; j < n_bas_; ++j){
			nuc2 = basis_[j].nuc_num;
			a2 = basis_[j].a;
			
			r2[0] = coords_[3*nuc2+0] + basis_[j].dx;
			r2[1] = coords_[3*nuc2+1] + basis_[j].dy;
			r2[2] = coords_[3*nuc2+2] + basis_[j].dz;
			
			su = a1 + a2;
			
			AB2 = mtrx::Norm2(r1 - r2);
			
			for(k = 0; k < 3; ++k){
				rp[k] = ( a1*r1[k] + a2*r2[k] )/su;
			}
			
			matS.Set(i, j, pow(mtrx::PI/su,1.5)*exp(-a1*a2*AB2/su)*pow(2*a1/mtrx::PI,0.75)*pow(2*a2/mtrx::PI,0.75));
			matS.Set(j, i, matS.Get(i, j));
			
			matH.Set(i, j, -(a1*a2)/su*matS.Get(i,j)*( 2*(a1*a2)/su*AB2-3 )); //KES addition
			for(k = 0; k < static_cast<int>(charges_.size()); ++k){
				CP2 = (coords_[3*k + 0] - rp[0])*(coords_[3*k + 0] - rp[0]) 
				+ (coords_[3*k + 1] - rp[1])*(coords_[3*k + 1] - rp[1])
				+ (coords_[3*k + 2] - rp[2])*(coords_[3*k + 2] - rp[2]);
				
				matH.Set( i, j, matH.Get(i, j) - ( (CP2 < ZERO_EPS) ? charges_[k]*4/su*sqrt(2./mtrx::PI)*pow(a1*a2, 0.75)*exp(-a1*a2*AB2/su)
					: charges_[k]*2/su*sqrt(2/(su*CP2))*pow(a1*a2, 0.75)*exp(-a1*a2*AB2/su)*erf(sqrt(su*CP2)) ) );
			}
			matH.Set(j, i, matH.Get(i, j));
		}
	}
	return {matS, matH};
}
/*
template <typename Container>
std::pair<int, double> Molecule<Container>::Choose(int num_to_choose) const{ //OPTIMIZATION REQUIRED!
	int num = 0;
	for(std::size_t i = 0; i < charges_.size(); ++i){
		for(std::size_t j = 0; j < basis_[i].size(); ++j){
			if(num == num_to_choose) return {i, basis_[i][j]};
			++num;
		}
	}
	return {static_cast<int>(charges_.size() - 1), basis_[n_bas_ - 1][n_bas_ - 1]};
}*/

template <typename Container>
double Molecule<Container>::Vnn(void) const {
    int a, b, nnuc = static_cast<int>(charges_.size());
    double v = 0, d;
    for(a = 0; a < nnuc; ++a){
        for(b = a + 1; b < nnuc; ++b){
			d = Dist(a, b);
            if(d > ZERO_EPS && charges_[a] != 0 && charges_[b] != 0) v += (charges_[a])*(charges_[b])/(Dist(a, b));
        }
    }
    return v;
}

template <typename Container>
double Molecule<Container>::Ee_1e() const{
	auto SH = BuildSH();
	auto res = mtrx::VarTask(SH.second, SH.first);
	return res.eig_val.Get(0, 0) + Vnn();
}

template <typename Container>
double Molecule<Container>::Ee_RHF(bool use_prev_C){
	int i, j, p, q, k;
	auto SH = BuildSH();
	mtrx::Matrix<double> matH = SH.second;
	mtrx::Matrix<double> matS = SH.first;
	/*auto single_e_res = mtrx::VarTask(SH.second, SH.first);
	mtrx::Matrix<double> matC = single_e_res.eig_vec;*/
	mtrx::Matrix<double> matC;
	if(use_prev_C && hasC_){
		matC = this->matC_;
	}
	else{
		auto single_e_res = mtrx::VarTask(SH.second, SH.first);
		matC = single_e_res.eig_vec;
	}
	mtrx::Matrix<double> matP = BuildP_RHF(matC);
	mtrx::Matrix<double> matP1(n_bas_, n_bas_);
	mtrx::Matrix<double> matF(n_bas_, n_bas_);
	
	double a1, a2, a3, a4, su1, su2, AB2, CD2, PQ2, d, next_integral, value_to_add, max_dP;
	Container rp(3);
	Container rq(3);
	Container r1(3);
	Container r2(3);
	Container r3(3);
	Container r4(3);
	int nuc1, nuc2, nuc3, nuc4;
	
	for(i = 0; i < n_bas_; ++i){
		nuc1 = basis_[i].nuc_num;
		a1 = basis_[i].a;
		
		r1[0] = coords_[3*nuc1+0] + basis_[i].dx;
		r1[1] = coords_[3*nuc1+1] + basis_[i].dy;
		r1[2] = coords_[3*nuc1+2] + basis_[i].dz;
		
		for(j = i; j < n_bas_; ++j){
			nuc2 = basis_[j].nuc_num;
			a2 = basis_[j].a;
			
			su1 = a1 + a2;
			
			r2[0] = coords_[3*nuc2+0] + basis_[j].dx;
			r2[1] = coords_[3*nuc2+1] + basis_[j].dy;
			r2[2] = coords_[3*nuc2+2] + basis_[j].dz;
			
			AB2 = mtrx::Norm2(r1 - r2);
			
			for(k = 0; k < 3; ++k){
				rp[k] = ( a1*r1[k] + a2*r2[k] )/su1;
			}
			
			for(p = 0; p < n_bas_; ++p){
				nuc3 = basis_[p].nuc_num;
				a3 = basis_[p].a;
				
				r3[0] = coords_[3*nuc3+0] + basis_[p].dx;
				r3[1] = coords_[3*nuc3+1] + basis_[p].dy;
				r3[2] = coords_[3*nuc3+2] + basis_[p].dz;
				
				for(q = p; q < n_bas_; ++q){
					nuc4 = basis_[q].nuc_num;
					a4 = basis_[q].a;
					
					su2 = a3 + a4;
					
					r4[0] = coords_[3*nuc4+0] + basis_[q].dx;
					r4[1] = coords_[3*nuc4+1] + basis_[q].dy;
					r4[2] = coords_[3*nuc4+2] + basis_[q].dz;
					
					CD2 = mtrx::Norm2(r3 - r4);
					
					for(k = 0; k < 3; ++k){
						rq[k] = ( a3*r3[k] + a4*r4[k] )/su2;
					}
					
					PQ2 = mtrx::Norm2(rp - rq);
					
					d = 1./(4*su1) + 1./(4*su2);
					
					/*if(PQ2==0) return 16/(su1*su2)/sqrt(PI*(su1+su2))*pow(a1*a2*a3*a4,0.75)*exp(-a1*a2*AB2/su1-a3*a4*CD2/su2);
					return 8/(su1*su2*sqrt((su1+su2)*(PQ2/(4*d))))*pow(a1*a2*a3*a4,0.75)*exp(-a1*a2*AB2/su1-a3*a4*CD2/su2)*erf(sqrt(PQ2/(4*d)));*/
					
					if(std::abs(PQ2) < ZERO_EPS) next_integral = 16/(su1*su2)/sqrt(mtrx::PI*(su1+su2))*pow(a1*a2*a3*a4,0.75)*exp(-a1*a2*AB2/su1-a3*a4*CD2/su2);
					else next_integral = 8/(su1*su2*sqrt((su1+su2)*(PQ2/(4*d))))*pow(a1*a2*a3*a4,0.75)*exp(-a1*a2*AB2/su1-a3*a4*CD2/su2)*erf(sqrt(PQ2/(4*d)));
					
					integrals_.Set(i, j, p, q, next_integral);
				}
			}
		}
	}
	
	for(k = 0; k < SCF_MAX_ITERS; ++k){
		for(i = 0; i < n_bas_; ++i){
			for(j = i; j < n_bas_; ++j){
				value_to_add = 0;
				for(p = 0; p < n_bas_; ++p){
					for(q = 0; q < n_bas_; ++q){
						value_to_add += matP.Get(p, q)*(2*integrals_.Get(i,j,p,q) - integrals_.Get(i,q,p,j));
					}
				}
				matF.Set(i, j, matH.Get(i, j) + value_to_add);
				matF.Set(j, i, matF.Get(i, j));
			}
		}
		auto res = mtrx::VarTask(matF, matS);
		matC = res.eig_vec;
		matP1 = BuildP_RHF(matC);
		max_dP = 0;
		for(i = 0; i < n_bas_; ++i){
			for(j = 0; j < n_bas_; ++j){
				if(std::abs(matP.Get(i, j) - matP1.Get(i, j)) > max_dP){
					max_dP = std::abs(matP.Get(i, j) - matP1.Get(i, j));
				}
			}
		}
		matP = matP1;
		//std::cout << max_dP << ", iter " << k << std::endl;
		if(max_dP < SCF_EPS) break;
	}
	if(k == SCF_MAX_ITERS) /*std::cout << "NO CONVERGENCE!!!" << std::endl*/return 0;
	else {
		this->matC_ = matC;
		this->hasC_ = true;
	}
	
	double Ee = 0;
	for(i = 0; i < n_bas_; ++i){
		for(j = 0; j < n_bas_; ++j){
			Ee += 2*matP.Get(i, j)*matH.Get(i, j);
		}
	}
	for(i = 0; i < n_bas_; ++i){
		for(j = 0; j < n_bas_; ++j){
			for(p = 0; p < n_bas_; ++p){
				for(q = 0; q < n_bas_; ++q){
					Ee += matP.Get(i, j)*matP.Get(p, q)*(2*integrals_.Get(i,j,p,q) - integrals_.Get(i,q,p,j));
				}
			}
		}
	}
	Ee += Vnn();
	integrals_.Clear();
	//std::cout << "Ee(RHF) = " << Ee << ", SCF convergence: " << (k != SCF_MAX_ITERS) << std::endl;
	return Ee;
}

template <typename Container>
double Molecule<Container>::Ee_UHF(bool use_prev_C){
	int i, j, p, q, k;
	//std::cout << "here" << std::endl;
	auto SH = BuildSH();
	//std::cout << "there" << std::endl;
	mtrx::Matrix<double> &matH = SH.second;
	mtrx::Matrix<double> &matS = SH.first;
	//mtrx::Print(matS);
	/*auto single_e_res = mtrx::VarTask(SH.second, SH.first);
	mtrx::Matrix<double> matC = single_e_res.eig_vec;*/
	mtrx::Matrix<double> matCa;
	mtrx::Matrix<double> matCb;
	if(use_prev_C && hasCab_){
		matCa = this->matCa_;
		matCb = this->matCb_;
	}
	else{
		//std::cout << "here" << std::endl;
		auto single_e_res = mtrx::VarTask(matH, matS);
		//std::cout << "here" << std::endl;
		matCa = single_e_res.eig_vec;
		matCb = matCa;
		mtrx::Matrix<double> matH1(n_bas_, n_bas_);
		for(i = 0; i < n_bas_; ++i){
			for(j = i; j < n_bas_; ++j){
				//matCb.Set(i, j, matCb.Get(i, j) + 0.1*matCb.Get(i, j));
				matH1.Set(i, j, matH.Get(i, j) + mtrx::GenRand(0., 0.1)*matH.Get(i, j));
				matH1.Set(j, i, matH1.Get(i, j));
				/*mtrx::Print(matCa);
				mtrx::Print(matCb);*/
			}
		}
		matCb = (mtrx::VarTask(matH1, SH.first)).eig_vec;
	}
	
	mtrx::Matrix<double> matPa = BuildP_UHF(matCa, 'a');
	mtrx::Matrix<double> matPb = BuildP_UHF(matCb, 'b');
	mtrx::Matrix<double> matPa1(n_bas_, n_bas_);
	mtrx::Matrix<double> matPb1(n_bas_, n_bas_);
	mtrx::Matrix<double> matFa(n_bas_, n_bas_);
	mtrx::Matrix<double> matFb(n_bas_, n_bas_);
	
	double a1, a2, a3, a4, su1, su2, AB2, CD2, PQ2, d, next_integral, integral_k, integral_j, value_to_add_a, value_to_add_b, max_dP;
	Container rp(3);
	Container rq(3);
	Container r1(3);
	Container r2(3);
	Container r3(3);
	Container r4(3);
	int nuc1, nuc2, nuc3, nuc4;
	
	for(i = 0; i < n_bas_; ++i){
		nuc1 = basis_[i].nuc_num;
		a1 = basis_[i].a;
		
		r1[0] = coords_[3*nuc1+0] + basis_[i].dx;
		r1[1] = coords_[3*nuc1+1] + basis_[i].dy;
		r1[2] = coords_[3*nuc1+2] + basis_[i].dz;
		
		for(j = i; j < n_bas_; ++j){
			nuc2 = basis_[j].nuc_num;
			a2 = basis_[j].a;
			
			su1 = a1 + a2;
			
			r2[0] = coords_[3*nuc2+0] + basis_[j].dx;
			r2[1] = coords_[3*nuc2+1] + basis_[j].dy;
			r2[2] = coords_[3*nuc2+2] + basis_[j].dz;
			
			AB2 = mtrx::Norm2(r1 - r2);
			
			for(k = 0; k < 3; ++k){
				rp[k] = ( a1*r1[k] + a2*r2[k] )/su1;
			}
			
			for(p = 0; p < n_bas_; ++p){
				nuc3 = basis_[p].nuc_num;
				a3 = basis_[p].a;
				
				r3[0] = coords_[3*nuc3+0] + basis_[p].dx;
				r3[1] = coords_[3*nuc3+1] + basis_[p].dy;
				r3[2] = coords_[3*nuc3+2] + basis_[p].dz;
				
				for(q = p; q < n_bas_; ++q){
					nuc4 = basis_[q].nuc_num;
					a4 = basis_[q].a;
					
					su2 = a3 + a4;
					
					r4[0] = coords_[3*nuc4+0] + basis_[q].dx;
					r4[1] = coords_[3*nuc4+1] + basis_[q].dy;
					r4[2] = coords_[3*nuc4+2] + basis_[q].dz;
					
					CD2 = mtrx::Norm2(r3 - r4);
					
					for(k = 0; k < 3; ++k){
						rq[k] = ( a3*r3[k] + a4*r4[k] )/su2;
					}
					
					PQ2 = mtrx::Norm2(rp - rq);
					
					d = 1./(4*su1) + 1./(4*su2);
					
					/*if(PQ2==0) return 16/(su1*su2)/sqrt(PI*(su1+su2))*pow(a1*a2*a3*a4,0.75)*exp(-a1*a2*AB2/su1-a3*a4*CD2/su2);
					return 8/(su1*su2*sqrt((su1+su2)*(PQ2/(4*d))))*pow(a1*a2*a3*a4,0.75)*exp(-a1*a2*AB2/su1-a3*a4*CD2/su2)*erf(sqrt(PQ2/(4*d)));*/
					
					if(std::abs(PQ2) < ZERO_EPS) next_integral = 16/(su1*su2)/sqrt(mtrx::PI*(su1+su2))*pow(a1*a2*a3*a4,0.75)*exp(-a1*a2*AB2/su1-a3*a4*CD2/su2);
					else next_integral = 8/(su1*su2*sqrt((su1+su2)*(PQ2/(4*d))))*pow(a1*a2*a3*a4,0.75)*exp(-a1*a2*AB2/su1-a3*a4*CD2/su2)*erf(sqrt(PQ2/(4*d)));
					
					integrals_.Set(i, j, p, q, next_integral);
				}
			}
		}
	}
	for(k = 0; k < SCF_MAX_ITERS; ++k){
		for(i = 0; i < n_bas_; ++i){
			for(j = i; j < n_bas_; ++j){
				value_to_add_a = 0;
				value_to_add_b = 0;
				for(p = 0; p < n_bas_; ++p){
					for(q = 0; q < n_bas_; ++q){
						integral_j = integrals_.Get(i,j,p,q);
						integral_k = integrals_.Get(i,q,p,j);
						//value_to_add_a += matPa.Get(p, q)*(integral_j - integral_k) + matPb.Get(p, q)*integral_j;
						//value_to_add_b += matPb.Get(p, q)*(integral_j - integral_k) + matPa.Get(p, q)*integral_j;
						value_to_add_a += (matPa.Get(p, q) + matPb.Get(p, q))*integral_j - matPa.Get(p, q)*integral_k;
						value_to_add_b += (matPa.Get(p, q) + matPb.Get(p, q))*integral_j - matPb.Get(p, q)*integral_k;
					}
				}
				matFa.Set(i, j, matH.Get(i, j) + value_to_add_a);
				matFa.Set(j, i, matFa.Get(i, j));
				matFb.Set(i, j, matH.Get(i, j) + value_to_add_b);
				matFb.Set(j, i, matFb.Get(i, j));
			}
		}
		auto res_a = mtrx::VarTask(matFa, matS);
		auto res_b = mtrx::VarTask(matFb, matS);
		matCa = res_a.eig_vec;
		matCb = res_b.eig_vec;
		matPa1 = BuildP_UHF(matCa, 'a');
		matPb1 = BuildP_UHF(matCb, 'b');
		max_dP = 0;
		for(i = 0; i < n_bas_; ++i){
			for(j = 0; j < n_bas_; ++j){
				if(std::abs(matPa.Get(i, j) - matPa1.Get(i, j)) > max_dP){
					max_dP = std::abs(matPa.Get(i, j) - matPa1.Get(i, j));
				}
			}
		}
		for(i = 0; i < n_bas_; ++i){
			for(j = 0; j < n_bas_; ++j){
				if(std::abs(matPb.Get(i, j) - matPb1.Get(i, j)) > max_dP){
					max_dP = std::abs(matPb.Get(i, j) - matPb1.Get(i, j));
				}
			}
		}
		matPa = matPa1;
		matPb = matPb1;
		//std::cout << max_dP << ", iter " << k + 1 << std::endl;
		if(max_dP < SCF_EPS) break;
	}
	if(k == SCF_MAX_ITERS) /*std::cout << "NO CONVERGENCE!!!" << std::endl*/return 0;
	else {
		this->matCa_ = matCa;
		this->matCb_ = matCb;
		this->hasCab_ = true;
	}
	
	double Ee = 0;
	for(i = 0; i < n_bas_; ++i){
		for(j = 0; j < n_bas_; ++j){
			Ee += (matPa.Get(i, j) + matPb.Get(i, j))*matH.Get(i, j);
		}
	}
	for(i = 0; i < n_bas_; ++i){
		for(j = 0; j < n_bas_; ++j){
			for(p = 0; p < n_bas_; ++p){
				for(q = 0; q < n_bas_; ++q){
					//Ee += matP.Get(i, j)*matP.Get(p, q)*(2*integrals_.Get(i,j,p,q) - integrals_.Get(i,q,p,j));
					//Ee += (matPa.Get(i, j)*matPa.Get(p, q) + matPb.Get(i, j)*matPb.Get(p, q))*(integrals_.Get(i,j,p,q) - integrals_.Get(i,q,p,j));
					Ee += 0.5*(matPa.Get(i, j) + matPb.Get(i, j))*(matPa.Get(p, q) + matPb.Get(p, q))*integrals_.Get(i, j, p, q);
				}
			}
		}
	}
	for(i = 0; i < n_bas_; ++i){
		for(j = 0; j < n_bas_; ++j){
			for(p = 0; p < n_bas_; ++p){
				for(q = 0; q < n_bas_; ++q){
					//Ee += matP.Get(i, j)*matP.Get(p, q)*(2*integrals_.Get(i,j,p,q) - integrals_.Get(i,q,p,j));
					//Ee += (matPa.Get(i, j)*matPa.Get(p, q) + matPb.Get(i, j)*matPb.Get(p, q))*(integrals_.Get(i,j,p,q) - integrals_.Get(i,q,p,j));
					Ee -= 0.5*(matPa.Get(i, j)*matPa.Get(p, q) + matPb.Get(i, j)*matPb.Get(p, q))*integrals_.Get(i, q, p, j);
				}
			}
		}
	}
	Ee += Vnn();
	integrals_.Clear();
	/*mtrx::Print(matH);
	mtrx::Print(matFa);
	mtrx::Print(matFb);*/
	//std::cout << "Ee(UHF) = " << Ee << ", SCF convergence: " << (k != SCF_MAX_ITERS) << std::endl;
	return Ee;
}

template <typename Container>
void Molecule<Container>::Print() const {
	for(std::size_t i = 0; i < charges_.size(); ++i){
		std::cout << charges_[i] << " " << coords_[3*i + 0] << " " << coords_[3*i + 1] << " " << coords_[3*i + 2] << std::endl;
	}
}

template <typename Container>
double Ee_1e(Molecule<Container> &m){
	auto SH = m.BuildSH();
	auto res = mtrx::VarTask(SH.second, SH.first);
	std::cout << res.eig_val.Get(0, 0) + m.Vnn() << std::endl;
	return res.eig_val.Get(0, 0) + m.Vnn();
}

}



