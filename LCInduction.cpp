#include<iostream>
#include<eigen3/Eigen/Dense>
#include<fstream>
#include <string>
#include<cmath> 

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Matriz{
	private:
	int n;
	public:
	MatrixXd K;
	Matriz(int size, double L, double C) : n(size), K(MatrixXd::Zero(n,n)){
		for (int i = 0; i < n; ++i) {
        K(i, i) = 2.0 * 1/(L*C); // Diagonal
        if (i > 0) {
            K(i, i - 1) = -1/(L*C); // Elemento inferior izquierdo
            K(i - 1, i) = -1/(L*C); // Elemento superior derecho
        }
    }
	}
	MatrixXd AutoVec(){
	Eigen::EigenSolver<MatrixXd> es(K);
	MatrixXd Vec = es.eigenvectors().real();
	return Vec; 
	}
	
	VectorXd AutoVal(){
	Eigen::EigenSolver<MatrixXd> es(K);
	VectorXd Val = es.eigenvalues().real();
  return Val;
	}

	void PrintM(){
	std::cout<<K<<std::endl;
	}
};

/* Tareas pendientes:
 * Factorización LU.
 * Cargas de los N circuitos.
 * Corrientes de los N circuitos.
 */
int main(){
	int size = 0;
	double M = 1;
	double L;
	double C;
	double t;
	double delta = 1000;
	MatrixXd K;
	std::cout<<"la cantidad de circuitos que desea acoplar"<<std::endl;
	std::cin>>size;
		std::cout<<"Introduzca la inductancia del sistema de circuitos" <<std::endl;
		std::cin>>L;
		std::cout<<"Introduzca la capacitancia del sistema de circuitos" <<std::endl;
		std::cin>>C;
	Matriz Mat(size,L,C);
	Mat.PrintM();
	MatrixXd Vecd = Mat.AutoVec();
	VectorXd Vald = Mat.AutoVal();
	for(int i=0; i<size; i++){	
		std::ofstream datos("datos"+ std::to_string(i) +".dat");
		 t=0.0;
			for(int j=0; j<delta ; j++){
				for(int k=0; k<size; k++){
					datos << Vecd.col(i)[k]*std::cos(Vald[i]*t) << "\t";
				}
			t = t + 0.01;
			datos << "\n";
			}

		datos.close();
	}






}


