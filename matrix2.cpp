#include <iostream>
#include <eigen3/Eigen/Dense>
#include <chrono>

using namespace Eigen;
using namespace std;
using namespace std::chrono;

class CircuitLC {
private:
    MatrixXd A;
    VectorXd b;
    int n;

public:
    CircuitLC(int size) : n(size), A(MatrixXd::Random(size, size)), b(VectorXd::Random(size)) {
        // Para asegurar que el sistema tenga una solución única, haremos la matriz A diagonalmente dominante
        for (int i = 0; i < n; i++) {
            A(i, i) += n; // Aumentamos el valor de la diagonal principal
        }
    }

    VectorXd gaussSeidel(double tol = 1e-6, int maxIter = 10000) {
        VectorXd x = VectorXd::Zero(n); // Solución inicial en ceros
        VectorXd x_old = x;

        for (int iter = 0; iter < maxIter; ++iter) {
            for (int i = 0; i < n; ++i) {
                double sum = 0;
                for (int j = 0; j < n; ++j) {
                    if (j != i) {
                        sum += A(i, j) * x[j];
                    }
                }
                x[i] = (b[i] - sum) / A(i, i);
            }
            if ((x - x_old).norm() < tol) {
                break;
            }
            x_old = x;
        }

        return x;
    }

    VectorXd jacobi(double tol = 1e-6, int maxIter = 10000) {
        VectorXd x = VectorXd::Zero(n); // Solución inicial en ceros
        VectorXd x_new = x;

        for (int iter = 0; iter < maxIter; ++iter) {
            for (int i = 0; i < n; ++i) {
                double sum = 0;
                for (int j = 0; j < n; ++j) {
                    if (j != i) {
                        sum += A(i, j) * x[j];
                    }
                }
                x_new[i] = (b[i] - sum) / A(i, i);
            }
            if ((x_new - x).norm() < tol) {
                break;
            }
            x = x_new;
        }

        return x;
    }

    // Método LU utilizando Eigen
    VectorXd solveLU() {
        PartialPivLU<MatrixXd> lu_decomp(A); // Descomposición LU con pivoteo parcial
        VectorXd x = lu_decomp.solve(b);     // Resolver el sistema LUx = b
        return x;
    }

    void measureExecutionTime() {
        auto start = high_resolution_clock::now();
        VectorXd solutionGS = gaussSeidel();
        auto end = high_resolution_clock::now();
        auto durationGS = duration_cast<microseconds>(end - start);
        cout << "Gauss-Seidel solution: " << solutionGS.transpose() << endl;
        cout << "Gauss-Seidel execution time: " << durationGS.count() << " microseconds" << endl;

        start = high_resolution_clock::now();
        VectorXd solutionJacobi = jacobi();
        end = high_resolution_clock::now();
        auto durationJacobi = duration_cast<microseconds>(end - start);
        cout << "Jacobi solution: " << solutionJacobi.transpose() << endl;
        cout << "Jacobi execution time: " << durationJacobi.count() << " microseconds" << endl;

        // Medir el tiempo de ejecución del método LU
        start = high_resolution_clock::now();
        VectorXd solutionLU = solveLU();
        end = high_resolution_clock::now();
        auto durationLU = duration_cast<microseconds>(end - start);
        cout << "LU solution: " << solutionLU.transpose() << endl;
        cout << "LU execution time: " << durationLU.count() << " microseconds" << endl;
    }

    void printSystem() {
        cout << "Matrix A:\n" << A << endl;
        cout << "Vector b:\n" << b.transpose() << endl;
    }
};

int main() {
    int size;
    cout << "Ingrese el tamaño de la matriz n*n: ";
    cin >> size;

    CircuitLC circuit(size);
    circuit.printSystem();
    circuit.measureExecutionTime();

    return 0;
}
