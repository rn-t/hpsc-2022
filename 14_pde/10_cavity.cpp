#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

/*
not_final_report_mpiから結構コピーしています。
*/
int main(int argc, char **argv){
    int nx = 41;
    int ny = 41;
    int nt = 500;
    int nit = 50;
    double dx = 2 / static_cast<double>(nx - 1);
    double dy = 2 / static_cast<double>(ny - 1);
    double dt = 0.01;
    double rho = 1;
    double nu = 0.02;

    double x[nx];
    double y[ny];

    for(int i=0; i<nx; i++) x[i] = static_cast<double>(i * dx);
    for(int i=0; i<ny; i++) y[i] = static_cast<double>(i * dy);

    double X[nx][ny];
    double Y[nx][ny];

    double u[ny][nx] = {0};
    double v[ny][nx] = {0};
    double p[ny][nx] = {0};
    double b[ny][nx] = {0};

    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny ; j++){
            X[i][j] = x[j];
        }
    }
    
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny ; j++){
            Y[i][j] = y[i];
        }
    }

    for (int n = 0; n < nt; n++){
        for (int j = 1; j < ny-1; j++){
            for (int i=1; i < nx-1 ;i++){
                b[j][i] = rho * (1 /dt *
                        ((u[j][i+1] - u[j][i-1]) / (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) -
                        pow(((u[j][i+1] - u[j][i-1]) / (2 * dx)),2) - 2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy) *
                         (v[j][i+1] -v[j][i-1]) / (2 * dx)) - pow(((v[j+1][i] - v[j-1][i]) / (2 * dy)), 2));
            }
            if(j == 1){
            for(int i = 1; i < nx-1; i++) cout << b[j][i] << " ";
            cout << "\n";
            }
        }
    
        for(int it = 0; it < nit; it++){
            double pn[ny][nx];
            for (int i = 0; i<ny;i++){
                for (int j = 0; j<nx; j++){
                    pn[i][j] = p[i][j];
                }
            }
            for (int j = 1; j < ny-1; j++){
                for (int i = 1; i < nx-1; i++){
                    p[j][i] = (pow(dy, 2) * (pn[j][i+1] + pn[j][i-1]) +
                            pow(dx, 2) * (pn[j+1][i] + pn[j-1][i]) -
                            b[j][i] * pow(dx, 2) * pow(dy, 2))
                            / (2 * (pow(dx, 2) + pow(dy, 2)));
                }
            }
            for (int j = 0; j < ny; j++) p[j][nx-1] = p[j][nx-2];
            for (int i = 0; i < nx; i++) p[0][i] = p[1][i];
            for (int j = 0; j < ny; j++) p[j][0] = p[j][1];
            for (int i = 0; i < nx; i++) p[ny-1][i] = 0;

        }
        
        double un[ny][nx];
        for (int i = 0; i<ny;i++){
            for (int j = 0; j<nx; j++){
                un[i][j] = u[i][j];
            }
        }
        double vn[ny][nx];
        for (int i = 0; i<ny;i++){
            for (int j = 0; j<nx; j++){
                vn[i][j] = v[i][j];
            }
        }

        for (int j = 1; j < ny-1; j++){
            for (int i = 1; i < nx-1; i++){
                u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i-1])
                                - un[j][i] * dt / dy * (un[j][i] - un[j-1][i])
                                - dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1])
                                + nu * dt / pow(dx, 2) * (un[j][i+1] - 2 * un[j][i] + un[j][i-1])
                                + nu * dt / pow(dy, 2) * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]);
                v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i-1])
                                - vn[j][i] * dt / dy * (vn[j][i] - vn[j-1][i])
                                - dt / (2 * rho * dx) * (p[j+1][i] - p[j-1][i])
                                + nu * dt / pow(dx, 2) * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1])
                                + nu * dt / pow(dy, 2) * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]);
            }
        }
        for (int i = 0; i < nx; i++) u[0][i] = 0;
        for (int j = 0; j < ny; j++) u[j][0] = 0;
        for (int j = 0; j < ny; j++) u[j][nx-1] = 0;
        for (int i = 0; i < nx; i++) u[ny-1][i] = 1;
        
        for (int i = 0; i < nx; i++) v[0][i] = 0;
        for (int i = 0; i < nx; i++) v[ny-1][i] = 0;
        for (int i = 0; i < ny; i++) v[i][0] = 0;
        for (int i = 0; i < ny; i++) v[i][nx-1] = 0;
    
    /*
        for (int i = 0; i < ny; i++){
            std::cout << "[";
            for (int j = 0; j < nx ; j++){
                std::cout << p[i][j] << ", ";
            }
            std::cout << "]\n";
        }
        std::cout << "\n";
    */
    }
}
