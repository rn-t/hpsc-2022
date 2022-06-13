#include <iostream>
#include <vector>
#include <cmath>

std::vector<std::vector<double>> build_up_b(
    double rho,
    double dt,
    double dx, 
    double dy, 
    std::vector<std::vector<double>> u,
    std::vector<std::vector<double>> v
    ) 
{
    //bの定義と初期化
    std::vector<std::vector<double>> b(u.size(), std::vector<double>(u[0].size(), 0));

    for (int i = 1; i < b.size() - 1; i++){
        for (int j = 1; j < b[0].size() - 1; j++){
            b[i][j] = (rho * (1 / dt * ((u[i][j+1] - u[i][j-1]) / (2 * dx) +
                                        (v[i+1][j] - v[i-1][j] / (2 * dy))) -
                        std::pow(((u[i][j+1] - u[i][j-1]) / (2 * dx)), 2) -
                        2 * ((u[i+1][j] - u[i-1][j]) / (2 * dy) *
                                (v[i][j+1] - v[i][j-1]) / (2 * dx)) -
                        std::pow(((v[i+1][j] - v[i-1][j]) / (2 * dy)), 2)));
        }
    }

    //jが最大値の場合j+1が配列外参照になるため書き換える
    for (int i = 1; i < b.size() - 1; i++){
        int j = b[0].size() - 1;
            b[i][j] = (rho * (1 / dt * ((u[i][0] - u[i][j-1]) / (2 * dx) +
                                        (v[i+1][j] - v[i-1][j] / (2 * dy))) -
                        std::pow(((u[i][0] - u[i][j-1]) / (2 * dx)), 2) -
                        2 * ((u[i+1][j] - u[i-1][j]) / (2 * dy) *
                                (v[i][0] - v[i][j-1]) / (2 * dx)) -
                        std::pow(((v[i+1][j] - v[i-1][j]) / (2 * dy)), 2)));
    }

    //jが最小値の場合j-1が配列外参照になるため書き換える
    for (int i = 1; i < b.size() - 1; i++){
        int j = 0;
        int jmax = b[0].size() - 1;
            b[i][j] = (rho * (1 / dt * ((u[i][j+1] - u[i][jmax]) / (2 * dx) +
                                        (v[i+1][j] - v[i-1][j] / (2 * dy))) -
                        std::pow(((u[i][j+1] - u[i][jmax]) / (2 * dx)), 2) -
                        2 * ((u[i+1][j] - u[i-1][j]) / (2 * dy) *
                                (v[i][j+1] - v[i][jmax]) / (2 * dx)) -
                        std::pow(((v[i+1][j] - v[i-1][j]) / (2 * dy)), 2)));
    }  
    return b;
}

//pythonではグローバル変数となっているbやnitを引数として追加
std::vector<std::vector<double>> pressure_poisson_periodic(
    std::vector<std::vector<double>> p,
    double dx, 
    double dy,
    int nit,
    std::vector<std::vector<double>> b
    )
{
    std::vector<std::vector<double>> pn;
    for (int q = 0; q < nit; q++){
        pn = p;
    
        for (int i = 1; i < p.size() - 1; i++){
            for (int j = 1; j < p[0].size(); j++){
                p[i][j] = (((pn[i][j+1] + pn[i][j-1]) * std::pow(dy, 2) +
                            (pn[i+1][j] + pn[i-1][j]) * std::pow(dx, 2)) /
                           (2 * (std::pow(dx, 2) + std::pow(dy, 2))) -
                           (std::pow(dx, 2) * std::pow(dy, 2)) * b[i][j]);
            }
        }

        //jが最大値の場合j+1が配列外参照になるため書き換える
        for (int i = 1; i < p.size() - 1; i++){
            int j = p[0].size() - 1;
                p[i][j] = (((pn[i][0] + pn[i][j-1]) * std::pow(dy, 2) +
                            (pn[i+1][j] + pn[i-1][j]) * std::pow(dx, 2)) /
                           (2 * (std::pow(dx, 2) + std::pow(dy, 2))) -
                           (std::pow(dx, 2) * std::pow(dy, 2)) * b[i][j]);
        }
        
        //jが最小値の場合j-1が配列外参照になるため書き換える
        for (int i = 1; i < p.size() - 1; i++){
            int j = 0;
            int jmax = p[0].size() - 1;
                p[i][j] = (((pn[i][j+1] + pn[i][jmax]) * std::pow(dy, 2) +
                            (pn[i+1][j] + pn[i-1][j]) * std::pow(dx, 2)) /
                           (2 * (std::pow(dx, 2) + std::pow(dy, 2))) -
                           (std::pow(dx, 2) * std::pow(dy, 2)) * b[i][j]);
        }

    }
    return p;
}

//2dのベクトルの合計を出す関数
double vec2d_sum(std::vector<std::vector<double>> input_vec){
    double output_scalar = 0;
    for (int i = 0; i < input_vec.size(); i++){
        for (int j = 0; j < input_vec[0].size(); j++){
            output_scalar += input_vec[i][j];
        }
    }
    return output_scalar;
}


int main(){
    int nx = 41;
    int ny = 41;
    int nt = 10;
    int nit = 50;
    int c = 1;
    double dx = 2 / static_cast<double>(nx - 1);
    double dy = 2 / static_cast<double>(ny - 1);
    
    //xの初期化
    std::vector<double> x;
    for (int i = 0; i < nx; i++) x.push_back(static_cast<double>(i * dx));

    //yの初期化
    std::vector<double> y;
    for (int i = 0; i < ny; i++) y.push_back(static_cast<double>(i * dy));

    //Xの初期化(numpy.meshgridの代替)
    std::vector<std::vector<double>> X(nx, std::vector<double>(ny));
    for (int j = 0; j < ny; j++) X[j] = x;

    //Yの初期化(numpy.meshgridの代替)
    std::vector<std::vector<double>> Y(nx, std::vector<double>(ny));
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny ; j++){
            Y[i][j] = y[i];
        }
    }

    double rho = 1;
    double nu = 0.1;
    double F = 1;
    double dt = 0.01;

    std::vector<std::vector<double>> u(nx, std::vector<double>(ny, 0));
    std::vector<std::vector<double>> un(nx, std::vector<double>(ny, 0));
    
    std::vector<std::vector<double>> v(nx, std::vector<double>(ny, 0));
    std::vector<std::vector<double>> vn(nx, std::vector<double>(ny, 0));

    std::vector<std::vector<double>> p(nx, std::vector<double>(ny, 1));
    std::vector<std::vector<double>> pn(nx, std::vector<double>(ny, 1));

    std::vector<std::vector<double>> b(nx, std::vector<double>(ny, 0));

    double udiff = 1;
    int stepcount = 0;

    while (udiff > 0.001){
        un = u;
        vn = v;

        b = build_up_b(rho, dt, dx, dy, u, v);
        p = pressure_poisson_periodic(p, dx, dy, nit, b);

        for (int i = 1; i < u.size() - 1; i++){
            for (int j = 1; j < u[0].size() - 1; j++){
                u[i][j] = (un[i][j] -
                           un[i][j] * dt / dx *
                          (un[i][j] - un[i][j-1]) -
                           vn[i][j] * dt / dy *
                          (un[i][j] - un[i-1][j]) -
                           dt / (2 * rho * dx) *
                          (p[i][j+1] - p[i][j-1]) + 
                           nu * (dt / std::pow(dx, 2) *
                          (un[i][j+1] - 2 * un[i][j] + un[i][j-1]) +
                           dt / std::pow(dy, 2) *
                          (un[i+1][j] - 2 * un[i][j] + un[i-1][j])) +
                           F * dt);
            }
        }
        
        for (int i = 1; i < v.size() - 1; i++){
            for (int j = 1; j < v[0].size() - 1; j++){
                v[i][j] = (vn[i][j] -
                           un[i][j] * dt / dx *
                          (vn[i][j] - vn[i][j-1]) -
                           vn[i][j] * dt / dy *
                          (vn[i][j] - vn[i-1][j]) -
                           dt / (2 * rho * dy) *
                          (p[i+1][j] - p[i-1][j]) + 
                           nu * (dt / std::pow(dx, 2) *
                          (vn[i][j+1] - 2 * vn[i][j] + vn[i][j-1]) +
                           dt / std::pow(dy, 2) *
                          (vn[i+1][j] - 2 * vn[i][j] + vn[i-1][j])));
            }
        }
        
        //jが最大の場合
        for (int i = 1; i < u.size() - 1; i++){
            int j = u[0].size() - 1;
                u[i][j] = (un[i][j] -
                           un[i][j] * dt / dx *
                          (un[i][j] - un[i][j-1]) -
                           vn[i][j] * dt / dy *
                          (un[i][j] - un[i-1][j]) -
                           dt / (2 * rho * dx) *
                          (p[i][0] - p[i][j-1]) + 
                           nu * (dt / std::pow(dx, 2) *
                          (un[i][0] - 2 * un[i][j] + un[i][j-1]) +
                           dt / std::pow(dy, 2) *
                          (un[i+1][j] - 2 * un[i][j] + un[i-1][j])) +
                           F * dt);
        }
        
        //jが最小値の場合
        for (int i = 1; i < u.size() - 1; i++){
            int j = 0;
            int jmax = u[0].size() - 1;
                u[i][j] = (un[i][j] -
                           un[i][j] * dt / dx *
                          (un[i][j] - un[i][jmax]) -
                           vn[i][j] * dt / dy *
                          (un[i][j] - un[i-1][j]) -
                           dt / (2 * rho * dx) *
                          (p[i][0] - p[i][jmax]) + 
                           nu * (dt / std::pow(dx, 2) *
                          (un[i][0] - 2 * un[i][j] + un[i][jmax]) +
                           dt / std::pow(dy, 2) *
                          (un[i+1][j] - 2 * un[i][j] + un[i-1][j])) +
                           F * dt);
        }

        //jが最大の場合
        for (int i = 1; i < v.size() - 1; i++){
            int j = v[0].size() - 1;
                v[i][j] = (vn[i][j] -
                           un[i][j] * dt / dx *
                          (vn[i][j] - vn[i][j-1]) -
                           vn[i][j] * dt / dy *
                          (vn[i][j] - vn[i-1][j]) -
                           dt / (2 * rho * dy) *
                          (p[i+1][j] - p[i-1][j]) + 
                           nu * (dt / std::pow(dx, 2) *
                          (vn[i][0] - 2 * vn[i][j] + vn[i][j-1]) +
                           dt / std::pow(dy, 2) *
                          (vn[i+1][j] - 2 * vn[i][j] + vn[i-1][j])));
        }

        //jが最小の場合
        for (int i = 1; i < v.size() - 1; i++){
            int j = 0;
            int jmax = v[0].size() - 1;
                v[i][j] = (vn[i][j] -
                           un[i][j] * dt / dx *
                          (vn[i][j] - vn[i][jmax]) -
                           vn[i][j] * dt / dy *
                          (vn[i][j] - vn[i-1][j]) -
                           dt / (2 * rho * dy) *
                          (p[i+1][j] - p[i-1][j]) + 
                           nu * (dt / std::pow(dx, 2) *
                          (vn[i][0] - 2 * vn[i][j] + vn[i][jmax]) +
                           dt / std::pow(dy, 2) *
                          (vn[i+1][j] - 2 * vn[i][j] + vn[i-1][j])));
        }

        for (int j = 0; j < u[0].size(); j++) u[0][j] = 0;
        for (int j = 0; j < u[0].size(); j++) u[u.size()-1][j] = 0;

        for (int j = 0; j < v[0].size(); j++) v[0][j] = 0;
        for (int j = 0; j < v[0].size(); j++) v[v.size()-1][j] = 0;
    
        udiff = (vec2d_sum(u) - vec2d_sum(un)) / vec2d_sum(u);
        
        stepcount++;
    }
    
    std::cout << "stepcount = " << stepcount << "\n";
    
    for (int i = 0; i < nx; i++){
        std::cout << "[";
        for (int j = 0; j < ny ; j++){
            std::cout << u[i][j] << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << "\n";
    
    
    //確認用コード
    /*
    for (int i = 0; i < nx; i++){
        std::cout << "[";
        for (int j = 0; j < ny ; j++){
            std::cout << u[i][j] << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << "\n";
    */    


}