#include <iostream>
#include <vector>
#include <cmath>
#include "mpi.h"

std::vector<std::vector<double>> build_up_b(
    double rho,
    double dt,
    double dx, 
    double dy, 
    std::vector<std::vector<double>> u,
    std::vector<std::vector<double>> v, 
    int rank, 
    int size
    ) 
{
    //bの定義と初期化
    std::vector<std::vector<double>> b(u.size(), std::vector<double>(u[0].size(), 0));

    int start = 1;
    int end = b.size() -1;    
    if (rank == 0){
        start = 2;
    }
    if (rank == size - 1){
        end -= 1;
    }

    for (int i = start; i < end; i++){
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
    for (int i = start; i < end; i++){
        int j = b[0].size() - 1;
            b[i][j] = (rho * (1 / dt * ((u[i][0] - u[i][j-1]) / (2 * dx) +
                                        (v[i+1][j] - v[i-1][j] / (2 * dy))) -
                        std::pow(((u[i][0] - u[i][j-1]) / (2 * dx)), 2) -
                        2 * ((u[i+1][j] - u[i-1][j]) / (2 * dy) *
                                (v[i][0] - v[i][j-1]) / (2 * dx)) -
                        std::pow(((v[i+1][j] - v[i-1][j]) / (2 * dy)), 2)));
    }

    //jが最小値の場合j-1が配列外参照になるため書き換える
    for (int i = start; i < end; i++){
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
    std::vector<std::vector<double>> b,
    int rank, 
    int size
    )
{
    std::vector<std::vector<double>> pn;

    int start = 1;
    int end = b.size() -1;    
    if (rank == 0){
        start = 2;
    }
    if (rank == size - 1){
        end -= 1;
    }

    for (int q = 0; q < nit; q++){
        pn = p;
    
        for (int i = start; i < end; i++){
            for (int j = 1; j < p[0].size(); j++){
                p[i][j] = (((pn[i][j+1] + pn[i][j-1]) * std::pow(dy, 2) +
                            (pn[i+1][j] + pn[i-1][j]) * std::pow(dx, 2)) /
                           (2 * (std::pow(dx, 2) + std::pow(dy, 2))) -
                           (std::pow(dx, 2) * std::pow(dy, 2)) * b[i][j]);
            }
        }

        //jが最大値の場合j+1が配列外参照になるため書き換える
        for (int i = start; i < end; i++){
            int j = p[0].size() - 1;
                p[i][j] = (((pn[i][0] + pn[i][j-1]) * std::pow(dy, 2) +
                            (pn[i+1][j] + pn[i-1][j]) * std::pow(dx, 2)) /
                           (2 * (std::pow(dx, 2) + std::pow(dy, 2))) -
                           (std::pow(dx, 2) * std::pow(dy, 2)) * b[i][j]);
        }
        
        //jが最小値の場合j-1が配列外参照になるため書き換える
        for (int i = start; i < end; i++){
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


int main(int argc, char **argv){
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int nx = 41;
    int ny = 41;
    int nt = 10;
    int nit = 50;
    int c = 1;
    double dx = 2 / static_cast<double>(nx - 1);
    double dy = 2 / static_cast<double>(ny - 1);
    
    if (nx < size){
        if (nx <= rank){
            //xで分割するため、xよりrankが多い場合に利用できないことから終了させる。
            MPI_Finalize();
        }
        size = nx;
    }

    int true_nx = nx;
    int nx_start;
    if (nx % size == 0){
        nx = nx / size;
        nx_start = rank * nx;
    }else{
        nx = nx / size;
        nx_start = rank * nx;
        if (rank < (true_nx % size)){
            nx++;
            nx_start += rank;
        }else{
            nx_start += (true_nx % size);
        }
    }
    //後で使うので、nx_startを全てのrankに渡しておく
    int all_nx_start[size + 1];
    MPI_Allgather(&nx_start, 1, MPI_INT, all_nx_start, 1, MPI_INT, MPI_COMM_WORLD);
    all_nx_start[size] = true_nx;

    //std::cout << rank << "," << nx_start << "\n";
    
    //X, Yは計算には使わないので、rank0にのみ定義する
    if (rank == 0){
        //xの初期化
        std::vector<double> x;
        for (int i = 0; i < true_nx; i++) x.push_back(static_cast<double>(i * dx));

        //yの初期化
        std::vector<double> y;
        for (int i = 0; i < ny; i++) y.push_back(static_cast<double>(i * dy));

        //Xの初期化(numpy.meshgridの代替)
        std::vector<std::vector<double>> X(true_nx, std::vector<double>(ny));
        for (int j = 0; j < ny; j++) X[j] = x;

        //Yの初期化(numpy.meshgridの代替)
        std::vector<std::vector<double>> Y(true_nx, std::vector<double>(ny));
        for (int i = 0; i < true_nx; i++){
            for (int j = 0; j < ny ; j++){
                Y[i][j] = y[i];
            }
        }
    }

    double rho = 1;
    double nu = 0.1;
    double F = 1;
    double dt = 0.01;

    std::vector<std::vector<double>> u(nx+2, std::vector<double>(ny, 0));
    std::vector<std::vector<double>> un(nx+2, std::vector<double>(ny, 0));
    
    std::vector<std::vector<double>> v(nx+2, std::vector<double>(ny, 0));
    std::vector<std::vector<double>> vn(nx+2, std::vector<double>(ny, 0));

    std::vector<std::vector<double>> p(nx+2, std::vector<double>(ny, 1));
    std::vector<std::vector<double>> pn(nx+2, std::vector<double>(ny, 1));

    std::vector<std::vector<double>> b(nx+2, std::vector<double>(ny, 0));

    //通信用の変数を定義
    int rank_back;
    if(rank == 0){
        rank_back = size - 1;
    }else{
        rank_back = rank - 1;
    }

    int rank_forward;
    if(rank == size - 1){
        rank_forward = 0;
    }else{
        rank_forward = rank + 1;
    }

    double send_to_back[ny];
    double send_to_forward[ny];

    double recv_from_back[ny];
    double recv_from_forward[ny];

    std::vector<double> temp(ny);
    
    double udiff = 1;
    int stepcount = 0;

    while (udiff > 0.001){
        
        //uを送信
        for (int j = 0; j < ny; j++) send_to_back[j] = u[1][j];
        for (int j = 0; j < ny; j++) send_to_forward[j] = u[u.size()-2][j];
            
        MPI_Send(send_to_back, ny, MPI_DOUBLE, rank_back, 0, MPI_COMM_WORLD);
        MPI_Recv(recv_from_forward, ny, MPI_DOUBLE, rank_forward, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < ny; j++) temp[j] = recv_from_forward[j];
        u[u.size()-1] = temp;

        MPI_Send(send_to_forward, ny, MPI_DOUBLE, rank_forward, 1, MPI_COMM_WORLD);
        MPI_Recv(recv_from_back, ny, MPI_DOUBLE, rank_back, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < ny; j++) temp[j] = recv_from_back[j];
        u[0] = temp;
        
        //vを送信
        for (int j = 0; j < ny; j++) send_to_back[j] = v[1][j];
        for (int j = 0; j < ny; j++) send_to_forward[j] = v[v.size()-2][j];
            
        MPI_Send(send_to_back, ny, MPI_DOUBLE, rank_back, 0, MPI_COMM_WORLD);
        MPI_Recv(recv_from_forward, ny, MPI_DOUBLE, rank_forward, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < ny; j++) temp[j] = recv_from_forward[j];
        v[v.size()-1] = temp;

        MPI_Send(send_to_forward, ny, MPI_DOUBLE, rank_forward, 1, MPI_COMM_WORLD);
        MPI_Recv(recv_from_back, ny, MPI_DOUBLE, rank_back, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < ny; j++) temp[j] = recv_from_back[j];
        v[0] = temp;
        
        //pを送信
        for (int j = 0; j < ny; j++) send_to_back[j] = p[1][j];
        for (int j = 0; j < ny; j++) send_to_forward[j] = p[p.size()-2][j];
            
        MPI_Send(send_to_back, ny, MPI_DOUBLE, rank_back, 0, MPI_COMM_WORLD);
        MPI_Recv(recv_from_forward, ny, MPI_DOUBLE, rank_forward, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < ny; j++) temp[j] = recv_from_forward[j];
        p[p.size()-1] = temp;

        MPI_Send(send_to_forward, ny, MPI_DOUBLE, rank_forward, 1, MPI_COMM_WORLD);
        MPI_Recv(recv_from_back, ny, MPI_DOUBLE, rank_back, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < ny; j++) temp[j] = recv_from_back[j];
        p[0] = temp;

        //bを送信
        for (int j = 0; j < ny; j++) send_to_back[j] = b[1][j];
        for (int j = 0; j < ny; j++) send_to_forward[j] = b[b.size()-2][j];
            
        MPI_Send(send_to_back, ny, MPI_DOUBLE, rank_back, 0, MPI_COMM_WORLD);
        MPI_Recv(recv_from_forward, ny, MPI_DOUBLE, rank_forward, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < ny; j++) temp[j] = recv_from_forward[j];
        b[b.size()-1] = temp;

        MPI_Send(send_to_forward, ny, MPI_DOUBLE, rank_forward, 1, MPI_COMM_WORLD);
        MPI_Recv(recv_from_back, ny, MPI_DOUBLE, rank_back, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < ny; j++) temp[j] = recv_from_back[j];
        b[0] = temp;

        un = u;
        vn = v;

        b = build_up_b(rho, dt, dx, dy, u, v, rank, size);
        p = pressure_poisson_periodic(p, dx, dy, nit, b, rank, size);

        int start = 1;
        int end = u.size() -1;    
        if (rank == 0){
            start = 2;
        }
        if (rank == size - 1){
            end -= 1;
        }
        
        for (int i = start; i < end; i++){
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
        
        
        for (int i = start; i < end; i++){
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
        for (int i = start; i < end; i++){
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
        for (int i = start; i < end; i++){
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
        for (int i = start; i < end; i++){
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
        for (int i = start; i < end; i++){
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

        //分割する前に端となる部分のみ0にする
        if (rank == 0){
            for (int j = 0; j < u[0].size(); j++) u[1][j] = 0;
            for (int j = 0; j < v[0].size(); j++) v[1][j] = 0;
        }
        
        if (rank == size - 1){
            for (int j = 0; j < u[0].size(); j++) u[u.size()-2][j] = 0;
            for (int j = 0; j < v[0].size(); j++) v[v.size()-2][j] = 0;
        }

        //udiffの計算がずれないように、計算用に端に追加された部分を0にする。
        for (int j = 0; j < u[0].size(); j++) u[0][j] = 0;
        for (int j = 0; j < u[0].size(); j++) u[u.size()-1][j] = 0;
        for (int j = 0; j < v[0].size(); j++) v[0][j] = 0;
        for (int j = 0; j < v[0].size(); j++) v[v.size()-1][j] = 0;
        
        //unの状態も前のstepから計算用に追加された部分をコピーして変化しているため、0にしておく。
        for (int j = 0; j < u[0].size(); j++) un[0][j] = 0;
        for (int j = 0; j < u[0].size(); j++) un[u.size()-1][j] = 0;
    
        double local_usum = vec2d_sum(u);
        double local_unsum = vec2d_sum(un);
        double usum, unsum;

        //allreduceを用いて、usum, unsumを全て加算したものを計算して全てのrankに渡す。
        MPI_Allreduce(&local_usum, &usum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);      
        MPI_Allreduce(&local_unsum, &unsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        udiff = (usum - unsum) / usum;
        
        //udiff確認用コード
        /*
        if(rank == 0){
        std::cout << "udiff = " << udiff << "\n";
        }
        */
        
        stepcount++;
    }
    
    if (rank == 0){
        std::cout << "stepcount = " << stepcount << "\n";
    }

    //u, vの最初と最後の要素を削除
    u.erase(u.begin());
    u.pop_back();
    v.erase(v.begin());
    v.pop_back();
    
    for (int i = 1; i < size; i++){
        double temp_send[ny];
        double temp_recv[ny];
        std::vector<double> temp_vec(ny);
        for (int j = 0; j < all_nx_start[i+1] - all_nx_start[i]; j++){
            
            //uを送信
            if (rank == i){ 
                for (int k = 0; k < ny; k++) temp_send[k] = u[j][k];
                MPI_Send(temp_send, ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }else if(rank == 0){
                MPI_Recv(temp_recv, ny, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int k = 0; k < ny; k++) temp_vec[k] = temp_recv[k];
                u.push_back(temp_vec);
            }
            //vを送信
            if (rank == i){ 
                for (int k = 0; k < ny; k++) temp_send[k] = v[j][k];
                MPI_Send(temp_send, ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }else if(rank == 0){
                MPI_Recv(temp_recv, ny, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int k = 0; k < ny; k++) temp_vec[k] = temp_recv[k];
                v.push_back(temp_vec);
            }
        }

    }

    if (rank == 0){
        for (int i = 0; i < v.size(); i++){
            std::cout << "[";
            for (int j = 0; j < ny ; j++){
                std::cout << v[i][j] << ", ";
            }
            std::cout << "]\n";
        }
        std::cout << "\n";
    }

    
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