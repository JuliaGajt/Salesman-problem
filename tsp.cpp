#include <iostream>
#include <utility>
#include <vector>
#include <limits>
#include "tsp.hpp"

double get_forbidden_cost(){
    double INF = std::numeric_limits<double>::infinity();
    return INF;
}

void TSP_cost_matrix::rows_cols(matrix &cost_matrix1) {
    std::vector<int> pom1(cost_matrix1.size(), 0);
    std::vector<int> pom2(cost_matrix1.size(), 0);
    for (ulong i = 0; i < cost_matrix1.size(); i++) {
        pom1[i] = i + 1;
        pom2[i] = i + 1;
    }
    rows = pom1;
    cols = pom2;
}

void TSP_cost_matrix::reduce_all_rows(matrix &cost_matrix1, std::vector<double> min) {
    int length = cost_matrix1.size();
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < length; j++) {
            if (cost_matrix1[i][j] != INF) {
                cost_matrix1[i][j] -= min[i];
            }
        }
        low_bound += min[i];
    }
}

std::vector<double> TSP_cost_matrix::find_min(matrix cost_matrix1) {
    int length = cost_matrix1.size();
    vector min(length, 0);
    for (int i = 0; i < length; i++) {
        if (cost_matrix1[i][0] == INF) { min[i] = cost_matrix1[i][1]; }
        else { min[i] = cost_matrix1[i][0]; }
        for (int j = 0; j < length - 1; j++) {
            if (cost_matrix1[i][j + 1] != INF) {
                min[i] = std::min(min[i], cost_matrix1[i][j + 1]);
            }
        }
    }
    return min;
}

void TSP_cost_matrix::transform(matrix &cost_matrix1, int k, int l) {
    std::vector<std::vector<double>> pom(k, std::vector<double>(l, 0));
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < l; j++) {
            pom[j][i] = cost_matrix1[i][j];
        }
    }
    cost_matrix1 = pom;
}

void TSP_cost_matrix::reduce_all_cols(matrix &cost_matrix1) {
    int length = cost_matrix1.size();
    std::vector<std::vector<double>> pom(length, std::vector<double>(length, 0));
    pom = cost_matrix1;
    vector vec(length, 0);
    transform(pom, length, length);
    vec = find_min(pom);
    for (auto e : vec) {
        if (e != 0) {
            reduce_all_rows(pom, vec);
            break;
        }
    }
    transform(pom, length, length);
    cost_matrix1 = pom;
}


double TSP_cost_matrix::sum_min_val(matrix cost_matrix1, int i, int j) {
    int length = cost_matrix1.size();
    vector pom1 = cost_matrix1[i];
    transform(cost_matrix1, length, length);
    vector pom2 = cost_matrix1[j];
    vector min{pom1[0], pom2[0]};

    if (j == 0 || pom1[0] == INF) {
        if (pom1[1] != INF && 1 != j) { min[0] = pom1[1]; }
        else { min[0] = pom1[2]; }
    }
    if (i == 0 || pom2[0] == INF) {
        if (pom2[1] != INF && 1 != i) { min[1] = pom2[1]; }
        else { min[1] = pom2[2]; }
    }
    for (int k = 0; k < length; k++) {
        if (pom1[k] != INF && k != j) {
            min[0] = std::min(min[0], pom1[k]);
        } else { continue; }
    }
    for (int k = 0; k < length; k++) {
        if (pom2[k] != INF && k != i) {
            min[1] = std::min(min[1], pom2[k]);\

        } else { continue; }
    }
    double min1;
    min1 = min[0] + min[1];
    return min1;
}

int TSP_cost_matrix::if_visited(int i, int j) {
    int a = 0;
    for (auto l : visited) {
        if (l[0] == rows[i] || l[1] == cols[j]) {
            a = 1;
            break;
        }
    }
    return a;
}

std::vector<int> TSP_cost_matrix::find_next_path(matrix cost_matrix1) {
    register_t1 points_val;
    std::vector<int> pom{0, 0};
    int length = cost_matrix1.size();
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < length; j++) {
            if (cost_matrix1[i][j] == 0) {
                int b = if_visited(i, j);
                if (b == 0) {
                    std::vector<int> a{i, j};
                    points_val[sum_min_val(cost_matrix1, i, j)] = a;
                }
            }
        }
    }
    vector vec;
    for (auto &e : points_val) {
        vec.push_back(e.first);
    }
    double max = 0;
    for (unsigned long long i = 0; i < (vec.size() - 1); i++) {
        if (vec[i] >= vec[i + 1]) { max = vec[i]; }
        else { max = vec[i + 1]; }
    }
    pom = points_val[max];
    int a1 = pom[0];
    int a2 = pom[1];
    std::vector<int> pom1{rows[a1], cols[a2]};
    visited.push_back(pom1);
    return pom;
}

void TSP_cost_matrix::del_point(matrix &cost_matrix1, int i, int j) {
    int length = cost_matrix1.size();
    int a1 = rows[i];
    int a2 = cols[j];
    std::vector<int> row_p(length, 0);
    std::vector<int> col_p(length, 0);
    for (int l = 0; l < length; l++) {
        row_p[l] = rows[l];
        col_p[l] = cols[l];
        if (l >= i) {
            row_p[l] = rows[l + 1];
        }
        if (l >= j) {
            col_p[l] = cols[l + 1];
        }
    }
    rows = row_p;
    cols = col_p;
    rows[rows.size()] = get_forbidden_cost();
    cols[cols.size()] = get_forbidden_cost();
    for (int o = 0; o < length; o++) {
        if (rows[o] == a2) {
            for (int p = 0; p < length; p++) {
                if (cols[p] == a1) {
                    cost_matrix1[o][p] = get_forbidden_cost();
                    break;
                }
            }
        } else { continue; }
    }
}

void TSP_cost_matrix::delete_row_col(matrix &cost_matrix1, int i, int j) {
    ulong length = cost_matrix1.size() - 1;
    std::vector<std::vector<double>> pom(length, std::vector<double>(length, 0));
    for (ulong k = 0; k < length; k++) {
        for (ulong l = 0; l < length; l++) {
            if (k >= i) { pom[k][l] = cost_matrix1[k + 1][l]; }
            else { pom[k][l] = cost_matrix1[k][l]; }
        }
    }
    for (ulong k = 0; k < length; k++) {
        for (ulong l = 0; l < length; l++) {
            if (l < j) {
                pom[k][l] = pom[k][l];
            } else {
                if (l >= j) {
                    pom[k][l] = pom[k][l + 1];
                    if (l >= (pom.size() - 1) && k >= i) { pom[k][l] = cost_matrix1[k + 1][l + 1]; }
                    else {
                        if (l >= (pom.size() - 1)) { pom[k][l] = cost_matrix1[k][l + 1]; }
                    }
                }
            }
        }
    }
    cost_matrix1 = pom;
    del_point(cost_matrix1, i, j);
}



std::vector<std::vector<int>> TSP_cost_matrix::findBestPath(matrix cost_matrix1) {
    rows_cols(cost_matrix1);
    std::vector<std::vector<double>> a = cost_matrix1;
    int length = cost_matrix1.size();
    for (int i = 0; i < length; i++) {
        vector min = find_min(a);
        reduce_all_rows(a, min);
        reduce_all_cols(a);
        std::vector<int> vec = find_next_path(a);
        delete_row_col(a, vec[0], vec[1]);
    }
    return visited;
}

std::vector<int> tsp(std::vector<std::vector<double>> cost_matrix1) {
    TSP_cost_matrix matrix;
    int length = cost_matrix1.size();
    std::vector<std::vector<int>> vec = matrix.findBestPath(std::move(cost_matrix1));
    std::vector<int> fin{vec[0][0],vec[0][1]};
    int k=1;
    while(fin.size()<(length+1)){
        for(int i = 1; i< length; i++){
            if(vec[i][0]==fin[k]){
                fin.push_back(vec[i][1]);
                k++;
            }
        }
    }
    return fin;
}
