//
// Created by Julia on 16.11.2019.
//

#ifndef UNTITLED1_TSP_HPP
#define UNTITLED1_TSP_HPP

#include <iostream>
#include <utility>
#include <vector>
#include <limits>
#include <map>

double get_forbidden_cost();

class TSP_cost_matrix {
public:

    // INICJALIZACJA ZMIENNYCH I STAŁYCH
    using vector = std::vector<double>;
    using matrix = std::vector<std::vector<double>>;
    using ulong = unsigned long long;
    double INF = get_forbidden_cost();
    double low_bound = 0;
    using register_t1 = std::map<double, std::vector<int>>;
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<std::vector<int>> visited;


    void rows_cols(matrix &cost_matrix1);

    void reduce_all_rows(matrix &cost_matrix1, std::vector<double> min);    // REDUKACJA WARTOŚCI W RZĘDACH

    vector find_min(matrix cost_matrix1);                                   // ZNAJDYWANIE MINIMUM W MACIERZY

    void reduce_all_cols(matrix &cost_matrix1);                             // REDUKCJA W KOLUMNACH

    static void transform(matrix &cost_matrix1, int k, int l);              // TRANSPONOWANIE MACIERZY

    double sum_min_val(matrix cost_matrix1, int i, int j);                  // SUMA MINIMALNYCH WARTOŚCI (DROGA)

    int if_visited(int i, int j);                                           // SPRAWDZANIE  CZY DANY PKT BYŁ ODWIEDZONY

    std::vector<int> find_next_path(matrix cost_matrix1);                   // ZNAJDYWANIE DROGI

    void delete_row_col(matrix &cost_matrix1, int i, int j);                // USUWANIE KOLUMN I PUNKTÓW
    void del_point(matrix &cost_matrix1, int i, int j);

    std::vector<std::vector<int>> findBestPath(matrix cost_matrix1);

};

std::vector<int> tsp(std::vector<std::vector<double>> cost_matrix1);

#endif //UNTITLED1_TSP_HPP
