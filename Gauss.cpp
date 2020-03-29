#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

void forward_elimination(double matrix[], double f[], int var_links[], int N);
void back_substitution(double matrix[], double f[], int var_links[], double X[], int N);

void solveSystem(double matrix[3*3], double f[3], double &a, double &b, double &c) {
    int N = 3;

    double * X; //Искомые переменные (не переставленные)
    int * var_links; //Отображает переставленные переменные на не переставленные

    X = (double*) malloc(sizeof(double) * N);
    var_links = (int*) malloc(sizeof(int) * N);
    for (int i = 0; i < N; ++i)
        var_links[i] = i;


    forward_elimination(matrix, f, var_links, N);

    back_substitution(matrix, f, var_links, X, N);
    a = X[0];
    b = X[1];
    c = X[2];

    free(X);
    free(var_links);
}


void back_substitution(double matrix[], double f[], int var_links[], double X[], int N) {
    double sum;
    for (int i = N-1; i >= 0; --i) {
        sum = f[i];
        for (int j = i + 1; j < N; ++j) {
            sum -= matrix[N * i + j] * X[var_links[j]];
        }
        X[var_links[i]] = sum;
    }
}


void forward_elimination(double matrix[], double f[], int var_links[], int N) {
    double max_val; // Максимальное значение в строке
    int max_num; // Номер элемента с максимальным значением в столбце
    double mul; // переменная для хранения разных множителей
    int tmp;

    for (int i = 0; i < N; ++i) {
        max_val = abs((double)matrix[(N + 1) * i]);
        max_num = i;
        for (int j = i + 1; j < N; ++j) {
            if ( abs(matrix[N * i + j]) > abs(max_val) ) {
                max_val = abs(matrix[N * i + j]);
                max_num = j;
            }
        }
        // printf("|\nmax: %d\ni: %d\n|", max_num, i);/////////////////////

        tmp = var_links[max_num]; // Меняем местами переменные
        var_links[max_num] = var_links[i];
        var_links[i] = tmp;

        for (int j = 0; j < N; ++j) { // Меняем местами столбцы
            max_val = matrix[j * N + max_num];
            matrix[j * N + max_num] = matrix[j * N + i];
            matrix[j * N + i] = max_val;
        }

        mul = 1 / matrix[(N + 1) * i]; // Разделим один раз, чтобы потом умножать на 1/a[i][i]
        matrix[(N + 1) * i] = 1;
        f[i] *= mul;
        for (int j = i + 1; j < N; ++j) {
            matrix[N * i + j] *= mul;
        }

        for (int k = i + 1; k < N; ++k) { // Обнуляем нижнюю часть столбца
            // k - номер строки
            // j - номер столбца
            mul = matrix[N * k + i];
            matrix[N * k + i] = 0;
            for (int j = i + 1; j < N; ++j)
                matrix[N * k + j] -= matrix[N * i + j] * mul;
            f[k] -= f[i] * mul;
        }

    }
}
