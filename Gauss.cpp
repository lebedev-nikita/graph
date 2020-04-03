#include <cstdio>
#include <cstdlib>
// #include <cmath>
#include <utility>

#define N 3

using namespace std;

void inline forward_elimination(double matrix[N][N], double f[N], int var_links[N]);
void inline back_substitution(double matrix[N][N], double f[N],
                              int var_links[N], double X[N]);

void solveSystem(double matrix[N][N], double f[N], double &a, double &b, double &c) {

    double X[N]; //Искомые переменные (не переставленные)
    int var_links[N] = {0,1,2}; //Отображает переставленные переменные на не переставленные

    forward_elimination(matrix, f, var_links);

    back_substitution(matrix, f, var_links, X);

    a = X[0];
    b = X[1];
    c = X[2];
}


void inline back_substitution(double matrix[N][N], double f[N],
                              int var_links[N], double X[N])
{
    double sum;
    int i,j;
    for (i = N-1; i >= 0; --i) {
        sum = f[i];
        for (j = i + 1; j < N; ++j) {
            sum -= matrix[i][j] * X[var_links[j]];
        }
        X[var_links[i]] = sum;
    }
}


void inline forward_elimination(double matrix[N][N], double f[N], int var_links[N])
{
    double max_val; // Максимальное значение в строке
    int max_num; // Номер элемента с максимальным значением в столбце
    double mul; // переменная для хранения разных множителей
    int tmp;
    int i,j,k;

    for (i = 0; i < N; ++i) {
        max_val = abs(matrix[i][i]);
        max_num = i;
        for (j = i + 1; j < N; ++j) {
            if ( abs(matrix[i][j]) > abs(max_val) ) {
                max_val = abs(matrix[i][j]);
                max_num = j;
            }
        }

        swap(var_links[max_num],var_links[i]); // Меняем местами переменные

        for (j = 0; j < N; ++j) { // Меняем местами столбцы
            swap(matrix[j][max_num], matrix[j][i]);
        }

        mul = 1 / matrix[i][i]; // Разделим один раз, чтобы потом умножать на 1/a[i][i]
        matrix[i][i] = 1;
        f[i] *= mul;
        for (j = i + 1; j < N; ++j) {
            matrix[i][j] *= mul;
        }

        for (k = i + 1; k < N; ++k) { // Обнуляем нижнюю часть столбца
            // k - номер строки
            // j - номер столбца
            mul = matrix[k][i];
            matrix[k][i] = 0;
            for (j = i + 1; j < N; ++j)
                matrix[k][j] -= matrix[i][j] * mul;
            f[k] -= f[i] * mul;
        }

    }
}
