#include <stdio.h>
#include <stdlib.h>

#define TRIVIAL_MULTIPLICATION_BOUND 8

//Пример реализации алгоритма штрассена

void trivial_mul(int** left, int** right, int** result, int n) {
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            result[row][col] = 0;
            for (int mid = 0; mid < n; mid++) {
                result[row][col] += left[row][mid] * right[mid][col];
            }
        }
    }
}

void strassen_mul_2x2(int** lb, int** rb, int** result) {
    int a = lb[0][0];
    int b = lb[0][1];
    int c = lb[1][0];
    int d = lb[1][1];

    int e = rb[0][0];
    int f = rb[0][1];
    int g = rb[1][0];
    int h = rb[1][1];

    int p1 = (a + d) * (e + h);
    int p2 = (c - a) * (f + h);
    int p3 = (b - d) * (g + h);
    int p4 = d * (g - e);
    int p5 = a * (f - h);
    int p6 = c + d;
    int p7 = a + b;

    result[0][0] = p1 + p4 - p5 + p3;
    result[0][1] = p5 + p7;
    result[1][0] = p4 + p6;
    result[1][1] = p1 + p2 - p3 - p6;
}

void strassen_mul(int** left, int** right, int** result, int n) {
    if (n <= TRIVIAL_MULTIPLICATION_BOUND) {
        trivial_mul(left, right, result, n);
        return;
    }

    // Разделяем матрицы на 2x2 блоки
    int new_size = n / 2;
    int** lb[2], ** rb[2], ** res[2];

    for (int i = 0; i < 2; i++) {
        lb[i] = malloc(2 * sizeof(int*));
        rb[i] = malloc(2 * sizeof(int*));
        res[i] = malloc(2 * sizeof(int*));
        for (int j = 0; j < 2; j++) {
            lb[i][j] = malloc(new_size * sizeof(int));
            rb[i][j] = malloc(new_size * sizeof(int));
            res[i][j] = malloc(new_size * sizeof(int));
        }
    }

    // Заполняем блоки
    for (int i = 0; i < new_size; i++) {
        for (int j = 0; j < new_size; j++) {
            lb[0][0][i * new_size + j] = left[i][j];
            lb[0][1][i * new_size + j] = left[i][j + new_size];
            lb[1][0][i * new_size + j] = left[i + new_size][j];
            lb[1][1][i * new_size + j] = left[i + new_size][j + new_size];

            rb[0][0][i * new_size + j] = right[i][j];
            rb[0][1][i * new_size + j] = right[i][j + new_size];
            rb[1][0][i * new_size + j] = right[i + new_size][j];
            rb[1][1][i * new_size + j] = right[i + new_size][j + new_size];
        }
    }

    // Рекурсивно умножаем
    strassen_mul_2x2(lb, rb, res);

    // Собираем результат
    for (int i = 0; i < new_size; i++) {
        for (int j = 0; j < new_size; j++) {
            result[i][j] = res[0][0][i * new_size + j];
            result[i][j + new_size] = res[0][1][i * new_size + j];
            result[i + new_size][j] = res[1][0][i * new_size + j];
            result[i + new_size][j + new_size] = res[1][1][i * new_size + j];
        }
    }

    // Освобождаем память
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            free(lb[i][j]);
            free(rb[i][j]);
            free(res[i][j]);
        }
        free(lb[i]);
        free(rb[i]);
        free(res[i]);
    }
}

int main() {
    int n = 10; // Размер матриц
    int** left, ** right, ** result;

    // Выделяем память для матриц
    left = malloc(n * sizeof(int*));
    right = malloc(n * sizeof(int*));
    result = malloc(n * sizeof(int*));

    for (int i = 0; i < n; i++) {
        left[i] = malloc(n * sizeof(int));
        right[i] = malloc(n * sizeof(int));
        result[i] = malloc(n * sizeof(int));

        for (int j = 0; j < n; j++) {
            left[i][j] = i + 1; // Пример заполнения
            right[i][j] = j + 11; // Пример заполнения
        }
    }

    strassen_mul(left, right, result, n);

    // Вывод результата
    printf("Result:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", result[i][j]);
        }
        printf("\n");
    }

    // Освобождаем память
    for (int i = 0; i < n; i++) {
        free(left[i]);
        free(right[i]);
        free(result[i]);
    }

    free(left);
    free(right);
    free(result);

    return 0;
}


//Approximating matrix multiplication and low-rank approximation

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 120
#define N 69
#define O 40
#define C 69

void print_matrix(double matrix[M][N], int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void random_seed() {
    srand(0); // Установка начального значения для генератора случайных чисел
}

void generate_random_matrix(double matrix[M][N], int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = ((double)rand() / RAND_MAX) * 10;
        }
    }
}

void matrix_multiply(double A[M][N], double B[N][O], double result[M][O]) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < O; j++) {
            result[i][j] = 0;
            for (int k = 0; k < N; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

double vector_norm(double vector[], int size) {
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += vector[i] * vector[i];
    }
    return sqrt(sum);
}

double random_choice(double* array, int size) {
    double rand_value = (double)rand() / RAND_MAX;
    double cumulative_sum = 0;
    for (int i = 0; i < size; i++) {
        cumulative_sum += array[i];
        if (rand_value <= cumulative_sum) {
            return array[i];
        }
    }
    return array[size - 1]; // На случай, если ничего не выбрано
}

int main() {
    random_seed();

    double A[M][N];
    double B[N][O];
    double K[M][C];
    double R[C][O];

    generate_random_matrix(A, M, N);
    generate_random_matrix(B, N, O);

    printf("Matrix A:\n");
    print_matrix(A, M, N);

    printf("Matrix B:\n");
    print_matrix(B, N, O);

    for (int t = 0; t < C; t++) {
        int i_t = rand() % N;

        double random_array[N];
        for (int i = 0; i < N; i++) {
            random_array[i] = (double)rand() / RAND_MAX;
        }

        double sum = 0;
        for (int i = 0; i < N; i++) {
            sum += random_array[i];
        }

        double p_i[N];
        for (int i = 0; i < N; i++) {
            p_i[i] = random_array[i] / sum;
        }

        double choice_p_i = random_choice(p_i, N);

        for (int i = 0; i < M; i++) {
            K[i][t] = A[i][i_t] / sqrt(C * choice_p_i);
        }

        choice_p_i = random_choice(p_i, N);

        for (int j = 0; j < O; j++) {
            R[t][j] = B[i_t][j] / sqrt(C * choice_p_i);
        }
    }

    printf("Matrix C:\n");
    print_matrix(K, M, C);

    printf("Matrix R:\n");
    print_matrix(R, C, O);

    double AB[M][O];
    double CR[M][O];

    matrix_multiply(A, B, AB);
    matrix_multiply(K, R, CR);

    double E_sum = 0;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < O; j++) {
            E_sum += (AB[i][j] - CR[i][j]) * (AB[i][j] - CR[i][j]);
        }
    }

    double E = sqrt(E_sum);

    printf("E: %f\n", E);

    return 0;
}





//PCA


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

void center_data(double* X, double* X_centered, int num_samples, int num_features) {
    double* mean = (double*)malloc(num_features * sizeof(double));
    for (int j = 0; j < num_features; j++) {
        mean[j] = 0.0;
        for (int i = 0; i < num_samples; i++) {
            mean[j] += X[i * num_features + j];
        }
        mean[j] /= num_samples;
    }

    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_features; j++) {
            X_centered[i * num_features + j] = X[i * num_features + j] - mean[j];
        }
    }

    free(mean);
}

void compute_covariance_matrix(double* X_centered, double* covariance_matrix, int num_samples, int num_features) {
    for (int i = 0; i < num_features; i++) {
        for (int j = 0; j < num_features; j++) {
            covariance_matrix[i * num_features + j] = 0.0;
            for (int k = 0; k < num_samples; k++) {
                covariance_matrix[i * num_features + j] += X_centered[k * num_features + i] * X_centered[k * num_features + j];
            }
            covariance_matrix[i * num_features + j] /= (num_samples - 1);
        }
    }
}

// Примерные функции для нахождения собственных значений и векторов

void eigen_decomposition(double* covariance_matrix, double* eigenvalues, double* eigenvectors, int n) {
    // Копируем матрицу ковариации, так как LAPACK может изменять входную матрицу.
    double* A = (double*)malloc(n * n * sizeof(double));
    for (int i = 0; i < n * n; i++) {
        A[i] = covariance_matrix[i];
    }

    // Выделяем память для собственных векторов
    double* work = (double*)malloc(6 * n * sizeof(double)); // Размер рабочей области должен быть достаточным
    int* ipiv = (int*)malloc(n * sizeof(int));
    int info;

    // Вычисляем собственные значения и собственные вектора с помощью LAPACK (dsyev)
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, A, n, eigenvalues); //'V' for eigenvectors, 'U' for upper triangular form

    if (info != 0) {
        fprintf(stderr, "Ошибка при вычислении собственных значений и векторов: %d\n", info);
        exit(1);
    }

    // Копируем собственные векторы в выходной массив
    for (int i = 0; i < n * n; i++) {
        eigenvectors[i] = A[i];
    }

    free(A);
    free(work);
    free(ipiv);
}

void sort_eigenpairs(double* eigenvalues, double* eigenvectors, int n) {
    // Сортировка собственных значений и векторов по убыванию
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (eigenvalues[i] < eigenvalues[j]) {
                double temp_val = eigenvalues[i];
                eigenvalues[i] = eigenvalues[j];
                eigenvalues[j] = temp_val;

                for (int k = 0; k < n; k++) {
                    double temp_vec = eigenvectors[k * n + i];
                    eigenvectors[k * n + i] = eigenvectors[k * n + j];
                    eigenvectors[k * n + j] = temp_vec;
                }
            }
        }
    }
}

void pca(double* X, int num_samples, int num_features, int num_components, double* X_reduced, double* eigenvalues) {
    double* X_centered = (double*)malloc(num_samples * num_features * sizeof(double));
    center_data(X, X_centered, num_samples, num_features);

    double* covariance_matrix = (double*)malloc(num_features * num_features * sizeof(double));
    compute_covariance_matrix(X_centered, covariance_matrix, num_samples, num_features);

    double* eigenvectors = (double*)malloc(num_features * num_features * sizeof(double));
    eigen_decomposition(covariance_matrix, eigenvalues, eigenvectors, num_features);

    sort_eigenpairs(eigenvalues, eigenvectors, num_features);

    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_components; j++) {
            X_reduced[i * num_components + j] = 0.0;
            for (int k = 0; k < num_features; k++) {
                X_reduced[i * num_components + j] += X_centered[i * num_features + k] * eigenvectors[k * num_features + j];
            }
        }
    }

    free(X_centered);
    free(covariance_matrix);
    free(eigenvectors);
}

int main() {
    // Пример использования
    int num_samples = 100;
    int num_features = 5;
    int num_components = 2;

    double* X = (double*)malloc(num_samples * num_features * sizeof(double));

    // Заполнение массива X случайными данными
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_features; j++) {
            X[i * num_features + j] = rand() / (double)RAND_MAX;
        }
    }

    double* X_reduced = (double*)malloc(num_samples * num_components * sizeof(double));
    double* eigenvalues = (double*)malloc(num_components * sizeof(double));

    pca(X, num_samples, num_features, num_components, X_reduced, eigenvalues);

    printf("Преобразованные данные:\n");
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_components; j++) {
            printf("%f ", X_reduced[i * num_components + j]);
        }
        printf("\n");
    }

    printf("Собственные значения:\n");
    for (int i = 0; i < num_components; i++) {
        printf("%f ", eigenvalues[i]);
    }
    printf("\n");

    free(X);
    free(X_reduced);
    free(eigenvalues);

    return 0;
}