#include <stdio.h>
#include <iostream>
#include <ctime>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <chrono>
#include <conio.h>

using namespace std::chrono;
using namespace std;

// количество строк в исходной квадратной матрице
const int MATRIX_SIZE = 1500;

/// Функция InitMatrix() заполняет переданную в качестве 
/// параметра квадратную матрицу случайными значениями
/// matrix - исходная матрица СЛАУ
void InitMatrix(double** matrix)
{
	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		matrix[i] = new double[MATRIX_SIZE + 1];
	}

	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		for (int j = 0; j <= MATRIX_SIZE; ++j)
		{
			matrix[i][j] = rand() % 2500 + 1;
		}
	}
}

/// Функция SerialGaussMethod() решает СЛАУ методом Гаусса 
/// matrix - исходная матрица коэффиициентов уравнений, входящих в СЛАУ,
/// последний столбей матрицы - значения правых частей уравнений
/// rows - количество строк в исходной матрице
/// result - массив ответов СЛАУ
double SerialGaussMethod(double **matrix, const int rows, double* result)
{
	int k;
	double koef;

	// -2- добавление строк кода для измерения времени прямого хода метода Гаусса
	high_resolution_clock::time_point t1, t2;
	t1 = high_resolution_clock::now();

	// прямой ход метода Гаусса
	for (k = 0; k < rows; ++k)
	{
		//
		for (int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];

			for (int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}

	t2 = high_resolution_clock::now();
	duration<double> duration = (t2 - t1);
	std::cout << "Время выполнения прямого хода метода Гаусса последовательно: " << duration.count() << " секунд\n" << std::endl;

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];

		//
		for (int j = k + 1; j < rows; ++j)
		{
			result[k] -= matrix[k][j] * result[j];
		}

		result[k] /= matrix[k][k];
	}
	return duration.count();
}


// -3- создание функции ParallelGaussMethod и введение параллелизма

/// Функция ParallelGaussMethod() решает СЛАУ методом Гаусса 
/// matrix - исходная матрица коэффиициентов уравнений, входящих в СЛАУ,
/// последний столбей матрицы - значения правых частей уравнений
/// rows - количество строк в исходной матрице
/// result - массив ответов СЛАУ
double ParallelGaussMethod(double **matrix, const int rows, double* result)
{
	//int k;
	//double koef;

	high_resolution_clock::time_point t1, t2;
	t1 = high_resolution_clock::now();

	// прямой ход метода Гаусса
	for (int k = 0; k < rows; ++k)
	{
		//

		cilk_for(int i = k + 1; i < rows; ++i)
		{	
			double koef = -matrix[i][k] / matrix[k][k];

			for (int j = k; j <= rows; ++j)
			{				
				matrix[i][j] += koef * matrix[k][j];
			}
			
		}
		
	}

	t2 = high_resolution_clock::now();
	duration<double> duration = (t2 - t1);
	std::cout << "Время выполнения прямого хода метода Гаусса параллельно: " << duration.count() << " секунд\n" << std::endl;

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (int k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];
		cilk::reducer_opadd<double>tmp(result[k]);
		cilk_for (int j = k + 1; j < rows; ++j)
		{
			tmp -= matrix[k][j] * result[j];
		}
		result[k] = tmp.get_value()/matrix[k][k];
		//result[k] /= matrix[k][k];
	}
	
	return duration.count();
	}

/// Функция CopyMatr() делает копию матрицы (одна для параллельного метода, другая - для последовательного) 
/// s_mat - инициализированная последовательная матрица
/// p_mat - параллельная матрица (в которую копируем значения)
/// razm - размерность матрицы
void CopyMatr(double **s_mat, double **p_mat, int razm)
{
	for (int i = 0; i<razm; i++)
	{
		for (int j = 0; j <= razm; j++)
		{
			p_mat[i][j] = s_mat[i][j];
		}
	}
}

int main()
{
	srand((unsigned)time(0));
	setlocale(LC_ALL, "RUS");
	int i;

	// кол-во строк в матрице, приводимой в качестве примера
	//const int test_matrix_lines = 4;

	// кол-во строк в матрице для выполнения заданий
	const int test_matrix_lines = MATRIX_SIZE;

	double **s_test_matrix = new double*[test_matrix_lines]; // для SerialGaussMethod
	double **p_test_matrix = new double*[test_matrix_lines]; // для ParallelGaussMethod

	// цикл по строкам
	for (i = 0; i < test_matrix_lines; ++i)
	{
		// (test_matrix_lines + 1)- количество столбцов в тестовой матрице,
		// последний столбец матрицы отведен под правые части уравнений, входящих в СЛАУ
		s_test_matrix[i] = new double[test_matrix_lines + 1];
		p_test_matrix[i] = new double[test_matrix_lines + 1];
	}

	// массив решений СЛАУ для Serial метода
	double *serial_result = new double[test_matrix_lines];

	// массив решений СЛАУ для Parallel метода
	double *parallel_result = new double[test_matrix_lines];

	/*
	// инициализация тестовой матрицы
	s_test_matrix[0][0] = 2; s_test_matrix[0][1] = 5;  s_test_matrix[0][2] = 4;  s_test_matrix[0][3] = 1;  s_test_matrix[0][4] = 20;
	s_test_matrix[1][0] = 1; s_test_matrix[1][1] = 3;  s_test_matrix[1][2] = 2;  s_test_matrix[1][3] = 1;  s_test_matrix[1][4] = 11;
	s_test_matrix[2][0] = 2; s_test_matrix[2][1] = 10; s_test_matrix[2][2] = 9;  s_test_matrix[2][3] = 7;  s_test_matrix[2][4] = 40;
	s_test_matrix[3][0] = 3; s_test_matrix[3][1] = 8;  s_test_matrix[3][2] = 9;  s_test_matrix[3][3] = 2;  s_test_matrix[3][4] = 37;
	*/

	// -2- инициализация матрицы с количеством строк MATRIX_SIZE
	InitMatrix(s_test_matrix);

	// создание копии матрицы
	CopyMatr(s_test_matrix, p_test_matrix, test_matrix_lines);

	// выполнение последовательной и параллельной реализации
	double serT = SerialGaussMethod(s_test_matrix, test_matrix_lines, serial_result); // 
	double parT = ParallelGaussMethod(p_test_matrix, test_matrix_lines, parallel_result); //

	for (i = 0; i < test_matrix_lines; ++i)
	{
		delete[]s_test_matrix[i];
		delete[]p_test_matrix[i];
	}

	// печать результатов
	// решения СЛАУ
	/*printf("Результаты решения:\n\n");

	for (i = 0; i < test_matrix_lines; ++i)
	{
		printf("x(%d) = %lf        x(%d) = %lf\n", i, serial_result[i], i, parallel_result[i]);
	}*/

	// разница между параллельным и последовательным решением
	/*printf("Результаты решения:\n");

	for (i = 0; i < test_matrix_lines; ++i)
	{
		printf("dif(%d) = %f\n", i, serial_result[i]-parallel_result[i]);
	}*/
	
	printf("\nУскорение = %f\n", serT/parT);
	
	delete[]parallel_result;
	delete[]serial_result;
	_getch();
	return 0;
}
