#include<iostream>
#include<fstream>
#include<omp.h>
#include<string>
#include<math.h>
#include<locale.h>

//матрица иее копия для сравнения последовтаельного и параллельного алгоритмов
double **mat, **mat2;
int N;//размер матрицы

void printMatrix(double** mat, int n)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			std::cout << mat[i][j] << "  ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

//вернут true при успешном вводе, false при ошибке
bool input(std::string filename)
{

	std::ifstream fin;
	fin.open(filename);
	if (!fin.is_open())
		return false;
	fin >> N;
	mat = new double*[N];
	mat2 = new double*[N];
	for (int i = 0; i < N; i++)
	{
		mat[i] = new double[N];
		mat2[i] = new double[N];
		for (int j = 0; j < N; j++)
			mat[i][j] = mat2[i][j] = 0;
		int pi;
		fin >> pi;
		for (int k = 0; k < pi; k++)
		{
			int j, a;
			fin >> j >> a;
			mat[i][j - 1] = mat2[i][j - 1] = a;
		}
	}
	fin.close();
	return true;
}

//определитель матрицы
//не параллельный
double det(double** mat, int N)
{
	double det = 1;
	const double eps = 1E-8;
	int perc = 0;
	std::cout << "Не параллельный " << 0 << "%";
	for (int i = 0; i < N; i++)
	{
		if (N > 100)
			if (i > (perc + 1)*N / 100)
			{
				perc++;
				std::cout << "\rНе параллельный " << perc << "%";
			}

		if (abs(mat[i][i]) < eps)
		{
			int k = i + 1;
			for (; k < N; k++)
				if (abs(mat[k][i] >= eps))
				{
					break;
				}
			if (k == N)
			{
				std::cout << std::endl;
				return 0;
			}
			/*for (int j = 0; j < N; j++)
		{
			auto tmp = mat[i][j];
			mat[i][j] = mat[k][j];
			mat[k][j] = tmp;
		}*/
			auto tmp = mat[i];
			mat[i] = mat[k];
			mat[k] = tmp;
			det *= -1;
		}
		double diag = mat[i][i];
		for (int j = i; j < N; j++)
			mat[i][j] /= diag;
		for (int k = i + 1; k < N; k++)
		{
			double coeff = mat[k][i];
			if (coeff)
				for (int j = i; j < N; j++)
					mat[k][j] -= coeff * mat[i][j];
		}
		det *= diag;

	}
	std::cout << std::endl;
	return det;
}

//определитель матрицы
//параллельный
double det_par(double** mat, int N)
{
	double det = 1;
	const double eps = 1E-8;
	int perc = 0;
	std::cout << "Параллельный " << 0 << "%";
	for (int i = 0; i < N; i++)
	{
		if (N > 100)
			if (i > (perc + 1)*N / 100)
			{
				perc++;
				std::cout << "\rПараллельный " << perc << "%";
			}

		if (abs(mat[i][i]) < eps)
		{
			int k = i + 1;
			for (; k < N; k++)
				if (abs(mat[k][i] >= eps))
				{
					break;
				}
			if (k == N)
			{
				std::cout << std::endl;
				return 0;
			}

			//#pragma omp parallel for 
			//			for (int j = 0; j < N; j++)
			//			{
			//				auto tmp = mat[i][j];
			//				mat[i][j] = mat[k][j];
			//				mat[k][j] = tmp;
			//			}
			auto tmp = mat[i];
			mat[i] = mat[k];
			mat[k] = tmp;
			det *= -1;
		}
		double diag = mat[i][i];
#pragma omp parallel 
		{

#pragma omp for
			for (int j = i; j < N; j++)
				mat[i][j] /= diag;

#pragma omp barrier

#pragma omp for 
			for (int k = i + 1; k < N; k++)
			{
				if (mat[k][i])
				{
//#pragma omp parallel for
					for (int j = i + 1; j < N; j++)
						mat[k][j] -= mat[k][i] * mat[i][j];
				}


				mat[k][i] = 0;
			}
		}
		det *= diag;

	}
	std::cout << std::endl;
	return det;
}


bool needPrintMatrix(char needPrint)
{
	return  (needPrint == 'Y' || needPrint == 'y' || needPrint == 'н' || needPrint == 'Н');
}


int main()
{
	setlocale(LC_ALL, "Rus");
	std::string file;
	std::cout << "Название файла" << std::endl;
	std::cin >> file;
	char needPrint;
	std::cout << "Выводить матрицы? Y/N" << std::endl;
	std::cin >> needPrint;
	if (input(file))
	{
		if (needPrintMatrix(needPrint))
		{
			std::cout << "Введенная матрица" << std::endl;
			printMatrix(mat, N);
		}


		auto start = omp_get_wtime();
		double d = det(mat, N);
		auto end = omp_get_wtime();
		std::cout << "det = " << d << std::endl;
		std::cout << "t=" << end - start << std::endl;
		std::cout << std::endl;

		start = omp_get_wtime();
		double d2 = det_par(mat2, N);
		end = omp_get_wtime();
		std::cout << "det_par = " << d2 << std::endl;
		std::cout << "t=" << end - start << std::endl;
		std::cout << std::endl;

		if (needPrintMatrix(needPrint))
		{
			std::cout << "После расчета " << std::endl;
			printMatrix(mat, N);
			printMatrix(mat2, N);
		}
	}
	else
		std::cout << "Не найден файл" << std::endl;


	system("pause");
	return 0;
}