#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

// Оператор для вычитания векторов
vector<float> operator-(const vector<float>& a, const vector<float>& b) {
	if (a.size() != b.size())
		throw("a.size() != b.size()");
	vector<float> c(a.size());
	for (int i = 0; i < a.size(); ++i)
		c[i] = a[i] - b[i];
	return c;
}
// Оператор для суммы векторов
vector<float> operator+(const vector<float>& a, const vector<float>& b) {
	if (a.size() != b.size())
		throw("a.size() != b.size()");
	vector<float> c(a.size());
	for (int i = 0; i < a.size(); ++i)
		c[i] = a[i] + b[i];
	return c;
}
// Оператор для умножения векторов на число
vector<float> operator*(const vector<float>& a, const float& b) {
	vector<float> c(a.size());
	for (int i = 0; i < a.size(); ++i)
		c[i] = a[i] * b;
	return c;
}
// Оператор для деления вектора на число
vector<float> operator/(const vector<float>& a, const float& b) {
	vector<float> c(a.size());
	if (b != 0) {
		for (int i = 0; i < a.size(); ++i)
			c[i] = a[i] / b;
	}
	return c;
}
class GauseClass {
public:
	vector<vector<float>> MainMatrix;				// Главная матрица
	vector<vector<float>> OriginalMatrix;			// Оригинальная матрица
	vector<vector<float>> SingleMatrix;				// Еденичная матрица
	vector<float> Vector_B;							// Результат
	vector<float> OriginalVector_B;					// Cтолбец свободных членов
	vector<float> ResidualVector;					// Вектор невязки
	int size;										// Размерность матрицы
	int PermutationsNumber;							// Число перестановок
	bool type;										// Тип матрицы, вырожденная или невырожденная
	GauseClass() {
		PermutationsNumber = 0;
		ReadInputFile("input.txt");
		OriginalMatrix = MainMatrix;
		OriginalVector_B = Vector_B;
		CreateSingleMatrix();
		OutputMainElement();
		InverseMatrix();
		TriangleMatrix();
		type = ChechTypeMatrix();
		if (type) {
			FindResult();
			FindResidualVector();
		}
		SaveResultFile("result.txt");
	}
private:
	// Чтение данных из файла
	void ReadInputFile(string FileName) {
		ifstream input(FileName);
		if (input.is_open())
		{
			input >> size;
			for (int i = 0; i < size; ++i) {
				vector<float> temp;
				for (int q = 0; q < size; ++q) {
					temp.push_back(0);
					input >> temp[q];
				}
				MainMatrix.push_back(temp);
			}
			for (int i = 0; i < size; ++i) {
				Vector_B.push_back(0);
				input >> Vector_B[i];
			}
		}
		else {
			cout << "input.txt not found";
		}
		input.close();
	}
	// Сохранение результатов в файл
	void SaveResultFile(string FileName) {
		ofstream output(FileName);
		if (output.is_open()) {
			output << "Исходная матрица:" << endl;
			for (int i = 0; i < size; ++i) {
				for (int q = 0; q < size; ++q) {
					output << " " << setprecision(4) << OriginalMatrix[i][q] << "\t";
				}
				output << endl;
			}
			output << endl << "Столбец свободных членов:" << endl;
			for (int i = 0; i < size; ++i) {
				output << " " << OriginalVector_B[i] << endl;
			}
			output << endl << "Det(A) = " << FindDeterminant() << "\t";
			if (type) {
				output << "Матрица является невырожденной" << endl;
				output << endl;
				output << "Результат:" << "\t\tНевязки:" << endl;
				for (int i = 0; i < size; ++i) {
					output << "X[" << i + 1 << "] = " << Vector_B[i] << "\t\t" << ResidualVector[i] << endl;
				}
				output << endl << "Обратная матрица:" << endl;
				for (int i = 0; i < size; ++i) {
					for (int q = 0; q < size; ++q) {
						output << " " << setprecision(4) << SingleMatrix[i][q] << "\t";
					}
					output << endl;
				}
			}
			else {
				output << "Матрица является вырожденной" << endl;
			}
		}
		else {
			cout << "Error save results!";
		}
	}
	// Заполнение еденичной матрицы
	void CreateSingleMatrix() {
		for (int i = 0; i < size; ++i) {
			vector<float> temp;
			for (int q = 0; q < size; ++q) {
				temp.push_back(0);
			}
			SingleMatrix.push_back(temp);
			SingleMatrix[i][i] = 1;
		}
	}
	// Преобразование матрицы к треугольному виду
	void TriangleMatrix() {
		float temp_num;
		for (int q = 0; q < size-1; ++q) {
			for (int i = q+1; i < size; ++i) {
				vector<float> temp = MainMatrix[q];
				temp_num = MainMatrix[i][q] / MainMatrix[q][q];
				Vector_B[i] = Vector_B[i] - Vector_B[q] * temp_num;
				MainMatrix[i] = MainMatrix[i] - temp * temp_num;
			}
		}
	}
	// Поиск результата
	void FindResult() {
		float sum;
		for (int i = size-1; i >= 0; --i) {
			sum = Vector_B[i];
			for (int q = i; q < size-1; ++q) {
				sum -= MainMatrix[i][q+1] * Vector_B[q+1];
			}
			Vector_B[i] = sum / MainMatrix[i][i];
		}
	}
	// Алгоритм выбора главного элемента
	void OutputMainElement() {
		int index;
		bool NeedPermutation = 0;
		for (int i = 0; i < size; ++i) {
			float max = abs(MainMatrix[i][i]);
			for (int q = i; q < size; ++q) {
				if (abs(MainMatrix[q][i]) > max) {
					max = abs(MainMatrix[q][i]);
					index = q;
					NeedPermutation = true;
				}
			}
			if (NeedPermutation) {
				// Перестановка строк основной матрицы
				vector<float> temp = MainMatrix[i];
				MainMatrix[i] = MainMatrix[index];
				MainMatrix[index] = temp;
				PermutationsNumber++;
				// Перестановка cтолбца свободных членов
				float temp_B = Vector_B[i];
				Vector_B[i] = Vector_B[index];
				Vector_B[index] = temp_B;
				// Перестановка еденичной матрицы
				temp = SingleMatrix[i];
				SingleMatrix[i] = SingleMatrix[index];
				SingleMatrix[index] = temp;
				
				NeedPermutation = false;
			}
		}
	}
	// Поиск обратной матрицы
	void InverseMatrix() {
		vector<vector<float>> tempOriginMatrix = MainMatrix;
		float temp_num;
		// Прямой ход (Зануление нижнего левого угла)
		for (int q = 0; q < size - 1; ++q) {
			SingleMatrix[q] = SingleMatrix[q] / tempOriginMatrix[q][q];
			tempOriginMatrix[q] = tempOriginMatrix[q] / tempOriginMatrix[q][q];
			for (int i = q + 1; i < size; ++i) {
				// Преобразование основной матрицы
				vector<float> temp = tempOriginMatrix[q];
				temp_num = tempOriginMatrix[i][q] / tempOriginMatrix[q][q];
				tempOriginMatrix[i] = tempOriginMatrix[i] - temp * temp_num;
				// Преоброзование еденичной матрицы
				temp = SingleMatrix[q];
				SingleMatrix[i] = SingleMatrix[i] - temp * temp_num;
			}
		}
		// Обратный ход (Зануление верхнего правого угла)
		vector<vector<float>> some = tempOriginMatrix;
		for (int q = size - 1; q >= 0; --q) {
			SingleMatrix[q] = SingleMatrix[q] / some[q][q];
			tempOriginMatrix[q] = tempOriginMatrix[q] / some[q][q];
			for (int i = q - 1; i >= 0; --i) {
				// Преобразование основной матрицы
				vector<float> temp = tempOriginMatrix[q];
				temp_num = tempOriginMatrix[i][q] / tempOriginMatrix[q][q];
				tempOriginMatrix[i] = tempOriginMatrix[i] - temp * temp_num;
				// Преоброзование еденичной матрицы
				temp = SingleMatrix[q];
				SingleMatrix[i] = SingleMatrix[i] - temp * temp_num;
			}
		}
	}
	// Определить тип матрицы, true - невырожденная, false - вырожденная
	bool ChechTypeMatrix() {
		bool answer;
		if (FindDeterminant() == 0) {
			answer = false;
		}
		else {
			answer = true;
		}
		return answer;
	}
	// Поиск вектора невязки
	void FindResidualVector() {
		ResidualVector = OriginalVector_B;
		for (int i = 0; i < size; ++i) {
			float sum = 0;
			for (int q = 0; q < size; ++q) {
				sum += OriginalMatrix[i][q] * Vector_B[q];
			}
			ResidualVector[i] -= sum;
			sum = 0;
		}
	}
public:
	// Поиск определителя
	float FindDeterminant() {
		float answer = MainMatrix[0][0];
		for (int i = 1; i < size; ++i) {
			answer *= pow(-1, PermutationsNumber) * MainMatrix[i][i];
		}
		return answer;
	}
	// Вывод матрицы в консоль
	void PrintMatrix(vector<vector<float>> &Matrix) {
		for (int i = 0; i < size; ++i) {
			for (int q = 0; q < size; ++q) {
				cout << " " << setprecision(4) << Matrix[i][q] << "\t";
			}
			cout << endl;
		}
	}
	// Вывод вектора в консоль
	void PrintVector(vector<float> &vec) {
		for (int i = 0; i < size; ++i) {
			cout << " " << vec[i] << endl;
		}
	}
};