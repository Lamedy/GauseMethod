#include "GauseClass.h"

int main() {
	setlocale(LC_ALL, "RUS");
	GauseClass A;
	cout << "Исходная матрица:" << endl;
	A.PrintMatrix(A.OriginalMatrix);
	cout << endl << "Столбец свободных членов:" << endl;
	A.PrintVector(A.OriginalVector_B);
	cout << endl << "Det(A) = " << A.FindDeterminant() << "\t";
	if (A.type) {
		cout << "Матрица является невырожденной" << endl;
		cout << endl;
		cout << "Результат:" << "\t\tНевязки:" << endl;
		for (int i = 0; i < A.size; ++i) {
			cout << "X[" << i + 1 << "] = " << A.Vector_B[i] << "\t\t" << A.ResidualVector[i] << endl;
		}
		cout << endl << "Обратная матрица:" << endl;
		A.PrintMatrix(A.SingleMatrix);
	}
	else {
		cout << "Матрица является вырожденной" << endl;
	}
}
