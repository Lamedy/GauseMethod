#include "GauseClass.h"

int main() {
	setlocale(LC_ALL, "RUS");
	GauseClass A;
	cout << "�������� �������:" << endl;
	A.PrintMatrix(A.OriginalMatrix);
	cout << endl << "������� ��������� ������:" << endl;
	A.PrintVector(A.OriginalVector_B);
	cout << endl << "Det(A) = " << A.FindDeterminant() << "\t";
	if (A.type) {
		cout << "������� �������� �������������" << endl;
		cout << endl;
		cout << "���������:" << "\t\t�������:" << endl;
		for (int i = 0; i < A.size; ++i) {
			cout << "X[" << i + 1 << "] = " << A.Vector_B[i] << "\t\t" << A.ResidualVector[i] << endl;
		}
		cout << endl << "�������� �������:" << endl;
		A.PrintMatrix(A.SingleMatrix);
	}
	else {
		cout << "������� �������� �����������" << endl;
	}
}