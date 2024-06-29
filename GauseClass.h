#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

// �������� ��� ��������� ��������
vector<float> operator-(const vector<float>& a, const vector<float>& b) {
	if (a.size() != b.size())
		throw("a.size() != b.size()");
	vector<float> c(a.size());
	for (int i = 0; i < a.size(); ++i)
		c[i] = a[i] - b[i];
	return c;
}
// �������� ��� ����� ��������
vector<float> operator+(const vector<float>& a, const vector<float>& b) {
	if (a.size() != b.size())
		throw("a.size() != b.size()");
	vector<float> c(a.size());
	for (int i = 0; i < a.size(); ++i)
		c[i] = a[i] + b[i];
	return c;
}
// �������� ��� ��������� �������� �� �����
vector<float> operator*(const vector<float>& a, const float& b) {
	vector<float> c(a.size());
	for (int i = 0; i < a.size(); ++i)
		c[i] = a[i] * b;
	return c;
}
// �������� ��� ������� ������� �� �����
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
	vector<vector<float>> MainMatrix;				// ������� �������
	vector<vector<float>> OriginalMatrix;			// ������������ �������
	vector<vector<float>> SingleMatrix;				// ��������� �������
	vector<float> Vector_B;							// ���������
	vector<float> OriginalVector_B;					// C������ ��������� ������
	vector<float> ResidualVector;					// ������ �������
	int size;										// ����������� �������
	int PermutationsNumber;							// ����� ������������
	bool type;										// ��� �������, ����������� ��� �������������
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
	// ������ ������ �� �����
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
	// ���������� ����������� � ����
	void SaveResultFile(string FileName) {
		ofstream output(FileName);
		if (output.is_open()) {
			output << "�������� �������:" << endl;
			for (int i = 0; i < size; ++i) {
				for (int q = 0; q < size; ++q) {
					output << " " << setprecision(4) << OriginalMatrix[i][q] << "\t";
				}
				output << endl;
			}
			output << endl << "������� ��������� ������:" << endl;
			for (int i = 0; i < size; ++i) {
				output << " " << OriginalVector_B[i] << endl;
			}
			output << endl << "Det(A) = " << FindDeterminant() << "\t";
			if (type) {
				output << "������� �������� �������������" << endl;
				output << endl;
				output << "���������:" << "\t\t�������:" << endl;
				for (int i = 0; i < size; ++i) {
					output << "X[" << i + 1 << "] = " << Vector_B[i] << "\t\t" << ResidualVector[i] << endl;
				}
				output << endl << "�������� �������:" << endl;
				for (int i = 0; i < size; ++i) {
					for (int q = 0; q < size; ++q) {
						output << " " << setprecision(4) << SingleMatrix[i][q] << "\t";
					}
					output << endl;
				}
			}
			else {
				output << "������� �������� �����������" << endl;
			}
		}
		else {
			cout << "Error save results!";
		}
	}
	// ���������� ��������� �������
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
	// �������������� ������� � ������������ ����
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
	// ����� ����������
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
	// �������� ������ �������� ��������
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
				// ������������ ����� �������� �������
				vector<float> temp = MainMatrix[i];
				MainMatrix[i] = MainMatrix[index];
				MainMatrix[index] = temp;
				PermutationsNumber++;
				// ������������ c������ ��������� ������
				float temp_B = Vector_B[i];
				Vector_B[i] = Vector_B[index];
				Vector_B[index] = temp_B;
				// ������������ ��������� �������
				temp = SingleMatrix[i];
				SingleMatrix[i] = SingleMatrix[index];
				SingleMatrix[index] = temp;
				
				NeedPermutation = false;
			}
		}
	}
	// ����� �������� �������
	void InverseMatrix() {
		vector<vector<float>> tempOriginMatrix = MainMatrix;
		float temp_num;
		// ������ ��� (��������� ������� ������ ����)
		for (int q = 0; q < size - 1; ++q) {
			SingleMatrix[q] = SingleMatrix[q] / tempOriginMatrix[q][q];
			tempOriginMatrix[q] = tempOriginMatrix[q] / tempOriginMatrix[q][q];
			for (int i = q + 1; i < size; ++i) {
				// �������������� �������� �������
				vector<float> temp = tempOriginMatrix[q];
				temp_num = tempOriginMatrix[i][q] / tempOriginMatrix[q][q];
				tempOriginMatrix[i] = tempOriginMatrix[i] - temp * temp_num;
				// �������������� ��������� �������
				temp = SingleMatrix[q];
				SingleMatrix[i] = SingleMatrix[i] - temp * temp_num;
			}
		}
		// �������� ��� (��������� �������� ������� ����)
		vector<vector<float>> some = tempOriginMatrix;
		for (int q = size - 1; q >= 0; --q) {
			SingleMatrix[q] = SingleMatrix[q] / some[q][q];
			tempOriginMatrix[q] = tempOriginMatrix[q] / some[q][q];
			for (int i = q - 1; i >= 0; --i) {
				// �������������� �������� �������
				vector<float> temp = tempOriginMatrix[q];
				temp_num = tempOriginMatrix[i][q] / tempOriginMatrix[q][q];
				tempOriginMatrix[i] = tempOriginMatrix[i] - temp * temp_num;
				// �������������� ��������� �������
				temp = SingleMatrix[q];
				SingleMatrix[i] = SingleMatrix[i] - temp * temp_num;
			}
		}
	}
	// ���������� ��� �������, true - �������������, false - �����������
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
	// ����� ������� �������
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
	// ����� ������������
	float FindDeterminant() {
		float answer = MainMatrix[0][0];
		for (int i = 1; i < size; ++i) {
			answer *= pow(-1, PermutationsNumber) * MainMatrix[i][i];
		}
		return answer;
	}
	// ����� ������� � �������
	void PrintMatrix(vector<vector<float>> &Matrix) {
		for (int i = 0; i < size; ++i) {
			for (int q = 0; q < size; ++q) {
				cout << " " << setprecision(4) << Matrix[i][q] << "\t";
			}
			cout << endl;
		}
	}
	// ����� ������� � �������
	void PrintVector(vector<float> &vec) {
		for (int i = 0; i < size; ++i) {
			cout << " " << vec[i] << endl;
		}
	}
};