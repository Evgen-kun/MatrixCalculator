#include "Matrix.h"
#include <iostream>
using namespace std;

class Matrix
{
public:
	string name;
	float a[50][50], e[50][50], b[50][1], x[50][1];
	int n, coef = 1;
	float det;

	void set__N(int n) {
		this->n = n;
	}

	int get__N() {
		return n;
	}

	float get__det() {
		return det;
	}

	float getB_ByIndex(int i) {
		return b[i][0];
	}

	void set__B(Matrix l) {
		b = l.x;
	}

	void input() {
		cout << "Enter the dimension of the matrix: ";
		cin >> n;

		cout << "\nSolve a system of equations by the Gaussian elimination(1)\nFind determinant(2)\n";
		cout<<"Find inverse matrix(3)\nLU decomposition(4)\n";
		cout << "Solve a system of equations by the LU decomposition(5)\n" << endl;
		int m;
		cin >> m;
		switch (m)
		{
		case 1:
			{
				input_A();
				input_B();
				cout << endl << "Source matrix" << endl;
				output();
				gauss();
				output__det();
				output();
				break;
			}

		case 2:
			{
				input_A();
				for (int i = 0; i < n; i++) {
					b[i][0] = 0;
				}
				cout << endl << "Source matrix" << endl;
				output__A();
				determinant();
				output__det();
				break;
			}
			
		case 3:
			{
				input_A();
				cout << endl << "Source matrix" << endl;
				output__A();
				input_E();
				inverse_matrix();
				break;
			}

		case 4:
			{
				input_A();
				cout << endl << "Source matrix" << endl;
				output__A();
				//input_E();
				Matrix l, u;
				LU_decomposition(l, u);
				break;
			}

		case 5:
			{
				input_A();
				cout << endl << "Source matrix" << endl;
				output__A();
				Matrix l, u;
				LU_decomposition(l, u);
				l.input_B();
				for(int i = 0; i < n; i++) {
					float s = 0;
					for(int j = 0; j < i; j++) {
						s += l.a[i][j] * l.getB_ByIndex(j);
					}
					l.x[i][0] = l.getB_ByIndex(i) - s;
				}
			}
		}
		
	}

	void input_A() {
		cout << endl << "Input the matrix: " << endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cin >> a[i][j];
			}
		}
	}

	void input_B() {
		cout << endl << "Enter a column of constant terms: " << endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < 1; j++) {
				cin >> b[i][j];
			}
		}
	}

	void input_E() {
		for(int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				e[i][j] = (i == j) ? 1 : 0;
			}
		}
	}

	void gauss() {

		det = 1;

		if(check_1()) return;

		triangular_matrix(0);

		find__determinant();

		if (check_2()) return;

		find__answer();

		output__answer();
		cout << endl;
		
	}

	void determinant() {
		det = 1;

		if (check_1()) return;

		triangular_matrix(0);

		find__determinant();

		if (det == 0) {
			cout << "Singular matrix" << endl << endl;
			return;
		}
	}

	void inverse_matrix() {

		det = 1;
		
		if (check_1()) return;

		triangular_matrix(1);

		find__determinant();

		if (det == 0) {
			cout<<"Singular matrix" << endl << endl;
			return;
		}
		
		find__inverse_matrix();

		output__A();
		output__inverse_matrix();
	}

	bool check_1() {
		float s1 = 0, s2 = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				s1 += abs(a[j][i]);
				s2 += abs(a[i][j]);
			}
			if (s1 == 0 || s2 == 0) {
				det = 0;
				cout << "Singular matrix. Infinitely many solutions" << endl << endl;
				return true;
			}
			s1 = 0;
			s2 = 0;
		}
		return false;
	}

	bool check_2() {
		if(det == 0) {

			if (a[n - 1][n - 1] == 0 && b[n - 1][0] != 0) {
				cout << "No solutions" << endl << endl;
			}

			if (a[n - 1][n - 1] == 0 && b[n - 1][0] == 0) {
				cout << "Singular matrix. Infinitely many solutions" << endl << endl;
			}
			return true;
		}
		return false;
	}

	void triangular_matrix(int v) {
		int r = 0;
		
		for (int k = 0; k < n - 1 + r; k++) {
			if (v == 0) output();
			if (v == 1) output__A();
			cout << endl;
			float min = a[0 + k - r][0 + k];
			int min_index = 0 + k - r;

			if (min == 0) {
				int t = 0;
				while (t < n - k + r) {
					if (a[0 + k + t - r][0 + k] != 0) {
						min = a[0 + k + t - r][0 + k];
						min_index = 0 + k + t - r;
						break;
					}
					t++;
				}
				
				if(min == 0) {
					r++;
					continue;
				}
			}
		
			for (int i = 0 + k - r; i < n; i++) {
				if (abs(a[i][0 + k]) < min && a[i][0 + k] != 0) {
					min = a[i][0 + k];
					min_index = i;
				}
			}

			if (min_index != 0 + k - r) {
				float t;
				for (int j = 0 + k; j < n; j++) {
					t = a[0 + k - r][j];
					a[0 + k - r][j] = a[min_index][j];
					a[min_index][j] = t;
				}

				switch (v)
				{
				case 0:
					{
						t = b[0 + k - r][0];
						b[0 + k - r][0] = b[min_index][0];
						b[min_index][0] = t;
						break;
					}

				case 1:
					{
						for (int j = 0; j < n; j++) {
							t = e[0 + k - r][j];
							e[0 + k - r][j] = e[min_index][j];
							e[min_index][j] = t;
						}
						break;
					}
				}
				
				coef *= -1;
			}

			if (a[0 + k - r][0 + k] < 0) {
				for (int j = 0; j < n; j++) {
					a[0 + k - r][j] *= -1.0;
				}

				switch (v)
				{
				case 0:
					{
						b[0 + k - r][0] *= -1.0;
						break;
					}

				case 1:
					{
						for (int j = 0; j < n; j++) {
							e[0 + k - r][j] *= -1.0;
						}
						break;
					}
				}
				
				coef *= -1;
			}

			float above = a[0 + k - r][0 + k], below, temp;
		
			for(int i = 1 + k - r; i < n; i++) {
				below = a[i][0 + k];
				if(below == 0) continue;
				for (int j = 0; j < n; j++) {
					temp = above * a[i][j] - below * a[0 + k - r][j];
					a[i][j] = temp;
				}

				switch (v)
				{
				case 0:
					{
						temp =  above * b[i][0] - below * b[0 + k - r][0];
						b[i][0] = temp;
						break;
					}

				case 1:
					{
						for (int j = 0; j < n; j++) {
							temp = above * e[i][j] - below * e[0 + k - r][j];
							e[i][j] = temp;
						}
						break;
					}
				}
				
				coef *= above;
			}
			
		}
	}

	

	void find__determinant() {
		for(int i = 0; i < n; i++) {
			det *= a[i][i];
		}
		det /= coef;
	}

	void output__A() {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++)
				cout << a[i][j] << "\t";
			cout << endl;
		}
		cout << endl;
	}

	void output__B() {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < 1; j++)
				cout << b[i][j] << "\t";
			cout << endl;
		}
		cout << endl;
	}

	void output__inverse_matrix() {
		cout << "Inverse matrix" << endl << endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++)
				cout << e[i][j] << "\t";
			cout << endl;
		}
		cout << endl;
	}

	void find__inverse_matrix() {
		int k = 1;

		for(int j = n-1; j>=0; j--) {
			float below = a[j][j];
			for(int i = n-1-k; i>=0; i--) {
				float above = a[i][j];
				for (int t = 0; t < n; t++) {
					a[i][t] *= below;
					e[i][t] *= below;
				}
				
				a[i][j] = a[i][j] - above * a[j][j];
				
				for (int t = 0; t < n; t++) {
					e[i][t] = e[i][t] - above * e[j][t];
				}
			}
			k++;
		}

		for(int i = 0; i<n; i++) {
			for(int j = 0; j<n; j++) {
				e[i][j] /= a[i][i];
			}
			a[i][i] /= a[i][i];
		}
		
	}

	void LU_decomposition(Matrix l, Matrix u) {
		l.set__N(n);
		u.set__N(n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
			l.a[i][j] = u.a[i][j] = 0;
		}
		
		
		float sum = 0;
		for (int i = 0; i < n; i++)
		{
			for (int j = i; j < n; j++)
			{
				sum = 0;
				for (int k = 0; k < i; k++)
					sum += l.a[i][k] * u.a[k][j];
				u.a[i][j] = a[i][j] - sum;
			}
			for (int j = i + 1; j < n; j++)
			{
				sum = 0;
				for (int k = 0; k < i; k++)
					sum += l.a[j][k] * u.a[k][i];
				l.a[j][i] = (a[j][i] - sum) / u.a[i][i];
			}
		}

		for(int i = 0; i < n; i++) {
			l.a[i][i] = 1;
		}

		cout << "L-matrix" << endl;
		l.output__A();
		cout << "U-matrix" << endl;
		u.output__A();

		LU_determinant(u);
	}

	void find__answer() {
		for(int i = n - 1; i >= 0; i--) {
			x[i][0] = (b[i][0] - sum(i)) / a[i][i];
		}
	}

	void output__answer() {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < 1; j++)
				cout << "X" << i+1 << "= " << x[i][j] << "\t";
			cout << endl;
		}
		cout << endl;
	}

	void output__det() {
		cout << "determinant = " << get__det() << endl << endl;
	}

	void output() {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++)
				cout << a[i][j] << "\t";
			cout << "|" << b[i][0];
			cout << endl;
		}
		cout << endl;
	}

	void LU_determinant(Matrix u) {
		float s = 0;
		for(int i = 0; i < n; i++) {
			s += u.a[i][i];
		}
		cout << "LU determinant = " << s << endl << endl;
	}

private:
	float sum(int i) {
		float s = 0;
		if (i < n - 1) {
			for (int j = 1; j < n - i; j++) {
				s += a[i][i + j] * x[i + j][0];
			}
			return s;
		}	
		else {
			return 0;
		}
	}
};