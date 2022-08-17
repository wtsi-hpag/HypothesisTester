#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>
//A highly simplified, specialised and streamlined interface for matrices -- I could have used another implementation but I think this lightweight mimicry helps my case
class squareMatrix
{
	public:
		int Dimension;
		squareMatrix()
		{
			Dimension = 0;
		}
		squareMatrix(int dim)
		{
			Dimension = dim;
			Data = std::vector<std::vector<double>>(dim,std::vector<double>(dim,0.0));
		};
	
		
		double & operator()(unsigned int i, unsigned int j)
		{
			return Data[i][j];
		}
		const double & operator()(unsigned int i, unsigned int j) const
		{
			return Data[i][j];
		}

		std::string Display()
		{
			std::ostringstream s;
			for (int i = 0; i < Dimension; ++i)
			{
				if (i == 0 || i == Dimension -1)
					s<<"| ";
				else
					s<< "( ";

				for (int j = 0; j < Dimension; ++j)
				{
					s<< std::left << std::setw(15) << Data[i][j] << " ";
				}
				if (i == 0 || i == Dimension -1)
					s<<"|\n";
				else
					s<< ")\n";
			}
			return s.str();
		}

		squareMatrix operator*(const squareMatrix & rhs)
		{
			if (rhs.Dimension != Dimension)
			{
				printf("Cannot compute product of different dimension square matrices");
				exit(1);
			}

			squareMatrix out(Dimension);
			for (int i = 0; i < Dimension; ++i)
			{
				for (int j = 0; j < Dimension; ++j)
				{
					for (int k = 0; k < Dimension; ++k)
					{
						out(i,j) += Data[i][k] * rhs(k,j);
					}
				}
			}
			return out;
		}


		double log_LU_Determinant()
		{
			squareMatrix P = Identity(Dimension);
			squareMatrix L = Identity(Dimension);
			squareMatrix A = *this;
			// double p = 1;
			for (int i = 0; i < Dimension-1; ++i)
			{
				// std::cout << "LU step " << i << std::endl;
				double base = A(i,i);
				if (abs(base) < 1e-10)
				{
					squareMatrix newP = Identity(Dimension);
					for (int q = i+1; q < Dimension; ++q)
					{
						if (abs(A(q,i)) > 1e-10)
						{
							newP(i,i) = 0;
							newP(q,q) = 0;
							newP(q,i) = 1;
							newP(i,q) = 1;
							break;
						}
					}
					A = newP * A;
					// L = newP * L;
					P = newP * P;
					base = A(i,i);
				}
				// squareMatrix ell = Identity(Dimension);
				for (int q = i + 1; q < Dimension; ++q)
				{
					L(q,i) = A(q,i)/base;
					double ell = - A(q,i)/base;
					A(q,i) = 0;
					for (int z = i + 1; z < Dimension; ++z)
					{
						A(q,z) += ell * A(i,z);
					}
				}
				// A = ell * A;		


				// std::cout << i << "\n" << A.Display÷ß) << std::endl;
				if (A.isUpperTriangular() && L.isLowerTriangular())
				{
					break;
				}
				
			}
			double p = 0;
			for (int i = 0; i < Dimension; ++i)
			{
				p += std::log(abs((double)A(i,i)));
			}
			return p;
		}

		squareMatrix Transpose()
		{
			squareMatrix out(Dimension);
			for (int i = 0; i < Dimension; ++i)
			{
				out(i,i) = Data[i][i];
				for (int j = i+1; j < Dimension; ++j)
				{
					out(i,j) = Data[j][i];
					out(j,i) = Data[i][j];
				}
			}
			return out;
		}
		static squareMatrix Random(int dim,int bounds)
		{
			squareMatrix s(dim);
			for (int i = 0; i < dim; ++i)
			{
				for (int j = 0; j < dim; ++j)
				{
					s(i,j) = (rand() - RAND_MAX/2) % bounds;
					s(i,j) /= 3;
				}
			}
			return s;
		}

		static squareMatrix Identity(int dim)
		{
			squareMatrix out(dim);
			for (int i = 0; i < dim; ++i)
			{
				out(i,i) = 1;
			}
			return out;
		}

		bool isLowerTriangular()
		{
			for (int i = 0; i < Dimension; ++i)
			{
				for (int j = i+1; j < Dimension; ++j)
				{
					if (abs(Data[i][j]) > 0)
					{
						return false;
					}
				}
			}
			return true;
		}
		bool isUpperTriangular()
		{
			for (int i = 0; i < Dimension; ++i)
			{
				for (int j = i+1; j < Dimension; ++j)
				{
					if (abs(Data[j][i]) > 0)
					{
						return false;
					}
				}
			}
			return true;
		}

	private:
		
		std::vector<std::vector<double>> Data;


};
