#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>

//A highly simplified, specialised and streamlined interface for matrices -- I could have used another implementation but I think this lightweight mimicry helps my case
//Though it can be used for other purposes within the Hypothesis subclasses (i.e. data binning), the squareMatrix class is designed entirely around the computation of the LU form, and hence the log-determinant of the matrix which is needed for the GAI.
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
	
		//!Access to matrix is through (i,j) instead of [i][j], as this allows me to protect the internal structure of the arrays		
		double & operator()(unsigned int i, unsigned int j)
		{
			return Data[i][j];
		}
		//! If squareMatrix is ever const, need this to mimic the non-const one
		const double & operator()(unsigned int i, unsigned int j) const
		{
			return Data[i][j];
		}

		//! Shouldn't be used much in real life, but useful to be able to printout the matrix as a diagnostic.
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

		//What this is here for: computing the logarithm of the determinant via the LU decomposition method
		//The bits that are commented out are those necessary if you actually want the LU decomposition - since we just want the determinant they're only useful for diagnostics and checking the decomposition was successful 
		double log_LU_Determinant()
		{
			// squareMatrix P = Identity(Dimension);
			squareMatrix L = Identity(Dimension);
			squareMatrix A = *this;
			for (int i = 0; i < Dimension-1; ++i)
			{
				std::cout << i << "/"<< Dimension-1 << std::endl;
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
					L = newP * L;
					// P = newP * P;
					base = A(i,i);
				}
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
