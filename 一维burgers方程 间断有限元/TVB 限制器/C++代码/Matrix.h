#ifndef _MATRIX_H_
#define _MATRIX_H_
#pragma once

#include <vector>
#include <cassert>

using namespace std;
/////////////////////////////////////////////////
//////////
/////////////////////////////////////////////////
template <typename T>
class CMatrix
{
public:
	CMatrix();
	~CMatrix();
	void Resize(unsigned int row,unsigned int col);
	inline T& operator ()(unsigned int row,unsigned int col);
private:
	vector<vector<T> > m_Matrix;
	unsigned int m_Col;
	unsigned int m_Row;
};
template <typename T>
CMatrix<T>::CMatrix()
{
	m_Col = 0;
	m_Row = 0;
}
template <typename T>
CMatrix<T>::~CMatrix()
{};
template <typename T>
T& CMatrix<T>::operator ()(unsigned int row,unsigned int col)
{
	return m_Matrix[row][col];
}
template <typename T>
void CMatrix<T>::Resize(unsigned int row, unsigned int col)
{
	m_Row = row;
	m_Col = col;

	m_Matrix.resize(m_Row);
	for (unsigned int i = 0; i != m_Row; i++)
		m_Matrix[i].resize(m_Col);
}
///////////////////////////////
#endif