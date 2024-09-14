/// ------------------------------------------
/// @file Matrix.h
///
/// @brief Header/Source file for matrix object
///
/// @note Must implement all functions upon definition due to template format
/// ------------------------------------------
#pragma once

#include <stdexcept>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>

/// @brief Templated class for storing, acsessing and performing operations on a matrix of values
template <typename T>
class Matrix
{
    public:
        ///--------------------------------------------------------
        /// @brief Constructor for a matrix object
        ///
        /// @note Memory will only be allocated with no construction taking place
        ///
        /// @tparam T type to store in the matrix, type must be copy
        /// constructable and have all maths operations (+-*/) implemented
        /// with at least itself for any calculations to be performed
        ///
        /// @param rows number of rows to allocate
        /// @param cols number of columns to allocate
        ///
        /// @throws std::invalid_argument if rows/cols < 1
        Matrix(const size_t& rows, const size_t& cols)
        {
            if (rows < 1 or cols < 1)
            {
                throw std::invalid_argument("Cols/Rows of a matrix must be above 0");
            }

            m_cols = cols;
            m_rows = rows;

            m_data = new T[m_cols * m_rows];
        };

        /// @brief Constructor using
        /// @param matData array of arrays to pipe into matrix
        Matrix(const std::initializer_list<std::initializer_list<T>>& matData)
        {
            if (matData.size() == 0)
            {
                throw std::invalid_argument("Cols/Rows of a matrix must be above 0");
            }

            // Check all rows for non zero and identical length
            size_t colLen = matData.begin()->size();
            for (auto row : matData)
            {
                if (row.size() == 0)
                {
                    throw std::invalid_argument("Cols/Rows of a matrix must be above 0");
                }

                if (row.size() != colLen)
                {
                    throw std::invalid_argument("Columns must all be of the same length");
                }
            }

            m_cols = colLen;
            m_rows = matData.size();

            m_data = new T[m_cols * m_rows];
            size_t row = 0;
            for (auto rowData : matData)
            {
                size_t col = 0;
                for (T val : rowData)
                {
                    set(row, col++, val);
                }
                row++;
            }
        };

        ///--------------------------------------------------------
        /// @brief Copy constructor
        ///
        /// @tparam T type stored by matrix
        ///
        /// @param mat reference to copied matrix
        Matrix<T>(Matrix<T> const& mat)
        {
            m_cols = mat.getColCount();
            m_rows = mat.getRowCount();

            m_data = new T[m_cols * m_rows];
            memcpy(m_data, mat.get_data(), m_rows * m_cols * sizeof(T));
        };

        ///--------------------------------------------------------
        /// @brief Destructor
        ~Matrix()
        {
            delete[] m_data;
        };

        /// @brief Assignment operator
        /// @param mat matrix object being assigned from
        /// @return reference to assigned matrix
        Matrix<T>& operator=(Matrix<T> const& mat)
        {
            m_cols = mat.getColCount();
            m_rows = mat.getRowCount();

            // assigned so data will already be present, delete old array to resize
            delete[] m_data;
            m_data = new T[m_cols * m_rows];
            memcpy(m_data, mat.get_data(), m_rows * m_cols * sizeof(T));
            return *this;
        }

        ///--------------------------------------------------------
        /// @brief Operator overload of +, implements matrix addition
        ///
        /// @param mat reference to rval matrix
        ///
        /// @return result of summed matricies
        Matrix<T> operator+(Matrix<T> const& mat)
        {
            if (mat.getColCount() != m_cols or mat.getRowCount() != m_rows)
            {
                throw std::invalid_argument("Matrix addition requires matricies of same dimensions");
            }

            Matrix<T> outMat(m_rows, m_cols);

            // Can sum the memory regions as both have same dimensions
            for (size_t i = 0; i < m_cols * m_rows; i++)
            {
                outMat.get_data()[i] = m_data[i] + mat.get_data()[i];
            }

            return outMat;
        }

        ///--------------------------------------------------------
        /// @brief Operator overload of -, implements matrix subtraction
        ///
        /// @param mat reference to rval matrix
        ///
        /// @return result of subtracted matricies
        Matrix<T> operator-(Matrix<T> const& mat)
        {
            if (mat.getColCount() != m_cols or mat.getRowCount() != m_rows)
            {
                throw std::invalid_argument("Matrix subtraction requires matricies of same dimensions");
            }

            Matrix<T> outMat(m_rows, m_cols);

            // Can sub the memory regions as both have same dimensions
            for (size_t i = 0; i < m_cols * m_rows; i++)
            {
                outMat.get_data()[i] = m_data[i] - mat.get_data()[i];
            }

            return outMat;
        }

        ///--------------------------------------------------------
        /// @brief Operator overload of *, implements matrix dot product
        ///
        /// @param mat reference to rval matrix
        ///
        /// @return result of dot product matricies
        Matrix<T> operator*(Matrix<T> const& mat)
        {
            if (mat.getColCount() != m_cols or mat.getRowCount() != m_rows)
            {
                throw std::invalid_argument("Dot product requires matricies of same dimensions");
            }

            Matrix<T> outMat(m_rows, m_cols);

            // Can product the memory regions as both have same dimensions
            for (size_t i = 0; i < m_cols * m_rows; i++)
            {
                outMat.get_data()[i] = m_data[i] * mat.get_data()[i];
            }

            return outMat;
        }

        /// @brief Operator overload of /, implements number division of matrix
        ///
        /// @param num to divide matrix by
        ///
        /// @return multiplied matrix
        Matrix<T> operator/(T const& num)
        {
            Matrix<T> outMat(m_rows, m_cols);

            for (size_t i = 0; i < m_rows * m_cols; i++)
            {
                outMat.get_data()[i] = m_data[i] / num;
            }

            return outMat;
        };

        ///--------------------------------------------------------
        /// @brief Operator overload of %, implements matrix cross product
        ///
        /// @param mat reference to rval matrix
        ///
        /// @return result of dot product matricies
        Matrix<T> operator%(Matrix<T> const& mat)
        {
            if (m_cols != mat.getRowCount())
            {
                throw std::invalid_argument("Cross product requires matricies of the dimensions: (m,p) % (p,n)");
            }

            Matrix<T> outMat(m_rows, mat.getColCount());

            for (size_t i = 0; i < outMat.getRowCount(); i++)
            {
                for (size_t j = 0; j < outMat.getColCount(); j++)
                {
                    T sum = 0;
                    for (size_t idx = 0; idx < m_cols; idx++)
                    {
                        T product = get(i, idx) * mat.get(idx, j);
                        sum += product;
                    }
                    outMat.set(i,j, sum);
                }
            }

            return outMat;
        };

        ///--------------------------------------------------------
        /// @brief Operator overload of ==, compares two matricies
        ///
        /// @param mat rval mat to compare
        ///
        /// @return are matricies equal in dimension and content?
        bool operator==(Matrix<T> const& mat)
        {
            if (mat.getColCount() != m_cols or mat.getRowCount() != m_rows)
            {
                return false;
            }

            for (size_t i = 0; i < m_cols * m_rows; i++)
            {
                if (m_data[i] != mat.get_data()[i])
                {
                    std::cout << "failed: (" << i / m_cols << "," << i % m_rows << "), " <<
                    m_data[i] << " != " << mat.get_data()[i] << std::endl;
                    return false;
                }
            }

            return true;
        };

        ///--------------------------------------------------------
        /// @brief Operator overload of !=, compares two matricies
        ///
        /// @param mat rval mat to compare
        ///
        /// @return are matricies not equal in dimension or content?
        bool operator!=(Matrix<T> const& mat)
        {
            return !(*this == mat);
        };

        ///--------------------------------------------------------
        /// @brief Operator overload of *, implements number multiplication of matrix
        ///
        /// @note all numbers cast to double to ensure floating point compatability
        ///
        /// @param num to multiply matrix by
        ///
        /// @return multiplied matrix
        Matrix<T> operator*(double const& num)
        {
            Matrix<T> outMat(m_rows, m_cols);

            for (size_t i = 0; i < m_rows * m_cols; i++)
            {
                outMat.get_data()[i] = m_data[i] * num;
            }

            return outMat;
        };

        ///--------------------------------------------------------
        /// @brief Gets the value at the row col position
        ///
        /// @param row to get value from
        /// @param col to get value from
        ///
        /// @returns value at given location
        T get(const size_t& row, const size_t& col) const
        {
            return m_data[_trans_coord(row, col)];
        };

        /// @brief Sets the value at the row col position
        ///
        /// @param row to set value at
        /// @param col to set value at
        /// @param val to set coordinate to
        void set(const size_t& row, const size_t& col, const T& val)
        {
            m_data[_trans_coord(row, col)] = val;
        };

        ///--------------------------------------------------------
        /// @brief Returns a reference of the internal array for the matrix
        ///
        /// @tparam T type of matrix
        ///
        /// @returns reference of internal array
        T* get_data() const
        {
            return m_data;
        };

        ///--------------------------------------------------------
        /// @brief Get the number of rows in the matrix
        ///
        /// @return number of rows in the matrix
        size_t getRowCount() const
        {
            return m_rows;
        };

        ///--------------------------------------------------------
        /// @brief Get the number of columns in the matrix
        ///
        /// @return number of columns in the matrix
        size_t getColCount() const
        {
            return m_cols;
        };

        ///--------------------------------------------------------
        /// @brief Create the transpose of the matrix
        ///
        /// @return the transposed form of the matrix
        Matrix<T> transpose() const
        {
            // Create matrix with transposed dimensions
            Matrix<T> transposeMat(m_cols, m_rows);

            for (size_t i = 0; i < m_rows; i++)
            {
                for (size_t j = 0; j < m_cols; j++)
                {
                    transposeMat.set(j, i, get(i,j));
                }
            }

            return transposeMat;
        };

        ///--------------------------------------------------------
        /// @brief Returns the sub matrix defined by exculding the row and col given
        ///
        /// @param row to exclude when creating new matrix
        /// @param col to exclude when creating new matrix
        ///
        /// @returns matrix of (m-1,n-1) size with given row/col excluded
        Matrix<T> createSubMatrix(const size_t& row, const size_t& col) const
        {
            Matrix<T> outMat(m_rows-1, m_cols-1);

            size_t rowSkip = 0;
            for (size_t i = 0; i < m_rows; i++)
            {
                if (rowSkip == 0 and i == row)
                {
                    rowSkip = 1;
                    continue;
                }

                size_t colSkip = 0;
                for (size_t j = 0; j < m_cols; j++)
                {
                    if (colSkip == 0 and j == col)
                    {
                        colSkip = 1;
                        continue;
                    }

                    outMat.set(i - rowSkip, j - colSkip, get(i, j));
                }
            }

            return outMat;
        };

        ///--------------------------------------------------------
        /// @brief Finds the minor of the matrix at point (i,j)
        ///
        /// @param i row to find minor for
        /// @param j col to find minor for
        ///
        /// @returns value of minor at (i,j)
        T minor(const size_t& i, const size_t& j) const
        {
            return createSubMatrix(i,j).determinant();
        };

        ///--------------------------------------------------------
        /// @brief Finds the cofactor of the matrix at point (i,j)
        ///
        /// @param i row to find cofactor for
        /// @param j col to find cofactor for
        ///
        /// @returns value of cofactor at (i,j)
        T cofactor(const size_t& i, const size_t& j) const
        {
            return minor(i,j) * pow(-1, i + j);
        };

        ///--------------------------------------------------------
        /// @brief Calculates the determinant for the matrix
        ///
        /// @returns value of the determinant for the matrix
        T determinant() const
        {
            if (m_cols != m_rows)
            {
                throw std::invalid_argument("Matrix must be square to have a determinant");
            }

            if (m_cols == 2)
            {
                return get(0,0) * get(1,1) - get(1,0) * get(0,1);
            }
            else if (m_cols == 1)
            {
                return get(0,0);
            }

            size_t workingRow = _find_zeros_row();

            T det = 0;
            for (size_t j = 0; j < m_cols; j++)
            {
                if (get(workingRow, j) != 0)
                {
                    T res = get(workingRow, j) * cofactor(workingRow, j);
                    det += res;
                }
            }

            return det;
        };

        ///--------------------------------------------------------
        /// @brief Calculates the adjoint matrix
        ///
        /// @returns the adjoint matrix of the matrix
        Matrix<T> adjoint()
        {
            Matrix<T> outMat(m_rows, m_cols);

            for (size_t i = 0; i < m_rows; i++)
            {
                for (size_t j = 0; j < m_cols; j++)
                {
                    outMat.set(i,j, cofactor(i,j));
                }
            }

            return outMat.transpose();
        };

        ///--------------------------------------------------------
        /// @brief Calculate the inverse matrix
        ///
        /// @return the inverse matrix
        Matrix<T> inverse()
        {
            T det = determinant();
            if (det == 0)
            {
                throw std::invalid_argument("Matrix determinant is zero, no inverse exists");
            }

            return adjoint() / det;
        };

        ///--------------------------------------------------------
        /// @brief Returns a matrix with each value of the input reciprocated
        ///
        /// @return reciprocated matrix
        Matrix<T> reciprocal()
        {
            Matrix<T> outMat(m_rows, m_cols);

            for (size_t i = 0; i < m_rows * m_cols; i++)
            {
                outMat.get_data()[i] = 1 / m_data[i];
            }

            return outMat;
        };

        ///--------------------------------------------------------
        /// @brief Creates an identity matrix of size len
        ///
        /// @param len side length of the identity matrix
        ///
        /// @return identity matrix of requested size
        static Matrix<T> identity(const size_t& len)
        {
            Matrix<T> id(len, len);
            for (size_t i = 0; i < len; i++)
            {
                for (size_t j = 0; j < len; j++)
                {
                    if (i==j)
                    {
                        id.set(i,j, (T) 1);
                    }
                    else
                    {
                        id.set(i,j, (T) 0);
                    }
                }
            }

            return id;
        }

    private:
        /// @brief Stores all matrix values, one dimensional to exploit memory adjacency benifits
        T* m_data;

        /// @brief the number of columns in the matrix
        size_t m_cols;

        /// @brief the number of rows in the matrix
        size_t m_rows;

        ///--------------------------------------------------------
        /// @brief Translates a coordinate to the index location of the value
        ///
        /// @param row coord to translate
        /// @param col coord to translate
        ///
        /// @returns index to offset memory pointer to
        ///
        /// @throws std::invalid_argument if row/col location is out of bounds
        size_t _trans_coord(const size_t& row, const size_t& col) const
        {
            if (!_check_bounds(row, col))
            {
               std::string err = _gen_coord_err_string(row, col);
               throw std::invalid_argument(err.c_str());
            }

           return row * m_cols + col;
        };

        ///--------------------------------------------------------
        /// @brief Finds the row containing the most zeros
        ///
        /// @note if no zeros are found, row 1 (0) will be returned
        ///
        /// @returns index of row with most zeros
        size_t _find_zeros_row() const
        {
            size_t highestZerosCount = 0;
            size_t zerosRow = 0;

            for (size_t i = 0; i < m_rows; i++)
            {
                size_t zeroCount = 0;
                for (size_t j = 0; j < m_cols; j++)
                {
                    if (get(i,j) == 0)
                    {
                        zeroCount++;
                    }
                }

                if (zeroCount > highestZerosCount)
                {
                    highestZerosCount = zeroCount;
                    zerosRow = i;
                }
            }

            return zerosRow;
        };

        ///--------------------------------------------------------
        /// @brief Determines if a coordinate is out of bounds for this matrix
        ///
        /// @param row of coordinate (i in common notation)
        /// @param col of coordinate (j in common notation)
        ///
        /// @returns is coordinate in bounds of this array?
        bool _check_bounds(const size_t& row, const size_t& col) const
        {
            // True if both coords are in bounds
            return (row < m_rows) and (col < m_cols);
        };

        ///--------------------------------------------------------
        /// @brief Generates the string needed for the error reporting for a bad coord
        ///
        /// @param row offending coord row
        /// @param col offending coord col
        /// @return string of error message
        std::string _gen_coord_err_string(const size_t& row, const size_t& col) const
        {
            std::stringstream err;
            err << "Bad coordinate, (" << row << "," << col << ") is not within the bounds of (" << m_rows - 1 << "," << m_cols - 1 << ")";
            return err.str();
        };
};

///--------------------------------------------------------
/// @brief Overload of <<, used to convert matrix into an output stream
///
/// @param os output stream
/// @param mat matrix to push to output stream
///
/// @return output stream
template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat)
{
    for (size_t i = 0; i < mat.getRowCount(); i++)
    {
        for (size_t j = 0; j < mat.getColCount(); j++)
        {
            os << mat.get(i,j);
            if (j+1 != mat.getColCount())
            {
                os << ", ";
            }
        }

        if (i+1 != mat.getRowCount())
        {
            os << std::endl;
        }
    }

    return os;
}
