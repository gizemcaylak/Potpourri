#ifndef MATRIX_MATH
#define MATRIX_MATH

#include "Matrix.h"

#include <map>



namespace MatrixMath
{   
    
    /** @brief Performs context sensitive element wise modifications on input matrixes and returns a new matrix
     * 
     *  @param m1 The first input matrix
     *  @param m2 The second input matrix, must at least one dimension must equal the first while 
     *  the other is size one, or both must be size one
     *  @return A new matrix as a result of the element wise modification
     */
    template<class T>
    static Matrix<T> elementWiseOperation(const Matrix<T> & m1, const Matrix<T> & m2, T (*foo)(const T &, const T &));
    template<class T>
    static T add(const T & v1, const T & v2);
    template<class T>
    static T subtract(const T & v1, const T & v2);
    template<class T>
    static T multiply(const T & v1, const T & v2);
    template<class T>
    static T divide(const T & v1, const T & v2);
    
    template<class T>
    static Matrix<T> transpose(const Matrix<T> & m);
    
    template<class T>
    static Matrix<T> MatrixProduct(const Matrix<T> & m1, const Matrix<T> & m2);
    
    
    /** @brief Makes a matrix such that every rows value for a given column is unique
     *  
     *  Makes unique by dropping all rows with a key equivalent to one that already exists. Where the key is defined
     *  as the vector of values made up from the passed columns in each row
     * 
     *  @param m The input matrix
     *  @param c The vector of columns to base uniqueness on
     * 
     *  @return The new matrix with unique rows based on the specified columns
     */
    template<class T>
    static Matrix<T> makeRowsUniqueOnColumns(const Matrix<T> & m, vector<int> c);
    
    //REFORMULATE TO WORK WITH VECTORS OF FUNCTION POINTERS
    template<class T>
    static Matrix<T> removeRowsWith(const Matrix<T> & m, const vector<int> & c, const vector<T> & v, bool (*foo1)(const bool, const bool), bool (*foo2)(const T &, const T &));
    
    template<class T>
    static bool EQUALS(const T & t1, const T & t2);
    static bool AND(const bool b1, const bool b2);
    static bool OR(const bool b1, const bool b2);
    
    template<class T>
    static pair<Matrix<T>, Matrix<T> > permuteColumns(const Matrix<T> & m1, const Matrix<T> & m2 );
    
/*-----------------------Implementation-----------------------------*/
    


    template<class T>
    Matrix<T> elementWiseOperation(const Matrix<T> & m1, const Matrix<T> & m2, T (*foo)(const T &, const T &))
    {
        //Extract both matrix dimensions
        int r1 = m1.dim(0);  int r2 = m2.dim(0); int c1 = m1.dim(1);int c2 = m2.dim(1);
        
        //If they are equal sized
        Matrix<T> n_m(r1, c1);
        if((r1 == r2 && c1 == c2) )
        {
            for(int i=0; i<r1; i++)
                for(int j=0; j<c1; j++)
                    n_m.a(i,j) = foo(m1.e(i,j) , m2.e(i,j));

            return n_m;
        }
        //Or if the second one is a single element
        else if((r2 == 1 && c2 == 1))
        {
            for(int i=0; i<r1; i++)
                for(int j=0; j<c1; j++)
                    n_m.a(i,j) = foo(m1.e(i,j) , m2.e(0,0));   
                    
            return n_m;
        }
        //Or if the second has the same number of columns
        else if(r2 == 1 && c2 == c1)
        {
            for(int i=0; i<r1; i++)
                for(int j=0; j<c1; j++)
                    n_m.a(i,j) = foo(m1.e(i,j) , m2.e(0,j));   
                    
            return n_m;
        }
        //Or if the second has the same number of rows
        else if(c2 == 1 && r2 == r1)
        {
            for(int i=0; i<r1; i++)
                for(int j=0; j<c1; j++)
                    n_m.a(i,j) = foo(m1.e(i,j) , m2.e(i,0));   
                    
            return n_m;
        }
        else
        {
            cout<<endl<<"DIMENSION MISMATCH"<<endl;
            return Matrix<T>();
        }
    }
    template<class T>
    T add(const T & v1, const T & v2)
    {
        return v1 + v2;
    }
    template<class T>
    T subtract(const T & v1, const T & v2)
    {
        return v1 - v2;
    }
    template<class T>
    T multiply(const T & v1, const T & v2)
    {
        return v1 * v2;
    }
    template<class T>
    T divide(const T & v1, const T & v2)
    {
        return v1 / v2;
    }
        

    template<class T>
    Matrix<T> transpose(const Matrix<T> & m)
    {
        int r = m.dim(0);
        int c = m.dim(1);
        Matrix<T> n_m(c,r);

        for(int i=0; i<c; i++)
            for(int j=0; j<r; j++)
                n_m.a(i,j)=m.a(j,i);

        return n_m;
    }
    
    template<class T>
    Matrix<T> MatrixProduct(const Matrix<T> & m1, const Matrix<T> & m2)
    {
        int C=0; int A=0; int B=0;
        int r1 = m1.dim(0);  int r2 = m2.dim(0); int c1 = m1.dim(1);int c2 = m2.dim(1);

        if(c1 == r2)
        {
            C = c1;
            A = r1;
            B = c2;
        }
        else
        {
            cout<<endl<<"INNER DIMENSION INNEQUALITY"<<endl;
            return Matrix<T>();
        }

        Matrix<T> n_m(A,B);

        //perform multiplication
        for(int i=0; i<A; i++)
        {
            for(int k=0; k<C; k++)
            {
                for(int j=0; j<B; j++)
                {
                    n_m.a(i,j) +=  m1.e(i,k) * m2.e(k,j);
                }
            }
        }

        return n_m;
    }
 
    template<class T>
    Matrix<T> makeRowsUniqueOnColumns(const Matrix<T> & m, vector<int> c)
    {
        Matrix<T> retVal;
        
        map<vector<T>, char> uniqueChecker;
        
        for(int i=0; i<m.dim(0); ++i)
        {
            
            //Generate Key
            vector<T> tempKey;
            for(int j=0; j<c.size(); ++j)
                tempKey.push_back(m.a(i,c[j]));
            
           // if(uniqueChecker.count(m.e(i, c)) == 0)
            if(uniqueChecker.count(tempKey) == 0)
            {
                retVal.addRow(m.getRow(i));
                //uniqueChecker.insert(make_pair(m.e(i, c), 1));
                uniqueChecker.insert(make_pair(tempKey, 1));
                
            }
        }
        return retVal;
    }
    
    template<class T>
    Matrix<T> removeRowsWith(const Matrix<T> & m, const vector<int> & c, const vector<T> & v, bool (*foo1)(const bool, const bool), bool (*foo2)(const T &, const T &))
    {
        if(c.size() != v.size())
        {
            cout<<"Remove Row Vector Size Missmatch"<<endl;
            return Matrix<T>();
        }
        
        
        Matrix<T> retVal;
        
        //Over each row in the matrix
        for(int i=0; i<m.dim(0); ++i)
        {
            bool rowBool = foo2(m.a(i, c[0]), v[0]);
             
            for(int j=1; j<c.size(); ++j)
            {
                rowBool = foo1(rowBool, foo2(m.a(i, c[j]), v[j]));
            }
            
            if(!rowBool)
                retVal.addRow(m.getRow(i));
            
        }
        return retVal;
    }
    

    
    template<class T>
    static bool EQUALS(const T & t1, const T & t2)
    {
        return t1 == t2;
    }
    static bool AND(const bool b1, const bool b2)
    {
        return b1 && b2;
    }
    static bool OR(const bool b1, const bool b2)
    {
        return b1 || b2;
    }

    template<class T>
    static pair<Matrix<T>, Matrix<T> > permuteColumns(const Matrix<T> & m1, const Matrix<T> & m2 )
    {

        vector<int> ordering(m1.dim(1) + m2.dim(1),0);
        //ordering.resize(m1.dim(1) + m2.dim(1));
        for(int i=0; i< ordering.size(); ++i)
            ordering[i] = i;

        random_shuffle(ordering.begin(), ordering.end());
        
        pair<Matrix<T>, Matrix<T> > retVal;
  
        for(int i=0; i<ordering.size(); ++i)
        {
            if(i<m1.dim(1))
            {
                if(ordering[i]<m1.dim(1))
                {
                    retVal.first.addColumn(m1.getColumn(ordering[i]));
                }
                else
                {
                    retVal.first.addColumn(m2.getColumn(ordering[i]-m1.dim(1)));
                }
                
            }
            else
            {
                if(ordering[i]<m1.dim(1))
                {
                    retVal.second.addColumn(m1.getColumn(ordering[i]));
                }
                else
                {
                    retVal.second.addColumn(m2.getColumn(ordering[i]-m1.dim(1)));
                }
            }
        }

        cout<<retVal.first.dim(0)<<" "<<retVal.first.dim(1)<<" "<<retVal.second.dim(0)<<" "<<retVal.second.dim(1)<<endl;
        return retVal;
    }
}




#endif //MATRIX_MATH