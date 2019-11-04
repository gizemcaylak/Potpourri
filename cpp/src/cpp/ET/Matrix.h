#ifndef MATRIX_H
#define MATRIX_H

#include<fstream>
#include<iostream>
#include<vector>
#include<algorithm>
#include<string>
#include<sstream>
#include<utility>
#include<omp.h>

using namespace std;

template<class T>
class Matrix;

template<class T>
ostream& operator<< (ostream &out, Matrix<T> & M);

/** @class cmp
 *  @brief Simple class for use in order by column
 */
template<class T>
class cmp 
{
    int byRow;
    bool asNumeric;
    bool ascending;
public:
    cmp(int r, bool n, bool a) : byRow(r), asNumeric(n), ascending(a) {}

    bool operator()(const vector<T> & v1 , const vector<T> & v2)
    {
        
        
        switch(asNumeric)
        {
            case true:
                switch(ascending)
                {
                    case true:
                        return stof(v1[byRow]) < stof(v2[byRow]);
                        break;
                    case false:
                        return stof(v1[byRow]) > stof(v2[byRow]);
                        break;
                }
                break;
            case false:
                switch(ascending)
                {
                    case true:
                        return v1[byRow]< v2[byRow];
                        break;
                    case false:
                        return v1[byRow]> v2[byRow];
                        break;
                }
                break;
        }
            
    }
};



/** @class Matrix
 *  @brief Two-dimensional templetized data structure. (y, x) (Rows, Columns).
 *
 *  Implemented as a vector of vectors of the type T.
 */
template<class T>
class Matrix
{
    public:
        /** @brief Constructs an empty Matrix.
         *
         *  Calls the empty constructors for a vector of vectors for the
         *  type T. Sets the columns and rows values to 0.
         */
        Matrix();

        /** @brief Constructs a rows by columns Matrix
         *
         *  Calls the default constructors for vectors of vectors.
         *  of the type T. Sets the columns and rows values accordingly.
         *  @param rows_in The initial number of rows to build.
         *  @param columns_in The initial number of columns to build.
         */
        Matrix(unsigned int rows_in, unsigned int columns_in);

        /** @brief Constructs a copy of a Matrix.
         *
         *  @param r A reference to the Matrix being copied.
         */
        Matrix(const Matrix & r);

        /** @brief Constructs a numerical Matrix from a formatted string.
         *
         *  Only use this to build numerical (decimal values work) matrices. Input should
         *  be in the form "# # ; # #" where # is a number and ; denotes a new row.
         *  @param num_string The formatted string input.
         *
         */
        Matrix(string num_string);

        /** @brief Destructor for a Matrix object
         *  @see clear
         */
        ~Matrix();

        /** @brief Clears a Matrix
         *
         *  Deletes all data in the Matrix and sets t
         *  rows and columns variables to 0.
         */
        void clear();

        
        friend ostream& operator<< <> (ostream &out, Matrix<T> & M);
        
        
        /** @brief Load a Matrix from a file
         *
         *  The file should be in a csv format with spaces as separators and new lines
         *  denoting a new row. This function can read numbers strings and chars, but
         *  the data type must be constant and it must be a rectangular Matrix.
         *
         *  First clears the current Matrix then reads the file and sets the size.
         * @param filename The name and path of the file to load from
         */
        bool loadFromFile(string filename, int rows);

        /** @brief Adds a single row to the Matrix.
         *
         *  Increases the row size of the Matrix by one, by adding a row of
         *  size columns default constructed types T. the row is added to the bottom.
         */
        void addRow();

        /** @brief Adds a single column to the Matrix.
         *
         *  Increases the column size of the Matrix by one, by adding a column of
         *  size rows default constructed types T. The column is added to the right.
         */
        void addColumn();

        /** @brief Adds a row of specified values to the Matrix.
         *
         *  Increases the row size of the Matrix by adding values specified by a vector.
         *  The vector must be equal in size to the rows and of the same type T.
         *  @param new_row The vector of type T to add to the Matrix.
         */
        void addRow(vector<T> new_row);

        /** @brief Adds a column of specified values to the Matrix.
         *
         *  Increases the column size of the Matrix by adding values specified by a vector.
         *  The vector must be equal in size to the columns and of the same type T.
         *  @param new_row The vector of type T to add to the Matrix.
         */
        void addColumn(vector<T> new_column);

        void setColumnLabels(vector<string> V);
        const vector<string> & getColumnLabels()const;
        
        /*access*/
        Matrix<T> get_sub(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)const;
        const vector<T> & getRow(unsigned int r)const;
        vector<T> & getRow(unsigned int r);
        vector<T> getRow(unsigned int r, vector<int> columns)const; 
        vector<T> getColumn(unsigned int c)const;
        const T & a(int y, int x)const;
        T & a(int y, int x);
        
        bool replaceRow(int r, const vector<T> & newRow);
        
        /*logic and info*/
        bool empty()const;
        int size()const;
        int dim(bool dim_choice)const;        

        Matrix orderRows(vector<int> v)const;
        Matrix mat_sort_row()const;
        Matrix mat_sort_column()const;
        
        //void sortByColumn( int column, const cmp & (*foo)(int)  );
        void sortByColumn(int column, bool asNumeric, bool ascending);
        //vector<int> sortByColumn(int column, bool asNumeric, bool ascending);
        
        bool sortHelper(const vector<T> & v1, const vector<T> & v2/*vector<T>::iterator it1, vector<T>::iterator it2*/);
        void sortt(int column);

        /*printing*/
        void print()const;
        void printToFile(ofstream & ofs, string delimiter)const;

    private:
        vector<vector<T> > d_; /**< Underlying vector of vectors storage structure for type T. */
        unsigned int rows_; /**< Holds the current number of rows in an unsigned int. */
        unsigned int columns_; /**< Holds the current number of columns in an unsigned int. */
        
        vector<string> columnLabels_;
};



/*IMPLEMENTATION---------------------------------------------------*/

/*begin constructor destructor-------------------------------------------*/
template<class T>
Matrix<T>::Matrix()
{
    d_ = vector<vector<T> >(0, vector<T>(0, T()));
    rows_=0;
    columns_=0;
}

template<class T>
Matrix<T>::Matrix(unsigned int rows_in, unsigned int columns_in)
{
    //d = vector<vector<T> >(rows_in, vector<T>(columns_in, T()));
    d_ = vector<vector<T> >(rows_in, vector<T>());

    //#pragma omp parallel for
    for(unsigned int i=0; i<rows_in; i++)
        d_[i] = vector<T>(columns_in, T());

    rows_ = rows_in;
    columns_ = columns_in;
}

template<class T>
Matrix<T>::Matrix(const Matrix & r)
{
    d_ = vector<vector<T> >(r.rows_, vector<T>());

    //#pragma omp parallel for
    for(unsigned int i=0; i<r.rows_; i++)
        d_[i] = vector<T>(r.columns_, T());

    for(unsigned int x=0; x<d_.size(); x++)
        for(unsigned int y=0; y<d_[0].size(); y++)
            d_[x][y] = r.d_[x][y];

    rows_ = r.rows_;
    columns_ = r.columns_;
}

template<class T>
Matrix<T>::Matrix(string num_string)
{
    stringstream ss;
    ss<<num_string;

    rows_ = 0;
    columns_ = 0;


    string temp_v;
    vector<T> current_row;

    while(!ss.eof())
    {
        ss>>temp_v;
        if(temp_v == ";")
        {
            if(!current_row.empty())
                addRow(current_row);

            current_row.clear();
        }
        else
        {
            stringstream ss2;
            ss2<<temp_v;
            T num;
            ss2>>num;

            current_row.push_back(num);
        }
    }
    if(!current_row.empty())
        addRow(current_row);
}

template<class T>
Matrix<T>::~Matrix()
{
    clear();
}

template<class T>
void Matrix<T>::clear()
{
    for(unsigned int i=0; i<d_.size(); i++)
        d_[i].clear();

    d_.clear();
    rows_ = 0;
    columns_ = 0;
    
    columnLabels_ = vector<string>();
}
/*end constructor destructor---------------------------------------*/

template<class T>
ostream& operator<< (ostream &out, Matrix<T> & M)
{
    for(int x=0; x<M.d_.size(); x++)
    {
        for(int y=0; y<M.d_[0].size(); y++)
        {
            out<<M.d_[x][y]<<'\t';
        }
        out<<endl;
    } 
    
     return(out);
}

/*begin file loading-----------------------------------------------*/
template<class T>
bool Matrix<T>::loadFromFile(string filename, int rows)
{
    ifstream FILE(filename.c_str());
    if(!FILE.is_open())
    {
        cout<<"FILE NOT OPENED"<<endl;
        return false;
    }

    clear();

    //Read all rows
    if(rows_ == -1)    
    {
        while(!FILE.eof())
        {
            stringstream ss;
            vector<T> temp_vec;
            string temp_line;

            getline(FILE, temp_line);
            ss<<temp_line;

            while(!ss.eof())
            {
                T temp_d;
                ss>>temp_d;
               // if(ss.eof())
               //     break;

                temp_vec.push_back(temp_d);
            }
           // if(FILE.eof())
           //     return true;
            addRow(temp_vec);
        }
    }
    else
    {
        for(int i=0; i<rows_; ++i)
        {
            stringstream ss;
            vector<T> temp_vec;
            string temp_line;

            getline(FILE, temp_line);
            ss<<temp_line;

            while(!ss.eof())
            {
                T temp_d;
                ss>>temp_d;
                //if(ss.eof())
                  //  break;

                temp_vec.push_back(temp_d);
            }
            if(FILE.eof())
                return true;
            addRow(temp_vec);
        }
    }

    FILE.close();
    return true;
}
/*end file loading-------------------------------------------------*/

/*begin dimension manipulation-------------------------------------*/
template<class T>
void Matrix<T>::addRow()
{
    d_.push_back(vector<T>(columns_));
    ++rows_;
}

template<class T>
void Matrix<T>::addColumn()
{
    for(unsigned int i=0; i<rows_; i++)
        d_[i].push_back(T());
    ++columns_;
}

template<class T>
void  Matrix<T>::addRow(vector<T> new_row)
{
    if(empty())
    {
        columns_ = new_row.size();
        rows_ = 1;

        d_.push_back(vector<T>(new_row.size(), T()));
        for(unsigned int i=0; i<new_row.size(); i++)
            d_[0][i]=new_row[i];
    }
    else if(new_row.size() == 0)
    {
        return;

    }
    else if(new_row.size() != columns_)
    {
        cout<<"NEW ROW DIMENSION MISMATCH "<<rows_<<" SIZE "<<new_row.size()<<" TO "<<dim(1)<<endl;
        return;
    }
    else
    {
        addRow();
        for(unsigned int i=0; i<new_row.size(); i++)
            d_[rows_-1][i]=new_row[i];
    }
}

template<class T>
void  Matrix<T>::addColumn(vector<T> new_column)
{
    if(empty())
    {
        rows_ = new_column.size();
        columns_ = 1;

        for(unsigned int i=0; i<new_column.size(); i++)
        {
            d_.push_back(vector<T>(1,T()));
            d_[i][0] = new_column[i];
        }
    }
    else if(new_column.size() != rows_)
    {
        cout<<"NEW COLUMN DIMENSION MISMATCH"<<endl;
        return;
    }
    else
    {
        addColumn();
        for(unsigned int i=0; i<new_column.size(); i++)
            d_[i][columns_-1]=new_column[i];
    }

}
/*end dimension manipulation---------------------------------------*/
template<class T>
void Matrix<T>::setColumnLabels(vector<string> V)
{
    columnLabels_ = V;   
}

template<class T>
const vector<string> & Matrix<T>::getColumnLabels()const
{
    return columnLabels_;
}

/*begin access-----------------------------------------------------*/
template<class T>
Matrix<T> Matrix<T>::get_sub(unsigned int y1, unsigned int y2, unsigned int x1, unsigned int x2)const
{
    int n_r=y2-y1+1;
    int n_c=x2-x1+1;
    if(y2>rows_-1 || x2>columns_-1 || y1<0 || x1<0 || n_r<0 || n_c<0)
    {
        cout<<"SUB_MATRIX DIMENSIONS INVALID"<<endl;
        return Matrix<T>();
    }
    Matrix<T> n_m(n_r, n_c);

    //extract sub Matrix
    for(unsigned int j=0; j<n_r; j++)
    {
         for(unsigned int i=0; i<n_c; i++)
         {
                n_m.d_[j][i] = d_[y1+j][x1+i];
         }
    }

    return n_m;
}

template<class T>
const vector<T> & Matrix<T>::getRow(unsigned int r)const
{
    return d_.at(r);
}

template<class T>
vector<T> & Matrix<T>::getRow(unsigned int r)
{
    return d_[r];
}

template<class T>
vector<T> Matrix<T>::getRow(unsigned int r, vector<int> columns)const
{
    vector<T> retVal;
    for(int i=0; i<columns.size(); ++i)
        retVal.push_back(a(r,columns[i]));
    
    return retVal;
}

template<class T>
vector<T> Matrix<T>::getColumn(unsigned int c)const
{
    vector<T> n_m(rows_, T());

    for(int i=0; i<rows_; i++)
        n_m[i]=d_[i][c];

    return n_m;
}

/*end access-------------------------------------------------------*/

/*begin direct manipulation----------------------------------------*/
template<class T>
T & Matrix<T>::a(int y, int x)
{
    return d_[y][x];
}

template<class T>
const T & Matrix<T>::a(int y, int x)const
{
    return d_[y][x];
}

/*end direct manipulation------------------------------------------*/

template<class T>
bool Matrix<T>::replaceRow(int r, const vector<T> & newRow)
{
    if(newRow.size() != dim(1))
    {
        cout<<"ROW SIZE MISMATCH "<<dim(1)<<" "<<newRow.size()<<endl;
        return false;
    }
    for(int i=0; i<newRow.size(); ++i)
        a(r, i) = newRow[i];
    
}


/*begin logic and info---------------------------------------------*/
template<class T>
bool Matrix<T>::empty()const
{
    if(rows_==0 && columns_==0)
        return true;
    else
        return false;
}

template<class T>
int Matrix<T>::size()const
{
    return columns_*rows_;
}

template<class T>
int Matrix<T>::dim(bool dim_choice)const
{
    if(dim_choice==0)
        return rows_;
    else
        return columns_;
}

/*end logic and info-----------------------------------------------*/

/*begin math-------------------------------------------------------*/




template<class T>
Matrix<T> Matrix<T>::orderRows(vector<int> v)const
{
    Matrix<T> retVal;
//cout<<v.size()<<endl;
    for(int i=0; i<v.size(); ++i)
    {
        retVal.addRow((*this).getRow(v[i]));
    }


    return retVal;
}

template<class T>
Matrix<T> Matrix<T>::mat_sort_row()const
{
    Matrix<T> n_m(rows_,columns_);
    for(int i=0; i<rows_; i++)
        for(int j=0; j<columns_; j++)
            n_m.d_[i][j]=d_[i][j];

    sort(n_m.d_.begin(), n_m.d_.end());

    return n_m;
}

template<class T>
Matrix<T> Matrix<T>::mat_sort_column()const
{
    Matrix<T> n_m(rows_,columns_);
    for(int i=0; i<rows_; i++)
        for(int j=0; j<columns_; j++)
            n_m.d_[i][j]=d_[i][j];

    n_m.mat_transpose() = n_m.mat_transpose() ;
    sort(n_m.d_.begin(), n_m.d_.end());
    //sort(n_m.d.rbegin(), n_m.d.rend()); //reverse sort
    n_m.mat_transpose() = n_m.mat_transpose();

    return n_m;
}

/*
template<class T>
void Matrix<T>::sortByColumn(int column, const cmp & (*foo)(int) )
{
    cmp<T> c = foo;
    
    sort(d.begin(), d.end(), );//cmp<T>(column));//cmp<T>(column));
}
*/

template<class T>
void Matrix<T>::sortByColumn(int column, bool asNumeric, bool ascending)
{
    
    stable_sort(d_.begin(), d_.end(), cmp<T>(column, asNumeric, ascending));//cmp<T>(column));
}


/*end math---------------------------------------------------------*/

/*begin printing---------------------------------------------------*/
template<class T>
void Matrix<T>::print()const
{
    for(int x=0; x<d_.size(); x++)
    {
        for(int y=0; y<d_[0].size(); y++)
        {
            cout<<d_[x][y]<<" ";
        }
        cout<<endl;
    }
}

template<class T>
void Matrix<T>::printToFile(ofstream & ofs , string delimiter)const
{
    for(int x=0; x<d_.size(); x++)
    {
        for(int y=0; y<d_[0].size(); y++)
        {
            ofs<<d_[x][y]<<delimiter;
        }
        ofs<<endl;
    }
}

/*end printing-----------------------------------------------------*/

/*begin helpers----------------------------------------------------*/
/*end helpers------------------------------------------------------*/

#endif // Matrix_H
