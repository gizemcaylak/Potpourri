/**
 * @author Tyler Cowman
 * 
 * This class adds some additional functionality for keeping track of
 * and changing the current working directory. To enable compiling for 
 * a Windows machine, this "should" be the only file that needs modification 
 * by adding a flag to swithc between the Linux system calls and Windows.
 */


#ifndef FILEPATH_H
#define FILEPATH_H

#include <string.h>
#include <iostream>
#include <stdio.h>
//#include<windows.h>
#include <unistd.h>
#include <sys/stat.h>

#include<vector>

using namespace std;
class FilePath
{
    public:
        FilePath();
        FilePath(const FilePath & cpy);

        void setFilePath(string filePath);
        string getFilePath()const;


        //navigation
        void moveUp();
        void moveDown();
        void moveDown(string newDirectory);

        void makeActive();
        void returnActive();

    protected:
    private:

        void structFromString(string filePath);

        //Private variables
        char PRE_ACTIVE_DIRECTORY[FILENAME_MAX];

        string path_;
        vector<int> directoryBreaks_;
        int depth_;


};

#endif // FILEPATH_H
