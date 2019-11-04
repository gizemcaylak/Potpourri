#include "FilePath.h"

FilePath::FilePath()
{
    char CREATION_PATH[FILENAME_MAX];
    if(getcwd(CREATION_PATH, sizeof(CREATION_PATH))!=NULL)
    {
        path_ = string(CREATION_PATH);
        structFromString(path_);
    }
    else
    {
        cout<<"FILEPATH ERROR"<<endl;
    }
}

FilePath::FilePath(const FilePath & cpy)
{

    strncpy(PRE_ACTIVE_DIRECTORY, cpy.PRE_ACTIVE_DIRECTORY, sizeof(cpy.PRE_ACTIVE_DIRECTORY));

    path_ = cpy.path_;
    directoryBreaks_ = cpy.directoryBreaks_;
    depth_ = cpy.depth_;
}

void FilePath::setFilePath(string filePath)
{
   
    if(getcwd(PRE_ACTIVE_DIRECTORY, sizeof(PRE_ACTIVE_DIRECTORY)) != NULL)
    {
        //if(SetCurrentDirectory(path_in.c_str()))
        if(chdir(filePath.c_str())==0)
            structFromString(filePath);
        else
            cout<<"FILE PATH INVALID "<<filePath<<endl;

        if(!chdir(PRE_ACTIVE_DIRECTORY)==0)
            cout<<"RETURN DIRECTORY INVALID"<<endl;
    }
}

string FilePath::getFilePath()const
{
//cout<<string(path_.begin(), path_.begin()+directoryBreaks_[depth_])<<endl;
    return string(path_.begin(), path_.begin()+directoryBreaks_[depth_]);
}

void FilePath::moveUp()
{
    if(depth_>=1)
        depth_ = depth_-1;
}

void FilePath::moveDown()
{
    if(depth_<(int)directoryBreaks_.size()-1)
        depth_ = depth_+1;
}

void FilePath::moveDown(string new_directory)
{

    makeActive();
    //CreateDirectory(new_directory.c_str(), NULL);
    mkdir(new_directory.c_str(), 0777);

    //remove path to old subdirectories

    if(depth_ < (int)directoryBreaks_.size()-1)
    {
        path_.erase(path_.begin()+directoryBreaks_[depth_], path_.end());

        int delete_last = directoryBreaks_.size()-depth_;
        for(int i=0; i<delete_last; ++i)
            {directoryBreaks_.pop_back();}
    }



    //append new directory
    path_.append("/" + new_directory);
    directoryBreaks_.push_back(path_.size());
    depth_=directoryBreaks_.size()-1;


    returnActive();

}


void FilePath::makeActive()
{
    if(getcwd(PRE_ACTIVE_DIRECTORY, sizeof(PRE_ACTIVE_DIRECTORY)) != NULL)
        if(!chdir(getFilePath().c_str())==0)
            cout<<"MAKE ACTIVE DIRECTORY INVALID"<<endl;
}

void FilePath::returnActive()
{

    if(!chdir(PRE_ACTIVE_DIRECTORY)==0)
        cout<<"RETURN DIRECTORY INVALID"<<endl;
}

void FilePath::structFromString(string path_in)
{
    //int t_begin=0;
    for(unsigned int i=0; i<path_in.size(); ++i)
    //Find the path segments
    if(path_in[i] == '/')
    {
        directoryBreaks_.push_back(i);
    }
    directoryBreaks_.push_back(path_in.size());
    depth_=directoryBreaks_.size()-1;
}
