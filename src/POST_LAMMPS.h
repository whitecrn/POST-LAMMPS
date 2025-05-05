#ifndef POST_LAMMPS_H_
#define POST_LAMMPS_H_
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <sstream>
#include <algorithm>

using namespace std;

class LAMMPS_DUMP
{
private:
struct Atom
{
    long step;long id;int type;double x;double y;double z;
    Atom(long step_=0,long id_=0,int type_=0,double x_=0,double y_=0,double z_=0):
        step(step_),id(id_),type(type_),x(x_),y(y_),z(z_){}
    
    bool operator<(const Atom& atom)
    {
        if (step!=atom.step) return step<atom.step;
        return id<atom.id;
    }
};

vector<vector<long>> atom_number; //atom number
vector<vector<Atom>> frames; // informations of each step
vector<vector<double>> boundary; // 0-1:xlo,xhi  2-3:ylo,yhi  4-5:zlo,zhi  6-8:xy,xz,yz

string file_name; //dump file name.
long step;
double potim;

public:
    // about constructor:
    LAMMPS_DUMP(string str="dump.atom",double p=0.1); //constructor
    ~LAMMPS_DUMP(); //destructor

    // about overload:
    LAMMPS_DUMP(const LAMMPS_DUMP &dump); //copy constructor
    LAMMPS_DUMP & operator=(const LAMMPS_DUMP &);

    // about information:
    long show_atom_number(long s,long type) const {return atom_number[s][type];} //return atom number of i step
    long show_step() const {return step;} //return total step;
    const string & show_file_name() const {return file_name;} //return file_name;

    // change private:
    void dt(double dt) {potim=dt;}

    // about calculations:
    double min_image(double,double) const;
    vector<pair<double,double>> msd(int) const;
    vector<pair<double,double>> rdf(int ct,int ne,int pix=100,int interval=1,double cutoff_r=8.0) const;
};

#endif