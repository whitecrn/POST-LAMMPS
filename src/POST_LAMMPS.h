#ifndef POST_LAMMPS_H_
#define POST_LAMMPS_H_
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <stdexcept>

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

    vector<vector<long>> atom_number; // atom number
    vector<vector<Atom>> frames; // informations of each step
    vector<vector<double>> boundary; // 0:a 1:b 2:c 3-5:xy,xz,yz

    string file_name; //dump file name.
    long step;
    double potim;

public:
    // constructor:
    LAMMPS_DUMP(string str="dump.atom",double p=0.1); //constructor
    ~LAMMPS_DUMP(); //destructor

    // overload:
    LAMMPS_DUMP(const LAMMPS_DUMP &dump); //copy constructor
    LAMMPS_DUMP & operator=(const LAMMPS_DUMP &); // = overload

    // information:
    long show_atom_number(long s,long type) const {return atom_number[s][type];} //return atom number of i step
    long show_step() const {return step;} //return total step;
    const string & show_file_name() const {return file_name;} //return file_name;

    // change private:
    void dt(double dt) {potim=dt;}

    // about calculations:
    double Min_image(double,double) const;
    vector<pair<double,double>> Msd(int) const;
    vector<pair<double,double>> Rdf(int ct=1,int ne=1,int pix=100,int interval=1,double cutoff_r=8.0) const;
    vector<pair<double,double>> Selected_Rdf(int ct=1,int ne=1,int pix=100,int interval=1,double cutoff_r=8.0,double z_min=0.0,double z_max=0.0,double volume_z_min=0.0,double volume_z_max=0.0) const;
    vector<vector<double>> Velocity(double mesh,int delta_t) const;
    vector<double> Flux(int type,double mesh) const;
    vector<double> Bond_Lenth(int ct=1,int ne=1,int pix=100,int interval=1,double cutoff_r=2.5) const;
    vector<double> Coordination(int ct=1,int ne=1,int pix=100,int interval=1,double cutoff_r=2.5) const;
    vector<long> Discharge_Speed(int type=1,double min_z=0.0) const;
};

//class LAMMPS_DATA
//{
//private:

//struct Atom
//{
//    long id;int type;double x;double y;double z;
//    Atom(long id_=0,int type_=0,double x_=0,double y_=0,double z_=0):
//        id(id_),type(type_),x(x_),y(y_),z(z_){}
    
//    bool operator<(const Atom& atom)
//    {
//        if (id!=atom.id) return id<atom.id;
//    }
//};

//vector<long> atom_number(1,0);
//vector<string> atom_type;
//vector<Atom> R;
//vector<double> boundary;

//string file_name;

//public:
    
    // constructor
    //LAMMPS_DATA(string str="lammps.data");
    //~LAMMPS_DATA();

    // about overload:
    //LAMMPS_DATA(const LAMMPS_DATA &data); //copy constructor
    //LAMMPS_DATA & operator=(const LAMMPS_DATA &); // = overload

    // information:
    //long show_atom_number(long type) const {return atom_number[type];} //return atom number of i step
    //const string & show_file_name() const {return file_name;} //return file_name;

    // crystal operation:
    //LAMMPS_DATA & Construct_interface(const LAMMPS_DATA& other_interface,int direction=3,double vacuum=20.0,double interface_vacuum=2.5) const;
    //LAMMPS_DATA & Replicate(int x1=1,int x2=1,int x3=1) const;
    //LAMMPS_DATA & Period(int dirction=3) const;

    // calculation
    //vector<pair<double,double>> & Rdf(int ct=1,int ne=1,int pix=100,double cutoff_r=8.0) const;
    //vector<pair<double,double>> & Xrd(double ray,double start=0.0,double end=90.0) const;
//};
#endif