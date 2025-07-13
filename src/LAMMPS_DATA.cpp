#include "POST_LAMMPS.h"

using namespace std;

LAMMPS_DATA::LAMMPS_DATA(string str):file_name(str)
{
    ifstream file;
    string line;
    file.open(file_name);

    if (!file.is_open())
    {
        cout << "Can't open " << file_name << ". " << endl;
        exit(1);
    }

    while(getline(file,line))
    {
        vector<double> x_b(9,0.0);
        for (int i=0;i<2;i++) getline(file,line);
        vector
        {long N;file >> N;atom_number[0]=N;}
        getline(file,line);getline(file,line);
        for (int i=0;i<3;i++)
        {
            getline(file,line);
            stringstream ss(line);
            ss >> x_b[i*2] >> x_b[i*2+1];
            boundary[i]=x_b[i*2+1]-x_b[i*2];
        }

        getline(file,line);getline(file,line);getline(file,line);

        for (int i=1;;i++)
        {
            getline(file,line);
            stringstream ss(line);
            int id;double mass;string temp;string element;
            ss >> id >> mass >> temp >> element;
            if (id!=i)
            {
                cout << "atom id error, exit." << endl;
                exit(1);
            }
            
        }
    }
}