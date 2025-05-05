#include "POST_LAMMPS.h"

using namespace std;

LAMMPS_DUMP::LAMMPS_DUMP(string str,double p):file_name(str),step(0),potim(p)
{
    ifstream file;
    string line;
    file.open(file_name);

    if (!file.is_open())
    {
        cout << "Can't open " << file_name << " . " << endl;
        exit(1); 
    }

    while(getline(file,line))
    {
        vector<double> x_b(9,0.0);
        vector<double> a(6,0.0);
        vector<long> current_atom_number(1);
        vector<Atom> current_frame;
        for (int i=0;i<2;i++) getline(file,line);
        {long N;file >> N;current_atom_number[0]=(N);}
        getline(file,line);
        getline(file,line);
        for (int i=0;i<3;i++) 
        {
            getline(file,line);
            stringstream ss(line);
            ss >> x_b[i*2] >> x_b[i*2+1];
            a[i]=x_b[i*2+1]-x_b[i*2];
        }

        
        for (int i=0;i<current_atom_number[0];i++)
        {
            Atom atom;
            getline(file,line);
            stringstream ss(line);
            ss >> atom.id >> atom.type >> atom.x >> atom.y >> atom.z;
            atom.x-=x_b[0];
            atom.y-=x_b[2];
            atom.z-=x_b[4];
            if (atom.type>=static_cast<long>(current_atom_number.size())) current_atom_number.resize(atom.type+1);
            atom.step=step;
            current_atom_number[atom.type]+=1;
            current_frame.push_back(atom);
        }
        sort(current_frame.begin(),current_frame.end());
        frames.push_back(current_frame);
        atom_number.push_back(current_atom_number);
        boundary.push_back(a);
        step++;
    }
}

LAMMPS_DUMP::~LAMMPS_DUMP()
{
    atom_number.clear();
    frames.clear();
    boundary.clear();
    step=0;
    file_name='\0';
}

LAMMPS_DUMP::LAMMPS_DUMP(const LAMMPS_DUMP &dump)
{
    file_name=dump.file_name;
    atom_number=dump.atom_number;
    frames=dump.frames;
    boundary=dump.boundary;
    step=dump.step;
    potim=dump.potim;
}

LAMMPS_DUMP & LAMMPS_DUMP::operator=(const LAMMPS_DUMP &dump)
{
    if (this==&dump) return *this;

    atom_number.clear();
    frames.clear();
    boundary.clear();
    step=0;
    file_name='\0';

    file_name=dump.file_name;
    atom_number=dump.atom_number;
    frames=dump.frames;
    boundary=dump.boundary;
    step=dump.step;
    potim=dump.potim;

    return *this;
}


vector<pair<double,double>> LAMMPS_DUMP::msd(int select_type) const
{
    long real_step=static_cast<long>(round(static_cast<double>(step)/2));
    vector<pair<double,double>> msd_data(real_step,make_pair(0.0,0.0)); //initialize msd
    vector<double> total_msd(real_step,0.0);
    for (int delta_t=0;delta_t<real_step;delta_t++)
    {
        double n=0;
        for (int t=0;t<real_step-1;t++)
        {
            if (select_type>static_cast<int>(atom_number[t].size()) || atom_number[t][select_type]<=0)
            {
                cout << "Can not find type: " << select_type << " or type number is 0 or less!" << endl;
                exit(1);
            }

            for (int atom=0;atom<atom_number[t][0];atom++)
                if (frames[t][atom].type==select_type)
                    {
                        double dx=frames[t+delta_t][atom].x-frames[t][atom].x;dx=min_image(dx,boundary[t][0]);
                        double dy=frames[t+delta_t][atom].y-frames[t][atom].y;dy=min_image(dy,boundary[t][1]);
                        double dz=frames[t+delta_t][atom].z-frames[t][atom].z;dz=min_image(dz,boundary[t][2]);
                        double dr=(dx*dx+dy*dy+dz*dz);
                        total_msd[delta_t]+=dr;
                        n++;
                    }
                else continue;
        }
        msd_data[delta_t].first=delta_t*potim; 
        msd_data[delta_t].second=total_msd[delta_t]/n;
    }
    return msd_data;
}

vector<pair<double,double>> LAMMPS_DUMP::rdf(int ct,int ne,int pix,int interval,double cutoff_r) const
{
    const double PI=3.1415926535;
    int size_rdf=static_cast<int>(round(pix*cutoff_r));
    double bin=1/static_cast<double>(pix);
    vector<pair<double,double>> rdf(size_rdf,make_pair(0.0,0.0));
    vector<double> total_rdf(size_rdf,0.0);
    if (interval>step)
    {
        cout << "Error: unvaild interval!" << endl;
        exit(1);
    }

    for (long s=0;s<step;s+=interval)
    {
        if (s>=static_cast<long>(frames.size()) || frames[s].empty()) continue; 
        for (int center=0;center<atom_number[s][0];center++)
        {
            for (int neighbor=0;neighbor<atom_number[s][0];neighbor++)
            {
                if (frames[s][center].type==ct && frames[s][neighbor].type==ne)
                {
                    double dx=frames[s][center].x-frames[s][neighbor].x;dx=min_image(dx,boundary[s][0]);
                    double dy=frames[s][center].y-frames[s][neighbor].y;dy=min_image(dy,boundary[s][1]);
                    double dz=frames[s][center].z-frames[s][neighbor].z;dz=min_image(dz,boundary[s][2]);
                    double dr=sqrt(dx*dx+dy*dy+dz*dz);
                    if (dr<cutoff_r)
                    {
                        int r_group=static_cast<int>(round(dr/bin));
                        total_rdf[r_group]+=1.0;
                    }
                }
            }
        }
    }
    for (int r=1;r<size_rdf;r++)
    {
        double shell_volume=4.0*PI*r*r*bin;
        rdf[r].first=r*bin;
        if (atom_number[0][ne]>0) rdf[r].second=total_rdf[r]/(shell_volume*atom_number[0][ct]*atom_number[0][ne]*step/interval);
    }
    return rdf;
}

double LAMMPS_DUMP::min_image(double a,double boundary) const
{
    if (boundary==0.0) return a;
    return a-boundary*round(a/boundary);
}