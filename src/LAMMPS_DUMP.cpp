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
        getline(file,line);getline(file,line);
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


vector<pair<double,double>> LAMMPS_DUMP::Msd(int select_type) const
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
                        double dx=frames[t+delta_t][atom].x-frames[t][atom].x;dx=Min_image(dx,boundary[t][0]);
                        double dy=frames[t+delta_t][atom].y-frames[t][atom].y;dy=Min_image(dy,boundary[t][1]);
                        double dz=frames[t+delta_t][atom].z-frames[t][atom].z;dz=Min_image(dz,boundary[t][2]);
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

vector<pair<double,double>> LAMMPS_DUMP::Rdf(int ct,int ne,int pix,int interval,double cutoff_r) const
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
                if (frames[s][center].type==ct && frames[s][neighbor].type==ne && center!=neighbor)
                {
                    double dx=frames[s][center].x-frames[s][neighbor].x;dx=Min_image(dx,boundary[s][0]);
                    double dy=frames[s][center].y-frames[s][neighbor].y;dy=Min_image(dy,boundary[s][1]);
                    double dz=frames[s][center].z-frames[s][neighbor].z;dz=Min_image(dz,boundary[s][2]);
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
        double r_mid=(r+0.5)*bin;
        double shell_volume=4.0*PI*r_mid*r_mid*bin;
        double rho=atom_number[0][ne]/(boundary[0][0]*boundary[0][1]*boundary[0][2]);
        rdf[r].first=r*bin;
        if (atom_number[0][ne]>0) rdf[r].second=total_rdf[r]/(shell_volume*atom_number[0][ct]*(step/interval)*rho);
    }
    return rdf;
}

vector<pair<double,double>> LAMMPS_DUMP::Selected_Rdf(int ct,int ne,int pix,int interval,double cutoff_r,double z_min,double z_max,double volume_z_min,double volume_z_max) const
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

    int N_A=0;
    int N_B=0;
    for (long s=0;s<step;s+=interval)
    {
        for (int neighbor=0;neighbor<atom_number[s][0];neighbor++)
        {
            if (frames[s][neighbor].type==ne && frames[s][neighbor].z<=volume_z_max && frames[s][neighbor].z>=volume_z_min) N_B++;
        }
    }

    for (long s=0;s<step;s+=interval)
    {
        if (s>=static_cast<long>(frames.size()) || frames[s].empty()) continue; 
        for (int center=0;center<atom_number[s][0];center++)
        {
            if (frames[s][center].z<=z_max && frames[s][center].z>=z_min)
            {
                N_A++;
                for (int neighbor=0;neighbor<atom_number[s][0];neighbor++)
                {
                    if (frames[s][center].type==ct && frames[s][neighbor].type==ne && center!=neighbor)
                    {
                        double dx=frames[s][center].x-frames[s][neighbor].x;dx=Min_image(dx,boundary[s][0]);
                        double dy=frames[s][center].y-frames[s][neighbor].y;dy=Min_image(dy,boundary[s][1]);
                        double dz=frames[s][center].z-frames[s][neighbor].z;dz=Min_image(dz,boundary[s][2]);
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
    }

    N_A=N_A/(step/interval);
    N_B=N_B/(step/interval);

    for (int r=1;r<size_rdf;r++)
    {
        double r_mid=(r+0.5)*bin;
        double shell_volume=4.0*PI*r_mid*r_mid*bin;
        double rho=N_B/(boundary[0][0]*boundary[0][1]*(volume_z_max-volume_z_min));
        rdf[r].first=r*bin;
        if (atom_number[0][ne]>0) rdf[r].second=total_rdf[r]/(shell_volume*N_A*(step/interval)*rho);
    }
    return rdf;  
}

vector<vector<double>> LAMMPS_DUMP::Velocity(double mesh,int delta_t) const
{
    long Nx=static_cast<long>(round(boundary[0][0]/mesh))+1;
    long Nz=static_cast<long>(round(boundary[0][2]/mesh))+1;

    vector<vector<double>> velocity(Nx,vector<double>(Nz,0.0));
    vector<vector<long>> n(Nx,vector<long>(Nz,0));
    for (int i=0;i<atom_number[0][0];i++)
    {
        if (frames[0][i].type==1)
        {
            for (int t=0;t<200;t++)
            {
                int x=static_cast<int>(round(frames[t][i].x/mesh));
                int z=static_cast<int>(round(frames[t][i].z/mesh));
                n[x][z]++;
                double d=frames[t+delta_t][i].z-frames[t][i].z;Min_image(d,boundary[t][2]);
                if (abs(d)>1.0) 
                {
                    velocity[x][z]+=d/(delta_t*potim);
                    n[x][z]++;
                }
            }
        }
        else continue;
    }

    for (int x=0;x<Nx;x++)
    {
        for (int z=0;z<Nz;z++) 
        {
            if (n[x][z]>0) velocity[x][z]/=n[x][z];
        }
    }
    return velocity;
}

vector<double> LAMMPS_DUMP::Flux(int type,double mesh) const
{
    int Nz=static_cast<int>(round(boundary[0][2]/mesh))+1;
    vector<double> flux(Nz,0.0);
    double distance=1.5;

    for (int i=0;i<atom_number[0][0];i++)
    {
        if (frames[0][i].type==type)
        {
            int z_start=static_cast<int>(round(frames[0][i].z/mesh));
            for (int t=1;t<step;t++)
            {
                int z_end=static_cast<int>(round(frames[t][i].z/mesh));
                int d=z_end-z_start;

                if (d>distance/mesh) 
                {
                    for (int p=z_start;p<z_end;p++) flux[p]+=1.0;
                    z_start=static_cast<int>(round(frames[t][i].z/mesh));
                }
                else if (d>0 && d<distance/mesh) continue;
                else if (d<0 && d>-distance/mesh) continue;
                else if (d<-distance/mesh) 
                {
                    for (int p=z_start;p<z_end;p--) flux[p]-=1.0;
                    z_start=static_cast<int>(round(frames[t][i].z/mesh));
                }                
            }
        }
    }

    for (int z=0;z<Nz;z++)
    {
        double area=boundary[0][0]*boundary[0][1];
        flux[z]=flux[z]/(area*step);
    }
    return flux;
}

vector<double> LAMMPS_DUMP::Bond_Lenth(int ct,int ne,int pix,int interval,double cutoff_r) const
{
    long Nz=static_cast<int>(round(boundary[0][2]*pix))+1;
    vector<double> bond_lenth(Nz,0.0);
    vector<long> n(Nz,0);
    cout << frames[0][2].x << endl;
    
    for (int i=0;i<step-interval;i+=interval)
    {
        for (int center=0;center<atom_number[0][0];center++)
        {
            bool judge=false;
            if (frames[0][center].type==ct)
            {
                double min=9999.9;
                for (int neighbor=0;neighbor<atom_number[0][0];neighbor++)
                {
                    if (frames[0][neighbor].type==ne && center!=neighbor)
                    {
                        double x=frames[i][center].x-frames[i][neighbor].x;x=Min_image(x,boundary[0][0]);
                        if (abs(x)>cutoff_r) continue;
                        double y=frames[i][center].y-frames[i][neighbor].y;y=Min_image(y,boundary[0][1]);
                        if (abs(y)>cutoff_r) continue;
                        double z=frames[i][center].z-frames[i][neighbor].z;z=Min_image(z,boundary[0][2]);
                        if (abs(z)>cutoff_r) continue;
                        double r=sqrt(x*x+y*y+z*z);
                        if (r>cutoff_r) continue;
                        else if (r<=cutoff_r && r<min) 
                        {
                            min=r;
                            judge=true;
                        }
                    }
                }
            if (judge)
            {
                long site=static_cast<long>(round(pix*frames[i][center].z));
                bond_lenth[site]+=min;
                n[site]++;
            }
            }
        }
    }
    for (int i=0;i<Nz;i++)
    {
        if (n[i]!=0) bond_lenth[i]/=n[i];
    }
    return bond_lenth;
}

vector<double> LAMMPS_DUMP::Coordination(int ct,int ne,int pix,int interval,double cutoff_r) const
{
    long Nz=static_cast<int>(round(boundary[0][2]*pix))+1;
    vector<double> coordination(Nz,0.0);
    vector<long> n(Nz,0);

    for (int i=0;i<step-interval;i+=interval)
    {
        for (int center=0;center<atom_number[0][0];center++)
        {
            if (frames[0][center].type==ct)
            {
                long site=static_cast<long>(round(pix*frames[i][center].z));
                for (int neighbor=0;neighbor<atom_number[0][0];neighbor++)
                {
                    if (frames[0][neighbor].type==ne && neighbor!=center)
                    {
                        double x=frames[i][center].x-frames[i][neighbor].x;x=Min_image(x,boundary[0][0]);
                        if (abs(x)>cutoff_r) continue;
                        double y=frames[i][center].y-frames[i][neighbor].y;y=Min_image(y,boundary[0][1]);
                        if (abs(y)>cutoff_r) continue;
                        double z=frames[i][center].z-frames[i][neighbor].z;z=Min_image(z,boundary[0][2]);
                        if (abs(z)>cutoff_r) continue;
                        double r=sqrt(x*x+y*y+z*z);
                        if (r<=cutoff_r) coordination[site]++;                      
                    }
                }
                n[site]++;
            }
        }
    }

    for (int i=0;i<Nz;i++)
    {
        if (n[i]!=0) coordination[i]/=n[i];
    }
    return coordination;

}

vector<long> LAMMPS_DUMP::Discharge_Speed(int type,double min_z) const
{
    vector<long> discharge_speed(step,0);
    for (int i=0;i<step;i++)
    {
        for (int j=0;j<atom_number[0][0];j++)
        {
            if (frames[i][j].type==type && frames[i][j].z>min_z) discharge_speed[i]++;
        }
    }

    return discharge_speed;
}

double LAMMPS_DUMP::Min_image(double a,double boundary) const
{
    if (boundary==0.0) return a;
    return a-boundary*round(a/boundary);
}