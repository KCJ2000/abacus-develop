#include"../H_Ewald_pw.h"
#include"gtest/gtest.h"
#include"gmock/gmock.h"
#include<iostream>
#include "module_cell/unitcell.h"
#include "module_base/mathzone.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_base/vector3.h"


/*********************************************************
 *  unit test of class H_Ewald_pw and related functions.
**********************************************************/


/**
 * -Tested functions:
 *  -H_Ewald_pw::
*/


//copy from vdw_test.cpp
//add zv init in construct_ucell

double H_Ewald_pw::alpha=0.0;
int H_Ewald_pw::mxr = 10;//for test 10 is large enough
double H_Ewald_pw::ewald_energy=0.0;

struct atomtype_
{
    std::string atomname;
    int Z_a;//the add part
    std::vector<std::vector<double>> coordinate;
};

struct stru_
{
    std::vector<double> cell;
    std::vector<atomtype_> all_type;
};

void construct_ucell(stru_ &stru, UnitCell &ucell)
{
    std::vector<atomtype_> coord = stru.all_type;
    ucell.a1 = ModuleBase::Vector3<double>(stru.cell[0], stru.cell[1], stru.cell[2]);
    ucell.a2 = ModuleBase::Vector3<double>(stru.cell[3], stru.cell[4], stru.cell[5]);
    ucell.a3 = ModuleBase::Vector3<double>(stru.cell[6], stru.cell[7], stru.cell[8]);
    ucell.latvec.e11 = stru.cell[0]; ucell.latvec.e12 = stru.cell[1]; ucell.latvec.e13 = stru.cell[2];
    ucell.latvec.e21 = stru.cell[3]; ucell.latvec.e22 = stru.cell[4]; ucell.latvec.e23 = stru.cell[5];
    ucell.latvec.e31 = stru.cell[6]; ucell.latvec.e32 = stru.cell[7]; ucell.latvec.e33 = stru.cell[8];
    ucell.GT = ucell.latvec.Inverse();
    ucell.G = ucell.GT.Transpose();//reciprocal space
    ucell.lat0 = 10.2;
    ucell.ntype = stru.all_type.size();
    ucell.atoms = new Atom[ucell.ntype];
    ucell.nat = 0;
    int nmax = 0;
    for (int i = 0; i < coord.size(); i++)
    {
        ucell.atoms[i].ncpp.zv = coord[i].Z_a;//the add part
        ucell.atoms[i].label = coord[i].atomname;
        ucell.atoms[i].ncpp.psd = coord[i].atomname;
        ucell.atoms[i].na = coord[i].coordinate.size();
        ucell.atoms[i].tau = new ModuleBase::Vector3<double>[ucell.atoms[i].na];
        ucell.atoms[i].taud = new ModuleBase::Vector3<double>[ucell.atoms[i].na];
        for (int j = 0; j < ucell.atoms[i].na; j++)
        {
            std::vector<double> this_atom = coord[i].coordinate[j];
            ucell.atoms[i].tau[j] = ModuleBase::Vector3<double>(this_atom[0], this_atom[1], this_atom[2]);
            ModuleBase::Mathzone::Cartesian_to_Direct(ucell.atoms[i].tau[j].x,
                                                      ucell.atoms[i].tau[j].y,
                                                      ucell.atoms[i].tau[j].z,
                                                      ucell.a1.x,
                                                      ucell.a1.y,
                                                      ucell.a1.z,
                                                      ucell.a2.x,
                                                      ucell.a2.y,
                                                      ucell.a2.z,
                                                      ucell.a3.x,
                                                      ucell.a3.y,
                                                      ucell.a3.z,
                                                      ucell.atoms[i].taud[j].x,
                                                      ucell.atoms[i].taud[j].y,
                                                      ucell.atoms[i].taud[j].z);
        }
        ucell.nat += ucell.atoms[i].na;
        if (ucell.atoms[i].na > nmax)
        {
            nmax = ucell.atoms[i].na;
        }

    }

    ucell.omega = abs(ucell.latvec.Det()) * ucell.lat0 * ucell.lat0 * ucell.lat0;
    ucell.itia2iat.create(ucell.ntype, nmax);
	int iat=0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			ucell.itia2iat(it, ia) = iat;
            ++iat;
		}	
	}
}

void Clear_ucell(UnitCell &ucell)
{
    for (int i = 0; i < ucell.ntype; i++)
    {
        delete[] ucell.atoms[i].tau;
        delete[] ucell.atoms[i].taud;
    }
    delete[] ucell.atoms;
}

void construct_pw_Basis(ModulePW::PW_Basis* rho_basis){
    //this construction has no physics explanation
    //just a numerical test
    rho_basis->ggecut = 1.0;
    rho_basis->ig_gge0 = 1;
    rho_basis->npw = 3;
    double k[3]={1.0,1.0,1.0};
    rho_basis->gg = k;
}

class EwaldTest: public testing::Test
{
    protected:
    UnitCell ucell;
    ModulePW::PW_Basis* rho_basis;
    

    void SetUp(){
        stru_ structure{std::vector<double>{0.5,0.5,0.0,0.5,0.0,0.5,0.0,0.5,0.5},
                        std::vector<atomtype_>{atomtype_{"Si",4,
                                                         std::vector<std::vector<double>>{
                                                             {0., 0., 0.},
                                                             {0.3, 0.25, 0.25}
                                                         }}}};
        construct_ucell(structure,ucell);
        construct_pw_Basis(rho_basis);
    }
    
    void TearDown(){
        Clear_ucell(ucell);
    }
};

TEST_F(EwaldTest, compute_ewald_test){

}

TEST_F(EwaldTest,rgen_test){
    double rmax = 1.0;
    int nrm = 0;
    
    ModuleBase::Vector3<double> *r = new ModuleBase::Vector3<double>[H_Ewald_pw::mxr];
    
    double *r2 = new double[H_Ewald_pw::mxr];
    int* irr = new int[H_Ewald_pw::mxr];
    ModuleBase::Vector3<double> dtau(1.0,1.0,1.0);
    H_Ewald_pw::rgen(
        dtau,
        rmax,
        irr,
        ucell.latvec,
        ucell.G,
        r,
        r2,
        nrm);
    
}