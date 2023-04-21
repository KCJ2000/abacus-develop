#ifndef H_EWALD_PW_H
#define H_EWALD_PW_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
#include "module_hamilt_pw/hamilt_pwdft/stress_func.h"

class H_Ewald_pw 
{
	public:

	// the calculation of ewald force is essential to Stress_Func and Forces
    // And the processes are the same
    // So here use friend class to give Stress_Function and Forces the power to get ewald force in this class.
    // But such strategy is not flexible because of the use of the static varible and funcitons.
    // The static varible does not belong to any object.In a calculation,the ewld is unique.
    // When ABACUS wants to deal with several systems at a time, this class will be a barrier.
    // kongfanhan - 23.4.17
    friend class Stress_Func<double, psi::DEVICE_CPU>; // Ewald stress
    friend class Stress_Func<double, psi::DEVICE_GPU>; // Ewald stress
    friend class Forces<double, psi::DEVICE_CPU>; // Ewald forces
    friend class Forces<double, psi::DEVICE_GPU>; // Ewald forces
	friend class Force_LCAO; // Ewald forces

    // stardard part , in this class ,these two functions have no contains.
    H_Ewald_pw();
    ~H_Ewald_pw();

	// the Ewald energy
    static double ewald_energy;

	// compute the Ewald energy
    static void compute_ewald(const UnitCell &cell, ModulePW::PW_Basis* rho_basis);

	private:

    static void rgen(
        const ModuleBase::Vector3<double> &dtau,
        const double &rmax,
        int *irr,
        const ModuleBase::Matrix3 &at,
        const ModuleBase::Matrix3 &bg,
        ModuleBase::Vector3<double> *r,
        double *r2,
        int  &nrm
    );

	// the coefficient of ewald method
	static double alpha;
    static int mxr;

};

#endif //ewald energy
