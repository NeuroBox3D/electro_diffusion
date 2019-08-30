/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Andreas Vogel
 * Creation date: <= 2013
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


#include "electric_circuit.h"

#include <algorithm>                                                             // for sort
#include <cmath>                                                                 // for abs, fabs
#include <iostream>                                                              // for operator<<

#include "common/error.h"                                                        // for UG_THROW
#include "common/math/math_vector_matrix/math_vector_functions.h"                // for VecTwoNorm
#include "lib_algebra/common/operations_mat/operations_mat.h"                    // for MatMult
#include "lib_algebra/small_algebra/small_matrix/densematrix_inverse.h"          // for Invert
#include "lib_algebra/small_algebra/small_matrix/densematrix_operations.h"       // for MatMult


using namespace std;

namespace ug {
namespace nernst_planck {


void ElectricCircuit::add_capacitor(size_t from, size_t to, number capacitance)
{
	vCap.push_back(Capacitor(from, to, capacitance));
	use_node(from); use_node(to);
}

void ElectricCircuit::add_resistor(size_t from, size_t to, number resistance)
{
	vRes.push_back(Resistor(from, to, resistance));
	use_node(from); use_node(to);
}

void ElectricCircuit::add_voltage_source(size_t from, size_t to, number voltage)
{
	vVolt.push_back(VoltageSource(from, to, voltage));
	use_node(from); use_node(to);
}

void ElectricCircuit::add_current_source(size_t from, size_t to, number current)
{
	vCurr.push_back(CurrentSource(from, to, current));
	use_node(from); use_node(to);
}

void ElectricCircuit::add_initial_solution(size_t node_id, number sol)
{
	inits.push_back(std::pair<size_t, number>(node_id-1, sol));
}

void ElectricCircuit::use_node(size_t n)
{
	if (n >= vUsedNodes.size()) vUsedNodes.resize(n+1,false);
	vUsedNodes[n] = true;
}

bool ElectricCircuit::all_nodes_used() const
{
	for (size_t n = 0; n < vUsedNodes.size(); ++n)
		if (!vUsedNodes[n])
			return false;
	return true;
}

std::string ElectricCircuit::unused_nodes() const
{
	std::stringstream ss;
	for (size_t n = 0; n < vUsedNodes.size(); ++n)
		if (!vUsedNodes[n])
		{
			if (!ss.str().empty()) ss<< ", ";
			ss << n+1;
		}

	return ss.str();
}

template <typename T>
void ElectricCircuit::fill_incidence_matrix(Matrix& A, size_t N, const T& vComp)
{
	const size_t N_c = vComp.size();
	A.resize(N,N);
	A = 0.0;
	for (size_t k = 0; k < N_c; ++k)
	{
		if (vComp[k].from > 0) A(vComp[k].from - 1, k) = 1;
		if (vComp[k].to > 0)   A(vComp[k].to   - 1, k) = -1;
	}
}

#ifndef NDEBUG
void ElectricCircuit::calculate_condition(Matrix A)
{
	// find biggest EV
	number delta = 1.0;
	size_t steps = 0;
	Vector x, x_new;
	x.resize(A.num_cols());
	x_new.resize(A.num_cols());
	x = 1.0/A.num_cols();

	while (delta > 1e-8 && steps < 1000)
	{
		MatMult(x_new, 1.0, A, x);
		x_new = x_new * (1.0 / VecTwoNorm(x_new));
		delta = VecTwoNorm(x_new-x);
		x = x_new;
		steps++;
	}
	MatMult(x_new, 1.0, A, x);
	number lambda_max = VecTwoNorm(x_new)  / VecTwoNorm(x);

	// find smallest EV
	if (!Invert(A))
	{
		UG_THROW("ElectricCircuit: Cannot invert impl. Euler matrix,"
							   " maybe time step too small ");
	}
	delta = 1.0;
	steps = 0;
	x = 1.0/A.num_cols();

	while (delta > 1e-8 && steps < 1000)
	{
		MatMult(x_new, 1.0, A, x);
		x_new = x_new * (1.0 / VecTwoNorm(x_new));
		delta = VecTwoNorm(x_new-x);
		x = x_new;
	}
	MatMult(x_new, 1.0, A, x);
	number lambda_min = VecTwoNorm(x) / VecTwoNorm(x_new);

	std::cout << std::endl << "matrix condition is: " << lambda_max / lambda_min << "." << std::endl;
}
#endif

void ElectricCircuit::check_and_calculate_initial_solution()
{
// check given initial conditions for completeness and consistency
	// find number of purely algebraic (non-differential) equations
	size_t NA = vVolt.size();

	vector<size_t> vZeroLines;
	for (size_t i=0; i<AC.num_rows(); i++)
	{
		bool bZeroLine = true;
		for (size_t j=0; j<AC.num_cols(); j++)
		{
			if (AC(i,j))
			{
				bZeroLine = false;
				break;
			}
		}
		if (bZeroLine) vZeroLines.push_back(i);
	}
	NA += vZeroLines.size();
//std::cout << "zeroLines: ";
//for (size_t i = 0; i < vZeroLines.size(); i++) std::cout << vZeroLines[i] << ", ";
//std::cout << std::endl;

	bool capAtGround = false;
	for (size_t i=0; i<vCap.size(); i++)
	{
		if (vCap[i].from == 0 || vCap[i].to == 0)
		{
			capAtGround = true;
			break;
		}
	}
	if (!capAtGround) NA++; // non-zero capacity matrix lines are lin.dep.

	if (inits.size() < N-NA)
		UG_THROW("ElectricCircuit: Not enough initial conditions supplied!\n"
				"You need at least "<< N-NA <<", but only provided " << inits.size() << ".");

//std::cout << "C = " << M << std::endl;
//std::cout << "G = " << B << std::endl;

	// construct row modification matrix
	Matrix S;
	S.resize(NA, N);
	S = 0.0;
	size_t currZRow = 0;
	vector<size_t>::iterator itZL = vZeroLines.begin();
	// line swaps
	for (size_t i=0; i<vVolt.size(); i++)
		S(NA-vVolt.size()+i, N-vVolt.size()+i) = 1;
	for (size_t i=0; i<N-vVolt.size(); i++)
	{
		if (itZL != vZeroLines.end() && i == *itZL)
		{
			S(currZRow,i) = 1;
			currZRow++;
			itZL++;
		}
		else if (!capAtGround)
		{
			S(NA-vVolt.size()-1,i) = 1;
		}
	}
//std::cout << "S = " << S << std::endl;
	// multiplying S*C gives a NAxN matrix with only zeros
	// multiplying S*G will give the algebraic matrix to be "inverted" to calculate S(i-Gu)=0

	// try solving with given initial values
	std::cout << "#algebraic terms: " << NA << std::endl;
	std::cout << "#free variables: " << N-NA << std::endl;
	std::cout << "#given variables: " << inits.size() << std::endl << std::endl;
	if (inits.size() > N-NA)
	{
		UG_THROW("Too many initial conditions supplied!\n"
				 "The system only has " << N-NA << " free variables, you provided " << inits.size() << "!");
	}
	// construct matrix to be inverted
	Matrix invertMe;
	invertMe = S * B;
//std::cout << "S*G = " << invertMe << std::endl;
	std::sort (inits.begin(), inits.end(), sortFunction());

	Matrix T;
	T.resize(N, NA);
	T = 0.0;
	size_t currCol = 0;
	Vector add2rhs;
	add2rhs.resize(N);
	add2rhs = 0.0;
	vector<pair<size_t, number> >::iterator itKV = inits.begin();
	// line swaps
	for (size_t i=0; i<N; i++)
	{
		if (i == itKV->first)
		{
			add2rhs[i] = itKV->second;
			itKV++;
		}
		else
		{
			T(i, currCol) = 1;
			currCol++;
		}
	}
	//std::cout << "T = " << T << std::endl;
	Vector rhs;
	rhs.resize(NA);
	rhs = S * b - invertMe * add2rhs;

	invertMe = invertMe * T;

//std::cout << "S*G*T = " << invertMe << std::endl;
	if(!Invert(invertMe))
		UG_THROW("ElectricCircuit: Algebraic matrix is not invertible.\n"
				 "Maybe the provided start values are contradictory!");

	u0 = T * invertMe * rhs + add2rhs;
//std::cout << "u0 = " << u0 << std::endl;
}

void ElectricCircuit::init()
{
	if (num_nodes() < 1)
		UG_THROW("ElectricCircuit: No nodes given.")

	// size of independent potentials (node 0 is always ground)
	NN = num_nodes() - 1;

	// number of unknowns (= #nodes + #VoltSources)
	N = NN + vVolt.size();

	if (!all_nodes_used())
		UG_THROW("ElectricCircuit: Given "<<NN+1<<" nodes [0, "<<NN<<"] in the circuit "
				 "the following nodes are unused: "<<unused_nodes());


	fill_incidence_matrix(AC, NN, vCap);
	fill_incidence_matrix(AR, NN, vRes);
	fill_incidence_matrix(AI, NN, vCurr);
	fill_incidence_matrix(AU, NN, vVolt);

	// fill capacitive nodal equations matrix
	M.resize(N, N);
	M = 0.0;
	for (size_t i = 0; i < NN; ++i)
		for (size_t j = 0; j < NN; ++j)
			for (size_t k = 0; k < vCap.size(); ++k)
				M(i,j) += AC(i,k) * vCap[k].C * AC(j,k);

	// fill resistive nodal equations matrix
	B.resize(N, N);
	B = 0.0;
	for (size_t i = 0; i < NN; ++i)
		for (size_t j = 0; j < NN; ++j)
			for (size_t k = 0; k < vRes.size(); ++k)
				B(i,j) += AR(i,k) * vRes[k].G * AR(j,k);

	// add currents over voltage sources to nodal equations (matrix part)
	for (size_t i = 0; i < N; ++i)
		for (size_t k = 0; k < vVolt.size(); ++k)
			B(i,k+NN) = abs(AU(i,k));

	// add voltage source equations (matrix part)
	for (size_t k = 0; k < vVolt.size(); ++k)
		for (size_t i = 0; i < N; ++i)
			B(k+NN, i) = abs(AU(i,k));

	b.resize(N);
	b = 0.0;

	// add voltage source equations (rhs part)
	for (size_t k = 0; k < vVolt.size(); ++k)
		b[k+NN] = vVolt[k].U;

	// add current sources
	for (size_t i = 0; i < NN; ++i)
		for (size_t k = 0; k < vCurr.size(); ++k)
			b[i] -= AI(i,k) * vCurr[k].I;

	check_and_calculate_initial_solution();

	F_write_to_file = 0;
	t=0;
}

std::vector<number> ElectricCircuit::solve_stationary()
{
	Vector x;
	x.resize(N);
	InverseMatMult(x, 1.0, B, b);

	std::vector<number> sol(NN);
	for(size_t i = 0; i < NN; ++i){
		sol[i] = x[i];
	}
	return sol;
}

// YES! We don't need an electric circuit solver per-se, but something that we can use as boundary condition!
void ElectricCircuit::update_rhs(std::vector<number> rhs)
{
//std::cout<<"";// << "(*!*)The input string was: "<< std::endl;

	for (size_t i=0; i<rhs.size(); i++)
	{
		if (i < b.size()) b[i] = rhs[i];
		else {UG_THROW("Given rhs vector is too long.");}
	}

}

void ElectricCircuit::write_to_file(bool wtf)
{
	if (wtf)
	{
		F_write_to_file = 1;
		std::cout << "[***From ElectricCircuit Plugin***]: The solution for the circuit "
					 "are written to the output file!" << std::endl;
		f = fopen("circuit_output.txt","w+");
	}
	else
	{
		F_write_to_file = 0;
	}
}

void ElectricCircuit::solve_euler(number delta_t, number h)
{
/*
	// corrections for better conditioning
	if (t == 0)
	{
		for (size_t i = NN; i < N; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				B(i,j) *= 1e-8;
				B(j,i) *= 1e-8;
			}
			b[i] *= 1e-8;
			u0[i] /= 1e-8;
		}
	}
*/
	Vector u, tmp;
	u = u0;
	tmp.resize(N);

	// output initial solution
	// invert system matrix
	Matrix inv;
	inv.resize(N,N);
	inv = h*B + M;

//std::cout << inv << std::endl;
#ifndef NDEBUG
	if (t==0) calculate_condition(inv);
#endif

	if (!Invert(inv))
	{
		UG_THROW("ElectricCircuit: Cannot invert impl. Euler matrix,"
							   " maybe time step too small ");
	}
	//std::cout << "t: " << t << std::endl;
	//std::cout << "B: " << B << std::endl;
	//std::cout << "M: " << M << std::endl;

	number relative_time = 0;

	while (relative_time < delta_t)
	{
		// compute next time point
		if (relative_time+2*h > delta_t)
		{
			// check if out of bounds, if yes:
			// set to final time point and adjust step size accordingly
			h = (relative_time+h > delta_t) ? delta_t - relative_time : (delta_t-relative_time) / (2-1e-8);
			relative_time += h;
			inv = h*B + M;
#ifndef NDEBUG
	if (t==0) calculate_condition(inv);
#endif
			if (!Invert(inv))
			{
				UG_THROW("ElectricCircuit: Cannot invert impl. Euler matrix,"
						" maybe time step too small ");
			}
		}
		else relative_time += h;

		// perform backward Euler step
		MatMultAdd(tmp, h, b, 1.0, M, u);
		MatMult(u, 1.0, inv, tmp);
	}

	u0=u;
	t += relative_time;
	if (F_write_to_file == 1)
	{
			fprintf(f, "%g ", t);
			for (size_t i = 0; i < N; ++i) fprintf(f, "%g ", u[i]);
			fprintf(f,"\n");
	}
}

void ElectricCircuit::solve_trapezoid(number tn, number h)
{
	// initial solution is zero everywhere
	number t = 0;
	Vector u, tmp, b_old;
	u = u0;
	tmp.resize(N);
	b_old = b;	//TODO:  this will have to be changed, when real time-dep. currents are added

	// invert system matrix
	Matrix inv;
	inv.resize(N,N);
	inv = 0.5*h*B + M;
	if(!Invert(inv)) UG_THROW("ElectricCircuit: Cannot invert trapezoid scheme matrix,"
							  " maybe time step too small ");

	while (fabs(t-tn) > 1e-8*tn)
	{
	   // compute next time point
	   t = t + h;

	   // check if out of bounds, if yes:
	   // set to final time point and adjust step size accordingly
	   if (t > tn)
	   {
		   h -= t - tn;
		   t = tn;
		   inv =  0.5*h*B + M;
		   if (!Invert(inv))
			   UG_THROW("ElectricCircuit: Cannot invert trapezoid scheme matrix,"
						" maybe time step too small ");
	   }

	   // perform trapezoidal scheme step
	   MatMultAdd(tmp, 0.5*h, b+b_old, 1.0, M - 0.5*h*B, u);
	   MatMult(u, 1.0, inv, tmp);

	   b_old = b;

	   // output of current solution
	   if (F_write_to_file == 1)
	   {
		   for (size_t i = 0; i < NN; ++i) fprintf(f,"%f ",u[i]);
		   fprintf(f,"\n");
	   }
	   else std::cout << "time "<< t << ": " << u << std::endl;
	}
}

std::vector<number> ElectricCircuit::get_solution()
{
	std::vector<number> sol(N);
	for (size_t i = 0; i<N; i++)
		sol[i]=u0[i];

	return sol;
}

std::vector<number> ElectricCircuit::get_rhs()
{
	std::vector<number> rhs(N);

	for (size_t i = 0; i<N; i++)
		rhs[i]=b[i];

	return rhs;
}

} // namespace nernst_planck
} // namespace ug
