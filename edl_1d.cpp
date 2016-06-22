/*
 * edl_1d.cpp
 *
 *  Created in: 03-2013
 *      Author: mbreit
 */

#include "edl_1d.h"

// std includes
#include <iostream>
#include <sstream>

namespace ug {
namespace nernst_planck {


void EDLSimulation::set_geom_specs(number a, number d, number A, number h)
{
	this->a = a;
	this->d = d;
	this->A = A;

	// compute number of unknowns
	size_t N1 = 20;
	size_t N2a = floor(20e-9/h);
	size_t N2b = floor(d/h);
	size_t N2c = floor(20e-9/h);
	size_t N3 = 50;
	size_t N = 1+N1+N2a+N2b+N2c+N3;

	// construct grid according to numbers of unknowns N_i
	grid.resize(N);
	for (size_t i=0; i<=N1; i++) grid[i] = (a-20e-9)/N1*i;
	for (size_t i=1; i<=N2a; i++) grid[N1+i] = (a-20e-9) + 20e-9/N2a*i;
	for (size_t i=1; i<=N2b; i++) grid[N1+N2a+i] = a + d/N2b*i;
	for (size_t i=1; i<=N2c; i++) grid[N1+N2a+N2b+i] = a+d + 20e-9/N2c*i;
	for (size_t i=1; i<=N3; i++) grid[N1+N2a+N2b+N2c+i] = a+d+20e-9 + A/N3*i;

	/*
	std::cout << "grid: " << std::endl;
	for (size_t i=0; i<N; i++) std::cout << grid[i] << ", ";
	std::cout << std::endl;
	*/

	// construct vector of unknowns
	this->N1 = N1+N2a;
	this->N2 = N2b;
	this->N3 = N3+N2c;
	u.resize(4*(this->N1+this->N3+2) + this->N1+this->N2+this->N3+1);
}



// 0:K, 1:Cl, 2:Na, 3:A, 4:Phi
void EDLSimulation::set_bnd_cond(std::vector<number> inner, std::vector<number> outer)
{
	// set intelligent initial solution for Newton iteration

	// dendrite: constant inner concentrations and potential
	for (size_t i=0; i<=N1; i++)
	{
		for (size_t j=0; j<4; j++) u[i+j*(N1+N3+2)] = inner[j];
		u[i+4*(N1+N3+2)] = inner[4];
	}

	// membrane: linear potential
	for (size_t i=N1+1; i<=N1+N2-1; i++) u[i+4*(N1+N3+2)] = inner[4] + (outer[4]-inner[4]) / N2 * (i-N1);

	// outer medium: constant outer concentration and potential
	for (size_t i=N1+1; i<=N1+1+N3; i++)
	{
		for (size_t j=0; j<4; j++) u[i+j*(N1+N3+2)] = outer[j];
		u[i+N2-1+4*(N1+N3+2)] = outer[4];
	}
}


void EDLSimulation::set_constants(std::vector<number> val, std::vector<number> perm)
{
	this->valency = val;
	this->eps_r1 = perm[0];
	this->eps_r2 = perm[1];
	this->eps_r3 = perm[2];
}


void EDLSimulation::set_verbosity_level(size_t lvl)
{
	m_verbLvl = lvl;
}


void EDLSimulation::set_min_reduction(number mr)
{
	m_minRed = mr;
}


void EDLSimulation::set_min_defect(number md)
{
	m_minDefect = md;
}


void EDLSimulation::set_max_iter(size_t mi)
{
	m_maxIt = mi;
}


void EDLSimulation::set_output_file(std::string file)
{
	m_outFile = file;
}


void EDLSimulation::compute_solution()
{
	// number of unknowns
	size_t N = 4*(N1+N3+2) + N1+N2+N3+1;

	// open filestreams for output
	std::vector<std::ofstream*> outStreams;
	for (size_t i = 0; i < 9; i++)
	{
		std::ostringstream oss;
		oss << m_outFile << i;
		outStreams.push_back(new std::ofstream(oss.str().c_str(), std::ofstream::out));
		if (!(*outStreams[i]).is_open())
		{UG_THROW("Opening the desired outFile " << m_outFile << " did not succeed.");}
	}

	// initialize old solution (u_old) and newton solution iterate (u_n)
	Vector u_n = u;

	// declare defect norms of solution iterate and initial solution
	number dn2, d0n2;

	// declare Jacobian and defect
	Matrix J_inv;
	J_inv.resize(N,N);
	Vector d;
	d.resize(N);

	// init defect (and its norm)
	defect(d, u);
	d0n2 = VecTwoNorm(d);

	// output
	if (m_verbLvl >= 1)
	{
		std::cout << "    ---- newton solver ----" << std::endl;
		std::cout << "    initial defect norm: " << d0n2 << std::endl;
	}

	// loop (Newton steps)
	size_t iter = 0;
	do
	{
		// update last solution
		u = u_n;

		// assemble and invert Jacobian
		jacobian(J_inv, u);

		if (!Invert(J_inv))
		{
			UG_THROW("Newton failed, as Jacobian was not invertible.");
		}

		// calculate next Newton iterate u_n = u - J^-1 * d
		MatMultAdd(u_n, 1.0, u, -1.0, J_inv, d);

		// compute new defect and norm
		defect(d, u_n);
		dn2 = VecTwoNorm(d);

		// increment counter
		iter++;

		// output
		if (m_verbLvl >= 1) std::cout << "    " << iter << "  reduction: " << dn2/d0n2 << std::endl;
	}
	while (dn2/d0n2 > m_minRed && dn2 > m_minDefect && iter < m_maxIt);
	// loop is performed as long as defect has not been sufficiently reduced
	// and no more than maxIter iterations have been performed

	// check convergence
	if (dn2/d0n2 <= m_minRed || dn2 <= m_minDefect)
	{
		u = u_n;
		if (m_verbLvl >= 1) std::cout << "    ---- iteration converged ----" << std::endl;
	}
	else
	{
		UG_THROW("Newton iteration did not succeed. "
				 "Maximum iterations reached, but defect reduction only " << dn2/d0n2);
	}

	// output solution to file (treat O1 separately)
	for (size_t j=0; j<5; j++)
	{
		for (size_t i=0; i<=N1+N2+N3; i++)
		{
			if (j<4)
			{
				if (i<=N1) *outStreams[j] << grid[i] << "\t" << u[i+j*(N1+N3+2)] << std::endl;
				else if (i>=N1+N2) *outStreams[j+4] << grid[i] << "\t" << u[i-N2+1+j*(N1+N3+2)] << std::endl;
			}
			else *outStreams[8] << grid[i] << "\t" << u[i+4*(N1+N3+2)] << std::endl;
		}
	}

	// close output streams
	for (size_t i = outStreams.size()-1; i < outStreams.size(); i--)
	{
		(*outStreams[i]).close();
		delete outStreams[i];
		outStreams.pop_back();
	}
}



void EDLSimulation::calc_EDL_ions()
{
	// Riemann integral (trapezoidal rule)
	std::vector<number> I = std::vector<number>(4,0.0);
	std::vector<number> bulk = std::vector<number>(4);
	for (size_t j=0; j<4; j++) bulk[j] = u[j*(N1+N3+2)];

	const number pi = 3.14159265;

	for (size_t i=0; i<N1; i++)
	{
		for (size_t j=0; j<4; j++)
		{
			I[j] += pi * (grid[i+1]-grid[i]) *
					(  grid[i] * (u[i+j*(N1+N3+2)] - bulk[j])
					 + grid[i+1] * (u[i+1+j*(N1+N3+2)] - bulk[j]));
		}
	}

	// save result
	// open filestreams for output
	std::vector<std::ofstream*> outStreams;
	for (size_t i = 0; i < 4; i++)
	{
		std::ostringstream oss;
		oss << m_outFile << "EDL_integral_" << i;
		outStreams.push_back(new std::ofstream(oss.str().c_str(), std::ofstream::app));
		if (!(*outStreams[i]).is_open())
		{UG_THROW("Opening the desired outFile " << m_outFile << " did not succeed.");}
	}
	// compile a list of integral values for various potential values
	for (size_t i=0; i<4; i++) *outStreams[i] << u[4*(N1+N3+2)] << "\t" << I[i] << std::endl;

	// close output streams
	for (size_t i = outStreams.size()-1; i < outStreams.size(); i--)
	{
		(*outStreams[i]).close();
		delete outStreams[i];
		outStreams.pop_back();
	}
}



void EDLSimulation::defect(Vector& d, Vector& u)
{
	/*
	number h1 = this->a / N1;
	number h2 = this->d / N2;
	number h3 = this->A / N3;
	*/
	number h1, h2, c1, c2, c3;
// dendritic parts
	// ion species
	d[0] = d[N1+N3+2] = d[2*(N1+N3+2)] = d[3*(N1+N3+2)] = 0.0;
	for (size_t i = 1; i<=N1; i++)
	{
		for (size_t j=0; j<4; j++)
		{
			d[i+j*(N1+N3+2)] = u[i+j*(N1+N3+2)]-u[i-1+j*(N1+N3+2)]
				   + valency[j]*F/R/T*u[i+j*(N1+N3+2)]*(u[i+4*(N1+N3+2)]-u[i-1+4*(N1+N3+2)]);
		}
	}
	// potential
	d[4*(N1+N3+2)] = 0.0;
	for (size_t i = 1; i<N1; i++)
	{
		// irregular grid schmock
		h1 = grid[i]-grid[i-1];
		h2 = grid[i+1]-grid[i];
		c1 = (h1+h2)/(2*h1);
		c2 = (h1+h2)*(h1+h2)/(2*h1*h2);
		c3 = (h1+h2)/(2*h2);

		d[i+4*(N1+N3+2)] = c1*u[i-1+4*(N1+N3+2)] - c2*u[i+4*(N1+N3+2)] + c3*u[i+1+4*(N1+N3+2)]
						   + h2 / grid[i] * (u[i+4*(N1+N3+2)]-u[i-1+4*(N1+N3+2)]);
		for (size_t j=0; j<4; j++) d[i+4*(N1+N3+2)] += h1*h2*valency[j]*F/eps_0/eps_r1*u[i+j*(N1+N3+2)];
	}

// interface condition for potential (inner): eps_3*phi'_3 - eps_1*phi'_1 + rho_f/eps_0 = 0
	h1 = grid[N1]-grid[N1-1];
	h2 = grid[N1+1]-grid[N1];

	d[N1+4*(N1+N3+2)] = eps_r1/h1*u[N1-1+4*(N1+N3+2)] - (eps_r1/h1+eps_r3/h2)*u[N1+4*(N1+N3+2)]
						+ (eps_r3/h2)*u[N1+1+4*(N1+N3+2)] + rho_f/eps_0;

// membrane parts
	// inner
	for (size_t i = N1+1; i<N1+N2; i++)
	{
		// irregular grid schmock
		h1 = grid[i]-grid[i-1];
		h2 = grid[i+1]-grid[i];
		c1 = (h1+h2)/(2*h1);
		c2 = (h1+h2)*(h1+h2)/(2*h1*h2);
		c3 = (h1+h2)/(2*h2);

		d[i+4*(N1+N3+2)] = c1*u[i-1+4*(N1+N3+2)] - c2*u[i+4*(N1+N3+2)] + c3*u[i+1+4*(N1+N3+2)]
										   + h2 / grid[i] * (u[i+4*(N1+N3+2)]-u[i-1+4*(N1+N3+2)]);
	}
	// left bnd
	/*
	h1 = grid[N1+1]-grid[N1];
	h2 = grid[N1+2]-grid[N1+1];
	c1 = (h1+h2)/(2*h1);
	c2 = (h1+h2)*(h1+h2)/(2*h1*h2);
	c3 = (h1+h2)/(2*h2);

	d[N1+4*(N1+N3+2)] = c1*u[N1+4*(N1+N3+2)] - c2*u[N1+1+4*(N1+N3+2)] + c3*u[N1+2+4*(N1+N3+2)]
						+ h2 / grid[N1] * (u[N1+1+4*(N1+N3+2)]-u[N1+4*(N1+N3+2)]);
	// right bnd
	h1 = grid[N1+N2-1]-grid[N1+N2-2];
	h2 = grid[N1+N2]-grid[N1+N2-1];
	c1 = (h1+h2)/(2*h1);
	c2 = (h1+h2)*(h1+h2)/(2*h1*h2);
	c3 = (h1+h2)/(2*h2);

	d[N1+N2+4*(N1+N3+2)] = c1*u[N1+N2-2+4*(N1+N3+2)] - c2*u[N1+N2-1+4*(N1+N3+2)] + c3*u[N1+N2+4*(N1+N3+2)]
									+ h1 / grid[N1+N2] * (u[N1+N2+4*(N1+N3+2)]-u[N1+N2-1+4*(N1+N3+2)]);
	*/

// interface condition for potential (inner): eps_2*phi'_2 - eps_3*phi'_3 + rho_f/eps_0 = 0
	h1 = grid[N1+N2]-grid[N1+N2-1];
	h2 = grid[N1+N2+1]-grid[N1+N2];

	d[N1+N2+4*(N1+N3+2)] = eps_r3/h1*u[N1+N2-1+4*(N1+N3+2)] - (eps_r3/h1+eps_r2/h2)*u[N1+N2+4*(N1+N3+2)]
						+ (eps_r2/h2)*u[N1+N2+1+4*(N1+N3+2)] + rho_f/eps_0;

// outer parts
	// ion species
	d[N1+N3+2-1] = d[2*(N1+N3+2)-1] = d[3*(N1+N3+2)-1] = d[4*(N1+N3+2)-1] = 0.0;
	for (size_t i = N1+1; i<N1+1+N3; i++)
	{
		for (size_t j=0; j<4; j++)
		{
			d[i+j*(N1+N3+2)] = u[i+1+j*(N1+N3+2)]-u[i+j*(N1+N3+2)]
				   + valency[j]*F/R/T*u[i+j*(N1+N3+2)]*(u[i+(N2-1)+1+4*(N1+N3+2)]-u[i+(N2-1)+4*(N1+N3+2)]);
		}
	}
	// potential
	d[4*(N1+N3+2)+N1+N2+N3] = 0.0;
	for (size_t i = N1+N2+1; i<N1+N2+N3; i++)
	{
		h1 = grid[i]-grid[i-1];
		h2 = grid[i+1]-grid[i];
		c1 = (h1+h2)/(2*h1);
		c2 = (h1+h2)*(h1+h2)/(2*h1*h2);
		c3 = (h1+h2)/(2*h2);

		d[i+4*(N1+N3+2)] = c1*u[i-1+4*(N1+N3+2)] - c2*u[i+4*(N1+N3+2)] + c3*u[i+1+4*(N1+N3+2)]
						   + h1 / grid[i] * (u[i+1+4*(N1+N3+2)]-u[i+4*(N1+N3+2)]);
		for (size_t j=0; j<4; j++) d[i+4*(N1+N3+2)] += h1*h2*valency[j]*F/eps_0/eps_r2*u[i-N2+1+j*(N1+N3+2)];
	}

	//std::cout << d << "  [" << d.size() << "]" << std::endl;
}



void EDLSimulation::jacobian(Matrix& J, Vector& u)
{
	number h1, h2, c1, c2, c3;

	J = 0.0;

// dendritic parts
	// ion species
	J(0,0) = J(N1+N3+2,N1+N3+2) = J(2*(N1+N3+2),2*(N1+N3+2)) = J(3*(N1+N3+2),3*(N1+N3+2)) = 1.0;
	for (size_t i = 1; i<=N1; i++)
	{
		for (size_t j=0; j<4; j++)
		{
			J(i+j*(N1+N3+2),i+j*(N1+N3+2)) = 1.0 + valency[j]*F/R/T*(u[i+4*(N1+N3+2)]-u[i-1+4*(N1+N3+2)]);
			J(i+j*(N1+N3+2),i-1+j*(N1+N3+2)) = -1.0;
			J(i+j*(N1+N3+2),i+4*(N1+N3+2)) = valency[j]*F/R/T*u[i+j*(N1+N3+2)];
			J(i+j*(N1+N3+2),i-1+4*(N1+N3+2)) = -valency[j]*F/R/T*u[i+j*(N1+N3+2)];
		}
	}
	// potential
	J(4*(N1+N3+2),4*(N1+N3+2)) = 1.0;
	for (size_t i = 1; i<N1; i++)
	{
		h1 = grid[i]-grid[i-1];
		h2 = grid[i+1]-grid[i];
		c1 = (h1+h2)/(2*h1);
		c2 = (h1+h2)*(h1+h2)/(2*h1*h2);
		c3 = (h1+h2)/(2*h2);

		J(i+4*(N1+N3+2),i+4*(N1+N3+2)) = -c2 + h2 / grid[i];
		J(i+4*(N1+N3+2),i-1+4*(N1+N3+2)) = c1 - h2 / grid[i];
		J(i+4*(N1+N3+2),i+1+4*(N1+N3+2)) = c3;
		for (size_t j=0; j<4; j++) J(i+4*(N1+N3+2),i+j*(N1+N3+2)) = h1*h2*valency[j]*F/eps_0/eps_r1;
	}

// interface condition for potential (inner): eps_1*phi'_1 - eps_2*phi'_2 = 0
	h1 = grid[N1]-grid[N1-1];
	h2 = grid[N1+1]-grid[N1];

	J(N1+4*(N1+N3+2),N1-1+4*(N1+N3+2)) = eps_r1/h1;
	J(N1+4*(N1+N3+2),N1+4*(N1+N3+2)) = -(eps_r1/h1+eps_r3/h2);
	J(N1+4*(N1+N3+2),N1+1+4*(N1+N3+2)) = eps_r3/h2;

// membrane parts
	// inner
	for (size_t i = N1+1; i<N1+N2; i++)
	{
		h1 = grid[i]-grid[i-1];
		h2 = grid[i+1]-grid[i];
		c1 = (h1+h2)/(2*h1);
		c2 = (h1+h2)*(h1+h2)/(2*h1*h2);
		c3 = (h1+h2)/(2*h2);

		J(i+4*(N1+N3+2),i+4*(N1+N3+2)) = -c2 + h2 / grid[i];
		J(i+4*(N1+N3+2),i-1+4*(N1+N3+2)) = c1 - h2 / grid[i];
		J(i+4*(N1+N3+2),i+1+4*(N1+N3+2)) = c3;
	}
	/*
	// left bnd
	h1 = grid[N1+1]-grid[N1];
	h2 = grid[N1+2]-grid[N1+1];
	c1 = (h1+h2)/(2*h1);
	c2 = (h1+h2)*(h1+h2)/(2*h1*h2);
	c3 = (h1+h2)/(2*h2);

	J(N1+4*(N1+N3+2),N1+2+4*(N1+N3+2)) = c3;
	J(N1+4*(N1+N3+2),N1+1+4*(N1+N3+2)) = -c2 + h2 / grid[N1];
	J(N1+4*(N1+N3+2),N1+4*(N1+N3+2)) = c1 - h2 / grid[N1];

	// right bnd
	h1 = grid[N1+N2-1]-grid[N1+N2-2];
	h2 = grid[N1+N2]-grid[N1+N2-1];
	c1 = (h1+h2)/(2*h1);
	c2 = (h1+h2)*(h1+h2)/(2*h1*h2);
	c3 = (h1+h2)/(2*h2);

	J(N1+N2+4*(N1+N3+2), N1+N2-2+4*(N1+N3+2)) = c1;
	J(N1+N2+4*(N1+N3+2), N1+N2-1+4*(N1+N3+2)) = -c2 - h1 / grid[N1+N2];
	J(N1+N2+4*(N1+N3+2), N1+N2+4*(N1+N3+2)) = c3 + h1 / grid[N1+N2];
	*/

// interface condition for potential (inner): eps_1*phi'_1 - eps_2*phi'_2 = 0
	h1 = grid[N1+N2]-grid[N1+N2-1];
	h2 = grid[N1+N2+1]-grid[N1+N2];

	J(N1+N2+4*(N1+N3+2),N1+N2-1+4*(N1+N3+2)) = eps_r3/h1;
	J(N1+N2+4*(N1+N3+2),N1+N2+4*(N1+N3+2)) = -(eps_r3/h1+eps_r2/h2);
	J(N1+N2+4*(N1+N3+2),N1+N2+1+4*(N1+N3+2)) = eps_r2/h2;

// outer parts
	// ion species
	J(N1+N3+2-1,N1+N3+2-1) = J(2*(N1+N3+2)-1,2*(N1+N3+2)-1) = J(3*(N1+N3+2)-1,3*(N1+N3+2)-1) = J(4*(N1+N3+2)-1,4*(N1+N3+2)-1) = 1.0;
	for (size_t i = N1+1; i<N1+1+N3; i++)
	{
		for (size_t j=0; j<4; j++)
		{
			J(i+j*(N1+N3+2),i+j*(N1+N3+2)) = -1.0 + valency[j]*F/R/T * (u[i+(N2-1)+1+4*(N1+N3+2)]-u[i+(N2-1)+4*(N1+N3+2)]);
			J(i+j*(N1+N3+2),i+1+j*(N1+N3+2)) = 1.0;
			J(i+j*(N1+N3+2),i+(N2-1)+4*(N1+N3+2)) = -valency[j]*F/R/T*u[i+j*(N1+N3+2)];
			J(i+j*(N1+N3+2),i+(N2-1)+1+4*(N1+N3+2)) = valency[j]*F/R/T*u[i+j*(N1+N3+2)];
		}
	}
	// potential
	J(4*(N1+N3+2)+N1+N2+N3,4*(N1+N3+2)+N1+N2+N3) = 1.0;
	for (size_t i = N1+N2+1; i<N1+N2+N3; i++)
	{
		h1 = grid[i]-grid[i-1];
		h2 = grid[i+1]-grid[i];
		c1 = (h1+h2)/(2*h1);
		c2 = (h1+h2)*(h1+h2)/(2*h1*h2);
		c3 = (h1+h2)/(2*h2);

		J(i+4*(N1+N3+2),i+4*(N1+N3+2)) = -c2 - h1 / grid[i];
		J(i+4*(N1+N3+2),i-1+4*(N1+N3+2)) = c1;
		J(i+4*(N1+N3+2),i+1+4*(N1+N3+2)) = c3 + h1 / grid[i];
		for (size_t j=0; j<4; j++) J(i+4*(N1+N3+2),i-N2+1+j*(N1+N3+2)) = h1*h2*valency[j]*F/eps_0/eps_r2;
	}
	/*
	size_t N = 4*(N1+N3+2) + N1+N2+N3+1;
	std::cout << "[ ";
	for (size_t i=0; i<N; i++)
	{
		for (size_t j=0; j<N; j++) std::cout << J(i,j) << " ";
		if (i != N-1) std::cout << std::endl << "  ";
	}
	std::cout << " ]" << std::endl;
	*/
}



} // namespace nernst_planck
} // namespace ug
