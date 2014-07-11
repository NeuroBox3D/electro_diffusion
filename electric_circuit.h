// created by Andreas Vogel

#include "common/ug_config.h"
#include "common/error.h"
#include <string.h>
#include <utility>	// for std::pair
#include <vector>
#include <queue>    // std::queue

#include "lib_algebra/cpu_algebra_types.h"
#include "bindings/lua/lua_user_data.h"		// for the current function
#include "common/util/smart_pointer.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace nernst_planck{


/**
 * \brief This class assembles and solves equations arising from electrical circuit analysis.
 *
 * \note It has initially been created to calculate boundary conditions for a fully dimensional
 * Poisson-Nernst-Planck simulation. This plan has been dropped and so this class is but
 * a bit of sediment at the bottom of the nernst_planck_plugin.
 *
 * \author Andreas Vogel, Dinu Patirniche, Markus Breit
 * \date 11-10-2013
 */



class ElectricCircuit
{
	typedef DenseVector< VariableArray1<number> > Vector;
	typedef DenseMatrix< VariableArray2<number> > Matrix;

	public:
		/// adds a capacitor to the circuit
		void add_capacitor(size_t from, size_t to, number capacitance);

		/// adds a resistor to the circuit
		void add_resistor(size_t from, size_t to, number resistance);

		/// adds a voltage source to the circuit
		void add_voltage_source(size_t from, size_t to, number voltage);

		/// adds a current source to the circuit
		void add_current_source(size_t from, size_t to, number current);

		/// sets the initial solution at a node
		void add_initial_solution(size_t node_id, number sol);

	protected:
		/// circuit component base type
		struct Component
		{
			Component(size_t from_, size_t to_) : from(from_), to(to_) {}
			size_t from;
			size_t to;
		};

		/// capacitor type
		struct Capacitor : public Component
		{
			Capacitor(size_t from, size_t to, number capacitance)
				: Component(from, to), C(capacitance) {}
			number C; ///< capacitance
		};

		/// resistor type
		struct Resistor : public Component
		{
			Resistor(size_t from, size_t to, number resistance)
				: Component(from, to), R(resistance), G(1.0/R) {}
			number R; ///< resistance
			number G; ///< conductance
		};

		/// voltage source type
		struct VoltageSource : public Component
		{
			VoltageSource(size_t from, size_t to, number voltage)
				: Component(from, to), U(voltage) {}
			number U; ///< voltage
		};

		/// current source type
		struct CurrentSource : public Component
		{
			CurrentSource(size_t from, size_t to, number current)
				: Component(from, to), I(current) {}
			number I; ///< current
		};

	protected:
		std::vector<Capacitor> vCap;		///< contains all capacitors of the circuit
		std::vector<Resistor> vRes;			///< contains all resistors of the circuit
		std::vector<VoltageSource> vVolt;	///< contains all voltage sources of the circuit
		std::vector<CurrentSource> vCurr;	///< contains all current sources of the circuit

	protected:
		/// number of nodes in the circuit
		size_t num_nodes() const {return vUsedNodes.size();}

		/// sets the number of nodes to be used
		void use_node(size_t n);

		/// checks whether all nodes are used
		bool all_nodes_used() const;

		/// collects unused nodes as string
		std::string unused_nodes() const;

		/// fill incidence matrix for any component type T
		template <typename T>
		void fill_incidence_matrix(Matrix& A, size_t N, const T& vComp);

	protected:
		size_t NN; ///< #nodes (without ground)
		size_t N ; ///< #unknowns (= #nodes + #Current @ VoltSource)

		Matrix AC;	///< incidence matrix capacitors
		Matrix AR;	///< incidence matrix resistors
		Matrix AI;	///< incidence matrix current sources
		Matrix AU;	///< incidence matrix voltage sources

		// original
		Matrix M;	///< system matrix for time-dep. parts of the equation
		Matrix B;	///< system matrix for time-indep. parts
		Vector b;	///< rhs

		std::vector<bool> vUsedNodes;
		std::vector<std::pair<size_t, number> > inits;	///< given initial solutions
		Vector u0;										///< all initial solutions
		number t; 										///< current time

		// flags
		int F_write_to_file;
		FILE *f;

	protected:
		struct sortFunction
		{
			bool operator() (const pair<size_t, number>& p1, const pair<size_t, number>& p2)
			{
				return (p1.first < p2.first);
			}
		};

#ifndef NDEBUG
		/// calculates the condition of a matrix (only for debugging purposes)
		void calculate_condition(Matrix A);
#endif

		/// checks whether given start values are consistent with the algebraic constraints
		/// and calculates the rest of the start values
		void check_and_calculate_initial_solution();

	public:
		/// inits the system
		void init();

	public:
		/// calculates the stationary solution
		std::vector<number> solve_stationary();

	public:
		/// set a new righthand side
		void update_rhs(std::vector<number> rhs);

		/// set whether output is to be written to file
		void write_to_file(bool wtf);

		/// closes the output file
		void close_file() {fclose(f);}

		/// solve time-dependent system using the backward Euler method
		void solve_euler(number delta_t, number h);

		/// solve time-dependent system using the trapezoid rule
		void solve_trapezoid(number tn, number h);

		/// returns the current solution
		std::vector<number> get_solution();

		/// returns the current righthand side
		std::vector<number> get_rhs();
};


} // namespace nernst_planck
} // namespace ug
