/*
 * morpho_gen.h
 *
 *  Created on: 07.09.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__MORPHO_GEN_H
#define UG__PLUGINS__NERNST_PLANCK__MORPHO_GEN_H

#include <cstddef>                                               // for size_t
#include <string>                                                // for string
#include <vector>                                                // for vector

#include "common/types.h"                                        // for number
#include "common/math/ugmath_types.h"                            // for vector3
#include "lib_algebra/small_algebra/small_matrix/densevector.h"  // for DenseVector
#include "lib_algebra/small_algebra/storage/variable_array.h"    // for Variable...
#include "lib_grid/grid/grid.h"                                  // for Grid
#include "lib_grid/attachments/attachment_pipe.hpp"              // for Attachme...
#include "lib_grid/common_attachments.h"                         // for APosition
#include "lib_grid/refinement/projectors/projection_handler.h"   // for ProjectionHandler
#include "lib_grid/tools/selector_grid.h"                        // for Selector
#include "lib_grid/tools/subset_handler_grid.h"                  // for SubsetHandler


namespace ug {

// forward declarations
class Vertex;
class Edge;
class Face;
template <typename TStorage> class DenseMatrix;


namespace nernst_planck {


class MorphoGen
{
	public:
		typedef Grid::VertexAttachmentAccessor<APosition> AAPosition;
		typedef Grid::FaceAttachmentAccessor<ANormal> AANormal;

	public:
		MorphoGen();
		virtual ~MorphoGen() {};

		void set_num_neck_filaments(size_t nFil);
		void set_num_filaments(size_t nFil);
		void set_fil_anisotropic(bool filAniso);
		void set_with_refinement_projector(bool withRP);
		void set_seed(size_t seed);
		void set_randomized(bool rand);
		void set_membrane_envelope_radius(number mem_env_rad);
		void set_filament_envelope_radius(number fil_env_rad);
		void set_resolution(size_t nRimVrt);

		/// creates a 3d spine morphology with 1d extensions
		void create_dendrite(const std::string& filename);

		/// creates a simple 2d long dendrite (usable for 1d/3d param optimization)
		void create_dendrite_2d(const std::string& filename, number radius);

		/// creates a simple 1d long dendrite (usable for 1d/3d param optimization)
		void create_dendrite_1d(const std::string& filename);


	protected:
		void create_shaft();
		void graft_spine();
		void make_neck_filaments();
		void make_spherical_filaments();
		void create_membrane_and_envelopes();
		void close_ends();
		void create_bounding_box();
		void tetrahedralize();
		void fix_dendritic_ends();
		void create_interfaces();
		void create_extensions();

	private:
		void create_circle(const vector3& center, const vector3& axis, number radius, size_t numRimVertices);

		void smooth_branching_point(number alpha, int numIterations, const vector3& spine_anchor);
		bool is_cut_by_spine_shaft
		(
			Face* f,
			const vector3& center,
			number radius,
			std::vector<Edge*>& outCutEdges,
			std::vector<size_t>& outCutEdgeIndices,
			std::vector<vector3>& outCutPositions
		);
		bool is_outside_spine_shaft(Vertex* v, const vector3& center, number radius);

		void defect_for_filament_distribution
		(
			const DenseVector<VariableArray1<number> >& bndPos,
			const DenseVector<VariableArray1<number> >& solNew,
			const DenseVector<VariableArray1<number> >& solOld,
			DenseVector<VariableArray1<number> >& defOut,
			number dt
		) const;

		void jacobian_for_filament_distribution
		(
			const DenseVector<VariableArray1<number> >& bndPos,
			const DenseVector<VariableArray1<number> >& solNew,
			DenseMatrix<VariableArray2<number> >& jacOut,
			number dt
		) const;

		void distribute_filaments(const DenseVector<VariableArray1<number> >& bndPos, DenseVector<VariableArray1<number> >& posOut);

		struct MinDistCalculator
		{
			void calculate_minDist
			(
				number& minDistOut,
				size_t& minIndOut,
				size_t ind,
				const std::vector<vector3>& vPoints,
				const std::vector<size_t>& vInd,
				size_t nPts
			);

			number neck_radius;
			number head_radius;
			number filament_envelope_radius;
			number filament_width;
			number neck_end_z;
			number head_center_z;
			number head_opening_theta;
		};

		void create_envelope(const std::vector<int>& vExtrudeSI, number offset, int newVolSI, const std::vector<int>& vNewFrontSI);
		void create_envelope(const std::vector<int>& vExtrudeSI, number offset, int newVolSI)
		{
			size_t sz = vExtrudeSI.size();
			std::vector<int> vFront(sz);
			for (size_t i = 0; i < sz; ++i)
				vFront[i] = newVolSI;
			create_envelope(vExtrudeSI, offset, newVolSI, vFront);
		}

		template <typename TIterator>
		void replace_subsets
		(
			const TIterator& begin,
			const TIterator& end,
			const std::vector<int>& oldSI,
			const std::vector<int>& newSI,
			bool deselect = false
		);

		void tetrahedralize_selection(int tetSI, const std::vector<vector3>* pvHoles = NULL);

		void find_ring_at_dendritic_ends
		(
			int si,
			int sideMultiplier,
			std::vector<Vertex*>& vrts,
			std::vector<Edge*>& edges
		);

		void mark_straight_edge_connection(Vertex* v0, Vertex* v1);
		void create_mark_straight_edge_connection(Vertex* v0, Vertex* v1);

	private:
		// sort index array using value array
		struct IndexCompare
		{
			IndexCompare(const std::vector<number>& _v) : v(_v) {}
			bool operator()(const size_t& a, const size_t& b) {return v[a] < v[b];}

			const std::vector<number>& v;
		};

	private:
		Grid 			m_grid;
		SubsetHandler	m_sh;
		Selector		m_sel;
		AAPosition		m_aaPos;
		AANormal		m_aaNorm;
		vector3			m_pivot;

		// projectors
		SubsetHandler		m_shProj;	// for projectors only
		ProjectionHandler	m_projHandler;

		// temporary objects
		std::vector<Edge*> m_tmpSpineEdges;
		std::vector<vector3> m_tmpFilPos;
		number m_tmpHeadHeight;

		// morphology params
		bool m_bFilAnisotropic;
		bool m_bWithRefinementProjector;

		number DENDRITE_LENGTH;
		number DENDRITE_RADIUS;
		number MEMBRANE_RADIUS;
		number m_membraneEnvelopeRadius;
		size_t m_dendriteRimVertices;

		number SPINE_NECK_RADIUS;
		number SPINE_NECK_LENGTH;
		number SPINE_HEAD_RADIUS;

		number FILAMENT_WIDTH;
		number FILAMENT_RIM_VERTICES;
		number m_filamentEnvelopeRadius;

		number EXTENSION_LENGTH;
		number EXTENSION_COMPARTMENT_LENGTH;

		number BOX_MARGIN;

		number SHORT_EDGES_FACTOR;

		size_t NECK_FILAMENTS;
		size_t NUM_FILAMENTS;

		number TETRAHEDRALIZATION_QUALITY;

		int INNER_SI;
		int OUTER_SI;
		int MEM_SI;
		int FIL_NECK_SI;
		int MEM_INNER_BND_SI;
		int MEM_OUTER_BND_SI;
		int MEM_NOFLUX_BND_SI;
		int SURF_CH_BND;
		int NOFLUX_BND;
		int DIRI_BND;
		int PSD_INNER_SI;
		int PSD_OUTER_SI;
		int EXT_LEFT_SI;
		int EXT_RIGHT_SI;
		int INTF_LEFT_CONSTRD_SI;
		int INTF_RIGHT_CONSTRD_SI;
		int INTF_LEFT_NODE1D_SI;
		int INTF_RIGHT_NODE1D_SI;
		int INTF_LEFT_NODEHD_SI;
		int INTF_RIGHT_NODEHD_SI;
		int EXT_LEFT_BND_SI;
		int EXT_RIGHT_BND_SI;
		int USELESS_SI;
		int INNER_FRONT_TMP_SI;
		int OUTER_FRONT_TMP_SI;
		int ENVELOPE_END_TMP_SI;
};



} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__MORPHO_GEN_H
