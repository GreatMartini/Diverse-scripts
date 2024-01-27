#include "RefineClasses.h"
#include "TimeInterpolateOperator.h"
#include "TimeRefinementIntegrator.h"
#include "RefineSchedule.h"
#include "RefineAlgorithm.h"
#include "StandardRefineTransactionFactory.h"

#include "RefineTimeTransaction.h"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "MainRestartData.h"
#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"
#include "Commons.h"
#include "SlicerDataWriter.h"
#include "SphereDataWriter.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "IntegrateDataWriter.h"
#include "PointDataWriter.h"
#include "ExternalInitialData.h"
#include "TimeInterpolator.h"

using namespace std;
using namespace SAMRAI;

#define DIMENSIONS 3



class Problem : 
   public mesh::StandardTagAndInitStrategy,
   public xfer::RefinePatchStrategy,
   public xfer::CoarsenPatchStrategy,
   public algs::TimeRefinementLevelStrategy
{
public:
	/*
	 * Constructor of the problem.
	 */
	Problem(
		const string& object_name,
		const tbox::Dimension& dim,
		std::shared_ptr<tbox::Database>& input_db,
		std::shared_ptr<geom::CartesianGridGeometry >& grid_geom, 
	   	std::shared_ptr<hier::PatchHierarchy >& patch_hierarchy, 
	   	MainRestartData& mrd,
		const double dt, 
		const bool init_from_files,
		const int console_output,
		const int timer_output,
		const int mesh_output_period,
		const vector<string> full_mesh_writer_variables,
		std::shared_ptr<appu::VisItDataWriter>& data_writer,
		const vector<int> slicer_output_period,
		const vector<set<string> > sliceVariables,
		vector<std::shared_ptr<SlicerDataWriter > >sliceWriters,
		const vector<int> sphere_output_period,
		const vector<set<string> > sphereVariables,
		vector<std::shared_ptr<SphereDataWriter > > sphereWriters,
		const vector<bool> slicer_analysis_output,
		const vector<bool> sphere_analysis_output,
		const vector<int> integration_output_period,
		const vector<set<string> > integralVariables,
		vector<std::shared_ptr<IntegrateDataWriter > > integrateDataWriters,
		const vector<bool> integration_analysis_output,
		const vector<int> point_output_period,
		const vector<set<string> > pointVariables,
		vector<std::shared_ptr<PointDataWriter > > pointDataWriters,
		const vector<bool> point_analysis_output);
  
	/*
	 * Destructor.
	 */
	~Problem();

	/*
	 * Block of subcycling inherited methods.
	 */
	void initializeLevelIntegrator(
	      const std::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_alg);

	double getLevelDt(
	      const std::shared_ptr<hier::PatchLevel>& level,
	      const double dt_time,
	      const bool initial_time);

	double getMaxFinerLevelDt(
	      const int finer_level_number,
	      const double coarse_dt,
	      const hier::IntVector& ratio_to_coarser);

	double advanceLevel(
	      const std::shared_ptr<hier::PatchLevel>& level,
	      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
	      const double current_time,
	      const double new_time,
	      const bool first_step,
	      const bool last_step,
	      const bool regrid_advance = false);

	void standardLevelSynchronization(
	      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
	      const int coarsest_level,
	      const int finest_level,
	      const double sync_time,
	      const std::vector<double>& old_times);

	void synchronizeNewLevels(
	      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
	      const int coarsest_level,
	      const int finest_level,
	      const double sync_time,
	      const bool initial_time);

	void resetTimeDependentData(
	      const std::shared_ptr<hier::PatchLevel>& level,
	      const double new_time,
	      const bool can_be_refined);

	void resetDataToPreadvanceState(
	      const std::shared_ptr<hier::PatchLevel>& level);

	bool usingRefinedTimestepping() const
	{
	      return d_refinedTimeStepping;
	} 

	/*
	 * Register the current iteration in the class.
	 */
	void registerIteration(int iter) {
		simPlat_iteration = iter;
	}

	/*
	 * Initialize the data from a given level.
   	 */
	virtual void initializeLevelData(
		const std::shared_ptr<hier::PatchHierarchy >& hierarchy ,
		const int level_number ,
		const double init_data_time ,
		const bool can_be_refined ,
		const bool initial_time ,
		const std::shared_ptr<hier::PatchLevel >& old_level=std::shared_ptr<hier::PatchLevel>() ,
		const bool allocate_data = true);
	void postInit();

	/*
	 * Reset the hierarchy-dependent internal information.
	 */
	virtual void resetHierarchyConfiguration(
		const std::shared_ptr<hier::PatchHierarchy >& new_hierarchy ,
		int coarsest_level ,
		int finest_level);

	/*
	 * Checks the finalization conditions              
	 */
	bool checkFinalization(
		const double simPlat_time, 
		const double simPlat_dt);

	/*
	 * This method sets the physical boundary conditions.
	*/
	void setPhysicalBoundaryConditions(
		hier::Patch& patch,
		const double fill_time,
		const hier::IntVector& ghost_width_to_fill);
	/*
	 * Set up external plotter to plot internal data from this class.        
	 * Tell the plotter about the refinement ratios.  Register variables     
	 * appropriate for plotting.                                            
	 */
	int setupPlotterMesh(appu::VisItDataWriter &plotter ) const;
	int setupSlicePlotter(vector<std::shared_ptr<SlicerDataWriter> > &plotters) const;
	int setupSpherePlotter(vector<std::shared_ptr<SphereDataWriter> > &plotters) const;
	int setupIntegralPlotter(vector<std::shared_ptr<IntegrateDataWriter> > &plotters) const;
	int setupPointPlotter(vector<std::shared_ptr<PointDataWriter> > &plotters) const;


	/*
	 * Map data on a patch. This mapping is done only at the begining of the simulation.
	 */
	void mapDataOnPatch(const double time, const bool initial_time, const int ln, const std::shared_ptr< hier::PatchLevel >& level);

	/*
	 * Sets the limit for the checkstencil routine
	 */
	void setStencilLimits(std::shared_ptr< hier::Patch > patch, int i, int j, int k, int v) const;

	/*
	 * Checks if the point has a stencil width
	 */
	void checkStencil(std::shared_ptr< hier::Patch > patch, int i, int j, int k, int v) const;

	/*
	 * Flood-Fill algorithm
	 */
	void floodfill(std::shared_ptr< hier::Patch > patch, int i, int j, int k, int pred, int seg) const;


	/*
	 * FOV correction for AMR
	 */
	void correctFOVS(const std::shared_ptr< hier::PatchLevel >& level);



	/*
	 * Interphase mapping. Calculates the FOV and its variables.
	 */
	void interphaseMapping(const double time, const bool initial_time, const int ln, const std::shared_ptr< hier::PatchLevel >& level, const int remesh);




	/*
	 * Initialize data on a patch. This initialization is done only at the begining of the simulation.
	 */
	void initializeDataOnPatch(
		hier::Patch& patch,
		const double time,
       		const bool initial_time);

	/*
	 *  Cell tagging routine - tag cells that require refinement based on a provided condition. 
	 */
	void applyGradientDetector(
	   	const std::shared_ptr< hier::PatchHierarchy >& hierarchy, 
	   	const int level_number,
	   	const double time, 
	   	const int tag_index,
	   	const bool initial_time,
	   	const bool uses_richardson_extrapolation_too);

	/*
	* Return maximum stencil width needed for user-defined
	* data interpolation operations.  Default is to return
	* zero, assuming no user-defined operations provided.
	*/
	hier::IntVector getRefineOpStencilWidth(const tbox::Dimension &dim) const
	{
		return hier::IntVector::getZero(dim);
	}

	/*
	* Pre- and post-processing routines for implementing user-defined
	* spatial interpolation routines applied to variables.  The 
	* interpolation routines are used in the MOL AMR algorithm
	* for filling patch ghost cells before advancing data on a level
	* and after regridding a level to fill portions of the new level
	* from some coarser level.  These routines are called automatically
	* from within patch boundary filling schedules; thus, some concrete
	* function matching these signatures must be provided in the user's
	* patch model.  However, the routines only need to perform some 
	* operations when "USER_DEFINED_REFINE" is given as the interpolation 
	* method for some variable when the patch model registers variables
	* with the MOL integration algorithm, typically.  If the 
	* user does not provide operations that refine such variables in either 
	* of these routines, then they will not be refined.
	*
	* The order in which these operations are used in each patch 
	* boundary filling schedule is:
	* 
	* - \b (1) {Call user's preprocessRefine() routine.}
	* - \b (2) {Refine all variables with standard interpolation operators.}
	* - \b (3) {Call user's postprocessRefine() routine.}
	* 
	* 
	* Also, user routines that implement these functions must use 
	* data corresponding to the d_scratch context on both coarse and
	* fine patches.
	*/
	virtual void preprocessRefine(
		hier::Patch& fine,
		const hier::Patch& coarse,
                const hier::Box& fine_box,
                const hier::IntVector& ratio)
	{
		NULL_USE(fine);
		NULL_USE(coarse);
		NULL_USE(fine_box);
		NULL_USE(ratio);
	}
	virtual void postprocessRefine(
		hier::Patch& fine,
                const hier::Patch& coarse,
                const hier::Box& fine_box,
                const hier::IntVector& ratio)
	{
		NULL_USE(fine);
		NULL_USE(coarse);
		NULL_USE(fine_box);
		NULL_USE(ratio);
	}



	/*
	* Return maximum stencil width needed for user-defined
	* data coarsen operations.  Default is to return
	* zero, assuming no user-defined operations provided.
	*/
	hier::IntVector getCoarsenOpStencilWidth( const tbox::Dimension &dim ) const {
		return hier::IntVector::getZero(dim);
	}

	/*
	* Pre- and post-processing routines for implementing user-defined
	* spatial coarsening routines applied to variables.  The coarsening 
	* routines are used in the MOL AMR algorithm synchronizing 
	* coarse and fine levels when they have been integrated to the same
	* point.  These routines are called automatically from within the 
	* data synchronization coarsen schedules; thus, some concrete
	* function matching these signatures must be provided in the user's
	* patch model.  However, the routines only need to perform some
	* operations when "USER_DEFINED_COARSEN" is given as the coarsening
	* method for some variable when the patch model registers variables
	* with the MOL level integration algorithm, typically.  If the
	* user does not provide operations that coarsen such variables in either
	* of these routines, then they will not be coarsened.
	*
	* The order in which these operations are used in each coarsening
	* schedule is:
	* 
	* - \b (1) {Call user's preprocessCoarsen() routine.}
	* - \b (2) {Coarsen all variables with standard coarsening operators.}
	* - \b (3) {Call user's postprocessCoarsen() routine.}
	* 
	*
	* Also, user routines that implement these functions must use
	* corresponding to the d_new context on both coarse and fine patches
	* for time-dependent quantities.
	*/
	virtual void preprocessCoarsen(
		hier::Patch& coarse,
                const hier::Patch& fine,
                const hier::Box& coarse_box,
                const hier::IntVector& ratio)
	{
		NULL_USE(fine);
		NULL_USE(coarse);
		NULL_USE(coarse_box);
		NULL_USE(ratio);
	}
	virtual void postprocessCoarsen(
		hier::Patch& coarse,
                const hier::Patch& fine,
                const hier::Box& coarse_box,
                const hier::IntVector& ratio)
	{
		NULL_USE(fine);
		NULL_USE(coarse);
		NULL_USE(coarse_box);
		NULL_USE(ratio);
	}

	/*
	 * Computes the dt to be used
	 */
	double computeDt() const;



	/*
	 * Checks if the point has to be stalled
	 */
	bool checkStalled(std::shared_ptr< hier::Patch > patch, int i, int j, int k, int v) const;



	/*
	 * Saves restart information relevant in Problem class.
	 */
	void putToRestart(MainRestartData& mrd);

	/*
	 * Loads restart information relevant in Problem class.
	 */
	void getFromRestart(MainRestartData& mrd);

	/*
	 * Allocates internal variable memory
	 */
	void allocateAfterRestart();

private:	 
	//Variables for the refine and coarsen algorithms

	std::shared_ptr< xfer::RefineAlgorithm > d_bdry_fill_init;
	std::shared_ptr< xfer::RefineAlgorithm > d_bdry_post_coarsen;
	std::vector< std::shared_ptr< xfer::RefineSchedule > > d_bdry_sched_postCoarsen;
	std::shared_ptr<TimeInterpolator> time_interpolate_operator_mesh1;
	std::shared_ptr< xfer::RefineAlgorithm > d_bdry_fill_advance1;
	std::vector< std::shared_ptr< xfer::RefineSchedule > > d_bdry_sched_advance1;
	std::shared_ptr< xfer::RefineAlgorithm > d_bdry_fill_advance6;
	std::vector< std::shared_ptr< xfer::RefineSchedule > > d_bdry_sched_advance6;
	std::shared_ptr< xfer::RefineAlgorithm > d_bdry_fill_advance11;
	std::vector< std::shared_ptr< xfer::RefineSchedule > > d_bdry_sched_advance11;
	std::shared_ptr< xfer::RefineAlgorithm > d_bdry_fill_advance16;
	std::vector< std::shared_ptr< xfer::RefineSchedule > > d_bdry_sched_advance16;
	std::shared_ptr< xfer::RefineAlgorithm > d_bdry_fill_analysis1;
	std::vector< std::shared_ptr< xfer::RefineSchedule > > d_bdry_sched_analysis1;

	std::shared_ptr< xfer::RefineAlgorithm > d_mapping_fill;
	std::shared_ptr< xfer::RefineAlgorithm > d_tagging_fill;
	std::shared_ptr< xfer::CoarsenAlgorithm > d_coarsen_algorithm;
	std::vector< std::shared_ptr< xfer::CoarsenSchedule > > d_coarsen_schedule;

	std::shared_ptr< xfer::RefineAlgorithm > d_fill_new_level;

	//Object name
	std::string d_object_name;

	//Pointers to the grid geometry and the patch hierarchy
   	std::shared_ptr<geom::CartesianGridGeometry > d_grid_geometry;
	std::shared_ptr<hier::PatchHierarchy > d_patch_hierarchy;

	//Identifiers of the fields and auxiliary fields
	int d_Phi_id, d_gammac_xx_id, d_gammac_xy_id, d_gammac_xz_id, d_gammac_yx_id, d_gammac_yy_id, d_gammac_yz_id, d_gammac_zx_id, d_gammac_zy_id, d_gammac_zz_id, d_A_xx_id, d_A_xy_id, d_A_xz_id, d_A_yy_id, d_A_yz_id, d_A_zz_id, d_trK_id, d_alpha_id, d_beta_x_id, d_beta_y_id, d_beta_z_id, d_Gam_x_id, d_Gam_y_id, d_Gam_z_id, d_Rscalar_id, d_HamCon_id, d_MomCon_x_id, d_MomCon_y_id, d_MomCon_z_id, d_rk1gammac_xx_id, d_rk1gammac_xy_id, d_rk1gammac_xz_id, d_rk1gammac_yy_id, d_rk1gammac_yx_id, d_rk1gammac_yz_id, d_rk1gammac_zz_id, d_rk1gammac_zx_id, d_rk1gammac_zy_id, d_rk1A_xx_id, d_rk1A_xy_id, d_rk1A_xz_id, d_rk1A_yy_id, d_rk1A_yz_id, d_rk1A_zz_id, d_rk1Phi_id, d_rk1trK_id, d_rk1beta_x_id, d_rk1beta_y_id, d_rk1beta_z_id, d_rk1alpha_id, d_rk1Gam_x_id, d_rk1Gam_y_id, d_rk1Gam_z_id, d_rk2gammac_xx_id, d_rk2gammac_xy_id, d_rk2gammac_xz_id, d_rk2gammac_yy_id, d_rk2gammac_yx_id, d_rk2gammac_yz_id, d_rk2gammac_zz_id, d_rk2gammac_zx_id, d_rk2gammac_zy_id, d_rk2A_xx_id, d_rk2A_xy_id, d_rk2A_xz_id, d_rk2A_yy_id, d_rk2A_yz_id, d_rk2A_zz_id, d_rk2Phi_id, d_rk2trK_id, d_rk2beta_x_id, d_rk2beta_y_id, d_rk2beta_z_id, d_rk2alpha_id, d_rk2Gam_x_id, d_rk2Gam_y_id, d_rk2Gam_z_id, d_rk3gammac_xx_id, d_rk3gammac_xy_id, d_rk3gammac_xz_id, d_rk3gammac_yy_id, d_rk3gammac_yx_id, d_rk3gammac_yz_id, d_rk3gammac_zz_id, d_rk3gammac_zx_id, d_rk3gammac_zy_id, d_rk3A_xx_id, d_rk3A_xy_id, d_rk3A_xz_id, d_rk3A_yy_id, d_rk3A_yz_id, d_rk3A_zz_id, d_rk3Phi_id, d_rk3trK_id, d_rk3beta_x_id, d_rk3beta_y_id, d_rk3beta_z_id, d_rk3alpha_id, d_rk3Gam_x_id, d_rk3Gam_y_id, d_rk3Gam_z_id, d_d_i_trK_beta_x_beta_y_beta_z_alpha_Gam_x_Gam_y_Gam_z_gammac_xx_gammac_xy_gammac_xz_gammac_yy_gammac_yx_gammac_yz_gammac_zz_gammac_zx_gammac_zy_A_xx_A_xy_A_xz_A_yy_A_yz_A_zz_Phi_id, d_d_j_trK_beta_x_beta_y_beta_z_alpha_Gam_x_Gam_y_Gam_z_gammac_xx_gammac_xy_gammac_xz_gammac_yy_gammac_yx_gammac_yz_gammac_zz_gammac_zx_gammac_zy_A_xx_A_xy_A_xz_A_yy_A_yz_A_zz_Phi_id, d_d_k_trK_beta_x_beta_y_beta_z_alpha_Gam_x_Gam_y_Gam_z_gammac_xx_gammac_xy_gammac_xz_gammac_yy_gammac_yx_gammac_yz_gammac_zz_gammac_zx_gammac_zy_A_xx_A_xy_A_xz_A_yy_A_yz_A_zz_Phi_id, d_stalled_1_id, d_Phi_p_id, d_gammac_xx_p_id, d_gammac_xy_p_id, d_gammac_xz_p_id, d_gammac_yx_p_id, d_gammac_yy_p_id, d_gammac_yz_p_id, d_gammac_zx_p_id, d_gammac_zy_p_id, d_gammac_zz_p_id, d_A_xx_p_id, d_A_xy_p_id, d_A_xz_p_id, d_A_yy_p_id, d_A_yz_p_id, d_A_zz_p_id, d_trK_p_id, d_alpha_p_id, d_beta_x_p_id, d_beta_y_p_id, d_beta_z_p_id, d_Gam_x_p_id, d_Gam_y_p_id, d_Gam_z_p_id;

	//Parameter variables
	double dissipation_factor_A_xx;
	double dissipation_factor_A_zz;
	double dissipation_factor_A_xy;
	double dissipation_factor_A_xz;
	double dissipation_factor_Phi;
	double yc;
	double dissipation_factor_trK;
	double lambda_z;
	double tend;
	double dissipation_factor_gammac_zx;
	double dissipation_factor_gammac_xx;
	double dissipation_factor_gammac_zz;
	double dissipation_factor_gammac_zy;
	double dissipation_factor_gammac_xz;
	double dissipation_factor_gammac_xy;
	double dissipation_factor_A_yy;
	double dissipation_factor_A_yz;
	double xc;
	double dissipation_factor_Gam_z;
	double dissipation_factor_Gam_y;
	double dissipation_factor_Gam_x;
	double a0;
	double rc;
	double dissipation_factor_beta_x;
	double dissipation_factor_beta_y;
	double lambda_r;
	double dissipation_factor_beta_z;
	double dissipation_factor_alpha;
	double dissipation_factor_gammac_yy;
	double dissipation_factor_gammac_yx;
	double dissipation_factor_gammac_yz;


	//mapping fields
	int d_nonSync_regridding_tag_id, d_interior_regridding_value_id, d_interior_i_id, d_interior_j_id, d_interior_k_id, d_FOV_1_id, d_FOV_xLower_id, d_FOV_xUpper_id, d_FOV_yLower_id, d_FOV_yUpper_id, d_FOV_zLower_id, d_FOV_zUpper_id;

	//Stencils of the discretization method variable
	int d_ghost_width, d_regionMinThickness;

	//initial dt
	double initial_dt;

   	const tbox::Dimension d_dim;

	//Subcycling variable
   	bool d_refinedTimeStepping;

	//Current iteration
	int simPlat_iteration;

	//Initialization from files
	bool d_init_from_restart;

	//List of particle variables


	//Variables to dump
	vector<set<string> > d_sliceVariables;
	vector<set<string> > d_sphereVariables;

	vector<set<string> > d_integralVariables;
	int d_mask_id;
	vector<set<string> > d_pointVariables;


	//FileWriter
	vector<std::shared_ptr<SlicerDataWriter > > d_sliceWriters;
	vector<std::shared_ptr<SphereDataWriter > > d_sphereWriters;
	vector<int> d_slicer_output_period;
	vector<int> d_sphere_output_period;
	vector<int> next_slice_dump_iteration;
	vector<int> next_sphere_dump_iteration;
	vector<bool> analysis_slice_dump;
	vector<bool> analysis_sphere_dump;
	int viz_mesh_dump_interval;
	int next_mesh_dump_iteration;
	set<string> d_full_mesh_writer_variables;
	std::shared_ptr<appu::VisItDataWriter > d_visit_data_writer;
	vector<std::shared_ptr<IntegrateDataWriter > > d_integrateDataWriters;
	vector<int> d_integration_output_period;
	vector<int> next_integration_dump_iteration;
	vector<std::shared_ptr<PointDataWriter > > d_pointDataWriters;
	vector<int> d_point_output_period;
	vector<int> next_point_dump_iteration;
	vector<bool> analysis_integration_dump;
	vector<bool> analysis_point_dump;


	vector<int> current_iteration;
	vector<int> bo_substep_iteration;

	//Console output variables
	int d_output_interval;
	int next_console_output;
	int d_timer_output_interval;
	int next_timer_output;

	//regridding options
	std::shared_ptr<tbox::Database> regridding_db;
	bool d_regridding;
	std::string d_regridding_field, d_regridding_type, d_regridding_field_shadow;
	double d_regridding_threshold, d_regridding_compressionFactor, d_regridding_mOffset, d_regridding_error, d_regridding_buffer;
	int d_regridding_min_level, d_regridding_max_level;

	//Gets the coarser patch that includes the box.
	const std::shared_ptr<hier::Patch >& getCoarserPatch(
		const std::shared_ptr< hier::PatchLevel >& level,
		const hier::Box interior, 
		const hier::IntVector ratio);

	static bool Equals(double d1, double d2);
	static inline int GetExpoBase2(double d);
};


