#pragma once
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include"struct_data.h"

class ConvexSolver
{
public:
	ConvexSolver();
	~ConvexSolver() {}

public:
	void calculate_both_side(std::vector<std::vector<float>>& data, RoadSide* const road_lr);

	void calculate_node_mid(const std::vector<Eigen::Vector2f>& pl,
		const std::vector<Eigen::Vector2f>& pr,
		std::vector<Eigen::Vector2f>* const nodes_mid);

	std::vector<int> subsection(const std::vector<Eigen::Vector2f>& pl,
		const std::vector<Eigen::Vector2f>& nodes_mid);

	void confirm_data_range(int border_left,
		int border_right,
		const Eigen::Vector2f& node_query,
		const std::vector<Eigen::Vector2f>& nodes_mid,
		const std::vector<Eigen::Vector2f>& pl,
		const std::vector<Eigen::Vector2f>& pr,
		std::vector<Eigen::Vector2f>* const base_data);

	void sphere_flipping(const std::vector<Eigen::Vector2f>& sub_data,
		const Eigen::Vector2f& node_query,
		std::vector<MappingPoint>* const mapping);

	void graham_scan_algorithm(std::vector<MappingPoint>& mapping,
		std::vector<Eigen::Vector2f>* const skonvexhull_points);

	void hgraham_scan_algorithm(std::vector<Eigen::Vector2f>& skonvexhull_points,
		DoublePointSet* const double_point_set);

	FAbMatrix compress_hulls(const DoublePointSet& hk_Inside_points,
		const Eigen::Vector2f& node_query);

	StepResult next_step(int border_left,
		int border_right,
		const std::vector<int>& sub_index,
		FAbMatrix& fab_matrix,
		Eigen::Vector2f& node_query,
		const std::vector<Eigen::Vector2f>& nodes_mid,
		int i);

	FAbMatrix generate_convexhulls(std::vector<MappingPoint>& space_mapping, Eigen::Vector2f& center_node_query);

	void add_space(const std::vector<int>& sub_index,
		const std::vector<Eigen::Vector2f>& pl,
		const std::vector<Eigen::Vector2f>& pr,
		std::vector<FAbMatrix>* const fab_set);
private:
	std::vector<MappingPoint> _sconvex_node_set;
	std::vector<Eigen::Vector2f> _hconvex_node_set;
};
