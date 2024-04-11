#pragma once
#include <vector>
#include "eigen/Eigen/Core"
#include "eigen/Eigen/Eigenvalues"
#include "eigen/Eigen/StdVector"
#include "eigen/Eigen/Dense"

struct RoadSide
{
	std::vector<Eigen::Vector2f> pl;
	std::vector<Eigen::Vector2f> pr;
};

struct MappingPoint
{
	Eigen::Vector2f original_point;
	Eigen::Vector2f project_point;
};

struct DoublePointSet
{
	std::vector<Eigen::Vector2f> hconvexhull_points;
	std::vector<Eigen::Vector2f> inside_peaks;
};

struct FAbMatrix
{
	Eigen::MatrixXf fa;
	Eigen::VectorXf fb;
};

struct StepResult
{
	bool circle_flag;
	Eigen::Vector2f node_query;
};

struct Line
{
	Eigen::Vector2f frontnode;
	Eigen::Vector2f backnode;
};
#pragma once
