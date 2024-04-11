#include"road_data.h"
#include"calculate_convex.h"
#include <iostream>

int main() {
	//Load road point cloud data
	std::vector<std::vector<float>> data = load_data();

	//Filter the initial data
	std::vector<std::vector<float>>::iterator itend = std::unique(data.begin(), data.end());
	data.erase(itend, data.end());    

	//Print road data information
	print_vector(data);

	ConvexSolver convex_solver;

	//Get road centerline information
	RoadSide road_left_right;
	convex_solver.calculate_both_side(data, &road_left_right);
	std::vector<Eigen::Vector2f> point_left = road_left_right.pl;
	std::vector<Eigen::Vector2f> point_right = road_left_right.pr;
	std::vector<Eigen::Vector2f> nodes_mid;
	convex_solver.calculate_node_mid(point_left, point_right, &nodes_mid);

	//Segment the road point cloud
	std::vector<int> sub_index = convex_solver.subsection(point_left, nodes_mid);

	//Initializes the query coordinates : node_query
	int index = 1;
	Eigen::Vector2f node_query = nodes_mid[index];

	//Store the road segment matrix set and the segment point matrix set
	std::vector<FAbMatrix> fab_set;

	for (int i = 0; i < sub_index.size(); i++)
	{
		int border_left, border_right;
		while (true)
		{
			if (i == 0)
			{
				border_left = 0;
				border_right = sub_index[i] - 1;
			}
			else if (i < sub_index.size() && i>0)
			{
				border_left = sub_index[i - 1];
				border_right = sub_index[i] - 1;
			}
			else
			{
				border_left = sub_index[i - 1];
				border_right = data.size();
			}

			//Determine the range of points to be projected
			std::vector<Eigen::Vector2f> sub_data;
			convex_solver.confirm_data_range(border_left, border_right, node_query, nodes_mid, point_left, point_right, &sub_data);

			//Sphere Flipping projection
			std::vector<MappingPoint> mapping;
			convex_solver.sphere_flipping(sub_data, node_query, &mapping);

			//Debug
			std::cout << "--------------------------------------------" << std::endl;
			std::cout << "i = " << i << std::endl;
			for (int t = 0; t < mapping.size(); t++)
			{
				std::cout << mapping[t].project_point[0] << "  " << mapping[t].project_point[1] << std::endl;
			}
			std::cout << "--------------------------------------------" << std::endl;

			/*The fast convex hull algorithm is used to calculate the star convex hull S 
			    by reflectionand output the set of star convex hull points*/
			std::vector<Eigen::Vector2f> skonvexhull_points;
			convex_solver.graham_scan_algorithm(mapping, &skonvexhull_points);

			//A new convex hull H is obtained from the vertex of the star convex hull
			DoublePointSet double_point_set;
			convex_solver.hgraham_scan_algorithm(skonvexhull_points, &double_point_set);

			//Compress convex hull H to output matrix information after compression
			FAbMatrix fab_matrix = convex_solver.compress_hulls(double_point_set, node_query);

			//Stores target matrix information
			fab_set.push_back(fab_matrix);

			//Full coverage of segmented road convex hull extension
			StepResult step_result = convex_solver.next_step(border_left, border_right, sub_index, fab_matrix, node_query, nodes_mid, i);

			node_query = step_result.node_query;
			if (step_result.circle_flag)
				break;
		}
	}

	//Convex hull is generated at the segment to prevent gaps
	convex_solver.add_space(sub_index, point_left, point_right, &fab_set);

	return 0;
}