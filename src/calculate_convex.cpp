#include"calculate_convex.h"
#include<iostream>

ConvexSolver::ConvexSolver() {
	_sconvex_node_set.resize(100);
	_hconvex_node_set.resize(1000);
}
//Get the dot aggregation point_left and point_right for the left and right sides of the lane
void ConvexSolver::calculate_both_side(std::vector<std::vector<float>>& data, RoadSide* const road_point)
{
	std::vector<Eigen::Vector2f> point_left, point_right;
	Eigen::Vector2f temp_point_left, temp_point_right;
	for (int i = 0; i < data.size(); i++)
	{
		temp_point_left[0] = data[i][0];
		temp_point_left[1] = data[i][1];
		temp_point_right[0] = data[i][2];
		temp_point_right[1] = data[i][3];
		point_left.push_back(temp_point_left);
		point_right.push_back(temp_point_right);
	}
	road_point->pl = point_left;
	road_point->pr = point_right;
}

//Get road centerline information
void ConvexSolver::calculate_node_mid(const std::vector<Eigen::Vector2f>& point_left,
	const std::vector<Eigen::Vector2f>& point_right,
	std::vector<Eigen::Vector2f>* const nodes_mid)
{
	Eigen::Vector2f node_mid;
	for (int i = 0; i < point_left.size(); i++)
	{
		node_mid = (point_left[i] + point_right[i]) / 2;
		nodes_mid->push_back(node_mid);
	}
}

//Segment the road point cloud
std::vector<int> ConvexSolver::subsection(const std::vector<Eigen::Vector2f>& point_left,
	const std::vector<Eigen::Vector2f>& nodes_mid)
{
	float sub_value1 = 0.3;
	float sub_value2 = 4.0;
	std::vector<int> sub_index;
	int i = 0;
	while (i < point_left.size() - 1)
	{
		float delta_Y = fabs(point_left[i][1] - point_left[i + 1][1]);
		float delta_y = fabs(nodes_mid[i][1] - nodes_mid[i + 1][1]);
		if (delta_Y > sub_value2 || delta_y > sub_value1)
		{
			sub_index.push_back(i + 1);
		}
		i++;
	}
	return sub_index;
}

//Determine the range of points to be projected, segment short sections, and expand them
void ConvexSolver::confirm_data_range(int border_left,
	int border_right,
	const Eigen::Vector2f& node_query,
	const std::vector<Eigen::Vector2f>& nodes_mid,
	const std::vector<Eigen::Vector2f>& point_left,
	const std::vector<Eigen::Vector2f>& point_right,
	std::vector<Eigen::Vector2f>* const base_data)
{
	float s1 = 3.0, s2 = 4.0;
	float d1;
	if (border_right != nodes_mid.size())
	{
		d1 = (nodes_mid[border_right] - nodes_mid[border_left]).norm();
	}
	else
	{
		d1 = (nodes_mid[border_right - 1] - nodes_mid[border_left]).norm();
	}

	if (d1 < s1)
	{
		for (int i = 0; i < nodes_mid.size(); i++)
		{
			float d2 = (nodes_mid[i] - node_query).norm();
			if (d2 < s2)
			{
				base_data->push_back(point_left[i]);
				base_data->push_back(point_right[i]);
			}
		}
	}
	else if (d1 < s1 >= 0 && border_right == nodes_mid.size())
	{
		for (int k = border_left; k < nodes_mid.size(); k++)
		{
			base_data->push_back(point_left[k]);
			base_data->push_back(point_right[k]);
		}
	}
	else
	{
		for (int j = border_left; j <= border_right; j++)
		{
			base_data->push_back(point_left[j]);
			base_data->push_back(point_right[j]);
		}
	}
}

//Sphere Flipping projection
void ConvexSolver::sphere_flipping(const std::vector<Eigen::Vector2f>& sub_data,
	const Eigen::Vector2f& node_query,
	std::vector<MappingPoint>* const mapping)
{
	float r = 25.0;
	MappingPoint temp_point;

	for (int j = 0; j < sub_data.size(); j++)
	{
		temp_point.project_point = node_query + (sub_data[j] - node_query) * ((2 * r -
			(sub_data[j] - node_query).norm()) / (sub_data[j] - node_query).norm());
		temp_point.original_point = sub_data[j];
		mapping->push_back(temp_point);
	}

}

//Compute the cross product of the vector
float convex_cross(Eigen::Vector2f& Pa, Eigen::Vector2f& Pb, Eigen::Vector2f& Pc)
{
	return (Pb[0] - Pa[0]) * (Pc[1] - Pa[1]) - (Pc[0] - Pa[0]) * (Pb[1] - Pa[1]);
}


//Fast convex hull algorithm and output the parent node of the mapping point
void ConvexSolver::graham_scan_algorithm(std::vector<MappingPoint>& mapping,
	std::vector<Eigen::Vector2f>* const skonvexhull_points)
{

	std::vector<Eigen::Vector2f> KonvexHullPoints;
	std::sort(mapping.begin(), mapping.end(), [&](MappingPoint& Pa, MappingPoint& Pb) {
		if (Pa.project_point[1] == Pb.project_point[1])
		return Pa.project_point[0] < Pb.original_point[0];
		else if (Pa.project_point[1] == Pb.project_point[1] && Pa.project_point[0] == Pb.original_point[0])
			return false;
		else
			return Pa.project_point[1] < Pb.project_point[1]; });
	_sconvex_node_set[0].project_point = mapping[0].project_point;
	_sconvex_node_set[0].original_point = mapping[0].original_point;

	std::sort(mapping.begin(), mapping.end(), [this](MappingPoint& Pa, MappingPoint& Pb) {
		float m = convex_cross(_sconvex_node_set[0].project_point, Pa.project_point, Pb.project_point);
	if (m > 0) {
		return true;
	}
	else if (m == 0 && (_sconvex_node_set[0].project_point - Pa.project_point).norm() -
		(_sconvex_node_set[0].project_point - Pb.project_point).norm() <= 0) {
		return true;
	}
	else {
		return false;
	}
		});
	_sconvex_node_set[1].project_point = mapping[1].project_point;
	_sconvex_node_set[1].original_point = mapping[1].original_point;

	int top = 1;
	for (int i = 2; i < mapping.size(); i++)
	{
		while (convex_cross(_sconvex_node_set[top - 1].project_point, _sconvex_node_set[top].project_point, mapping[i].project_point) < 0)
		{
			top--;		//Turn right, the middle point does not meet the requirement to eliminate
		}
		_sconvex_node_set[++top].project_point = mapping[i].project_point;
		_sconvex_node_set[top].original_point = mapping[i].original_point;
	}
	for (int j = 0; j <= top; j++)
	{
		skonvexhull_points->push_back(_sconvex_node_set[j].original_point);
	}
}


//A new convex hull H is obtained from the vertex of the star convex hull
void ConvexSolver::hgraham_scan_algorithm(std::vector<Eigen::Vector2f>& skonvexhull_points,
	DoublePointSet* const double_point_set)
{
	std::vector<Eigen::Vector2f> hconvexhull_point, inside_peak;
	std::sort(skonvexhull_points.begin(), skonvexhull_points.end(), [&](Eigen::Vector2f& p1, Eigen::Vector2f& p2) {
		if (p1[1] == p2[1])
		return p1[0] < p2[0];
		else
			return p1[1] < p2[1]; });

	_hconvex_node_set[0] = skonvexhull_points[0];
	std::sort(skonvexhull_points.begin() + 1, skonvexhull_points.end(), [this](Eigen::Vector2f& p1, Eigen::Vector2f& p2) {
		float n = convex_cross(_hconvex_node_set[0], p1, p2);
	if (n > 0) {
		return true;
	}
	else if (n == 0 && (_hconvex_node_set[0] - p1).norm() - (_hconvex_node_set[0] - p2).norm() <= 0) {
		return true;
	}
	else {
		return false;
	}});
	_hconvex_node_set[1] = skonvexhull_points[1];

	int top = 1;
	for (int i = 2; i < skonvexhull_points.size(); i++)
	{
		while (convex_cross(_hconvex_node_set[top - 1], _hconvex_node_set[top], skonvexhull_points[i]) < 0)
		{
			inside_peak.push_back(_hconvex_node_set[top]);
			top--;
		}
		_hconvex_node_set[++top] = skonvexhull_points[i];
	}
	for (int j = 0; j <= top; j++)
	{
		hconvexhull_point.push_back(_hconvex_node_set[j]);
	}
	double_point_set->hconvexhull_points = hconvexhull_point;
	double_point_set->inside_peaks = inside_peak;
}


//Compress convex hull H to output matrix information after compression
float h_cross(Eigen::Vector2f& p1, Eigen::Vector2f& p2)
{
	float x1, x2, y1, y2;
	x1 = p1[0];
	y1 = p1[1];
	x2 = p2[0];
	y2 = p2[1];
	return x1 * y2 - x2 * y1;
}

FAbMatrix ConvexSolver::compress_hulls(const DoublePointSet& hk_inside_points,
	const Eigen::Vector2f& node_query)
{
	//边缘集合
	std::vector<Eigen::Vector2f> hconvexhull_points = hk_inside_points.hconvexhull_points;
	std::vector<Eigen::Vector2f> inside_peaks = hk_inside_points.inside_peaks;
	Line line;
	std::vector<Line> line_set;
	for (int i = 0; i < hconvexhull_points.size() - 1; i++)
	{
		line.frontnode = hconvexhull_points[i];
		line.backnode = hconvexhull_points[i + 1];
		line_set.push_back(line);

		if (i == hconvexhull_points.size() - 2)
		{
			line.frontnode = hconvexhull_points[i + 1];
			line.backnode = hconvexhull_points[0];
			line_set.push_back(line);
		}
	}

	//Edge normal vector set
	std::vector<Eigen::Vector2f> norm_vector_set;
	for (int j = 0; j < line_set.size(); j++)
	{
		Eigen::Vector2f aq, ab, temp_vector, norm_vector;
		aq = node_query - line_set[j].frontnode;
		ab = line_set[j].backnode - line_set[j].frontnode;
		temp_vector = aq - (ab / ab.norm()) * (ab.dot(aq)) / ab.norm();
		norm_vector = temp_vector / temp_vector.norm();
		norm_vector_set.push_back(norm_vector);
	}

	//The triangle is defined, the convex polygon is segmented, the convex hull is segmented into the triangle set, 
	//and the position relation between the point and the triangle is judged
	std::vector<std::vector<float>> distance(line_set.size(), std::vector<float>(1, 0));
	for (int i = 0; i < inside_peaks.size(); i++)
	{
		for (int j = 0; j < line_set.size(); j++)
		{
			Eigen::Vector2f pa = node_query - inside_peaks[i];
			Eigen::Vector2f pb = line_set[j].frontnode - inside_peaks[i];
			Eigen::Vector2f pc = line_set[j].backnode - inside_peaks[i];
			float t1 = h_cross(pa, pb);
			float t2 = h_cross(pb, pc);
			float t3 = h_cross(pc, pa);
			if ((t1 > 0 && t2 > 0 && t3 > 0) || (t1 < 0 && t2 < 0 && t3 < 0))	//The vertices are inside the triangle
			{
				float temp_distance = (inside_peaks[i] - line_set[j].frontnode).dot(norm_vector_set[j]);
				distance[j].push_back(temp_distance);
				break;
			}
		}
	}

	int length_distance = distance.size();
	while (length_distance != line_set.size())
	{
		std::vector<float> a(4, 0.0);
		distance.push_back(a);
		length_distance = distance.size();
	}

	//Compression vector
	std::vector<Eigen::Vector2f> move_vector_set;
	for (int j = 0; j < distance.size(); j++)
	{
		std::sort(distance[j].begin(), distance[j].end(), [&](float& a, float& b) {return a > b; });
		float max_distance = distance[j][0];
		Eigen::Vector2f move_vector = norm_vector_set[j] * max_distance;
		move_vector_set.push_back(move_vector);
	}

	//The compressed linear representation matrix fax = fb;
	std::vector<Eigen::RowVector4f> compress_line_set(line_set.size(), Eigen::RowVector4f());
	for (int i = 0; i < line_set.size(); i++)
	{
		compress_line_set[i](0) = line_set[i].frontnode[0] + move_vector_set[i][0];
		compress_line_set[i](1) = line_set[i].frontnode[1] + move_vector_set[i][1];
		compress_line_set[i](2) = line_set[i].backnode[0] + move_vector_set[i][0];
		compress_line_set[i](3) = line_set[i].backnode[1] + move_vector_set[i][1];
	}

	FAbMatrix fab;
	Eigen::Matrix<float, Eigen::Dynamic, 2> a_matrix(compress_line_set.size(), 2);
	Eigen::VectorXf b_matrix(compress_line_set.size(), 1);
	for (int j = 0; j < compress_line_set.size(); j++)
	{
		float x1 = compress_line_set[j](0);
		float y1 = compress_line_set[j](1);
		float x2 = compress_line_set[j](2);
		float y2 = compress_line_set[j](3);
		a_matrix(j, 0) = y1 - y2;
		a_matrix(j, 1) = -(x1 - x2);
		b_matrix(j, 0) = ((y1 - y2) * x1 - (x1 - x2) * y1);
	}
	fab.fa = a_matrix;
	fab.fb = b_matrix;
	std::cout << "------------------------" << std::endl;

	std::cout << a_matrix.matrix() << std::endl;
	std::cout << b_matrix.matrix() << std::endl;

	std::cout << "------------------------" << std::endl;
	return fab;
}


//Take the next step to extend the convex hull
//分段道路凸包扩展全覆盖;判断Checknode是否在凸包内
bool check_in(FAbMatrix& fab_matrix, Eigen::Vector2f& node_query, Eigen::Vector2f& check_node)
{
	bool flag;
	Eigen::VectorXf node_result = fab_matrix.fa * node_query - fab_matrix.fb;
	Eigen::VectorXf test_result = fab_matrix.fa * check_node - fab_matrix.fb;

	for (int i = 0; i < test_result.size(); i++)
	{
		if (fabs(test_result[i]) < 1e-10)
		{
			test_result[i] = 0;
		}
	}
	Eigen::VectorXf m = node_result.cwiseProduct(test_result);
	if (m.minCoeff() >= 0) {
		flag = true;
	}
	else {
		flag = false;
	}
	return flag;
}
StepResult ConvexSolver::next_step(int border_left,
	int border_right,
	const std::vector<int>& sub_index,
	FAbMatrix& fab_matrix,
	Eigen::Vector2f& node_query,
	const std::vector<Eigen::Vector2f>& nodes_mid,
	int i)
{
	//Loop ends when check_node is inside the convex hull; There are less segmented road points, need to discuss to determine the check_node
	//If the segment is long, the penultimate point of the segment is taken as check_node
	Eigen::Vector2f check_node;
	StepResult step_result;
	step_result.circle_flag = false;
	if (border_right - border_left > 0 && border_right - border_left <= 3)
	{
		check_node = nodes_mid[border_right - 1];
	}
	else if (border_right - border_left == 0 && i != sub_index.size())
	{
		step_result.node_query = nodes_mid[sub_index[i] + 1];
		step_result.circle_flag = true;
	}
	else if (border_right - border_left == 0 && i == sub_index.size())
	{
		step_result.circle_flag = true;
	}
	else
	{
		check_node = nodes_mid[border_right - 2];
	}

	//check_node in convex hull, circleflag equal to true; The segment extension ends
	if (!step_result.circle_flag == true && i != sub_index.size())
	{
		if (check_in(fab_matrix, node_query, check_node))
		{
			step_result.node_query = nodes_mid[sub_index[i] + 1];
			step_result.circle_flag = true;
		}
	}
	else if (!step_result.circle_flag == true && i == sub_index.size())
	{
		if (check_in(fab_matrix, node_query, check_node))
		{
			step_result.circle_flag = true;
		}
	}

	//Checknode is not in the convex hull, and walk forward from the end point of the segment 
	//to determine the query coordinate node_query for the next extension to generate the convex hull
	Eigen::Vector2f jump_node;
	if (step_result.circle_flag == false)
	{
		if (border_right != nodes_mid.size())
			for (int k = border_right; k > border_left; k--)
			{
				jump_node = nodes_mid[k];
				if (check_in(fab_matrix, node_query, jump_node))
				{
					if (k + 1 < sub_index[i])
					{
						jump_node = nodes_mid[k + 1];
					}
					else
					{
						jump_node = nodes_mid[sub_index[i]];
					}
					step_result.node_query = jump_node;
					break;
				}
			}
		else
		{
			for (int j = border_right - 1; j > border_left; j--)
			{
				jump_node = nodes_mid[j];
				if (check_in(fab_matrix, node_query, jump_node))
				{
					if (j + 1 < nodes_mid.size())
					{
						jump_node = nodes_mid[j + 1];
					}
					else
					{
						jump_node = nodes_mid[nodes_mid.size()];
					}
					step_result.node_query = jump_node;
					break;
				}
			}

		}
	}
	return step_result;
}

//Convex hull is generated at the segment to prevent gaps
Eigen::Vector2f calculate_average_xy(const std::vector<Eigen::Vector2f>& space_data)
{
	Eigen::Vector2f center_node;
	float x, y;
	float sum_x = 0, sum_y = 0;
	int s = space_data.size();
	for (int i = 0; i < space_data.size(); i++)
	{
		sum_x += space_data[i][0];
		sum_y += space_data[i][1];
	}
	x = sum_x / s;
	y = sum_y / s;
	center_node[0] = x;
	center_node[1] = y;
	return center_node;
}

//Obtain the object matrix corresponding to the convex hull
FAbMatrix ConvexSolver::generate_convexhulls(std::vector<MappingPoint>& space_mapping, Eigen::Vector2f& center_node_query)
{
	std::vector<Eigen::Vector2f> space_skonvexhull_points;
	graham_scan_algorithm(space_mapping, &space_skonvexhull_points);

	DoublePointSet space_hk_inside_points;
	hgraham_scan_algorithm(space_skonvexhull_points, &space_hk_inside_points);

	FAbMatrix space_fab_matrix = compress_hulls(space_hk_inside_points, center_node_query);

	return space_fab_matrix;
}

void ConvexSolver::add_space(const std::vector<int>& sub_index,
	const std::vector<Eigen::Vector2f>& pl,
	const std::vector<Eigen::Vector2f>& pr,
	std::vector<FAbMatrix>* const fab_set)
{
	int number = 5;
	int border_left, border_right;
	for (int j = 0; j < sub_index.size(); j++)
	{												//很重要Space_data的作用域有效范围，若其在for循环之外，Space_data的数据一直在积累
		std::vector<Eigen::Vector2f> space_data;	//Space_data置于for循环之内，每次循环都会重置更新
		if (j == sub_index.size() - 1)
		{
			if (fabs(sub_index[j] - pl.size()) < number)
			{
				border_left = sub_index[j] - number;
				border_right = pl.size();
			}
			else
			{
				border_left = sub_index[j] - number;
				border_right = sub_index[j] + number;
			}
		}
		else
		{
			if (sub_index[j] < number)
			{
				border_left = 0;
				border_right = sub_index[j] + number;
			}
			else
			{
				border_left = sub_index[j] - number;
				border_right = sub_index[j] + number;
			}
		}

		//Gets Space data according to the defined range
		for (int i = border_left; i <= border_right; i++)
		{
			space_data.push_back(pl[i]);
			space_data.push_back(pr[i]);
		}

		//The node query takes the geometric center of the box fetching point set as the query coordinate point
		Eigen::Vector2f CenterNode_query = calculate_average_xy(space_data);

		//Mapping point set
		std::vector<MappingPoint> space_mapping;
		sphere_flipping(space_data, CenterNode_query, &space_mapping);

		//Debug
		std::cout << "##------------------------------##" << std::endl;
		std::cout << "j = " << j << std::endl;
		for (int t = 0; t < space_mapping.size(); t++)
		{
			std::cout << space_mapping[t].project_point[0] << "  " << space_mapping[t].project_point[1] << std::endl;
		}
		std::cout << "##------------------------------##" << std::endl;

		FAbMatrix space_fab_matrix = generate_convexhulls(space_mapping, CenterNode_query);

		fab_set->push_back(space_fab_matrix);
	}
}