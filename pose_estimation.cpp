// pose_estimation.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "get_pose.h"
#include "Eigen\Dense"
#include <unsupported/Eigen/Polynomials>

using namespace cv;
using namespace std;

int threeLines(cv::Mat& R, cv::Mat& tvec, const double& f, const cv::Mat& olines, const cv::Mat& ilines)
{
	//如果直线数不等于3，返回-1
	/*if (olines.size() < 3 || ilines.size() != 3)
		return(-1);*/
	//Mat数据不是double型，返回-2
	if (ilines.depth() !=  CV_64F)
		return(-2);

	//保存olines和ilines中直线参数
	Eigen::Vector3d V0[3], P0[3], v0[3], p0[3], V1[3];
	for (int j = 0; j < 3; j++){
		const double* data_i = ilines.ptr<double>(j);
		const double* data_o = olines.ptr<double>(j);
		
		for (int i = 0; i < 3; i++){
			V0[j](i) = data_o[i];
			v0[j](i) = data_i[i];
		}
		for (int i = 3; i < 6; i++){
			P0[j](i - 3) = data_o[i];
			p0[j](i - 3) = data_i[i];
		}
	}
	double d01 = v0[0](0)*p0[0](1) - v0[0](1)*p0[0](0);
	double d02 = v0[1](0)*p0[1](1) - v0[1](1)*p0[1](0);
	double d03 = v0[2](0)*p0[2](1) - v0[2](1)*p0[2](0);
	Eigen::Vector3d N[3];
	for (int i = 0; i < 3; i++){
		N[i] = (v0[i].cross(p0[i])) / f;
	}


	double A12, B12, A13, B13, C13;
	A12 = V0[0].dot(V0[1]);
	B12 = sqrt(1 - A12*A12);
	A13 = V0[0].dot(V0[2]);
	B13 = (V0[1].dot(V0[2]) - A13*A12) / B12;
	Eigen::Vector3d temp = V0[0].cross(V0[1]);
	C13 = temp.dot(V0[2]) / sqrt(temp.dot(temp));

	V1[0] = { 1, 0, 0 };
	V1[1] = { A12, B12, 0 };
	V1[2] = { A13, B13, C13 };

	Eigen::Matrix<double,3,3> M0, M1, Rm, Rv;
	Mat m_M0(3, 3, CV_64F), m_M1(3, 3, CV_64F), m_Rm(3, 3, CV_64F);
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			m_M0.at<double>(i, j) = V0[j](i);
			m_M1.at<double>(i, j) = V1[j](i);
		}
	}
	/*M0 << V0[0](0), V0[1](0), V0[2](0),
		V0[0](1), V0[1](1), V0[2](1),
		V0[0](2), V0[1](2), V0[2](2);

	M1 << V1[0](0), V1[1](0), V1[2](0),
		V1[0](1), V1[1](1), V1[2](1),
		V1[0](2), V1[1](2), V1[2](2);*/
	Mat I_3 = Mat::eye(3, 3, CV_64F);
	invert(m_M0, m_Rm, CV_SVD);
	m_Rm = m_M1*m_Rm;

	Mat svd_u, svd_vt, svd_w;
	SVD::compute(m_Rm, svd_w, svd_u, svd_vt, SVD::FULL_UV);

	m_Rm = svd_u * svd_vt;
	
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			Rm(i, j) = m_Rm.at<double>(i, j);
		}
	}

	M0 << sqrt(f*f + d01*d01), 0, 0,
		0, f, d01,
		0, -d01, f;

	M1 << v0[0](0), -v0[0](1), 0,
		v0[0](1), v0[0](0), 0,
		0, 0, 1;

	Rv = M1*M0 / M0(0, 0);
//	M1 = Rv*Rv.inverse();
	Eigen::Vector3d v1[3], p1[3], p2[3];
	double a2[3], b2[3], c2[3], d2[3], length;
	for (int i = 0; i < 3; i++){
		v1[i] = Rv.inverse()*v0[i];
		p1[i] = Rv.inverse()*p0[i];
		p2[i] = p1[i] * f / p1[i](2);

		a2[i] = v1[i](0) - p2[i](0)*v1[i](2) / f;
		b2[i] = v1[i](1) - p2[i](1)*v1[i](2) / f;
		length = sqrt(a2[i] * a2[i] + b2[i] * b2[i]);
		a2[i] = a2[i] / length;
		b2[i] = b2[i] / length;

		d2[i] = a2[i] * p2[i](1) - b2[i] * p2[i](0);
		c2[i] = d2[i] / f;
	}
	
	double miu[16], fai[16];
	double a22 = a2[1], b22 = b2[1], c22 = c2[1];
	double a23 = a2[2], b23 = b2[2], c23 = c2[2];
	miu[1] = B12*a22 * C13*b23;
	miu[2] = B12*a22 * B13*b23 - B12*b22*B13*a23;
	miu[3] = B12*b22 * C13*a23;
	miu[4] = A12*c22*B13*a23 - B12*a22*A13*c23;
	miu[5] = -A12*c22*C13*a23;
	miu[6] = -B12*a22*C13*c23;
	miu[7] = B12*c22*B13*a23 - B12*a22*B13*c23;
	miu[8] = -B12*c22*C13*a23;
	miu[9] = A12*b22*B13*a23 - B12*a22*A13*b23;
	miu[10] = -A12*b22*C13*a23;
	miu[11] = B12*c22*C13*b23 - B12*b22*C13*c23;
	miu[12] = B12*B13*(c22*b23 - b22*c23);
	miu[13] = A12*C13*(b22*b23 + c22*c23);
	miu[14] = A12*B13*(b22*b23 + c22*c23) - A13*B12*(c22*c23 + b22*b23);
	miu[15] = A12*A13*(c22*b23 - b22*c23);

	fai[1] = miu[1] * miu[1] + miu[6] * miu[6];
	fai[2] = 2 * (miu[1] * miu[2] + miu[6] * miu[7]);
	fai[3] = miu[2] * miu[2] + 2 * (miu[1] * miu[3] + miu[6] * miu[8]) + miu[7] * miu[7] - miu[11] * miu[11];
	fai[4] = 2 * (miu[2] * miu[3] + miu[7] * miu[8] - miu[11] * miu[12]);
	fai[5] = miu[3] * miu[3] + miu[8] * miu[8] - miu[12] * miu[12];
	fai[6] = 2 * (miu[1] * miu[4] + miu[6] * miu[9]);
	fai[7] = 2 * (miu[1] * miu[5] + miu[2] * miu[4] + miu[6] * miu[10] + miu[7] * miu[9] - miu[11] * miu[13]);
	fai[8] = 2 * (miu[2] * miu[5] + miu[3] * miu[4] + miu[7] * miu[10] + miu[8] * miu[9] - miu[11] * miu[14] - miu[12] * miu[13]);
	fai[9] = 2 * (miu[3] * miu[5] + miu[8] * miu[10] - miu[12] * miu[14]);
	fai[10] = miu[4] * miu[4] + miu[9] * miu[9] - miu[13] * miu[13];
	fai[11] = 2 * (miu[4] * miu[5] + miu[9] * miu[10] - miu[11] * miu[15] - miu[13] * miu[14]);
	fai[12] = miu[5] * miu[5] + miu[10] * miu[10] - miu[14] * miu[14] - 2 * miu[12] * miu[15];
	fai[13] = -2 * miu[13] * miu[15];
	fai[14] = -2 * miu[14] * miu[15];
	fai[15] = -miu[15] * miu[15];

	Eigen::VectorXd delta(9);
	delta[0] = fai[1] + fai[6] + fai[10] + fai[13] + fai[15];
	delta[1] = 2 * (fai[2] + fai[7] + fai[11] + fai[14]);
	delta[2] = 4 * (-fai[1] + fai[3] + fai[8] + fai[12] + fai[15]) + 2 * (-fai[6] + fai[13]);
	delta[3] = 2 * (fai[11] - fai[7]) + 6 * (fai[14] - fai[2]) + 8 * (fai[4] + fai[9]);
	delta[4] = 6 * fai[1] - 8 * fai[3] + 16 * fai[5] - 2 * fai[10] + 8 * fai[12] + 6 * fai[15];
	delta[5] = 6 * fai[2] - 8 * fai[4] - 2 * fai[7] + 8 * fai[9] - 2 * fai[11] + 6 * fai[14];
	delta[6] = -4 * fai[1] + 4 * fai[3] + 2 * fai[6] - 4 * fai[8] + 4 * fai[12] - 2 * fai[13] + 4 * fai[15];
	delta[7] = -fai[2] + 2 * fai[7] - 2 * fai[11] + 2 * fai[14];
	delta[8] = fai[1] - fai[6] + fai[10] - fai[13] + fai[15];

	//生成一个多项式解算器
	Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
	//解多项式
	solver.compute(delta);
	//导出根
	const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType & r = solver.roots();
	
	std::vector<double> real_root;
	solver.realRoots(real_root, 0.0001);

	std::vector<double> cosa, sina, cosb, sinb;
	Eigen::Matrix<double, 3, 3> Rab, Ro2i;
	for (int i = 0; i < real_root.size(); i++){
		double t = real_root[i];
		double sigma1, sigma2, sigma3;
		cosa.push_back((1 - t*t) / (1 + t*t));
		sina.push_back(2 * t / (1 + t*t));
		sigma1 = miu[1] * cosa[i] * cosa[i] + miu[2] * cosa[i] * sina[i] + miu[3] * sina[i] * sina[i] + miu[4] * cosa[i] + miu[5] * sina[i];
		sigma2 = miu[6] * cosa[i] * cosa[i] + miu[7] * cosa[i] * sina[i] + miu[8] * sina[i] * sina[i] + miu[9] * cosa[i] + miu[10] * sina[i];
		sigma3 = miu[11] * cosa[i] * sina[i] + miu[12] * sina[i] * sina[i] + miu[13] * cosa[i] + miu[14] * sina[i] + miu[15];
		cosb.push_back(sigma1 / sigma3);
		sinb.push_back(sigma2 / sigma3);

		Rab << cosb[i], sinb[i]*sina[i], sinb[i]*cosa[i],
			0, cosa[i], -sina[i],
			-sinb[i], cosb[i]*sina[i], cosb[i]*cosa[i];

		/*temp = Rab*V1[0];
		double error = temp(0)*b2[0] - temp(1)*a2[0] + temp(2)*c2[0];
		double sigma[7] = { 0 };
		sigma[1] = B12*c22*sina[i] + A12*b22;
		sigma[2] = B12*b22*sina[i] - A12*c22;
		sigma[3] = B12*a22*cosa[i];
		sigma[4] = C13*c23*cosa[i] + B13*c23*sina[i] + A13*b23;
		sigma[5] = C13*b23*cosa[i] + B13*b23*sina[i] - A13*c23;
		sigma[6] = B13*a23*cosa[i] - C13*a23*sina[i];
		error = (sigma[3] * sigma[5] - sigma[2] * sigma[6]) / (sigma[1] * sigma[5] - sigma[2] * sigma[4]) - cosb[i];
		error = cosb[i] * (B12*c22*sina[i] + A12*b22) + sinb[i] * (B12*b22*sina[i] - A12*c22) - B12*a22*cosa[i];*/

		/*temp = Rab*V1[1];
		error = temp(0)*b2[1] - temp(1)*a2[1] + temp(2)*c2[1];
		temp = Rab*V1[2];
		error = temp(0)*b2[2] - temp(1)*a2[2] + temp(2)*c2[2];*/

		Ro2i = Rv*Rab*Rm;

		//求T
		Eigen::Matrix<double, 3, 3> temp_matrix;
		Eigen::Vector3d P3[3], temp_T;
		for (int j = 0; j < 3; j++){
			P3[j] = Ro2i*P0[j];
		}
		temp << -P3[0].dot(N[0]), -P3[1].dot(N[1]), -P3[2].dot(N[2]);
		temp_matrix << N[0](0), N[0](1), N[0](2),
			N[1](0), N[1](1), N[1](2),
			N[2](0), N[2](1), N[2](2);
		temp_T = temp_matrix.inverse()*temp;

		//导出
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				R.at<double>(j+3*i, k) = Ro2i(j, k);
			}
			tvec.at<double>(i, j) = temp_T(j);
		}
	}

	return 0;
}

int _tmain(int argc, _TCHAR* argv[])
{
	Mat R(24, 3, CV_64F);
	Mat T(8, 3, CV_64F);
	double f = 4.037229;
	Mat ilines(3, 6, CV_64F, Scalar(0));
	Mat olines(3, 6, CV_64F, Scalar(0));

	olines.at<double>(0, 0) = 1;
//	olines.at<double>(0, 3) = 1;
	olines.at<double>(1, 1) = 1;
	double* data = olines.ptr<double>(2);
	data[0] = -1 / sqrt(2);
	data[1] = 1 / sqrt(2);
	data[4] = 81;

	data = ilines.ptr<double>(0);
	data[0] = -36 / sqrt(68377);
	data[1] = 259 / sqrt(68377);
	data[3] = 166 * 3.69 / 1000;
	data[4] = -134 * 3.69 / 1000;
	data[5] = f;
	data = ilines.ptr<double>(1);
	data[0] = 257 / sqrt(67493);
	data[1] = 38 / sqrt(67493);
	data[3] = 166 * 3.69 / 1000;
	data[4] = -134 * 3.69 / 1000;
	data[5] = f;
	data = ilines.ptr<double>(2);
	data[0] = 221 / sqrt(137050);
	data[1] = 297 / sqrt(137050);
	data[3] = -91 * 3.69 / 1000;
	data[4] = -172 * 3.69 / 1000;
	data[5] = f;

	threeLines(R, T, f, olines, ilines); 
	for (int i = 0; i < 24; i++){
		for (int j = 0; j < 3; j++){
			cout << R.at<double>(i, j)<<" ";
		}
		cout << endl;
	}
	system("pause");
	return 0;
}

