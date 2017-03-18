#include"opencv2/opencv.hpp"

using namespace cv;

/*用三条直线计算位置姿态。olines是物体坐标系下直线的表示，为3*6的矩阵，每行的6个参数表示一条直线。前三个为直线的方向向量，后三个表示直线上一个点的三维坐标。
ilines是图像中对应直线在摄像机坐标系下的表示，其中各参数的意义与olines中相同。
Mat中数据需要是double型，为24*3矩阵，tvec为8*3，double型*/
int threeLines(cv::Mat& R, cv::Mat& tvec, const double& f, const cv::Mat& olines, const cv::Mat& ilines);