# pose-estimation-with-three-lines
A C/C++ implementation of Determination of the Attitude of 3-D Objects from a Single Perspective View<br>
this project include  the function :<br>
int threeLines(cv::Mat& R, cv::Mat& tvec, const double& f, const cv::Mat& olines, const cv::Mat& ilines);<br>
input: olines:Coordinates of the three lines in the object coordinate system;<br>
       ilines:Coordinates of the three lines corresponding to the "olines" respectively in the image coordinate system;<br>
       f:focus of the camera;<br>
output:R:Rotation matrix;<br>
       tvec:Translation vector.
