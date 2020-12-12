#include "stdafx.h"
#include "common.h"
#include "surf_structures.h"
#include "surf_detection.h"

using namespace std;

Mat selectImage(int flag = -1) {
	char fname[MAX_PATH];
	Mat src;
	openFileDlg(fname);
	src = imread(fname, flag);
	return src;
}

int main() 
{

	cout << "Hello World!\n";

	//Mat_<uchar> img(3, 3);
	//img.setTo(1);
	//img(0, 0) = 2;
	//Mat_<int> J = surf::computeIntegralImage(img);

	//surf::Subregion region{ 0, 1, 2, 2 };
	//int sum = surf::computeSubregionSum(J, region);

	//cout << "Sum is " << sum << endl;

	//for (int i = 0; i < 4; i++)
	//{
	//	for (int j = 0; j < 4; j++)
	//	{
	//		cout << J(i, j) << "  ";
	//	}

	//	cout << endl;
	//}

	Mat_<uchar> img = selectImage(0);

	Mat_<int> J = surf::computeIntegralImage(img);

	//surf::BlobResponse response = surf::computeBlobResponseMap(J, 9);

	//imshow("img", img);
	//imshow("response", response.blobResponseMap);
	//waitKey(0);

	surf::SURFOctaves octaves = surf::generateOctaves(J);
	vector<surf::SURFPoint> keypoints =  surf::nonMaximumSuppression(octaves);

	cout << keypoints.size();

	return 0;
}