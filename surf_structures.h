#pragma once
#include <vector>

namespace surf
{
	struct SURFPoint
	{
		// location
		int x;
		int y;

		// scale at which the point was detected
		float scale;

		// strength of detected feature
		float metric;

		// orientation in radians
		float orientation;

		// Sign of Laplacian, -1 or 1
		int signOfLaplacian;

		bool operator < (const SURFPoint& o) const {
			return metric > o.metric;
		}
	};

	struct BlobResponse
	{
		Mat_<double> blobResponseMap;
		Mat_<int> trace;
	};

	struct Subregion
	{
		int x;
		int y;

		int height;
		int width;
	};

	typedef std::vector<float> SURFDescriptor;

	typedef std::vector<BlobResponse> SURFOctaves;

	typedef Mat_<int> IntegralImage;

	// Constants

	static const int FILTER_SIZE = 9;
	static const int NEXT_FILTER = 6;
	static const int NUM_SCALE_LEVELS = 4;
	static const int NEIGHBORHOOD_SIZE = 3;
	static const int NUM_OCTAVES = 4;
}