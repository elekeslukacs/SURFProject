#pragma once

#include "surf_structures.h"

namespace surf
{
	// Function that computes integral image
	IntegralImage computeIntegralImage(Mat_<uchar> img);

	int computeSubregionSum
	(
		IntegralImage J,
		Subregion subRegion
	);

	BlobResponse computeBlobResponseMap
	(
		IntegralImage J,
		int filterSize,
		int metricThreshold = 1000
	);


	SURFOctaves generateOctaves
	(
		IntegralImage J,
		int octaveNumber = NUM_OCTAVES
	);


	std::vector<SURFPoint> nonMaximumSuppression
	(
		SURFOctaves octaves
	);


	std::vector<SURFPoint> eliminateMarginKeypoints
	(
		std::vector<SURFPoint> keypoints,
		int height,
		int width
	);


}
