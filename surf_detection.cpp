#include "stdafx.h"
#include "common.h"
#include "surf_structures.h"
#include "surf_detection.h"

namespace surf
{

	IntegralImage computeIntegralImage(Mat_<uchar> img)
	{
		int height = img.rows + 1;
		int width = img.cols + 1;

		IntegralImage J(height, width);
		J.setTo(0);

		for (int i = 1; i < height; i++)
		{
			for (int j = 1; j < width; j++)
			{
				J(i, j) = J(i, j - 1) + J(i - 1, j) + img(i - 1, j - 1) - J(i - 1, j - 1);
			}
		}

		return J;
	}

	int computeSubregionSum
	(
		IntegralImage J,
		Subregion subRegion
	)
	{
		int x = subRegion.x;
		int y = subRegion.y;
		int height = subRegion.height;
		int width = subRegion.width;

		int sum = 0;

		sum = J(x + height, y + width) + J(x, y)
			- J(x, y + width) - J(x + height, y);

		return sum;
	}

	int computeDxxSum
	(
		IntegralImage J, 
		int filterSize, 
		int x, 
		int y
	) 
	{
		int L0 = filterSize / 3;
		int middle = filterSize / 2;

		int r1 = (filterSize - (2 * L0 - 1)) / 2;

		int height = 2 * L0 - 1;
		int width = L0;

		int white1_x = x - middle + r1;
		int white1_y = y - middle + 0;

		int black_x = x - middle + r1;
		int black_y = y - middle + L0;

		int white2_x = x - middle + r1;
		int white2_y = y - middle + 2 * L0;

		Subregion white1{ white1_x, white1_y, height, width };
		Subregion black{ black_x, black_y, height, width };
		Subregion white2{ white2_x, white2_y, height, width };

		int dxx_sum = computeSubregionSum(J, white1)
			- 2 * computeSubregionSum(J, black)
			+ computeSubregionSum(J, white2);

		return dxx_sum;
	}

	int computeDyySum
	(
		IntegralImage J, 
		int filterSize, 
		int x, 
		int y
	) 
	{
		int L0 = filterSize / 3;
		int middle = filterSize / 2;

		int c1 = (filterSize - (2 * L0 - 1)) / 2;

		int height = L0;
		int width = 2 * L0 - 1;

		int white1_x = x - middle + 0;
		int white1_y = y - middle + c1;

		int black_x = x - middle + L0;
		int black_y = y - middle + c1;

		int white2_x = x - middle + 2 * L0;
		int white2_y = y - middle + c1;

		Subregion white1{ white1_x, white1_y, height, width };
		Subregion black{ black_x, black_y, height, width };
		Subregion white2{ white2_x, white2_y, height, width };

		int dyy_sum = computeSubregionSum(J, white1)
			- 2 * computeSubregionSum(J, black)
			+ computeSubregionSum(J, white2);

		return dyy_sum;
	}

	int computeDxySum
	(
		IntegralImage J, 
		int filterSize, 
		int x, 
		int y
	)
	{
		int L0 = filterSize / 3;
		int middle = filterSize / 2;

		int r1 = (filterSize - (2 * L0 + 1)) / 2;
		int r2 = r1 + L0 + 1;

		int height = L0;
		int width = L0;

		int white1_x = x - middle + r1;
		int white1_y = y - middle + r1;

		int black1_x = x - middle + r2;
		int black1_y = y - middle + r1;

		int white2_x = x - middle + r2;
		int white2_y = y - middle + r2;

		int black2_x = x - middle + r1;
		int black2_y = y - middle + r2;

		Subregion white1{ white1_x, white1_y, height, width };
		Subregion black1{ black1_x, black1_y, height, width };
		Subregion white2{ white2_x, white2_y, height, width };
		Subregion black2{ black2_x, black2_y, height, width };

		int dxy_sum = computeSubregionSum(J, white1)
			- computeSubregionSum(J, black1)
			+ computeSubregionSum(J, white2)
			- computeSubregionSum(J, black2);

		return dxy_sum;
	}

	BlobResponse computeBlobResponseMap
	(
		IntegralImage J,
		int filterSize,
		int metricThreshold
	) 
	{
		float w = 0.9f;
		int middle = filterSize / 2;
		int L0 = filterSize / 3;

		int height = J.rows - 1;
		int width = J.cols - 1;

		BlobResponse response;

		Mat_<double> blobResponseMap(height, width);
		blobResponseMap.setTo(0);

		Mat_<int> traceSign(height, width);
		blobResponseMap.setTo(0);

		Mat_<double> dxx(height, width);
		dxx.setTo(0);

		Mat_<double> dyy(height, width);
		dyy.setTo(0);

		Mat_<double> dxy(height, width);
		dxy.setTo(0);

		for (int i = middle; i < height - middle; i++)
		{
			for (int j = middle; j < width - middle; j++)
			{
				int dxx_sum = computeDxxSum(J, filterSize, i, j);
				int dyy_sum = computeDyySum(J, filterSize, i, j);
				int dxy_sum = computeDxySum(J, filterSize, i, j);

				double det_response = (double) dxx_sum * dyy_sum - (w * (double)dxy_sum * dxy_sum);
				det_response /= ((double)L0 * L0 * L0 * L0);

				if (det_response >= metricThreshold) 
				{
					blobResponseMap(i, j) = det_response;
				}

				double trace = (double) dxx_sum + dyy_sum;
				if (trace > 0)
				{
					traceSign(i, j) = 1;
				}
				else
				{
					traceSign(i, j) = -1;
				}
			}
		}

		response.blobResponseMap = blobResponseMap;
		response.trace = traceSign;

		return response;
	}


	SURFOctaves generateOctaves
	(
		IntegralImage J,
		int octaveNumber
	) 
	{
		SURFOctaves octaves;

		int height = J.rows - 1;
		int width = J.cols - 1;

		int filter_size = FILTER_SIZE;
		int next_filter = NEXT_FILTER;

		for (int o = 0; o < octaveNumber; o++)
		{
			for (int i = 0; i < NUM_SCALE_LEVELS; i++)
			{
				int current_filter = filter_size + i * next_filter;
				BlobResponse response = computeBlobResponseMap(J, current_filter);
				octaves.push_back(response);
			}

			filter_size += next_filter;
			next_filter *= 2;
		}

		return octaves;
	}

	BlobResponse getBlobResponseFromOctaves
	(
		int octave_number,
		int scale_number,
		SURFOctaves octaves
	)
	{
		int i = octave_number * NUM_SCALE_LEVELS + scale_number;
		BlobResponse response = octaves[i];
		return response;
	}

	std::vector<SURFPoint> nonMaximumSuppression
	(
		SURFOctaves octaves
	)
	{
		std::vector<SURFPoint> keypoints;

		int height = octaves[0].blobResponseMap.rows;
		int width = octaves[0].blobResponseMap.cols;

		Mat_<uchar> non_maxima(height, width);
		non_maxima.setTo(0);

		int octave_number = octaves.size() / NUM_OCTAVES;
		
		// neighborhood size
		int n = NEIGHBORHOOD_SIZE;
		int middle = n / 2;
		int valid_scales = NUM_SCALE_LEVELS - 1;

		for (int o = 0; o < octave_number; o++)
		{
			for (int s = 1; s < valid_scales; s++)
			{
				BlobResponse previous = getBlobResponseFromOctaves(o, s - 1, octaves);
				BlobResponse current = getBlobResponseFromOctaves(o, s, octaves);
				BlobResponse next = getBlobResponseFromOctaves(o, s + 1, octaves);

				for (int i = middle; i < height - middle; i++)
				{
					for (int j = middle; j < width - middle; j++)
					{
						if (current.blobResponseMap(i, j) > 0)
						{
							int crop_i = i - middle;
							int crop_j = j - middle;

							bool isMax = true;

							for (int u = 0; u < n; u++)
							{
								for (int v = 0; v < n; v++)
								{
									if (current.blobResponseMap(crop_i + u, crop_j + v) > current.blobResponseMap(i, j))
									{
										isMax = false;
										break;
									}

									if (previous.blobResponseMap(crop_i + u, crop_j + v) > current.blobResponseMap(i, j))
									{
										isMax = false;
										break;
									}

									if (next.blobResponseMap(crop_i + u, crop_j + v) > current.blobResponseMap(i, j))
									{
										isMax = false;
										break;
									}
								}

								if (!isMax)
								{
									break;
								}
							}

							if (isMax)
							{
								non_maxima(i, j) = 255;

								SURFPoint keypoint;
								keypoint.x = i;
								keypoint.y = j;
								keypoint.scale = ((1 << o) * s + 1) * 0.4;
								keypoint.metric = current.blobResponseMap(i, j);
								keypoint.signOfLaplacian = current.trace(i, j);
								keypoint.orientation = 0.0f;

								keypoints.push_back(keypoint);
							}
						}
					}
				}
			}
		}

		imshow("keypoints", non_maxima);
		waitKey(0);

		std::sort(keypoints.begin(), keypoints.end());
		return keypoints;
	}

	std::vector<SURFPoint> eliminateMarginKeypoints
	(
		std::vector<SURFPoint> keypoints,
		int height,
		int width
	)
	{
		std::vector<SURFPoint> valid_keypoints;

		for (int i = 0; i < keypoints.size(); i++)
		{
			SURFPoint current_point = keypoints[i];
			float scale = current_point.scale;
			int x = current_point.x;
			int y = current_point.y;

			int middle = scale;

			int corner_x = x - 10 * scale - middle;
			int corner_y = y - 10 * scale - middle;
			if ((corner_x < 1 || corner_x > height) || (corner_y < 1 || corner_y > width))
				continue;

			corner_x = x - 10 * scale - middle;
			corner_y = y + 10 * scale + middle;
			if ((corner_x < 1 || corner_x > height) || (corner_y < 1 || corner_y > width))
				continue;

			corner_x = x + 10 * scale + middle;
			corner_y = y - 10 * scale - middle;
			if ((corner_x < 1 || corner_x > height) || (corner_y < 1 || corner_y > width))
				continue;

			corner_x = x + 10 * scale + middle;
			corner_y = y + 10 * scale + middle;
			if ((corner_x < 1 || corner_x > height) || (corner_y < 1 || corner_y > width))
				continue;

			valid_keypoints.push_back(current_point);
		}

		return valid_keypoints;
	}


}