#pragma once

struct DeepClassifierMatch {
	DeepClassifierMatch() :
		queryIdx(0),
		trainIdx(0),
		distance(0.0) {

	}

	DeepClassifierMatch(int _queryIdx, int _trainIdx, int _distance) :
		queryIdx(_queryIdx),
		trainIdx(_trainIdx),
		distance(_distance) {

	}

	int queryIdx;
	int trainIdx;
	float distance;
};
