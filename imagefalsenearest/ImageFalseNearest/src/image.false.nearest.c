/*
 * Copyright 2016 Martha Dais Ferreira, Rodrigo Fernandes de Mello
 * 
 * This file is part of ImageFalseNearest.
 *
 * ImageFalseNearest is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * ImageFalseNearest is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * should have received a copy of the GNU General Public License
 * along with ImageFalseNearest.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * See the file "COPYING" for the text of the license.
 * 
 * Contact: 
 * 	Martha D. Ferreira: daismf@icmc.usp.br
 * 	Rodrigo Fernandes de Mello: mello@icmc.usp.br
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>

/* 
 * Functions to sort the closest neighbors
 */
struct IndexedInteger {
  double value;
  int index;
} IndexedInteger;

/* 
 * Add indexes information in the input array
 */
void addIndices(struct IndexedInteger *array, size_t num) {
  int i;
  for (i = 0; i < num; ++i) {
    array[i].index = i;
  }
}

/* 
 * Sort the input array based on the values
 */
int comp(const void * a, const void * b) {
  struct IndexedInteger *i1, *i2;
  i1 = (struct IndexedInteger*) a;
  i2 = (struct IndexedInteger*) b;
  if (i1->value > i2->value) return  1;
  if (i1->value < i2->value) return -1;
  return 0;
}

/* 
 * Convert the input vector in a matrix
 */
double** toMatrix(double *euclidean, int nrows, int ncols) {
	int i, j;
	double **matrix = (double **) malloc(nrows * sizeof(double *));

	for (i = 0; i < nrows; i++) {
		matrix[i] = (double *) malloc(ncols * sizeof(double));
	}

	for (i=0; i < nrows; i++) {
		for (j=0; j < ncols; j++) {
			matrix[i][j] = euclidean[i+(j*(nrows-1))];
		}
	}
	return matrix;
}

/* 
 * Free matrix memory
 */
void freeMatrix(double **matrix, int nrows) {
	int i;
	for (i = 0; i < nrows; i++) free(matrix[i]);
	free(matrix);
}

/* 
 * Compute Euclidian distance between two vectors
 */
double squared_difference(double *vector1, double *vector2, int ncols) {
	int i;
	double sum = 0.0;

	for (i = 0; i < ncols; i++) {
		sum += pow(vector1[i] - vector2[i], 2.0);
	}
	return sqrt(sum);
}

/* 
 * Compute the modified FNN considering a phase space obtained from an image
 */
void ImageFalseNearest(double *in_vectors, 		// phase space with nvectors points, each vector with ncols elements
			int *in_nvectors, 		// number of vectors in the phase space
			int *in_ncols, 			// length of each vector in phase space
			double *in_euclidean, 		// euclidean distances between vectors
			int *in_endT, 			// time lag
			int *in_nn, 			// number of nearby neighbors to be considered
			double *in_rt, 			// radius for the neighborhood
			double *out, 			// fraction of false nearests
			int *out2) {                    // total number of nearests

	int i, j, nn, id, endT, nvectors, nrows, ncols, num, denum;
	double **vectors, **euclidean; 
	double dst;
	int *ids;
	double rt;

	// Receiving input variables
	endT = *in_endT;
	rt = *in_rt;
	nn=*in_nn;
	nvectors=*in_nvectors;
	ncols = *in_ncols;
	vectors=toMatrix(in_vectors, nvectors, ncols);
	euclidean=toMatrix(in_euclidean, nvectors, nvectors);

	// Algorithm
	nrows = nvectors - endT;
	num = denum = 0;

	// To keep neighbors ids
	ids = (int*) malloc(nn * sizeof(int));

	// gets all candidate neighbors for a vector in phase space
	struct IndexedInteger *possiblesNeighbors = (struct IndexedInteger *) malloc(nvectors * sizeof(struct IndexedInteger));
	
	// for each point or vector in phase space
	for (i = 0; i < nrows; i++) { 

		// Find nearest vectors according to variable 'nn' (nn=1)
		memcpy(possiblesNeighbors, euclidean[i], sizeof(double) * nvectors);
		addIndices(possiblesNeighbors, nvectors);		
		qsort(possiblesNeighbors, nrows, sizeof(struct IndexedInteger),comp);
		
		// Retrieving the nearby neighbors for the current vector under analysis
		id = 0;
		j = 0;
		while (id < nn) {
			if (i != possiblesNeighbors[j].index) {
				ids[id++] = possiblesNeighbors[j].index;
			}
			j++;
		}
		
		// for each neighbour
		for(j = 0; j < nn; j++) {
			dst = squared_difference(vectors[i], vectors[ids[j]], ncols);
			dst = (dst + squared_difference(vectors[i+endT], vectors[ids[j]+endT], ncols)) / dst;

			// if ratio > rt, we have another false neighbour
			if (dst>rt) num++; 
		} 

		// add number of neighbours founds to total number of neighbours
		denum += nn; 
	} 

	// free memories
	free(possiblesNeighbors);
	free(ids);
	freeMatrix(vectors, nvectors);
	freeMatrix(euclidean, nvectors);

	// this is the fraction of false neighbours
	*out = denum != 0 ? num/(denum*1.0) : -1; 

	// this is the total number of neighbours
	*out2 = denum; 
}
