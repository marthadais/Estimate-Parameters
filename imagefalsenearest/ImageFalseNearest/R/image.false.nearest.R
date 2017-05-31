# Copyright 2016 Martha Dais Ferreira, Rodrigo Fernandes de Mello
# 
# This file is part of ImageFalseNearest.
# 
# ImageFalseNearest is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# ImageFalseNearest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ImageFalseNearest.  If not, see <http://www.gnu.org/licenses/>.
# 
# See the file "COPYING" for the text of the license.
# 
# Contact: 
# 	Martha D. Ferreira: daismf@icmc.usp.br
# 	Rodrigo Fernandes de Mello: mello@icmc.usp.br


#  Compute the image embedding according to the mask dimension
embedd.image <- function(image.matrix, row, col) {
	vectors = NULL
	for (i in 1:(nrow(image.matrix)-row+1)) {
		for (j in 1:(ncol(image.matrix)-col+1)) {
			vectors = rbind(vectors, as.vector(image.matrix[i:(i+row-1), j:(j+col-1)]))
		}
	}
	vectors
}

# Compute the modified FNN for all possible codomain obtained from linear transformation
assess.number.of.kernels <- function(image.matrix, mask.nrows=10, mask.ncols=10, end.T=1, nkernels=1:5, rt=10, nn=1, plot.opt=T) {
	counter = 1
	results = list()
	heatmap = matrix(0, nrow=length(nkernels), ncol=end.T)

    # compute the embedding according to the input mask dimension
	vectors = embedd.image(image.matrix, mask.nrows, mask.ncols)
	vectors.nrow = nrow(vectors)
	vectors.ncol = ncol(vectors)

	# For all codomain dimensions possible
	for (k in 1:length(nkernels)) {
		nkernel = nkernels[k]

		# Compute a random linear transformation
		A = matrix(rnorm(mean=0, sd=1, n=vectors.ncol*nkernel), nrow=nkernel)
		A = (A - mean(A)) / sd(A)
		another.space = t(A%*%t(vectors))

		# Compute the euclidian distance between vectors in a codomain
		dist = as.matrix(dist(another.space))

		# Convert all results in a vector to call C funtion
		another.space <- (another.space - matrix(rep(apply(another.space, 2, min),nrow(another.space)),ncol=ncol(another.space),byrow=T)) / matrix(rep(apply(apply(another.space, 2, range), 2, diff), nrow(another.space)), ncol=ncol(another.space), byrow=T)
		another.space[which(is.nan(another.space))] = 0

		res <- rep(0, end.T)
		res2 <- rep(0, end.T)

		# For all possible number of neighbors
		for(i in 1:end.T) {
			# Call C function to compute the modified FNN
			a <- .C("ImageFalseNearest", vectors=as.double(another.space), 
				nvectors=as.integer(nrow(another.space)), 
				ncols=as.integer(ncol(another.space)),
				euclidean=as.double(dist),
				endT=as.integer(i), nn=as.integer(nn), rt=as.double(rt), 
				out=as.double(res[i]), out2=as.integer(res2[i]))
			res[i] = a[["out"]]
			res2[i]= a[["out2"]]
		}
		# Organize the returned information
		res <- rbind(res, res2)
		res[res==(-1)] <- NA
		rownames(res) <- c("fraction", "total")
		colnames(res) <- paste("m", 1:end.T, sep="")

		results[[k]] = res
		heatmap[counter,] = res[1,]
		counter = counter + 1
	}

	ret = list()
	ret$results = results
	ret$heatmap = heatmap

	ret
}

# Compute the modified FNN for all possible phase spaces obtained from an image
image.false.nearest <- function(image.matrix, mask.nrows=10, mask.ncols=10, end.T=1, rt=10, nn=1, plot.opt=T) {
	results = list()

	# For all mask row dimensions possible
	for (row in 1:mask.nrows) {
		results[[row]] = list()
		
		# For all mask colunm dimensions possible
		for (col in 1:mask.ncols) {
			# compute the embedding
			vectors = embedd.image(image.matrix, row, col)
			vectors.nrow = nrow(vectors)
			vectors.ncol = ncol(vectors)

			# Compute the euclidian distance between vectors in a phase space
			dist = as.matrix(dist(vectors))

			# Convert all results in a vector to call C funtion
			vectors <- (vectors - matrix(rep(apply(vectors, 2, min),nrow(vectors)),ncol=vectors.ncol,byrow=T)) / matrix(rep(apply(apply(vectors, 2, range), 2, diff), nrow(vectors)), ncol=vectors.ncol, byrow=T)
			vectors[which(is.nan(vectors))] = 0

			res <- rep(0, end.T)
			res2 <- rep(0, end.T)

			# For all possible number of neighbors
			for(i in 1:end.T) {
			    # Call C function to compute the modified FNN
				a <- .C("ImageFalseNearest", vectors=as.double(vectors), 
					nvectors=as.integer(nrow(vectors)), 
					ncols=as.integer(vectors.ncol),
					euclidean=as.double(dist),
					endT=as.integer(i), 
					nn=as.integer(nn), rt=as.double(rt), 
					out=as.double(res[i]), out2=as.integer(res2[i]))
				res[i] = a[["out"]]
				res2[i] = a[["out2"]]
			}
			# Organize the returned information
			res <- rbind(res, res2)
			res[res==(-1)] <- NA
			rownames(res) <- c("fraction", "total")
			colnames(res) <- paste("m", 1:end.T, sep="")

			results[[row]][[col]] = res
		}
	}
	
	# Plot results
	if (plot.opt==T) {
		plot.false.nearest(results)
	}

	results
}

# Plot results for mask evaluation
plot.false.nearest <- function(results, type="max") {

	heatmap = matrix(0, nrow=length(results), ncol=length(results[[1]]))

	for (rowId in 1:length(results)) {
		for (colId in 1:length(results[[rowId]])) {
			value = 0
			if (type == "max") {
				value = max(results[[rowId]][[colId]][1,])
			} else if (type == "min") {
				value = min(results[[rowId]][[colId]][1,])
			} else if (type == "ave") {
				value = mean(results[[rowId]][[colId]][1,])
			}
			heatmap[rowId, colId] = value
		}
	}
	image(x=1:nrow(heatmap), y=1:ncol(heatmap), z=heatmap, 
	      xlim=c(1,nrow(heatmap)), ylim=c(1,ncol(heatmap)), 
	      	xlab = "number of rows", ylab ="number of cols")
}

# Plot results for number of convolutional units evaluation
plot.assess.number.of.kernels <- function(ret, type="max") {

	image(x=1:nrow(ret$heatmap), y=1:ncol(ret$heatmap), z=ret$heatmap, 
	      xlim=c(1,nrow(ret$heatmap)), ylim=c(1,ncol(ret$heatmap)), 
	      	xlab = "number of kernels", ylab ="fraction of false nearest neighbors")
}
