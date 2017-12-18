# Copyright 2017 Martha Dais Ferreira, Rodrigo Fernandes de Mello
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

require(png)
require(jpeg)
require(pixmap)
require(ImageFalseNearest)

#
# Choosing best mask size considering all color channels
#  Input:
#     folder_initial: folder with the FNN results
#                  k: number of best mask sizes to be returned
#
#  Output:
#	   ret$histogram: histogram computed over the FNN result
#	    ret$maskList: mask sizes sorted
# 	   ret$maskKList: k best mask sizes
#
analise_mask_MR <- function(folder_initial, k=5){
	
	#Reading FNN result
	normMR = as.matrix(read.table(paste(folder_initial,"meanResultsAllColors.dat",sep="")))
	
	maxRow = nrow(normMR)
	maxCol = ncol(normMR)
	
	#computing histogram
	histogramRow = matrix(0,maxRow,maxCol)
	for(rowId in 1:nrow(normMR)){
		bestM = which.min(normMR[rowId,])
		histogramRow[rowId,bestM] = histogramRow[rowId,bestM] + 1
	}
			
	histogramCol = matrix(0,maxRow,maxCol)
	for(colId in 1:ncol(normMR)){
		bestM = which.min(normMR[,colId])
		histogramCol[bestM,colId] = histogramCol[bestM,colId] + 1
	}
	
	histogram = histogramRow + histogramCol
	
	#Saving histogram
	nameFile = paste(folder_initial,"histogramNormMR.dat",sep="")
	write.table(histogram,nameFile,row.names=F,col.names=F)
	
	
	#Converting histogram matrix to a list
	count = 1
	maskList = c()
	for(i in 1:maxCol){
		row = count
		col = i
		while(col < (maxCol+1)){
			if((histogram[row,col] > 0) && ((row != 1) && (col != 1)))maskList = rbind(maskList,cbind(row,col,histogram[row,col]))
			row = row+1
			col = col+1
		}
		if(row != col){
			row = i
			col = count
			while(row < (maxRow+1)){
				if(histogram[row,col] > 0)maskList = rbind(maskList,cbind(row,col,histogram[row,col]))
				row = row+1
				col = col+1
			}
		}
	}
	
	#Sorting vector
	ids = sort.list(maskList[,3],dec=T)	
	
	#output
	ret = list()
	ret$histogram = histogram
	ret$maskList = maskList[ids,]
	ret$maskKList = maskList[ids[1:k],]
	
	ret
}

#
# FNN to Mask Evaluation
# Input:
#      dataset: the name of the dataset
#     listfile: a file containing paths to dataset images (PNG)
#       maxRow: maximum number to evaluate the row mask size
#       maxCol: maximum number to evaluate the column mask size
#            k: number of the best mask sizes to be returned
#    threshold: threshold to stop the images analysis
#         endT: number of neighbors to be considered in FNN
#
# Output:
#    ret
#            ret$datset: dataset name
#    ret$folder_initial: folder create to output files
#         ret$bestMasks: best mask sizes found
#
mask_evaluation <- function(dataset, listfile, maxRow=0, maxCol=0, k=5, threshold=2e-2, endT=1){
	
	cat("-------Compute FNN for mask estimation.-------\n")
	
	#Creating folders to output
	folder_initial = paste(dataset,"/fnnMaskEstimate/",sep="")
	system(paste("mkdir ",dataset,sep=""))
	system(paste("mkdir ",folder_initial,sep=""))

	#Reading files
	Files = read.table(listfile)
	sizeFiles = length(Files[[1]])
	cat("Analysing ",sizeFiles," files\n")
	
	#Taking images dimension
	channels=1
	# img = readPNG(as.character(Files[[1]][1])) #png
	# img = as.matrix(readJPEG(as.character(Files[[1]][1]))) #jpg
	img = read.pnm(file = as.character(Files[[1]][1]),cellres=1)@grey #pgm
	if(length(dim(img))>2){
		channels=3
		img = img[,,1]
	}
	if((maxRow==0) || (maxCol==0)){
		maxRow = min(floor(nrow(img)/2),floor(ncol(img)/2))
		maxCol = maxRow
		cat("Maximum Mask Selected: ",maxRow,"x",maxCol,"\n")
	}
	
	#Preparing variables
	mean_normMR = matrix(0,maxRow,maxCol)
	mean_normFH = matrix(0,maxRow,maxCol)
	
	#Computing for three color channels
	for(color in 1:channels){
		
		#1-red 2-green 3-blue
		if(color == 1) color_name = "RED"
		if(color == 2) color_name = "GREEN"
		if(color == 3) color_name = "BLUE"
		if(channels == 1) color_name = "GRAY"
		cat("\nAnalysing color ",color_name,"\n")
		
		#Creating folder for color
		folder = paste(folder_initial,color_name,"/",sep="")
		system(paste("mkdir ",folder,sep=""))

		#Randomize examples
		ids = sample(length(Files[[1]]))
		files = Files[[1]][ids]
		
		#Preparing variables
		sumResults = matrix(0, nrow=maxRow, ncol=maxCol)
		histogramRowAll = matrix(0,maxRow,maxCol)
		histogramColAll = matrix(0,maxRow,maxCol)
		histogramFull = matrix(0,maxRow,maxCol)
		histogramFull_before = histogramFull

		frobNorm = c()
		curr_norm = 1
		i=1
		
		#Images analisys
		while ((curr_norm > threshold) && (i < sizeFiles)){
			if(channels==3){
				# img = readPNG(as.character(files[i]))[,,color] #png
				# img = as.matrix(readJPEG(as.character(files[i]))) #jpg
				img = read.pnm(file = as.character(files[i]),cellres=1)@grey #pgm
			} else {
				# img = readPNG(as.character(files[i])) #png
				# img = as.matrix(readJPEG(as.character(files[i]))) #jpg
				img = read.pnm(file = as.character(files[i]),cellres=1)@grey #pgm
			}
			
			#Taking file name for output
			fileFirstname = unlist(strsplit(as.character(files[i]), "/",fixed = TRUE))
			nameId = length(fileFirstname)
			fileFirstname = unlist(strsplit(fileFirstname[nameId], ".",fixed = TRUE))[1]
			cat("file[",i,"]=",fileFirstname,"\n")
			
			#Appling FNN for mask evaluation in a single image considering the current color
			cat("\tCalculating mask size...\n")
			resultsFNN = image.false.nearest(img,mask.nrows=maxRow, mask.ncols=maxCol, end.T=endT, plot.opt=F)
			
			#Creating FNN matrix output
			fractionMap = matrix(0, nrow=length(resultsFNN), ncol=length(resultsFNN[[1]]))
			for (rowId in 1:length(resultsFNN)) {
				for (colId in 1:length(resultsFNN[[rowId]])) {
					fractionMap[rowId, colId] = mean(resultsFNN[[rowId]][[colId]][1,])
				}
			}
			
			#Saving results for each image considering the current color
			#nameFile = paste(folder,fileFirstname,".dat",sep="")
			#write.table(fractionMap,nameFile,row.names=F,col.names=F)
			
			#Summing results
			sumResults = sumResults + fractionMap
			
			#Computing histogram for each image considering the current color
			histogramRow = matrix(0,maxRow,maxCol)
			for(rowId in 1:nrow(fractionMap)){
				bestM = which.min(fractionMap[rowId,])
				histogramRow[rowId,bestM] = histogramRow[rowId,bestM] + 1
			}
			
			histogramCol = matrix(0,maxRow,maxCol)
			for(colId in 1:ncol(fractionMap)){
				bestM = which.min(fractionMap[,colId])
				histogramCol[bestM,colId] = histogramCol[bestM,colId] + 1
			}
			
			histogramRowAll = histogramRowAll + histogramRow
			histogramColAll = histogramColAll + histogramCol
			histogramFull = histogramRowAll + histogramColAll
			
			#Saving results
			meanResults = sumResults/i
			write.table(meanResults,paste(folder,"outputMeanResults.dat",sep=""),row.names=F,col.names=F)
			write.table(histogramRowAll,paste(folder,"outputRowHistogram.dat",sep=""),row.names=F,col.names=F)
			write.table(histogramColAll,paste(folder,"outputColHistogram.dat",sep=""),row.names=F,col.names=F)
			write.table(histogramFull,paste(folder,"outputFullHistogram.dat",sep=""),row.names=F,col.names=F)
			
			#Computing Frobenius by Mean FNN
			histNorm = histogramFull/max(histogramFull)
			A = histogramFull_before - histNorm 
			curr_norm = norm(A,type="F")
			frobNorm = c(frobNorm, curr_norm)
			cat("\t Frobenius Norm: ",curr_norm,"\n")
			
			histogramFull_before = histNorm
			i=i+1
			
		}
	
		#Average FNN results considering the color under analysis 
                if(max(meanResults) != 0) 
                        normMR = meanResults/max(meanResults) 
                if(max(histogramFull) != 0) 
                        normFH = histogramFull/max(histogramFull) 	
		
		#Summing results
		mean_normMR = mean_normMR + normMR
		mean_normFH = mean_normFH + normFH
		
	}

	#Computing the mean of all
	mean_normMR = mean_normMR/channels
	mean_normFH = mean_normFH/channels
	
	#Saving Files
	nameFile = paste(folder_initial,"meanResultsAllColors.dat",sep="")
	write.table(mean_normMR,nameFile,row.names=F,col.names=F)
	nameFile = paste(folder_initial,"histogramAllColors.dat",sep="")
	write.table(mean_normFH,nameFile,row.names=F,col.names=F)
	
	#Computing k best mask sizes
	maskRes = analise_mask_MR(folder_initial, k)
	
	#Printing and saving
	nameFile = paste(folder_initial,"bestMasks.dat",sep="")
	write.table(maskRes$maskKList[,1:2],nameFile,row.names=F,col.names=F)
	cat("\nBest Mask:\n")
	print(maskRes$maskKList[,1:2])
	
	#output
	ret = list()
	ret$dataset = dataset
	ret$folder_initial = folder_initial
	ret$bestMasks = maskRes
	
	ret
}

