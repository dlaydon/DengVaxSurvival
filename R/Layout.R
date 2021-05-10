
MakeLayoutMatrix = function(NoCols, NoRows, AddTitle = TRUE, AddColNames = TRUE, AddRowNames = TRUE, BY_ROW = FALSE)
{
	if (!AddTitle)  			## just plots
	{
		LayoutMatrix = matrix(1:(NoCols*NoRows), ncol = NoCols, byrow = BY_ROW) 
		
	}	else if (!AddColNames)			## plots with Overall Plot Title
	{
		LayoutMatrix = matrix(1:(NoCols*NoRows) + 1, ncol = NoCols, byrow = BY_ROW) ### plus 1 as need 1st plot to be titles
		LayoutMatrix = rbind (rep(1, NoCols), LayoutMatrix)
		
	} 	else if (!AddRowNames)			## plots with Overall Plot Title and colnames
	{
		LayoutMatrix = matrix(1:(NoCols*NoRows) + 2, ncol = NoCols, byrow = BY_ROW) ### plus 2 as need 1st plot to be titles, 2nd plot needs to be plot name. 
		LayoutMatrix = cbind (rep(2, NoRows), LayoutMatrix)
		LayoutMatrix = rbind (rep(1, NoCols + 1), LayoutMatrix)
		
	}	else							## plots with Overall Plot Title and colnames and rownames
	{
		LayoutMatrix = matrix(1:(NoCols*(NoRows + 1)) + 2, ncol = NoCols, byrow = BY_ROW) ### plus 2 as need 1st plot to be titles, 2nd plot needs to be plot name. (NoRows + 1) as each colname will be a row.  
		LayoutMatrix = cbind (rep(2, (NoRows + 1)), LayoutMatrix)
		LayoutMatrix = rbind (rep(1, NoCols + 1), LayoutMatrix)
	}
	
	return(LayoutMatrix)
}

MakeHeightsVecForLayout = function(NoRows, AddTitle = TRUE, AddColNames = TRUE, AddRowNames = TRUE)
{
	if (!AddTitle)  			## just plots
	{
		Heights = NULL
		
	}	else if (!AddColNames)			## plots with Overall Plot Title
	{
		Heights = c(1,rep(4, NoRows))
		
	} 	else if (!AddRowNames)			## plots with Overall Plot Title and colnames
	{
		Heights = c(1, rep(4, NoRows))
		
	}	else							## plots with Overall Plot Title and colnames and rownames
	{
		Heights = c(1.75, 0.5, rep(4, NoRows))	
	}
	
	return(Heights)
}

MakeWidthsVecForLayout = function(NoCols, AddTitle = TRUE, AddColNames = TRUE, AddRowNames = TRUE, MarginToIndividualPlotRatio = 4)
{
	if (!AddTitle)  			## just plots
	{
		Widths = NULL
		
	}	else if (!AddColNames)			## plots with Overall Plot Title
	{
		Widths = NULL
		
	} 	else if (!AddRowNames)			## plots with Overall Plot Title and colnames
	{
		Widths = c(1, rep(MarginToIndividualPlotRatio, NoCols))
		
	}	else							## plots with Overall Plot Title and colnames and rownames
	{
		Widths = c(1, rep(MarginToIndividualPlotRatio, NoCols))
	}
	
	return(Widths)
}

MakeLayoutMatrix_2ndAttempt 		= function(NoCols, 	NoRows, 														AddTitle = TRUE, AddColNames = TRUE, 	AddRowNames = TRUE, BY_ROW = FALSE, RemoveCorner = TRUE)
{
	## Start by making a matrix. Then add space for colnames, rownames and full titles as necessary. 
	LayoutMatrix = matrix(1:(NoRows*NoCols), ncol = NoCols, byrow = BY_ROW)
	if (AddRowNames)
	{
		LayoutMatrix = LayoutMatrix + 1
		LayoutMatrix = cbind(rep(1, dim(LayoutMatrix)[1]), 		LayoutMatrix)
	}
	if (AddColNames) 
	{
		LayoutMatrix = LayoutMatrix + 1
		LayoutMatrix = rbind(rep(1, dim(LayoutMatrix)[2]), 		LayoutMatrix)
	}
	if (AddTitle) 
	{
		LayoutMatrix = LayoutMatrix + 1
		LayoutMatrix = rbind(rep(1, dim(LayoutMatrix)[2]), 		LayoutMatrix)
	}
	
	if (RemoveCorner)
		if (AddRowNames & AddColNames) ## remove weird "corner that ruins everyting in rest of plotting code (specifically position of text/rownames/colnames etc.) 
		{
			if (AddTitle) 	LayoutMatrix[2,1] = max(LayoutMatrix) + 1 	else  	### remove the corner (that is below the title) hence 2,1
				LayoutMatrix[1,1] = max(LayoutMatrix) + 1 			### remove the corner hence 1,1
			
		}
	return(LayoutMatrix)
}
MakeHeightsVecForLayout_2ndAttempt 	= function(			NoRows, Ratio_Plots_To_Title = 4, Ratio_Plots_To_ColNames = 2, 	AddTitle = TRUE, AddColNames = TRUE)
{
	Heights = rep(1, NoRows)	
	if (AddColNames) 	Heights = c(1/Ratio_Plots_To_ColNames	, Heights)
	if (AddTitle)  		Heights = c(1/Ratio_Plots_To_Title		, Heights)
	return(Heights)	
}
MakeWidthsVecForLayout_2ndAttempt 	= function(NoCols, 	Ratio_Plots_To_RowNames = 2, 																			AddRowNames = TRUE)
{
	Widths = rep(1, NoCols)
	if (AddRowNames) Widths = c(1/Ratio_Plots_To_RowNames, Widths)
	return(Widths)	
}
SetUpMultiPlot 						= function(NoCols = NULL, NoRows = NULL, ColNames = NULL, RowNames = NULL, MultiPlot_Title = NULL, AddTitle = TRUE, AddColNames = TRUE, AddRowNames = TRUE, BY_ROW = FALSE, RemoveCorner = TRUE, 
		Ratio_Plots_To_Title = 4, Ratio_Plots_To_ColNames = 2, Ratio_Plots_To_RowNames = 2, ColName_Cex = 1.4, RowName_Cex = 1.4, MultiPlot_Title_Cex = 2, OG_Mar = OrigMAR, modelrun = ModelRun)
{
	if (is.null(NoCols)) NoCols = length(ColNames)
	if (is.null(NoRows)) NoRows = length(RowNames)
	if (!is.null(MultiPlot_Title)) AddTitle 	= TRUE
	if (!is.null(ColNames		)) AddColNames 	= TRUE
	if (!is.null(RowNames		)) AddRowNames 	= TRUE
	
	LayoutMatrix 	= MakeLayoutMatrix_2ndAttempt(NoCols = NoCols, NoRows = NoRows, AddTitle = AddTitle, AddColNames = AddColNames, AddRowNames = AddRowNames, BY_ROW = BY_ROW, RemoveCorner = RemoveCorner)
	HEIGHTS 		= MakeHeightsVecForLayout_2ndAttempt(NoRows = NoRows, Ratio_Plots_To_Title = Ratio_Plots_To_Title, Ratio_Plots_To_ColNames = Ratio_Plots_To_ColNames, AddTitle = AddTitle, AddColNames = AddColNames)
	WIDTHS			= MakeWidthsVecForLayout_2ndAttempt	(NoCols = NoCols, Ratio_Plots_To_RowNames = Ratio_Plots_To_RowNames, AddRowNames = AddRowNames)
	layout(LayoutMatrix, heights = HEIGHTS, widths = WIDTHS)
	
	if (AddTitle)
	{	
		par(mar = rep(0, 4)) ### without this you often get the annoying "Error in plot.new() : figure margins too large" error. Set margins to zero for the title, then reset them to OrigMAR (global variable defined at start of script). 
		plot.new()
		if (is.null(MultiPlot_Title))
		{
			if (ShortenOutputString) MultiPlot_Title = AbbreviateOutputString(OutputSubDir, modelrun) else MultiPlot_Title = OutputSubDir
			#MultiPlot_Title = paste0(MultiPlot_Title, "\n", AGs_String, Countries_String, DiseaseSeverities_String, Phases_String)
		}
		text(0.5, 0.5, MultiPlot_Title , cex = MultiPlot_Title_Cex, font = 2)
	}
	if (AddColNames)
	{
		par(mar = rep(0, 4))
		plot(NA, ylim = 0:1, xlim = c(0,length(ColNames)) , yaxs='i', xaxs='i', xaxt = "n", yaxt = "n", bty="n")
		text(1:length(ColNames) - 0.5 , 0.5, ColNames, cex = ColName_Cex, font = 2)
	}
	if (AddRowNames)
	{
		par(mar = rep(0, 4))
		plot(NA, xlim = 0:1, ylim = c(0,length(RowNames)), yaxs='i', xaxs='i', xaxt = "n", yaxt = "n", bty="n")
		text(0.5, (1:length(RowNames) - 0.5), rev(RowNames), cex = RowName_Cex, font = 2, srt = 90) ### don't know why you have to reverse the rownames, but it otherwise the names come out in the wrong order. 
	}
	
	# reset par ("mar")
	par(mar = OG_Mar) ### reset default margins or S_Curves will look horrid. 
}



#MakeWidthsVecForLayout_ORIG = function(NoCols, AddTitle = TRUE, AddColNames = TRUE, AddRowNames = TRUE)
#{
#	if (!AddTitle)  			## just plots
#	{
#		Widths = NULL
#		
#	}	else if (!AddColNames)			## plots with Overall Plot Title
#	{
#		Widths = NULL
#		
#	} 	else if (!AddRowNames)			## plots with Overall Plot Title and colnames
#	{
#		Widths = c(1, rep(4, NoCols))
#		
#	}	else							## plots with Overall Plot Title and colnames and rownames
#	{
#		Widths = c(1, rep(4, NoCols))
#	}
#	
#	return(Widths)
#}

