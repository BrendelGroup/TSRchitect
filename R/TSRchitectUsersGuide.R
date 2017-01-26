#' @title \strong{TSRchitectUsersGuide}
#' @export

TSRchitectUsersGuide <- function(view=TRUE)
#	Find and optionally view TSRchitect User's Guide
{
	f <- system.file("doc","TSRchitectUsersGuide.pdf",package="TSRchitect")
	if(view) {
		if(.Platform$OS.type == "windows") 
			shell.exec(f)
		else
			system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
	}
	return(f)
}
