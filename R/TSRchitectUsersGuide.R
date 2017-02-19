#' @title \strong{TSRchitectUsersGuide}
#' @description Opens the TSRchitect User's Guide in a pdf viewer
#' on the user's system.
#'
#' @param view (logical) if TRUE (default) the User's Guide
#' is opened in the local pdf viewer.
#' If FALSE then the full path to the User's Guide is
#' returned.
#'
#' @return if view=FALSE, then the full path to the User's Guide
#' is returned.
#'
#' @examples
#' myPath <- TSRchitectUsersGuide(view=FALSE)
#'
#' @export


TSRchitectUsersGuide <- function(view=TRUE)
#	Find and optionally view TSRchitect User's Guide
{
	f <- system.file("doc", "TSRchitectUsersGuide.pdf", package="TSRchitect")
	if(view) {
		if(.Platform$OS.type == "windows")
			shell.exec(f)
		else
			system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
	}
	return(f)
}
