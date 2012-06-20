\name{getFileSepChar}
\alias{getFileSepChar}
\title{
getFileSepChar
}
\description{
Get the separator character for an existing readable text file.
}

\usage{
getFileSepChar(File)
}

\arguments{
  \item{File}{The name of the existing readable text file.}
}
\value{
  The character that separates data items in the file is returned. NULL is returned if the file cannot be read or if the algorithm does not recognize a separator character.
}
\author{ Niles Oien, CU CHANGE Lab  <Niles.Oien@colorado.edu> }
\examples{
 \dontrun{
	sepChar <- getFileSepChar( "/some/directory/structure/asciiFile.name" )
 }
}
\keyword{ hplot }
\keyword{ manip }
