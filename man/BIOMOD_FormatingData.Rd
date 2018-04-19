%\renewcommand{\usage}{\pkg{#1}}

\name{BIOMOD_FormatingData}
\alias{BIOMOD_FormatingData}

\title{ Initialize the datasets for usage in \pkg{biomod2} }

\description{
This function rearranges the user's input data to make sure they can be used within \pkg{biomod2}.
The function allows to select pseudo-absences or background data in the case that true absences data are not available, or to add pseudo-absence data to an existing set of absence  (see details).
}
\usage{
BIOMOD_FormatingData(resp.var,
                     expl.var,
                     resp.xy = NULL,
                     resp.name = NULL,
                     eval.resp.var = NULL,
                     eval.expl.var = NULL,
                     eval.resp.xy = NULL,
                     PA.nb.rep = 0,
                     PA.nb.absences = 1000,
                     PA.strategy = 'random',
                     PA.dist.min = 0,
                     PA.dist.max = NULL,
                     PA.sre.quant = 0.025,
                     PA.table = NULL,
                     na.rm = TRUE)
}

\arguments{

  \item{resp.var}{ a vector, \code{\link[sp]{SpatialPointsDataFrame}} (or \code{\link[sp]{SpatialPoints}} if you work with \sQuote{only presences} data) containing species data (a single species) in binary format (ones for presences, zeros for true absences and NA for indeterminate ) that will be \bold{used to build the species distribution models}. }

  \item{expl.var}{ a \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} or \code{\link[raster:stack]{RasterStack}} containing your explanatory variables that will be \bold{used to build your models}.}

  \item{resp.xy}{ optional 2 columns \code{matrix} containing the X and Y coordinates of resp.var (only consider if resp.var is a vector) that will be \bold{used to build your models}.}

  \item{eval.resp.var}{ a vector, \code{\link[sp]{SpatialPointsDataFrame}} your species data (a single species) in binary format (ones for presences, zeros for true absences and NA for indeterminate ) that will be \bold{used to evaluate the models with independent data (or past data for instance)}.}

  \item{eval.expl.var}{a \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} or \code{\link[raster:stack]{RasterStack}} containing your explanatory variables that will be \bold{used to evaluate the models with independent data (or past data for instance)}.}

  \item{eval.resp.xy}{optional 2 columns \code{matrix} containing the X and Y coordinates of resp.var (only consider if resp.var is a vector) that will be \bold{used to evaluate the modelswith independent data (or past data for instance)}.}

  \item{resp.name}{ response variable name (character). The species name. }

  \item{PA.nb.rep}{ number of required Pseudo Absences selection (if needed). 0 by Default.}

  \item{PA.nb.absences}{ number of pseudo-absence selected for each repetition (when PA.nb.rep > 0) of the selection (true absences included)}

  \item{PA.strategy}{ strategy for selecting the Pseudo Absences (must be \sQuote{random}, \sQuote{sre}, \sQuote{disk} or \sQuote{user.defined})}

  \item{PA.dist.min}{minimal distance to presences for \sQuote{disk} Pseudo Absences selection (in meters if the explanatory is a not projected raster (+proj=longlat) and in map units (typically also meters) when it is projected or when explanatory variables are stored within table )}

  \item{PA.dist.max}{maximal distance to presences for \sQuote{disk} Pseudo Absences selection(in meters if the explanatory is a not projected raster (+proj=longlat) and in map units (typically also meters) when it is projected or when explanatory variables are stored within table ) }

  \item{PA.sre.quant}{quantile used for \sQuote{sre} Pseudo Absences selection}

  \item{PA.table}{a \code{matrix} (or a \code{data.frame}) having as many rows than \code{resp.var} values. Each column corresponds to a Pseudo-absences selection. It contains \code{TRUE} or \code{FALSE} indicating which values of \code{resp.var} will be considered to build models. It must be used with \sQuote{user.defined} \code{PA.strategy}.}

  \item{na.rm}{logical, if TRUE, all points having one or several missing value for environmental data will be removed from the analysis}
}

\details{
  This function homogenizes the initial data for making sure the modelling exercise will be completed with all the required data. It supports different kind of inputs.

  IMPORTANT: When the explanatory data are given in \code{rasterLayer} or \code{rasterStack} objects, \pkg{biomod2} will be extract the variables onto the XY coordinates of the presence (and absence is any) vector. Be sure to give the XY coordinates (\sQuote{resp.xy}) in the same projection system than the raster objects. Same for the evaluation data in the case some sort of independent (or past) data are available (\sQuote{eval.resp.xy}).
  When the explanatory variables are given in \code{\link[sp]{SpatialPointsDataFrame}}, the same requirements are asked than for the raster objects. The XY coordinates must be given to make sure biomod2 can extract the explanatory variables onto the presence (absence) data
  When the explanatory variables are stored in a data.frame, make sure there are in the same order than the response variable. \pkg{biomod2} will simply merge the datasets without considering the XY coordinates.

  When both presence and absence data are available, and there is enough absences: set sQuote{PA.nb.rep} to 0. No pseudo-absence will be extracted.

  When no true absences are given or when there are not numerous enough. It's advise to make several pseudo absences selections. That way the influence of the pseudo-absence selection could then be estimated later on. If the user do not want to run several repetition, make sure to select a relatively high number pseudo-absence. Make sure the number of pseudo-absence data is not higher than the maximum number of potential pseudo-absence (e.g. do not select 10,000 pseudo-absence when the rasterStack or data.frame do not contain more than 2000 pixels or rows).

  \enumerate{
    \item{\bold{Response variable encoding}}{

      \code{BIOMOD_FormatingData} concerns a single species at a time so \code{resp.var} must be a uni-dimensional object.

      Response variable must be a \code{vector} or  a one column \code{data.frame}/\code{matrix}/\code{\link[sp]{SpatialPointsDataFrame}} ( \code{\link[sp]{SpatialPoints}} are also allowed if you work with \sQuote{only presences} data) object.
      As most of \pkg{biomod2} models need Presences AND Absences data, the response variable must contain some absences (if there are not, make sure to select pseudo-absence). In the input \code{resp.var} argument, the data should be coded in the following way :
      \itemize{
        \item{Presences : 1}
        \item{True Absences : 0 (if any)}
        \item{No Information : NA (if any, might latter be used for pseudo-absence)}
      }

      If \code{resp.var} is a non-spatial object (\code{vector}, \code{matrix}/\code{data.frame}) and that some models requiring spatial data are being used (e.g. MAXENT.Phillips) and/or pseudo absences spatially dependent (i.e 'disk'), make sure to give the XY coordinates of the sites/rows (\sQuote{resp.xy}).

    }

    \item{\bold{Explanatory variables encoding}}{

      Explanatory variables must be stored together in a multi-dimensional object. It may be a \code{matrix},  a \code{data.frame}, a \code{\link[sp]{SpatialPointsDataFrame}} or a \code{rasterStack} object. Factorial variables are allowed here even if that can lead to some models omissions.
    }

    \item{\bold{Evaluation Data}}{

      If you have data enough, we strongly recommend to split your dataset into 2 part : one for training/calibrating and testing the models and another to evaluate it. If you do it, fill the \code{eval.resp.var}, \code{eval.expl.var} and optionally the \code{eval.resp.xy} arguments with this data. The advantage of working with a specific dataset for evaluating your models is that you will be able to evaluate more properly your \sQuote{ensemble modeled} models. That being said, this argument is optional and you may prefer only to test (kind of evaluation) your models only with a \sQuote{cross-validation} procedure (see Models function). The best practice is to use one set of data for training/calibrating, one set of testing and one for evaluating. The calibration and testing of the data can be done automatically in \pkg{biomod2} in the Models function. The dataset for evaluation must be entered in \code{BIOMOD_FormatingData}.

    }

    \item{\bold{Pseudo Absences selection}}{

      The \code{PA.xxx}'s arguments let you parameterize your pseudo absences selection if you want some. It's an optional step.

      Pseudo absences will be selected within the \sQuote{background data} and might be constrained by a defined \sQuote{strategy}.

      \enumerate{
        \item{background data}{

        \sQuote{Background data} represents data there is no information whether the species of interest occurs or not. It is defined by the \sQuote{No Information} data of your \code{resp.var} if you give some. If not, (i.e Only presences data or all cells with a define presence or absence state) the background will be take into your \code{expl.var} object if it's a \code{RasterStack}.

        }

        \item{strategy}{

          The strategy allows to constrain the choice of pseudo-absence within the \sQuote{background data}.
          3 ways are currently implemented to select the pseudo-absences candidate cells (\code{PA.strategy} argument):
          \itemize{
            \item{\sQuote{random}: all cell of initial background are Pseudo absences candidates. The choice is made randomly given the number of pseudo-absence to select \code{PA.nb.absences}.}

            \item{\sQuote{disk}: you may define a minimal (\code{PA.dist.min}), respectively a maximal (\code{PA.dist.max}) distance to presences points for selecting your pseudo absences candidates. That may be useful if you don't want to select pseudo-absences too close to your presences (same niche and to avoid pseudo-replication), respectively too far from your presences (localized sampling strategy). }

            \item{\sQuote{sre}: Pseudo absences candidates have to be selected in condition that differs from a defined proportion (\code{PA.sre.quant}) of presences data. It forces pseudo absences to be selected outside of the broadly defined environmental conditions for the species. It means that a surface range envelop model (sre, similar the BIOCLIM) is first carried out (using the specified quantile) on the species of interest, and then the pseudo-absence data are extracted outside of this envelop. This particular case may lead to over optimistic models evaluations.}

            \item{\sQuote{user.defined}: In this case, pseudo absences selection should have been done in a previous step. This pseudo absences have to be reference into a well formatted \code{data.frame} (e.g. \code{PA.table} argument)}
          }

        }

      }



    }
  } % end enumerate

} %end detail section

\value{
  A \code{'data.formatted.Biomod.object'} for \code{\link[biomod2]{BIOMOD_Modeling}}.
  It is strongly advised to check whether this formatted data corresponds to what was expected. A summary is easily printed by simply tipping the name of the object. A generic plot function is also available to display the different dataset in the geographic space.
}

\author{ Wilfried Thuiller, Damien Georges }

\seealso{ \code{\link{BIOMOD_Modeling}}}
\examples{

# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd",
                     package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio7.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio11.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio12.grd",
                             package="biomod2"))
# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

myBiomodData
plot(myBiomodData)

}

\keyword{ models }
\keyword{ datasets }
