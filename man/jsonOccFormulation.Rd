% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbifDownload.R
\name{jsonOccFormulation}
\alias{jsonOccFormulation}
\title{Create a JSON Query for GBIF's Occurrence API}
\usage{
jsonOccFormulation(
  ...,
  user = NULL,
  email = NULL,
  sendNotification = TRUE,
  geodistance.values = NULL
)
}
\arguments{
\item{...}{A set of queries to combine to create for GBIF's occurrence
download API (see details below).}

\item{user}{A character string containing a valid username associated with a
GBIF account. If \code{NULL} then \code{Sys.getenv('GBIF_USER')} and
\code{getOptions('gbif_user')} are searched for a value instead.}

\item{email}{A character string containing a valid email associated with a
GBIF account. If \code{NULL} then \code{Sys.getenv('GBIF_EMAIL')} and
\code{getOptions('gbif_email')} are searched for a value instead.}

\item{sendNotification}{A scalar logical that, if \code{TRUE}, adds a portion
to the query that requests that an email be sent to the account specified in
\code{email} when the download is ready.}

\item{geodistance.values}{Argument is only used when using a 'geodistance'
query (see details below).}
}
\value{
A \code{jsonlite} object containing a query to use in GBIF's
occurrence download API. This object may have the \code{"taxonFrame"}
attribute defined if a call to \code{\link{taxonLookup}} was made. If so,
this attribute will be set to the data frame resulting from the taxonomic
query.
}
\description{
A function to create a JSON query for the
\href{https://techdocs.gbif.org/en/data-use/api-downloads.html}{occurrence download API}
of \href{https://www.gbif.org/}{GBIF}.
}
\details{
Queries to be converted to JSON are specified using named
parameters. Any named parameter that shares a name with one of the parameters
listed in GBIF's
\href{https://techdocs.gbif.org/en/openapi/v1/occurrence#/Searching\%20occurrences/searchOccurrence}{occurrence search API}
will be interpreted as a filter on the relevant records. For example, having
an argument entitled \code{COUNTRY} with value \code{"NO"} will result a
query that will only download records that have an ISO-3166-1 country code of
NO (that is all record registered in Norway). If the vector provided for this
argument has a length greater than one, then all records that match any of
the composite values will be returned. These query types correspond to the
'equals' and 'in' predicates described in the
\href{https://techdocs.gbif.org/en/data-use/api-downloads.html}{occurrence download API document}
respectively. Other predicates can be specified by adding a suffix to the
named argument. For example:
\describe{
\item{\code{key.min}}{Will set a query that will only return records that
have a value for the \code{key} that are greater than or equal to the value
provided in the argument. For example, supplying \code{YEAR.min = 1960} as
an argument will result in a query that only returns records that have been
collected from the year 1960 onward. There is also the suffix \code{.xmin}
that can be used to create a query that will retrieve records only if they
have a relevant value for \code{key} that are greater than the value
provided in the argument.}
\item{\code{key.max}}{Will set a query that will only return records that
have a value for the \code{key} that are less than or equal to the value
provided in the argument. For example, supplying \code{YEAR.max = 1990} as
an argument will result in a query that only returns records that have been
collected in the year 1990 or earlier onward. There is also the suffix
\code{.xmax} that can be used to create a query that will retrieve records
only if they have a relevant value for \code{key} that are less than the
value provided in the argument.}
\item{\code{key.range}}{Will set a query that will only return records that
have a value for the \code{key} that are greater than or equal to the
minimum value and less than or equal to the maximum value provided in the
vector of values provided for this argument. For example, supplying
\code{ELEVATION.range = c(400, 600)} as an argument will result in a query
that only returns records that have been collected at an elevation of
between 400 and 600 metres above sea level.}
\item{\code{key.like}}{Will set a query that will only return records that
have a value for the \code{key} that match a pattern set in the value given
for the argument. Details on how the string provided for the pattern
matching is interpreted in the specification for the 'like' predicate in the
\href{https://techdocs.gbif.org/en/data-use/api-downloads.html}{occurrence download API document}.}
\item{\code{key.notnull}}{If the value of this argument is \code{TRUE} then
it will result in a query that will only include records that have a
non-empty value for \code{key}. Similarly a value of \code{FALSE} will
result in a query that will only include records that have an empty value
for \code{key}.}
\item{\code{key.not}}{This presence of this suffix represents the logical
negation for all components of the suffix that comes before it. For example,
supplying \code{DATASET_KEY.not = "4fa7b334-ce0d-4e88-aaae-2e0c138d049e"}
will result in a query that includes records from all datasets except that
with the provided dataset key. Similarly, supplying
\code{ELEVATION.range.not = c(400, 600)} as an argument will result in a
query that retrieves all records that do not fall within the elevational
band that lies between 400 and 600 metres above sea level.}
}
This function will accept arguments that are named either using the uppercase
notation used in the examples provided here or using the camel case notation
used in the \code{\link[rgbif]{pred}} function.

By default the arguments \code{HAS_GEOSPATIAL_ISSUE = FALSE},
\code{HAS_COORDINATE = TRUE}, \code{OCCURRENCE_STATUS = "PRESENT"}, and
\code{BASIS_OF_RECORD.not = c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")} are
added to the list of arguments provided by the user. These are commonly
applied defaults and are added to avoid accidental requests for erroneous
records. If the user does not want these default to be applied then they can
override this by providing \code{NULL} as the relevant argument. For example,
supplying \code{OCCURRENCE_STATUS = NULL} will result in a query that
requests records that represent both presences and absences.

If an argument with the name '\code{geometry}' appears in the list of query
arguments with a value represented by an object of class \code{sfc} (or that
inherits from it), then a query will be produced that only includes records
that fall within the polygon specified in the object. Similarly, if an
argument with the name '\code{geodistance}' appears in the list of query
arguments with a value with class that inherits from \code{sfc} than a query
will be produced that only included records that fall within a distance
(specified by the \code{geodistance.values} argument) of the coordinates (as
resulting from the \code{\link[sf]{st_coordinates}} function) that make up
the relevant spatial feature.

If any arguments are passed to the function with names that match those in
\code{\link{taxonLookup}} and/or \code{\link[rgbif]{name_lookup}} functions,
then a call will be made to \code{taxonLookup} using those arguments. The
resulting data frame will be added as a \code{"taxonFrame"} attribute to
the output. The keys associated with the finest taxonomic resolution in the
taxonomic lookup query will then be used as keys to refine the occurrence
record search. For example, providing the arguments
\code{species = "Aeshna grandis"} and \code{order = "Odonata"}, will result
in a data frame returned in the \code{"taxonFrame"} attribute of the output
that will include the column \code{"speciesKey"} representing the keys in
GBIF's taxonomic database that are possible entries for the species
Aeshna grandis in the order Odonata. These keys values will be used in a
query as though they were passed to the function using the argument notation
\code{SPECIES_KEY = vectorOfSpeciesKeys} and therefore ensuring that only
records corresponding to this species are requested.
}
\seealso{
\code{\link[rgbif]{name_lookup}}, \code{\link{taxonLookup}},
\code{\link{gbifOccDownload}}, \code{\link[sf]{st_sfc}},
\code{\link[sf]{st_coordinates}}
}
\author{
Joseph D. Chipperfield, \email{joechip90@googlemail.com}
}
