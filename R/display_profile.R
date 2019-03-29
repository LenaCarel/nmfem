#' Display 3D profiles
#'
#' This function display profiles of 3 dimensions (day, hour, number of observations). It has been created to display
#' profiles from the \code{nmfem} package data.
#'
#' @param profile a vector or a matrix row containing the profile to display. The day/hour data are contained in the column names.
#' @param numclient logical. Whether the first value of the row is an identifier.
#' @param color color of the display. Possibilities are the ones provided by \url{colorbrewer2.org}.
#' @param language in which language the day/hour names are written. For now, the possibilities are "en" for english and "fr" for french.
#' @param theme A theme to use. The only valid values are "" and "dark".
#' @return Creates a 3D-heatmap displayed in the Viewer tab.
#' @examples
#' display_profile(travelers[sample(nrow(travelers),1), ], numclient = TRUE)
#' @importFrom plyr mapvalues
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @import dplyr
#' @importFrom d3heatmap d3heatmap
#' @export
#'

display_profile <- function(profile, numclient = FALSE, color = "Blues", language = "en", theme = "dark"){

  heure <- NULL
  numjour <- NULL
  temps <- NULL

  if(numclient){
    profile <- profile[ ,-1] # Delete card ID
  }

  profile <- as.data.frame(profile)
  profile <- tidyr::gather(
    data = profile,
    key = temps,
    value = n
  )

  profile$jour  <- substr(profile$temps, 1, nchar(profile$temps) - 2)
  profile$heure <- substr(profile$temps, nchar(profile$temps) - 1, nchar(profile$temps))

  profile <- profile[ ,-1]
  profile <- tidyr::spread(
    data  = profile,
    key   = heure,
    value = n
  )

  if(toupper(language) == "FR"){
    profile$numjour = plyr::mapvalues(profile$jour, from = c("lundi", "mardi", "mercredi",
                                                             "jeudi", "vendredi", "samedi",
                                                             "dimanche"),
                                      to = c(1:7)
    ) %>% as.numeric()
  }else{
    profile$numjour = plyr::mapvalues(profile$jour, from = c("Monday", "Tuesday", "Wednesday",
                                                             "Thursday", "Friday", "Saturday",
                                                             "Sunday"),
                                      to = c(1:7)
    ) %>% as.numeric()
  }

  profile <- profile[order(profile$numjour), ]
  profile <- subset(profile, select = -numjour)

  profile[is.na(profile)] <- 0

  rownames(profile) <- profile$jour
  profile <- profile[ ,-1]

  d3heatmap::d3heatmap(profile,
                       dendrogram = "none", colors = color,
                       show_grid = 3, theme = theme)
}

