#' Preprocesses a Cytotrol flowFrame object
#'
#' Our goal here is to use swap the marker names and the channel names within a
#' \code{flowFrame} object to ensure that the \code{flowFrame} objects across
#' centers can be merged into a single \code{flowSet}.
#'
#' We also preprocess the marker names to strip out any additional information
#' added to the marker name. For instance, NHLBI uses "IgD V500", which we reduce
#' to "IgD".
#'
#' @param flow_frame the \code{flowFrame} object to preprocess
#' @param markers_keep a character vector containing the markers to keep
#' @return the updated \code{flowFrame} object containing only the markers of
#' interest
preprocess_flowframe <- function(flow_frame, markers_keep) {
  if (missing(markers_keep)) {
    stop("The marker to keep must be specified.")
  }

  fr_rownames <- rownames(parameters(flow_frame)@data)
  
  # Preprocesses each of the columns in the flow_frame
  for (j in seq_len(ncol(flow_frame))) {
    marker_idx <- paste0(fr_rownames[j], "S")
    channel_idx <- paste0(fr_rownames[j], "N")

    marker <- flow_frame@description[[marker_idx]]
    channel <- flow_frame@description[[channel_idx]]

    # In the case the marker name is given, we swap the marker and channel
    # names.
    if (!is.null(marker) && channel != "<FITC-A>") {
      # Converts the marker names to a common name
      marker <- marker_conversion(marker)

      # Updates the channel information in the flow_frame with the marker
      flow_frame@description[[channel_idx]] <- marker
      flow_frame@parameters@data$name[j] <- marker

      # Updates the marker information in the flow_frame with the channel
      flow_frame@description[[marker_idx]] <- channel
      flow_frame@parameters@data$desc[j] <- channel
    } else if (is.null(marker) && channel == "FITC-A") {
      marker <- "Live"
      channel <- "FITC-A"
      flow_frame@description[[channel_idx]] <- marker
      flow_frame@parameters@data$name[j] <- marker

      # Updates the marker information in the flow_frame with the channel
      flow_frame@description[[marker_idx]] <- channel
      flow_frame@parameters@data$desc[j] <- channel
    }
  }
  colnames(exprs(flow_frame)) <- colnames(flow_frame)
  
  # Subset to markers of interest
  flow_frame <- flow_frame[, markers_keep]

  # The pData for the parameters are (sometimes?) resulting in class "AsIs"
  # rather than "character", which is causing errors in the gating.
  # We fix this here.
  pData(parameters(flow_frame))$name <- as.character(pData(parameters(flow_frame))$name)
  pData(parameters(flow_frame))$desc <- as.character(pData(parameters(flow_frame))$desc)

  flow_frame
}

#' Converts the Cytotrol marker names to a common name
#'
#' For the following list of marker names, we manually update the names so
#' that they are standard across centers.
marker_conversion <- Vectorize(function(marker) {
  # If marker name contains additional info, remove everything after the
  # space. (e.g., "IgD V500" to "IgD")
  marker <- strsplit(marker, " ")[[1]][1]

  if (marker == "19") {
    marker <- "CD19"
  } else if (marker %in% c("LIVE", "LIVE_GREEN", "Live/Dead", "live", "Live/green")) {
    marker <- "Live"
  } else if (marker == "IGD") {
    marker <- "IgD"
  } else if (marker %in% c("HLA", "HLADR", "HLA-DR")) {
    marker <- "HLADR"
  } else if (marker == "CD197") {
    marker <- "CCR7"
  } else if (marker == "CD194") {
    marker <- "CCR4"
  } else if (marker == "CD11C") {
    marker <- "CD11c"
  } else if (marker %in% c("CD3CD19CD20", "CD3+19+20", "CD3_CD19_CD20",
                           "CD3+CD19+CD20+", "CD3+CD19+CD20", "CD3+19+20", "CD3/19/20")) {
    marker <- "Lineage"
  } else if (marker == "CD196") {
    marker <- "CCR6"
  } else if (marker == "CD183") {
    marker <- "CXCR3"
  }

  marker
})
