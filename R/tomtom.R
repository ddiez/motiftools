getTomtomAlphabet <- function(x) {
  tmp <- lapply(xml_attrs(xml_find_all(x, "model/alphabet/letter")), function(l) {
    data.frame(
      id = l["id"],
      symbol = l["symbol"],
      aliases = ifelse(is.na(l["aliases"]), l["symbol"], l["aliases"]),
      color = l["colour"],
      equals = ifelse(is.na(l["equals"]), l["symbol"], l["equals"]),
      name = l["name"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, tmp)
}

getTomtomMotifProbabilities <- function(x, from = c("queries", "targets")) {
  path <- paste0(from, "/motif")
  lapply(xml_find_all(x, path), function(motif) {
    m <- xml_attrs(xml_find_all(motif, "pos"))
    m <- do.call(rbind, m)
    mode(m) <- "numeric"
    rownames(m) <- 1:nrow(m)
    t(as.matrix(m))
  })
}

getTomtomMotifInfo <- function(x, from = c("queries", "targets")) {
  path <- paste0(from, "/motif")
  tmp <- do.call(rbind, xml_attrs(xml_find_all(x, path)))  
  data.frame(
    db = tmp[, "db"],
    motif_id = tmp[, "id"],
    alt = tmp[, "alt"],
    width = as.numeric(tmp[, "length"]),
    nsites = as.numeric(tmp[, "nsites"]),
    e_value = as.numeric(tmp[, "evalue"]),
    stringsAsFactors = FALSE
  )
}

getTomtomMotifMatches <- function(x) {
  tmp <- lapply(xml_find_all(x, "matches/query"), function(motif) {
    query_id <- xml_attr(motif, "idx")
    tmp <- xml_attrs(xml_find_all(motif, "target"))
    tmp <- do.call(rbind, tmp)
    data.frame(query_id = query_id,
               target_id = tmp[, "idx"],
               offset = as.integer(tmp[, "off"]),
               p_value = as.numeric(tmp[, "pv"]),
               e_value = as.numeric(tmp[, "ev"]),
               q_value = as.numeric(tmp[, "qv"]),
               stringsAsFactors = FALSE
    )
  })
  tmp <- do.call(rbind, tmp)
  rownames(tmp) <- NULL
  tmp
}

#' readTOMTOM
#' 
#' Read results from tomtom.
#'
#' @param file Tomtom file in XML format.
#' @param description description for the experiment.
#'
#' @return MotifCompareResult
#' @export
#'
#' @examples
#' NULL
readTOMTOM <- function(file, description = NULL) {
  doc <- read_xml(file)
  root <- xml_root(doc)
  
  alpha_info <- getTomtomAlphabet(root)
  query_prob <- getTomtomMotifProbabilities(root, "queries")
  target_prob <- getTomtomMotifProbabilities(root, "targets")
  
  query_info <- getTomtomMotifInfo(root, "queries")
  target_info <- getTomtomMotifInfo(root, "targets")
  
  prob_info <- list(
    query = query_prob,
    target = target_prob
  )
  
  match_info <- getTomtomMotifMatches(root)
  
  alphabetData <- AnnotatedDataFrame(alpha_info)
  
  new("MotifCompareResult",
      info = list(
        tool = "tomtom",
        description = description,
        nquery = nrow(query_info),
        ntarget = nrow(target_info)
      ),
      alphabet = alpha_info,
      query = query_info,
      target = target_info,
      probabilities = prob_info,
      matches = match_info
  )
}

#' plotMotifMatches
#' 
#' Plot a matrix with rows and columns representing motifs in the query and target databases, and the fill
#' color representing one of the three statistics (p_value, e_value, q_value) measuring the significance of
#' the similarity between the motifs.
#'
#' @param x a MotifCompareResult object.
#' @param fill the statistic to plot. One of p_value, e_value, q_value (default: p_value).
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
plotMotifMatches <- function(x, fill = c("p_value", "e_value", "q_value")) {
  x <- x@matches
  ggplot(x, aes_string(x = "query_id", y = "target_id", fill = fill)) + 
    geom_tile() + 
    viridis::scale_fill_viridis(guide = guide_legend(), direction = -1) + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(
      aspect.ratio = 1,
      axis.ticks = element_blank()
    )
}

