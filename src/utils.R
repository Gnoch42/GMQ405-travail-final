# =============================================================================
# src/utils.R
# -----------------------------------------------------------------------------
# Petits utilitaires partagés (téléchargement robuste, opérateurs pratiques).
# =============================================================================

# Opérateur "valeur par défaut si NULL" (pratique pour lire le YAML).
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Téléchargement ROBUSTE d'un gros fichier (curl : reprise + réessais)
#'
#' download.file() de base est peu fiable pour les gros fichiers (~1 Go) des
#' serveurs gouvernementaux : délai court, troncature silencieuse. On délègue à
#' `curl` avec reprise (-C -), réessais et suivi des redirections, puis on
#' vérifie que la taille finale correspond à l'en-tête Content-Length du serveur.
#'
#' @param url  URL à télécharger.
#' @param dest Chemin de destination local.
#' @return (invisible) le chemin `dest` ; s'arrête en erreur si le fichier est
#'   incomplet ou si curl échoue.
robust.download <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  status <- system2("curl", c("-fL", "-C", "-", "--retry", "4", "--retry-delay", "5",
                              "-o", shQuote(dest), shQuote(url)))
  if (status != 0 || !file.exists(dest))
    stop("Échec du téléchargement : ", url, " (code curl ", status, ")")
  # Vérification de complétude via l'en-tête Content-Length.
  hdr <- suppressWarnings(system2("curl", c("-sIL", shQuote(url)), stdout = TRUE))
  cl  <- grep("(?i)content-length", hdr, value = TRUE, perl = TRUE)
  if (length(cl) > 0) {
    expected <- as.numeric(sub("\\D+", "", tail(cl, 1)))
    if (!is.na(expected) && file.info(dest)$size < expected)
      stop("Fichier incomplet (", file.info(dest)$size, " < ", expected, " o) : ", url)
  }
  invisible(dest)
}
