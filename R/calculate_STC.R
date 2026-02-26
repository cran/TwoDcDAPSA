#' Calculate TwoDcDAPSA: Swelling–Tenderness Contrast (STC) Score and Quartile Groups
#'
#' Computes STC (PC2) as a loading-weighted combination of standardized residuals
#' for Pain, Patient Global, SJC, and TJC after adjusting each for cDAPSA
#' via a natural spline model. STC emphasizes swelling relative to tenderness
#' under the default PC2 loadings.
#'
#' Includes input "tuning": coerces character columns to numeric (warning if NAs introduced)
#' and checks for out-of-range values (SJC 0-66, TJC 0-68, Pain/Patient Global 0-10)
#' with configurable handling. If \code{cDAPSA} is not provided, it is computed as
#' \code{SJC + TJC + Pain + Patient_Global}. If \code{cDAPSA} is provided, it is
#' verified against this sum (within \code{cdapsa_tolerance}); any discrepancy results
#' in an error.
#'
#' The function returns both 4-level quartiles (\code{STC_quartile}: Q1–Q4) and a
#' 3-level grouped factor (\code{STC_quartile_combine}) that combines Q1 and Q4
#' (levels: \code{Q1&4}, \code{Q2}, \code{Q3}).
#'
#' @param data A data.frame/tibble with the required columns.
#' @param cohort_id Name of the cohort id column.
#' @param cDAPSA Optional. Name of the cDAPSA column. If \code{NULL} (default),
#'   cDAPSA is computed as \code{SJC + TJC + Pain + Patient_Global}.
#'   If non-\code{NULL}, the provided column is verified to equal that sum
#'   within \code{cdapsa_tolerance}; otherwise an error is thrown.
#' @param Pain Name of the Pain column (0-10).
#' @param Patient_Global Name of the Patient Global column (0-10).
#' @param SJC Name of the Swollen Joint Count column (0-66).
#' @param TJC Name of the Tender Joint Count column (0-68).
#' @param oob_action What to do when an input is out of its valid range
#'   (SJC 0-66, TJC 0-68, Pain/Patient Global 0-10). One of:
#'   \code{"stop"} (error), \code{"na"} (keep rows but set STC/Quartile to NA),
#'   or \code{"drop"} (remove rows). Default is \code{"stop"}.
#' @param cdapsa_tolerance Numeric tolerance for comparing provided cDAPSA to
#'   the computed sum; default \code{1e-8}.
#' @param center_scale List of centers/scales used to standardize inputs.
#' @param ns_knots Numeric vector of interior knots for the spline on standardized cDAPSA.
#' @param ns_boundary_knots Numeric vector of boundary knots for the spline on standardized cDAPSA.
#' @param coef_list Named list of regression coefficients (intercept + 4 spline basis) for each component.
#' @param stc_loadings Numeric loadings (length 4) for Pain, Patient Global, SJC, TJC residuals.
#' @param stc_cutoffs Numeric vector of 5 cut points to define 4 quartile bins (include.lowest=TRUE).
#' @param resid_center_scale List with `center` and `scale` vectors for standardizing residuals.
#'
#' @return A tibble with \code{cohort_id}, \code{STC}, \code{STC_quartile},
#'   and \code{STC_quartile_combine}.
#'
#' @seealso \code{\link{calculate_PJC}} for the PROs-Joint Contrast (PJC; PC1).
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # Minimal example WITHOUT a cDAPSA column (it will be computed as SJC+TJC+Pain+PG)
#' df1 <- data.frame(
#'   id = 1:3,
#'   pain = c(4, 6, 8),
#'   pg   = c(3, 7, 9),
#'   sjc  = c(1, 3, 5),
#'   tjc  = c(0, 2, 4)
#' )
#' calculate_STC(
#'   df1,
#'   cohort_id = "id",
#'   cDAPSA = NULL,
#'   Pain = "pain",
#'   Patient_Global = "pg",
#'   SJC = "sjc",
#'   TJC = "tjc",
#'   oob_action = "na"
#' )
#'
#' # Example WITH a consistent cDAPSA column (verified against the sum)
#' df2 <- transform(df1, cdapsa = pain + pg + sjc + tjc)
#' calculate_STC(
#'   df2,
#'   cohort_id = "id",
#'   cDAPSA = "cdapsa",
#'   Pain = "pain",
#'   Patient_Global = "pg",
#'   SJC = "sjc",
#'   TJC = "tjc"
#' )
#'


calculate_STC <- function(
    data,
    cohort_id = "cohort_id",
    cDAPSA = NULL,
    Pain = "Pain",
    Patient_Global = "Patient_Global",
    SJC = "SJC",
    TJC = "TJC",
    oob_action = c("stop", "na", "drop"),
    cdapsa_tolerance = 1e-8,
    center_scale = list(
      Pain = c(center = 4.303191, scale = 2.798819),
      Patient_Global = c(center = 4.795213, scale = 2.791098),
      SJC = c(center = 3.783245, scale = 4.707089),
      TJC = c(center = 5.194149, scale = 7.371234),
      cDAPSA = c(center = 18.0758, scale = 14.03964)
    ),
    ns_knots = c(-0.7176679, -0.2190796, 0.4219626),
    ns_boundary_knots = c(-1.287483, 5.265392),
    coef_list = list(
      Pain = c(-1.48889, 1.93539, 2.25211, 3.35687, 2.68578),
      Patient_Global = c(-1.72890, 2.13640, 2.35881, 3.95251, 2.66605),
      SJC = c(-0.76905, 0.47397, 1.95020, 4.45945, 5.98404),
      TJC = c(-0.74115, 0.27891, 2.50892, 4.68559, 6.15326)
    ),
    stc_loadings = c(-0.05043539, -0.01476495, 0.73876801, -0.67190780),
    stc_cutoffs  = c(-Inf, -0.3881481, -0.1027692, 0.4411657, Inf),
    resid_center_scale = list(
      center = c(
        Pain = 1.155879e-15,
        "Patient Global" = 9.679019e-16,
        "Swollen Joint Count" = -2.764596e-15,
        "Tender Joint Count" = -3.534933e-15
      ),
      scale = c(
        Pain = 0.6478511,
        "Patient Global" = 0.6282206,
        "Swollen Joint Count" = 0.5895540,
        "Tender Joint Count" = 0.3902453
      )
    )
) {
  oob_action <- match.arg(oob_action)

  # -- helpers -------------------------------------------------------------
  coerce_if_char <- function(x, nm) {
    if (is.character(x)) {
      y <- suppressWarnings(as.numeric(x))
      n_bad <- sum(is.na(y) & !is.na(x))
      if (n_bad > 0) warning(sprintf("%s: NAs introduced by coercion for %d rows.", nm, n_bad), call. = FALSE)
      return(y)
    }
    x
  }

  # -- basic validation ----------------------------------------------------
  if (length(stc_cutoffs) != 5L || any(!is.finite(stc_cutoffs[2:4])) ||
      is.unsorted(stc_cutoffs, strictly = TRUE)) {
    stop("`stc_cutoffs` must be 5 strictly increasing values; outer values may be -Inf/Inf.", call. = FALSE)
  }

  # Columns we definitely need present
  base_required <- c(cohort_id, Pain, Patient_Global, SJC, TJC)
  missing_cols <- base_required[!base_required %in% names(data)]
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # -- tuning: coerce characters -> numeric -------------------------------
  dat <- data
  dat[[Pain]]            <- coerce_if_char(dat[[Pain]], Pain)
  dat[[Patient_Global]]  <- coerce_if_char(dat[[Patient_Global]], Patient_Global)
  dat[[SJC]]             <- coerce_if_char(dat[[SJC]], SJC)
  dat[[TJC]]             <- coerce_if_char(dat[[TJC]], TJC)
  if (!is.null(cDAPSA) && cDAPSA %in% names(dat)) {
    dat[[cDAPSA]] <- coerce_if_char(dat[[cDAPSA]], cDAPSA)
  }

  # -- Out-of-range handling (vectorised) ---------------------------------
  ranges  <- list(Pain = c(0,10), Patient_Global = c(0,10), SJC = c(0,66), TJC = c(0,68))
  col_map <- c(Pain = Pain, Patient_Global = Patient_Global, SJC = SJC, TJC = TJC)

  oob_mask   <- rep(FALSE, nrow(dat))
  oob_detail <- list()

  for (nm in names(col_map)) {
    vname <- col_map[[nm]]
    v     <- dat[[vname]]
    if (length(v) == 0L) next
    bad <- !is.na(v) & (v < ranges[[nm]][1] | v > ranges[[nm]][2])
    if (any(bad)) {
      oob_detail[[vname]] <- which(bad)
      oob_mask <- oob_mask | bad
    }
  }

  if (any(oob_mask)) {
    msgs <- vapply(names(oob_detail), function(vn) {
      idx <- oob_detail[[vn]]
      sprintf("%s is out of range (%d rows; %s)",
              vn, length(idx), paste(utils::head(idx, 5), collapse = ", "))
    }, character(1))

    if (oob_action == "stop") {
      stop(paste(c("Out-of-range values detected:", paste0(" - ", msgs)), collapse = "\n"), call. = FALSE)
    } else if (oob_action == "na") {
      warning(paste(c("Out-of-range values detected; STC set to NA for affected rows:",
                      paste0(" - ", msgs)), collapse = "\n"), call. = FALSE)
      # keep rows; blank out after computing STC
    } else if (oob_action == "drop") {
      warning(paste(c("Out-of-range rows dropped from output:",
                      paste0(" - ", msgs)), collapse = "\n"), call. = FALSE)
      dat <- dat[!oob_mask, , drop = FALSE]
      oob_mask <- rep(FALSE, nrow(dat))  # reset after dropping
    }
  }

  # -- cDAPSA presence/consistency ----------------------------------------
  cdapsa_sum <- dat[[SJC]] + dat[[TJC]] + dat[[Pain]] + dat[[Patient_Global]]

  if (is.null(cDAPSA)) {
    cDAPSA <- "..cdapsa_calc.."
    dat[[cDAPSA]] <- cdapsa_sum
  } else {
    if (!(cDAPSA %in% names(dat))) {
      stop(sprintf("Column '%s' (cDAPSA) not found in data.", cDAPSA), call. = FALSE)
    }
    both_ok <- !is.na(dat[[cDAPSA]]) & !is.na(cdapsa_sum)
    diffs   <- abs(dat[[cDAPSA]][both_ok] - cdapsa_sum[both_ok])
    if (length(diffs) && any(diffs > cdapsa_tolerance)) {
      bad_idx <- which(both_ok)[which(diffs > cdapsa_tolerance)]
      msg <- sprintf(
        "Provided cDAPSA differs from (SJC+TJC+Pain+Patient_Global) in %d row(s); e.g., indices: %s.",
        length(bad_idx), paste(utils::head(bad_idx, 5), collapse = ", ")
      )
      stop(msg, call. = FALSE)
    }
  }

  # -- NA warning (after coercion; include cDAPSA if present) -------------
  check_vars <- unique(c(cohort_id, cDAPSA, Pain, Patient_Global, SJC, TJC))
  col_na <- vapply(dat[, check_vars, drop = FALSE], function(x) any(is.na(x)), logical(1))
  if (any(col_na)) {
    warning(paste("Missing values detected in columns:", paste(names(col_na)[col_na], collapse = ", ")), call. = FALSE)
  }

  if (nrow(dat) == 0L) {
    out <- dat[, cohort_id, drop = FALSE]
    out$STC <- numeric(0)
    out$STC_quartile <- factor(character(0), levels = c("Q1","Q2","Q3","Q4"))
    out$STC_quartile_combine <- factor(character(0), levels = c("Q1&4","Q2","Q3"))
    return(out)
  }

  # -- Step 1: standardize inputs -----------------------------------------
  dat <- dplyr::mutate(
    dat,
    cDAPSA_std         = ( .data[[cDAPSA]]         - center_scale$cDAPSA["center"])         / center_scale$cDAPSA["scale"],
    Pain_std           = ( .data[[Pain]]           - center_scale$Pain["center"])           / center_scale$Pain["scale"],
    Patient_Global_std = ( .data[[Patient_Global]] - center_scale$Patient_Global["center"]) / center_scale$Patient_Global["scale"],
    SJC_std            = ( .data[[SJC]]            - center_scale$SJC["center"])            / center_scale$SJC["scale"],
    TJC_std            = ( .data[[TJC]]            - center_scale$TJC["center"])            / center_scale$TJC["scale"]
  )

  # -- Step 2: spline basis for standardized cDAPSA -----------------------
  ns_cDAPSA <- splines::ns(dat$cDAPSA_std, df = 4,
                           knots = ns_knots,
                           Boundary.knots = ns_boundary_knots)
  ns_cDAPSA <- as.data.frame(ns_cDAPSA)
  names(ns_cDAPSA) <- paste0("ns", 1:4)
  dat <- dplyr::bind_cols(dat, ns_cDAPSA)

  # -- Step 3: predicted values for each component ------------------------
  for (comp in c("Pain", "Patient_Global", "SJC", "TJC")) {
    coefs <- coef_list[[comp]]
    dat[[paste0(comp, "_yhat")]] <- coefs[1] + coefs[2]*dat$ns1 + coefs[3]*dat$ns2 + coefs[4]*dat$ns3 + coefs[5]*dat$ns4
  }

  # -- Step 4: residuals ---------------------------------------------------
  dat <- dplyr::mutate(
    dat,
    Pain_resid            = .data$Pain_std           - .data$Pain_yhat,
    Patient_Global_resid  = .data$Patient_Global_std - .data$Patient_Global_yhat,
    SJC_resid             = .data$SJC_std            - .data$SJC_yhat,
    TJC_resid             = .data$TJC_std            - .data$TJC_yhat
  )

  # -- Step 5: standardize residuals --------------------------------------
  dat <- dplyr::mutate(
    dat,
    Pain_resid_std = ( .data$Pain_resid - resid_center_scale$center["Pain"]) / resid_center_scale$scale["Pain"],
    Patient_Global_resid_std =
      ( .data$Patient_Global_resid - resid_center_scale$center["Patient Global"]) / resid_center_scale$scale["Patient Global"],
    SJC_resid_std  = ( .data$SJC_resid  - resid_center_scale$center["Swollen Joint Count"]) / resid_center_scale$scale["Swollen Joint Count"],
    TJC_resid_std  = ( .data$TJC_resid  - resid_center_scale$center["Tender Joint Count"])  / resid_center_scale$scale["Tender Joint Count"]
  )

  # -- Step 6: STC + quartiles + combined quartiles -----------------------
  dat <- dplyr::mutate(
    dat,
    STC = .data$Pain_resid_std * stc_loadings[1] +
      .data$Patient_Global_resid_std * stc_loadings[2] +
      .data$SJC_resid_std * stc_loadings[3] +
      .data$TJC_resid_std * stc_loadings[4],
    STC_quartile = cut(.data$STC, breaks = stc_cutoffs, include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4")),
    STC_quartile_combine = dplyr::case_when(
      is.na(.data$STC_quartile)              ~ NA_character_,
      .data$STC_quartile %in% c("Q1","Q4")   ~ "Q1&4",
      TRUE                                   ~ as.character(.data$STC_quartile)
    ),
    STC_quartile_combine = factor(.data$STC_quartile_combine, levels = c("Q1&4","Q2","Q3"))
  )

  # -- blank out STC for OOB rows if requested ----------------------------
  if (oob_action == "na" && any(oob_mask)) {
    dat$STC[oob_mask] <- NA_real_
    dat$STC_quartile[oob_mask] <- NA
    dat$STC_quartile_combine[oob_mask] <- NA
  }

  dplyr::select(dat, dplyr::all_of(c(cohort_id, "STC", "STC_quartile", "STC_quartile_combine")))
}


