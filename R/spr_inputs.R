#' Length-at-age
#'
#' Convert a user-specified length-at-age specification into an age-specific vector
#' spanning the full age range \code{ages = a_min:a_max}.
#'
#' Supported inputs:
#' \itemize{
#'   \item \code{type = "vb"}: von Bertalanffy growth with parameters \code{Linf}, \code{k}, \code{t0}.
#'   \item \code{type = "vector"}: external length-at-age vector \code{L} with \code{length(L) == length(ages)}.
#' }
#'
#' Von Bertalanffy form:
#' \deqn{L_a = L_\infty \left(1 - \exp\left[-k(a - t_0)\right]\right)}
#'
#' @param ages Numeric vector of ages (e.g., \code{a_min:a_max}).
#' @param len_at_age A named list specifying length-at-age. Required element \code{type}.
#'   For \code{type="vb"}, provide \code{Linf}, \code{k}, and \code{t0}. For \code{type="vector"}, provide \code{L}.
#'
#' @return Numeric vector \code{L_a} with \code{length(ages)}.
#'
#' @keywords internal
calc_len_at_age <- function(ages, len_at_age) {
  if (!is.numeric(ages) || length(ages) < 2) stop("ages must be a numeric vector with length >= 2.")
  if (!is.list(len_at_age) || is.null(len_at_age$type)) stop("len_at_age must be a list with element 'type'.")

  if (len_at_age$type == "vb") {
    req <- c("Linf", "k", "t0")
    miss <- req[!req %in% names(len_at_age)]
    if (length(miss) > 0) stop("len_at_age$type='vb' requires: ", paste(miss, collapse = ", "))
    Linf <- as.numeric(len_at_age$Linf)
    k    <- as.numeric(len_at_age$k)
    t0   <- as.numeric(len_at_age$t0)
    if (any(!is.finite(c(Linf, k, t0)))) stop("len_at_age vb parameters must be finite.")
    L <- Linf * (1 - exp(-k * (ages - t0)))
    if (any(!is.finite(L)) || any(L < 0)) stop("Computed length-at-age contains non-finite or negative values.")
    return(as.numeric(L))
  }

  if (len_at_age$type == "vector") {
    if (is.null(len_at_age$L)) stop("len_at_age$type='vector' requires len_at_age$L.")
    L <- as.numeric(len_at_age$L)
    if (length(L) != length(ages)) stop("len_at_age$L must have length equal to length(ages).")
    if (any(!is.finite(L)) || any(L < 0)) stop("len_at_age$L must be finite and nonnegative.")
    return(as.numeric(L))
  }

  stop("Unknown len_at_age$type. Use 'vb' or 'vector'.")
}

#' Length-at-age
#'
#' Build a length-at-age vector from a length specification.
#'
#' This is a user-facing wrapper around an internal helper used by
#' [spr_input()]. It is useful for sensitivity analyses where the user wants to
#' modify the length-at-age specification of an existing [spr_input()] object.
#'
#' @inheritParams calc_len_at_age
#' @return Numeric vector of length \code{length(ages)}.
#' @export
len_at_age <- function(ages, len_at_age) {
  calc_len_at_age(ages = ages, len_at_age = len_at_age)
}

#' Weight-at-age
#'
#' Convert a user-specified weight-at-age specification into an age-specific vector
#' spanning the full age range \code{ages = a_min:a_max}.
#'
#' Supported inputs:
#' \itemize{
#'   \item \code{type = "wl"}: derive weight-at-age from length-at-age using \code{alpha} and \code{beta}.
#'   \item \code{type = "vector"}: external weight-at-age vector \code{W} with \code{length(W) == length(ages)}.
#' }
#'
#' Weight-length form:
#' \deqn{w_a = \alpha L_a^\beta}
#'
#' @param ages Numeric vector of ages (e.g., \code{a_min:a_max}).
#' @param wt_at_age A named list specifying weight-at-age. Required element \code{type}.
#'   For \code{type="wl"}, provide \code{alpha} and \code{beta}. For \code{type="vector"}, provide \code{W}.
#' @param len_at_age Optional numeric vector \code{L_a} with \code{length(ages)}. Required when \code{wt_at_age$type = "wl"}.
#'
#' @return Numeric vector \code{w_a} with \code{length(ages)}.
#'
#' @keywords internal
calc_wt_at_age <- function(ages, wt_at_age, len_at_age = NULL) {
  if (!is.numeric(ages) || length(ages) < 2) stop("ages must be a numeric vector with length >= 2.")
  if (!is.list(wt_at_age) || is.null(wt_at_age$type)) stop("wt_at_age must be a list with element 'type'.")

  if (wt_at_age$type == "wl") {
    req <- c("alpha", "beta")
    miss <- req[!req %in% names(wt_at_age)]
    if (length(miss) > 0) stop("wt_at_age$type='wl' requires: ", paste(miss, collapse = ", "))
    if (is.null(len_at_age)) stop("wt_at_age$type='wl' requires len_at_age to be provided.")
    L <- as.numeric(len_at_age)
    if (length(L) != length(ages)) stop("len_at_age must have length equal to length(ages).")
    if (any(!is.finite(L)) || any(L < 0)) stop("len_at_age must be finite and nonnegative.")
    alpha <- as.numeric(wt_at_age$alpha)
    beta  <- as.numeric(wt_at_age$beta)
    if (any(!is.finite(c(alpha, beta)))) stop("wt_at_age wl parameters must be finite.")
    W <- alpha * L^beta
    if (any(!is.finite(W)) || any(W < 0)) stop("Computed weight-at-age contains non-finite or negative values.")
    return(as.numeric(W))
  }

  if (wt_at_age$type == "vector") {
    if (is.null(wt_at_age$W)) stop("wt_at_age$type='vector' requires wt_at_age$W.")
    W <- as.numeric(wt_at_age$W)
    if (length(W) != length(ages)) stop("wt_at_age$W must have length equal to length(ages).")
    if (any(!is.finite(W)) || any(W < 0)) stop("wt_at_age$W must be finite and nonnegative.")
    return(as.numeric(W))
  }

  stop("Unknown wt_at_age$type. Use 'wl' or 'vector'.")
}

#' Weight-at-age
#'
#' Build a weight-at-age vector from a weight specification.
#'
#' This is a user-facing wrapper around an internal helper used by
#' [spr_input()]. It is useful for sensitivity analyses where the user wants to
#' modify the weight specification of an existing [spr_input()] object.
#'
#' Use \code{len_at_age} when \code{wt_at_age$type == "wl"}.
#'
#' @inheritParams calc_wt_at_age
#' @return Numeric vector of length \code{length(ages)}.
#' @export
wt_at_age <- function(ages, wt_at_age, len_at_age = NULL) {
  calc_wt_at_age(ages = ages, wt_at_age = wt_at_age, len_at_age = len_at_age)
}

#' Maturity-at-age
#'
#' Convert a user-specified maturity specification into an age-specific vector
#' spanning the full age range \code{ages = a_min:a_max}.
#'
#' Supported inputs:
#' \itemize{
#'   \item \code{type = "vector"}: external maturity vector \code{m} with \code{length(m) == length(ages)}.
#'   \item \code{type = "logistic"}: logistic curve with parameters \code{a50} and \code{delta}.
#' }
#'
#' Logistic form:
#' \deqn{m_a = \frac{1}{1 + \exp\left(\frac{a - a_{50}}{\delta}\right)}}
#'
#' @param ages Numeric vector of ages (e.g., \code{a_min:a_max}).
#' @param maturity A named list specifying maturity. Required element \code{type}.
#'   For \code{type="vector"}, provide \code{m}. For \code{type="logistic"}, provide \code{a50} and \code{delta}.
#'
#' @return Numeric vector \code{maturity_at_age} with \code{length(ages)}.
#'
#' @keywords internal
calc_maturity_at_age <- function(ages, maturity) {
  stopifnot(is.numeric(ages), length(ages) >= 2, is.list(maturity), !is.null(maturity$type))

  if (maturity$type == "vector") {
    if (is.null(maturity$m)) stop("For maturity$type = 'vector', provide maturity$m.")
    m <- as.numeric(maturity$m)
    if (length(m) != length(ages)) stop("maturity$m must have length equal to length(ages).")
    if (any(!is.finite(m))) stop("maturity$m contains non-finite values.")
    if (any(m < 0 | m > 1, na.rm = TRUE)) warning("maturity$m had values outside [0,1].")
    return(pmin(1, pmax(0, m)))
  }

  if (maturity$type == "logistic") {
    if (is.null(maturity$a50) || is.null(maturity$delta)) {
      stop("For maturity$type = 'logistic', provide maturity$a50 and maturity$delta.")
    }
    a50 <- as.numeric(maturity$a50)
    delta <- as.numeric(maturity$delta)
    if (!is.finite(a50) || !is.finite(delta)) stop("maturity$a50 and maturity$delta must be finite.")
    if (abs(delta) < 1e-12) stop("maturity$delta is too close to zero.")
    if (delta > 0) warning("maturity$delta > 0 yields a decreasing maturity ogive with age. User beware!")
    m <- 1 / (1 + exp((ages - a50) / delta))
    return(pmin(1, pmax(0, as.numeric(m))))
  }

  stop("Unknown maturity$type. Use 'vector' or 'logistic'.")
}

#' Maturity-at-age
#'
#' Build a maturity-at-age vector from a maturity specification.
#'
#' This is a user-facing wrapper around an internal helper used by
#' [spr_input()]. It is useful for sensitivity analyses where the user wants to
#' modify the maturity specification of an existing [spr_input()] object.
#'
#' It returns maturity bounded to \eqn{[0,1]}.
#'
#' @inheritParams calc_maturity_at_age
#' @return Numeric vector of length \code{length(ages)} with values in \eqn{[0,1]}.
#' @export
maturity_at_age <- function(ages, maturity) {
  calc_maturity_at_age(ages = ages, maturity = maturity)
}

#' Fishery selectivity-at-age
#'
#' Convert a user-specified selectivity specification into an age-specific vector
#' spanning the full age range \code{ages = a_min:a_max}. Selectivity is
#' bounded to \eqn{[0, 1]} and may optionally be scaled so that \code{max(s_a) = 1}.
#'
#' Supported inputs:
#' \itemize{
#'   \item \code{type = "vector"}: external selectivity vector \code{s} with \code{length(s) == length(ages)}.
#'   \item \code{type = "knife_edge"}: step function with age-at-full-selection \code{a_c}.
#'   \item \code{type = "logistic"}: logistic curve with parameters \code{a50} and \code{delta}.
#' }
#'
#' Knife-edge form:
#' \deqn{s_a = 0 \text{ for } a < a_c,\quad s_a = 1 \text{ for } a \ge a_c}
#'
#' Logistic form:
#' \deqn{s_a = \frac{1}{1 + \exp\left(\frac{a - a_{50}}{\delta}\right)}}
#'
#' @param ages Numeric vector of ages (e.g., \code{a_min:a_max}).
#' @param selectivity A named list specifying selectivity. Required element \code{type}.
#'   For \code{type="vector"}, provide \code{s}. For \code{type="knife_edge"}, provide \code{a_c}.
#'   For \code{type="logistic"}, provide \code{a50} and \code{delta}.
#' @param scale_max Logical. If \code{TRUE}, scales selectivity so that \code{max(s_a) = 1} when possible.
#'
#' @return Numeric vector \code{selectivity_at_age} with \code{length(ages)}.
#'
#' @keywords internal
calc_selectivity_at_age <- function(ages, selectivity, scale_max = TRUE) {
  stopifnot(is.numeric(ages), length(ages) >= 2, is.list(selectivity), !is.null(selectivity$type))

  if (selectivity$type == "vector") {
    if (is.null(selectivity$s)) stop("For selectivity$type = 'vector', provide selectivity$s.")
    s <- as.numeric(selectivity$s)
    if (length(s) != length(ages)) stop("selectivity$s must have length equal to length(ages).")
    if (any(!is.finite(s))) stop("selectivity$s contains non-finite values.")
    if (any(s < 0, na.rm = TRUE)) warning("selectivity$s had values < 0.")
    s <- pmax(0, s)
    if (scale_max) {
      smax <- max(s, na.rm = TRUE)
      if (smax > 0) s <- s / smax
    }
    return(pmin(1, as.numeric(s)))
  }

  if (selectivity$type == "knife_edge") {
    if (is.null(selectivity$a_c)) stop("For selectivity$type = 'knife_edge', provide selectivity$a_c.")
    a_c <- as.numeric(selectivity$a_c)
    if (!is.finite(a_c)) stop("selectivity$a_c must be finite.")
    return(as.numeric(ifelse(ages >= a_c, 1, 0)))
  }

  if (selectivity$type == "logistic") {
    if (is.null(selectivity$a50) || is.null(selectivity$delta)) {
      stop("For selectivity$type = 'logistic', provide selectivity$a50 and selectivity$delta.")
    }
    a50 <- as.numeric(selectivity$a50)
    delta <- as.numeric(selectivity$delta)
    if (!is.finite(a50) || !is.finite(delta)) stop("selectivity$a50 and selectivity$delta must be finite.")
    if (abs(delta) < 1e-12) stop("selectivity$delta is too close to zero.")
    s <- 1 / (1 + exp((ages - a50) / delta))
    s <- pmin(1, pmax(0, as.numeric(s)))
    if (scale_max) {
      smax <- max(s, na.rm = TRUE)
      if (smax > 0) s <- s / smax
      s <- pmin(1, s)
    }
    return(as.numeric(s))
  }

  stop("Unknown selectivity$type. Use 'vector', 'knife_edge', or 'logistic'.")
}

#' Fishery selectivity-at-age
#'
#' Build a selectivity-at-age vector from a selectivity specification.
#'
#' This is a user-facing wrapper around an internal helper used by
#' [spr_input()]. It is useful for sensitivity analyses where the user wants to
#' modify the fishery selectivity specification of an existing [spr_input()]
#' object.
#'
#' It returns selectivity bounded to \eqn{[0,1]} and can optionally scale so that
#' \code{max(s_a) = 1}.
#'
#' @inheritParams calc_selectivity_at_age
#' @return Numeric vector of length \code{length(ages)} with values in \eqn{[0,1]}.
#'   When \code{selectivity$type == "fleets"}, the returned vector is the effective
#'   selectivity-at-age.
#' @export
selectivity_at_age <- function(ages, selectivity, scale_max = TRUE) {
  calc_selectivity_at_age(ages = ages, selectivity = selectivity, scale_max = scale_max)
}

#' Construct SPR input object
#'
#' Validates and expand inputs needed for per-recruit SPR/YPR calculations for
#' one or multiple species.
#'
#' For a single-species analysis, provide a single species specification via
#' `species` (can be named or unnamed). For multispecies analyses, provide a
#' named list of species specifications, where each element defines life history
#' and fishery inputs for one species.
#'
#' Each species specification must include:
#' \itemize{
#'   \item \code{len_at_age}: list with \code{type="vb"} (Linf, k, t0) or \code{type="vector"} (L).
#'   \item \code{wt_at_age}: list with \code{type="wl"} (alpha, beta) or \code{type="vector"} (W).
#'   \item \code{maturity}: list with \code{type="logistic"} (a50, delta) or \code{type="vector"} (m).
#'   \item \code{selectivity}: list with \code{type="logistic"} (a50, delta), \code{type="knife_edge"} (a_c),
#'     or \code{type="vector"} (s).
#'   \item \code{M}: numeric scalar or vector of length \code{length(ages)}.
#' }
#'
#' Rules for size-at-age:
#' \itemize{
#'   \item If \code{wt_at_age$type = "vector"}, the provided weight-at-age vector is used directly in the SPR model.
#'     In this case, \code{len_at_age} is not used to determine size-at-age for the SPR calculations, and a message
#'     is emitted to make this explicit.
#'   \item If \code{len_at_age$type = "vb"}, then \code{wt_at_age$type} must be \code{"wl"} to ensure consistency
#'     between growth and weight conversion.
#' }
#'
#' @param ages Numeric vector of ages (e.g., \code{a_min:a_max}).
#' @param species A named or unnamed list for single species specification, or a
#'   named list of species for multispecies specifications. See details.
#' @param rec_prop Optional numeric vector of recruitment proportions for multispecies analyses.
#'   Must be the same length as the number of species and sum to 1. Defaults to equal proportions.
#' @param R0 Numeric scalar recruitment scaling used for per-recruit
#'   calculations (default=1, or equal proportions if inputting multiple
#'   species). If defining R0 for multiple species, R0 proportions must be a
#'   vector equal to the number of species and sum to 1.
#' @param use_plus_group Logical. If \code{TRUE} (default), the terminal age will be treated as a plus group in pop dy equations.
#' @param scale_selex Logical. If \code{TRUE} (default), scales selectivity so max equals 1 when possible.
#'
#' @return A named list with expanded, validated inputs for running SPR functions:
#' \itemize{
#'   \item \code{ages}, \code{R0}, \code{use_plus_group}
#'   \item \code{species}: named list of per-species expanded vectors and original specifications
#'   \item \code{rec_prop}: recruitment proportions used for multispecies calculations
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' ages <- 3:52
#' inp1 <- spr_input(
#'   ages = ages,
#'   species = list(
#'   BS = list(
#'     len_at_age   = list(type="vb", Linf=52.6, k=0.062, t0=0.2),
#'     wt_at_age    = list(type="wl", alpha=8.13e-6, beta=3.177),
#'     maturity     = list(type="logistic", a50=20, delta=-2),
#'     selectivity  = list(type="knife_edge", a_c=10),
#'     M = 0.042
#'     )
#' )
#'
#' inp2 <- spr_input(
#'   ages = ages,
#'   species = list(
#'     RE = list(
#'       len_at_age  = list(type="vb", Linf=54.4, k=0.102, t0=-0.36),
#'       wt_at_age   = list(type="wl", alpha=1.143e-5, beta=3.098),
#'       maturity    = list(type="vector", m=rep(0, length(ages))), # placeholder
#'       selectivity = list(type="vector", s=rep(1, length(ages))),
#'       M = 0.042
#'     ),
#'     BS = list(
#'       len_at_age  = list(type="vb", Linf=52.6, k=0.062, t0=0.2),
#'       wt_at_age   = list(type="wl", alpha=8.13e-6, beta=3.177),
#'       maturity    = list(type="vector", m=rep(0, length(ages))), # placeholder
#'       selectivity = list(type="vector", s=rep(1, length(ages))),
#'       M = 0.042
#'     )
#'   ),
#'   rec_prop = c(RE=0.5, BS=0.5)
#' )
#' }
spr_input <- function(
    ages,
    species,
    rec_prop = NULL,
    R0 = 1,
    use_plus_group = TRUE,
    scale_selex = TRUE
) {
  if (!is.numeric(ages) || length(ages) < 2) stop("ages must be a numeric vector with length >= 2.")
  if (any(!is.finite(ages))) stop("ages contains non-finite values.")
  if (any(diff(ages) <= 0)) stop("ages must be strictly increasing (e.g., a_min:a_max).")

  if (!is.numeric(R0) || length(R0) != 1 || !is.finite(R0) || R0 <= 0) stop("R0 must be a single positive number.")
  if (!isTRUE(use_plus_group) && !isFALSE(use_plus_group)) stop("use_plus_group must be TRUE or FALSE.")
  if (!isTRUE(scale_selex) && !isFALSE(scale_selex)) stop("scale_selex must be TRUE or FALSE.")

  # normalize species input: single list -> named list of length 1
  if (!is.list(species)) stop("species must be a list (single species) or a named list of species lists.")
  is_single_species <- !is.null(species$len_at_age) || !is.null(species$wt_at_age) || !is.null(species$maturity)
  if (is_single_species) {
    species <- list(sp1 = species)
  } else {
    if (is.null(names(species)) || any(names(species) == "")) {
      stop("For multispecies, provide a *named* list of species (e.g., list(RE=..., BS=...)).")
    }
  }

  nsp <- length(species)

  # recruitment proportions
  if (is.null(rec_prop)) {
    rec_prop <- rep(1 / nsp, nsp)
    names(rec_prop) <- names(species)
  } else {
    if (!is.numeric(rec_prop) || any(!is.finite(rec_prop))) stop("rec_prop must be numeric and finite.")
    if (length(rec_prop) != nsp) stop("rec_prop length must equal the number of species.")
    if (any(rec_prop < 0)) stop("rec_prop must be nonnegative.")
    if (is.null(names(rec_prop))) names(rec_prop) <- names(species)
    s <- sum(rec_prop)
    if (abs(s - 1) > 1e-8) stop("rec_prop must sum to 1.")
  }

  # validate and expand one species at a time
  expand_one_species <- function(sp, sp_name) {
    req <- c("len_at_age", "wt_at_age", "maturity", "selectivity", "M")
    miss <- req[!req %in% names(sp)]
    if (length(miss) > 0) stop("Species '", sp_name, "' is missing: ", paste(miss, collapse = ", "))

    # M: scalar or vector across ages
    M <- sp$M
    if (!is.numeric(M) || any(!is.finite(M))) stop("Species '", sp_name, "': M must be numeric and finite.")
    if (length(M) == 1) M_at_age <- rep(as.numeric(M), length(ages))
    else if (length(M) == length(ages)) M_at_age <- as.numeric(M)
    else stop("Species '", sp_name, "': M must be a scalar or length(ages).")
    if (any(M_at_age < 0)) stop("Species '", sp_name, "': M must be nonnegative.")

    # If len_at_age is vb, require wt_at_age is wl
    if (is.list(sp$len_at_age) && !is.null(sp$len_at_age$type) &&
        sp$len_at_age$type == "vb") {
      if (!is.list(sp$wt_at_age) || is.null(sp$wt_at_age$type) || sp$wt_at_age$type != "wl") {
        stop("Species '", sp_name,
             "': len_at_age$type='vb' requires wt_at_age$type='wl' (alpha, beta) for consistency.")
      }
    }

    # Expand weight first if vector override requested
    if (is.list(sp$wt_at_age) && !is.null(sp$wt_at_age$type) && sp$wt_at_age$type == "vector") {
      W <- calc_wt_at_age(ages, sp$wt_at_age, len_at_age = NULL)
      message("Species '", sp_name,
              "': wt_at_age$type='vector' provided, using weight-at-age directly. ",
              "len_at_age will not be used to parameterize size-at-age in the SPR model.")
      L <- rep(NA_real_, length(ages))
    } else {
      L <- calc_len_at_age(ages, sp$len_at_age)
      W <- calc_wt_at_age(ages, sp$wt_at_age, len_at_age = L)
    }

    mat <- calc_maturity_at_age(ages, sp$maturity)
    sel <- calc_selectivity_at_age(ages, sp$selectivity, scale_max = scale_selex)

    list(
      name = sp_name,


      R0_base = as.numeric(R0), # store baseline and scaled recruitment
      R0 = as.numeric(R0 * rec_prop[sp_name]),  # per-species recruitment scaling for multispecies

      M_at_age = M_at_age,
      len_at_age = as.numeric(L),
      wt_at_age = as.numeric(W),
      maturity_at_age = as.numeric(mat),
      selex_at_age = as.numeric(sel),

      # store maturity parameterization for sensitivity analyses
      maturity_delta = if (is.list(sp$maturity) && identical(sp$maturity$type, "logistic")) {
        as.numeric(sp$maturity$delta)
      } else {
        NA_real_
      },

      maturity_spec = sp$maturity,

      spec = sp
    )
  }

  out_species <- lapply(names(species), function(nm) expand_one_species(species[[nm]], nm))
  names(out_species) <- names(species)

  list(
    ages = as.numeric(ages),
    R0 = as.numeric(R0),
    use_plus_group = isTRUE(use_plus_group),
    rec_prop = rec_prop,
    species = out_species
  )
}

#' Plot SPR input assumptions
#'
#' Visualize inputs supplied to spr_input() object for checking assumptions and
#' reporting (base vs sensitivity, or across species). By default, makes:
#'  - maturity and selectivity on the same panel
#'  - weight-at-age
#' Optional panels include length-at-age and M-at-age.
#'
#' @param x An object returned by spr_input().
#' @param species Character vector of species names to plot. Defaults to all in x$species.
#' @param compare Optional named list of additional spr_input objects (e.g.,
#'   sensitivities). Each element name is used as the scenario label in
#'   legends/facets. Example: compare = list(Base = inp_base, Sens1 = inp_sens)
#' @param panels Character vector specifying which panels to draw. Any of:
#'   c("mat_selex","wt","len","M")
#' @param overlay_maturity_selex Logical. If TRUE, maturity and selectivity are overlayed
#'   on a single plot (recommended).
#' @param return_data Logical. If TRUE, returns a list with plots and the tidy data used.
#'
#' @return By default, returns a named list of ggplot objects. If return_data=TRUE,
#'   returns list(plots=..., data=...).
#' @export
plot_spr_inputs <- function(
    x,
    species = NULL,
    compare = NULL,
    panels = c("mat_selex", "wt", "len", "M"),
    overlay_maturity_selex = TRUE,
    return_data = FALSE
  ) {
    if (!is.list(x) || is.null(x$ages) || is.null(x$species)) {
      stop("x must be an object returned by spr_input().")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2.")
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install tidyr.")
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr.")
    if (!requireNamespace("scales", quietly = TRUE)) stop("Please install scales.")

    facet_by_scenario <- !is.null(compare)

    if (facet_by_scenario) {
      scenario_list <- c(list(Base = x), compare)
      if (is.null(names(scenario_list)) || any(names(scenario_list) == "")) {
        stop("compare must be a named list, e.g., list(Sens1=..., Sens2=...).")
      }
    } else {
      scenario_list <- list(x)
      names(scenario_list) <- "scenario"
    }

    all_species <- names(x$species)
    if (is.null(species)) species <- all_species
    if (any(!species %in% all_species)) {
      stop("species must be a subset of names(x$species): ", paste(all_species, collapse = ", "))
    }

    build_df_one <- function(inp, scenario_nm) {
      ages <- inp$ages
      out <- lapply(species, function(sp) {
        spdat <- inp$species[[sp]]
        data.frame(
          scenario = scenario_nm,
          species  = sp,
          age      = ages,
          maturity = spdat$maturity_at_age,
          selex    = spdat$selex_at_age,
          weight   = spdat$wt_at_age,
          length   = spdat$len_at_age,
          M        = spdat$M_at_age,
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, out)
    }

    if (facet_by_scenario) {
      df <- do.call(rbind, Map(build_df_one, scenario_list, names(scenario_list)))
    } else {
      df <- build_df_one(x, scenario_nm = "single")
    }

    plots <- list()

    facet_if_needed <- function(p) {
      if (facet_by_scenario) p + ggplot2::facet_wrap(~ scenario, ncol = 1) else p
    }

    single_species_mode <- length(species) == 1

    # helper to deal with multispp cases where one or more inputs are shared
    tol <- 1e-12

    collapse_shared_wide <- function(df_wide, value_col, shared_label) {
      # Within each scenario-age, if value is identical across species, keep one row
      # and relabel species as "Shared ...", so the plot can show it once.
      df_wide |>
        dplyr::group_by(.data$scenario, .data$age) |>
        dplyr::mutate(
          .range = max(.data[[value_col]], na.rm = TRUE) - min(.data[[value_col]], na.rm = TRUE),
          .shared = .data$.range <= tol
        ) |>
        dplyr::ungroup() |>
        dplyr::group_by(.data$scenario, .data$age) |>
        dplyr::arrange(.data$species, .by_group = TRUE) |>
        dplyr::filter(!.data$.shared | dplyr::row_number() == 1) |>
        dplyr::ungroup() |>
        dplyr::mutate(species = ifelse(.data$.shared, shared_label, .data$species)) |>
        dplyr::select(-.data$.range, -.data$.shared)
    }

    # scenario, species, age, quantity, and value
    collapse_shared_long <- function(df_long, shared_prefix = "Shared") {
      out <- lapply(unique(df_long$quantity), function(q) {
        sub <- df_long[df_long$quantity == q, , drop = FALSE]
        collapse_shared_wide(
          df_wide = sub,
          value_col = "value",
          shared_label = paste0(shared_prefix, " ", q)
        )
      })
      do.call(rbind, out)
    }

    make_color_scale <- function(df_plot, shared_prefix = "Shared") {
      lev <- unique(df_plot$species)
      shared <- lev[grepl(paste0("^", shared_prefix), lev)]
      nonshared <- setdiff(lev, shared)

      pal <- scales::hue_pal()(max(1, length(nonshared)))
      vals <- setNames(pal[seq_along(nonshared)], nonshared)
      if (length(shared) > 0) {
        vals <- c(vals, setNames(rep("black", length(shared)), shared))
      }

      ggplot2::scale_color_manual(values = vals, breaks = lev, drop = FALSE)
    }

    # ---- Panel 1: maturity + selectivity ----
    if ("mat_selex" %in% panels) {
      if (isTRUE(overlay_maturity_selex)) {
        df_long <- tidyr::pivot_longer(
          df,
          cols = c("maturity", "selex"),
          names_to = "quantity",
          values_to = "value"
        )

        if (isTRUE(single_species_mode)) {
          p <- ggplot2::ggplot(
            df_long,
            ggplot2::aes(x = age, y = value, linetype = quantity, group = quantity)
          ) +
            ggplot2::geom_line(linewidth = 1) +
            ggplot2::scale_linetype_manual(values = c(maturity = 1, selex = 3)) +
            ggplot2::labs(x = "Age", y = "Maturity or selectivity", linetype = NULL) +
            ggplot2::coord_cartesian(ylim = c(0, 1))
          p <- facet_if_needed(p)
          plots$mat_selex <- p + theme_tier4()
        } else {
          # collapse shared lines independently for maturity and selex
          df_long2 <- collapse_shared_long(df_long, shared_prefix = "Shared")

          p <- ggplot2::ggplot(
            df_long2,
            ggplot2::aes(x = age, y = value, color = species, linetype = quantity)
          ) +
            ggplot2::geom_line(linewidth = 1) +
            ggplot2::scale_linetype_manual(values = c(maturity = 1, selex = 3)) +
            ggplot2::labs(
              x = "Age",
              y = "Maturity or selectivity",
              color = "Species",
              linetype = NULL
            ) +
            ggplot2::coord_cartesian(ylim = c(0, 1)) +
            make_color_scale(df_long2, shared_prefix = "Shared")

          p <- facet_if_needed(p)
          plots$mat_selex <- p + scale_color_tier4() + scale_fill_tier4() + theme_tier4()
        }
      } else {
        # separate plots for maturity and selectivity
        if (isTRUE(single_species_mode)) {
          p_mat <- ggplot2::ggplot(df, ggplot2::aes(age, maturity)) +
            ggplot2::geom_line(linewidth = 1) +
            ggplot2::labs(x = "Age", y = "Maturity") +
            ggplot2::coord_cartesian(ylim = c(0, 1))

          p_sel <- ggplot2::ggplot(df, ggplot2::aes(age, selex)) +
            ggplot2::geom_line(linewidth = 1) +
            ggplot2::labs(x = "Age", y = "Selectivity") +
            ggplot2::coord_cartesian(ylim = c(0, 1))

          p_mat <- facet_if_needed(p_mat)
          p_sel <- facet_if_needed(p_sel)
          plots$maturity <- p_mat + theme_tier4()
          plots$selectivity <- p_sel + theme_tier4()
        } else {
          df_m <- collapse_shared_wide(df, "maturity", "Shared maturity")
          df_s <- collapse_shared_wide(df, "selex", "Shared selectivity")

          p_mat <- ggplot2::ggplot(df_m, ggplot2::aes(age, maturity, color = species)) +
            ggplot2::geom_line(linewidth = 1) +
            ggplot2::labs(x = "Age", y = "Maturity", color = "Species") +
            ggplot2::coord_cartesian(ylim = c(0, 1)) +
            make_color_scale(df_m, shared_prefix = "Shared")

          p_sel <- ggplot2::ggplot(df_s, ggplot2::aes(age, selex, color = species)) +
            ggplot2::geom_line(linewidth = 1) +
            ggplot2::labs(x = "Age", y = "Selectivity", color = "Species") +
            ggplot2::coord_cartesian(ylim = c(0, 1)) +
            make_color_scale(df_s, shared_prefix = "Shared")

          p_mat <- facet_if_needed(p_mat)
          p_sel <- facet_if_needed(p_sel)
          plots$maturity <- p_mat + scale_color_tier4() + scale_fill_tier4() + theme_tier4()
          plots$selectivity <- p_sel + scale_color_tier4() + scale_fill_tier4() + theme_tier4()
        }
      }
    }

    # ---- Panel 2: weight-at-age ----
    if ("wt" %in% panels) {
      if (isTRUE(single_species_mode)) {
        p <- ggplot2::ggplot(df, ggplot2::aes(age, weight)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::labs(x = "Age", y = "Weight-at-age") +
          theme_tier4()
        if (facet_by_scenario) p <- p + ggplot2::facet_wrap(~ scenario, ncol = 1, scales = "free_y")
        plots$wt <- p
      } else {
        df_w <- collapse_shared_wide(df, "weight", "Shared weight-at-age")
        p <- ggplot2::ggplot(df_w, ggplot2::aes(age, weight, color = species)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::labs(x = "Age", y = "Weight-at-age", color = "Species") +
          make_color_scale(df_w, shared_prefix = "Shared") +
          scale_color_tier4() +
          scale_fill_tier4() +
          theme_tier4()
        if (facet_by_scenario) p <- p + ggplot2::facet_wrap(~ scenario, ncol = 1, scales = "free_y")
        plots$wt <- p
      }
    }

    # ---- Optional Panel 3: length-at-age ----
    if ("len" %in% panels) {
      if (isTRUE(single_species_mode)) {
        p <- ggplot2::ggplot(df, ggplot2::aes(age, length)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::labs(x = "Age", y = "Length-at-age") +
          theme_tier4()
        if (facet_by_scenario) p <- p + ggplot2::facet_wrap(~ scenario, ncol = 1, scales = "free_y")
        plots$len <- p
      } else {
        df_l <- collapse_shared_wide(df, "length", "Shared length-at-age")
        p <- ggplot2::ggplot(df_l, ggplot2::aes(age, length, color = species)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::labs(x = "Age", y = "Length-at-age", color = "Species") +
          make_color_scale(df_l, shared_prefix = "Shared") +
          scale_color_tier4() +
          scale_fill_tier4() +
          theme_tier4()
        if (facet_by_scenario) p <- p + ggplot2::facet_wrap(~ scenario, ncol = 1, scales = "free_y")
        plots$len <- p
      }
    }

    # ---- Optional Panel 4: M-at-age -----
    if ("M" %in% panels) {
      if (isTRUE(single_species_mode)) {
        p <- ggplot2::ggplot(df, ggplot2::aes(age, M)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::labs(x = "Age", y = "Natural mortality (M)") +
          theme_tier4()
        if (facet_by_scenario) p <- p + ggplot2::facet_wrap(~ scenario, ncol = 1, scales = "free_y")
        plots$M <- p
      } else {
        df_M <- collapse_shared_wide(df, "M", "Shared M")
        p <- ggplot2::ggplot(df_M, ggplot2::aes(age, M, color = species)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::labs(x = "Age", y = "Natural mortality (M)", color = "Species") +
          make_color_scale(df_M, shared_prefix = "Shared") +
          scale_color_tier4() +
          scale_fill_tier4() +
          theme_tier4()
        if (facet_by_scenario) p <- p + ggplot2::facet_wrap(~ scenario, ncol = 1, scales = "free_y")
        plots$M <- p
      }
    }

    if (return_data) return(list(plots = plots, data = df))
    plots
  }


#' Arrange SPR input plots into a report-ready grid
#'
#' Combine list output from \code{plot_spr_inputs()} into a multi-panel figure
#' using \pkg{patchwork}. Because \pkg{patchwork} is in Suggests, this function
#' will only run when it is installed.
#'
#' @details
#' Typical usage:
#' \preformatted{
#'   inp_plots <- plot_spr_inputs(inp, panels = c("mat_selex","wt","len","M"))
#'   grid_spr_inputs(inp_plots, order = c("mat_selex","wt","len","M"), ncol = 1)
#' }
#'
#' @param plots A named list of ggplot objects returned by \code{plot_spr_inputs()}.
#' @param order Character vector giving the plot names (in \code{plots}) in the order to arrange.
#'   Defaults to \code{names(plots)}.
#' @param ncol Number of columns in the grid.
#' @param nrow Optional number of rows. If \code{NULL}, patchwork chooses automatically.
#' @param guides Passed to \code{patchwork::plot_layout(guides = ...)}. Default \code{"collect"}
#'   collects legends when possible.
#' @param title Optional figure title (character). Added with \code{patchwork::plot_annotation()}.
#' @param caption Optional figure caption (character). Added with \code{patchwork::plot_annotation()}.
#'
#' @return A patchwork object (a ggplot-like object) that can be printed or saved.
#' @export
grid_spr_inputs <- function(
    plots,
    order = NULL,
    ncol = 1,
    nrow = NULL,
    guides = c("collect", "keep"),
    title = NULL,
    caption = NULL
) {
  guides <- match.arg(guides)

  if (!is.list(plots) || length(plots) == 0) {
    stop("grid_spr_inputs(): 'plots' must be a non-empty list returned by plot_spr_inputs().")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("grid_spr_inputs(): Please install 'patchwork' to arrange plots (it is in Suggests).")
  }

  if (is.null(order)) order <- names(plots)
  if (is.null(order) || any(order == "")) order <- names(plots)

  missing <- setdiff(order, names(plots))
  if (length(missing) > 0) {
    stop("grid_spr_inputs(): requested plot(s) not found: ", paste(missing, collapse = ", "))
  }

  plot_list <- plots[order]

  pw <- patchwork::wrap_plots(plotlist = plot_list, ncol = ncol, nrow = nrow) +
    patchwork::plot_layout(guides = guides)

  if (!is.null(title) || !is.null(caption)) {
    pw <- pw + patchwork::plot_annotation(title = title, caption = caption)
  }

  pw
}
